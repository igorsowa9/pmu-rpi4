#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include "uldaq.h"
#include "utility.h"
#include <time.h>
#include <sys/timeb.h>
#include <sys/time.h>
#include "MQTTClient.h"

#define MAX_DEV_COUNT  100
#define MAX_STR_LENGTH 64
#define MAX_SCAN_OPTIONS_LENGTH 256

#define M_PI 3.14159265358979323846
//#define M_E 2.71828182845904523536

#define ADDRESS     "tcp://localhost:1883"
#define CLIENTID    "ExampleClientPub"
#define TOPIC       "pmu_topic"
#define TEST_PAYL   "Hello World!"
#define QOS         0
#define TIMEOUT     10000L
volatile MQTTClient_deliveryToken deliveredtoken;

void connlost(void *context, char *cause)
{
    printf("\nConnection lost\n");
    printf("     cause: %s\n", cause);
}
void send_phasor(double f, double abs, double ph_deg, float ptt, MQTTClient client, MQTTClient_message pubmsg, MQTTClient_deliveryToken token) {
        char str[50];
        sprintf(str, "%f\t%f\t%f\t%f", abs, ph_deg, f, ptt);
        printf("sending: %s\n", str);
        pubmsg.payload = str;
        pubmsg.payloadlen = strlen(str);//12; //sizeof(double); //strlen(PAYLOAD);
        //MQTTClient_publishMessage(client, TOPIC, &pubmsg, &token);
}
int main(void)
{
	int descriptorIndex = 0;
	DaqDeviceDescriptor devDescriptors[MAX_DEV_COUNT];
	DaqDeviceInterface interfaceType = ANY_IFC;
	DaqDeviceHandle daqDeviceHandle = 0;
	unsigned int numDevs = MAX_DEV_COUNT;

	// set some variables that are used to acquire data
	int lowChan = 0;
	int highChan = 1;
	AiInputMode inputMode;
	Range range;
	int samplesPerChannel = 64; // 100 would be desired, but for Fs=5kHz daq feeds in 64-samples batches - if it is not 2^x then works bad somehow!!!
	// cont.: 32 seems minimum for desired frequencies (test with printing totalcount only) e.g. for 10kHz, 32 i.e. batch every 3.2ms ~ 300 frames/second with 10kHz
	// for 100 kHz 32 samples in a package works with 1 channel
	int samplesPerChannel_checked = 64;
	double rate = 40000;//40k/38400;
	ScanOption scanOptions = (ScanOption) (SO_DEFAULTIO | SO_CONTINUOUS);
	AInScanFlag flags = AINSCAN_FF_DEFAULT;

	int hasAI = 0;
	int hasPacer = 0;
	int numberOfChannels = 0;
	int index = 0;

	char inputModeStr[MAX_STR_LENGTH];
	char rangeStr[MAX_STR_LENGTH];
	char scanOptionsStr[MAX_SCAN_OPTIONS_LENGTH];

	int chanCount = 0;
	double* buffer = NULL;
	UlError err = ERR_NO_ERROR;

	int i = 0;
	int __attribute__((unused)) ret;
	char c;

	// Get descriptors for all of the available DAQ devices
	err = ulGetDaqDeviceInventory(interfaceType, devDescriptors, &numDevs);

	if (err != ERR_NO_ERROR)
		goto end;

	// verify at least one DAQ device is detected
	if (numDevs == 0)
	{
		printf("No DAQ device is detected\n");
		goto end;
	}

	printf("Found %d DAQ device(s)\n", numDevs);
	for (i = 0; i < (int) numDevs; i++)
		printf("  [%d] %s: (%s)\n", i, devDescriptors[i].productName, devDescriptors[i].uniqueId);

	if(numDevs > 1)
		descriptorIndex = selectDAQDevice(numDevs);

	// get a handle to the DAQ device associated with the first descriptor
	daqDeviceHandle = ulCreateDaqDevice(devDescriptors[descriptorIndex]);

	if (daqDeviceHandle == 0)
	{
		printf ("\nUnable to create a handle to the specified DAQ device\n");
		goto end;
	}

	// verify the specified device supports analog input
	err = getDevInfoHasAi(daqDeviceHandle, &hasAI);
	if (!hasAI)
	{
		printf("\nThe specified DAQ device does not support analog input\n");
		goto end;
	}

	// verify the specified device supports hardware pacing for analog input
	err = getAiInfoHasPacer(daqDeviceHandle, &hasPacer);
	if (!hasPacer)
	{
		printf("\nThe specified DAQ device does not support hardware paced analog input\n");
		goto end;
	}

	printf("\nConnecting to device %s - please wait ...\n", devDescriptors[descriptorIndex].devString);

	// establish a connection to the DAQ device
	err = ulConnectDaqDevice(daqDeviceHandle);

	if (err != ERR_NO_ERROR)
		goto end;

	// get the first supported analog input mode
	err = getAiInfoFirstSupportedInputMode(daqDeviceHandle, &numberOfChannels, &inputMode, inputModeStr);

	if (highChan >= numberOfChannels)
		highChan = numberOfChannels - 1;

	chanCount = highChan - lowChan + 1;

	// allocate a buffer to receive the data
	buffer = (double*) malloc(chanCount * samplesPerChannel * sizeof(double));

	if(buffer == NULL)
	{
		printf("\nOut of memory, unable to create scan buffer\n");
		goto end;
	}

	// get the first supported analog input range
	err = getAiInfoFirstSupportedRange(daqDeviceHandle, inputMode, &range, rangeStr);

	ConvertScanOptionsToString(scanOptions, scanOptionsStr);

	printf("\n%s ready\n", devDescriptors[descriptorIndex].devString);
	printf("    Function demonstrated: ulAInscan()\n");
	printf("    Channels: %d - %d\n", lowChan, highChan);
	printf("    Input mode: %s\n", inputModeStr);
	printf("    Range: %s\n", rangeStr);
	printf("    Samples per channel: %d\n", samplesPerChannel);
	printf("    Rate: %f\n", rate);
	printf("    Scan options: %s\n", scanOptionsStr);
	printf("\nHit ENTER to continue\n");

	ret = scanf("%c", &c);

	// clear the display
	ret = system("clear");

	// MQTT:
	MQTTClient client;
    	MQTTClient_connectOptions conn_opts = MQTTClient_connectOptions_initializer;
    	MQTTClient_message pubmsg = MQTTClient_message_initializer;
    	MQTTClient_deliveryToken token;
    	int rc;
    	MQTTClient_create(&client, ADDRESS, CLIENTID,
        	MQTTCLIENT_PERSISTENCE_NONE, NULL);
    	conn_opts.keepAliveInterval = 20;
    	conn_opts.cleansession = 1;
    	MQTTClient_setCallbacks(client, NULL, connlost, NULL, NULL);
    	if ((rc = MQTTClient_connect(client, &conn_opts)) != MQTTCLIENT_SUCCESS)
    	{
        	printf("Failed to connect, return code %d\n", rc);
        	exit(EXIT_FAILURE);
    	}
    	//pubmsg.payload = PAYLOAD;
    	//pubmsg.payloadlen = strlen(PAYLOAD);
        pubmsg.qos = QOS;
	pubmsg.retained = 0;
	deliveredtoken = 0;

	// integration for 2 channels and 2xmore frequency
        int int_rate = rate;
        int n_channel = highChan - lowChan + 1;
	int signal_channel = 1;

	/// ############################
        // prepare for phasor estimation
	int debug_sliding = 0;
	int debug_tagging = 0;

	int periods_for_calc = 3;
	int f_nom = 50;
        long long int mem_totalcount = 0;
        int M = periods_for_calc*rate/f_nom;
	int ki, k, mi, m, xi, cursor;

	int samples_in_period = rate/f_nom; // 50k/50=1000
	int sliding_in_samples = samplesPerChannel_checked; // when samplesPerChannel, sliding is every new batch of samples
	int phasor_counter = 0;

	int x_buffer_size = n_channel*((periods_for_calc)*samples_in_period + sliding_in_samples); // n_channel*(3*1000+16)=2*3032 - for all samples
	int x_buffer_ch1_size = (periods_for_calc)*samples_in_period + sliding_in_samples;

	double* x_buffer = malloc(x_buffer_size*sizeof(double)); // array[400]
	double* x_buffer_ch1 = malloc(x_buffer_ch1_size*sizeof(double));
	// other signals separately as _s2 etc.

	//double x_buffer[400];
	for(xi=0;xi<x_buffer_size;xi++) { // just some initial non-zero values in memory
	        x_buffer[xi] = xi+0.2;
        }
	printf("Check the size of the Wm1, Wm2!! It should be: %d\n", M);
        double complex Wm1[5][2400], Wm2[5][2400]; //2400/2304
        printf("%f + i*%f\n", creal(cexp(-I)), cimag(cexp(-I)));

        for (ki=0; ki<5; ki++){  // for k 1 to 5 i.e. 3+/-2, not universal, depending on M etc.
                k = ki+1;
                for(m=0; m<M; m++){
                        Wm1[ki][m] = cexp(-2*I*M_PI*k*m/M);
                        Wm2[ki][m] = cexp(2*I*M_PI*k*m/M);
                }
        }
	printf("test1");

	// assigning timetags to samples $$$$$$$$$$$$$$$$$$$$$$
        float dt; // in microseconds
        float itt_compens; // compensation for all tags: trigger 1us, trigger dt/2
        float v_pps = 0.5; // it is about 3.3V
        float itt_last, itt_first, itt_ultimate, itt_pps, itt_last_ch1, ptt_ch1;
        int samplesAcq = samplesPerChannel_checked * n_channel; // 16*2 for 2 channels and 100kS
        float real_rate;

        int pps_detect, pps_xi, first_pps_xi, pps_idx, first_pps_detect, first_pps_detect_diff;
        long long int pps_idx_sum = 0;
	long long int first_pps_shift = 0;
        long long int totalCount_updated, totalCount_updated_from_firstPPS;
        float last_pps_sample;
	first_pps_detect = 0;
	first_pps_detect_diff = 0; // flag for avoiding DIFF missing samples after first pps trigger
	int reset_at_pps = 0;
	int pps_count = 0;
	// $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

	// DFT
	double Xk_row_sum_real, Xk_row_sum_imag;
        double complex Xk[5];
        double complex Xk_H1, Xk_H2, Xk_H3; // hann in dft domain

	// Sliding
	double complex Xk_slide[5]; // for storing new Xk for 5 k of interest
        double complex into0, slide; // for storing temporary stuff during slides
        int n, sui, n_m;
        double re_into0, im_into0;
	float phasor_abs;

	// Interpolation
	double B, alpha, delta_bin, time_win, delta_f, f_estim, A_estim, ph_estim;
        int ki_max, k_max, esign, ki_max_esign;
        double Xk_Habs[3], Xk_Hmax;
	B = 1199.5; //1199.5/1151.5 // sum of Hann windowing for 2400 i.e. 3*40k/50 (check in matlab: sum(hann(2400)) )

	// #############################################################
	// start the acquisition #######################################
	err = ulAInScan(daqDeviceHandle, lowChan, highChan, inputMode, range, samplesPerChannel, &rate, scanOptions, flags, buffer);

	struct timeval timer_usec;
	long long int timestamp_usec_start, timestamp_usec_stop;
        long long int timestamp_usec; /* timestamp in microsecond */
        long long int timestamp_usec_prev;
	long long int total_count;

	// calculated dt and compensation based on the real rate ~500035.74 for 2ch/50kS
        real_rate = 40000;//rate;
        dt = 1/(real_rate*n_channel) * 1000000;
        itt_compens = 0.0; // no trigger so no predictable compensation (?)

	printf("real_rate: %f dt: %f", real_rate, dt);
//	exit(0);
	if(err == ERR_NO_ERROR)
	{
		ScanStatus status;
		TransferStatus transferStatus;

		// get the initial status of the acquisition
		ulAInScanStatus(daqDeviceHandle, &status, &transferStatus);

		/* starting timestamp. */
                if (!gettimeofday(&timer_usec, NULL)) {
                timestamp_usec_start = ((long long int) timer_usec.tv_sec) * 1000000ll + (long long int) timer_usec.tv_usec;
                }
                else {
                timestamp_usec_start = -1;
                }
                printf("\nstart tag: %lld", timestamp_usec_start);

		while(status == SS_RUNNING && err == ERR_NO_ERROR && !enter_press() )
		{
			// get the current status of the acquisition
			err = ulAInScanStatus(daqDeviceHandle, &status, &transferStatus);

			// reset the cursor to the top of the display and
			// show the termination message
			//resetCursor();
			//printf("Hit 'Enter' to terminate the process\n\n");
			//printf("Active DAQ device: %s (%s)\n\n", devDescriptors[descriptorIndex].productName, devDescriptors[descriptorIndex].uniqueId);
			//printf("actual scan rate = %f\n\n", rate);

			index = transferStatus.currentIndex;
			/*if (transferStatus.currentTotalCount < 4000000 && transferStatus.currentTotalCount > 0 ) {
				printf("\n\ncurrentScanCount =  %-10llu", transferStatus.currentScanCount);
				printf("\ncurrentTotalCount = %-10llu", transferStatus.currentTotalCount);
				//printf("\ncurrentIndex =      %-10d \n", index);
			}
			printf("\nPre-Current daq buffer:\n");
                        for(xi=0;xi<n_channel*samplesPerChannel_checked;xi++) {
	                       printf("%+-10.6f ", buffer[xi]);
        	                if((xi+1)%10==0){printf("\n");}
                	}*/

			// DETECT FIRST PPS
			for(xi=0; xi<samplesAcq; xi++) {
				//printf("%+-10.6f ", buffer[xi]);
				if (first_pps_detect==0 && xi==0 && last_pps_sample<v_pps && buffer[0]>=v_pps) {
 	                               	printf("<-PPS!\t");
                                       	first_pps_detect = 1;
					//pps_count++;
                             		//pps_detect = 1;
					first_pps_xi = xi;
					//pps_xi = xi;
					pps_idx_sum = transferStatus.currentTotalCount - (samplesAcq-first_pps_xi);
                                	first_pps_shift = transferStatus.currentTotalCount - (samplesAcq-first_pps_xi);
				}
                                if (first_pps_detect==0 && xi>=n_channel && xi%n_channel==0  && buffer[xi-n_channel]<v_pps && buffer[xi]>=v_pps) {
                                       	printf("<-PPS!\t");
                                       	first_pps_detect = 1;
					//pps_count++;
					//pps_detect = 1;
                                       	first_pps_xi = xi;
					//pps_xi = xi;
					pps_idx_sum = transferStatus.currentTotalCount - (samplesAcq-first_pps_xi);
					first_pps_shift = transferStatus.currentTotalCount - (samplesAcq-first_pps_xi);
				}
				//if ((xi+1)%10==0){printf("\n");}
			}
			//if (first_pps_detect == 1) {exit(0);}
			//continue;

			if(index >= 0)
			{
				// display the data
				if(debug_sliding==1){
				for (i = 0; i < chanCount; i++)
				{
                                        //printf("nowy sample%+-10.6f\n", buffer[index]);
                                        //printf("b0: %+-10.6f\n", buffer[0]);
					//printf("b1: %+-10.6f\n", buffer[1]);
                                        //printf("b2: %+-10.6f\n", buffer[2]);
                                        //printf("b3: %+-10.6f\n", buffer[3]);
                                        //printf("b4: %+-10.6f\n", buffer[4]);
					//printf("%+-10.6f\n", buffer[index+3]);
                                        printf("size buffer: %d\n", (int)(sizeof(buffer))); // it's passed i.e. it's not the real size of the buffer
                                        printf("size buffer element: %d\n", (int)(sizeof(buffer[0])));
				        printf("size: %d\n", (int)(sizeof(buffer) / sizeof(buffer[0])));
                                }}

				if(first_pps_detect == 1 && transferStatus.currentTotalCount != mem_totalcount && (int)transferStatus.currentScanCount % (int)samplesPerChannel_checked == 0){
					if ((int)transferStatus.currentTotalCount-(int)mem_totalcount != n_channel*samplesPerChannel_checked) {
						if (first_pps_detect_diff == 0) {
							printf("\nDIFF due to the first pps trigger.");
							first_pps_detect_diff = 1;
							mem_totalcount = transferStatus.currentTotalCount;
						} else {
							printf("\nDIFF: %d", (int)transferStatus.currentTotalCount-(int)mem_totalcount);
							printf("\n\nMISSED SAMPLES? Check Scan/Total Counts vs. samplesPerChannel_checked, samplesPerChannel_checked, sliding_in_samples !!! \n");
							exit(0);
						}
					}
					/*if ((float)((int)transferStatus.currentTotalCount-(int)mem_totalcount)/(float)n_channel != samplesPerChannel_checked) {
						if (first_pps_detect_diff == 0) {
                                                        printf("\nDIFF due to the first pps trigger.");
                                                        first_pps_detect_diff = 1;
							mem_totalcount = transferStatus.currentTotalCount;
						} else {
							printf("\nDIFF: %d", (int)transferStatus.currentTotalCount-(int)mem_totalcount);
							printf("\n\nWRONG sampled per channel settings? Check Scan/Total Counts vs. samplesPerChannel_checked !!! \n");
							exit(0);
						}
					}*/
					mem_totalcount = transferStatus.currentTotalCount; // in order to avoid repeating for the same counter

					if(debug_sliding==1){
						printf("\n\ncurrentScanCount =  %-10llu", transferStatus.currentScanCount);
	   	                                printf("\ncurrentTotalCount = %-10llu", transferStatus.currentTotalCount);
						printf("\nCurrent daq buffer:\n");
						for(xi=0;xi<n_channel*samplesPerChannel_checked;xi++) {
							printf("%+-10.6f ", buffer[xi]);
						if((xi+1)%10==0){printf("\n");}
						}
					}
					//continue;
					// move samples to make space for new samples - the "space" has the value of the old samples
					for(xi=0; xi<x_buffer_size-n_channel*samplesPerChannel_checked; xi++) {
                                                x_buffer[xi] = x_buffer[xi+n_channel*samplesPerChannel_checked];
                                        }

					// copy new ones to the end
					memcpy(x_buffer+(x_buffer_size-n_channel*samplesPerChannel_checked), buffer, n_channel*samplesPerChannel_checked*sizeof(double));

					// current buffer for calculations
					if(debug_sliding==1){
						printf("\n READY x buffer: %d samples or n_channels*(3*int_rate/50+sliding_in_samples):\n", n_channel*(periods_for_calc*int_rate/f_nom+sliding_in_samples));
						for(xi=0;xi<x_buffer_size;xi++) {
							printf("%+-10.6f ", x_buffer[xi]);
							if((xi+1)%10==0){printf("\n");}
						}
					}

					// ASSIGNING ITT to the samples etc.
					totalCount_updated_from_firstPPS = transferStatus.currentTotalCount - first_pps_shift; // for counting total
					totalCount_updated = transferStatus.currentTotalCount - pps_idx_sum; // for counting every second and timetags ptt
					if(debug_tagging==1){
						printf("\ntotalCount_updated: %d dt: %f real_rate: %f itt_compens: %f", (int)totalCount_updated, dt, real_rate, itt_compens);
					}
                                	itt_last =  (totalCount_updated-1)*dt + itt_compens;
                                	itt_first = round((totalCount_updated-1)*dt + itt_compens - (samplesAcq-1)*dt);

                               		if(debug_tagging==1) {printf("\nCurrent daq buffer: (compensated itt from %dr to %f)\n", (int)itt_first, itt_last);}
                               		for(xi=0; xi<samplesAcq; xi++) {
                                       		if(debug_tagging==1) {printf("%+-10.6f ", buffer[xi]);}
		                        	// detection of PPS
						if (xi==0 && last_pps_sample<v_pps && buffer[0]>=v_pps) {
        	                                       	if(debug_tagging==1) {printf("<-PPS(last)!\t"); printf("last pps s: %f", last_pps_sample);}
							pps_detect = 1;
							pps_xi = xi;
                       	                		pps_count++;
						}
                               	        	if (xi>=n_channel && xi%n_channel==0  && buffer[xi-n_channel]<v_pps && buffer[xi]>=v_pps) {
                                       	        	if(debug_tagging==1) {printf("<-PPS!\t");}
							pps_detect = 1;
							pps_xi = xi;
							pps_count++;
                                       		}
                                       		if (debug_tagging==1 && (xi+1)%10==0){printf("\n");}
	                               	}
					last_pps_sample = buffer[samplesAcq-n_channel];

					if (pps_detect == 1) {
                                        	//printf("\npps_xi=%d", pps_xi);
                                        	pps_idx = totalCount_updated-(samplesAcq-pps_xi-1);
                                        	//printf("\npps_idx: %d", pps_idx);
                                        	pps_idx_sum = pps_idx_sum + pps_idx - 1;

                                        	itt_ultimate = round((pps_idx-1)*dt + itt_compens);
                                        	//printf("\nitt_ultimate: %d", (int)itt_ultimate);
                                    	    	itt_pps = round(itt_compens);
                                    	    	itt_last = (samplesAcq-pps_xi-1)*dt + itt_compens;

                                        	printf("\nUpdated itt. itt_first(same): %d , itt_ultimate: %d , itt_pps: %d, itt_last: %f", (int)itt_first, (int)itt_ultimate, (int)itt_pps, itt_last );
					}
					pps_detect = 0;

					if(totalCount_updated_from_firstPPS < periods_for_calc*samples_in_period*n_channel) { // 3*rate/50
						//printf("\nWaiting. Channel count so far: %-10llu. After first PPS: %d. PPS count: %d. Updated Total count: %d",
							//transferStatus.currentScanCount, (int)totalCount_updated_from_firstPPS, pps_count, (int)totalCount_updated);
						printf("\npps_idx_sum: %lld first_pps_shift: %lld", pps_idx_sum, first_pps_shift);
						continue;
					}
					else if((totalCount_updated_from_firstPPS >= periods_for_calc*samples_in_period*n_channel) && (totalCount_updated_from_firstPPS < (periods_for_calc*samples_in_period+sliding_in_samples)*n_channel)) {
						// minimum amount for first DFT
						//printf("\nFirst DFT. Channel count so far: %-10llu. After first PPS: %d. PPS count: %d. Updated Total count: %d",
							//transferStatus.currentScanCount, (int)totalCount_updated_from_firstPPS, pps_count, (int)totalCount_updated);

						if(debug_sliding==1){
							printf("\nNewest %d samples: n_channel*(3*int_rate/50). Oldest excess is cut off.\n", n_channel*periods_for_calc*int_rate/f_nom);
							for(xi = n_channel*sliding_in_samples; xi < x_buffer_size; xi++) {
	        	                                        printf("%+-10.6f ", x_buffer[xi]);
        	        	                                if((xi+1)%10==0){printf("\n");}
                	        	                }
						}

						// copy, create new array with only one signal memcpy(x_buffer+(x_buffer_size-n_channel*samplesPerChannel_checked), buffer, n_channel*samplesPerChannel_checked*sizeof(double));
						printf("\nCh1 samples (whole buffer, also old ones) copied to x_buffer_ch1:\n");
						for(xi=0; xi<x_buffer_ch1_size; xi++) {
							//printf("\nxi= %d , to be copied: %+-10.6f ", xi, x_buffer[1+n_channel*xi]);
							memcpy(x_buffer_ch1+xi, x_buffer+n_channel*xi+1, sizeof(double));
						}
						if(debug_sliding==1){
							for(xi=0; xi<x_buffer_ch1_size; xi++) {
                                                        	printf("%+-10.6f ", x_buffer_ch1[xi]);
                                                        	if((xi+1)%10==0){printf("\n");}
                                                	}
							printf("\nCh1 samples (with old discarded):\n");
							for(xi=sliding_in_samples; xi<x_buffer_ch1_size; xi++) {
                                                        	printf("%+-10.6f ", x_buffer_ch1[xi]);
                                                        	if((xi+1)%10==0){printf("\n");}
                                                	}
						}
						// derive phasor timetab from itt_last from whole buffer (all channels)
						itt_last_ch1 = itt_last;
						ptt_ch1 = round(itt_last_ch1 - (M/2-0.5) * dt * n_channel); // average of two middle tags from ch1
						if(debug_tagging==1){
							printf("\n itt_last (all-channels-buffer): %f", itt_last);
							printf("\n itt_last (ch1-channel-buffer): %f", itt_last_ch1);
							printf("\n ptt_ch1 (round): %f", ptt_ch1);
						}
						// continue as before but with x_buffer_ch1

					        for (ki=0; ki<5; ki++){
							Xk_row_sum_real = 0;
							Xk_row_sum_imag = 0;

               						for(mi=0; mi<M; mi++){
           							Xk_row_sum_real = Xk_row_sum_real + x_buffer_ch1[sliding_in_samples + mi]*creal(Wm1[ki][mi]); //sliding_in_samples  in order to remove excess from buffer (which is for sliding later)
								Xk_row_sum_imag = Xk_row_sum_imag + x_buffer_ch1[sliding_in_samples + mi]*cimag(Wm1[ki][mi]);
               						}
							Xk[ki] = Xk_row_sum_real + Xk_row_sum_imag * I;
							//printf("\nXk%d = %.1f% + .1fi", ki+1, creal(Xk[ki]), cimag(Xk[ki])); // prints each Xk of interest
      						}
						phasor_counter++;
						//printf("\n End of Xk calculations for all k of interest. Xk3=%f %fj\nPhasor counter=%d\n", creal(Xk[2]), cimag(Xk[2]), phasor_counter);

					}
					else if(totalCount_updated_from_firstPPS >= (periods_for_calc*samples_in_period+sliding_in_samples)*n_channel) {
						// sufficient amount for sliding i.e. minimum + one period
						// copy, create new array with only one signal memcpy(x_buffer+(x_buffer_size-n_channel*samplesPerChannel_checked), buffer, n_channel*samplesPerChannel_checked*sizeof(double));
                                                for(xi=0; xi<x_buffer_ch1_size; xi++) {
                                                        //printf("\nxi= %d , to be copied: %+-10.6f ", xi, x_buffer[1+n_channel*xi]);
                                                        memcpy(x_buffer_ch1+xi, x_buffer+n_channel*xi+1, sizeof(double));
                                                }
						if(debug_sliding==1) {
							printf("\nCh1 samples (whole buffer, also old ones) copied to x_buffer_ch1:\n");
                	                                for(xi=0; xi<x_buffer_ch1_size; xi++) {
                        	                                printf("%+-10.6f ", x_buffer_ch1[xi]);
                                	                        if((xi+1)%10==0){printf("\n");}
                                        	        }
                 	                               printf("\nCh1 samples (with old discarded):\n");
                        	                        for(xi=sliding_in_samples; xi<x_buffer_ch1_size; xi++) {
                                	                        printf("%+-10.6f ", x_buffer_ch1[xi]);
                                        	                if((xi+1)%10==0){printf("\n");}
                                                	}
						}

						// derive phasor timetag from itt_last from whole buffer (all channels)
	                        	        itt_last_ch1 = itt_last;
                	                        ptt_ch1 = round(itt_last_ch1 - (M/2-0.5) * dt * n_channel); // average of two middle tags from ch1
						if (reset_at_pps==0) {
							ptt_ch1 = ptt_ch1 + (pps_count-1)*1000000;
						}
						if(debug_tagging==1){
							printf("\n itt_last (all-channels-buffer): %f", itt_last);
							printf("\n itt_last (ch1-channel-buffer): %f", itt_last_ch1);
                        	                	printf("\n ptt_ch1: %f", ptt_ch1);
                                                }
						// continue as before but with x_buffer_ch1

						if(debug_sliding==1/*transferStatus.currentScanCount>10000*/){
                                                	printf("\nNewest %d samples of ch1: (3*int_rate/50). Oldest excess is cut off.\n", periods_for_calc*int_rate/f_nom);
                                                	for(xi = sliding_in_samples; xi < x_buffer_ch1_size; xi++) {
                                                        	printf("%+-10.6f ", x_buffer_ch1[xi]);
                                                        	if((xi+1)%10==0){printf("\n");}
                                               		}
						}

						for(ki=0;ki<5;ki++) {
							k=ki+1;
							into0 = Xk[ki];
							slide = 0;
							for(sui=0; sui<sliding_in_samples; sui++){
								n = M+sui;
								n_m = n%M;

								re_into0 = creal(into0) + creal(Wm1[ki][n_m]) * (-x_buffer_ch1[n-M] + x_buffer_ch1[n]);
								im_into0 = cimag(into0) + cimag(Wm1[ki][n_m]) * (-x_buffer_ch1[n-M] + x_buffer_ch1[n]);

								into0 = re_into0 + im_into0*I;
								slide = Wm2[ki][n_m+1] * into0;

								if(debug_sliding==1)
								{
									printf("\nki=%d n=%d n_m=%d Wm1[%d,%d]= %f + i*%f",ki,n,n_m,ki,n_m, creal(Wm1[ki][n_m]), cimag(Wm1[ki][n_m]));
									printf("\ninto0:%f %fi\tWm1:%f %fi\t samples sum:%f", creal(into0), cimag(into0), creal(Wm1[ki][n_m]), cimag(Wm1[ki][n_m]), (-x_buffer_ch1[n-M] + x_buffer_ch1[n]) );
									printf("\n%d %d remove(%d): %f add(%d):%f ... into0_abs:%f  into0_deg:%f  slide_abs:%f  slide_deg:%f",n_m, n, n-M, x_buffer_ch1[n-M], n, x_buffer_ch1[n],
										cabs(into0),
										carg(into0)*180/M_PI,
										cabs(slide),
                                                                        	carg(slide)*180/M_PI);
								}
							}
							Xk[ki] = slide;
						}
						if(debug_sliding==1)
						{
							printf("\nXk_slide: after sliding:");
                                                	for (ki=0;ki<5;ki++) {
                                                        	printf("\nXk%d = %.1f% + .1fi", ki+1, creal(Xk[ki]), cimag(Xk[ki]));
                                                	}
						}
						phasor_counter++;
					}

                                Xk_H1 = -0.25*creal(Xk[0]) + 0.5*creal(Xk[1]) - 0.25*creal(Xk[2]) + (-0.25*cimag(Xk[0]) + 0.5*cimag(Xk[1]) - 0.25*cimag(Xk[2]))*I;
                                Xk_H2 = -0.25*creal(Xk[1]) + 0.5*creal(Xk[2]) - 0.25*creal(Xk[3]) + (-0.25*cimag(Xk[1]) + 0.5*cimag(Xk[2]) - 0.25*cimag(Xk[3]))*I;
                                Xk_H3 = -0.25*creal(Xk[2]) + 0.5*creal(Xk[3]) - 0.25*creal(Xk[4]) + (-0.25*cimag(Xk[2]) + 0.5*cimag(Xk[3]) - 0.25*cimag(Xk[4]))*I;

                                time_win = (float)periods_for_calc / (float)f_nom;

                                Xk_Habs[0] = cabsf(Xk_H1)/B;
                                Xk_Habs[1] = cabsf(Xk_H2)/B;
                                Xk_Habs[2] = cabsf(Xk_H3)/B;

				Xk_Hmax = Xk_Habs[1];
                                ki_max = 1;
                                k_max = ki_max + 2;
                                if (Xk_Habs[0]>Xk_Habs[2]) {
  	                        	esign = -1;
                                } else {
                                        esign = 1;
                                }
                                ki_max_esign = ki_max + esign;
                                alpha = Xk_Habs[ki_max] / Xk_Habs[ki_max_esign];
                                delta_bin = esign * (2-alpha) / (1+alpha);
                                delta_f = 1/time_win;

                                f_estim = (k_max + delta_bin) * delta_f;
                                A_estim = 2*Xk_Hmax*(M_PI*delta_bin*(1-delta_bin*delta_bin))/sin(M_PI*delta_bin);
                                ph_estim = cargf(Xk_H2)-M_PI*delta_bin;

				printf("\nSending f,A,ph,ptt as MQTT message.\n");
                                send_phasor(f_estim, A_estim, ph_estim*180/M_PI, ptt_ch1, client, pubmsg, token);

				}
			//last_pps_sample = buffer[samplesAcq-n_channel];
			//usleep(1000);
			}
		total_count = transferStatus.currentTotalCount;
		}

		// stop the acquisition if it is still running
		if (status == SS_RUNNING && err == ERR_NO_ERROR)
		{
			err = ulAInScanStop(daqDeviceHandle);
		}
	}

	// disconnect from the DAQ device
	ulDisconnectDaqDevice(daqDeviceHandle);
	// end from MQTT
	MQTTClient_disconnect(client, 10000);
	MQTTClient_destroy(&client);
end:
	/* end timestamp. */
        if (!gettimeofday(&timer_usec, NULL)) {
        timestamp_usec_stop = ((long long int) timer_usec.tv_sec) * 1000000ll + (long long int) timer_usec.tv_usec;
        }
        else {
        timestamp_usec_stop = -1;
        }
	printf("\nstop-start timetags: %lld Total count: %-10llu", timestamp_usec_stop-timestamp_usec_start, total_count);
	printf("\nRate: %0.f Hz. One period (50Hz) has ~%0.f samples. Acquisition (should be) every: %d samples i.e. every %0.fus", rate, rate/50, samplesPerChannel_checked, samplesPerChannel_checked/rate*1000000);
	printf("\nResulting average acquisition: total time/total acquisitions = %0.2f\n", (float)(timestamp_usec_stop-timestamp_usec_start)/((float)total_count/(float)samplesPerChannel_checked));

	// release the handle to the DAQ device
	if(daqDeviceHandle)
		ulReleaseDaqDevice(daqDeviceHandle);

	// release the scan buffer
	if(buffer)
		free(buffer);

	if(err != ERR_NO_ERROR)
	{
		char errMsg[ERR_MSG_LEN];
		ulGetErrMsg(err, errMsg);
		printf("Error Code: %d \n", err);
		printf("Error Message: %s \n", errMsg);
	}

	return 0;
}
