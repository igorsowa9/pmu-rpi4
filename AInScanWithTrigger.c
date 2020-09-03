/*
    Wrapper call demonstrated:        ulAInSetTrigger()

    Purpose:                          Setup an external trigger

    Demonstration:                    Uses the first available trigger type to
                                      set up an external trigger that is used
                                      to start a scan

    Steps:
    1. Call ulGetDaqDeviceInventory() to get the list of available DAQ devices
    2. Call ulCreateDaqDevice() to to get a handle for the first DAQ device
    3. Verify the DAQ device has an analog input subsystem
    3. Verify the analog input subsystem has a hardware pacer
    4. Get the supported trigger types
    5. Get the supported queue types
    6. Fill the channel array
    7. Call ulConnectDaqDevice() to establish a UL connection to the DAQ device
    8. Call ulAInSetTrigger to set the external trigger
    9. Call ulAInScan() to start the scan of A/D input channels
    10. Call ulAiScanStatus to check the status of the background operation
    11. Display the data for each channel
    12. Call ulAinScanStop() to stop the background operation
    13. Call ulDisconnectDaqDevice and ulReleaseDaqDevice() before exiting the process
*/

#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <time.h>
#include <sys/timeb.h>
#include <sys/time.h>
#include <math.h>
#include "uldaq.h"
#include "utility.h"
#include "math.h"

#define MAX_DEV_COUNT  100
#define MAX_STR_LENGTH 64
#define MAX_SCAN_OPTIONS_LENGTH 256

#define M_PI 3.14159265358979323846

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
	int samplesPerChannel = 16; //64
	double rate = 50000; //5-100k
	ScanOption scanOptions = (ScanOption) (SO_DEFAULTIO | SO_CONTINUOUS | SO_EXTTRIGGER);
	//ScanOption scanOptions = (ScanOption) (SO_DEFAULTIO | SO_CONTINUOUS | SO_EXTTRIGGER);
	AInScanFlag flags = AINSCAN_FF_DEFAULT;

	int hasAI = 0;
	int hasPacer = 0;
	int numberOfChannels = 0;
	TriggerType triggerType;
	int index = 0;

	char inputModeStr[MAX_STR_LENGTH];
	char rangeStr[MAX_STR_LENGTH];
	char triggerTypeStr[MAX_STR_LENGTH];
	char scanOptionsStr[MAX_SCAN_OPTIONS_LENGTH];

	// allocate a buffer to receive the data
	int chanCount = 0;
	double *buffer = NULL;
	UlError err = ERR_NO_ERROR;

	int __attribute__ ((unused))ret;
	char c;
	int i = 0;

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

	// establish a connection to the device
	err = ulConnectDaqDevice(daqDeviceHandle);

	if (err != ERR_NO_ERROR)
	{
		goto end;
	}

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

	// get the first supported trigger type (this returns a digital trigger type)
	getAiInfoFirstTriggerType(daqDeviceHandle, &triggerType, triggerTypeStr);

	ConvertScanOptionsToString(scanOptions, scanOptionsStr);

	printf("\n%s ready\n", devDescriptors[descriptorIndex].devString);
	printf("    Function demonstrated: ulAInSetTrigger()\n");
	printf("    Channels: %d - %d\n", lowChan, highChan);
	printf("    Input mode: %s\n", inputModeStr);
	printf("    Range: %s\n", rangeStr);
	printf("    Samples per channel: %d\n", samplesPerChannel);
	printf("    Rate: %f\n", rate);
	printf("    Scan options: %s\n", scanOptionsStr);
	printf("    Trigger type: %s\n", triggerTypeStr);
	printf("\nHit ENTER to continue\n");

	ret = scanf("%c", &c);

	// clear the display
	ret = system("clear");

	//	set the trigger
	//
	// this example uses the default values for setting the trigger so there is no need to call this function ...
	// if you want to change the trigger type (or any other trigger parameter), uncomment this function call and
	// change the trigger type (or any other parameter)
	//err = ulAInSetTrigger( daqDeviceHandle, triggerType, 0, 0.0, 0.0, 0);

	// variables to measure time
	struct timeval timer_usec;
        long long int timestamp_usec_start, timestamp_usec_stop;
        long long int timestamp_usec; /* timestamp in microsecond */
        long long int timestamp_usec_prev;
        long long int total_count;
	int is_first_timetag = 0;

	// other variables
	int mem_totalcount = 0;
	int retrigger = 0;
	int n_channel = highChan - lowChan + 1;

	// assigning timetags to samples
	float dt; // in microseconds
	float itt_compens; // compensation for all tags: trigger 1us, trigger dt/2
	float v_pps = 1.5; // it is about 3.3V
	float itt_last, itt_first, itt_ultimate, itt_pps;
	int samplesAcq = samplesPerChannel * n_channel; // 16*2 for 2 channels and 100kS
	float real_rate;

	int pps_detect, pps_xi, pps_idx;
	long long int pps_idx_sum = 0;
	long long int totalCount_updated;
	float last_pps_sample;

	int xi;
	int f_nom = 50;
	int periods_for_calc = 3;
	int samples_in_period = rate/f_nom;
	int sliding_in_samples = samplesPerChannel;

	while(1==1) // for retriggering? (when by PPS)
	{
	// start the acquisition
	err = ulAInScan(daqDeviceHandle, lowChan, highChan, inputMode, range, samplesPerChannel, &rate, scanOptions, flags, buffer);

	// calculated dt and compensation based on the real rate ~500035.74 for 2ch/50kS
	real_rate = rate;
	dt = 1/(rate*n_channel) * 1000000;
	itt_compens = 1 + dt/2;

	if(err == ERR_NO_ERROR)
	{
		ScanStatus status;
		TransferStatus transferStatus;

		// get the initial status of the acquisition
		ulAInScanStatus(daqDeviceHandle, &status, &transferStatus);

		printf ("Hit 'Enter' to quit waiting for trigger\n\n");
		printf("Active DAQ device: %s (%s)\n\n", devDescriptors[descriptorIndex].productName, devDescriptors[descriptorIndex].uniqueId);
		printf ("Waiting for trigger ...\n");

		while(status == SS_RUNNING && err == ERR_NO_ERROR && !enter_press())
		{
			if(transferStatus.currentTotalCount > 4.004*n_channel*rate)
			{
			break;
			}
			if(is_first_timetag==0)
			{
				/* starting timestamp. */
		                if (!gettimeofday(&timer_usec, NULL)) {
               			timestamp_usec_start = ((long long int) timer_usec.tv_sec) * 1000000ll + (long long int) timer_usec.tv_usec;
                		}
                		else {
                		timestamp_usec_start = -1;
                		}
                		printf("\nstart tag: %lld", timestamp_usec_start);
				is_first_timetag = 1;
			}

			// get the current status of the acquisition
			err = ulAInScanStatus(daqDeviceHandle, &status, &transferStatus);

			printf("currentScanCount =  %-10llu buffer[0]: %+-10.6f \n", transferStatus.currentScanCount, buffer[0]);
                        //printf("currentTotalCount = %-10llu \n", transferStatus.currentTotalCount);

			index = transferStatus.currentIndex;
			if(err == ERR_NO_ERROR && index >= 0 &&
				transferStatus.currentTotalCount != mem_totalcount && transferStatus.currentTotalCount % samplesPerChannel == 0) // to avoid repetitions and to print only every cycle of acquisition (for all CHs)
			{
				mem_totalcount = transferStatus.currentTotalCount;
				//resetCursor();
				//printf("%-40s\n\n","Hit 'Enter' to terminate the process");
				//printf("Active DAQ device: %s (%s)\n\n", devDescriptors[descriptorIndex].productName, devDescriptors[descriptorIndex].uniqueId);
				//printf("actual scan rate = %f\n\n", rate);

				//printf("no-reps currentScanCount =  %-10llu \n", transferStatus.currentScanCount);
				//printf("currentTotalCount = %-10llu \n", transferStatus.currentTotalCount);
				//printf("currentIndex =      %-10d \n\n", index);

				/* Example of timestamp in microsecond. */
                                /*if (!gettimeofday(&timer_usec, NULL)) {
                                timestamp_usec = ((long long int) timer_usec.tv_sec) * 1000000ll + (long long int) timer_usec.tv_usec;
                                }
                                else {
                                timestamp_usec = -1;
                                }

				printf("\nno-reps currentTotalCount = %llu", transferStatus.currentTotalCount);
                                timestamp_usec_prev = timestamp_usec;
				*/
				totalCount_updated = transferStatus.currentTotalCount - pps_idx_sum;
				itt_last =  round((totalCount_updated-1)*dt + itt_compens);
				itt_first = round((totalCount_updated-1)*dt + itt_compens - (samplesAcq-1)*dt);

				printf("\nCurrent daq buffer: (compensated itt from %dr to %dr)\n", (int)itt_first, (int)itt_last);
				for(xi=0; xi<samplesAcq; xi++) {
                                        printf("%+-10.6f ", buffer[xi]);
					if (xi==0 && last_pps_sample<v_pps && buffer[0]>=v_pps) {
						printf("<-PPS!\t"); pps_detect = 1; pps_xi = xi;
                                        }
					if (xi>=n_channel && xi%n_channel==0  && buffer[xi-n_channel]<v_pps && buffer[xi]>=v_pps) {
						printf("<-PPS!\t"); pps_detect = 1; pps_xi = xi;
					}
                                        if ((xi+1)%10==0){printf("\n");}
				}
				if (pps_detect == 1) {
					printf("\npps_xi=%d", pps_xi);
					pps_idx = totalCount_updated-(samplesAcq-pps_xi-1);
					printf("\npps_idx: %d", pps_idx);
					pps_idx_sum = pps_idx_sum + pps_idx - 1;

					itt_ultimate = round((pps_idx-1)*dt + itt_compens);
					printf("\nitt_ultimate: %d", (int)itt_ultimate);
					itt_pps = round(itt_compens);
					itt_last = round((samplesAcq-pps_xi-1)*dt + itt_compens);

					printf("\nUpdated itt. itt_first(same): %d , itt_ultimate: %d , itt_pps: %d, itt_last: %d", (int)itt_first, (int)itt_ultimate, (int)itt_pps, (int)itt_last );
				}
				//last_pps_sample = buffer[samplesAcq-n_channel]; // saving last sample from channel1 for comparison for PPS
				//if () {
				//	printf("\nPPS! ")
				//}

				printf("\n-------------------------------\n");
				pps_detect = 0;

				continue;
				// ALGORITM FOR SLIDING STARTS HERE:

				if(transferStatus.currentScanCount < periods_for_calc*samples_in_period) { // 3*rate/50
                                        printf("\n\nWaiting. Channel 2 count so far: %-10llu", transferStatus.currentScanCount);
                                }
                                else if((transferStatus.currentScanCount >= periods_for_calc*samples_in_period) && (transferStatus.currentScanCount < periods_for_calc*samples_in_period+sliding_in_samples)) {
                                        // minimum amount for first DFT
                                        printf("\n\nFirst DFT. Channel 2 count so far: %-10llu", transferStatus.currentScanCount);
				}
                                else if(transferStatus.currentScanCount >= periods_for_calc*samples_in_period+sliding_in_samples) {
                                        // sufficient amount for sliding i.e. minimum + one period
                                        printf("\n\nSliding possible. Channel 2 count so far: %-10llu", transferStatus.currentScanCount);
					exit(0);
				}



				//usleep(1000);
			}
			last_pps_sample = buffer[samplesAcq-n_channel]; // saving last sample from channel1 for comparison for PPS
			usleep(1);
		total_count = transferStatus.currentTotalCount;
		}

		if (index < 0)
			printf("Trigger cancelled by user\n");

		// stop the acquisition if it is still running
		if (status == SS_RUNNING && err == ERR_NO_ERROR)
		{
			err = ulAInScanStop(daqDeviceHandle);
		}
	}
	}
	// disconnect from the DAQ device
	ulDisconnectDaqDevice(daqDeviceHandle);

end:

	/* end timestamp. */
        if (!gettimeofday(&timer_usec, NULL)) {
        timestamp_usec_stop = ((long long int) timer_usec.tv_sec) * 1000000ll + (long long int) timer_usec.tv_usec;
        }
        else {
        timestamp_usec_stop = -1;
        }
        printf("\nstop-start timetags: %lld Total count: %-10llu", timestamp_usec_stop-timestamp_usec_start, total_count);
        printf("\nRate: %0.f Hz. One period (50Hz) has ~%0.f samples. Acquisition (should be) every: %d samples i.e. every %0.fus", rate, rate/50, samplesPerChannel, samplesPerChannel/rate*1000000);
        printf("\nResulting average acquisition: total time/total acquisitions = %0.2f\n", (float)(timestamp_usec_stop-timestamp_usec_start)/((float)total_count/(float)samplesPerChannel));

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

