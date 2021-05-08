/* PMU calculation and MQTT publication configurations:*/

/* Connection settings */
#define LOWCH   4 // vPPS from RTDS GTAO
#define HIGHCH  5

/* algorithm/publishing */
#define REPRATE     200 // reporting rate

/* MQTT */
#define ADDRESS     "tcp://localhost:1883" /* remove/comment conn_opts.username, conn_opts.password */
//#define ADDRESS     "ssl://platone.eng.it:8883"
//#define ADDRESS     "platone.eng.it:1883"

#define CLIENT_KEY_FILE "/etc/ssl/certs/DST_Root_CA_X3.pem"
#define CLIENT_KEY_PASS ""
#define SERVER_KEY_FILE ""
#define CLIENT_PRIVATE_KEY_FILE ""
#define CAPATH      "/etc/ssl/certs/"

#define CLIENTID    "ExampleClientPub"
#define TOPIC       "platone/ddemo"
#define USERNAME    "ddemo"
#define PASSW       "Ddemo01!"

/* Connection settings */
#define LOWCH	4 // vPPS from RTDS GTAO
#define HIGHCH  5

/* Other settings */
// most of the settings in the algorithm file - be careful with changes
#define DT      12.5
#define SRATE   40000
#define BHANN   1199.5

/*
Currently connected:
ch0 - pin16 of RPI
ch1 - signal analog card 7
ch2 - PPS-GPS
ch3 - signal analog card 7
ch4* - signal analog card 6 (vPPS from RTDS)
ch5 - signal analog card 7

* - currently also connected to the trigger of DAQ

*/

