#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
// #include "hiredis.h"

// redisContext *c;
// redisReply *reply;


/*
  Driven by input and output on stdin/stdout.

  1. Start with command line parameters.
  2. Receive (read) stress from Erlang.
  3. Calculate strain (sleep) and save particles to DB.
  4. Send strain to Erlang server.
  5. Go to 2.
*/
int main(int argc, char* argv[])
{
	int     option = 0;
	float   x = 0.0, y = 0.0, z = 0.0;

/*
   On startup, check for database. That means we're being
   restarted. Open the database, get our last timestep (and
   its particles) and contine with subsequent timestep (recovery).

    struct timeval timeout = {1,500000}; //1.5 seconds
    c = redisConnectWithTimeout((char *)"127.0.0.1", 6379, timeout);
    if (c == NULL || c->err) {
        if (c) {
            printf("Connection error: %s\n", c->errstr);
            redisFree(c);
        } else {
            printf("Connection error: can't allocate redis context\n");
        }
        return(0);
    }
*/
    while ((option = getopt(argc, argv,"x:y:z:")) != -1) {
        switch (option) {
            case 'x' : 
             	x = strtod(optarg, NULL); 
                break;
            case 'y' : 
             	y = strtod(optarg, NULL); 
                break;
            case 'z' : 
             	z = strtod(optarg, NULL); 
                break;
             default: 
                exit(EXIT_FAILURE);
        }
    }

    float  stress_x, stress_y;
    float  strain_x, strain_y;

	scanf("%f %f", &stress_x, &stress_y);
    sleep(2);
	printf("strain %f %f\n", 1.0/stress_x, 1.0/stress_y);
    fflush(stdout);
    sleep(20);
    return(-22);
}