#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include "hiredis.h"

redisContext *c;
redisReply *reply;

int main(int argc, char* argv[])
{
	int     option = 0;
	float   x = 0.0, y = 0.0, z = 0.0;

/*
   On startup, check for database. That means we're being
   restarted. Open the database, get our last timestep (and
   its particles) and contine with subsequent timestep (recovery).
*/ 
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


	scanf("%f %f %f", &x, &y, &z);
	printf("%f %f %f\n", x, y, z);
    fflush(stdout);
    sleep(20);
    return(-22);
}