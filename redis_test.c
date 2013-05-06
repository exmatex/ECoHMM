/*
	ExMatEx RunTime Prototype v1.1 - Storage Benchmark Program
	Redis Test System
	Christopher Mitchell, LANL, 2013

    gcc -std=c99 redis_test.c -o redis_test -I/Users/a107908/redis-2.6.13/deps/hiredis -L/Users/a107908/redis-2.6.13/deps/hiredis -lhiredis
*/

//Includes
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "hiredis.h"

//Global Variables
int i = 0; //Loop counter var

redisContext *c;
redisReply *reply;

long numParticles; //How many particles to test with per cell
int numDoublesInParticle; //How many doubles are used to store a particle
long numCells; //How many cells in the "simulation"

//Main Function

int main(int argc, char *argv[]){
	//Get Parameters from command line
	if(argc != 4){
		printf("ERROR: Parameters Not Specified, syntax: >redis_test numParticles numDoublesInParticle numCells\n");
		printf("Exiting...\n");
		return 1;
	}

    //Ping the server - Optional sanity check on the connection
    reply = redisCommand(c, "PING");
    printf("PING: %s\n", reply->str);
    freeReplyObject(reply);
	
	numParticles = atol(argv[1]);
	numDoublesInParticle = atoi(argv[2]);
	numCells = atol(argv[3]);
	
	//Setup a connection to Redis
	struct timeval timeout = {1,500000}; //1.5 seconds
	c = redisConnectWithTimeout((char *)"127.0.0.1", 6379, timeout);
	if (c == NULL || c->err) {
    	if (c) {
            printf("Connection error: %s\n", c->errstr);
            redisFree(c);
        } else {
            printf("Connection error: can't allocate redis context\n");
        }
        exit(1);
    }
    
    //Ping the server - Optional sanity check on the connection
    reply = redisCommand(c, "PING");
    printf("PING: %s\n", reply->str);
    freeReplyObject(reply);
    
    //Setup an array based on command line parameters
    int numElements = numCells * numParticles * numDoublesInParticle;
    double * testDat = (double *) malloc(numElements * sizeof(double));
    printf("Generating Test Data...\n");
    for(i = 0; i < numElements; i++){
    	testDat[i] = (double)rand()/(double)RAND_MAX;
    	//printf("%u) %f\n", i, testDat[i]);
    }
    
    //Push the Data into Redis
    reply = redisCommand(c, "DEL testMDArray");
    freeReplyObject(reply);
    printf("Pushing Data into Redis...\n");
    for(i = 0; i < numElements; i++){
    	char tmp[sizeof(double)];
    	memcpy(&tmp, &testDat[i], sizeof(double));
    	reply = redisCommand(c, "RPUSH testMDArray %b", tmp, sizeof(double));
    	//printf("Insert [%d]: %s\n", i, reply->str);
    }
    
    //Pull the Data back from Redis, & Print it
    printf("Pull the data back from Redis & print it...\n");
    reply = redisCommand(c,"LRANGE testMDArray 0 -1"); //Get entry from element 0 -> end
    if (reply->type == REDIS_REPLY_ARRAY) {
        for (i = 0; i < reply->elements; i++) {
        	double tmp;
        	memcpy(&tmp, reply->element[i]->str, sizeof(double));
            //printf("%u) %f\n", i, tmp);
        }
    }
    freeReplyObject(reply);
    
    //Cleanup
    free (testDat);
}