#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>

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
    double  stress[50];

    for (int i=0; i<50; i++) {
	   fscanf(stdin, "%lf", &stress[i]);
    }

	printf("stress:\n");
    for (int i=0; i<50; i++) {
        printf("%2d: %lf\n", i, stress[i]);
    }
    return(-22);
}