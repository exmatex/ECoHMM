#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>

int main(int argc, char* argv[])
{
	int     option = 0;
	float   x = 0.0, 
			y = 0.0, 
			z = 0.0;

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


	//scanf("%f %f %f", &x, &y, &z);
	printf("%f %f %f\n", x, y, z);
}