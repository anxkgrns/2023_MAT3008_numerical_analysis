#include "nr.h"
#include <ctype.h>
#include <stdlib.h>
#include <time.h>


int a = -3;
int b = 4;
float m = 0.5;
float s = 1.5;

int sample[4] = {100,1000,10000,100000};

int main()
{
	long seed = time(NULL);
	char **filename = malloc(8 * sizeof(char *));
	filename[0] = "guass_sample_100.txt";
	filename[1] = "guass_sample_1000.txt";
	filename[2] = "guass_sample_10000.txt";
	filename[3] = "guass_sample_100000.txt";
	filename[4] = "uniform_sample_100.txt";
	filename[5] = "uniform_sample_1000.txt";
	filename[6] = "uniform_sample_10000.txt";
	filename[7] = "uniform_sample_100000.txt";

	for(int i =0; i<8; i++){
		FILE *fp = fopen(filename[i], "w");
		if (fp == NULL) {
			printf("Error opening file!\n");
			exit(1);
		}
		if(i<4){
			for(int j = 0; j<sample[i]; j++){
				fprintf(fp, "%f\n", gasdev(&seed)*s+m);
			}
		}else{
			for(int j = 0; j<sample[i-4]; j++){
				fprintf(fp, "%f\n", ran1(&seed)*(b-a)+a);
			}
		}
		fclose(fp);
	}

}