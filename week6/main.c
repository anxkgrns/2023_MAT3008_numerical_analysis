#include "nr.h"
#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>


// vector<float> a(3,0);
// vector<float> b(3,0);

int main()
{
	char** filename = malloc(3 * sizeof(char *));
	filename[0] = "fitdata1.dat";
	filename[1] = "fitdata2.dat";
	filename[2] = "fitdata3.dat";
	FILE *fp = NULL;



	for(int k =0; k<3; k++){
		fp = fopen(filename[k], "r");
		printf("\n		|%s|\n", filename[k]);
		printf("-----------------------------------------------\n");
		if (fp == NULL) {
			printf("Error opening file!\n");
			exit(1);
		}
		float data[4][100] = {0};
		float num;
		int line = 0, temp = 0;
		int data_size = 1;
		while(fscanf(fp, "%f", &num) != EOF){
			//printf("!!!\n");

			data[temp][line] = num;
			temp++;

			if(temp >= 4){
				line++;
				data_size++;
				temp = 0;
			}
		}

		float X_prime[100] = {0};
		for(int i = 0; i < data_size; i++){
			X_prime[i] = data[2][i];
		}

		float Y_prime[100] = {0};
		for(int i = 0; i < data_size; i++){
			Y_prime[i] = data[3][i];
		}

		float F[3][100] = {0};
		for(int i = 0; i < data_size; i++){
			for(int j = 0; j < 2; j++){
				F[j][i] = data[j][i];
			}
			F[2][i] = 1;
		}

		//F^T * F * a1 = F^T * X_prime

		// A1 = F^T * F
		float ** A1 = malloc(4 * sizeof(float *));
		for(int i = 0; i < 4; i++){
			A1[i] = malloc(4 * sizeof(float));
		}

		for(int i = 0; i < 3; i++){
			for(int j = 0; j < 3; j++){
				float temp = 0;
				for(int k = 0; k < data_size; k++){
					temp += F[j][k] * F[i][k];
				}
				A1[i+1][j+1] = temp;
			}
		}

		// b1 = F^T * X_prime
		float ** b1 = malloc(4 * sizeof(float *));
		for(int i = 0; i < 4; i++){
			b1[i] = malloc(2 * sizeof(float));
		}

		for(int i = 0; i < 3; i++){
			float temp = 0;
			for(int j = 0; j < data_size; j++){
				temp += F[i][j] * X_prime[j];
			}
			b1[i+1][1] = temp;
		}
		// A1 * a1 = b1
		gaussj(A1, 3, b1, 1);

		for(int i=1;i<4;i++){
			printf(" a%d : %.6f ", i, b1[i][1]);
			if(i!=3)printf("|");
		}
		printf("\n");

		//F^T * F * a2 = F^T * Y_prime 

		// A2 = F^T * F
		float ** A2 = malloc(4 * sizeof(float *));
		for(int i = 0; i < 4; i++){
			A2[i] = malloc(4 * sizeof(float));
		}

		for(int i = 0; i < 3; i++){
			for(int j = 0; j < 3; j++){
				float temp = 0;
				for(int k = 0; k < data_size; k++){
					temp += F[j][k] * F[i][k];
				}
				A2[i+1][j+1] = temp;
			}
		}

		// b2 = F^T * X_prime
		float ** b2 = malloc(4 * sizeof(float *));
		for(int i = 0; i < 4; i++){
			b2[i] = malloc(2 * sizeof(float));
		}

		for(int i = 0; i < 3; i++){
			float temp = 0;
			for(int j = 0; j < data_size; j++){
				temp += F[i][j] * Y_prime[j];
			}
			b2[i+1][1] = temp;
		}
		// A1 * a1 = b1
		gaussj(A2, 3, b2, 1);

		for(int i=1;i<4;i++){
			printf(" a%d : %2.6f ", i+3, b2[i][1]);
			if(i!=3)printf("|");
		}
		printf("\n");


		printf("-----------------------------------------------\n");

		fclose(fp);
	}

}