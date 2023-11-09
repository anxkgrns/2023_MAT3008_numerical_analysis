#include "nr.h"
#include <ctype.h>
#include <stdlib.h>

int m = 1;

void sol_guass(float ** A,float ** B, int n){
	printf("1. Guass-Jordan\n");
	gaussj(A, n, B, m);
	for(int i = 1; i <= n; i++) {
		printf("x%d = %f ", i, B[i][1]);
	}
	printf("\n\n");
}

void use_guass_lud(char* filename, float ** A,float ** B, int n){
	printf("%s - Inverse matrix\n", filename);
	float **A_1 = malloc((n+1) * sizeof(float *));
	float **B_1 = malloc((n+1) * sizeof(float *));
	for (int i = 1; i <= n; i++) {
		A_1[i] = malloc((n+1) * sizeof(float));
		B_1[i] = malloc((m+1) * sizeof(float));
	}
	for(int i = 1; i <= n; i++) {
		for(int j = 1; j <= n; j++) {
			A_1[i][j] = A[i][j];
		}
		for(int j = 1; j <= m; j++) {
			B_1[i][j] = B[i][j];
		}
	}


	gaussj(A, n, B, m);
	for(int i = 1; i <= n; i++) {
		for(int j = 1; j <= n; j++) {
			printf("%f    ",A[i][j]);
		}
		printf("\n");
	}
	printf("\n");

	float d = 0;
	int *indx = malloc((n+1) * sizeof(int));
	ludcmp(A_1, n, indx,&d);
	float determinant = 1;
	for(int i = 1; i <= n; i++) {
		determinant *= A_1[i][i];
	}

	printf("determinant :  %f\n\n",determinant);

	for(int i = 1; i <= n; i++) {
		free(A_1[i]);
		free(B_1[i]);
	}
	free(A_1);
	free(B_1);
	free(indx);
}

void sol_lud(float ** A,float ** B, int n){
	printf("2.LU Decomposition\n");

	float d = 0;
	int *indx = malloc((n+1) * sizeof(int));
	ludcmp(A, n, indx,&d);
	float * B_lu = malloc((n+1) * sizeof(float));
	for (int i = 1; i <= n; i++) {
		B_lu[i] = B[i][m];
	}
	lubksb(A, n, indx, B_lu);

	for(int i = 1; i <= n; i++) {
		printf("x%d = %f ", i, B_lu[i]);
	}
	free(indx);
	free(B_lu);
	printf("\n\n");
}

void use_lud(float ** A,float ** B, int n, float * x){
	float d = 0;
	int *indx = malloc((n+1) * sizeof(int));
	ludcmp(A, n, indx,&d);
	float * B_lu = malloc((n+1) * sizeof(float));
	for (int i = 1; i <= n; i++) {
		B_lu[i] = B[i][m];
	}
	lubksb(A, n, indx, B_lu);
	for(int i = 1; i <= n; i++) {
		x[i] = B_lu[i];
	}
	free(indx);
	free(B_lu);
}

void sol_svd(float ** A,float ** B, int n){
	printf("3. Singular Value Decomposition\n");

	float *w = malloc((n+1) * sizeof(float));
	float **v = malloc((n+1) * sizeof(float *));
	for (int i = 1; i <= n; i++) {
		v[i] = malloc((n+1) * sizeof(float));
	}
	svdcmp(A, n, n, w, v);
	float * x_sv = malloc((n+1) * sizeof(float));
	float * B_sv = malloc((n+1) * sizeof(float));
	for (int i = 1; i <= n; i++) {
		B_sv[i] = B[i][m];
	}
	svbksb(A, w, v, n, n, B_sv, x_sv);
	
	for(int i = 1; i <= n; i++) {
		printf("x%d = %f ", i, x_sv[i]);
	}
	free(w);
	free(B_sv);
	free(x_sv);
	for(int i = 1; i <= n; i++) {
		free(v[i]);
	}
	free(v);
	printf("\n\n");
}

void sol_mprove(char* filename, float ** A,float ** B, int n){
	printf("%s - iterative improvement\n", filename);
	float d = 0;
	int *indx = malloc((n+1) * sizeof(int));
	float * B_mp = malloc((n+1) * sizeof(float));
	for (int i = 1; i <= n; i++) {
		B_mp[i] = B[i][m];
	}
	float ** A_lu = malloc((n+1) * sizeof(float *));
	for (int i = 1; i <= n; i++) {
		A_lu[i] = malloc((n+1) * sizeof(float));
	}
	for(int i = 1; i <= n; i++) {
		for(int j = 1; j <= n; j++) {
			A_lu[i][j] = A[i][j];
		}
	}
	float * x_mp = malloc((n+1) * sizeof(float));
	use_lud(A_lu, B, n, x_mp);

	mprove(A, A_lu, n, indx, B_mp, x_mp);
	

	for(int i = 1; i <= n; i++) {
		printf("x%d = %f ", i, x_mp[i]);
	}
	free(indx);
	free(B_mp);
	printf("\n\n");
}

int main()
{

	FILE *fp = NULL;
	// #1
	int n; // n
	float c_1 = 0;
	int i = 1;
	int j = 1;
	char** filename = malloc(3 * sizeof(char *));
	filename[0] = "lineq1.dat";
	filename[1] = "lineq2.dat";
	filename[2] = "lineq3.dat";

	// #1
	for (int k = 0; k < 3; k++) {
		fp = fopen(filename[k], "r");
		if (fp == NULL) {
			fprintf(stderr, "파일 %s를 열 수 없습니다.\n", filename[k]);
			return 0;
		}
		fscanf(fp, "%d %d", &n, &n);
		float **A = malloc((n+1) * sizeof(float *));
		float **B = malloc((n+1) * sizeof(float *));
		for (int i = 1; i <= n; i++) {
			A[i] = malloc((n+1) * sizeof(float));
			B[i] = malloc((m+1) * sizeof(float));
		}
		i = 1;
		j = 1;
		while(fscanf(fp, "%f", &c_1) != EOF) {
			if(i == n+1) {
				B[j][m] = c_1;
				j++;
				continue;
			}
			A[i][j] = c_1;
			j++;
			if (j == n+1) {
				i++;
				j = 1;
			}
		}
		//copy A to A_copy1,2,3
		float **A_copy1 = malloc((n+1) * sizeof(float *));
		float **A_copy2 = malloc((n+1) * sizeof(float *));
		float **A_copy3 = malloc((n+1) * sizeof(float *));
		
		for (int i = 1; i <= n; i++) {
			A_copy1[i] = malloc((n+1) * sizeof(float));
			A_copy2[i] = malloc((n+1) * sizeof(float));
			A_copy3[i] = malloc((n+1) * sizeof(float));
		}
		for (i = 1; i <= n; i++) {
			for (j = 1; j <= n; j++) {
				A_copy1[i][j] = A[i][j];
				A_copy2[i][j] = A[i][j];
				A_copy3[i][j] = A[i][j];
			}
		}
		//copy B to B_copy1,2,3
		float **B_copy1 = malloc((n+1) * sizeof(float *));
		float **B_copy2 = malloc((n+1) * sizeof(float *));
		float **B_copy3 = malloc((n+1) * sizeof(float *));
		for (int i = 1; i <= n; i++) {
			B_copy1[i] = malloc((m+1) * sizeof(float));
			B_copy2[i] = malloc((m+1) * sizeof(float));
			B_copy3[i] = malloc((m+1) * sizeof(float));
		}
		for (i = 1; i <= n; i++) {
			for (j = 1; j <= m; j++) {
				B_copy1[i][j] = B[i][j];
				B_copy2[i][j] = B[i][j];
				B_copy3[i][j] = B[i][j];
			}
		}
		printf("----------------------------------------\n\n");
		printf("File: %s\n\n", filename[k]);
		sol_guass(A_copy1, B_copy1, n);
		sol_lud(A_copy2, B_copy2, n);
		sol_svd(A_copy3, B_copy3, n);
		for(int i = 1; i <= n; i++) {
			free(A[i]);
			free(A_copy1[i]);
			free(A_copy2[i]);
			free(A_copy3[i]);
			free(B[i]);
			free(B_copy1[i]);
			free(B_copy2[i]);
			free(B_copy3[i]);
		}
		free(A);
		free(A_copy1);
		free(A_copy2);
		free(A_copy3);
		free(B);
		free(B_copy1);
		free(B_copy2);
		free(B_copy3);
		if(k == 2) printf("----------------------------------------\n");
		fclose(fp);
	}

	// #2
	for (int k = 0; k < 3; k++) {
		printf("\n");
		fp = fopen(filename[k], "r");
		if (fp == NULL) {
			fprintf(stderr, "파일 %s를 열 수 없습니다.\n", filename[k]);
			return 0;
		}
		fscanf(fp, "%d %d", &n, &n);
		float **A = malloc((n+1) * sizeof(float *));
		float **B = malloc((n+1) * sizeof(float *));
		for (int i = 1; i <= n; i++) {
			A[i] = malloc((n+1) * sizeof(float));
			B[i] = malloc((m+1) * sizeof(float));
		}
		i = 1;
		j = 1;
		while(fscanf(fp, "%f", &c_1) != EOF) {
			if(i == n+1) {
				B[j][m] = c_1;
				j++;
				continue;
			}
			A[i][j] = c_1;
			j++;
			if (j == n+1) {
				i++;
				j = 1;
			}
		}
		//copy A to A_copy1,2,3
		float **A_copy1 = malloc((n+1) * sizeof(float *));
		
		for (int i = 1; i <= n; i++) {
			A_copy1[i] = malloc((n+1) * sizeof(float));
		}
		for (i = 1; i <= n; i++) {
			for (j = 1; j <= n; j++) {
				A_copy1[i][j] = A[i][j];
			}
		}
		//copy B to B_copy1,2,3
		float **B_copy1 = malloc((n+1) * sizeof(float *));
		for (int i = 1; i <= n; i++) {
			B_copy1[i] = malloc((m+1) * sizeof(float));
		}
		for (i = 1; i <= n; i++) {
			for (j = 1; j <= m; j++) {
				B_copy1[i][j] = B[i][j];
			}
		}
		sol_mprove(filename[k],A_copy1, B_copy1, n);
		for(int i = 1; i <= n; i++) {
			free(A[i]);
			free(A_copy1[i]);
			free(B[i]);
			free(B_copy1[i]);
		}
		free(A);
		free(A_copy1);
		free(B);
		free(B_copy1);
		if(k == 2) printf("----------------------------------------\n");
		fclose(fp);
	}

	// #3
	for (int k = 0; k < 3; k++) {
		printf("\n");
		fp = fopen(filename[k], "r");
		if (fp == NULL) {
			fprintf(stderr, "파일 %s를 열 수 없습니다.\n", filename[k]);
			return 0;
		}
		fscanf(fp, "%d %d", &n, &n);
		float **A = malloc((n+1) * sizeof(float *));
		float **B = malloc((n+1) * sizeof(float *));
		for (int i = 1; i <= n; i++) {
			A[i] = malloc((n+1) * sizeof(float));
			B[i] = malloc((m+1) * sizeof(float));
		}
		i = 1;
		j = 1;
		while(fscanf(fp, "%f", &c_1) != EOF) {
			if(i == n+1) {
				B[j][m] = c_1;
				j++;
				continue;
			}
			A[i][j] = c_1;
			j++;
			if (j == n+1) {
				i++;
				j = 1;
			}
		}
		//copy A to A_copy1,2,3
		float **A_copy1 = malloc((n+1) * sizeof(float *));
		
		for (int i = 1; i <= n; i++) {
			A_copy1[i] = malloc((n+1) * sizeof(float));
		}
		for (i = 1; i <= n; i++) {
			for (j = 1; j <= n; j++) {
				A_copy1[i][j] = A[i][j];
			}
		}
		//copy B to B_copy1,2,3
		float **B_copy1 = malloc((n+1) * sizeof(float *));
		for (int i = 1; i <= n; i++) {
			B_copy1[i] = malloc((m+1) * sizeof(float));
		}
		for (i = 1; i <= n; i++) {
			for (j = 1; j <= m; j++) {
				B_copy1[i][j] = B[i][j];
			}
		}
		use_guass_lud(filename[k],A_copy1, B_copy1, n);
		
		for(int i = 1; i <= n; i++) {
			free(A[i]);
			free(A_copy1[i]);
			free(B[i]);
			free(B_copy1[i]);
		}
		free(A);
		free(A_copy1);
		free(B);
		free(B_copy1);
		if(k == 2) printf("----------------------------------------\n");
		fclose(fp);
	}
	
	return 0;
}


