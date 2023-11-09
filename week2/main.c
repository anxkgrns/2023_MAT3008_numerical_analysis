#include "nr.h"

void newton(float a, float *b, float *c){
	*b = bessj0(a);
	*c = -bessj1(a);
}

float func(float x){
	return 2*x*x+5*x+3;
}

void newton2(float a, float *b, float *c){
	*b = func(a);
	*c = 4*a+5;
}



int main()
{
	float x1=1.0,x2=10.0;
	int n = 30;
	float xb1[100]={0,},xb2[100]={0,};
	int nb=0;
	int convergence_speed = 0;

	zbrak(bessj0,x1,x2,n,xb1,xb2,&nb);
	printf("%d개의 root 존재\n",nb);
	for (int i = 1; i <= nb; i++) {
        printf("영점 #%d: [%.2f, %.2f]\n", i, xb1[i], xb2[i]);
    }
	printf("\nBisection:\n");
	for (int j = 1; j <= nb ; j++){
		convergence_speed = 0;
		float root = rtbis(bessj0,xb1[j],xb2[j],0.000001,&convergence_speed);
		printf("root%d : %f / convergence speed : %d\n",j,root,convergence_speed);
	}
	printf("\nLinear interpolation\n");
	for (int j = 1; j <= nb ; j++){
		convergence_speed = 0;
		float root = rtflsp(bessj0,xb1[j],xb2[j],0.000001,&convergence_speed);
		printf("root%d : %f / convergence speed : %d\n",j,root,convergence_speed);
	}
	printf("\nSecant\n");
	for (int j = 1; j <= nb ; j++){
		convergence_speed = 0;
		float root = rtsec(bessj0,xb1[j],xb2[j],0.000001,&convergence_speed);
		printf("root%d : %f / convergence speed : %d\n",j,root,convergence_speed);
	}
	printf("\nNewton-Raphson\n");
	for (int j = 1; j <= nb ; j++){
		convergence_speed = 0;
		float root = rtnewt(newton,xb1[j],xb2[j],0.000001,&convergence_speed);
		printf("root%d : %f / convergence speed : %d\n",j,root,convergence_speed);
	}
	printf("\nNewton with bracketing\n");
	for (int j = 1; j <= nb ; j++){
		convergence_speed = 0;
		float root = rtsafe(newton,xb1[j],xb2[j],0.000001,&convergence_speed);
		printf("root%d : %f / convergence speed : %d\n",j,root,convergence_speed);
	}
	printf("\nMuller\n");
	for (int j = 1; j <= nb ; j++){
		convergence_speed = 0;
		float root = muller(bessj0,xb1[j],xb2[j],0.000001,&convergence_speed);
		printf("root%d : %f / convergence speed : %d\n",j,root,convergence_speed);
	}

	printf("\nWith function : 2*x*x+5*x+3\n");

	float x3=-10.0,x4=1.0;

	float xb3[100]={0,},xb4[100]={0,};
	nb=0;
	zbrak(func,x3,x4,n,xb3,xb4,&nb);

	printf("%d개의 root 존재\n",nb);
	for (int i = 1; i <= nb; i++) {
        printf("영점 #%d: [%.2f, %.2f]\n", i, xb3[i], xb4[i]);
    }
	printf("\nBisection:\n");
	for (int j = 1; j <= nb ; j++){
		convergence_speed = 0;
		float root = rtbis(func,xb3[j],xb4[j],0.000001,&convergence_speed);
		printf("root%d : %f / convergence speed : %d\n",j,root,convergence_speed);
	}
	printf("\nLinear interpolation\n");
	for (int j = 1; j <= nb ; j++){
		convergence_speed = 0;
		float root = rtflsp(func,xb3[j],xb4[j],0.000001,&convergence_speed);
		printf("root%d : %f / convergence speed : %d\n",j,root,convergence_speed);
	}
	printf("\nSecant\n");
	for (int j = 1; j <= nb ; j++){
		convergence_speed = 0;
		float root = rtsec(func,xb3[j],xb4[j],0.000001,&convergence_speed);
		printf("root%d : %f / convergence speed : %d\n",j,root,convergence_speed);
	}
	printf("\nNewton-Raphson\n");
	for (int j = 1; j <= nb ; j++){
		convergence_speed = 0;
		float root = rtnewt(newton2,xb3[j],xb4[j],0.000001,&convergence_speed);
		printf("root%d : %f / convergence speed : %d\n",j,root,convergence_speed);
	}
	printf("\nNewton with bracketing\n");
	for (int j = 1; j <= nb ; j++){
		convergence_speed = 0;
		float root = rtsafe(newton2,xb3[j],xb4[j],0.000001,&convergence_speed);
		printf("root%d : %f / convergence speed : %d\n",j,root,convergence_speed);
	}
	printf("\nMuller\n");
	for (int j = 1; j <= nb ; j++){
		convergence_speed = 0;
		float root = muller(func,xb3[j],xb4[j],0.000001,&convergence_speed);
		printf("root%d : %f / convergence speed : %d\n",j,root,convergence_speed);
	}
	return 0;
}


