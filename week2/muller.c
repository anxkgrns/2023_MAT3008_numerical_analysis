
#include <math.h>
#include<stdio.h>
#define MAXIT 30
int sgn(float x) {
    if (x > 0.0) {
        return 1;  // 양수인 경우 1 반환
    } else if (x < 0.0) {
        return -1; // 음수인 경우 -1 반환
    } else {
        return 0;  // 0인 경우 0 반환
    }
}


float muller(float (*func)(float), float x0, float x2, float xacc,int* convergence_speed)
{
	void nrerror(char error_text[]);
	int j;
	float f0,f1,f2;
	float p0,p1,p2,p3,a,b,c;
	float x1;
	x1 = (x0+x2)/2.0;
	p0 = x0;
	p1 = x1;
	p2 = x2;

	for (j=1;j<=MAXIT;j++) {
		(*convergence_speed)++;
		f0=(*func)(p0);
		f1=(*func)(p1);
		f2=(*func)(p2);

		c=f2;
		float mom = (p0-p2)*(p1-p2)*(p0-p1);
		b=(((p0-p2)*(p0-p2)*(f1-f2))-((p1-p2)*(p1-p2)*(f0-f2)))/mom;
		a=(((p1-p2)*(f0-f2))-((p0-p2)*(f1-f2)))/mom;
	
		p3=p2-(2*c)/(b+sgn(b)*sqrt(b*b-4*a*c));		

		if (fabs(p2-p3) < xacc) return p3;
		p0=p1;
		p1=p2;
		p2=p3;
	}
	nrerror("Maximum number of iterations exceeded in muller");
	return 0.0;
}
#undef MAXIT
