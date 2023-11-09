#include "nr.h"

extern int mm,n,mpt;
extern DP h,c2,anorm;
extern Vec_DP *x_p;

void NR::difeq(const int k, const int k1, const int k2, const int jsf,
	const int is1, const int isf, Vec_I_INT &indexv, Mat_O_DP &s,
	Mat_I_DP &y)
{
	DP temp,temp1,temp2;

	Vec_DP &x=*x_p;
	if (k == k1) {
		if ((n+mm & 1) != 0) {
			s[2][3+indexv[0]]=1.0;
			s[2][3+indexv[1]]=0.0;
			s[2][3+indexv[2]]=0.0;
			s[2][jsf]=y[0][0];
		} else {
			s[2][3+indexv[0]]=0.0;
			s[2][3+indexv[1]]=1.0;
			s[2][3+indexv[2]]=0.0;
			s[2][jsf]=y[1][0];
		}
	} else if (k > k2-1) {
		s[0][3+indexv[0]] = -(y[2][mpt-1]-c2)/(2.0*(mm+1.0));
		s[0][3+indexv[1]]=1.0;
		s[0][3+indexv[2]] = -y[0][mpt-1]/(2.0*(mm+1.0));
		s[0][jsf]=y[1][mpt-1]-(y[2][mpt-1]-c2)*y[0][mpt-1]/
			(2.0*(mm+1.0));
		s[1][3+indexv[0]]=1.0;
		s[1][3+indexv[1]]=0.0;
		s[1][3+indexv[2]]=0.0;
		s[1][jsf]=y[0][mpt-1]-anorm;
	} else {
		s[0][indexv[0]] = -1.0;
		s[0][indexv[1]] = -0.5*h;
		s[0][indexv[2]]=0.0;
		s[0][3+indexv[0]]=1.0;
		s[0][3+indexv[1]] = -0.5*h;
		s[0][3+indexv[2]]=0.0;
		temp1=x[k]+x[k-1];
		temp=h/(1.0-temp1*temp1*0.25);
		temp2=0.5*(y[2][k]+y[2][k-1])-c2*0.25*temp1*temp1;
		s[1][indexv[0]]=temp*temp2*0.5;
		s[1][indexv[1]] = -1.0-0.5*temp*(mm+1.0)*temp1;
		s[1][indexv[2]]=0.25*temp*(y[0][k]+y[0][k-1]);
		s[1][3+indexv[0]]=s[1][indexv[0]];
		s[1][3+indexv[1]]=2.0+s[1][indexv[1]];
		s[1][3+indexv[2]]=s[1][indexv[2]];
		s[2][indexv[0]]=0.0;
		s[2][indexv[1]]=0.0;
		s[2][indexv[2]] = -1.0;
		s[2][3+indexv[0]]=0.0;
		s[2][3+indexv[1]]=0.0;
		s[2][3+indexv[2]]=1.0;
		s[0][jsf]=y[0][k]-y[0][k-1]-0.5*h*(y[1][k]+y[1][k-1]);
		s[1][jsf]=y[1][k]-y[1][k-1]-temp*((x[k]+x[k-1])
			*0.5*(mm+1.0)*(y[1][k]+y[1][k-1])-temp2
			*0.5*(y[0][k]+y[0][k-1]));
		s[2][jsf]=y[2][k]-y[2][k-1];
	}
}
