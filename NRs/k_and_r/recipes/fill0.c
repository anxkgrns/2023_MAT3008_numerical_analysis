
void fill0(u,n)
double **u;
int n;
{
	int i,j;
	for (j=1;j<=n;j++)
		for (i=1;i<=n;i++)
			u[i][j]=0.0;
}
