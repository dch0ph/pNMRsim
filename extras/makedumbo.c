#include <stdio.h>
#include <math.h>

const double TWOPI=2*M_PI;

#define N 6

int main()
{
  const int steps=32;

  /* Original DUMBO-1 coefficients */
  const double as[7]={0.1506,0.0325,0.0189,0.0238,0.0107, 0.0038, -0.0013 };
  const double bs[7]={0.0000,0.1310,0.1947,0.0194,0.1124,-0.0456, 0.0869};

  double phis[steps];
  int i,j;

  for (i=0;i<steps;i++) {
	double phi=0.0;
	for (j=0;j<=N;j++)
	  phi+=as[j]*cos((TWOPI*j*i)/steps)+bs[j]*sin((TWOPI*j*i)/steps);

	phis[i]=360.0*phi;
  }

  for (i=0;i<steps;i++)
	printf("%g\n",phis[i]);
  for (i=0;i<steps;i++)
	printf("%g\n",phis[steps-i]+180.0);
     
  return 0;
}
