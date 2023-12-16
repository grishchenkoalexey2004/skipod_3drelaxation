#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <omp.h>
#define  Max(a,b) ((a)>(b)?(a):(b))

#define N (2*2*2*2*2*2*2+2)
double   maxeps = 0.1e-7;
int itmax = 50;
int i,j,k;

double eps;
double A [N][N][N];

// вариант 40
// автор: Грищенко Алексей 

int c_num_threads = 2;

void relax();
void init();
void verify(); 


int main(int an, char **as)
{
	int it;

	double start= omp_get_wtime();
	init();

	for(it=1; it<=itmax; it++)
	{
		eps = 0.;
		relax();
		printf( "it=%4i   eps=%f\n", it,eps);
		if (eps < maxeps) break;
	}

	verify();

	printf("Время выполнения: %f\n",omp_get_wtime()-start);

	return 0;
}


void init()
{ 
	#pragma omp parallel for private(k,j,i) num_threads(c_num_threads) collapse(3)
	for(k=0; k<=N-1; k++)
	for(j=0; j<=N-1; j++)
	for(i=0; i<=N-1; i++)
	{
		if(i==0 || i==N-1 || j==0 || j==N-1 || k==0 || k==N-1)
			A[i][j][k]= 0.;
		else A[i][j][k]= ( 4. + i + j + k) ;
	}
} 

void relax()
{

	#pragma omp parallel num_threads(c_num_threads)
	{
		int iam = omp_get_thread_num();
		for (int newi=1;newi<=N-2+c_num_threads-1;newi++){
			int i = newi-iam;
			#pragma omp for private(k,j)
			for(k=1; k<=N-2; k++)
			for(j=1; j<=N-2; j++)
			{
				if ((i>=1) && (i<=N-2)){
					A[i][j][k] = (A[i-1][j][k]+A[i+1][j][k])/2.;
				}
			}
		}
	}
	
	
	#pragma omp parallel num_threads(c_num_threads)
	{
		int iam = omp_get_thread_num();
		for (int newj=1;newj<=N-2+c_num_threads-1;newj++){
			int j = newj-iam;
			#pragma omp for private(k,i)
			for(k=1; k<=N-2; k++)
			for(i=1; i<=N-2; i++)
			{
				if ((j>=1) && (j<=N-2)){
					A[i][j][k] =(A[i][j-1][k]+A[i][j+1][k])/2.;
				}
			}

		}	
	}

	#pragma omp parallel num_threads(c_num_threads)
	{
		int iam = omp_get_thread_num();
		for(int newk=1; newk<=N-2+c_num_threads-1; newk++){
			int k = newk-iam;
			#pragma omp for private(i,j)
			for(j=1; j<=N-2; j++)
			for(i=1; i<=N-2; i++)
			{
				if ((k>=1) && (k<=N-2)){
					double e;
					e=A[i][j][k];
					A[i][j][k] = (A[i][j][k-1]+A[i][j][k+1])/2.;
					eps=Max(eps,fabs(e-A[i][j][k]));
					#pragma omp flush(eps)
				}
				
			}
		}
	}	
		
	
	
}

void verify()
{
	double s;

	s=0.;

	#pragma omp parallel for private(k,j,i) reduction(+:s) num_threads(c_num_threads) collapse(3)
	for(k=0; k<=N-1; k++)
	for(j=0; j<=N-1; j++)
	for(i=0; i<=N-1; i++)
	{
		s=s+A[i][j][k]*(i+1)*(j+1)*(k+1)/(N*N*N);
	}

	printf("  S = %f\n",s);

}
//S = 1172094.718571

