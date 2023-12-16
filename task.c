#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <omp.h>

// вариант 40
// автор: Грищенко Алексей
// параллелизация с помощью task 
// один из вариантов распаралелленной программы

#define  Max(a,b) ((a)>(b)?(a):(b))

#define N (2*2*2*2*2*2*2*2+2)

double   maxeps = 0.1e-7;
int itmax = 20;
int i,j,k;

double eps;
double start,end;

int n_threads;
int tid;

double A[N][N][N];

void relax();
void init();
double verify(); 


int main(int an, char **as)
{
	int it;
	double result; 

	start = omp_get_wtime();

	init();


	for(it=1; it<=itmax; it++)
	{
		eps = 0.;
		relax();
		if (eps < maxeps) break;
	}

	result = verify();

	end = omp_get_wtime();

	printf("%d %d %f %f\n",N,n_threads,end-start,result);
	return 0;
}


void init()
{ 
	#pragma omp parallel private(i,j,k)
	{
		#pragma omp single
		{
			for(i=0; i<=N-1; i++)
			#pragma omp task untied 
			{
				for(j=0; j<=N-1; j++)
				for(k=0; k<=N-1; k++)
				{
					if(i==0 || i==N-1 || j==0 || j==N-1 || k==0 || k==N-1)
						A[i][j][k]= 0.;
					else A[i][j][k]= ( 4. + i + j + k) ;
				}
			}
			
		}
		
	}
	
} 

void relax()
{


	#pragma omp parallel private(i,j,k)
	{
		if (omp_get_thread_num()==0){
			n_threads = omp_get_num_threads();
		}

		#pragma omp single
		{
			for(i=1; i<=N-2; i++){
				#pragma omp taskloop
				for(j=1; j<=N-2; j++)
				{
					for(k=1; k<=N-2; k++)
					{
						A[i][j][k] = (A[i-1][j][k]+A[i+1][j][k])/2.;
					}
				}

				#pragma omp taskwait
			}
			
		}



		#pragma omp single
		{
			for(j=1; j<=N-2; j++)
			{
				#pragma omp taskloop
				for(i=1; i<=N-2; i++)
				for(k=1; k<=N-2; k++)
				
				{
					A[i][j][k] =(A[i][j-1][k]+A[i][j+1][k])/2.;
				}

				#pragma omp taskwait
			}	
			
		}



		#pragma omp single
		{
			for(k=1; k<=N-2; k++)
			{	
				#pragma omp taskloop
				for(i=1; i<=N-2; i++)
				for(j=1; j<=N-2; j++)
				{	
					double e;
					e=A[i][j][k];
					A[i][j][k] = (A[i][j][k-1]+A[i][j][k+1])/2.;
					eps=Max(eps,fabs(e-A[i][j][k]));
				}
				#pragma omp taskwait
			}	
		}
		
		


	}
}


double verify()
{
	double s;

	s=0.;
	
	#pragma omp parallel private(i,j,k)
	{
		#pragma omp single
		{
			for(i=0; i<=N-1; i++)
			#pragma omp task untied
			for(j=0; j<=N-1; j++)
			for(k=0; k<=N-1; k++)
			{	
				#pragma omp atomic update
					s+=A[i][j][k]*(i+1)*(j+1)*(k+1)/(N*N*N);
				
			}

		}
	}
	
	return s;

}


