#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <omp.h>

// параллелизация с помощью for 
// один из вариантов распаралелленной программы


// Для загрузки фаила
// scp  var40_variations.c edu-cmc-skpod23-324-03@polus.hpc.cs.msu.ru:/home_edu/edu-cmc-skpod23-324/edu-cmc-skpod23-324-03/var40_variations.c


// Для компиляции фаила: 
// gcc -std=c99 -fopenmp var40_variations.c

// Для запуска фаила 
// mpisubmit.pl --stdout=results.txt ./a.out

#define  Max(a,b) ((a)>(b)?(a):(b))


double   maxeps = 0.1e-7;
int itmax = 35;
int i,j,k;

double eps;
double start,end;

double*** A;

int matr_dims[4] = {2*2*2*2*2+2,2*2*2*2*2*2+2,2*2*2*2*2*2*2+2,2*2*2*2*2*2*2*2+2};
int num_threads[12] = {1,2,3,4,6,8,10,20,40,80,120,160};


int cur_dim;
int c_num_threads;

void relax();
void init();
double verify(); 
void delete_data();
void create_data();


int main(int an, char **as)
{
	
	int it;
	double result; 




	for (int f =0;f<4;f++){

		cur_dim = matr_dims[f];

		create_data();
		for (int m =0; m<12;m++){


			c_num_threads = num_threads[m];

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
			
			printf("%d   %d   %f    %f\n",cur_dim,c_num_threads,end-start,result);

		}

		delete_data();
	}




	return 0;
}

void create_data(){
	A = (double ***) malloc(cur_dim * sizeof(double**));
    for(int l= 0; l < cur_dim; l++){
        A[l] = (double **)malloc(cur_dim * sizeof(double*));
		for (int q =0;q<cur_dim;q++){
			A[l][q] = (double *) malloc(cur_dim * sizeof(double));
		}
    }
}



void init()
{ 

	#pragma omp parallel for private(k,j,i) num_threads(c_num_threads) collapse(3)
	for(k=0; k<=cur_dim-1; k++)
	for(j=0; j<=cur_dim-1; j++)
	for(i=0; i<=cur_dim-1; i++)
	{
		if(i==0 || i==cur_dim-1 || j==0 || j==cur_dim-1 || k==0 || k==cur_dim-1)
			A[i][j][k]= 0.;
		else A[i][j][k]= ( 4. + i + j + k) ;
	}
} 

void delete_data(){
	for (int l=0;l<cur_dim;l++){
		for (int q=0;q<cur_dim;q++){
			free(A[l][q]);
		}
		free(A[l]);
	}
}

void relax()
{


	for(i=1; i<=cur_dim-2; i++)
	#pragma omp parallel for private(k,j) num_threads(c_num_threads)
	for(k=1; k<=cur_dim-2; k++)
	for(j=1; j<=cur_dim-2; j++)
	{
		A[i][j][k] = (A[i-1][j][k]+A[i+1][j][k])/2.;
	}
	
	

	for(j=1; j<=cur_dim-2; j++)
	#pragma omp parallel for private(k,i) num_threads(c_num_threads) 
	for(k=1; k<=cur_dim-2; k++)
	for(i=1; i<=cur_dim-2; i++)
	{
		A[i][j][k] =(A[i][j-1][k]+A[i][j+1][k])/2.;
	}



	
	for(k=1; k<=cur_dim-2; k++)	
	#pragma omp parallel for private(i,j) num_threads(c_num_threads)
	for(j=1; j<=cur_dim-2; j++)
	for(i=1; i<=cur_dim-2; i++)
	{	
		double e;
		e=A[i][j][k];
		A[i][j][k] = (A[i][j][k-1]+A[i][j][k+1])/2.;
		eps=Max(eps,fabs(e-A[i][j][k]));
	}

}


double verify()
{
	double s;

	s=0.;
	
	#pragma omp parallel for private(k,j,i) reduction(+:s) num_threads(c_num_threads) collapse(3)
	for(k=0; k<=cur_dim-1; k++)
	for(j=0; j<=cur_dim-1; j++)
	for(i=0; i<=cur_dim-1; i++)
	{
		s=s+A[i][j][k]*(i+1)*(j+1)*(k+1)/(cur_dim*cur_dim*cur_dim);
	}

	return s;

}
