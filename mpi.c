#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>

// вариант 40
// автор: Грищенко Алексей
// параллелизация с помощью mpi 


#define  Max(a,b) ((a)>(b)?(a):(b))

#define N (2*2*2*2*2*2*2*2+2)

double   maxeps = 0.1e-7;
int itmax = 20;
int i,j,k;




double A[N][N][N];

void relax();
void init();
double verify(); 
double get_max_time();

double get_max_time(double* timearr,int proc_num){

}


int main(int argc, char **argv)
{
	int it;
	double result; 
    double eps;
    double start,end,maxtime;
    int proc_count;
    int proc_id
    double *times = NULL;

    MPI_INIT(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD,&proc_count);
    MPI_Comm_rank(MPI_COMM_WORLD,&proc_id);

    start = MPI_Wtime()

    if (proc_id == 0){
        times = (double*) malloc(time)
    }

	init();


	for(it=1; it<=itmax; it++)
	{
		eps = 0.;
		relax();
		if (eps < maxeps) break;
	}

	result = verify();

    end = MPI_Wtime();
	times[proc_id] = end-start;

	MPI_Barrier();

	if (proc_id == 0)
	{
		maxtime = get_max_time(times,proc_count);		
		printf("%d %d %f %f\n",N,proc_count,maxtime,result);		
	}

    MPI_FINALIZE();
	return 0;
}


void init()
{ 
	for(i=0; i<=N-1; i++)
	for(j=0; j<=N-1; j++)
	for(k=0; k<=N-1; k++)

	{
		if(i==0 || i==N-1 || j==0 || j==N-1 || k==0 || k==N-1)
			A[i][j][k]= 0.;
		else A[i][j][k]= ( 4. + i + j + k) ;
	}
} 

void relax()
{

	
		for(i=1; i<=N-2; i++)
		for(j=1; j<=N-2; j++)
		for(k=1; k<=N-2; k++)
		{
			A[i][j][k] = (A[i-1][j][k]+A[i+1][j][k])/2.;
		}
		

		for(j=1; j<=N-2; j++)
		for(i=1; i<=N-2; i++)
		for(k=1; k<=N-2; k++)
		
		{
			A[i][j][k] =(A[i][j-1][k]+A[i][j+1][k])/2.;
		}

		for(k=1; k<=N-2; k++)	
		for(i=1; i<=N-2; i++)
		for(j=1; j<=N-2; j++)
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
	
	for(i=0; i<=N-1; i++)
	for(j=0; j<=N-1; j++)
	for(k=0; k<=N-1; k++)
	{
		s=s+A[i][j][k]*(i+1)*(j+1)*(k+1)/(N*N*N);
	}

	return s;

}


