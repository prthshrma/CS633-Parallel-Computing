#include <stdio.h>
#include <string.h>
#include "mpi.h"
#include<math.h>
#include<stdlib.h>

double multipleMPISends(int N,int size,int rank);
double mpiPack(int N,int size,int rank);
double mpiDerived(int N,int size,int rank);

int rank,row,col,N,stepSize,p,size,commRank_L,commRank_R,commRank_T,commRank_B;

int main( int argc, char *argv[])
{

 double time,sTime,eTime,maxTime=-1;
 N=sqrt(atoi(argv[1]));
 stepSize=atoi(argv[2]);
 MPI_Init(&argc,&argv);
 MPI_Comm_rank( MPI_COMM_WORLD, &rank);
 MPI_Comm_size(MPI_COMM_WORLD, &size);
 
 FILE *fptr;
 fptr = fopen("data.txt","a+");     //file for storing time 
 if(fptr == NULL)
 {
    printf("Error!");   
    exit(1);             
 }
if(fptr!=NULL)
{
 fseek (fptr, 0, SEEK_END);
 int s= ftell(fptr);

    if (0 == s&&rank==0) {
	 fprintf(fptr,"function\tsqrt(N)\ttime\tP\n");      //Name of columns in data.txt file
    }
}

if(rank==0)
{
	printf("P = %d, sqrt(N) = %d\n",size,N);
}

 time=multipleMPISends(N,size,rank);  	//Computation and Communication time Using multiple sends
 MPI_Reduce (&time, &maxTime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
 if(maxTime>0)
 {
	//printf("Send Receive Time is for N:%d is %lf\n",N,maxTime);
	fprintf(fptr,"multipleMPISends\t%d\t%lf\t%d\n",N,maxTime,size);
 }

 maxTime=-1;
 time=mpiPack(N,size,rank); 	//Computation and Communication time using mpiPack
 MPI_Reduce (&time, &maxTime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
 if(maxTime>0)
 {    
	//printf("MPI_Packed Time is %lf\n",maxTime);
	fprintf(fptr,"mpiPack\t%d\t%lf\t%d\n",N,maxTime,size);
 }

 maxTime=-1;
 time=mpiDerived(N,size,rank);    //Computation and Communication time using mpiDerived
 MPI_Reduce (&time, &maxTime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
 if(maxTime>0)
 {
	//printf("Derived vector Time is %lf\n",maxTime);
	fprintf(fptr,"mpiDerived\t%d\t%lf\t%d\n",N,maxTime,size);
 }

 MPI_Finalize();
 fclose(fptr);
return 0;
}//end of main()



//start of function definition multipleMPISends()
double multipleMPISends(int N,int size,int rank)
{

MPI_Status status;
double** buf = (double**)malloc(N*sizeof(double*));
double** t = (double**)malloc(N*sizeof(double*));
for(int i=0;i<N;i++)
{
	buf[i] = (double*)malloc(N*sizeof(double));
	t[i] = (double*)malloc(N*sizeof(double)); 	
}

for(int i=0;i<N;i++)
{
	for(int j=0;j<N;j++)
	{
		buf[i][j]=rand() % size ;     //random initialization of data
	}
}
p = sqrt((double)size);
row=rank/p;	//row of processGrid
col=rank%p;	//column of processGrid

double sTime,eTime,time;
sTime=MPI_Wtime();//timer started

for(int steps=0;steps<stepSize;steps++)
{
MPI_Barrier(MPI_COMM_WORLD);
if(row==0)
{
	if(col==0)
	{
	       commRank_R = 1;
	       commRank_B = p;
	       double recvRight[N],recvDown[N];
	       for(int r=0;r<N;r++)
	       {
	       	MPI_Send(&buf[r][N-1], 1, MPI_DOUBLE,commRank_R , rank, MPI_COMM_WORLD);
	       	MPI_Recv(&recvRight[r], 1, MPI_DOUBLE, commRank_R, commRank_R, MPI_COMM_WORLD, &status);
	       }
	       for(int c=0;c<N;c++)
	       {
	       	 MPI_Send(&buf[N-1][c], 1, MPI_DOUBLE,commRank_B , rank, MPI_COMM_WORLD);
	       	 MPI_Recv(&recvDown[c], 1, MPI_DOUBLE, commRank_B, commRank_B, MPI_COMM_WORLD, &status);
	       }

	       //local computation
	       for(int i=0;i<N;i++)
	       {
	       	for(int j=0;j<N;j++)
	       	{
	       		if(i==0)
	       		{
	       			if(j==0)
	       			{
	       				t[i][j]=(buf[i][j+1]+buf[i+1][j])/2;
	       			}
	       			else if(j==N-1)
	       			{
	       				t[i][j]=(buf[i][j-1]+buf[i+1][j]+recvRight[i])/3;
	       			}
	       			else
	       			{
	       				t[i][j]=(buf[i][j+1]+buf[i][j-1]+buf[i+1][j])/3;
	       			}
	       		}
	       		else if(i==N-1)
	       		{
	       			if(j==0)
	       			{
	       				t[i][j]=(buf[i][j+1]+buf[i-1][j]+recvDown[j])/3;
	       			}
	       			else if(j==N-1)
	       			{
	       				t[i][j]=(buf[i][j-1]+buf[i-1][j]+recvRight[i]+recvDown[j])/4;
	       			}
	       			else
	       			{
	       				t[i][j]=(buf[i][j+1]+buf[i][j-1]+buf[i-1][j]+recvDown[j])/4;
	       			}
	       		}
	       		else if(j==0)
	       		{
	       			t[i][j]=(buf[i-1][j]+buf[i+1][j]+buf[i][j+1])/3;
	       		}
	       		else if(j==N-1)
	       		{
	       			t[i][j]=(buf[i-1][j]+buf[i+1][j]+buf[i][j-1]+recvRight[i])/4;
	       		}
	       		else
	       		{
	       			t[i][j]=(buf[i-1][j]+buf[i+1][j]+buf[i][j-1]+buf[i][j+1])/4;
	       		}
	       	}
	       }
	}
	else if(col==p-1)
	{
	       commRank_L = rank-1;
	       commRank_B = rank+p;
	       double recvLeft[N],recvDown[N];
	       for(int r=0;r<N;r++)
	       {
	       	 MPI_Send(&buf[r][0], 1, MPI_DOUBLE,commRank_L , rank, MPI_COMM_WORLD);
	       	 MPI_Recv(&recvLeft[r], 1, MPI_DOUBLE, commRank_L, commRank_L, MPI_COMM_WORLD, &status);
	       }
	       for(int c=0;c<N;c++)
	       {
	       	 MPI_Send(&buf[N-1][c], 1, MPI_DOUBLE,commRank_B , rank, MPI_COMM_WORLD);
	       	 MPI_Recv(&recvDown[c], 1, MPI_DOUBLE, commRank_B, commRank_B, MPI_COMM_WORLD, &status);
	       }
	      
	       //local computation
	       for(int i=0;i<N;i++)
	       {
	       	for(int j=0;j<N;j++)
	       	{
	       		if(i==0)
	       		{
	       			if(j==0)
	       			{
	       				t[i][j]=(buf[i][j+1]+buf[i+1][j]+recvLeft[i])/3;
	       			}
	       			else if(j==N-1)
	       			{
	       				t[i][j]=(buf[i][j-1]+buf[i+1][j])/2;
	       			}
	       			else
	       			{
	       				t[i][j]=(buf[i][j+1]+buf[i][j-1]+buf[i+1][j])/3;
	       			}
	       		}
	       		else if(i==N-1)
	       		{
	       			if(j==0)
	       			{
	       				t[i][j]=(buf[i][j+1]+buf[i-1][j]+recvDown[j]+recvLeft[i])/4;
	       			}
	       			else if(j==N-1)
	       			{
	       				t[i][j]=(buf[i][j-1]+buf[i-1][j]+recvDown[j])/3;
	       			}
	       			else
	       			{
	       				t[i][j]=(buf[i][j+1]+buf[i][j-1]+buf[i-1][j]+recvDown[j])/4;
	       			}
	       		}
	       		else if(j==0)
	       		{
	       			t[i][j]=(buf[i-1][j]+buf[i+1][j]+buf[i][j+1]+recvLeft[i])/4;
	       		}
	       		else if(j==N-1)
	       		{
	       			t[i][j]=(buf[i-1][j]+buf[i+1][j]+buf[i][j-1])/3;
	       		}
	       		else
	       		{
	       			t[i][j]=(buf[i-1][j]+buf[i+1][j]+buf[i][j-1]+buf[i][j+1])/4;
	       		}
	       	}
	       }
	}
	else
        {
               commRank_L = rank-1;
               commRank_R = rank+1;
               commRank_B = rank+p;
               double recvLeft[N],recvDown[N],recvRight[N];
               for(int r=0;r<N;r++)
	       {
	       	 MPI_Send(&buf[r][N-1], 1, MPI_DOUBLE,commRank_R , rank, MPI_COMM_WORLD);
	       	 MPI_Recv(&recvRight[r], 1, MPI_DOUBLE, commRank_R, commRank_R, MPI_COMM_WORLD, &status);
	       }
	       for(int r=0;r<N;r++)
	       {
	       	 MPI_Send(&buf[r][0], 1, MPI_DOUBLE,commRank_L , rank, MPI_COMM_WORLD);
	       	 MPI_Recv(&recvLeft[r], 1, MPI_DOUBLE, commRank_L, commRank_L, MPI_COMM_WORLD, &status);
	       }
	       for(int c=0;c<N;c++)
	       {
	       	 MPI_Send(&buf[N-1][c], 1, MPI_DOUBLE,commRank_B , rank, MPI_COMM_WORLD);
	       	 MPI_Recv(&recvDown[c], 1, MPI_DOUBLE, commRank_B, commRank_B, MPI_COMM_WORLD, &status);
	       }

         	//local computation
		for(int i=0;i<N;i++)
	       {
	       	for(int j=0;j<N;j++)
	       	{
	       		if(i==0)
	       		{
	       			if(j==0)
	       			{
	       				t[i][j]=(buf[i][j+1]+buf[i+1][j]+recvLeft[i])/3;
	       			}
	       			else if(j==N-1)
	       			{
	       				t[i][j]=(buf[i][j-1]+buf[i+1][j]+recvRight[i])/3;
	       			}
	       			else
	       			{
	       				t[i][j]=(buf[i][j+1]+buf[i][j-1]+buf[i+1][j])/3;
	       			}
	       		}
	       		else if(i==N-1)
	       		{
	       			if(j==0)
	       			{
	       				t[i][j]=(buf[i][j+1]+buf[i-1][j]+recvDown[j]+recvLeft[i])/4;
	       			}
	       			else if(j==N-1)
	       			{
	       				t[i][j]=(buf[i][j-1]+buf[i-1][j]+recvDown[j]+recvRight[i])/4;
	       			}
	       			else
	       			{
	       				t[i][j]=(buf[i][j+1]+buf[i][j-1]+buf[i-1][j]+recvDown[j])/4;
	       			}
	       		}
	       		else if(j==0)
	       		{
	       			t[i][j]=(buf[i-1][j]+buf[i+1][j]+buf[i][j+1]+recvLeft[i])/4;
	       		}
	       		else if(j==N-1)
	       		{
	       			t[i][j]=(buf[i-1][j]+buf[i+1][j]+buf[i][j-1]+recvRight[i])/4;
	       		}
	       		else
	       		{
	       			t[i][j]=(buf[i-1][j]+buf[i+1][j]+buf[i][j-1]+buf[i][j+1])/4;
	       		}
	       	}
	       }       	
       	
                 
        }
} 
else if(row==p-1) 
{
        if(col==0)
        { 
               commRank_R = rank+1;
               commRank_T = rank-p;
	       double recvRight[N],recvTop[N];
	       for(int r=0;r<N;r++)
	       {
	       	 MPI_Send(&buf[r][N-1], 1, MPI_DOUBLE,commRank_R , rank, MPI_COMM_WORLD);
	       	 MPI_Recv(&recvRight[r], 1, MPI_DOUBLE, commRank_R, commRank_R, MPI_COMM_WORLD, &status);
	       }
	       for(int c=0;c<N;c++)
	       {
	       	 MPI_Send(&buf[0][c], 1, MPI_DOUBLE,commRank_T , rank, MPI_COMM_WORLD);
	       	 MPI_Recv(&recvTop[c], 1, MPI_DOUBLE, commRank_T, commRank_T, MPI_COMM_WORLD, &status);
	       }
		
		//local computation
	       for(int i=0;i<N;i++)
	       {
	       	for(int j=0;j<N;j++)
	       	{
	       		if(i==0)
	       		{
	       			if(j==0)
	       			{
	       				t[i][j]=(buf[i][j+1]+buf[i+1][j]+recvTop[j])/3;
	       			}
	       			else if(j==N-1)
	       			{
	       				t[i][j]=(buf[i][j-1]+buf[i+1][j]+recvRight[i]+recvTop[j])/4;
	       			}
	       			else
	       			{
	       				t[i][j]=(buf[i][j+1]+buf[i][j-1]+buf[i+1][j]+recvTop[j])/4;
	       			}
	       		}
	       		else if(i==N-1)
	       		{
	       			if(j==0)
	       			{
	       				t[i][j]=(buf[i][j+1]+buf[i-1][j])/2;
	       			}
	       			else if(j==N-1)
	       			{
	       				t[i][j]=(buf[i][j-1]+buf[i-1][j]+recvRight[i])/3;
	       			}
	       			else
	       			{
	       				t[i][j]=(buf[i][j+1]+buf[i][j-1]+buf[i-1][j])/3;
	       			}
	       		}
	       		else if(j==0)
	       		{
	       			t[i][j]=(buf[i-1][j]+buf[i+1][j]+buf[i][j+1])/3;
	       		}
	       		else if(j==N-1)
	       		{
	       			t[i][j]=(buf[i-1][j]+buf[i+1][j]+buf[i][j-1]+recvRight[i])/4;
	       		}
	       		else
	       		{
	       			t[i][j]=(buf[i-1][j]+buf[i+1][j]+buf[i][j-1]+buf[i][j+1])/4;
	       		}
	       	}
	       }
        }
        else if(col==p-1)
        { 
               commRank_L = rank-1;
               commRank_T = rank-p;
               double recvLeft[N],recvTop[N];
	       for(int r=0;r<N;r++)
	       {
	       	 MPI_Send(&buf[r][0], 1, MPI_DOUBLE,commRank_L , rank, MPI_COMM_WORLD);
	       	 MPI_Recv(&recvLeft[r], 1, MPI_DOUBLE, commRank_L, commRank_L, MPI_COMM_WORLD, &status);
	       }
	       for(int c=0;c<N;c++)
	       {
	       	 MPI_Send(&buf[0][c], 1, MPI_DOUBLE,commRank_T , rank, MPI_COMM_WORLD);
	       	 MPI_Recv(&recvTop[c], 1, MPI_DOUBLE, commRank_T, commRank_T, MPI_COMM_WORLD, &status);
	       }
	       
	       //local computation
	       for(int i=0;i<N;i++)
	       {
	       	for(int j=0;j<N;j++)
	       	{
	       		if(i==0)
	       		{
	       			if(j==0)
	       			{
	       				t[i][j]=(buf[i][j+1]+buf[i+1][j]+recvLeft[i]+recvTop[j])/4;
	       			}
	       			else if(j==N-1)
	       			{
	       				t[i][j]=(buf[i][j-1]+buf[i+1][j]+recvTop[j])/3;
	       			}
	       			else
	       			{
	       				t[i][j]=(buf[i][j+1]+buf[i][j-1]+buf[i+1][j]+recvTop[j])/4;
	       			}
	       		}
	       		else if(i==N-1)
	       		{
	       			if(j==0)
	       			{
	       				t[i][j]=(buf[i][j+1]+buf[i-1][j]+recvLeft[i])/3;
	       			}
	       			else if(j==N-1)
	       			{
	       				t[i][j]=(buf[i][j-1]+buf[i-1][j])/2;
	       			}
	       			else
	       			{
	       				t[i][j]=(buf[i][j+1]+buf[i][j-1]+buf[i-1][j])/3;
	       			}
	       		}
	       		else if(j==0)
	       		{
	       			t[i][j]=(buf[i-1][j]+buf[i+1][j]+buf[i][j+1]+recvLeft[i])/4;
	       		}
	       		else if(j==N-1)
	       		{
	       			t[i][j]=(buf[i-1][j]+buf[i+1][j]+buf[i][j-1])/3;
	       		}
	       		else
	       		{
	       			t[i][j]=(buf[i-1][j]+buf[i+1][j]+buf[i][j-1]+buf[i][j+1])/4;
	       		}
	       	}
	       }
        }
        else
        { 
               commRank_L = rank-1;
               commRank_R = rank+1;
               commRank_T = rank-p;
               double recvLeft[N],recvTop[N],recvRight[N];
               for(int r=0;r<N;r++)
	       {
	       	 MPI_Send(&buf[r][N-1], 1, MPI_DOUBLE,commRank_R , rank, MPI_COMM_WORLD);
	       	 MPI_Recv(&recvRight[r], 1, MPI_DOUBLE, commRank_R, commRank_R, MPI_COMM_WORLD, &status);
	       }
	       for(int r=0;r<N;r++)
	       {
	       	 MPI_Send(&buf[r][0], 1, MPI_DOUBLE,commRank_L , rank, MPI_COMM_WORLD);
	       	 MPI_Recv(&recvLeft[r], 1, MPI_DOUBLE, commRank_L, commRank_L, MPI_COMM_WORLD, &status);
	       }
	       for(int c=0;c<N;c++)
	       {
	       	 MPI_Send(&buf[0][c], 1, MPI_DOUBLE,commRank_T , rank, MPI_COMM_WORLD);
	       	 MPI_Recv(&recvTop[c], 1, MPI_DOUBLE, commRank_T, commRank_T, MPI_COMM_WORLD, &status);
	       }
       	
       	//local computation
       	for(int i=0;i<N;i++)
	       {
	       	for(int j=0;j<N;j++)
	       	{
	       		if(i==0)
	       		{
	       			if(j==0)
	       			{
	       				t[i][j]=(buf[i][j+1]+buf[i+1][j]+recvLeft[i]+recvTop[j])/4;
	       			}
	       			else if(j==N-1)
	       			{
	       				t[i][j]=(buf[i][j-1]+buf[i+1][j]+recvRight[i]+recvTop[j])/4;
	       			}
	       			else
	       			{
	       				t[i][j]=(buf[i][j+1]+buf[i][j-1]+buf[i+1][j]+recvTop[j])/4;
	       			}
	       		}
	       		else if(i==N-1)
	       		{
	       			if(j==0)
	       			{
	       				t[i][j]=(buf[i][j+1]+buf[i-1][j]+recvLeft[i])/3;
	       			}
	       			else if(j==N-1)
	       			{
	       				t[i][j]=(buf[i][j-1]+buf[i-1][j]+recvRight[i])/3;
	       			}
	       			else
	       			{
	       				t[i][j]=(buf[i][j+1]+buf[i][j-1]+buf[i-1][j])/3;
	       			}
	       		}
	       		else if(j==0)
	       		{
	       			t[i][j]=(buf[i-1][j]+buf[i+1][j]+buf[i][j+1]+recvLeft[i])/4;
	       		}
	       		else if(j==N-1)
	       		{
	       			t[i][j]=(buf[i-1][j]+buf[i+1][j]+buf[i][j-1]+recvRight[i])/4;
	       		}
	       		else
	       		{
	       			t[i][j]=(buf[i-1][j]+buf[i+1][j]+buf[i][j-1]+buf[i][j+1])/4;
	       		}
	       	}
	       }
       }
} 
else if(col==0)
{
		commRank_B = rank+p;
                commRank_R = rank+1;
                commRank_T = rank-p;
		double recvDown[N],recvTop[N],recvRight[N];
               for(int r=0;r<N;r++)
	       {
	       	 MPI_Send(&buf[r][N-1], 1, MPI_DOUBLE,commRank_R , rank, MPI_COMM_WORLD);
	       	 MPI_Recv(&recvRight[r], 1, MPI_DOUBLE, commRank_R, commRank_R, MPI_COMM_WORLD, &status);
	       }
	       for(int c=0;c<N;c++)
	       {
	       	 MPI_Send(&buf[N-1][c], 1, MPI_DOUBLE, commRank_B , rank, MPI_COMM_WORLD);
	       	 MPI_Recv(&recvDown[c], 1, MPI_DOUBLE, commRank_B, commRank_B, MPI_COMM_WORLD, &status);
	       }
	       for(int c=0;c<N;c++)
	       {
	       	 MPI_Send(&buf[0][c], 1, MPI_DOUBLE,commRank_T , rank, MPI_COMM_WORLD);
	       	 MPI_Recv(&recvTop[c], 1, MPI_DOUBLE, commRank_T, commRank_T, MPI_COMM_WORLD, &status);
	       }
	       
	       //local Computation
	       for(int i=0;i<N;i++)
	       {
	       	for(int j=0;j<N;j++)
	       	{
	       		if(i==0)
	       		{
	       			if(j==0)
	       			{
	       				t[i][j]=(buf[i][j+1]+buf[i+1][j]+recvTop[j])/3;
	       			}
	       			else if(j==N-1)
	       			{
	       				t[i][j]=(buf[i][j-1]+buf[i+1][j]+recvTop[j])/3;
	       			}
	       			else
	       			{
	       				t[i][j]=(buf[i][j+1]+buf[i][j-1]+buf[i+1][j]+recvTop[j])/4;
	       			}
	       		}
	       		else if(i==N-1)
	       		{
	       			if(j==0)
	       			{
	       				t[i][j]=(buf[i][j+1]+buf[i-1][j]+recvDown[j])/3;
	       			}
	       			else if(j==N-1)
	       			{
	       				t[i][j]=(buf[i][j-1]+buf[i-1][j]+recvRight[i]+recvDown[j])/4;
	       			}
	       			else
	       			{
	       				t[i][j]=(buf[i][j+1]+buf[i][j-1]+buf[i-1][j]+recvDown[j])/4;
	       			}
	       		}
	       		else if(j==0)
	       		{
	       			t[i][j]=(buf[i-1][j]+buf[i+1][j]+buf[i][j+1])/3;
	       		}
	       		else if(j==N-1)
	       		{
	       			t[i][j]=(buf[i-1][j]+buf[i+1][j]+buf[i][j-1]+recvRight[i])/4;
	       		}
	       		else
	       		{
	       			t[i][j]=(buf[i-1][j]+buf[i+1][j]+buf[i][j-1]+buf[i][j+1])/4;
	       		}
	       	}
	       }
}
else if(col==p-1)
{
		commRank_L = rank-1;
                commRank_B = rank+p;
                commRank_T = rank-p;
		double recvDown[N],recvTop[N],recvLeft[N];
               for(int r=0;r<N;r++)
	       {
	       	 MPI_Send(&buf[r][0], 1, MPI_DOUBLE, commRank_L , rank, MPI_COMM_WORLD);
	       	 MPI_Recv(&recvLeft[r], 1, MPI_DOUBLE, commRank_L, commRank_L, MPI_COMM_WORLD, &status);
	       }
	       for(int c=0;c<N;c++)
	       {
	       	 MPI_Send(&buf[N-1][c], 1, MPI_DOUBLE, commRank_B , rank, MPI_COMM_WORLD);
	       	 MPI_Recv(&recvDown[c], 1, MPI_DOUBLE, commRank_B, commRank_B, MPI_COMM_WORLD, &status);
	       }
	       for(int c=0;c<N;c++)
	       {
	       	 MPI_Send(&buf[0][c], 1, MPI_DOUBLE,commRank_T , rank, MPI_COMM_WORLD);
	       	 MPI_Recv(&recvTop[c], 1, MPI_DOUBLE, commRank_T, commRank_T, MPI_COMM_WORLD, &status);
	       }
	       
	       //local Computation
	       for(int i=0;i<N;i++)
	       {
	       	for(int j=0;j<N;j++)
	       	{
	       		if(i==0)
	       		{
	       			if(j==0)
	       			{
	       				t[i][j]=(buf[i][j+1]+buf[i+1][j]+recvTop[j]+recvLeft[i])/4;
	       			}
	       			else if(j==N-1)
	       			{
	       				t[i][j]=(buf[i][j-1]+buf[i+1][j]+recvTop[j])/3;
	       			}
	       			else
	       			{
	       				t[i][j]=(buf[i][j+1]+buf[i][j-1]+buf[i+1][j]+recvTop[j])/4;
	       			}
	       		}
	       		else if(i==N-1)
	       		{
	       			if(j==0)
	       			{
	       				t[i][j]=(buf[i][j+1]+buf[i-1][j]+recvDown[j]+recvLeft[i])/4;
	       			}
	       			else if(j==N-1)
	       			{
	       				t[i][j]=(buf[i][j-1]+buf[i-1][j]+recvDown[j])/3;
	       			}
	       			else
	       			{
	       				t[i][j]=(buf[i][j+1]+buf[i][j-1]+buf[i-1][j]+recvDown[j])/4;
	       			}
	       		}
	       		else if(j==0)
	       		{
	       			t[i][j]=(buf[i-1][j]+buf[i+1][j]+buf[i][j+1]+recvLeft[i])/4;
	       		}
	       		else if(j==N-1)
	       		{
	       			t[i][j]=(buf[i-1][j]+buf[i+1][j]+buf[i][j-1])/3;
	       		}
	       		else
	       		{
	       			t[i][j]=(buf[i-1][j]+buf[i+1][j]+buf[i][j-1]+buf[i][j+1])/4;
	       		}
	       	}
	       }
}
else
{
		commRank_L = rank-1;
                commRank_R = rank+1;
                commRank_T = rank-p;
		commRank_B = rank+p;

		double recvDown[N],recvTop[N],recvLeft[N],recvRight[N];
               for(int r=0;r<N;r++)
	       {
	       	 MPI_Send(&buf[r][N-1], 1, MPI_DOUBLE,commRank_R , rank, MPI_COMM_WORLD);
	       	 MPI_Recv(&recvRight[r], 1, MPI_DOUBLE, commRank_R, commRank_R, MPI_COMM_WORLD, &status);
	       }
               for(int r=0;r<N;r++)
	       {
	       	 MPI_Send(&buf[r][0], 1, MPI_DOUBLE, commRank_L , rank, MPI_COMM_WORLD);
	       	 MPI_Recv(&recvLeft[r], 1, MPI_DOUBLE, commRank_L, commRank_L, MPI_COMM_WORLD, &status);
	       }
	       for(int c=0;c<N;c++)
	       {
	       	 MPI_Send(&buf[N-1][c], 1, MPI_DOUBLE, commRank_B , rank, MPI_COMM_WORLD);
	       	 MPI_Recv(&recvDown[c], 1, MPI_DOUBLE, commRank_B, commRank_B, MPI_COMM_WORLD, &status);
	       }
	       for(int c=0;c<N;c++)
	       {
	       	 MPI_Send(&buf[0][c], 1, MPI_DOUBLE,commRank_T , rank, MPI_COMM_WORLD);
	       	 MPI_Recv(&recvTop[c], 1, MPI_DOUBLE, commRank_T, commRank_T, MPI_COMM_WORLD, &status);
	       }
	       
	       //local Computation
	       for(int i=0;i<N;i++)
	       {
	       	for(int j=0;j<N;j++)
	       	{
	       		if(i==0)
	       		{
	       			if(j==0)
	       			{
	       				t[i][j]=(buf[i][j+1]+buf[i+1][j]+recvTop[j]+recvLeft[i])/4;
	       			}
	       			else if(j==N-1)
	       			{
	       				t[i][j]=(buf[i][j-1]+buf[i+1][j]+recvTop[j]+recvRight[i])/4;
	       			}
	       			else
	       			{
	       				t[i][j]=(buf[i][j+1]+buf[i][j-1]+buf[i+1][j]+recvTop[j])/4;
	       			}
	       		}
	       		else if(i==N-1)
	       		{
	       			if(j==0)
	       			{
	       				t[i][j]=(buf[i][j+1]+buf[i-1][j]+recvDown[j]+recvLeft[i])/4;
	       			}
	       			else if(j==N-1)
	       			{
	       				t[i][j]=(buf[i][j-1]+buf[i-1][j]+recvDown[j]+recvRight[i])/4;
	       			}
	       			else
	       			{
	       				t[i][j]=(buf[i][j+1]+buf[i][j-1]+buf[i-1][j]+recvDown[j])/4;
	       			}
	       		}
	       		else if(j==0)
	       		{
	       			t[i][j]=(buf[i-1][j]+buf[i+1][j]+buf[i][j+1]+recvLeft[i])/4;
	       		}
	       		else if(j==N-1)
	       		{
	       			t[i][j]=(buf[i-1][j]+buf[i+1][j]+buf[i][j-1]+recvRight[i])/4;
	       		}
	       		else
	       		{
	       			t[i][j]=(buf[i-1][j]+buf[i+1][j]+buf[i][j-1]+buf[i][j+1])/4;
	       		}
	       	}
	       }
}

for(int i=0;i<N;i++)
{
	for(int j=0;j<N;j++)
		buf[i][j]=t[i][j];
}
} // end of loop steps

eTime=MPI_Wtime();   //timer stopped
time=eTime-sTime;
return time;

/*for(int i=0;i<N;i++)
{
	for(int j=0;j<N;j++)
		printf("%d  %lf ,   ",rank,buf[i][j]);
printf("\n");
}*/

}//end of multipleMPISends

double mpiPack(int N,int size,int rank)
{
MPI_Status status,status1,status2;
MPI_Request req1,req2,req;

double** buf = (double**)malloc(N*sizeof(double*));
double** t = (double**)malloc(N*sizeof(double*));
for(int i=0;i<N;i++)
{
	buf[i] = (double*)malloc(N*sizeof(double));
	t[i] = (double*)malloc(N*sizeof(double)); 	
}
for(int i=0;i<N;i++)
{
	for(int j=0;j<N;j++)
	{
		buf[i][j]=rand() % size ;   //random initialization of data
        }
}
p = sqrt((double)size);
row=rank/p;	//row of Process Grid 
col=rank%p;	//column of Process Grid

double sTime,eTime,time;
sTime=MPI_Wtime();     //timer started 

for(int steps=0;steps<stepSize;steps++)
{
MPI_Barrier(MPI_COMM_WORLD);
if(row==0)
{
	if(col==0)
	{
	       commRank_R = 1;
	       commRank_B = p;
	       int rpos=0,dpos=0;
	       double recvRight[N],recvDown[N];
	       double rightPack[N],downPack[N],rightUnpack[N],downUnpack[N];
	       for(int r=0;r<N;r++)
	       {
	         MPI_Pack(&buf[r][N-1],1,MPI_DOUBLE,rightPack,8*N,&rpos,MPI_COMM_WORLD);
	       	 MPI_Pack(&buf[N-1][r],1,MPI_DOUBLE,downPack,8*N,&dpos,MPI_COMM_WORLD);
	       }
		      
	       MPI_Send(rightPack, rpos, MPI_PACKED,commRank_R , rank, MPI_COMM_WORLD);
	       MPI_Send(downPack, dpos, MPI_PACKED,commRank_B , rank, MPI_COMM_WORLD);
      	       rpos=0; dpos=0;
	       MPI_Recv(rightUnpack, 8*N,MPI_PACKED, commRank_R, commRank_R, MPI_COMM_WORLD,&status);
	       MPI_Recv(downUnpack, 8*N, MPI_PACKED, commRank_B, commRank_B, MPI_COMM_WORLD,&status);

	       for(int c=0;c<N;c++)
	       {
	       	 MPI_Unpack(rightUnpack,8*N,&rpos,&recvRight[c],1,MPI_DOUBLE,MPI_COMM_WORLD);
	         MPI_Unpack(downUnpack,8*N,&dpos,&recvDown[c],1,MPI_DOUBLE,MPI_COMM_WORLD);
	       }
	       
	       //local computation
	       for(int i=0;i<N;i++)
	       {
	       	for(int j=0;j<N;j++)
	       	{
	       		if(i==0)
	       		{
	       			if(j==0)
	       			{
	       				t[i][j]=(buf[i][j+1]+buf[i+1][j])/2;
	       			}
	       			else if(j==N-1)
	       			{
	       				t[i][j]=(buf[i][j-1]+buf[i+1][j]+recvRight[i])/3;
	       			}
	       			else
	       			{
	       				t[i][j]=(buf[i][j+1]+buf[i][j-1]+buf[i+1][j])/3;
	       			}
	       		}
	       		else if(i==N-1)
	       		{
	       			if(j==0)
	       			{
	       				t[i][j]=(buf[i][j+1]+buf[i-1][j]+recvDown[j])/3;
	       			}
	       			else if(j==N-1)
	       			{
	       				t[i][j]=(buf[i][j-1]+buf[i-1][j]+recvRight[i]+recvDown[j])/4;
	       			}
	       			else
	       			{
	       				t[i][j]=(buf[i][j+1]+buf[i][j-1]+buf[i-1][j]+recvDown[j])/4;
	       			}
	       		}
	       		else if(j==0)
	       		{
	       			t[i][j]=(buf[i-1][j]+buf[i+1][j]+buf[i][j+1])/3;
	       		}
	       		else if(j==N-1)
	       		{
	       			t[i][j]=(buf[i-1][j]+buf[i+1][j]+buf[i][j-1]+recvRight[i])/4;
	       		}
	       		else
	       		{
	       			t[i][j]=(buf[i-1][j]+buf[i+1][j]+buf[i][j-1]+buf[i][j+1])/4;
	       		}
	       	}
	       }
	}
	else if(col==p-1)
	{
	       commRank_L = rank-1;
	       commRank_B = rank+p;
	       int lpos=0,dpos=0;
	       double recvLeft[N],recvDown[N];
	       double leftPack[N],downPack[N],leftUnpack[N],downUnpack[N];
	       for(int r=0;r<N;r++)
	       {
	       	 MPI_Pack(&buf[r][0],1,MPI_DOUBLE,leftPack,8*N,&lpos,MPI_COMM_WORLD);
	       	 MPI_Pack(&buf[N-1][r],1,MPI_DOUBLE,downPack,8*N,&dpos,MPI_COMM_WORLD);
	       }
	      
	       MPI_Send(leftPack, lpos, MPI_PACKED,commRank_L , rank, MPI_COMM_WORLD);
	       MPI_Send(downPack, dpos, MPI_PACKED,commRank_B , rank, MPI_COMM_WORLD);
               lpos=0; dpos=0;
	       MPI_Recv(leftUnpack, 8*N,MPI_PACKED, commRank_L, commRank_L, MPI_COMM_WORLD,&status);
	       MPI_Recv(downUnpack, 8*N, MPI_PACKED, commRank_B, commRank_B, MPI_COMM_WORLD,&status);
	
	       for(int c=0;c<N;c++)
	       {       	 
	       	MPI_Unpack(leftUnpack,8*N,&lpos,&recvLeft[c],1,MPI_DOUBLE,MPI_COMM_WORLD);
	         	MPI_Unpack(downUnpack,8*N,&dpos,&recvDown[c],1,MPI_DOUBLE,MPI_COMM_WORLD);
	       }
	       
	       //local computation
	       for(int i=0;i<N;i++)
	       {
	       	for(int j=0;j<N;j++)
	       	{
	       		if(i==0)
	       		{
	       			if(j==0)
	       			{
	       				t[i][j]=(buf[i][j+1]+buf[i+1][j]+recvLeft[i])/3;
	       			}
	       			else if(j==N-1)
	       			{
	       				t[i][j]=(buf[i][j-1]+buf[i+1][j])/2;
	       			}
	       			else
	       			{
	       				t[i][j]=(buf[i][j+1]+buf[i][j-1]+buf[i+1][j])/3;
	       			}
	       		}
	       		else if(i==N-1)
	       		{
	       			if(j==0)
	       			{
	       				t[i][j]=(buf[i][j+1]+buf[i-1][j]+recvDown[j]+recvLeft[i])/4;
	       			}
	       			else if(j==N-1)
	       			{
	       				t[i][j]=(buf[i][j-1]+buf[i-1][j]+recvDown[j])/3;
	       			}
	       			else
	       			{
	       				t[i][j]=(buf[i][j+1]+buf[i][j-1]+buf[i-1][j]+recvDown[j])/4;
	       			}
	       		}
	       		else if(j==0)
	       		{
	       			t[i][j]=(buf[i-1][j]+buf[i+1][j]+buf[i][j+1]+recvLeft[i])/4;
	       		}
	       		else if(j==N-1)
	       		{
	       			t[i][j]=(buf[i-1][j]+buf[i+1][j]+buf[i][j-1])/3;
	       		}
	       		else
	       		{
	       			t[i][j]=(buf[i-1][j]+buf[i+1][j]+buf[i][j-1]+buf[i][j+1])/4;
	       		}
	       	}
	       }
	}
	else
        {
               commRank_L = rank-1;
               commRank_R = rank+1;
               commRank_B = rank+p;
               int lpos=0,rpos=0,dpos=0,tpos=0;
               double recvLeft[N],recvDown[N],recvRight[N],recvTop[N];
               double rightPack[N],downPack[N],rightUnpack[N],downUnpack[N],leftPack[N],leftUnpack[N],topPack[N],topUnpack[N];
               
               for(int r=0;r<N;r++)
	       {
	       	 
	       	 MPI_Pack(&buf[r][N-1],1,MPI_DOUBLE,rightPack,8*N,&rpos,MPI_COMM_WORLD);
	       	 MPI_Pack(&buf[r][0],1,MPI_DOUBLE,leftPack,8*N,&lpos,MPI_COMM_WORLD);
	       	 MPI_Pack(&buf[N-1][r],1,MPI_DOUBLE,downPack,8*N,&dpos,MPI_COMM_WORLD);
	       	
	       }
	       MPI_Send(rightPack, rpos, MPI_PACKED,commRank_R , rank, MPI_COMM_WORLD);
	       MPI_Send(leftPack, lpos, MPI_PACKED,commRank_L , rank, MPI_COMM_WORLD);
	       MPI_Send(downPack, dpos, MPI_PACKED,commRank_B , rank, MPI_COMM_WORLD);
	       rpos=0; dpos=0; lpos=0; 
	       MPI_Recv(rightUnpack, 8*N,MPI_PACKED, commRank_R, commRank_R, MPI_COMM_WORLD,&status);
	       MPI_Recv(leftUnpack, 8*N, MPI_PACKED, commRank_L, commRank_L, MPI_COMM_WORLD,&status);
	       MPI_Recv(downUnpack, 8*N, MPI_PACKED, commRank_B, commRank_B, MPI_COMM_WORLD,&status);
 
	       for(int c=0;c<N;c++)
	       {
	       	 MPI_Unpack(rightUnpack,8*N,&rpos,&recvRight[c],1,MPI_DOUBLE,MPI_COMM_WORLD);
	       	 MPI_Unpack(leftUnpack,8*N,&lpos,&recvLeft[c],1,MPI_DOUBLE,MPI_COMM_WORLD);
	       	 MPI_Unpack(downUnpack,8*N,&dpos,&recvDown[c],1,MPI_DOUBLE,MPI_COMM_WORLD);
	       }
               
               //local computation
		for(int i=0;i<N;i++)
	       {
	       	for(int j=0;j<N;j++)
	       	{
	       		if(i==0)
	       		{
	       			if(j==0)
	       			{
	       				t[i][j]=(buf[i][j+1]+buf[i+1][j]+recvLeft[i])/3;
	       			}
	       			else if(j==N-1)
	       			{
	       				t[i][j]=(buf[i][j-1]+buf[i+1][j]+recvRight[i])/3;
	       			}
	       			else
	       			{
	       				t[i][j]=(buf[i][j+1]+buf[i][j-1]+buf[i+1][j])/3;
	       			}
	       		}
	       		else if(i==N-1)
	       		{
	       			if(j==0)
	       			{
	       				t[i][j]=(buf[i][j+1]+buf[i-1][j]+recvDown[j]+recvLeft[i])/4;
	       			}
	       			else if(j==N-1)
	       			{
	       				t[i][j]=(buf[i][j-1]+buf[i-1][j]+recvDown[j]+recvRight[i])/4;
	       			}
	       			else
	       			{
	       				t[i][j]=(buf[i][j+1]+buf[i][j-1]+buf[i-1][j]+recvDown[j])/4;
	       			}
	       		}
	       		else if(j==0)
	       		{
	       			t[i][j]=(buf[i-1][j]+buf[i+1][j]+buf[i][j+1]+recvLeft[i])/4;
	       		}
	       		else if(j==N-1)
	       		{
	       			t[i][j]=(buf[i-1][j]+buf[i+1][j]+buf[i][j-1]+recvRight[i])/4;
	       		}
	       		else
	       		{
	       			t[i][j]=(buf[i-1][j]+buf[i+1][j]+buf[i][j-1]+buf[i][j+1])/4;
	       		}
	       	}
	       }       	
       }
} 
else if(row==p-1) 
{
        if(col==0)
        { 
               commRank_R = rank+1;
               commRank_T = rank-p;
	       int tpos=0,rpos=0;
               double recvTop[N],recvRight[N];
               double rightPack[N],topPack[N],rightUnpack[N],topUnpack[N];
               
               for(int r=0;r<N;r++)
	       {
	       	 MPI_Pack(&buf[r][N-1],1,MPI_DOUBLE,rightPack,8*N,&rpos,MPI_COMM_WORLD);
	       	 MPI_Pack(&buf[0][r],1,MPI_DOUBLE,topPack,8*N,&tpos,MPI_COMM_WORLD);
	       }
	       MPI_Send(rightPack, rpos, MPI_PACKED,commRank_R , rank, MPI_COMM_WORLD);
	       MPI_Send(topPack, tpos, MPI_PACKED,commRank_T , rank, MPI_COMM_WORLD);
	       rpos=0; tpos=0; 
	       MPI_Recv(rightUnpack, 8*N,MPI_PACKED, commRank_R, commRank_R, MPI_COMM_WORLD,&status);
	       MPI_Recv(topUnpack, 8*N, MPI_PACKED, commRank_T, commRank_T, MPI_COMM_WORLD,&status);
	       for(int c=0;c<N;c++)
	       {
	       	 MPI_Unpack(rightUnpack,8*N,&rpos,&recvRight[c],1,MPI_DOUBLE,MPI_COMM_WORLD);
	       	 MPI_Unpack(topUnpack,8*N,&tpos,&recvTop[c],1,MPI_DOUBLE,MPI_COMM_WORLD);
	       }
	       
	       //local computation
	       for(int i=0;i<N;i++)
	       {
	       	for(int j=0;j<N;j++)
	       	{
	       		if(i==0)
	       		{
	       			if(j==0)
	       			{
	       				t[i][j]=(buf[i][j+1]+buf[i+1][j]+recvTop[j])/3;
	       			}
	       			else if(j==N-1)
	       			{
	       				t[i][j]=(buf[i][j-1]+buf[i+1][j]+recvRight[i]+recvTop[j])/4;
	       			}
	       			else
	       			{
	       				t[i][j]=(buf[i][j+1]+buf[i][j-1]+buf[i+1][j]+recvTop[j])/4;
	       			}
	       		}
	       		else if(i==N-1)
	       		{
	       			if(j==0)
	       			{
	       				t[i][j]=(buf[i][j+1]+buf[i-1][j])/2;
	       			}
	       			else if(j==N-1)
	       			{
	       				t[i][j]=(buf[i][j-1]+buf[i-1][j]+recvRight[i])/3;
	       			}
	       			else
	       			{
	       				t[i][j]=(buf[i][j+1]+buf[i][j-1]+buf[i-1][j])/3;
	       			}
	       		}
	       		else if(j==0)
	       		{
	       			t[i][j]=(buf[i-1][j]+buf[i+1][j]+buf[i][j+1])/3;
	       		}
	       		else if(j==N-1)
	       		{
	       			t[i][j]=(buf[i-1][j]+buf[i+1][j]+buf[i][j-1]+recvRight[i])/4;
	       		}
	       		else
	       		{
	       			t[i][j]=(buf[i-1][j]+buf[i+1][j]+buf[i][j-1]+buf[i][j+1])/4;
	       		}
	       	}
	       }
        }
        else if(col==p-1)
        { 
               commRank_L = rank-1;
               commRank_T = rank-p;
               int lpos=0,tpos=0;
               double recvLeft[N],recvTop[N];
               double leftPack[N],leftUnpack[N],topPack[N],topUnpack[N];
               
               for(int r=0;r<N;r++)
	       {
	       	MPI_Pack(&buf[r][0],1,MPI_DOUBLE,leftPack,8*N,&lpos,MPI_COMM_WORLD);
	       	MPI_Pack(&buf[0][r],1,MPI_DOUBLE,topPack,8*N,&tpos,MPI_COMM_WORLD);
	       }
	       MPI_Send(leftPack, lpos, MPI_PACKED,commRank_L , rank, MPI_COMM_WORLD);
	       MPI_Send(topPack, tpos, MPI_PACKED,commRank_T , rank, MPI_COMM_WORLD);
	       lpos=0; tpos=0; 
	       MPI_Recv(leftUnpack, 8*N,MPI_PACKED, commRank_L, commRank_L, MPI_COMM_WORLD,&status);
	       MPI_Recv(topUnpack, 8*N, MPI_PACKED, commRank_T, commRank_T, MPI_COMM_WORLD,&status); 
	       for(int c=0;c<N;c++)
	       {
	       	 MPI_Unpack(leftUnpack,8*N,&lpos,&recvLeft[c],1,MPI_DOUBLE,MPI_COMM_WORLD);
	       	 MPI_Unpack(topUnpack,8*N,&tpos,&recvTop[c],1,MPI_DOUBLE,MPI_COMM_WORLD);
	       }               
	       
	       //local computation
	       for(int i=0;i<N;i++)
	       {
	       	for(int j=0;j<N;j++)
	       	{
	       		if(i==0)
	       		{
	       			if(j==0)
	       			{
	       				t[i][j]=(buf[i][j+1]+buf[i+1][j]+recvLeft[i]+recvTop[j])/4;
	       			}
	       			else if(j==N-1)
	       			{
	       				t[i][j]=(buf[i][j-1]+buf[i+1][j]+recvTop[j])/3;
	       			}
	       			else
	       			{
	       				t[i][j]=(buf[i][j+1]+buf[i][j-1]+buf[i+1][j]+recvTop[j])/4;
	       			}
	       		}
	       		else if(i==N-1)
	       		{
	       			if(j==0)
	       			{
	       				t[i][j]=(buf[i][j+1]+buf[i-1][j]+recvLeft[i])/3;
	       			}
	       			else if(j==N-1)
	       			{
	       				t[i][j]=(buf[i][j-1]+buf[i-1][j])/2;
	       			}
	       			else
	       			{
	       				t[i][j]=(buf[i][j+1]+buf[i][j-1]+buf[i-1][j])/3;
	       			}
	       		}
	       		else if(j==0)
	       		{
	       			t[i][j]=(buf[i-1][j]+buf[i+1][j]+buf[i][j+1]+recvLeft[i])/4;
	       		}
	       		else if(j==N-1)
	       		{
	       			t[i][j]=(buf[i-1][j]+buf[i+1][j]+buf[i][j-1])/3;
	       		}
	       		else
	       		{
	       			t[i][j]=(buf[i-1][j]+buf[i+1][j]+buf[i][j-1]+buf[i][j+1])/4;
	       		}
	       	}
	       }
        }
        else
        { 
               commRank_L = rank-1;
               commRank_R = rank+1;
               commRank_T = rank-p;
               int lpos=0,rpos=0,tpos=0;
               double recvLeft[N],recvRight[N],recvTop[N];
               double rightPack[N],rightUnpack[N],leftPack[N],leftUnpack[N],topPack[N],topUnpack[N];
               
               for(int r=0;r<N;r++)
	       {
	       	 MPI_Pack(&buf[r][N-1],1,MPI_DOUBLE,rightPack,8*N,&rpos,MPI_COMM_WORLD);
	       	 MPI_Pack(&buf[r][0],1,MPI_DOUBLE,leftPack,8*N,&lpos,MPI_COMM_WORLD);
	       	 MPI_Pack(&buf[0][r],1,MPI_DOUBLE,topPack,8*N,&tpos,MPI_COMM_WORLD);
	       }
	       MPI_Send(rightPack, rpos, MPI_PACKED,commRank_R , rank, MPI_COMM_WORLD);
	       MPI_Send(leftPack, lpos, MPI_PACKED,commRank_L , rank, MPI_COMM_WORLD);
	       MPI_Send(topPack, tpos, MPI_PACKED,commRank_T , rank, MPI_COMM_WORLD);
	       rpos=0; tpos=0; lpos=0; 
	       MPI_Recv(rightUnpack, 8*N,MPI_PACKED, commRank_R, commRank_R, MPI_COMM_WORLD, &status);
	       MPI_Recv(leftUnpack, 8*N,MPI_PACKED, commRank_L, commRank_L, MPI_COMM_WORLD, &status);
	       MPI_Recv(topUnpack, 8*N, MPI_PACKED, commRank_T, commRank_T, MPI_COMM_WORLD, &status);
	        
	       for(int c=0;c<N;c++)
	       {	 
	       	 MPI_Unpack(rightUnpack,8*N,&rpos,&recvRight[c],1,MPI_DOUBLE,MPI_COMM_WORLD);
	       	 MPI_Unpack(leftUnpack,8*N,&lpos,&recvLeft[c],1,MPI_DOUBLE,MPI_COMM_WORLD);
	       	 MPI_Unpack(topUnpack,8*N,&tpos,&recvTop[c],1,MPI_DOUBLE,MPI_COMM_WORLD);
	       }
              
              //local computation
               for(int i=0;i<N;i++)
	       {
	       	for(int j=0;j<N;j++)
	       	{
	       		if(i==0)
	       		{
	       			if(j==0)
	       			{
	       				t[i][j]=(buf[i][j+1]+buf[i+1][j]+recvLeft[i]+recvTop[j])/4;
	       			}
	       			else if(j==N-1)
	       			{
	       				t[i][j]=(buf[i][j-1]+buf[i+1][j]+recvRight[i]+recvTop[j])/4;
	       			}
	       			else
	       			{
	       				t[i][j]=(buf[i][j+1]+buf[i][j-1]+buf[i+1][j]+recvTop[j])/4;
	       			}
	       		}
	       		else if(i==N-1)
	       		{
	       			if(j==0)
	       			{
	       				t[i][j]=(buf[i][j+1]+buf[i-1][j]+recvLeft[i])/3;
	       			}
	       			else if(j==N-1)
	       			{
	       				t[i][j]=(buf[i][j-1]+buf[i-1][j]+recvRight[i])/3;
	       			}
	       			else
	       			{
	       				t[i][j]=(buf[i][j+1]+buf[i][j-1]+buf[i-1][j])/3;
	       			}
	       		}
	       		else if(j==0)
	       		{
	       			t[i][j]=(buf[i-1][j]+buf[i+1][j]+buf[i][j+1]+recvLeft[i])/4;
	       		}
	       		else if(j==N-1)
	       		{
	       			t[i][j]=(buf[i-1][j]+buf[i+1][j]+buf[i][j-1]+recvRight[i])/4;
	       		}
	       		else
	       		{
	       			t[i][j]=(buf[i-1][j]+buf[i+1][j]+buf[i][j-1]+buf[i][j+1])/4;
	       		}
	       	}
	       }
       }
} 
else if(col==0)
{
		commRank_B = rank+p;
                commRank_R = rank+1;
                commRank_T = rank-p;
                
                int rpos=0,dpos=0,tpos=0;
               double recvDown[N],recvRight[N],recvTop[N];
               double rightPack[N],downPack[N],rightUnpack[N],downUnpack[N],topPack[N],topUnpack[N];
               
               for(int r=0;r<N;r++)
	       {
	       	 MPI_Pack(&buf[r][N-1],1,MPI_DOUBLE,rightPack,8*N,&rpos,MPI_COMM_WORLD);
	       	 MPI_Pack(&buf[N-1][r],1,MPI_DOUBLE,downPack,8*N,&dpos,MPI_COMM_WORLD);
	       	 MPI_Pack(&buf[0][r],1,MPI_DOUBLE,topPack,8*N,&tpos,MPI_COMM_WORLD);
	       }
	       MPI_Send(rightPack, rpos, MPI_PACKED,commRank_R , rank, MPI_COMM_WORLD);
	       MPI_Send(downPack, dpos, MPI_PACKED,commRank_B , rank, MPI_COMM_WORLD);
	       MPI_Send(topPack, tpos, MPI_PACKED,commRank_T , rank, MPI_COMM_WORLD);
	       rpos=0; dpos=0; tpos=0; 
	       MPI_Recv(rightUnpack, 8*N,MPI_PACKED, commRank_R, commRank_R, MPI_COMM_WORLD, &status);
	       MPI_Recv(downUnpack, 8*N, MPI_PACKED, commRank_B, commRank_B, MPI_COMM_WORLD, &status);
	       MPI_Recv(topUnpack, 8*N, MPI_PACKED, commRank_T, commRank_T, MPI_COMM_WORLD, &status);
	        
	       for(int c=0;c<N;c++)
	       {	 
	       	 MPI_Unpack(rightUnpack,8*N,&rpos,&recvRight[c],1,MPI_DOUBLE,MPI_COMM_WORLD);
	       	 MPI_Unpack(downUnpack,8*N,&dpos,&recvDown[c],1,MPI_DOUBLE,MPI_COMM_WORLD);
	       	 MPI_Unpack(topUnpack,8*N,&tpos,&recvTop[c],1,MPI_DOUBLE,MPI_COMM_WORLD);
	       }
	       
	       for(int i=0;i<N;i++)
	       {
	       	for(int j=0;j<N;j++)
	       	{
	       		if(i==0)
	       		{
	       			if(j==0)
	       			{
	       				t[i][j]=(buf[i][j+1]+buf[i+1][j]+recvTop[j])/3;
	       			}
	       			else if(j==N-1)
	       			{
	       				t[i][j]=(buf[i][j-1]+buf[i+1][j]+recvTop[j])/3;
	       			}
	       			else
	       			{
	       				t[i][j]=(buf[i][j+1]+buf[i][j-1]+buf[i+1][j]+recvTop[j])/4;
	       			}
	       		}
	       		else if(i==N-1)
	       		{
	       			if(j==0)
	       			{
	       				t[i][j]=(buf[i][j+1]+buf[i-1][j]+recvDown[j])/3;
	       			}
	       			else if(j==N-1)
	       			{
	       				t[i][j]=(buf[i][j-1]+buf[i-1][j]+recvRight[i]+recvDown[j])/4;
	       			}
	       			else
	       			{
	       				t[i][j]=(buf[i][j+1]+buf[i][j-1]+buf[i-1][j]+recvDown[j])/4;
	       			}
	       		}
	       		else if(j==0)
	       		{
	       			t[i][j]=(buf[i-1][j]+buf[i+1][j]+buf[i][j+1])/3;
	       		}
	       		else if(j==N-1)
	       		{
	       			t[i][j]=(buf[i-1][j]+buf[i+1][j]+buf[i][j-1]+recvRight[i])/4;
	       		}
	       		else
	       		{
	       			t[i][j]=(buf[i-1][j]+buf[i+1][j]+buf[i][j-1]+buf[i][j+1])/4;
	       		}
	       	}
	       }
}
else if(col==p-1)
{
		commRank_L = rank-1;
                commRank_B = rank+p;
                commRank_T = rank-p;
                int lpos=0,dpos=0,tpos=0;
               double recvLeft[N],recvDown[N],recvTop[N];
               double downPack[N],downUnpack[N],leftPack[N],leftUnpack[N],topPack[N],topUnpack[N];
               
               for(int r=0;r<N;r++)
	       {
	       	 MPI_Pack(&buf[r][0],1,MPI_DOUBLE,leftPack,8*N,&lpos,MPI_COMM_WORLD);
	       	 MPI_Pack(&buf[N-1][r],1,MPI_DOUBLE,downPack,8*N,&dpos,MPI_COMM_WORLD);
	       	 MPI_Pack(&buf[0][r],1,MPI_DOUBLE,topPack,8*N,&tpos,MPI_COMM_WORLD);
	       }
	       MPI_Send(leftPack, lpos, MPI_PACKED,commRank_L , rank, MPI_COMM_WORLD);
	       MPI_Send(downPack, dpos, MPI_PACKED,commRank_B , rank, MPI_COMM_WORLD);
	       MPI_Send(topPack, tpos, MPI_PACKED,commRank_T , rank, MPI_COMM_WORLD);
	       lpos=0; dpos=0; tpos=0; 
	       MPI_Recv(leftUnpack, 8*N,MPI_PACKED, commRank_L, commRank_L, MPI_COMM_WORLD, &status);
	       MPI_Recv(downUnpack, 8*N, MPI_PACKED, commRank_B, commRank_B, MPI_COMM_WORLD, &status);
	       MPI_Recv(topUnpack, 8*N, MPI_PACKED, commRank_T, commRank_T, MPI_COMM_WORLD, &status);
	        
	       for(int c=0;c<N;c++)
	       {	 
	       	 MPI_Unpack(leftUnpack,8*N,&lpos,&recvLeft[c],1,MPI_DOUBLE,MPI_COMM_WORLD);
	       	 MPI_Unpack(downUnpack,8*N,&dpos,&recvDown[c],1,MPI_DOUBLE,MPI_COMM_WORLD);
	       	 MPI_Unpack(topUnpack,8*N,&tpos,&recvTop[c],1,MPI_DOUBLE,MPI_COMM_WORLD);
	       }
               
               //local Computation       
	       for(int i=0;i<N;i++)
	       {
	       	for(int j=0;j<N;j++)
	       	{
	       		if(i==0)
	       		{
	       			if(j==0)
	       			{
	       				t[i][j]=(buf[i][j+1]+buf[i+1][j]+recvTop[j]+recvLeft[i])/4;
	       			}
	       			else if(j==N-1)
	       			{
	       				t[i][j]=(buf[i][j-1]+buf[i+1][j]+recvTop[j])/3;
	       			}
	       			else
	       			{
	       				t[i][j]=(buf[i][j+1]+buf[i][j-1]+buf[i+1][j]+recvTop[j])/4;
	       			}
	       		}
	       		else if(i==N-1)
	       		{
	       			if(j==0)
	       			{
	       				t[i][j]=(buf[i][j+1]+buf[i-1][j]+recvDown[j]+recvLeft[i])/4;
	       			}
	       			else if(j==N-1)
	       			{
	       				t[i][j]=(buf[i][j-1]+buf[i-1][j]+recvDown[j])/3;
	       			}
	       			else
	       			{
	       				t[i][j]=(buf[i][j+1]+buf[i][j-1]+buf[i-1][j]+recvDown[j])/4;
	       			}
	       		}
	       		else if(j==0)
	       		{
	       			t[i][j]=(buf[i-1][j]+buf[i+1][j]+buf[i][j+1]+recvLeft[i])/4;
	       		}
	       		else if(j==N-1)
	       		{
	       			t[i][j]=(buf[i-1][j]+buf[i+1][j]+buf[i][j-1])/3;
	       		}
	       		else
	       		{
	       			t[i][j]=(buf[i-1][j]+buf[i+1][j]+buf[i][j-1]+buf[i][j+1])/4;
	       		}
	       	}
	       }
}
else
{
		commRank_L = rank-1;
                commRank_R = rank+1;
                commRank_T = rank-p;
		commRank_B = rank+p;
		int lpos=0,rpos=0,dpos=0,tpos=0;
               double recvLeft[N],recvDown[N],recvRight[N],recvTop[N];
               double rightPack[N],downPack[N],rightUnpack[N],downUnpack[N],leftPack[N],leftUnpack[N],topPack[N],topUnpack[N];
               
               for(int r=0;r<N;r++)
	       {
	       	 MPI_Pack(&buf[r][N-1],1,MPI_DOUBLE,rightPack,8*N,&rpos,MPI_COMM_WORLD);
	       	 MPI_Pack(&buf[r][0],1,MPI_DOUBLE,leftPack,8*N,&lpos,MPI_COMM_WORLD);
	       	 MPI_Pack(&buf[N-1][r],1,MPI_DOUBLE,downPack,8*N,&dpos,MPI_COMM_WORLD);
	       	 MPI_Pack(&buf[0][r],1,MPI_DOUBLE,topPack,8*N,&tpos,MPI_COMM_WORLD);
	       }
	       MPI_Send(rightPack, rpos, MPI_PACKED,commRank_R , rank, MPI_COMM_WORLD);
	       MPI_Send(leftPack, lpos, MPI_PACKED,commRank_L , rank, MPI_COMM_WORLD);
	       MPI_Send(downPack, dpos, MPI_PACKED,commRank_B , rank, MPI_COMM_WORLD);
	       MPI_Send(topPack, tpos, MPI_PACKED,commRank_T , rank, MPI_COMM_WORLD);
	       rpos=0; dpos=0; lpos=0; tpos=0;  
	       MPI_Recv(rightUnpack, 8*N,MPI_PACKED, commRank_R, commRank_R, MPI_COMM_WORLD, &status);
	       MPI_Recv(leftUnpack, 8*N,MPI_PACKED, commRank_L, commRank_L, MPI_COMM_WORLD, &status);
	       MPI_Recv(downUnpack, 8*N, MPI_PACKED, commRank_B, commRank_B, MPI_COMM_WORLD, &status);
	       MPI_Recv(topUnpack, 8*N, MPI_PACKED, commRank_T, commRank_T, MPI_COMM_WORLD, &status);
	        
	       for(int c=0;c<N;c++)
	       { 
	       	 MPI_Unpack(rightUnpack,8*N,&rpos,&recvRight[c],1,MPI_DOUBLE,MPI_COMM_WORLD);
	       	 MPI_Unpack(leftUnpack,8*N,&lpos,&recvLeft[c],1,MPI_DOUBLE,MPI_COMM_WORLD);
	       	 MPI_Unpack(downUnpack,8*N,&dpos,&recvDown[c],1,MPI_DOUBLE,MPI_COMM_WORLD);
	       	 MPI_Unpack(topUnpack,8*N,&tpos,&recvTop[c],1,MPI_DOUBLE,MPI_COMM_WORLD);
	       }
		
		//local Computation               
	       for(int i=0;i<N;i++)
	       {
	       	for(int j=0;j<N;j++)
	       	{
	       		if(i==0)
	       		{
	       			if(j==0)
	       			{
	       				t[i][j]=(buf[i][j+1]+buf[i+1][j]+recvTop[j]+recvLeft[i])/4;
	       			}
	       			else if(j==N-1)
	       			{
	       				t[i][j]=(buf[i][j-1]+buf[i+1][j]+recvTop[j]+recvRight[i])/4;
	       			}
	       			else
	       			{
	       				t[i][j]=(buf[i][j+1]+buf[i][j-1]+buf[i+1][j]+recvTop[j])/4;
	       			}
	       		}
	       		else if(i==N-1)
	       		{
	       			if(j==0)
	       			{
	       				t[i][j]=(buf[i][j+1]+buf[i-1][j]+recvDown[j]+recvLeft[i])/4;
	       			}
	       			else if(j==N-1)
	       			{
	       				t[i][j]=(buf[i][j-1]+buf[i-1][j]+recvDown[j]+recvRight[i])/4;
	       			}
	       			else
	       			{
	       				t[i][j]=(buf[i][j+1]+buf[i][j-1]+buf[i-1][j]+recvDown[j])/4;
	       			}
	       		}
	       		else if(j==0)
	       		{
	       			t[i][j]=(buf[i-1][j]+buf[i+1][j]+buf[i][j+1]+recvLeft[i])/4;
	       		}
	       		else if(j==N-1)
	       		{
	       			t[i][j]=(buf[i-1][j]+buf[i+1][j]+buf[i][j-1]+recvRight[i])/4;
	       		}
	       		else
	       		{
	       			t[i][j]=(buf[i-1][j]+buf[i+1][j]+buf[i][j-1]+buf[i][j+1])/4;
	       		}
	       	}
	       }
}

for(int i=0;i<N;i++)
{
	for(int j=0;j<N;j++)
		buf[i][j]=t[i][j];
}
}  // end of loop steps

eTime=MPI_Wtime();     //timer stopped
time=eTime-sTime;
return time;
/*for(int i=0;i<N;i++)
{
	for(int j=0;j<N;j++)
		printf("%d  %lf ,   ",rank,buf[i][j]);
printf("\n");
}*/

}//end of mpiPack

//Start of mpiDerived Function
double mpiDerived(int N,int size,int rank)
{
MPI_Status status;
MPI_Datatype newvtypeRow,newvtypeCol;
int count=N,blocklen=1,strideR=N,strideC=1; 

MPI_Type_vector (count,blocklen, strideR, MPI_DOUBLE, &newvtypeRow);
MPI_Type_commit (&newvtypeRow);
MPI_Type_vector (count,blocklen, strideC, MPI_DOUBLE, &newvtypeCol);
MPI_Type_commit (&newvtypeCol);

double** t = (double**)malloc(N*sizeof(double*));

for (int i=0; i<N; i++) 
         t[i] = (double *)malloc(N * sizeof(double)); 

int len=sizeof(double *)*N+sizeof(double)*N*N;		//2D array Contiguous memory location
double *ptr,**buf;
buf = (double **)malloc(len);
ptr=(double*)(buf+N);
for(int i=0;i<N;i++)
{
buf[i]=(ptr+N*i);
} 

for(int i=0;i<N;i++)
{
	for(int j=0;j<N;j++)
	{
		buf[i][j]=rand() % size ;
        }
}
p = sqrt((double)size);
row=rank/p;		//row of Process Grid 
col=rank%p;		//column of Process Grid

double sTime,eTime,time;
sTime=MPI_Wtime();//timer started 

for(int steps=0;steps<stepSize;steps++)
{
MPI_Barrier(MPI_COMM_WORLD);
if(row==0)
{
	if(col==0)
	{
	       commRank_R = 1;
	       commRank_B = p;
	       double recvRight[N],recvDown[N];
	       
       	MPI_Send(&buf[0][N-1], 1, newvtypeRow,commRank_R , rank, MPI_COMM_WORLD);
       	MPI_Send(&buf[N-1][0], 1,newvtypeCol,commRank_B , rank, MPI_COMM_WORLD);
		MPI_Recv(recvRight, N, MPI_DOUBLE, commRank_R, commRank_R, MPI_COMM_WORLD, &status);
       	MPI_Recv(recvDown, N, MPI_DOUBLE, commRank_B, commRank_B, MPI_COMM_WORLD, &status);

 	       //local Computation
 	       for(int i=0;i<N;i++)
	       {
	       	for(int j=0;j<N;j++)
	       	{
	       		if(i==0)
	       		{
	       			if(j==0)
	       			{
	       				t[i][j]=(buf[i][j+1]+buf[i+1][j])/2;
	       			}
	       			else if(j==N-1)
	       			{
	       				t[i][j]=(buf[i][j-1]+buf[i+1][j]+recvRight[i])/3;
	       			}
	       			else
	       			{
	       				t[i][j]=(buf[i][j+1]+buf[i][j-1]+buf[i+1][j])/3;
	       			}
	       		}
	       		else if(i==N-1)
	       		{
	       			if(j==0)
	       			{
	       				t[i][j]=(buf[i][j+1]+buf[i-1][j]+recvDown[j])/3;
	       			}
	       			else if(j==N-1)
	       			{
	       				t[i][j]=(buf[i][j-1]+buf[i-1][j]+recvRight[i]+recvDown[j])/4;
	       			}
	       			else
	       			{
	       				t[i][j]=(buf[i][j+1]+buf[i][j-1]+buf[i-1][j]+recvDown[j])/4;
	       			}
	       		}
	       		else if(j==0)
	       		{
	       			t[i][j]=(buf[i-1][j]+buf[i+1][j]+buf[i][j+1])/3;
	       		}
	       		else if(j==N-1)
	       		{
	       			t[i][j]=(buf[i-1][j]+buf[i+1][j]+buf[i][j-1]+recvRight[i])/4;
	       		}
	       		else
	       		{
	       			t[i][j]=(buf[i-1][j]+buf[i+1][j]+buf[i][j-1]+buf[i][j+1])/4;
	       		}
	       	}
	       }
	}
	else if(col==p-1)
	{
	       commRank_L = rank-1;
	       commRank_B = rank+p;
	       double recvLeft[N],recvDown[N];
	       
       	 MPI_Send(&buf[0][0], 1, newvtypeRow, commRank_L , rank, MPI_COMM_WORLD);
       	 MPI_Recv(recvLeft, N, MPI_DOUBLE, commRank_L, commRank_L, MPI_COMM_WORLD, &status);
	       
       	 MPI_Send(&buf[N-1][0], 1, newvtypeCol,commRank_B , rank, MPI_COMM_WORLD);
       	 MPI_Recv(recvDown, N, MPI_DOUBLE, commRank_B, commRank_B, MPI_COMM_WORLD, &status);

	       //local Computation
	       for(int i=0;i<N;i++)
	       {
	       	for(int j=0;j<N;j++)
	       	{
	       		if(i==0)
	       		{
	       			if(j==0)
	       			{
	       				t[i][j]=(buf[i][j+1]+buf[i+1][j]+recvLeft[i])/3;
	       			}
	       			else if(j==N-1)
	       			{
	       				t[i][j]=(buf[i][j-1]+buf[i+1][j])/2;
	       			}
	       			else
	       			{
	       				t[i][j]=(buf[i][j+1]+buf[i][j-1]+buf[i+1][j])/3;
	       			}
	       		}
	       		else if(i==N-1)
	       		{
	       			if(j==0)
	       			{
	       				t[i][j]=(buf[i][j+1]+buf[i-1][j]+recvDown[j]+recvLeft[i])/4;
	       			}
	       			else if(j==N-1)
	       			{
	       				t[i][j]=(buf[i][j-1]+buf[i-1][j]+recvDown[j])/3;
	       			}
	       			else
	       			{
	       				t[i][j]=(buf[i][j+1]+buf[i][j-1]+buf[i-1][j]+recvDown[j])/4;
	       			}
	       		}
	       		else if(j==0)
	       		{
	       			t[i][j]=(buf[i-1][j]+buf[i+1][j]+buf[i][j+1]+recvLeft[i])/4;
	       		}
	       		else if(j==N-1)
	       		{
	       			t[i][j]=(buf[i-1][j]+buf[i+1][j]+buf[i][j-1])/3;
	       		}
	       		else
	       		{
	       			t[i][j]=(buf[i-1][j]+buf[i+1][j]+buf[i][j-1]+buf[i][j+1])/4;
	       		}
	       	}
	       }
	}
	else
        {
               commRank_L = rank-1;
               commRank_R = rank+1;
               commRank_B = rank+p;
               double recvLeft[N],recvDown[N],recvRight[N];
               
       	 MPI_Send(&buf[0][N-1], 1, newvtypeRow,commRank_R , rank, MPI_COMM_WORLD);
       	 MPI_Recv(recvRight, N, MPI_DOUBLE, commRank_R, commRank_R, MPI_COMM_WORLD, &status);

       	 MPI_Send(&buf[0][0], 1, newvtypeRow,commRank_L , rank, MPI_COMM_WORLD);
       	 MPI_Recv(recvLeft, N, MPI_DOUBLE, commRank_L, commRank_L, MPI_COMM_WORLD, &status);

       	 MPI_Send(&buf[N-1][0], 1, newvtypeCol,commRank_B , rank, MPI_COMM_WORLD);
       	 MPI_Recv(recvDown, N, MPI_DOUBLE, commRank_B, commRank_B, MPI_COMM_WORLD, &status);

       	//local Computation
       	for(int i=0;i<N;i++)
	       {
	       	for(int j=0;j<N;j++)
	       	{
	       		if(i==0)
	       		{
	       			if(j==0)
	       			{
	       				t[i][j]=(buf[i][j+1]+buf[i+1][j]+recvLeft[i])/3;
	       			}
	       			else if(j==N-1)
	       			{
	       				t[i][j]=(buf[i][j-1]+buf[i+1][j]+recvRight[i])/3;
	       			}
	       			else
	       			{
	       				t[i][j]=(buf[i][j+1]+buf[i][j-1]+buf[i+1][j])/3;
	       			}
	       		}
	       		else if(i==N-1)
	       		{
	       			if(j==0)
	       			{
	       				t[i][j]=(buf[i][j+1]+buf[i-1][j]+recvDown[j]+recvLeft[i])/4;
	       			}
	       			else if(j==N-1)
	       			{
	       				t[i][j]=(buf[i][j-1]+buf[i-1][j]+recvDown[j]+recvRight[i])/4;
	       			}
	       			else
	       			{
	       				t[i][j]=(buf[i][j+1]+buf[i][j-1]+buf[i-1][j]+recvDown[j])/4;
	       			}
	       		}
	       		else if(j==0)
	       		{
	       			t[i][j]=(buf[i-1][j]+buf[i+1][j]+buf[i][j+1]+recvLeft[i])/4;
	       		}
	       		else if(j==N-1)
	       		{
	       			t[i][j]=(buf[i-1][j]+buf[i+1][j]+buf[i][j-1]+recvRight[i])/4;
	       		}
	       		else
	       		{
	       			t[i][j]=(buf[i-1][j]+buf[i+1][j]+buf[i][j-1]+buf[i][j+1])/4;
	       		}
	       	}
	       }       	
       }
} 
else if(row==p-1) 
{
        if(col==0)
        { 
               commRank_R = rank+1;
               commRank_T = rank-p;
	       double recvRight[N],recvTop[N];
	       
        	 MPI_Send(&buf[0][N-1], 1, newvtypeRow,commRank_R , rank, MPI_COMM_WORLD);
       	 MPI_Recv(recvRight, N, MPI_DOUBLE, commRank_R, commRank_R, MPI_COMM_WORLD, &status);

       	 MPI_Send(&buf[0][0], 1, newvtypeCol,commRank_T , rank, MPI_COMM_WORLD);
       	 MPI_Recv(recvTop, N, MPI_DOUBLE, commRank_T, commRank_T, MPI_COMM_WORLD, &status);

		//local Computation
		for(int i=0;i<N;i++)
	       {
	       	for(int j=0;j<N;j++)
	       	{
	       		if(i==0)
	       		{
	       			if(j==0)
	       			{
	       				t[i][j]=(buf[i][j+1]+buf[i+1][j]+recvTop[j])/3;
	       			}
	       			else if(j==N-1)
	       			{
	       				t[i][j]=(buf[i][j-1]+buf[i+1][j]+recvRight[i]+recvTop[j])/4;
	       			}
	       			else
	       			{
	       				t[i][j]=(buf[i][j+1]+buf[i][j-1]+buf[i+1][j]+recvTop[j])/4;
	       			}
	       		}
	       		else if(i==N-1)
	       		{
	       			if(j==0)
	       			{
	       				t[i][j]=(buf[i][j+1]+buf[i-1][j])/2;
	       			}
	       			else if(j==N-1)
	       			{
	       				t[i][j]=(buf[i][j-1]+buf[i-1][j]+recvRight[i])/3;
	       			}
	       			else
	       			{
	       				t[i][j]=(buf[i][j+1]+buf[i][j-1]+buf[i-1][j])/3;
	       			}
	       		}
	       		else if(j==0)
	       		{
	       			t[i][j]=(buf[i-1][j]+buf[i+1][j]+buf[i][j+1])/3;
	       		}
	       		else if(j==N-1)
	       		{
	       			t[i][j]=(buf[i-1][j]+buf[i+1][j]+buf[i][j-1]+recvRight[i])/4;
	       		}
	       		else
	       		{
	       			t[i][j]=(buf[i-1][j]+buf[i+1][j]+buf[i][j-1]+buf[i][j+1])/4;
	       		}
	       	}
	       }
        }
        else if(col==p-1)
        { 
               commRank_L = rank-1;
               commRank_T = rank-p;
               double recvLeft[N],recvTop[N];
	       
       	 MPI_Send(&buf[0][0], 1, newvtypeRow,commRank_L , rank, MPI_COMM_WORLD);
       	 MPI_Recv(recvLeft, N, MPI_DOUBLE, commRank_L, commRank_L, MPI_COMM_WORLD, &status);

       	 MPI_Send(&buf[0][0], 1, newvtypeCol,commRank_T , rank, MPI_COMM_WORLD);
       	 MPI_Recv(recvTop, N, MPI_DOUBLE, commRank_T, commRank_T, MPI_COMM_WORLD, &status);
	    
	     //local Computation
	     for(int i=0;i<N;i++)
	        {
	       	for(int j=0;j<N;j++)
	       	{
	       		if(i==0)
	       		{
	       			if(j==0)
	       			{
	       				t[i][j]=(buf[i][j+1]+buf[i+1][j]+recvLeft[i]+recvTop[j])/4;
	       			}
	       			else if(j==N-1)
	       			{
	       				t[i][j]=(buf[i][j-1]+buf[i+1][j]+recvTop[j])/3;
	       			}
	       			else
	       			{
	       				t[i][j]=(buf[i][j+1]+buf[i][j-1]+buf[i+1][j]+recvTop[j])/4;
	       			}
	       		}
	       		else if(i==N-1)
	       		{
	       			if(j==0)
	       			{
	       				t[i][j]=(buf[i][j+1]+buf[i-1][j]+recvLeft[i])/3;
	       			}
	       			else if(j==N-1)
	       			{
	       				t[i][j]=(buf[i][j-1]+buf[i-1][j])/2;
	       			}
	       			else
	       			{
	       				t[i][j]=(buf[i][j+1]+buf[i][j-1]+buf[i-1][j])/3;
	       			}
	       		}
	       		else if(j==0)
	       		{
	       			t[i][j]=(buf[i-1][j]+buf[i+1][j]+buf[i][j+1]+recvLeft[i])/4;
	       		}
	       		else if(j==N-1)
	       		{
	       			t[i][j]=(buf[i-1][j]+buf[i+1][j]+buf[i][j-1])/3;
	       		}
	       		else
	       		{
	       			t[i][j]=(buf[i-1][j]+buf[i+1][j]+buf[i][j-1]+buf[i][j+1])/4;
	       		}
	       	}
	       }
        }
        else
        { 
               commRank_L = rank-1;
               commRank_R = rank+1;
               commRank_T = rank-p;
               double recvLeft[N],recvTop[N],recvRight[N];
              
       	 MPI_Send(&buf[0][N-1], 1, newvtypeRow,commRank_R , rank, MPI_COMM_WORLD);
       	 MPI_Recv(recvRight, N, MPI_DOUBLE, commRank_R, commRank_R, MPI_COMM_WORLD, &status);

       	 MPI_Send(&buf[0][0], 1, newvtypeRow,commRank_L , rank, MPI_COMM_WORLD);
       	 MPI_Recv(recvLeft, N, MPI_DOUBLE, commRank_L, commRank_L, MPI_COMM_WORLD, &status);

       	 MPI_Send(&buf[0][0], 1, newvtypeCol,commRank_T , rank, MPI_COMM_WORLD);
       	 MPI_Recv(recvTop, N, MPI_DOUBLE, commRank_T, commRank_T, MPI_COMM_WORLD, &status);

       	//local Computation
       	for(int i=0;i<N;i++)
	       {
	       	for(int j=0;j<N;j++)
	       	{
	       		if(i==0)
	       		{
	       			if(j==0)
	       			{
	       				t[i][j]=(buf[i][j+1]+buf[i+1][j]+recvLeft[i]+recvTop[j])/4;
	       			}
	       			else if(j==N-1)
	       			{
	       				t[i][j]=(buf[i][j-1]+buf[i+1][j]+recvRight[i]+recvTop[j])/4;
	       			}
	       			else
	       			{
	       				t[i][j]=(buf[i][j+1]+buf[i][j-1]+buf[i+1][j]+recvTop[j])/4;
	       			}
	       		}
	       		else if(i==N-1)
	       		{
	       			if(j==0)
	       			{
	       				t[i][j]=(buf[i][j+1]+buf[i-1][j]+recvLeft[i])/3;
	       			}
	       			else if(j==N-1)
	       			{
	       				t[i][j]=(buf[i][j-1]+buf[i-1][j]+recvRight[i])/3;
	       			}
	       			else
	       			{
	       				t[i][j]=(buf[i][j+1]+buf[i][j-1]+buf[i-1][j])/3;
	       			}
	       		}
	       		else if(j==0)
	       		{
	       			t[i][j]=(buf[i-1][j]+buf[i+1][j]+buf[i][j+1]+recvLeft[i])/4;
	       		}
	       		else if(j==N-1)
	       		{
	       			t[i][j]=(buf[i-1][j]+buf[i+1][j]+buf[i][j-1]+recvRight[i])/4;
	       		}
	       		else
	       		{
	       			t[i][j]=(buf[i-1][j]+buf[i+1][j]+buf[i][j-1]+buf[i][j+1])/4;
	       		}
	       	}
	       }
       }
} 
else if(col==0)
{
		commRank_B = rank+p;
                commRank_R = rank+1;
                commRank_T = rank-p;

		double recvDown[N],recvTop[N],recvRight[N];
               
       	 MPI_Send(&buf[0][N-1], 1, newvtypeRow,commRank_R , rank, MPI_COMM_WORLD);
       	 MPI_Recv(recvRight, N, MPI_DOUBLE, commRank_R, commRank_R, MPI_COMM_WORLD, &status);

       	 MPI_Send(&buf[N-1][0], 1, newvtypeCol, commRank_B , rank, MPI_COMM_WORLD);
       	 MPI_Recv(recvDown, N, MPI_DOUBLE, commRank_B, commRank_B, MPI_COMM_WORLD, &status);

       	 MPI_Send(&buf[0][0], 1, newvtypeCol,commRank_T , rank, MPI_COMM_WORLD);
       	 MPI_Recv(recvTop, N, MPI_DOUBLE, commRank_T, commRank_T, MPI_COMM_WORLD, &status);

		//local Computation       
	       for(int i=0;i<N;i++)
	       {
	       	for(int j=0;j<N;j++)
	       	{
	       		if(i==0)
	       		{
	       			if(j==0)
	       			{
	       				t[i][j]=(buf[i][j+1]+buf[i+1][j]+recvTop[j])/3;
	       			}
	       			else if(j==N-1)
	       			{
	       				t[i][j]=(buf[i][j-1]+buf[i+1][j]+recvTop[j])/3;
	       			}
	       			else
	       			{
	       				t[i][j]=(buf[i][j+1]+buf[i][j-1]+buf[i+1][j]+recvTop[j])/4;
	       			}
	       		}
	       		else if(i==N-1)
	       		{
	       			if(j==0)
	       			{
	       				t[i][j]=(buf[i][j+1]+buf[i-1][j]+recvDown[j])/3;
	       			}
	       			else if(j==N-1)
	       			{
	       				t[i][j]=(buf[i][j-1]+buf[i-1][j]+recvRight[i]+recvDown[j])/4;
	       			}
	       			else
	       			{
	       				t[i][j]=(buf[i][j+1]+buf[i][j-1]+buf[i-1][j]+recvDown[j])/4;
	       			}
	       		}
	       		else if(j==0)
	       		{
	       			t[i][j]=(buf[i-1][j]+buf[i+1][j]+buf[i][j+1])/3;
	       		}
	       		else if(j==N-1)
	       		{
	       			t[i][j]=(buf[i-1][j]+buf[i+1][j]+buf[i][j-1]+recvRight[i])/4;
	       		}
	       		else
	       		{
	       			t[i][j]=(buf[i-1][j]+buf[i+1][j]+buf[i][j-1]+buf[i][j+1])/4;
	       		}
	       	}
	       }
}
else if(col==p-1)
{
		commRank_L = rank-1;
                commRank_B = rank+p;
                commRank_T = rank-p;
		double recvDown[N],recvTop[N],recvLeft[N];
               
       	 MPI_Send(&buf[0][0], 1, newvtypeRow, commRank_L , rank, MPI_COMM_WORLD);
       	 MPI_Recv(recvLeft, N, MPI_DOUBLE, commRank_L, commRank_L, MPI_COMM_WORLD, &status);

       	 MPI_Send(&buf[N-1][0], 1, newvtypeCol, commRank_B , rank, MPI_COMM_WORLD);
       	 MPI_Recv(recvDown, N, MPI_DOUBLE, commRank_B, commRank_B, MPI_COMM_WORLD, &status);

       	 MPI_Send(&buf[0][0], 1, newvtypeCol,commRank_T , rank, MPI_COMM_WORLD);
       	 MPI_Recv(recvTop, N, MPI_DOUBLE, commRank_T, commRank_T, MPI_COMM_WORLD, &status);
		
		//local Computation
	       for(int i=0;i<N;i++)
	       {
	       	for(int j=0;j<N;j++)
	       	{
	       		if(i==0)
	       		{
	       			if(j==0)
	       			{
	       				t[i][j]=(buf[i][j+1]+buf[i+1][j]+recvTop[j]+recvLeft[i])/4;
	       			}
	       			else if(j==N-1)
	       			{
	       				t[i][j]=(buf[i][j-1]+buf[i+1][j]+recvTop[j])/3;
	       			}
	       			else
	       			{
	       				t[i][j]=(buf[i][j+1]+buf[i][j-1]+buf[i+1][j]+recvTop[j])/4;
	       			}
	       		}
	       		else if(i==N-1)
	       		{
	       			if(j==0)
	       			{
	       				t[i][j]=(buf[i][j+1]+buf[i-1][j]+recvDown[j]+recvLeft[i])/4;
	       			}
	       			else if(j==N-1)
	       			{
	       				t[i][j]=(buf[i][j-1]+buf[i-1][j]+recvDown[j])/3;
	       			}
	       			else
	       			{
	       				t[i][j]=(buf[i][j+1]+buf[i][j-1]+buf[i-1][j]+recvDown[j])/4;
	       			}
	       		}
	       		else if(j==0)
	       		{
	       			t[i][j]=(buf[i-1][j]+buf[i+1][j]+buf[i][j+1]+recvLeft[i])/4;
	       		}
	       		else if(j==N-1)
	       		{
	       			t[i][j]=(buf[i-1][j]+buf[i+1][j]+buf[i][j-1])/3;
	       		}
	       		else
	       		{
	       			t[i][j]=(buf[i-1][j]+buf[i+1][j]+buf[i][j-1]+buf[i][j+1])/4;
	       		}
	       	}
	       }
}
else
{
		commRank_L = rank-1;
                commRank_R = rank+1;
                commRank_T = rank-p;
		commRank_B = rank+p;

		double recvDown[N],recvTop[N],recvLeft[N],recvRight[N];
               
       	 MPI_Send(&buf[0][N-1], 1, newvtypeRow,commRank_R , rank, MPI_COMM_WORLD);
       	 MPI_Recv(recvRight, N, MPI_DOUBLE, commRank_R, commRank_R, MPI_COMM_WORLD, &status);

       	 MPI_Send(&buf[0][0], 1, newvtypeRow, commRank_L , rank, MPI_COMM_WORLD);
       	 MPI_Recv(recvLeft, N, MPI_DOUBLE, commRank_L, commRank_L, MPI_COMM_WORLD, &status);

       	 MPI_Send(&buf[N-1][0], 1, newvtypeCol, commRank_B , rank, MPI_COMM_WORLD);
       	 MPI_Recv(recvDown, N, MPI_DOUBLE, commRank_B, commRank_B, MPI_COMM_WORLD, &status);

       	 MPI_Send(&buf[0][0], 1, newvtypeCol,commRank_T , rank, MPI_COMM_WORLD);
       	 MPI_Recv(recvTop, N, MPI_DOUBLE, commRank_T, commRank_T, MPI_COMM_WORLD, &status);

	       //local Computation
	       for(int i=0;i<N;i++)
	       {
	       	for(int j=0;j<N;j++)
	       	{
	       		if(i==0)
	       		{
	       			if(j==0)
	       			{
	       				t[i][j]=(buf[i][j+1]+buf[i+1][j]+recvTop[j]+recvLeft[i])/4;
	       			}
	       			else if(j==N-1)
	       			{
	       				t[i][j]=(buf[i][j-1]+buf[i+1][j]+recvTop[j]+recvRight[i])/4;
	       			}
	       			else
	       			{
	       				t[i][j]=(buf[i][j+1]+buf[i][j-1]+buf[i+1][j]+recvTop[j])/4;
	       			}
	       		}
	       		else if(i==N-1)
	       		{
	       			if(j==0)
	       			{
	       				t[i][j]=(buf[i][j+1]+buf[i-1][j]+recvDown[j]+recvLeft[i])/4;
	       			}
	       			else if(j==N-1)
	       			{
	       				t[i][j]=(buf[i][j-1]+buf[i-1][j]+recvDown[j]+recvRight[i])/4;
	       			}
	       			else
	       			{
	       				t[i][j]=(buf[i][j+1]+buf[i][j-1]+buf[i-1][j]+recvDown[j])/4;
	       			}
	       		}
	       		else if(j==0)
	       		{
	       			t[i][j]=(buf[i-1][j]+buf[i+1][j]+buf[i][j+1]+recvLeft[i])/4;
	       		}
	       		else if(j==N-1)
	       		{
	       			t[i][j]=(buf[i-1][j]+buf[i+1][j]+buf[i][j-1]+recvRight[i])/4;
	       		}
	       		else
	       		{
	       			t[i][j]=(buf[i-1][j]+buf[i+1][j]+buf[i][j-1]+buf[i][j+1])/4;
	       		}
	       	}
	       }
}

for(int i=0;i<N;i++)
{
	for(int j=0;j<N;j++)
		buf[i][j]=t[i][j];
}
}
eTime=MPI_Wtime();   //timer stopped
time=eTime-sTime;
return time;
/*for(int i=0;i<N;i++)
{
	for(int j=0;j<N;j++)
		printf("%d  %lf ,   ",rank,buf[i][j]);
printf("\n");
}*/

}//end of mpiDerived

