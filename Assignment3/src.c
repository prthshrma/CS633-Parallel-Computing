#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include "mpi.h"
#include <math.h>
#include <float.h>

#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))


int main(int argc, char *argv[])
{
	int rank, numtasks;
    	int row=0,column=0;
	if(argc!=2)
	{
		printf("Sorry, You entered Command Line Argument Wrongly..!!\n");
		exit(0);
	}
    	// Setup
    	MPI_Init(&argc, &argv);
    	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    	MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
	
	if(rank==0) 
        {
                // Opening file, creating file pointer for counting rows and columns
                FILE* fp = fopen(argv[1], "r");  
                if (!fp)
                {
                        printf("Can't open file\n");
                        return 0;
                }

                char line[1024];    // to store a line from data

                // Reading data line by line
                while (fgets(line,1024,fp)) 
                {
                    column = 0;
                    row++;

                    // to skip first row as it contains column names
                    if (row == 1)
                        continue;

                    // Splitting the data
                    char* value = strtok(line, ", ");

                    while (value) {
                        value = strtok(NULL, ", ");
                        column++;
                    }
                }
                // Close the file
                fclose(fp);  
         }
	MPI_Bcast(&row,1 , MPI_INT, 0, MPI_COMM_WORLD);  //broadcasting row to every other process
	MPI_Bcast(&column, 1, MPI_INT, 0, MPI_COMM_WORLD);   //broadcasting column to every other process
       
    	
	float *data = (float*)malloc(sizeof(float)*(row-1)*(column)); //array to store data from file 
	int i=0,j=0;	
	//using 1 process to read the entire data from a file    	
    	if(rank==0) 
    	{
		// Opening file, creating file pointer
    		FILE* fp = fopen(argv[1], "r");  
    		if (!fp)
    		{
        		printf("Can't open file\n");
			return 0;
    		}
	
        	char line[1024];    // to store a line from data

		// Reading data line by line
		while (fgets(line,1024,fp)) 
		{
		    j= 0;
		    i++;
		    
		    // to skip first row as it contains column names
		    if (i == 1)
		        continue;

		    // Splitting the data
		    char* value = strtok(line, ", ");
	  
		    while (value) {
			data[(i-2)*column+j] = atof(value);   // converting string into float
		        value = strtok(NULL, ", ");
		        j++;
		    }
		}
		// Close the file
		fclose(fp);  
		
	}
	    	int scatterRow;			// number of rows each process get on distribution
		scatterRow = (row-1)/numtasks;   

		//receive array for each process
		float *recvBuffer = (float*)malloc(sizeof(float)*scatterRow*column);
		
		MPI_Barrier(MPI_COMM_WORLD);

		//Timer starts
		double sTime = MPI_Wtime();
		
		// Distributing data to every process
		MPI_Scatter (data, scatterRow*column, MPI_FLOAT, recvBuffer, scatterRow*column, MPI_FLOAT, 0, MPI_COMM_WORLD);	
		
		float *minPerYear = (float*)malloc(sizeof(float)*(column-2));  // to store per year minimum on root
		float *minYearPerProcess = (float*)malloc(sizeof(float)*(column-2));  // to store per year minimum on every process
		for(int i=0; i<(column-2); i++)
			minYearPerProcess[i] = recvBuffer[i+2];
		
		// every process finding minimum per year
		for(int i=1; i<scatterRow; i++)
		{
			for(int j=0; j<(column-2); j++)
			{
		        	minYearPerProcess[j] = MIN(minYearPerProcess[j],recvBuffer[i*column + j+2]);
			}
		}
		
		// Reducing every year Minimum on root process  
		MPI_Reduce(minYearPerProcess, minPerYear, 41, MPI_FLOAT, MPI_MIN, 0, MPI_COMM_WORLD);
		
		float globalMin = FLT_MAX;  // to store global min across all stations and all years

		//opening output file
		FILE* fout = fopen("output.txt","w+");
		if(!rank)
		{
			for(int j=0; j<(column-2); j++)
			{
		                minPerYear[j] = MIN(minPerYear[j],data[(row-2)*column + j+2]);
				globalMin = MIN(globalMin,minPerYear[j]);
				fprintf(fout,"%0.2f,",minPerYear[j]);	//Writing yearwise minimum across all stations
			}
			fprintf(fout,"\n%0.2f\n",globalMin);     //Writing Global Minimum
		}
		// Timer stops
		double eTime = MPI_Wtime();
	  		
	  	double time, maxtime;
	  	time = eTime-sTime; 
	  	MPI_Reduce (&time, &maxtime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);   //Reduce used to calculate maxtime across all process

		if(!rank)
		{
			printf("maxtime : %lf , global min : %0.2f\n",maxtime,globalMin);
			fprintf(fout,"%lf\n",maxtime);  //Writing maximum time across process
		}
		fclose(fout);  //closing of output file pointer

	    // finalize
	    MPI_Finalize(); 
	    return 0;  
}

