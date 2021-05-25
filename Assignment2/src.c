#include <stdio.h>
#include "mpi.h"
#include <stdlib.h>

double MPI_Bcast_optimized(int no_of_Groups,int D,int ppn,int nodes,int NodesPerGroup);
double MPI_Gather_optimized(int no_of_Groups,int D,int ppn,int nodes,int NodesPerGroup);
double MPI_Reduce_optimized(int no_of_Groups,int D,int ppn,int nodes,int NodesPerGroup);
double MPI_Alltoallv_optimized(int no_of_Groups,int D,int ppn,int nodes,int NodesPerGroup);

double MPI_Bcast_default(int no_of_Groups,int D,int ppn,int nodes,int NodesPerGroup);
double MPI_Gather_default(int no_of_Groups,int D,int ppn,int nodes,int NodesPerGroup);
double MPI_Reduce_default(int no_of_Groups,int D,int ppn,int nodes,int NodesPerGroup);
double MPI_Alltoallv_default(int no_of_Groups,int D,int ppn,int nodes,int NodesPerGroup);

int main(int argc, char *argv[]) 
{
  int size, rank, len, gRank;

  MPI_Init (&argc, &argv);
  MPI_Comm_size (MPI_COMM_WORLD, &size);
  MPI_Comm_rank (MPI_COMM_WORLD, &rank);

  MPI_Group W_group;
  MPI_Comm_group (MPI_COMM_WORLD, &W_group);
  
  int no_of_Groups=4;
  int D=atoi(argv[1]);
  int ppn=atoi(argv[2]);
  int nodes=size/ppn;
  int NodesPerGroup=nodes/no_of_Groups;
  
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
 	 fprintf(fptr,"c\tD\tnodes\tppn\tmode\ttime\n");      //Name of columns in data.txt file
     }
  }

  // Optimized Bcast
  double avgOptBcast = MPI_Bcast_optimized(no_of_Groups, D, ppn, nodes, NodesPerGroup);
  if(avgOptBcast>0)
  {
	printf("Bcast Opt : %lf \n",avgOptBcast);
  	int mode=1;
 	fprintf(fptr,"Bcast\t%d\t%d\t%d\t%d\t%lf\n",D,nodes,ppn,mode,avgOptBcast);
  }
  
  // Default Bcast
  double avgDefBcast = MPI_Bcast_default(no_of_Groups, D, ppn, nodes, NodesPerGroup);
  if(avgDefBcast>0)
  {
	printf("Bcast Def : %lf \n",avgDefBcast);
  	int mode=0;
 	fprintf(fptr,"Bcast\t%d\t%d\t%d\t%d\t%lf\n",D,nodes,ppn,mode,avgDefBcast);
  }
  
  // Optimized Gather
  double avgOptGather = MPI_Gather_optimized(no_of_Groups, D, ppn, nodes, NodesPerGroup);
  if(avgOptGather>0)
  {
	printf("Gather Opt : %lf \n",avgOptGather);
  	int mode=1;
 	fprintf(fptr,"Gather\t%d\t%d\t%d\t%d\t%lf\n",D,nodes,ppn,mode,avgOptGather);
  }
  
  //Default Gather
  double avgDefGather = MPI_Gather_default(no_of_Groups, D, ppn, nodes, NodesPerGroup);
  if(avgDefGather>0)
  {
	printf("Gather Def : %lf \n",avgDefGather);
  	int mode=0;
 	fprintf(fptr,"Gather\t%d\t%d\t%d\t%d\t%lf\n",D,nodes,ppn,mode,avgDefGather);
  }
  
  //Optimized Reduce
  double avgOptReduce = MPI_Reduce_optimized(no_of_Groups, D, ppn, nodes, NodesPerGroup);
  if(avgOptReduce>0)
  {
	printf("Reduce Opt : %lf \n",avgOptReduce);
  	int mode=1;
 	fprintf(fptr,"Reduce\t%d\t%d\t%d\t%d\t%lf\n",D,nodes,ppn,mode,avgOptReduce);
  }
  
  //Default Reduce
  double avgDefReduce = MPI_Reduce_default(no_of_Groups, D, ppn, nodes, NodesPerGroup);
  if(avgDefReduce>0)
  {
 	printf("Reduce Def : %lf \n",avgDefReduce);
  	int mode=0;
 	fprintf(fptr,"Reduce\t%d\t%d\t%d\t%d\t%lf\n",D,nodes,ppn,mode,avgDefReduce);
  }
  
  //Optimized Alltoallv
  double avgOptAlltoallv = MPI_Alltoallv_optimized(no_of_Groups, D, ppn, nodes, NodesPerGroup);
  if(avgOptAlltoallv>0)
  {
  	printf("Alltoallv Opt : %lf \n",avgOptAlltoallv);
  	int mode=1;
 	fprintf(fptr,"Alltoallv\t%d\t%d\t%d\t%d\t%lf\n",D,nodes,ppn,mode,avgOptAlltoallv);
  }
  
  //Default Alltoallv
  double avgDefAlltoallv = MPI_Alltoallv_default(no_of_Groups, D, ppn, nodes, NodesPerGroup);
  if(avgDefAlltoallv>0)
  {
  	printf("Alltoallv Def : %lf \n",avgDefAlltoallv);
  	int mode=0;
 	fprintf(fptr,"Alltoallv\t%d\t%d\t%d\t%d\t%lf\n",D,nodes,ppn,mode,avgDefAlltoallv);
  }

  
  MPI_Finalize();
  fclose(fptr);
  return 0;
}

double MPI_Bcast_optimized(int no_of_Groups,int D,int ppn,int nodes,int NodesPerGroup)
{
	int rank;
	MPI_Comm_rank (MPI_COMM_WORLD, &rank);
	int size;
	MPI_Comm_size( MPI_COMM_WORLD, &size );

	MPI_Group W_group;
  	MPI_Comm_group (MPI_COMM_WORLD, &W_group);
  	
  	double maxtime[5],avgTime=0;
  	int count = D*128;       // Converting into no. of doubles
  	
  	double *buf=(double*)malloc(sizeof(double)*count);
  	for (int i=0; i<count; i++)
  		buf[i] = i;
  		
  	for(int b=0; b<5; b++)
  	{
		double sTime = MPI_Wtime();
		
		// Optimization Type 2 for hostfile category 2 and 3
		int ranksComm1[nodes], i, j=-1;
  		for (i=0; i<nodes; i++)
  		{    ranksComm1[++j] = i;
  		}
  		int rankSizeComm1 = j+1;

  		// create a new group
  		MPI_Group newGroup1;
  		MPI_Group_incl (W_group, rankSizeComm1, ranksComm1, &newGroup1);
  
  		MPI_Comm newComm1;
  		MPI_Comm_create_group (MPI_COMM_WORLD, newGroup1, 123, &newComm1);
  
  		MPI_Comm newComm2;
  		MPI_Comm_split (MPI_COMM_WORLD, rank%nodes, rank, &newComm2); 
 		
 		if(rank<nodes)
 			MPI_Bcast(buf, count, MPI_DOUBLE, 0, newComm1);

  		MPI_Bcast(buf, count, MPI_DOUBLE, 0, newComm2);
  		
  		// Optimization Type 1 for hostfile category 1
  		/*
  		  int ranksComm1[no_of_Groups], i, j=-1;
		  for (i=0; i<size; i++)
		  {    if (i%(ppn*NodesPerGroup) == 0) ranksComm1[++j] = i;
		  }

		  int rankSizeComm1 = j+1;

		  // create a new group
		  MPI_Group newGroup1;
		  MPI_Group_incl (W_group, rankSizeComm1, ranksComm1, &newGroup1);
		  
		  MPI_Comm newComm1;
		  MPI_Comm_create_group (MPI_COMM_WORLD, newGroup1, 123, &newComm1);
	
		  int ranks[nodes];// i, 
		  j=-1;
		  for (i=0; i<size; i++)
		  {    if (i%ppn == 0) ranks[++j] = i;
		  }

		  int rankSize = j+1;

		  MPI_Group newGroup;
		  MPI_Group_incl (W_group, rankSize, ranks, &newGroup);
		  
		  MPI_Comm newComm;
		  MPI_Comm_create_group (MPI_COMM_WORLD, newGroup, 123, &newComm);

		  MPI_Comm newComm2;
		  if(rank%ppn==0)
		  	MPI_Comm_split (newComm, rank/(ppn*NodesPerGroup), rank, &newComm2);		
		 
		  MPI_Comm newComm3;
		  MPI_Comm_split (MPI_COMM_WORLD, rank/ppn, rank, &newComm3); 
		  
		  if(rank%(ppn*NodesPerGroup)==0)
		  	MPI_Bcast(buf, count, MPI_DOUBLE, 0, newComm1);

		  if(rank%ppn==0)
			MPI_Bcast(buf, count, MPI_DOUBLE, 0, newComm2);
		  
		  MPI_Bcast(buf, count, MPI_DOUBLE, 0, newComm3);
  		*/
  		
  		//Optimization type 3 for hostfile category 2 and 3
  		/*
		  MPI_Status status[nodes];
  		  MPI_Request request[nodes];  	
  		  if(rank==0)
		  {
		  	for (int r=1; r<nodes; r++)
		    		MPI_Isend (buf, count, MPI_DOUBLE, r, r, MPI_COMM_WORLD,&request[r-1]);
			MPI_Waitall(nodes-1,request,status);
		  }
		  else if(rank<nodes)
		  {
		  	MPI_Recv (buf, count, MPI_DOUBLE, 0, rank, MPI_COMM_WORLD, &status[0]);
		  }
		  
		  MPI_Comm newComm3;
		  MPI_Comm_split (MPI_COMM_WORLD, rank%nodes, rank, &newComm3); 
		    
		  
		  MPI_Bcast(buf, count, MPI_DOUBLE, 0, newComm3);
  		*/
  		
  		double eTime = MPI_Wtime();
  		
  		double time;
  		time = eTime-sTime; 
  		MPI_Reduce (&time, &maxtime[b], 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

  	}
  		if(!rank)
  		{
  			
  			for(int i=0;i<5;i++)
  			{
                	        //printf("%lf , ",maxtime[i]);
        	                avgTime += maxtime[i];
	                }
                //printf("\n");

  			avgTime = avgTime/5;
  			//printf("time is %lf \n",maxtime);
  		}	
  		
  		return avgTime;
}

double MPI_Gather_optimized(int no_of_Groups,int D,int ppn,int nodes,int NodesPerGroup)
{
	int rank;
    MPI_Comm_rank (MPI_COMM_WORLD, &rank);
    int size;
	MPI_Comm_size( MPI_COMM_WORLD, &size );

	MPI_Group W_group;
  	MPI_Comm_group (MPI_COMM_WORLD, &W_group);
  	
  	double maxtime[5],avgTime=0;
  	int count = D*128;
  	
  	double* recvMessage2=(double*)malloc(sizeof(double)*ppn*count);
  	double* recvMessage1=(double*)malloc(sizeof(double)*no_of_Groups*NodesPerGroup*ppn*count);
	double* recvMessage3=(double*)malloc(sizeof(double)*NodesPerGroup*ppn*count);
  	
  	double *buf=(double*)malloc(sizeof(double)*count);
  	for (int i=0; i<count; i++)
  		buf[i] = i;
  		
  	for(int b=0; b<5; b++)
  	{
		double sTime = MPI_Wtime();
		
		// Optimization Type 2 for hostfile category 2 and 3
		int ranksComm1[nodes], i, j=-1;
  		for (i=0; i<nodes; i++)
  		{    ranksComm1[++j] = i;
  		}
  		int rankSizeComm1 = j+1;

  		// create a new group
  		MPI_Group newGroup1;
  		MPI_Group_incl (W_group, rankSizeComm1, ranksComm1, &newGroup1);
  
  		MPI_Comm newComm1;
  		MPI_Comm_create_group (MPI_COMM_WORLD, newGroup1, 123, &newComm1);
  
  		MPI_Comm newComm2;
  		MPI_Comm_split (MPI_COMM_WORLD, rank%nodes, rank, &newComm2); 
 		
 		MPI_Gather(buf, count, MPI_DOUBLE, recvMessage2, count, MPI_DOUBLE, 0, newComm2);

  		if(rank<nodes)
			MPI_Gather(recvMessage2, ppn*count, MPI_DOUBLE, recvMessage1, ppn*count, MPI_DOUBLE, 0, newComm1);

		// Optimization Type 1 for hostfile category 1
		/*
		
		  int ranksComm1[no_of_Groups], i, j=-1;
		  for (i=0; i<size; i++)
		  {    if (i%(ppn*NodesPerGroup) == 0) ranksComm1[++j] = i;
		  }

		  int rankSizeComm1 = j+1;

		  // create a new group
		  MPI_Group newGroup1;
		  MPI_Group_incl (W_group, rankSizeComm1, ranksComm1, &newGroup1);
		  
		  MPI_Comm newComm1;
		  MPI_Comm_create_group (MPI_COMM_WORLD, newGroup1, 123, &newComm1);
	
		  int ranks[nodes];// i, 
		  j=-1;
		  for (i=0; i<size; i++)
		  {    if (i%ppn == 0) ranks[++j] = i;
		  }

		  int rankSize = j+1;

		  MPI_Group newGroup;
		  MPI_Group_incl (W_group, rankSize, ranks, &newGroup);
		  
		  MPI_Comm newComm;
		  MPI_Comm_create_group (MPI_COMM_WORLD, newGroup, 123, &newComm);

		  MPI_Comm newComm2;
		  if(rank%ppn==0)
		  	MPI_Comm_split (newComm, rank/(ppn*NodesPerGroup), rank, &newComm2);		
		 
		  MPI_Comm newComm3;
		  MPI_Comm_split (MPI_COMM_WORLD, rank/ppn, rank, &newComm3);
		
		MPI_Gather(buf, count, MPI_DOUBLE, recvMessage2, count, MPI_DOUBLE, 0, newComm3);

		 if(rank%ppn==0)
			MPI_Gather(recvMessage2, ppn*count, MPI_DOUBLE, recvMessage3, ppn*count, MPI_DOUBLE, 0, newComm2);
		 
		 if(rank%(ppn*NodesPerGroup)==0)
	MPI_Gather(recvMessage3, NodesPerGroup*ppn*count, MPI_DOUBLE, recvMessage1, NodesPerGroup*ppn*count, MPI_DOUBLE, 0, newComm1);

		*/
  		
  		double eTime = MPI_Wtime();
  		
  		double time;
  		time = eTime-sTime; 
  		MPI_Reduce (&time, &maxtime[b], 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

  	}
  		if(!rank)
  		{
  			
  			for(int i=0;i<5;i++)
  			{
                	        //printf("%lf , ",maxtime[i]);
        	                avgTime += maxtime[i];
	                }
                //printf("\n");

  			avgTime = avgTime/5;
  			//printf("time is %lf \n",maxtime);
  		}	
  		
  		return avgTime;
}

double MPI_Reduce_optimized(int no_of_Groups,int D,int ppn,int nodes,int NodesPerGroup)
{
	int rank;
        MPI_Comm_rank (MPI_COMM_WORLD, &rank);
        int size;
	MPI_Comm_size( MPI_COMM_WORLD, &size );

	MPI_Group W_group;
  	MPI_Comm_group (MPI_COMM_WORLD, &W_group);
  	
  	double maxtime[5],avgTime=0;
  	int count = D*128;
  	
	double* recvMax3=(double*)malloc(sizeof(double)*count);
  	double* recvMax2=(double*)malloc(sizeof(double)*count);
  	double* recvMax=(double*)malloc(sizeof(double)*count);

  	
  	double *buf=(double*)malloc(sizeof(double)*count);
  	for (int i=0; i<count; i++)
  		buf[i] = i;
  		
  	for(int b=0; b<5; b++)
  	{
		double sTime = MPI_Wtime();
		
		// Optimization Type 2 for hostfile category 2 and 3
		int ranksComm1[nodes], i, j=-1;
  		for (i=0; i<nodes; i++)
  		{    ranksComm1[++j] = i;
  		}
  		int rankSizeComm1 = j+1;

  		// create a new group
  		MPI_Group newGroup1;
  		MPI_Group_incl (W_group, rankSizeComm1, ranksComm1, &newGroup1);
  
  		MPI_Comm newComm1;
  		MPI_Comm_create_group (MPI_COMM_WORLD, newGroup1, 123, &newComm1);
  
  		MPI_Comm newComm2;
  		MPI_Comm_split (MPI_COMM_WORLD, rank%nodes, rank, &newComm2); 
 		
 		MPI_Reduce(buf, recvMax2, count, MPI_DOUBLE, MPI_MAX, 0, newComm2);    

  		if(rank<nodes)
			MPI_Reduce(recvMax2, recvMax, count, MPI_DOUBLE, MPI_MAX, 0, newComm1); 	

		// Optimization Type 1 for hostfile category 1
  		/*
		  int ranksComm1[no_of_Groups], i, j=-1;
		  for (i=0; i<size; i++)
		  {    if (i%(ppn*NodesPerGroup) == 0) ranksComm1[++j] = i;
		  }

		  int rankSizeComm1 = j+1;

		  // create a new group
		  MPI_Group newGroup1;
		  MPI_Group_incl (W_group, rankSizeComm1, ranksComm1, &newGroup1);
		  
		  MPI_Comm newComm1;
		  MPI_Comm_create_group (MPI_COMM_WORLD, newGroup1, 123, &newComm1);
	
		  int ranks[nodes];// i, 
		  j=-1;
		  for (i=0; i<size; i++)
		  {    if (i%ppn == 0) ranks[++j] = i;
		  }

		  int rankSize = j+1;

		  MPI_Group newGroup;
		  MPI_Group_incl (W_group, rankSize, ranks, &newGroup);
		  
		  MPI_Comm newComm;
		  MPI_Comm_create_group (MPI_COMM_WORLD, newGroup, 123, &newComm);

		  MPI_Comm newComm2;
		  if(rank%ppn==0)
		  	MPI_Comm_split (newComm, rank/(ppn*NodesPerGroup), rank, &newComm2);		
		 
		  MPI_Comm newComm3;
		  MPI_Comm_split (MPI_COMM_WORLD, rank/ppn, rank, &newComm3);
		
		  MPI_Reduce(buf, recvMax3, count, MPI_DOUBLE, MPI_MAX, 0, newComm3);    

		  if(rank%ppn==0)
			MPI_Reduce(recvMax3, recvMax2, count, MPI_DOUBLE, MPI_MAX, 0, newComm2);

		  if(rank%(ppn*NodesPerGroup)==0)
			MPI_Reduce(recvMax2, recvMax, count, MPI_DOUBLE, MPI_MAX, 0, newComm1);   	

		*/
  		
  		double eTime = MPI_Wtime();
  		
  		double time;
  		time = eTime-sTime; 
  		MPI_Reduce (&time, &maxtime[b], 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

  	}
  		if(!rank)
  		{
  			
  			for(int i=0;i<5;i++)
  			{
                	        //printf("%lf , ",maxtime[i]);
        	                avgTime += maxtime[i];
	                }
                //printf("\n");

  			avgTime = avgTime/5;
  			//printf("time is %lf \n",maxtime);
  		}	
  		
  		return avgTime;
}

double MPI_Alltoallv_optimized(int no_of_Groups,int D,int ppn,int nodes,int NodesPerGroup)
{
	int rank,size;
        MPI_Comm_rank (MPI_COMM_WORLD, &rank);
	MPI_Comm_size( MPI_COMM_WORLD, &size );
  	
  	double maxtime[5],avgTime=0;
  	int count = D*128;
  	
  	int sCount[size],sDispls[size],rCount[size],rDispls[size];    
    	long dis=0;
    	int elements=count/size;
    	for (int i=0; i<size; i++) 
    	{
      	 	 sCount[i] = elements;
      	 	 rCount[i] = elements;
      	 	 rDispls[i] = dis ;
      	 	 sDispls[i] = dis;
       	 dis+=elements;
    	}
   	double *sBuf=(double*)malloc(sizeof(double)*dis);

	double *rBuf=(double*)malloc(sizeof(double)*dis);

    	for(long i=0;i<dis;i++)
	{
		sBuf[i]=rank;
		rBuf[i]=-1;
	}
  	
  	for(int b=0; b<5; b++)
  	{
		double sTime = MPI_Wtime();
		for(int i=0;i<size;i++)
 			MPI_Scatterv (sBuf, sCount,rDispls, MPI_DOUBLE, rBuf+elements*i, elements, MPI_DOUBLE, i, MPI_COMM_WORLD);
  		double eTime = MPI_Wtime();
  		
  		double time;
  		time = eTime-sTime; 
  		MPI_Reduce (&time, &maxtime[b], 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

  	}
  		if(!rank)
  		{
  			
  			for(int i=0;i<5;i++)
  			{
                	        //printf("%lf , ",maxtime[i]);
        	                avgTime += maxtime[i];
	                }
                 	//printf("\n");

  			avgTime = avgTime/5;
  			//printf("time is %lf \n",maxtime);
  		}	
  		
  		return avgTime;
}

// Defaults


double MPI_Bcast_default(int no_of_Groups,int D,int ppn,int nodes,int NodesPerGroup)
{
	int rank;
        MPI_Comm_rank (MPI_COMM_WORLD, &rank);

	double maxtime[5],avgTime=0;
  	int count = D*128;
  	
	double *buf=(double*)malloc(sizeof(double)*count);
  	for (int i=0; i<count; i++)
  		buf[i] = i;
  		
  	for(int b=0; b<5; b++)
  	{
  		double sTime = MPI_Wtime();
  		MPI_Bcast(buf, count, MPI_DOUBLE, 0, MPI_COMM_WORLD);
 		double eTime = MPI_Wtime();
 		
 		double time;
  		time = eTime-sTime; 
  		MPI_Reduce (&time, &maxtime[b], 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  	}
  	
  	if(!rank)
  	{
  			
  		for(int i=0;i<5;i++)
  		{
                        //printf("%lf , ",maxtime[i]);
                        avgTime += maxtime[i];
                }
                //printf("\n");

  		avgTime = avgTime/5;
  		//printf("time is %lf \n",maxtime);
  	}	
  		
  	return avgTime;
}

double MPI_Gather_default(int no_of_Groups,int D,int ppn,int nodes,int NodesPerGroup)
{
	int rank;
        MPI_Comm_rank (MPI_COMM_WORLD, &rank);

	double maxtime[5],avgTime=0;
  	int count = D*128;
  	
  	double* recvMessage1=(double*)malloc(sizeof(double)*no_of_Groups*NodesPerGroup*ppn*count);
	double *buf=(double*)malloc(sizeof(double)*count);
  	for (int i=0; i<count; i++)
  		buf[i] = i;
  		
  	for(int b=0; b<5; b++)
  	{
  		double sTime = MPI_Wtime();
  		MPI_Gather(buf, count, MPI_DOUBLE, recvMessage1, count, MPI_DOUBLE, 0, MPI_COMM_WORLD);
 		double eTime = MPI_Wtime();
 		
 		double time;
  		time = eTime-sTime; 
  		MPI_Reduce (&time, &maxtime[b], 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  	}
  	
  	if(!rank)
  	{
  			
  		for(int i=0;i<5;i++)
  		{
			//printf("%lf , ",maxtime[i]);
			avgTime += maxtime[i];
		}
		//printf("\n");
  		avgTime = avgTime/5;
  		//printf("time is %lf \n",maxtime);
  	}	
  		
  	return avgTime;
}

double MPI_Reduce_default(int no_of_Groups,int D,int ppn,int nodes,int NodesPerGroup)
{
	int rank;
        MPI_Comm_rank (MPI_COMM_WORLD, &rank);

	double maxtime[5],avgTime=0;
  	int count = D*128;
  	
  	double *max=(double*)malloc(sizeof(double)*count);
	double *buf=(double*)malloc(sizeof(double)*count);
  	for (int i=0; i<count; i++)
  		buf[i] = i;
  		
  	for(int b=0; b<5; b++)
  	{
  		double sTime = MPI_Wtime();
  		MPI_Reduce(buf, max, count, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
 		double eTime = MPI_Wtime();
 		
 		double time;
  		time = eTime-sTime; 
  		MPI_Reduce (&time, &maxtime[b], 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  	}
  	
  	if(!rank)
  	{
  			
  		for(int i=0;i<5;i++)
		{
			//printf("%lf , ",maxtime[i]);
  			avgTime += maxtime[i];
  		}
		//printf("\n");
		avgTime = avgTime/5;
  		//printf("time is %lf \n",maxtime);
  	}	
  		
  	return avgTime;
}

double MPI_Alltoallv_default(int no_of_Groups,int D,int ppn,int nodes,int NodesPerGroup)
{
	int rank,size;
        MPI_Comm_rank (MPI_COMM_WORLD, &rank);
        MPI_Comm_size( MPI_COMM_WORLD, &size );
        
        int count = D*128;
    
 	int sCount[size],sDispls[size],rCount[size],rDispls[size];    
    	long dis=0;
    	int elements=count/size;
    	for (int i=0; i<size; i++) 
    	{
      	 	 sCount[i] = elements;
      	 	 rCount[i] = elements;
      	 	 rDispls[i] =dis ;
      	 	 sDispls[i] =dis;
       	 dis+=elements;
    	}
   	double *sBuf=(double*)malloc(sizeof(double)*dis);

	double *rBuf=(double*)malloc(sizeof(double)*dis);

    	for(long i=0;i<dis;i++)
	{
		sBuf[i]=rank;
		rBuf[i]=-1;
	}


	double maxtime[5],avgTime=0;
  		
  	for(int b=0; b<5; b++)
  	{
  		double sTime = MPI_Wtime();
  		MPI_Alltoallv( sBuf, sCount, sDispls, MPI_DOUBLE,rBuf, rCount, rDispls, MPI_DOUBLE, MPI_COMM_WORLD );
 		double eTime = MPI_Wtime();
 		
 		double time;
  		time = eTime-sTime; 
  		MPI_Reduce (&time, &maxtime[b], 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  	}
  	
  	if(!rank)
  	{
  			
  		for(int i=0;i<5;i++)
		{
			//printf("%lf , ",maxtime[i]);
  			avgTime += maxtime[i];
  		}
		//printf("\n");
		avgTime = avgTime/5;
  		//printf("time is %lf \n",maxtime);
  	}	
  		
  	return avgTime;
}

