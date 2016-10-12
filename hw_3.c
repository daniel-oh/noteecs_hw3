#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>
#include <assert.h>

void split_data(int data_size, int world_rank,
		      int world_size, int* subdata_start,
		      int* subdata_size) {
  if (world_size > data_size) {
    // checks to see if more procs than data 
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  *subdata_start = data_size / world_size * world_rank;
  *subdata_size = data_size / world_size;
  if (world_rank == world_size - 1) {
    // Last processor gets remainder
    *subdata_size += data_size % world_size;
  }
}

void initialize_data(double** data, int subdata_start, int sub_final, int n) {
  int i, j;
  for (i = subdata_start; i <= sub_final; i++) {
    for (j = 0; j < n; j++) {
      data[i][j] = i*sin(i) + j*cos(j) + sqrt(i + j);
    }
  }
}

double minmax(double x) {
  double min = 100;
  double max = -100;
  if (x < min) {
    min = x;
  }
  if (min > max) {
    max = min;
  }
  return max;
}

double foo(double x) {
  double y;
  int i;
  y = x;
  for (i = 1; i <= 10; i++) {
    y = y+ sin(x*i)/pow(2.0,i);
  }
  return y;
}

double math(double** data, int i, int j) {
  double z = 0;
  
  //  z += foo(data[i-1][j]);
  // z += foo(data[i+1][j]);
  //  z += foo(data[i][j-1]);
  //  z += foo(data[i][j+1]);
  z += foo(data[i][j]);;
  z /= 5;

  return minmax(z);
}


double sum(double** data, int subdata_start, int n,  int sub_final) {
  int i, j;
  double sum = 0;
  for (i = subdata_start; i <= sub_final; i++) {
    for (j = 0; j < n; j++) {
      sum += data[i][j];
    }
  }
  return sum;
}


void send_datadown(double * send, int world_rank, int size) {
  MPI_Send(send, size, MPI_DOUBLE, world_rank + 1, 0, MPI_COMM_WORLD);
}

void send_dataup(double * send, int world_rank, int size) {
  MPI_Send(send, size, MPI_DOUBLE, world_rank - 1, 0, MPI_COMM_WORLD);
}

void recv_up(double * recv, int world_rank, int size) {
  MPI_Recv(recv, size, MPI_DOUBLE, world_rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
}

void recv_down(double * recv, int world_rank, int size) {
  MPI_Recv(recv, size, MPI_DOUBLE, world_rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
}
	   
int main(int argc, char** argv) {
  if (argc != 3) {
    fprintf(stderr, "Usage: hw3 num_m num_n\n");
    exit(1);
  }

 int m = atoi(argv[1]);
 int n = atoi(argv[2]);

 MPI_Init(NULL, NULL);

 int world_rank;
 int world_size;
 MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
 MPI_Comm_size(MPI_COMM_WORLD, &world_size);

 double total_time = 0.0;
 double** data;
 double** subdata;
 double send[n];
 double recv[n];
 int data_size;
 int i, j, k;
 int subdata_start, subdata_size;
 data = (double**)malloc(sizeof(double*)*m);
 for (i = 0; i < m; i++) {
   data[i] = (double*)malloc(sizeof(double)*n);
 }
 /*
 for(i = 0; i < m; i++) {
   for(j = 0; j < n; j++) {
     data[i][j] = 0;
   }
   }*/
 
 MPI_Barrier(MPI_COMM_WORLD);
 total_time -= MPI_Wtime();

 split_data(m, world_rank, world_size, &subdata_start, &subdata_size);
 int sub_final = subdata_start + subdata_size - 1;
 printf("Process %d of procs=%d covering data indexed in %d - %d\n", world_rank, world_size, subdata_start, sub_final);
 /*
 //empty subdata array
 subdata = (double**)malloc(sizeof(double*)*subdata_size);
 for (j = 0; j < subdata_size; j++) {
   subdata[j] = (double*)malloc(sizeof(double)*n);
   }*/

 initialize_data(data, subdata_start, sub_final, n);
 MPI_Barrier(MPI_COMM_WORLD);
 /*
 for (j = subdata_start + 1; j <= sub_final - 1; j++) {
   for (k = 0 + 1; k < n - 1; k++) {
     subdata[j][k] = math(data, j, k);
     printf("The data at [%d, %d] is: %f\n", j, k, subdata[j][k]);
   }
   }*/

 
 /*
 for (k = 1; k < n; k++) {
   subdata[2][k] = math(data, 2, k);
   printf("The data at [%d, %d] is: %f\n", 2, k, subdata[2][k]);
   }*/
 
 /*
 //run the main part here
 if (world_size == 1) {
   //do one processor math
   for 
 }
 else {
   //do other shit
   if (world_rank == 0) {
     //do non repeatable math
   }
   else if (world_rank == world_size - 1) {
     //do go up math
   }
   else if ((world_size / 2) - 1 == world_rank) {
     //send down then recv up
   }
   else if 
 }
 */

 if (world_rank == 0) {
   //do math
   for (i = subdata_start + 1; i < sub_final; i++) {
     for (j = 1; j < n - 1; j++) {
       subdata[i][j] = math(data, i, j);
     }
   }
   for (i = subdata_start + 1; i < sub_final; i++) {
     for (j = 1; j < n - 1; j++) {
       data[i][j] = subdata[i][j];
     }
   }
 }
   
 /*if (world_rank % 2 == 0 ) {
   //send message
   //receive message
   //do math
 }
 if (world_rank % 2 ==
 */

 
 double xsum = sum(data, subdata_start, n, sub_final);
 double xtotal = 0;
 double sqxtotal = 0;
 printf("Rank %d has sum %f\n", world_rank, xsum);

 MPI_Reduce(&xsum, &xtotal, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

 if (world_rank == 0) {
   printf("This is root (Rank %d) with total sum %f\n", world_rank, xtotal);
   sqxtotal = pow(xtotal, 2);
   printf("And the square of this is: %f\n", sqxtotal);
 }
 
 MPI_Barrier(MPI_COMM_WORLD);
 total_time += MPI_Wtime();
 if (world_rank == 0) {
   printf("Data size = %d, Matrix Size = %d\n", m * n * (int)sizeof(double), m * n);
   printf("Time it took = %lf\n", total_time);
 }

 for (j = 0; j < subdata_size; j++) {
   free(subdata);
 }
 free(subdata);
    
 
 for (i = 0; i < m; i++) {
   free(data[i]);
 }
 free(data);

 MPI_Finalize();

 return 0;
}
