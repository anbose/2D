#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <time.h>
#include <mpi.h>


void send_boundary(const int N, const int count, const int bnum, const double *localrho, double *leftboundary, double *rightboundary){
    
    int world_size,world_rank;
    int leftpeer,rightpeer,tag_send,tag_recv;

    double *leftbuff, *rightbuff;

    MPI_Comm_size(MPI_COMM_WORLD,&world_size);
    MPI_Comm_rank(MPI_COMM_WORLD,&world_rank);

#ifdef DEBUG
    printf("DEBUG :  BoundaryExchange %d \n",world_rank);
#endif 

    leftpeer = (world_rank == 0) ? world_size - 1 : world_rank - 1;
    rightpeer = (world_rank == world_size - 1 ? 0 : world_rank + 1);
    tag_send = 0;
    tag_recv = tag_send;

    leftbuff = new double[bnum*N];
    rightbuff = new double[bnum*N];

    for (int i = 0; i < bnum*N; i++){
        leftbuff[i] = localrho[i];
        rightbuff[i] = localrho[count - bnum*N + i];
    }

    MPI_Sendrecv(leftbuff, bnum*N, MPI_DOUBLE, leftpeer, tag_send, rightboundary, bnum*N, MPI_DOUBLE, rightpeer, tag_recv, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    
    MPI_Barrier(MPI_COMM_WORLD);

    MPI_Sendrecv(rightbuff, bnum*N, MPI_DOUBLE, rightpeer, tag_send, leftboundary, bnum*N, MPI_DOUBLE, leftpeer, tag_recv, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    MPI_Barrier(MPI_COMM_WORLD);

    delete[] leftbuff;
    delete[] rightbuff;
}

void send_boundary_from_global(const int N, const int bnum, const int *firstindex, const double *globalrho, double *leftboundary, double *rightboundary){

    int world_size,world_rank;
    int leftpeer,rightpeer,tag_send,tag_recv;

    double *leftbuff, *rightbuff;

    MPI_Comm_size(MPI_COMM_WORLD,&world_size);
    MPI_Comm_rank(MPI_COMM_WORLD,&world_rank);

#ifdef DEBUG
    printf("DEBUG :  BoundaryExchangeGlobal %d \n",world_rank);
#endif

    leftpeer = (world_rank == 0) ? world_size - 1 : world_rank - 1;
    rightpeer = (world_rank == world_size - 1 ? 0 : world_rank + 1);
    tag_send = 0;
    tag_recv = tag_send;

    leftbuff = new double[bnum*N];
    rightbuff = new double[bnum*N];

    for (int i = 0; i < bnum*N; i++){
        leftbuff[i] = globalrho[firstindex[world_rank]*N+i];
        rightbuff[i] = globalrho[(firstindex[world_rank+1]-bnum)*N + i];
    }

    MPI_Sendrecv(leftbuff, bnum*N, MPI_DOUBLE, leftpeer, tag_send, rightboundary, bnum*N, MPI_DOUBLE, rightpeer, tag_recv, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    MPI_Barrier(MPI_COMM_WORLD);

    MPI_Sendrecv(rightbuff, bnum*N, MPI_DOUBLE, rightpeer, tag_send, leftboundary, bnum*N, MPI_DOUBLE, leftpeer, tag_recv, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    MPI_Barrier(MPI_COMM_WORLD);

    delete[] leftbuff;
    delete[] rightbuff;
}
