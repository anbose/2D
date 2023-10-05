#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <time.h>
#include <mpi.h>

void update_local_density(const int N, const int *firstindex, double *rhol, const double *lboundary, const double *rboundary, \
                            const double *f_xl, const double *f_yl, const double *Parameters){

    int world_rank;

    MPI_Comm_rank(MPI_COMM_WORLD,&world_rank);

#ifdef DEBUG
    printf("DEBUG :  DensityUpdate %d \n",world_rank);
#endif
    
    int ylim = (firstindex[world_rank+1] - firstindex[world_rank]);

    double *xFlux = new double[ylim*N];
    double *yFlux = new double[ylim*N];
    double *leftfluxboundary = new double[N];
    double *rightfluxboundary = new double[N];
    double *yFluxbounded = new double[(ylim+2)*N];

    double *delta_rho = new double[ylim*N];
    double *rho_bounded = new double[N*(ylim+4)];

    merge_boundary(N,ylim,2,rhol,rho_bounded,lboundary,rboundary);

    calculate_flux(N,ylim,xFlux,yFlux,rho_bounded,f_xl,f_yl,Parameters);

    MPI_Barrier(MPI_COMM_WORLD);

    send_boundary(N,ylim*N,1,yFlux,leftfluxboundary,rightfluxboundary);

    MPI_Barrier(MPI_COMM_WORLD);

    merge_boundary(N,ylim,1,yFlux,yFluxbounded,leftfluxboundary,rightfluxboundary);

    for(int j=1;j<=ylim;j++){
        for(int i=0;i<N;i++){
            
            int index_plus = (i==N-1?0:i+1);
            int index_minus = (i==0?N-1:i-1);

            delta_rho[(j-1)*N+i] = 0.5*Parameters[P_dx]*(xFlux[(j-1)*N+index_minus] - xFlux[(j-1)*N+index_plus]);
            delta_rho[(j-1)*N+i] += 0.5*Parameters[P_dy]*(yFluxbounded[(j-1)*N+i] - yFluxbounded[(j+1)*N+i]);

            rhol[(j-1)*N+i] += 0.5*Parameters[P_dt]*delta_rho[(j-1)*N+i];
        }
    }

    delete[] xFlux; delete[] yFlux;
    delete[] leftfluxboundary; delete[] rightfluxboundary; delete[] yFluxbounded;
    delete[] delta_rho; delete[] rho_bounded;
}
