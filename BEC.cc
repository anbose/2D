#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <time.h>
#include <mpi.h>

using namespace std;

int N, Nits, NSample;
double L, dx, dt, D0, M0, Teff, meanrho, rhoc;
int *firstindex, *counts, *dspls;
double *rho,*localrho,*leftboundary, *rightboundary;
double *force_x, *force_xl, *force_y, *force_yl;

double eps = 1.e-100;

int leftpeer,rightpeer,tag_send,tag_recv;

int nproc,myid;

#include "funcs.h"
#include "declarations.h"
#include "density_update.h"

int main(int argc, char* argv[]){

    MPI_Init(NULL,NULL);
    MPI_Comm_size(MPI_COMM_WORLD,&nproc);
    MPI_Comm_rank(MPI_COMM_WORLD,&myid);

    N = 128;
    L = 2;

    Nits = 2;
    dt = 1.e-4;
    NSample = 1;

    D0 = atof(argv[1]);
    M0 = 1.0;
    Teff = D0/M0;

    meanrho = 400.;
    rhoc = 1000.;

    set_parameters();
    if(myid==0)
	    printf("mean density : %f\n",meanrho);

    rho = new double[N*N];
    //localrho = new double[N*N];
    leftboundary = new double[2*N];
    rightboundary = new double[2*N];

    force_x = new double[N*N];
    //force_xl = new double[N*N];
    force_y = new double[N*N];
    //force_yl = new double[N*N];

    firstindex = new int[nproc+1];
    counts = new int[nproc];
    dspls = new int[nproc];


    initialize_global_rho(N,rho,Parameters); // set initial condition
    //read_file("data/rho_finalT1.bin",N);
    set_force(N,force_x,force_y);

    if(myid==0){
        printf("time %f : density at some point : %f \n", 0.0, rho[10]);
        save_rho("rho_initial.bin",N,rho);
    }

    ////////////////////////////////////////////////////// Distribute data to local variables //////////////////////////////////////////////////////

    dist_index(N,nproc,firstindex);

    MPI_Barrier(MPI_COMM_WORLD);

    for(int i=0;i<nproc;i++){
        counts[i] = N*(firstindex[i+1]-firstindex[i]);
        dspls[i] = N*firstindex[i];
    }

    localrho = new double[counts[myid]];
    force_xl = new double[counts[myid]];
    force_yl = new double[counts[myid]];

    MPI_Scatterv(rho,counts,dspls,MPI_DOUBLE,localrho,counts[myid],MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);

    MPI_Scatterv(force_x,counts,dspls,MPI_DOUBLE,force_xl,counts[myid],MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);

    MPI_Scatterv(force_y,counts,dspls,MPI_DOUBLE,force_yl,counts[myid],MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);

    ////////////////////////////////////////////////////// Distribute boundaries //////////////////////////////////////////////////////

    //send_boundary(N,counts[myid],localrho,leftboundary,rightboundary);

    ////////////////////////////////////////////////////// Update local rho //////////////////////////////////////////////////////

    int neg_events = 0;
    double global_sum = total_mass(N*N,rho,Parameters);
    //cout << "global sum is : " << global_sum << endl;
    double phys_time = 0.;
    bool global_neg;

    for(int t=1;t<=Nits;t++){

        send_boundary(N,counts[myid],2,localrho,leftboundary,rightboundary);

        /*if(t==1 && myid==0){
            printf("time %f : density at some point : %f \n", 0.0, rho[10]);
        }*/

        MPI_Barrier(MPI_COMM_WORLD);
        update_local_density(N,firstindex,localrho,leftboundary,rightboundary,force_xl,force_yl,Parameters);

        printf("node %d : local mass : %f \n",myid,total_mass(counts[myid],localrho,Parameters));
        
        bool cpt_neg = negative_density_check(counts[myid],localrho);
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Allreduce(&cpt_neg,&global_neg,1,MPI_C_BOOL,MPI_LOR,MPI_COMM_WORLD);

        if(global_neg){
            double local_sum = total_mass(counts[myid],rho,Parameters);
            MPI_Barrier(MPI_COMM_WORLD);
            MPI_Allreduce(&local_sum,&global_sum,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
            normalize_rho(counts[myid],global_sum,localrho,Parameters);
            neg_events += 1;
            //if(myid==0)
                //cout << "neg event found!" << endl;
        }

        MPI_Barrier(MPI_COMM_WORLD); 

        phys_time += Parameters[P_dt];

        if(t%NSample==0){
            MPI_Gatherv(localrho,counts[myid],MPI_DOUBLE,rho,counts,dspls,MPI_DOUBLE,0,MPI_COMM_WORLD);
            global_sum = total_mass(N*N,rho,Parameters);
            if(myid==0){
                printf("global sum is %f\n",global_sum);
                //printf("time %f : density at some point : %f \n", phys_time, localrho[10]);
            }
        }
    }

    ////////////////////////////////////////////////////// data handling //////////////////////////////////////////////////////

    MPI_Barrier(MPI_COMM_WORLD);

    MPI_Gatherv(localrho,counts[myid],MPI_DOUBLE,rho,counts,dspls,MPI_DOUBLE,0,MPI_COMM_WORLD);

    if(myid==0){
        global_sum = total_mass(N*N,rho,Parameters);
        printf("Parameters : \n");
        printf("Temp : %f, mean density : %f\n",Teff,meanrho);
        printf("Total mass : %f\n",global_sum);
        printf("negative density events : %f times of total run \n",(double)neg_events/Nits);
        save_rho("rho_final.bin",N,rho);
    }
    /*if(myid==1){
        printf("left vector for process %d is : \n", myid);
        for(int i=0;i<N;i++){
            printf("index %d : %f \n",i,localrho[i]);
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if(myid==0){
        printf("right vector for process %d is : \n", myid);
        for(int i=0;i<N;i++){
            printf("index %d : %f \n",i,rightboundary[i]);
        }
    }*/

    delete[] rho; delete[] localrho;
    delete[] leftboundary; delete[] rightboundary;

    delete[] firstindex; delete[] counts; delete[] dspls;
    delete[] force_x; delete[] force_xl; delete[] force_y; delete[] force_yl;

    MPI_Finalize();
    return 0;
}
