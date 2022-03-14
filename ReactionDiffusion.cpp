#include <iostream>
#include <iomanip>
#include <fstream>
#include <mpi.h>
#include "ReactionDiffusion.h"

#define RD ReactionDiffusion
#define f1(u,v) (eps*u*(1.0-u)* (u - (v+b)/a))
#define f2(u,v) (u*u*u - v)

RD::ReactionDiffusion  (double dt, int T, int Nx, int Ny, double a,
                        double b, double mu1, double mu2, double eps,
                        double dx, double dy, int start, int end, int partition) {
    this->dt   = dt;
    this->T    = T;
    this->Nx   = Nx;
    this->Ny   = Ny;
    this->a    = a;
    this->b    = b;
    this->mu1  = mu1;
    this->mu2  = mu2;
    this->eps  = eps;
    this->dx   = dx;
    this->dy   = dy;
    this->start = start;        // Off-by-one due to left padding  (Ny)
    this->end   = end;          // Off-by-one due to right padding ((sz-1)*Ny)
    this->partition = partition;
    this->sz = end-start + 2; // Add 2 columns for left and right edge

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // Allocate memory for the matrix for each process, initialise to zero
    try {
        U1 = new double[(partition+4)*Ny]();
        U2 = new double[(partition+4)*Ny]();
        V1 = new double[(partition+4)*Ny]();
        V2 = new double[(partition+4)*Ny]();
    } catch (std::bad_alloc& ex) {
        std::cout << "Out of memory!" << std::endl;
        std::exit(1);
    }
}

RD::~ReactionDiffusion() {
    if(U1) delete U1;
    if(U2) delete U2;
    if(V1) delete V1;
    if(V2) delete V2;
}

void RD::SetInitialConditions() {
    // Set boundary conditions for u, for each slice
    for (int i = 1; i<sz-1; ++i) {
        for (int j = Ny/2+1; j<Ny; ++j) {
            U1[j+Ny*i] = 1.0;
        }
    }

    // Set boundary conditions for v, for each slice (if applicable)
    if (start < Nx/2) {
        for (int i = 1; i<(end > Nx/2 ? Nx/2-start+1 : sz-1); ++i) {
            for (int j = 0; j<Ny; ++j) {
                V1[j+Ny*i] = a/2.0;
            }
        }
    }
}

void RD::TimeIntegrate() {
    int timeSteps = int(T/dt);
    // Timesteps must be calculated sequentially
    for (int k = 0; k<timeSteps; ++k) {
        TimeIntegrateSingle();
        // if (rank==0) std::cout << "Completed timestep " << k << std::endl;
    }
}

void RD::TimeIntegrateSingle() {
    double mu1_val  = mu1/(dx*dy);
    double mu2_val  = mu2/(dx*dy);

    // Copy to local scope, avoid referencing
    double ddt = dt;
    int Nxx = Nx;
    int Nyy = Ny;

    for (int i = 1; i < sz-1; ++i) {
        for (int j = 0; j < Nyy; ++j) {
            int indx = j+Nyy*i;
            double u1_val = U1[indx];
            double v1_val = V1[indx];
            int multiplier = 4-(start==0 && i==1)-(j==0)-(end==Nxx && i==Nxx-1)-(j==Nyy-1);
            // Iterate over u matrix
            U2[indx] = ddt * (mu1_val * (
                                (end==Nxx && i==Nxx-1 ? 0 : U1[indx+Nyy]) +
                                (start==0 && i==1     ? 0 : U1[indx-Nyy]) +
                                (j < Nyy-1 ? U1[indx+1] : 0)   +
                                (j         ? U1[indx-1] : 0)   -
                                (multiplier)*u1_val)           +
                                f1(u1_val, v1_val)) + u1_val;
            // Iterate over v matrix
            V2[indx] = ddt * (mu2_val * (
                                (end==Nxx && i==Nxx-1 ? 0 : V1[indx+Nyy]) +
                                (start==0 && i==1     ? 0 : V1[indx-Nyy]) +
                                (j < Nyy-1 ? V1[indx+1] : 0)   +
                                (j         ? V1[indx-1] : 0)   -
                                (multiplier)*v1_val)           +
                                f2(u1_val, v1_val)) + v1_val;
        }
    }

    // Barrier to synchronise processes
    MPI_Barrier(MPI_COMM_WORLD);

    // Save current time step for next iteration
    for (int i = Ny; i < ((sz-1)*Ny); ++i) {
        U1[i] = U2[i];
        V1[i] = V2[i];
    }

    // Barrier to synchronise processes
    MPI_Barrier(MPI_COMM_WORLD);

    MPI_Request request;
    // Send right edge (1), non-blocking
    if (rank < size-1)  MPI_Isend(&U2[((sz-2)*Nyy)],    Nyy, MPI_DOUBLE, rank+1, 2, MPI_COMM_WORLD, &request);
    // Send left edge (2), non-blocking
    if (rank > 0)       MPI_Isend(&U2[Nyy],             Nyy, MPI_DOUBLE, rank-1, 1, MPI_COMM_WORLD, &request);
    // Receive right edge (1), blocking
    if (rank < size-1)  MPI_Recv(&U1[((sz-1)*Nyy)],     Nyy, MPI_DOUBLE, rank+1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    // Receive left edge (2), blocking
    if (rank > 0)       MPI_Recv(&U1[0],                Nyy, MPI_DOUBLE, rank-1, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    // Send right edge (1), non-blocking
    if (rank < size-1)  MPI_Isend(&V2[((sz-2)*Nyy)],    Nyy, MPI_DOUBLE, rank+1, 2, MPI_COMM_WORLD, &request);
    // Send left edge (2), non-blocking
    if (rank > 0)       MPI_Isend(&V2[Nyy],             Nyy, MPI_DOUBLE, rank-1, 1, MPI_COMM_WORLD, &request);
    // Receive right edge (1), blocking
    if (rank < size-1)  MPI_Recv(&V1[((sz-1)*Nyy)],     Nyy, MPI_DOUBLE, rank+1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    // Receive left edge (2), blocking
    if (rank > 0)       MPI_Recv(&V1[0],                Nyy, MPI_DOUBLE, rank-1, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    // Barrier to synchronise processes
    MPI_Barrier(MPI_COMM_WORLD);
}

void RD::writeOutput() {
    // Gather entire grid into one process
    if (rank == 0) {
        U_output = new double[Nx*Ny*2]();
        V_output = new double[Nx*Ny*2]();
    }

    int sendSize  = (partition)*Ny;

    if (rank == 0) {
        MPI_Gather(&U1[Ny], sendSize, MPI_DOUBLE, U_output, sendSize, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    } else {
        MPI_Gather(&U1[Ny], sendSize, MPI_DOUBLE, NULL, 0, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }

    if (rank == 0) {
        MPI_Gather(&V1[Ny], sendSize, MPI_DOUBLE, V_output, sendSize, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    } else {
        MPI_Gather(&V1[Ny], sendSize, MPI_DOUBLE, NULL, 0, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }

    // Barrier to ensure that root thread (0) has received all data
    MPI_Barrier(MPI_COMM_WORLD);

    // Save data to output.txt
    if (rank == 0) {
        std::ofstream outputFile("output.txt");
        for (int j = 0; j < Ny; ++j) {
            for (int i = 0; i < Nx; ++i) {
                outputFile << std::setw(5) << i << std::setw(5) << j << std::setw(15) << U_output[j + Ny * i] << std::setw(15) << V_output[j + Ny * i] << std::endl;
            }
            outputFile << std::endl;
        }
        outputFile.close();
    }
}