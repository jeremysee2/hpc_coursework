#include <iostream>
#include <iomanip>
#include <fstream>
#include "ReactionDiffusion.h"

#define RD ReactionDiffusion
#define f1(u,v) (eps*u*(1.0-u)* (u - (v+b)/a))
#define f2(u,v) (u*u*u - v)

RD::ReactionDiffusion  (double dt, int T, int Nx, int Ny, double a,
                        double b, double mu1, double mu2, double eps,
                        double dx, double dy) {
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

    // Allocate memory for the matrix, initialise to zero
    try {
        U1 = new double[Nx*Ny]();
        U2 = new double[Nx*Ny]();
        V1 = new double[Nx*Ny]();
        V2 = new double[Nx*Ny]();
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
    // Set boundary conditions for u
    for (int i = 0; i<Nx; ++i) {
        for (int j = Ny/2+1; j<Ny; ++j) {
            U1[j+Ny*i] = 1.0;
        }
    }

    // Set boundary conditions for v
    for (int i = 0; i<(Nx/2); ++i) {
        for (int j = 0; j<Ny; ++j) {
            V1[j+Ny*i] = a/2.0;
        }
    }
}

void RD::TimeIntegrate() {
    int timeSteps = int(T/dt);
    // Cannot parallelise this, due to race condition
    for (int k = 0; k<timeSteps; ++k) {
        TimeIntegrateSingle();
    }
}

void RD::TimeIntegrateSingle() {
    double mu1_val  = mu1/(dx*dy);
    double mu2_val  = mu2/(dx*dy);

    // Copy to local scope, avoid referencing
    double ddt = dt;
    int Nxx = Nx;
    int Nyy = Ny;

    #pragma omp parallel for schedule(static) collapse(2) num_threads(1)
    for (int i = 0; i < Nxx; ++i) {
        for (int j = 0; j < Nyy; ++j) {
            int indx = j+Nyy*i;
            double u1_val = U1[indx];
            double v1_val = V1[indx];
            int multiplier = 4-(i==0)-(j==0)-(i==Nxx-1)-(j==Nyy-1);
            // Iterate over u matrix
            U2[indx] = ddt * (mu1_val * (
                                (i < Nxx-1 ? U1[indx+Nyy] : 0) +
                                (i         ? U1[indx-Nyy] : 0) +
                                (j < Nyy-1 ? U1[indx+1] : 0)   +
                                (j         ? U1[indx-1] : 0)   -
                                (multiplier)*u1_val)           +
                                f1(u1_val, v1_val)) + u1_val;
            // Iterate over v matrix
            V2[indx] = ddt * (mu2_val * (
                                (i < Nxx-1 ? V1[indx+Nyy] : 0) +
                                (i         ? V1[indx-Nyy] : 0) +
                                (j < Nyy-1 ? V1[indx+1] : 0)   +
                                (j         ? V1[indx-1] : 0)   -
                                (multiplier)*v1_val)           +
                                f2(u1_val, v1_val)) + v1_val;
        }
    }

    // Save current time step for next iteration
    int sz = Nx*Ny;
    #pragma omp parallel for schedule(static) num_threads(1)
    for (int i = 0; i < sz; ++i) {
        U1[i] = U2[i];
        V1[i] = V2[i];
    }
}

void RD::writeOutput() {
    std::ofstream outputFile("output.txt");
    for (int j = 0; j < Ny; ++j) {
        for (int i = 0; i < Nx; ++i) {
            outputFile << std::setw(5) << i << std::setw(5) << j << std::setw(15) << U1[j + Ny * i] << std::setw(15) << V1[j + Ny * i] << std::endl;
        }
        outputFile << std::endl;
    }
    outputFile.close();
}