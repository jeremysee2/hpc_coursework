#include <iostream>
#include <iomanip>
#include <fstream>
#include "ReactionDiffusion.h"

#define RD ReactionDiffusion

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
}

RD::~ReactionDiffusion() {
    if(U1) delete U1;
    if(U2) delete U2;
    if(V1) delete V1;
    if(V2) delete V2;
}

int RD::SetInitialConditions() {
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
    return 0;
}

int RD::SetParameters() {
    // Allocate memory for the matrix
    U1 = new double[Nx*Ny]();
    U2 = new double[Nx*Ny]();
    V1 = new double[Nx*Ny]();
    V2 = new double[Nx*Ny]();
    return 0;
}

int RD::TimeIntegrate() {
    #pragma omp parallel for
    for (int k = 0; k<int(T/dt); ++k) {
        TimeIntegrateSingle();
    }
    return 0;
}

int RD::TimeIntegrateSingle() {
    // Iterate over u matrix
    double h = dx;
    for (int i = 0; i < Nx; ++i) {
        for (int j = 0; j < Ny; ++j) {
            U2[j+Ny*i] = dt * (mu1 /(h*h) * (
                                (i < Nx-1 ? U1[j+Ny*(i+1)] : 0) +
                                (i        ? U1[j+Ny*(i-1)] : 0) +
                                (j < Ny-1 ? U1[(j+1)+Ny*i] : 0) +
                                (j        ? U1[(j-1)+Ny*i] : 0) -
                                (4-(i==0)-(j==0)-(i==Nx-1)-(j==Ny-1))*U1[j+Ny*i]) +
                                f1_(U1[j+Ny*i], V1[j+Ny*i])) + U1[j+Ny*i];
        }
    }


    // Iterate over v matrix
    for (int i = 0; i < Nx; ++i) {
        for (int j = 0; j < Ny; ++j) {
            V2[j+Ny*i] = dt * (mu1 /(h*h) * (
                                (i < Nx-1 ? V1[j+Ny*(i+1)] : 0) +
                                (i        ? V1[j+Ny*(i-1)] : 0) +
                                (j < Ny-1 ? V1[(j+1)+Ny*i] : 0) +
                                (j        ? V1[(j-1)+Ny*i] : 0) -
                                (4-(i==0)-(j==0)-(i==Nx-1)-(j==Ny-1))*V1[j+Ny*i]) +
                                f2_(U1[j+Ny*i], V1[j+Ny*i])) + V1[j+Ny*i];
        }
    }

    // Save current time step for next iteration
    for (int i = 0; i < Nx; ++i) {
        for (int j = 0; j < Ny; ++j) {
            U1[j+Ny*i] = U2[j+Ny*i];
            V1[j+Ny*i] = V2[j+Ny*i];
        }
    }
    return 0;
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

double RD::f1_(double u, double v) {
    return double (eps*u*(1.0-u)* (u - (v+b)/a));
}

double RD::f2_(double u, double v) {
    return double (u*u*u - v);
}