#include <iostream>
#include "ReactionDiffusion.h"
#include "cblas.h"

// #define F77NAME(x) x##_
// extern "C" {

// }

ReactionDiffusion::ReactionDiffusion (double dt, int T, int Nx, int Ny, double a,
                                      double b, double mu1, double mu2, double eps) {
    ReactionDiffusion::dt   = dt;
    ReactionDiffusion::T    = T;
    ReactionDiffusion::Nx   = Nx;
    ReactionDiffusion::Ny   = Ny;
    ReactionDiffusion::a    = a;
    ReactionDiffusion::b    = b;
    ReactionDiffusion::mu1  = mu1;
    ReactionDiffusion::mu2  = mu2;
    ReactionDiffusion::eps  = eps;
}

ReactionDiffusion::~ReactionDiffusion() {
    delete ReactionDiffusion::matrixU;
    delete ReactionDiffusion::matrixV;
}

int ReactionDiffusion::SetInitialConditions() {
    return 0;
}

int ReactionDiffusion::SetParameters() {
    // Allocate memory for the matrix
    ReactionDiffusion::matrixU = new double[ReactionDiffusion::Nx*ReactionDiffusion::Ny];
    ReactionDiffusion::matrixV = new double[ReactionDiffusion::Nx*ReactionDiffusion::Ny];
    return 0;
}

int ReactionDiffusion::TimeIntegrate() {
    // Iterate over u matrix
    for (int i = 0; i < Nx; ++i) {
        for (int j = 0; j < Ny; ++j) {

        }
    }

    // Iterate over v matrix
    for (int i = 0; i < Nx; ++i) {
        for (int j = 0; j < Ny; ++j) {
            
        }
    }
    return 0;
}

double ReactionDiffusion::f1_(double u, double v) {
    double res = ReactionDiffusion::eps*u*(1.0-u)*
                 (u - (v+ReactionDiffusion::b)/ReactionDiffusion::a);
    return res;
}

double ReactionDiffusion::f2_(double u, double v) {
    double res = u*u*u - v;
    return res;
}