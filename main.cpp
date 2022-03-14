/**
 * @file main.cpp
 *
 * High-Performance Computing Coursework.
 * Solves PDE using Barkley model, parallelised code with OpenMP.
 * 
 */
#include <iostream>
#include <cmath>
#include <mpi.h>
#include <boost/program_options.hpp>

#include "ReactionDiffusion.h"

using namespace std;
namespace po = boost::program_options;

/**
 * @brief Solves PDE using Barkley model.
 */
int main(int argc, char* argv[]) {
    // MPI Setup
    int sz, rank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &sz);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if(sz < 1) {
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    // Default parameters
    int Nx = 101, Ny = 101, T = 100, dx = 1, dy = 1;
    double dt = 0.001, a = 0.75, b = 0.06, eps = 50.0, mu1 = 5.0, mu2 = 0.0;

    if (rank == 0) {
        // Task 1: Command line input parser
        // Declare the supported options.
        po::options_description desc("Allowed options");
        desc.add_options()
            ("dt",  po::value<double>(),  "Time-step to use.")
            ("T",   po::value<int>(),     "Total integration time.")
            ("Nx",  po::value<int>(),     "Number of grid points in x")
            ("Ny",  po::value<int>(),     "Number of grid points in y")
            ("a",   po::value<double>(),  "Value of parameter a")
            ("b",   po::value<double>(),  "Value of parameter b")
            ("mu1", po::value<double>(),  "Value of parameter mu1")
            ("mu2", po::value<double>(),  "Value of parameter mu2")
            ("eps", po::value<double>(),  "Value of parameter epsilon")
        ;

        // Parse parameters from command line
        po::variables_map vm;
        po::store(po::parse_command_line(argc, argv, desc), vm);
        po::notify(vm);    

        if (vm.count("dt"))      dt  = vm["dt"].as<double>();
        if (vm.count("T"))       T   = vm["T"].as<int>();
        if (vm.count("Nx"))      Nx  = vm["Nx"].as<int>();
        if (vm.count("Ny"))      Ny  = vm["Ny"].as<int>();
        if (vm.count("a"))       a   = vm["a"].as<double>();
        if (vm.count("b"))       b   = vm["b"].as<double>();
        if (vm.count("mu1"))     mu1 = vm["mu1"].as<double>();
        if (vm.count("mu2"))     mu2 = vm["mu2"].as<double>();
        if (vm.count("eps"))     eps = vm["eps"].as<double>();

        // Display current parameters
        cout << "Arguments: " << dt << ", " << T << ", " << Nx
            << ", " << Ny << ", " << a << ", " << b << ", " << mu1
            << ", " << mu2 << ", " << eps << endl;

        if (dt <= 0.0 || T < 0.0 || Nx < 0 || Ny < 0 || a < 0.0 ||
            b < 0.0 || mu1 < 0.0 || mu2 < 0.0 || eps < 0.0) {
                throw invalid_argument("Invalid parameters!");
        }
    }

    // MPI: Broadcast parameters to all processes
    MPI_Bcast(&dt,  1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&T,   1, MPI_INT,    0, MPI_COMM_WORLD);
    MPI_Bcast(&Nx,  1, MPI_INT,    0, MPI_COMM_WORLD);
    MPI_Bcast(&Ny,  1, MPI_INT,    0, MPI_COMM_WORLD);
    MPI_Bcast(&a,   1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&b,   1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&mu1, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&mu2, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&eps, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
    // MPI: Calculate partition of grid to be calculated by each process
    int partition   = Nx / sz + 1;
    int myStart     = rank*partition;
    int myEnd       = min(Nx, myStart+partition);

    cout << "Starting with " << sz << " processes from process " << rank << endl;

    // Store parameters in ReactionDiffusion Object
    ReactionDiffusion reactor(dt, T, Nx, Ny, a, b, mu1, mu2, eps, dx, dy, myStart, myEnd);

    // Set initial conditions
    reactor.SetInitialConditions();

    // Perform integration
    reactor.TimeIntegrate();

    // File output, values of u,v at each gridpoint x,y
    reactor.writeOutput();

    // Gracefully exit MPI environment
    MPI_Finalize();
    return 0;
}
