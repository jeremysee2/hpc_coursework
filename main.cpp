/**
 * @file main.cpp
 *
 * High-Performance Computing Coursework.
 * Solves PDE using Barkley model, parallelised code with MPI.
 * 
 */
#include <iostream>
#include <iomanip>
#include <cmath>

#include "cblas.h"

using namespace std;
namespace po = boost::program_options;

// #define F77NAME(x) x##_
// extern "C" {

// }

/**
 * @brief Solves PDE using Barkley model.
 */
int main(int argc, char* argv[]) {
    // Task 1: Command line input parser
    // Declare the supported options.
    po::options_description desc("Allowed options");
    desc.add_options()
        ("dt",  po::value<double>(),  "Time-step to use.")
        ("T",   po::value<double>(),  "Total integration time.")
        ("Nx",  po::value<int>(),     "Number of grid points in x")
        ("Ny",  po::value<int>(),     "Number of grid points in y")
        ("a",   po::value<double>(),  "Value of parameter a")
        ("b",   po::value<double>(),  "Value of parameter b")
        ("mu1", po::value<double>(),  "Value of parameter mu1")
        ("mu2", po::value<double>(),  "Value of parameter mu2")
        ("eps", po::value<double>(),  "Value of parameter epsilon")
    ;

    // Default parameters
    int Nx = 101, Ny = 101, T = 100, dx = 1, dy = 1;
    double dt = 0.001, a = 0.75, b = 0.06, eps = 50.0, mu1 = 5.0, mu2 = 0.0;

    // Parse parameters from command line
    po::variables_map vm;
    po::store(po::parse_command_line(ac, av, desc), vm);
    po::notify(vm);    

    if (vm.count("dt"))      dt  = vm["dt"].as<double>();
    if (vm.count("T"))       T   = vm["Tt"].as<int>();
    if (vm.count("Nx"))      Nx  = vm["Nx"].as<int>();
    if (vm.count("Ny"))      Ny  = vm["Ny"].as<int>();
    if (vm.count("a"))       a   = vm["a"].as<double>();
    if (vm.count("b"))       b   = vm["b"].as<double>();
    if (vm.count("mu1"))     mu1 = vm["a"].as<double>();
    if (vm.count("mu2"))     mu2 = vm["b"].as<double>();
    if (vm.count("eps"))     eps = vm["a"].as<double>();

    // Display current parameters
    cout << "Arguments: " << dt << ", " << T << ", " << Nx
         << ", " << Ny << ", " << a << ", " << b << ", " << mu1
         << ", " << mu2 << ", " << eps << endl;

    // TODO: Calculate PDE solution

    // File output, values of u,v at each gridpoint x,y
    std::ofstream outputFile;
    outputFile.open("output.txt");
    outputFile.close();
}
