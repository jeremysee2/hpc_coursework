#ifndef REACTIONDIFFUSION_H
#define REACTIONDIFFUSION_H

class ReactionDiffusion {
    private:
        int Nx, Ny, T, dx, dy, start, end, sz, rank, size;
        double dt, a, b, eps, mu1, mu2;
        double* U1;
        double* V1;
        double* U2;
        double* V2;
        double* U_output;
        double* V_output;
    public:
        /**
         * @brief Set the Parameters required for the problem, initialise the 2D matrix.
         * 
         * @param dt    Time-step to use.
         * @param T     Total integration time.
         * @param Nx    Number of grid points in x
         * @param Ny    Number of grid points in y
         * @param a     Value of parameter a
         * @param b     Value of parameter b
         * @param mu1   Value of parameter mu1
         * @param mu2   Value of parameter mu2
         * @param eps   Value of parameter epsilon
         */
        ReactionDiffusion(double dt, int T, int Nx, int Ny, double a,
                          double b, double mu1, double mu2, double eps,
                          double dx, double dy, int start, int end);
        /// Destructor, free memory
        ~ReactionDiffusion();

        /**
         * @brief Set the Initial Conditions of the problem
         */
        void SetInitialConditions();

        /**
         * @brief Perform the full solution.
         */
        void TimeIntegrate();

        /**
         * @brief Perform one time step of the integration.
         */
        void TimeIntegrateSingle();

        /**
         * @brief Perform integration at boundary.
         */
        void TimeIntegrateBC();

        /**
         * @brief Getter method to output entire grid.
         */
        void writeOutput();
};

#endif