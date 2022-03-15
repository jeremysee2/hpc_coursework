#ifndef REACTIONDIFFUSION_H
#define REACTIONDIFFUSION_H

class ReactionDiffusion {
    private:
        int Nx, Ny, T, dx, dy;
        double dt, a, b, eps, mu1, mu2;
        double* U1;
        double* V1;
        double* U2;
        double* V2;
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
                          double dx = 1, double dy = 1);
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
         * @brief Getter method to output entire grid.
         */
        void writeOutput();
};

#endif