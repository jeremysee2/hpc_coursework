#ifndef REACTIONDIFFUSION_H
#define REACTIONDIFFUSION_H

class ReactionDiffusion {
    private:
        double f1_(double u, double v); /// Calculate parameter f1
        double f2_(double u, double v); /// Calculate parameter f2

        int Nx, Ny, T, dx, dy;
        double dt, a, b, eps, mu1, mu2;
        double* matrixU;
        double* matrixV;
    public:
        /// Constructor
        ReactionDiffusion(double dt, int T, int Nx, int Ny, double a,
                          double b, double mu1, double mu2, double eps);
        /// Destructor, free memory
        ~ReactionDiffusion();

        /**
         * @brief Set the Parameters required for the problem, initialise the 2D matrix.
         * 
         * @return int Error code.
         */
        int SetParameters();

        /**
         * @brief Set the Initial Conditions of the problem
         * 
         * @return int Error code.
         */
        int SetInitialConditions();

        /**
         * @brief Perform the full solution.
         * 
         * @return int Error code.
         */
        int TimeIntegrate();
};

#endif