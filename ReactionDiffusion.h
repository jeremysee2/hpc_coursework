#ifndef REACTIONDIFFUSION_H
#define REACTIONDIFFUSION_H

class ReactionDiffusion {
    private:

    public:
        // Constructor
        ReactionDiffusion();
        // Destructor
        ~ReactionDiffusion();
        int SetParameters();
        int SetInitialConditions();
        int TimeIntegrate();
};

#endif