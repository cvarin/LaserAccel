#ifndef _INCLUDE_LaserAccel_hpp
#define _INCLUDE_LaserAccel_hpp

#include "MotionEQs.hpp"
#include "LaserBeams.hpp"

/****************** Structures ************************************************/
struct SolverParams
{
   double t,tf;     // Time window
   double h;        // Starting stepsize 
   double eps_abs;  // Absolute precision requested
};

/****************** Accessible functions **************************************/
void RPLB_Axial_Trajectory(double P, double wo, double T, double dzo, 
                            double lambda, double phio, double q, double m,
                             double z0, double v0, SolverParams &sp);

# endif // End of file
