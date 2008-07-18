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
void RPLB_Axial_MultiParticle(double P, double wo, double T, double zf, 
                               double dzo, double lambda, double phio, double q, 
                                double m, int N, double z0, double cw,
                                 SolverParams sp);
void RPLB_Axial_Trajectory(double P, double wo, double T, double zf, double dzo, 
                            double lambda, double phio, double q, double m,
                             double z0, double v0, SolverParams sp);
void RPLB_3D_Trajectory(double P, double wo, double T, double zf, double dzo, 
                         double lambda, double phio, double q, double m,
                          double r0, double vr0, double z0, double vz0, 
                           SolverParams sp);
double *RPLB_Phase_Scan(double P, double wo, double T, double zf, double dzo, 
                         double lambda, double q, double m, double z0, 
                          double v0, SolverParams sp, int N);

# endif // End of file
