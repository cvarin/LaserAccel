#ifndef _INCLUDE_SODE1D_hpp
#define _INCLUDE_SODE1D_hpp

#include"Structures.hpp"

int func (double t, const double y[], double dydt[], void* params);

void SODE1D(double t,double t1,double h,double eps_abs,
             parameters params, double y[], motion &traj, int &count);

# endif // End of file
