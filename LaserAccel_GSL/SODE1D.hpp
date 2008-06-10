#ifndef _INCLUDE_SODE1D_hpp
#define _INCLUDE_SODE1D_hpp

#include <vector>

/****************** Structures ************************************************/
struct motion1D
{
   std::vector<double> time;
   std::vector<double> position;
   std::vector<double> velocity; 
};

/****************** Accessible functions ************************************************************/
int GSL_Example(double t, const double y[], double dydt[], void* params);
void SODE1D(double t, double t1, double h, 
             double eps_abs, double mu, double y[]);

# endif // End of file
