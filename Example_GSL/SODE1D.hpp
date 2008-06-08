#ifndef _INCLUDE_SODE1D_hpp
#define _INCLUDE_SODE1D_hpp

/****************** Structures ************************************************/
struct motion
{
   std::vector<double> time;
   std::vector<double> position;
   std::vector<double> velocity; 
};

/******************************************************************************/
struct GSL_ExampleParams
{
   double mu;
};

/****************** Accessible functions ************************************************************/
int GSL_Example(double t, const double y[], double dydt[], void* params);
void SODE1D(double t,double t1,double h,double eps_abs,
             GSL_ExampleParams params, double y[], motion &traj, int &count);

# endif // End of file
