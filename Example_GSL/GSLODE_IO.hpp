#ifndef _INCLUDE_WriteTo_hpp
#define _INCLUDE_WriteTo_hpp

using namespace std;

#include"SODE1D.hpp"

/******************************************************************************/
void Read_Input_File(GSL_ExampleParams *params, const char *filename);
void WriteToFile(vector<double> t, vector<double> p, 
                  vector<double> v, int count);
void WriteToGnuplot(void* trajectory, int count);

#endif // End of file
