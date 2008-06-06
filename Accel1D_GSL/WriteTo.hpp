#ifndef _INCLUDE_WriteTo_hpp
#define _INCLUDE_WriteTo_hpp

using namespace std;

void WriteToFile(vector<double> time, vector<double> pos,
     vector<double> vel, int count);

void WriteToGnuplot(void* trajectory, int count);

#endif // End of file
