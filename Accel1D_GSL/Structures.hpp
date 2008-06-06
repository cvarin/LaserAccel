#ifndef _INCLUDE_Structures_hpp
#define _INCLUDE_Structures_hpp

#include<vector>
             
struct motion
{
   std::vector<double> time;
   std::vector<double> position;
   std::vector<double> velocity; 
};
        
struct parameters
{
   double mu;
};

#endif // End of file
