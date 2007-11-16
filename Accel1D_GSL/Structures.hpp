#ifndef INC_Structures_hpp
#define INC_Structures_hpp

#include <vector>

             
struct Motion{
       std::vector<double> time;
       std::vector<double> position;
       std::vector<double> velocity; 
       };
        
struct Parameters{
             double mu;
             };
             
#endif // INC_Structures_hpp

