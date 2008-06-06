// Purpose =  reading a file and writing the parameters in the structure 'parameters'
//This Function is not used

#include <iostream>
#include <fstream>

#include "Structures.hpp"

using namespace std;

void GetParams(void* params){
     
     parameters *PtoParams = (parameters *)params; // Typecast from void* to parameters*
     FILE *ModelParams = fopen("ModelParams.ini", "r");
     
     if( ModelParams == NULL )
     {
     cout << "Error: Unable to open \'ModelParams.ini\'." << endl;
     }
     
     else
     {
      fscanf(ModelParams, "mu : %lf", &(PtoParams->mu));

      fclose(ModelParams);
      }
     
     }
