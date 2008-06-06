// Purpose =  reading a file and writing the parameters in the structure 'parameters'
//This Function is not used

#include <iostream>
#include <fstream>

#include "Structures.hpp"

using namespace std;

/******************************************************************************/
void GetParams(parameters *params)
{     
    FILE *ModelParams = fopen("ModelParams.ini", "r");
    
    if( ModelParams == NULL )
    {
        cout << "Error: Unable to open \'ModelParams.ini\'." << endl;
    }
    else fscanf(ModelParams, "mu : %lf", &params->mu);
    fclose(ModelParams);
}

/****************** End of file ***********************************************/
