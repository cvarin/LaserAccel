#include<fstream>
#include<iostream>
#include<vector>

#include "Constants.hpp"
#include "Gnuplot_i.hpp"
#include "SODE1D.hpp"

using namespace std;

/******************************************************************************/
/****************** External functions ****************************************/
/******************************************************************************/
void Read_Input_File(GSL_ExampleParams *params)
{     
    FILE *ModelParams = fopen("ModelParams.ini", "r");
    
    if( ModelParams == NULL )
    {
        cout << "Error: Unable to open \'ModelParams.ini\'." << endl;
    }
    else fscanf(ModelParams, "mu : %lf", &params->mu);
    fclose(ModelParams);
}

/******************************************************************************/
void WriteToFile(vector<double> t, vector<double> p, 
                  vector<double> v, int count)
{      
     ofstream OutFile("Solution.dat", ios::out);    
     for(int i=0 ; i<=count ; i++)
     {
        OutFile << t[i]  << "\t" 
                << p[i]  << "\t" 
                << v[i]  << "\n";
     }
     OutFile.close();
}

/******************************************************************************/
void WriteToGnuplot(void* trajectory, int count)
{     
     Gnuplot Plot;
     motion* PtoTraj = (motion*)trajectory;//Typecast from void* to motion*
     
     vector<double> TempTime = PtoTraj->time;
     vector<double> TempPosition = PtoTraj->position;
     vector<double> TempVelocity = PtoTraj->velocity;

     Plot.set_style("lines");
     Plot.set_xlabel("Temps");
     Plot.set_ylabel("Amplitude");
     Plot.plot_xy(TempTime,TempPosition,"Position");
     Plot.plot_xy(TempTime,TempVelocity,"Vitesse");
    
     cout << "Press [ENTER] or [CTRL-C] to quit." << endl;
     getchar();    
}

/****************** End of file ***********************************************/
