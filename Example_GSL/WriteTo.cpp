#include<fstream>
#include<iostream>
#include<vector>

#include "PhysConsts.hpp"
#include "Gnuplot_i.hpp"
#include "Structures.hpp"

using namespace std;

/******************************************************************************/
void WriteToFile(vector<double> time,vector<double> pos,
     vector<double> vel,int count){      
     
     ofstream OutFile("Solution.dat", ios::out);
           
     for(int i=0 ; i<=count ; i++){
          OutFile << time[i] << "\t" 
                  << pos[i]  << "\t" 
                  << vel[i]  << "\n";
          }
          
     OutFile.close();
}

/******************************************************************************/
void WriteToGnuplot(void* trajectory,int count){
     
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
