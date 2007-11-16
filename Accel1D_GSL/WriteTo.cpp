#include<fstream>
#include<iostream>
#include<vector>

#include "PhysConsts.hpp"
#include "Gnuplot_i.hpp"
#include "Structures.hpp"


void WriteToGnuplot(void* trajectory,int count){
     
     Gnuplot Plot;
     Motion* PtoTraj = (Motion*)trajectory;//Typecast from void* to motion*
     
     std::vector<double> TempTime = PtoTraj->time;
     std::vector<double> TempPosition = PtoTraj->position;
     std::vector<double> TempVelocity = PtoTraj->velocity;

     Plot.set_style("lines");
     Plot.set_xlabel("Temps");
     Plot.set_ylabel("Amplitude");
     Plot.plot_xy(TempTime,TempPosition,"Position");
     Plot.plot_xy(TempTime,TempVelocity,"Vitesse");
    
     std::cout << "Press [ENTER] or [CTRL-C] to quit." << std::endl;
     std::getchar();
     
     }

