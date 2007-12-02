#include<fstream>
#include<iostream>
#include<vector>

#include "PhysConsts.hpp"
#include "Gnuplot_i.hpp"
#include "SODE1D.hpp"

void WriteToGnuplot(Motion traj, int count)
{
    Gnuplot Plot;

    Plot.set_style("lines");
    Plot.set_xlabel("Temps");
    Plot.set_ylabel("Amplitude");
    Plot.plot_xy(traj.time,traj.position,"Position");
    Plot.plot_xy(traj.time,traj.velocity,"Vitesse");
    
    std::cout << "Press [ENTER] or [CTRL-C] to quit." << std::endl;
    std::getchar();
}

