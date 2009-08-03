#! /usr/bin/env python
# -*- coding: utf-8 -*-

######### Modules #########################################
import matplotlib
matplotlib.use('WXAgg')
import numpy
from matplotlib.pylab import *

#from matplotlib import rc
#rc('text', usetex=True)

######### Paramètres de simulation ########################
filename = "../output/balayage_puissance.dat"
W0 = 30.0
spotsize = 5.0e-6
#P_seuil = 0.1*(spotsize/0.8e-6)**4*(1.0/W0/0.511)
#P_seuil = 20
#Wmax0   = 1.17
P_seuil = 30
Wmax0   = 1.165

######### Les données sont lues dans le fichier ###########
print '\r Reading %s' %filename
data = load(filename);

P = data[:,0]
I = data[:,1]
Wmax  = data[:,2]
Wmin  = data[:,4]
Wmax_analytique = Wmax0*sqrt(240/pi*(P - P_seuil))

######### Crée le graphique ###############################
figure(figsize=(8,6))
get_current_fig_manager().window.SetPosition([300,200])
grid(True)
plot(P,Wmax,label=r'Maximum')
plot(P,Wmin,label=r'Minimum')
plot(P,Wmax_analytique,label=r'1D Model (%.2f,' %P_seuil + '%.2f)' %Wmax0)
axhline(y=-W0, linewidth=1.0, color='black',label='Initial')
xlabel(r'Laser power (TW)')
ylabel(r'Net energy gain $W - W_0$ (MeV)')
legend(loc = 'center right')

###########################################################
show()

######### End of script ###################################
