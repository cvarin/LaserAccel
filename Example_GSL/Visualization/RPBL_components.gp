#
# Script Gnuplot pour tracer les composantes de champ de
# faisceaux polarisés radialement.
#
################################################################################

datafile = './data/RPLB_components.dat'

set xlabel "Normalized radial coordinate (r_0/w_0)"
set ylabel "Normalized amplitude"

plot datafile using 1:2 with lines 3,\
     datafile using 1:3 with lines 1

################################################################################
pause -1 "Hit return to continue"
reset
