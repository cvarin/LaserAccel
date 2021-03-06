#
# Script Gnuplot pour tracer les composantes de champ de
# faisceaux polarisÚs radialement.
#
################################################################################

load "./Visualization/RPLB_components.path"

set term x11 font "arial,18" enhanced
set xlabel "Normalized radial coordinate ( r_0 / w_0 )"
set ylabel "Normalized amplitude"

plot datafile using 1:2 with lines 3,\
     datafile using 1:3 with lines 1

################################################################################
pause -1 "[Gnuplot message] Hit return to continue"
reset
