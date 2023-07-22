#!/usr/bin/env python
from string import *
import os, getopt, sys
import matplotlib.pyplot as plot
import numpy as np

############# run on chrysalis with python 3, python 2 on anvil does not work

ddt=[0.2, 0.1, 0.05, 0.02, 0.01, 0.005, 0.002, 0.001]


seriesP=[ 0.003077509496527265,0.002971564287359151,0.00263778452855153,0.001755425382169638,0.001118177876390719,0.0006099010453942045,0.0002608441075285435,0.0001272050274434055 ]
seriesPhy=[ 0.003086108128760258,0.002972611990412168,0.002636010412023117,0.00175533930861668,0.001118370217551115,0.0006091090534282601,0.0002605342125618298,0.0001270081037904339 ]
seriesV=[ 0.00693267677466352,0.006481375979017358,0.005289501660810096,0.001935662595967612,0.0009215228580399363,0.0007046734778502718,0.0002752273001902675,0.00017930304565655 ]
seriesEdry=[ 0.006129476690641761,0.007223143634597555,0.00708357263368673,0.002762966755376245,0.001217748855692548,0.0004386992155546833,0.001146164225601199,0.00035028576305705 ]
seriesEstar=[ 0.006212394002788509,0.006912952061895793,0.006839840418866516,0.002742073244445472,0.001390152505598283,0.0004116095426320242,0.001410724895259956,0.001009662491319665 ]




#make hor line
hline=np.array(seriesV)-np.array(seriesV)+0.0024


seriesPerf=np.asarray(ddt)*seriesP[7]*500

MS=10
FS=14

plot.loglog(ddt, seriesV,    label="const volume (V)",       color="limegreen", linewidth=3,marker='^',markersize=MS )
plot.loglog(ddt, seriesPhy,  label="const pressure HY",      color="black", linewidth=3,marker='s',markersize=MS )
plot.loglog(ddt, seriesP,    label="const pressure (P)",     color="red", linewidth=3,marker='^',markersize=MS )
plot.loglog(ddt, seriesEdry, label="EAM cpdry", color="orange", linewidth=3,marker='v',markersize=MS )
plot.loglog(ddt, seriesEstar,label="EAM cpstar", color="blue", linewidth=3,marker='p',markersize=MS )
plot.loglog(ddt, seriesPerf, label="1st order scaling", color="gray", linewidth=1.5 )
plot.loglog(ddt, hline, label="uncertainty P vs V", color="gray", linewidth=1.5, linestyle='--')

plot.tick_params(labelsize=FS)
plot.legend(loc='lower left',fontsize=FS)
plot.xlabel("time step", fontsize=FS)
plot.ylabel("normalized error", fontsize=FS)
plot.gca().invert_xaxis()
plot.ylim([1e-5,2e-2])

plot.savefig('p-v-convergence.pdf')






