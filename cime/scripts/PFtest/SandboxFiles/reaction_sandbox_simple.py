import sys
import os
try:
  pflotran_dir = os.environ['PFLOTRAN_DIR']
except KeyError:
  try:
    pflotran_dir = '../../'
  except KeyError:
    print('PFLOTRAN_DIR must point to PFLOTRAN installation directory and be defined in system environment variables.')
    sys.exit(1)
sys.path.append(pflotran_dir + '/src/python')
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import math
import pflotran as pft

path = []
path.append('.')

scale_string = 'linear'
#scale_string = 'log'

aq_labels = []
aq_labels.append('Aaq')
aq_labels.append('Baq')
aq_labels.append('Caq')
aq_labels.append('Daq')
aq_labels.append('Eaq')
aq_labels.append('Faq')

aq_colors = []
aq_colors.append('blue')
aq_colors.append('green')
aq_colors.append('red')
aq_colors.append('cyan')
aq_colors.append('magenta')
aq_colors.append('y')

im_labels = []
im_labels.append('Xim')
im_labels.append('Yim')

im_colors = []
im_colors.append('darkorange')
im_colors.append('navy')

f = plt.figure(figsize=(16,6))
f.suptitle("Simple Reaction Sandbox",fontsize=16)

# concentration profiles
plt.subplot(1,2,2)
plt.title('Concentration Profile')

files = pft.get_tec_filenames('pflotran',[2])
filenames = pft.get_full_paths(path,files)

plt.xlabel('X [m]')
plt.ylabel('Concentration [M]')

plt.yscale(scale_string)

maxval = -1.e20
minval = 1.e20
for ifile in range(len(filenames)):
  columns = [4,5,6,7,8,9]
  for icol in range(len(columns)):
    data = pft.Dataset(filenames[ifile],1,columns[icol])
    ydata = data.get_array('y')
    maxval = max(maxval,np.amax(ydata))
    minval = min(minval,np.amin(ydata))
    plt.plot(data.get_array('x'),ydata,
             label=aq_labels[icol],c=aq_colors[icol])
plt.ylim(0.95*minval,1.05*maxval)

#'best'         : 0, (only implemented for axis legends)
#'upper right'  : 1,
#'upper left'   : 2,
#'lower left'   : 3,
#'lower right'  : 4,
#'right'        : 5,
#'center left'  : 6,
#'center right' : 7,
#'lower center' : 8,
#'upper center' : 9,
#'center'       : 10,
plt.legend(loc=1)
# xx-small, x-small, small, medium, large, x-large, xx-large, 12, 14
plt.setp(plt.gca().get_legend().get_texts(),fontsize='small')
#      plt.setp(plt.gca().get_legend().get_texts(),linespacing=0.)
plt.setp(plt.gca().get_legend().get_frame().set_fill(False))
plt.setp(plt.gca().get_legend().draw_frame(False))
#        plt.gca().yaxis.get_major_formatter().set_powerlimits((-1,1))

plt.twinx()
plt.ylabel('Concentration [mol/m^3]')
#plt.ylim(0.,1.1e-5)
plt.yscale(scale_string)
for ifile in range(len(filenames)):
  columns = [10,11]
  for icol in range(len(columns)):
    data = pft.Dataset(filenames[ifile],1,columns[icol])
    ydata = data.get_array('y')
    maxval = max(maxval,np.amax(ydata))
    minval = min(minval,np.amin(ydata))
    plt.plot(data.get_array('x'),data.get_array('y'),
             label=im_labels[icol],c=im_colors[icol])
plt.ylim(0.95*minval,1.05*maxval)

plt.legend(loc=4)
# xx-small, x-small, small, medium, large, x-large, xx-large, 12, 14
plt.setp(plt.gca().get_legend().get_texts(),fontsize='small')
#      plt.setp(plt.gca().get_legend().get_texts(),linespacing=0.)
plt.setp(plt.gca().get_legend().get_frame().set_fill(False))
plt.setp(plt.gca().get_legend().draw_frame(False))
#        plt.gca().yaxis.get_major_formatter().set_powerlimits((-1,1))

# concentration breakthrough at observation point
plt.subplot(1,2,1)
plt.title('Concentration Breakthrough')

files = []
files.append('pflotran-obs-0.tec')
filenames = pft.get_full_paths(path,files)

plt.xlabel('Time [y]')
plt.ylabel('Concentration [M]')

#plt.xlim(0.,1.)
#plt.ylim(0.,1.)
#plt.grid(True)
plt.yscale(scale_string)

for ifile in range(len(filenames)):
  columns = [2,3,4,5,6,7]
  for icol in range(len(columns)):
    data = pft.Dataset(filenames[ifile],1,columns[icol])
    ydata = data.get_array('y')
    maxval = max(maxval,np.amax(ydata))
    minval = min(minval,np.amin(ydata))
    plt.plot(data.get_array('x'),data.get_array('y'),
             label=aq_labels[icol],c=aq_colors[icol])
plt.ylim(0.95*minval,1.05*maxval)

#'best'         : 0, (only implemented for axis legends)
#'upper right'  : 1,
#'upper left'   : 2,
#'lower left'   : 3,
#'lower right'  : 4,
#'right'        : 5,
#'center left'  : 6,
#'center right' : 7,
#'lower center' : 8,
#'upper center' : 9,
#'center'       : 10,
plt.legend(loc=4)
# xx-small, x-small, small, medium, large, x-large, xx-large, 12, 14
plt.setp(plt.gca().get_legend().get_texts(),fontsize='small')
#      plt.setp(plt.gca().get_legend().get_texts(),linespacing=0.)
plt.setp(plt.gca().get_legend().get_frame().set_fill(False))
plt.setp(plt.gca().get_legend().draw_frame(False))
#        plt.gca().yaxis.get_major_formatter().set_powerlimits((-1,1))

f.subplots_adjust(hspace=0.2,wspace=0.40,
                  bottom=.12,top=.85,
                  left=.08,right=.92)

plt.show()
