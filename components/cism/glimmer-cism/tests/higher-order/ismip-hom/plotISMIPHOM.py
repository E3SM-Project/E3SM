#!/usr/bin/env python
# This script plots the results of running ISMIP-HOM experiments using Glimmer.
# Before running this script, run runISMIPHOM.py to generate the results.
# See the accompanying README file for more information.
# To see all command line options run: python plotResults.py --help
# Written February 26, 2010 by Glen Granzow at the University of Montana.

import os
import matplotlib.figure
from matplotlib import pyplot
from optparse   import OptionParser
from runISMIPHOM import appendToList, defaultExperiments, defaultSizes
from math       import sqrt

# Output flags
printNaNwarnings = False
savePlotInFile = True
plotType = '.png'

# Lists of model classifications
fullStokesModels = ['aas1','aas2','cma1','fpa2','ghg1','jvj1','mmr1','oga1','rhi1','rhi3','spr1','ssu1','yko1']
lmlaModels = ['ahu1','ahu2','bds1','fpa1','mbr1','rhi2','tpa1']

# Function to read data files
# Returns a list of tuples [(x0,|v0|),(x1,|v1|),...,(xm,|vm|)]
def read(filename,experiment):
  if experiment in ['b','d','e','f']:
#   Read two numbers, x and vx, from each line in the file
    n = 2 
  else:
#   Read four numbers, x, y, vx, and vy, from each line in the file
    n = 4
    X = 0 # X=0 => plot v vs x;  X=1 => plot v vs y 
    Y = 1 - X
    VX,VY = 2,3
    target = 0.25

  inputfile = open(filename)
  data = list()
  for line in inputfile:
    row = map(float,line.strip().split()[:n])
    if reduce(lambda a,b: a or b,[x != x for x in row]):
      if printNaNwarnings:
        print 'WARNING: NaN in file', filename, line,
    else:
      data.append(tuple(row))
  inputfile.close()

  if experiment in ['b','d','e','f']:
    return data
  else:
    if target in [row[Y] for row in data]:
#     Extract the points with the desired (target) y
      return [(row[X],sqrt(row[VX]**2+row[VY]**2)) for row in data if row[Y]==target]
    else:
#     Interpolate to the desired (target) y value
      below = 0
      above = 1
      for row in data:
        y = row[Y]
        if  below < y < target: below = y
        if target < y < above:  above = y
#     Extract the bracketing data
      dataA = [(row[X],sqrt(row[VX]**2+row[VY]**2)) for row in data if row[Y]==above]
      dataB = [(row[X],sqrt(row[VX]**2+row[VY]**2)) for row in data if row[Y]==below]
      if len(dataA) != len(dataB):
        print 'WARNING: unequal number of x values in file', filename
      for (a,b) in zip(dataA,dataB):
        if a[0]!=b[0]: 
          print 'WARNING: the x values are not the same in file', filename
#     Return the interpolated values
      alpha = (target-below)/(above-below)
      return [(a[0],alpha*a[1]+(1-alpha)*b[1]) for (a,b) in zip(dataA,dataB)]

if __name__ == '__main__':

# Parse the command line arguments
  parser = OptionParser()
  parser.add_option('-e','--exp',dest='experiments',type='string',action='callback',callback=appendToList,help='Which ISMIP-HOM experiments to run')
  parser.add_option('-s','--size',dest='sizes',type='string',action='callback',callback=appendToList,help='Which domain sizes to run')
  parser.add_option('-p','--prefix',dest='prefix',default='glm1',help='Prefix to use for model output files (defaults to glm1)')
  parser.add_option('-t','--title',dest='subtitle',help='Subtitle to place on the created graph')
  parser.add_option('-l','--lmla',dest='lmla',action='store_true',help='Compare to lmla models instead of all partial stokes models')
  options, args = parser.parse_args()
# If the user didn't specify a list of experiments or domain sizes, run the whole suite
  if options.experiments == None: options.experiments = defaultExperiments
  if options.sizes == None: options.sizes = defaultSizes

# Loop over the experiments requested on the command line
  for experiment in options.experiments:
    print 'ISMIP-HOM', experiment.upper()

#   Create the figure on which the plot axes will be placed
    figure = pyplot.figure(subplotpars=matplotlib.figure.SubplotParams(top=.85,bottom=.15))
    figure.text(0.5,0.92,'ISMIP-HOM Experiment '+experiment.upper(),horizontalalignment='center',size='large')
    if options.subtitle:
      figure.text(0.5,0.89,options.subtitle,horizontalalignment='center',size='small')
    figure.text(0.5,0.1,'Normalized X coordinate',horizontalalignment='center',size='small')
    figure.text(0.06,0.5,'Ice Speed (m/a)',rotation='vertical',verticalalignment='center')
#   Create the (three column) legend
    prop = matplotlib.font_manager.FontProperties(size='x-small')
    Line2D = matplotlib.lines.Line2D([],[],color=(0,0,0))
    figure.legend([Line2D],['Model Output'],loc=(0.1,0.05),prop=prop).draw_frame(False)
    Line2D.set_linestyle(':')
    Line2D.set_color((1,0,0))
    Patch = matplotlib.patches.Patch(edgecolor=None,facecolor=(1,0,0),alpha=0.25)
    figure.legend([Line2D,Patch],['Full Stokes Mean','Full Stokes Std. Dev.'],loc=(0.3,0.02),prop=prop).draw_frame(False)
    Line2D.set_color((0,0,1))
    Patch.set_facecolor((0,0,1))
    figure.legend([Line2D,Patch],['First Order Mean','First Order Std. Dev.'],loc=(0.55,0.02),prop=prop).draw_frame(False)

#   Loop over the sizes requested on the command line
    for i, size in enumerate(map(int,options.sizes)):

      if experiment == 'f': 
        if size != 100 or len(options.sizes) > 1:
          print 'NOTE: Experiment f uses a domain size of 100 km only'
        size = 100

#     Create the plot axes for this domain size
      if len(options.sizes) == 1:
        axes = figure.add_subplot(111)
      else:
        axes = figure.add_subplot(2,3,i+1)
        for tick in axes.xaxis.get_major_ticks():
          tick.label1.set_fontsize('xx-small')
        for tick in axes.yaxis.get_major_ticks():
          tick.label1.set_fontsize('xx-small')
      axes.set_title('%d km' % size, size='medium')

#     Get the Glimmer output data
      filename = os.path.join('output',options.prefix+experiment+'%03d'%size+'.txt')
      glimmerData = read(filename,experiment)
#     The Glimmer data is on a staggered grid;
#     Interpolate to obtain the value at x=0 and x=1
#     using periodic boundary conditions
#      v = (glimmerData[0][1] + glimmerData[-1][1])/2
#      glimmerData = [(0.0,v)]+glimmerData+[(1.0,v)]

#     Plot the Glimmer data
      axes.plot([row[0] for row in glimmerData],
                [row[1] for row in glimmerData],color='black')

#     Get the data from other models for comparison
      firstOrder = 0
      fullStokes = 1
      count = [0,0]
      sumV  = [[0.0 for v in glimmerData],[0.0 for v in glimmerData]]
      sumV2 = [[0.0 for v in glimmerData],[0.0 for v in glimmerData]]
      for (path,directories,filenames) in os.walk('ismip_all'):
        for filename in filenames:
          modelName = filename[0:4]
          modelExperiment = filename[4]
          modelSize = filename[5:8]
          if (modelExperiment != experiment) or (int(modelSize) != size) \
            or (options.lmla and not modelName in lmlaModels + fullStokesModels):
              continue
          data = read(os.path.join(path,filename),experiment)
          if modelName in fullStokesModels:
            index = fullStokes
          else:
            index = firstOrder
          count[index] += 1
#         Interpolate onto the x values from the Glimmer model run
          for (i,target) in enumerate([row[0] for row in glimmerData]):
            below = -999
            above =  999
            for (j,x) in enumerate([row[0] for row in data]):
              if  below <  x <= target: b,below = j,x
              if target <= x <  above:  a,above = j,x
            if above == below:
              v = data[a][1]
            else:
              if below == -999: # Use the periodic boundary condition at x = 0
                xBelow = data[-1][0] - 1
                vBelow = data[-1][1]
              else:
                xBelow,vBelow = data[b]
              if above ==  999: # Use the periodic boundary condition at x = 1
                xAbove = data[0][0] + 1
                vAbove = data[0][1]
              else:
                xAbove,vAbove = data[a]
              if xAbove == xBelow:
                print 'Surprise!',above,below,xAbove,xBelow,vAbove,vBelow
                v = (vAbove+vBelow)/2
              else:
                alpha = (target-xBelow)/(xAbove-xBelow)
                v = alpha*vAbove + (1-alpha)*vBelow
            sumV [index][i] += v
            sumV2[index][i] += v*v

#     Calculate statistics of the other model results
      if sum(count) == 0:
        print 'To compare with other models you need to download the ISMIP-HOM results from: http://www.the-cryosphere.net/2/95/2008/tc-2-95-2008-supplement.zipand unzip to a directory named ismip_all.  The ismip_all directory must be in the directory from which you are running this script.'
      else:
#       Find the mean and standard deviation of the velocities at each x
        for index in (firstOrder,fullStokes):
          if count[index] == 0:
            continue
          mean = list()
          standardDeviation = list()
          for i in range(len(glimmerData)):
            mean.append(sumV[index][i]/count[index])
            standardDeviation.append(sqrt(sumV2[index][i]/count[index]-mean[-1]**2))

#         Plot the mean using a dotted line
          color = (index,0,1-index) # blue for first order (index=0); red for full Stokes (index=1)
          x = [row[0] for row in glimmerData]
          axes.plot(x,mean,':',color=color)

#         Plot a filled polygon showing the mean plus and minus one standard deviation
          meanMinusSD = [m-sd for (m,sd) in zip(mean,standardDeviation)]
          meanPlusSD  = [m+sd for (m,sd) in zip(mean,standardDeviation)]
          x = x + list(reversed(x))
          y = meanPlusSD + list(reversed(meanMinusSD))
          axes.fill(x,y,facecolor=color,edgecolor=color,alpha=0.25)

          if index == firstOrder:
#           Calculate some statistics comparing the Glimmer data with the other models
            error = [abs(v-m)/m for (v,m) in zip([row[1] for row in glimmerData],mean)]
            maximum = max(error)
            position = glimmerData[error.index(maximum)][0]
            total   = sum([e for e in error])
            compare = sum([(s/m) for (s,m) in zip(standardDeviation,mean)])
            n = len(glimmerData)
            print '\t'.join([str(size)+' km',str(total/n),str(compare/n),str(position)])

    if savePlotInFile:
      filename = os.path.join('output','ISMIP-HOM-'+experiment.upper()+'-'+options.prefix+plotType)
      print 'Writing:', filename
      pyplot.savefig(filename)

#   Experiment f should be run for one size (100 km) only
    if experiment == 'f': break

  if not savePlotInFile:
    pyplot.show()
