#!/usr/bin/env python

"""
This script plots the results of running ISMIP-HOM experiments using Glimmer.
Before running this script, run runISMIP_HOM.py to generate the results.
runISMIP_HOM.py generates a set of output files that will follow the pattern:
ismip-hom-?[-MOD].RESO.[pPROC.]out.nc, where `?` is a POSIX metacharacter,
[-MOD] is an optional user specified file modifier, RESO is the size of the
experiment, and [pPROC.] is the optional number of processors used to run the model
(prefixed by a p). This script will plot the set of experiments that follow the
above pattern. If more than one set exists, it will plot the last modified set,
unless a specific set is specified. See the accompanying README file for more
information.
"""

# To see all command line args.run: python plotResults.py --help
# Written February 26, 2010 by Glen Granzow at the University of Montana.
# Reconfigured by Joseph H Kennedy at ORNL on August 7, 2015 to work with the regression testing
#     NOTE: Did not adjust inner workings except where needed. This could
#     seriously use a clean-up (probably an entire rewrite).

import os
import sys
import glob
import argparse
import numpy as np

from math import sqrt, sin, cos, pi
from netCDF import *
from optparse import OptionParser

import matplotlib.figure
from matplotlib import pyplot

# Output flags
printNaNwarnings = False
savePlotInFile = True
plotType = '.png'


# The command line options
# ------------------------
parser = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

def unsigned_int(x):
    """
    Allows argparse to understand unsigned integers. 
    """
    x = int(x)
    if x < 1:
        raise argparse.ArgumentTypeError("This argument is an unsigned int type! Should be an integer greater than zero.")
    return x

parser.add_argument('-o', '--output-dir', default='./output',
        help="Directory containing the tests output files. Warning: if there is a " \
        +"path passed via the -f/--out-file option, this argument will be " \
        +"ignored.")

parser.add_argument('-f', '--output-file',
        help="One of the tests output files from the set you would like to plot. If a path is " \
        +"passed via this option, the -o/--output-dir option will be ignored.")

#FIXME: do we need a valid range here?
parser.add_argument('--sizes', nargs='*', type=unsigned_int, metavar='KM', 
        help="List (separated by spaces) the domain sizes to plot. "
            +"Note: sizes will only be applied to experiments a though e. Experiment f has only one size (100 km). ")

# Plot specific options
parser.add_argument('-s', '--subtitle',
        help="Subtitle to place on the created graph.")
parser.add_argument('-a', '--all-ps', action='store_true',
        help="Compare to all partial stokes models instead of just lmla models.")


# ===========================================================
# Define some variables and functions used in the main script
# ===========================================================

def get_file_pattern():
    """
    Get the search pattern for a single set of test output files. 
    """
    if args.output_file:
        out_d, out_f = os.path.split(args.output_file)
        if out_d:
            args.output_dir = out_d
            args.output_file = out_f
   
        if not out_f.endswith('.out.nc'):
            print('\nERROR: You must select an output netCDF file (*.out.nc) with the -f/--output-file option!')
            print('You selected: '+os.path.join(args.output_dir, args.output_file)+'\n')
            sys.exit(1)

    else:
        outpath = os.path.join(args.output_dir, '*.out.nc')
        matching = glob.glob(outpath)
        if len(matching) > 0:
            newest = max(matching, key=os.path.getmtime)
            
        else:
            print('\nERROR: NO *.out.nc FILES DETECTED!')
            print(  '==================================')
            print(  'Either specify a location to look for the test output')
            print(  'files with the -o/--output-dir option, or one of test output')
            print(  'files from the set you would like to plot with the -f/--output-file option.\n')
            sys.exit(1)

        args.output_file = os.path.basename(newest)
        
    experiment, mod, size, processors = split_file_name(args.output_file)

    pattern = os.path.join(args.output_dir, 'ismip-hom-?'+mod+'.????'+processors+'.out.nc')

    print('\nFinding results to plot. Using search pattern: '+pattern)
    print('   NOTE: there may be another set or sets of output files that do not match this search pattern. ' \
          +'Use the `-f/--output-file` to specify any output file within a set to plot that set. See the ' \
          +'README.md for more information.\n' )
    return pattern


def split_file_name(file_name):
    """
    Get the experiment, filename modifiers, size, and number of processors out of an out.nc filename.
    """
    experiment = ''
    mod = ''
    size = ''
    processors = ''

    file_details = file_name.replace('.out.nc','') .split('.')
    if len(file_details) > 2:
        processors = '.'+file_details[2]
    size = file_details[1]

    file_base = file_details[0].split('-')
    experiment = file_base[2]
    if len(file_base) > 3:
        mod = '-'+'-'.join(file_base[3:])

    return (experiment, mod, size, processors)

def get_experiments_and_sizes(pattern):
    """
    Get all the experiments and sizes for a single set of test output files. 
    """
    matching = glob.glob(pattern)
    
    sizes = set()
    experiments = set()
    
    # so we can ignore the size for f
    matching_f_tests = glob.glob(pattern.replace('-?','-f'))
    if matching_f_tests:
        experiments.add('f')
        matching = list(set(matching) - set(matching_f_tests))
    
    for match in matching:
        mt = os.path.basename(match)
        ex, mod, sz, pr = split_file_name(mt)
        sizes.add(sz)
        experiments.add(ex)
  
    return (list(experiments), list(sizes))


# Lists of model classifications
# Note: 'aas1' is skipped below -- see line ~165
fullStokesModels = ['aas1','aas2','cma1','fpa2','ghg1','jvj1','mmr1','oga1','rhi1','rhi3','spr1','ssu1','yko1']
lmlaModels = ['ahu1','ahu2','bds1','fpa1','mbr1','rhi2','tpa1']

# Classification used in Pattyn et al. 2008 paper, but not used by this script.  Listed here for reference only.
#nonFSModels = ['ahu1','ahu2','bds1','cma2','dpo1','fpa1','lpe1','mbr1','mtk1','oso1','rhi2','rhi4','rhi5','tpa1']

# Function to read data files
# Returns a list of tuples [(x0,|v0|),(x1,|v1|),...,(xm,|vm|)]
def read(filename,experiment):
    if experiment in ['b','d','e']:
        #   Read two numbers, x and vx, from each line in the file
        n = 2 
        VX = 1  # VX is in 1 position for the 2d tests
        X = 0
    elif experiment in ['f',]:
        # Read three numbers, x, y and vx,vy, from each line in the file
        n = 5
        VX,VY = 3,4  # VX is in position 3 for the F test
        X = 0
        Y = 1
        target = 0.0  # y-position in km to get the profile from
    elif experiment in ['f-elevation',]:  # special case for returning the elevation data from the 'f' experiment
        n = 3
        X,Y = 0,1
        # surface elevation is in position 2
        ELEV = 2
        target = 0.0  # y-position in km to get the profile from
    else:
        # Read four numbers, x, y, vx, and vy, from each line in the file
        n = 4
        X = 0 # X=0 => plot v vs x;  X=1 => plot v vs y 
        Y = 1 - X
        VX,VY = 2,3
        target = 0.25  # desired dimensionless y-position to use for the profile

    with  open(filename) as inputfile:
        data = list()
        for line in inputfile:
            row = map(float,line.strip().split()[:n])
            if reduce(lambda a,b: a or b,[x != x for x in row]):
                if printNaNwarnings:
                    print 'WARNING: NaN in file', filename, line,
            else:
                data.append(tuple(row))

    if experiment in ['b','d','e']:
        return [( row[X], abs(row[VX]) ) for row in data]  # use the absolute value of the x-velocity
    else:
        if target in [row[Y] for row in data]:
            # Extract the points with the desired (target) y
            if experiment in ['f-elevation',]:
                dataA = [(row[X], row[ELEV]) for row in data if row[Y]==target]
            else:
                dataA = [(row[X],sqrt(row[VX]**2+row[VY]**2)) for row in data if row[Y]==target]
            return dataA
        else:
            #print "Note: Plotting model output data along profile at dimensionless y-position "+str(target)+" requires interpolation.  "+filename
            # Interpolate to the desired (target) y value
            below = -100000.0
            above =  100000.0
            for row in data:
                y = row[Y]
                if  below < y < target: below = y
                if target < y < above:  above = y
            #print 'got above=',above,' below=',below
            # Extract the bracketing data
            if experiment in ['f-elevation',]:
                dataA = [(row[X],row[ELEV]) for row in data if row[Y]==above]
                dataB = [(row[X],row[ELEV]) for row in data if row[Y]==below]
            else:
                dataA = [(row[X],sqrt(row[VX]**2+row[VY]**2)) for row in data if row[Y]==above]
                dataB = [(row[X],sqrt(row[VX]**2+row[VY]**2)) for row in data if row[Y]==below]
            #if len(dataA) != len(dataB):
            #  print 'WARNING: unequal number of x values in file', filename
            for (a,b) in zip(dataA,dataB):
                if a[0]!=b[0]: 
                    print 'WARNING: the x values are not the same in file', filename
            # Return the interpolated values
            alpha = (target-below)/(above-below)
            return [(a[0],alpha*a[1]+(1-alpha)*b[1]) for (a,b) in zip(dataA,dataB)]



# =========================
# Actual script starts here
# =========================
def main():
    """
    Plot the ISMIP-HOM test results.
    """
    
    if args.all_ps:
         nonFSmodelType='All Partial Stokes'
    else:
         nonFSmodelType='First Order'
    print 'NOTE: The category being used for models approximating Full Stokes is: '+nonFSmodelType
    print 'For more information, see details of option -a by invoking:   python plotISMIP_HOM.py --help \n'


    # =========================================================
    #  First generate the standard ascii files defined by ISMIP-HOM
    #  (This used to be in the run script, but is more useful and natural here.)
    #  Note: It is not strictly necessary to generate these files
    #  just to plot the results, but this approach is left in 
    #  to use existing code, and on the off-chance it is ever
    #  necessary to generate these standard ascii ISMIP-HOM files...
    # =========================================================

    pattern = get_file_pattern() 

    experiments, sizes = get_experiments_and_sizes(pattern)

    if args.sizes and (set(args.sizes) <= set(map(int,sizes))):
        sizes = args.sizes
    elif args.sizes:
        print("ERROR: No output files found for that size!")
        print("   Run `ls "+pattern.replace('-?','-[a-e]')+"`\non the command line to see the set of output files. (Note, test f's size is ignored.)")
        sys.exit(1)

    # sort experiments and sizes
    sizes.sort()
    experiments.sort()

    # Loop over the experiments
    for experiment in experiments:

        # Loop over the sizes
        for size in map(int,sizes):
            try:
                # Extract the output data for comparison to the other models

                # NOTE: The script now assumes that uvel_extend and vvel_extend are ALWAYS present.
                # Those fields contain the ice velocity computed at the upper right corner of each grid cell.
                # They appear to be on the x1,y1 grid in their metadata but are actually on the x0,y0 grid.
                # The additional row/column include the first halo value past ewn/nsn.
                # That value is valid at both 0.0 and 1.0 on the non-dimensional coordinate system.
                # Matrix manipulations for each test case below are done to create a larger matrix that goes from 0.0 to 1.0, inclusive.
                # NOTE: The cases below are only writing [x,y,u,v] to the text file.  This is the minimum needed to compare to other models.
                # In the future, the other additional fields specified in section 4 of http://homepages.ulb.ac.be/~fpattyn/ismip/ismiphom.pdf
                # can be added.  wvel and the stresses are on the x1,y1 grid, so they would need to be interpolated to the x0,y0 grid
                # since we are using that as the coordinate system in the text files.

                # Open the netCDF file that was written by CISM

                # Standard filename format used by both scripts
                if experiment == 'f': 
                    size = 100
                
                out_file = pattern.replace('-?','-'+experiment).replace('????',str(size).zfill(4))
                netCDFfile = NetCDFFile(out_file,'r')
                if netCDF_module == 'Scientific.IO.NetCDF':
                    velscale = netCDFfile.variables['uvel_extend'].scale_factor
                else:
                    velscale = 1.0


                if experiment in ['f',]:
                    # Convert CISM output data to the rotated coord system used by the problem setup
                    alpha = -3.0 * pi/180  # defined in run script
                    if netCDF_module == 'Scientific.IO.NetCDF':
                        thkscale = netCDFfile.variables['thk'].scale_factor
                        wvelscale = netCDFfile.variables['wvel_ho'].scale_factor
                    else:
                        thkscale = 1.0
                        wvelscale = 1.0
                    
                    usurf = netCDFfile.variables['usurf'][-1,:,:] * thkscale  # get last time level
                    usurfStag = (usurf[1:,1:] + usurf[1:,:-1] + usurf[:-1,:-1] + usurf[:-1, :-1]) / 4.0
                    uvelS = netCDFfile.variables['uvel'][-1,0,:,:] * velscale  # top level of last time
                    vvelS = netCDFfile.variables['vvel'][-1,0,:,:] * velscale  # top level of last time
                    wvelS = netCDFfile.variables['wvel_ho'][-1,0,:,:] * wvelscale  # top level of last time
                    wvelStag = (wvelS[1:,1:] + wvelS[1:,:-1] + wvelS[:-1,:-1] + wvelS[:-1, :-1]) / 4.0
                    x0 = netCDFfile.variables['x0'][:]
                    y0 = netCDFfile.variables['y0'][:]
                    # calculate rotated xprime coordinates along the surface - xx, yy are used by code below to write the output
                    xx = x0 * cos(alpha) + (usurfStag[20,:]-7000.0) * sin(alpha)
                    xx = xx/1000.0 - 50.0
                    yy = y0/1000.0 - 50.0
                    # calculate rotated uvel/vvel at surface
                    uvelSprime =  uvelS[:,:] * cos(alpha) + wvelStag[:,:] * sin(alpha)
                    wvelSprime = -uvelS[:,:] * sin(alpha) + wvelStag[:,:] * cos(alpha)

                    nan = np.ones(uvelSprime.shape)*-999.0  # create a dummy matrix for uncalculated values.

                    # ===========================================
                    # optional bit of code to plot out vertical velocity on the rotated grid.
                    # This can be compared to Fig. 14b in the tc-2007-0019-sp3.pdf document
                    figure = pyplot.figure(subplotpars=matplotlib.figure.SubplotParams(top=.85,bottom=.15))
                    axes = figure.add_subplot(111)
                    pc = axes.pcolor(wvelSprime)
                    pyplot.colorbar(pc)
                    cntr=axes.contour(wvelSprime, [-0.5, -0.4, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4],colors='k')
                    axes.clabel(cntr)
                    pyplot.title('Compare to Fig. 14b of tc-2007-0019-sp3.pdf')
                    pyplot.draw()
                    pyplot.show()
                    # ===========================================

                else:  # all other tests
                    # Make x,y position arrays that can be used by all test cases.
                    #   Want x/y positions to include the periodic edge at both the beginning and end
                    xx = netCDFfile.variables['x0'][:]/(1000.0*float(size))
                    xx = np.concatenate(([0.0],xx,[1.0]))
                    yy = netCDFfile.variables['y0'][:]/(1000.0*float(size))
                    yy = np.concatenate(([0.0],yy,[1.0]))
                    if experiment in ('b','d'):
                        yy = (yy[len(yy)/2],)  # for the 2-d experiments, just use the middle y-index

                    # Figure out u,v since all experiments needs at least one of them (avoids duplicate code in each case below
                    #   Want to use last time level.  Most experiments should only have a single time level, but F may have many in the file.
                    #   Apparently some older versions of netCDF4 give an error when using the -1 dimension if the size is 1, hence this bit of seemingly unnecessary logic...
                    if netCDFfile.variables['uvel_extend'][:].shape[0] == 1:
                        t = 0
                    else:
                        t = -1
                    us = netCDFfile.variables['uvel_extend'][t,0,:,:] * velscale  # top level of last time
                    us = np.concatenate( (us[:,-1:], us), axis=1)  # copy the column at x=1.0 to x=0.0
                    us = np.concatenate( (us[-1:,:], us), axis=0)  # copy the row at y=1.0 to y=0.0
                    vs = netCDFfile.variables['vvel_extend'][t,0,:,:] * velscale  # top level of last time
                    vs = np.concatenate( (vs[:,-1:], vs), axis=1)  # copy the column at x=1.0 to x=0.0
                    vs = np.concatenate( (vs[-1:,:], vs), axis=0)  # copy the row at y=1.0 to y=0.0
                    ub = netCDFfile.variables['uvel_extend'][t,-1,:,:] * velscale  # bottom level of last time
                    ub = np.concatenate( (ub[:,-1:], ub), axis=1)  # copy the column at x=1.0 to x=0.0
                    ub = np.concatenate( (ub[-1:,:], ub), axis=0)  # copy the row at y=1.0 to y=0.0
                    vb = netCDFfile.variables['vvel_extend'][t,-1,:,:] * velscale  # bottom level of last time
                    vb = np.concatenate( (vb[:,-1:], vb), axis=1)  # copy the column at x=1.0 to x=0.0
                    vb = np.concatenate( (vb[-1:,:], vb), axis=0)  # copy the row at y=1.0 to y=0.0
                    #nan = ub*np.NaN  # create a dummy matrix for uncalculated values.
                    nan = np.ones(ub.shape)*-999.0  # create a dummy matrix for uncalculated values.

                # make arrays of the variables needed for each experiment
                # the extended-grid velocities have the periodic edge in the last x-position.  We also want it in the first x-position.
                # After building the 2-d array as needed for each variable from the raw file data, then build a list called 'data'.
                if experiment == 'a':
                    #  This is supposed to be: [('uvel',0),('vvel',0),('wvel',0),('tau_xz',-1),('tau_yz',-1),[deltap]]
                    data = (us, vs, nan, nan, nan, nan)
                elif experiment == 'b':
                    #  This is supposed to be: uvel(0), wvel(0), tau_xz(-1), deltaP
                    data = (us, nan, nan, nan)
                elif experiment == 'c':
                    #  This is supposed to be: [uvel',0),('vvel',0),('wvel',0),('uvel',-1),('vvel',-1),('tau_xz',-1),('tau_yz',-1), deltap]
                    data = (us, vs, nan, ub, vb, nan, nan, nan)
                elif experiment == 'd':
                    #  This is supposed to be:  [('uvel',0),('wvel',0),('uvel',-1),('tau_xz',-1), deltap]
                    data = (us, nan, nan, nan, nan)
                elif experiment == 'f':  # should be: x, y, zs, vx, vy, vz
                    #  This is supposed to be: [('usurf',None),('uvel',0),('vvel',0),('wvel',0)]
                    data = (nan, uvelSprime, vvelS, nan)

                # Write a "standard" ISMIP-HOM file (example file name: "cis1a020.txt") in the "output" subdirectory 
                ISMIP_HOMfilename = out_file.replace('.out.nc','.txt')
                ISMIP_HOMfile = open(ISMIP_HOMfilename,'w')
                for i, x in enumerate(xx):
                    for j, y in enumerate(yy):
                        if experiment in ('a','c','f'):  # include x and y positions
                            ISMIP_HOMfile.write('\t'.join(map(str,[x,y]+[v[j,i] for (v) in data]))+'\n')
                        else:  # only include x position
                            ISMIP_HOMfile.write('\t'.join(map(str,[x]+[v[j,i] for (v) in data]))+'\n')
                ISMIP_HOMfile.close()
                netCDFfile.close()
            except:
                print 'Warning: The CISM output file for experiment '+experiment+' at size '+str(size)+' could NOT be read successfully!'
                raise

    # =========================================================
    #  Now actually analyze the results.
    # =========================================================
    
    # Loop over the experiments requested on the command line
    for experiment in experiments:
        print 'ISMIP-HOM', experiment.upper()

        # Create the figure on which the plot axes will be placed
        figure = pyplot.figure(subplotpars=matplotlib.figure.SubplotParams(top=.85,bottom=.15))
        figure.text(0.5,0.92,'ISMIP-HOM Experiment '+experiment.upper(),horizontalalignment='center',size='large')
        if args.subtitle:
            figure.text(0.5,0.89,args.subtitle,horizontalalignment='center',size='small')
        figure.text(0.5,0.1,'Normalized X coordinate',horizontalalignment='center',size='small')
        figure.text(0.06,0.5,'Ice Speed (m/a)',rotation='vertical',verticalalignment='center')
        # Create the (three column) legend
        prop = matplotlib.font_manager.FontProperties(size='x-small')
        Line2D = matplotlib.lines.Line2D([],[],color=(0,0,0))
        figure.legend([Line2D],['Model Output'],loc=(0.1,0.05),prop=prop).draw_frame(False)
        Line2D.set_linestyle(':')
        Line2D.set_color((1,0,0))
        Patch = matplotlib.patches.Patch(edgecolor=None,facecolor=(1,0,0),alpha=0.25)
        figure.legend([Line2D,Patch],['Full Stokes Mean','Full Stokes Std. Dev.'],loc=(0.3,0.02),prop=prop).draw_frame(False)
        Line2D.set_color((0,0,1))
        Patch.set_facecolor((0,0,1))
        figure.legend([Line2D,Patch],[nonFSmodelType+' Mean',nonFSmodelType+' Std. Dev.'],loc=(0.55,0.02),prop=prop).draw_frame(False)

        # Loop over the sizes requested on the command line
        for i, size in enumerate(map(int,sizes)):
            try:
                if experiment == 'f': 
                    if size != 100 or len(sizes) > 1:
                        print 'NOTE: Experiment f uses a domain size of 100 km only'
                    size = 100
                
                out_file = pattern.replace('-?','-'+experiment).replace('????',str(size).zfill(4))

                # Create the plot axes for this domain size
                if len(sizes) == 1:
                    axes = figure.add_subplot(111)
                else:
                    axes = figure.add_subplot(2,3,i+1)
                    for tick in axes.xaxis.get_major_ticks():
                        tick.label1.set_fontsize('xx-small')
                    for tick in axes.yaxis.get_major_ticks():
                        tick.label1.set_fontsize('xx-small')
                axes.set_title('%d km' % size, size='medium')

                # Get the Glimmer output data
                glimmerData = read(out_file.replace('.out.nc','.txt'),experiment)
                # The Glimmer data is on a staggered grid;
                # Interpolate to obtain the value at x=0 and x=1
                # using periodic boundary conditions
                #v = (glimmerData[0][1] + glimmerData[-1][1])/2
                #glimmerData = [(0.0,v)]+glimmerData+[(1.0,v)]

                # Plot the Glimmer data
                axes.plot([row[0] for row in glimmerData],
                          [row[1] for row in glimmerData],color='black')

                # Get the data from other models for comparison
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
                        #print 'name, exp, size', modelName, modelExperiment, modelSize
                        if modelName == 'aas1':
                            # Skip the 'aas1' model because its output files in the tc-2007-0019-sp2.zip file do not follow the proper naming convention.  MJH 11/5/13
                            continue
                        if (modelExperiment != experiment) or (modelExperiment != 'f' and int(modelSize) != size) \
                        or (not args.all_ps and not modelName in lmlaModels + fullStokesModels):
                            continue # continue next loop iteration if not the size or not the experiment desired or if we just want FO comparison and this model is not FO or FS.
                        elif (modelExperiment == 'f'):
                            if (modelSize == '001'):
                                  continue # ignore the sliding version for now
                            if modelName == 'cma1':
                                  continue  # the cma1 'f' experiments made the x,y coords dimensionless instead of dimensional - ignore for convenience
                        print 'Using data from file:',os.path.join(path,filename)
                        data = read(os.path.join(path,filename),experiment)
                        if modelName in fullStokesModels:
                            index = fullStokes
                        else:
                            index = firstOrder
                        count[index] += 1

                        #axes.plot([row[0] for row in data], [row[1] for row in data] )   ## OPTIONAL: print out every individual model in its native x-coordinates.

                        # Interpolate onto the x values from the Glimmer model run
                        for (i,target) in enumerate([row[0] for row in glimmerData]):
                            below = -99999.0
                            above =  99999.0
                            for (j,x) in enumerate([row[0] for row in data]):
                                if  below <  x <= target: b,below = j,x
                                if target <= x <  above:  a,above = j,x
                            if above == below:
                                v = data[a][1]
                            else:
                                if below == -99999.0: # Use the periodic boundary condition at x = 0
                                    xBelow = data[-1][0] - 1
                                    vBelow = data[-1][1]
                                else:
                                    xBelow,vBelow = data[b]
                                
                                if above ==  99999.0: # Use the periodic boundary condition at x = 1
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

                # Calculate statistics of the other model results
                if sum(count) == 0:
                    print 'To compare with other models you need to download the ISMIP-HOM results from: http://www.the-cryosphere.net/2/95/2008/tc-2-95-2008-supplement.zip and unzip the contained file tc-2007-0019-sp2.zip into a directory named ismip_all.  The ismip_all directory must be in the directory from which you are running this script.'
                else:
                    # Find the mean and standard deviation of the velocities at each x
                    for index in (firstOrder,fullStokes):
                        if count[index] == 0:
                            continue
                        mean = list()
                        standardDeviation = list()
                        for i in range(len(glimmerData)):
                            mean.append(sumV[index][i]/count[index])
                            standardDeviation.append(sqrt(sumV2[index][i]/count[index]-mean[-1]**2))

                        # Plot the mean using a dotted line
                        color = (index,0,1-index) # blue for first order (index=0); red for full Stokes (index=1)
                        x = [row[0] for row in glimmerData]
                        axes.plot(x,mean,':',color=color)

                        # Plot a filled polygon showing the mean plus and minus one standard deviation
                        meanMinusSD = [m-sd for (m,sd) in zip(mean,standardDeviation)]
                        meanPlusSD  = [m+sd for (m,sd) in zip(mean,standardDeviation)]
                        x = x + list(reversed(x))
                        y = meanPlusSD + list(reversed(meanMinusSD))
                        axes.fill(x,y,facecolor=color,edgecolor=color,alpha=0.25)

                        if index == firstOrder:
                            # Calculate some statistics comparing the Glimmer data with the other models
                            pcterror = [100.0*abs(glimmer-others)/others for (glimmer,others) in zip([row[1] for row in glimmerData],mean)]
                            abserror = [abs(glimmer-others)        for (glimmer,others) in zip([row[1] for row in glimmerData],mean)]
                            maximum = max(pcterror)
                            position = glimmerData[pcterror.index(maximum)][0]
                            total   = sum([e for e in pcterror])
                            compare = sum([(s/m) for (s,m) in zip(standardDeviation,mean)])
                            n = len(glimmerData)
                            #print '\t'.join([str(size)+' km',str(total/n),str(compare/n),str(position)])
                            print 'Size='+str(size)+' km' 
                            print '  Mean percent error along flowline of CISM relative to mean of first-order models='+str(total/float(n))+'%'
                            print '  Mean COD (stdev/mean) along flowline of mean of first-order models (excluding CISM)='+str(compare/float(n)*100.0)+'%'
                            print '  Max. CISM percent error='+str(maximum)+'% at x-position '+str(position)
                            print '  Max. CISM absolute error='+str(max(abserror))+' m/yr at x-position '+str(glimmerData[abserror.index(max(abserror))][0])

            except:
                print "Error in analyzing/plotting experiment ",experiment," at size ",size," km"

        if savePlotInFile:
            plot_file = pattern.replace('-?','-'+experiment).replace('.????','').replace('.out.nc',plotType)
            print 'Writing:', plot_file
            pyplot.savefig(plot_file)

        # Experiment f can also have a surface profile plotted
        if experiment == 'f':
            # rather than getting the data from the text file, we are going to read it directly.
            # this is because the velocities and usrf are on different grids, so it is difficult to include them
            # both in the standard ISMIP-HOM text file format that has a single x,y coord. system
            size = 100
            out_file = pattern.replace('-?','-'+experiment).replace('????',str(size).zfill(4))
            netCDFfile = NetCDFFile(out_file,'r')
            if netCDF_module == 'Scientific.IO.NetCDF':
                thkscale = netCDFfile.variables['thk'].scale_factor
            else:
                thkscale = 1.0
            usurf = netCDFfile.variables['usurf'][-1,:,:] * thkscale  # get last time level
            x1 = netCDFfile.variables['x1'][:]
            y1 = netCDFfile.variables['y1'][:]

            #  Create the usurf figure
            ufigure = pyplot.figure(subplotpars=matplotlib.figure.SubplotParams(top=.85,bottom=.15))
            ufigure.text(0.5,0.92,'ISMIP-HOM Experiment F: Surface elevation',horizontalalignment='center',size='large')
            if args.subtitle:
                ufigure.text(0.5,0.89,args.subtitle,horizontalalignment='center',size='small')
            ufigure.text(0.5,0.1,'X coordinate',horizontalalignment='center',size='small')
            ufigure.text(0.06,0.5,'upper surface (m)',rotation='vertical',verticalalignment='center')
            # Create the (three column) legend
            prop = matplotlib.font_manager.FontProperties(size='x-small')
            Line2D = matplotlib.lines.Line2D([],[],color=(0,0,0))
            ufigure.legend([Line2D],['Model Output'],loc=(0.1,0.05),prop=prop).draw_frame(False)
            Line2D.set_linestyle(':')
            Line2D.set_color((1,0,0))
            Patch = matplotlib.patches.Patch(edgecolor=None,facecolor=(1,0,0),alpha=0.25)
            ufigure.legend([Line2D,Patch],['Full Stokes Mean','Full Stokes Std. Dev.'],loc=(0.3,0.02),prop=prop).draw_frame(False)
            Line2D.set_color((0,0,1))
            Patch.set_facecolor((0,0,1))
            ufigure.legend([Line2D,Patch],[nonFSmodelType+' Mean',nonFSmodelType+' Std. Dev.'],loc=(0.55,0.02),prop=prop).draw_frame(False)
            # Create the plot axes 
            axes2 = ufigure.add_subplot(111)
            axes2.set_title('%d km' % size, size='medium')

            # Convert CISM output data to the rotated coord system used by the problem setup
            alpha = -3.0 * pi/180  # defined in run script
            # use integer floor division operator to get an index close to the center  TODO should be interpolating if needed...
            yp = len(y1)//2
            # calculate rotated xprime, zprime coordinates along the surface (this is the coord. sys. used for this test case)
            xprime =  x1 * cos(alpha) + (usurf[yp,:]-7000.0) * sin(alpha)
            xprime = xprime/1000.0 - 50.0
            zprime = -x1 * sin(alpha) + (usurf[yp,:]-7000.0) * cos(alpha)
            # Plot CISM output
            axes2.plot(xprime, zprime, color='black')

            # create glimmerData so we can re-use the code from above
            glimmerData = list()
            for i in range(len(xprime)):
                glimmerData.append(tuple([xprime[i], zprime[i]]))

            # Now plot the other models - yucky code copied from above
            # Get the data from other models for comparison
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
                    #print 'name, exp, size', modelName, modelExperiment, modelSize
                    if modelName == 'aas1':
                        # Skip the 'aas1' model because its output files in the tc-2007-0019-sp2.zip file do not follow the proper naming convention.  MJH 11/5/13
                        continue
                    if (modelExperiment != experiment) or (modelExperiment != 'f' and int(modelSize) != size) \
                    or (not args.all_ps and not modelName in lmlaModels + fullStokesModels):
                        continue # continue next loop iteration if not the size or not the experiment desired or if we just want FO comparison and this model is not FO or FS.
                    elif (modelExperiment == 'f'):
                        if (modelSize == '001'):
                            continue # ignore the sliding version for now
                        if modelName == 'cma1':
                            continue  # the cma1 'f' experiments made the x,y coords dimensionless instead of dimensional - ignore for convenience
                    print 'Using data from file:',os.path.join(path,filename)
                    data = read(os.path.join(path,filename), experiment='f-elevation')
                    if modelName in fullStokesModels:
                        index = fullStokes
                    else:
                        index = firstOrder
                    count[index] += 1

                    #axes2.plot([row[0] for row in data], [row[1] for row in data] )   ## OPTIONAL: print out every individual model in its native x-coordinates.

                    # Interpolate onto the x values from the Glimmer model run
                    for (i,target) in enumerate([row[0] for row in glimmerData]):
                        below = -99999.0
                        above =  99999.0
                        for (j,x) in enumerate([row[0] for row in data]):
                            if  below <  x <= target: b,below = j,x
                            if target <= x <  above:  a,above = j,x
                        if above == below:
                            v = data[a][1]
                        else:
                            if below == -99999.0: # Use the periodic boundary condition at x = 0
                                xBelow = data[-1][0] - 1
                                vBelow = data[-1][1]
                            else:
                                xBelow,vBelow = data[b]
                            if above ==  99999.0: # Use the periodic boundary condition at x = 1
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

                    # Calculate statistics of the other model results
                    if sum(count) == 0:
                        print 'To compare with other models you need to download the ISMIP-HOM results from: http://www.the-cryosphere.net/2/95/2008/tc-2-95-2008-supplement.zip and unzip the contained file tc-2007-0019-sp2.zip into a directory named ismip_all.  The ismip_all directory must be in the directory from which you are running this script.'
                    else:
                        # Find the mean and standard deviation of the velocities at each x
                        for index in (firstOrder,fullStokes):
                            if count[index] == 0:
                                continue
                            mean = list()
                            standardDeviation = list()
                            for i in range(len(glimmerData)):
                                mean.append(sumV[index][i]/count[index])
                                standardDeviation.append(sqrt(sumV2[index][i]/count[index]-mean[-1]**2))

                            # Plot the mean using a dotted line
                            color = (index,0,1-index) # blue for first order (index=0); red for full Stokes (index=1)
                            x = [row[0] for row in glimmerData]
                            axes2.plot(x,mean,':',color=color)

                            # Plot a filled polygon showing the mean plus and minus one standard deviation
                            meanMinusSD = [m-sd for (m,sd) in zip(mean,standardDeviation)]
                            meanPlusSD  = [m+sd for (m,sd) in zip(mean,standardDeviation)]
                            x = x + list(reversed(x))
                            y = meanPlusSD + list(reversed(meanMinusSD))
                            axes2.fill(x,y,facecolor=color,edgecolor=color,alpha=0.25)

            if savePlotInFile:
                plot_file = pattern.replace('-?','-'+experiment+'-SurfaceElevation').replace('.????','').replace('.out.nc',plotType)
                print 'Writing:', plot_file
                pyplot.savefig(plot_file)

        # Experiment f should be run for one size (100 km) only
        if experiment == 'f': 
              break

    if not savePlotInFile:
        pyplot.show()


# Run only if this is being run as a script.
if __name__=='__main__':
    
    # get the command line arguments
    args = parser.parse_args()
    
    # run the script
    sys.exit(main())
