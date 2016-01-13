#!/usr/bin/env python

"""
This script plots the results of the Ross Ice Shelf Experiments.
The plots compare the model output with measured velocities from the 
Ross Ice Shelf Geophysical and Glaciological Survey 1973-78 (RIGGS).
"""

# Written April 5, 2010 by Glen Granzow at the University of Montana.
# Reconfigured by Joseph H Kennedy at ORNL on August 7, 2015 to work with the regression testing
#     NOTE: Did not adjust inner workings except where needed.

import os
import sys
import glob
import numpy

from netCDF import *
from matplotlib import pyplot, colors

import argparse
parser = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

# The command line options
# ------------------------
parser.add_argument('-o', '--output-dir', default='./output',
        help="Directory containing the tests output files. Warning: if there is a" \
        +"path passed via the -f/--out-file option, this argument will be" \
        +"ignored.")

parser.add_argument('-f', '--output-file',
        help="The tests output file you would like to plot. If a path is" \
        +"passed via this option, the -o/--output-dir option will be ignored.")

parser.add_argument('-c', '--colormap', default='jet',
        help='Specify a color map')
parser.add_argument('-m', '--use-mask', action='store_true',
        help='Mask the grounded region')
parser.add_argument('-n', '--ncontours', type=int, default=0,
        help='Plot n filled contours instead of a continuous color plot')
parser.add_argument('-v','--vmax', type=float, default=0.0,
        help='Set upper limit of color scale')


# ===========================================================
# Define some variables and functions used in the main script
# ===========================================================

def get_in_file():
    if args.output_file:
        out_d, out_f = os.path.split(args.output_file)
        if out_d:
            args.output_dir = out_d
            args.output_file = out_f
    
        print("\nUsing "+os.path.join(args.output_dir, args.output_file)+"\n")

    else:
        outpath = os.path.join(args.output_dir, '*.out.nc')
        matching = glob.glob(outpath)
        if len(matching) == 1:
            newest = matching[0]
            print("\nUsing "+newest+"\n")
    
        elif len(matching) > 1:
            newest = max(matching, key=os.path.getmtime)
            print("\nWARNING: MULTIPLE *.out.nc FILES DETECTED!")
            print(  "==========================================")
            print(  "Ploting the most recently modified file in the output directory:")
            print(  "    "+newest)
            print(  "To plot another file, specify it with the -f/--outfile option.\n")
            
        else:
            print("\nERROR: NO *.out.nc FILES DETECTED!")
            print(  "==================================")
            print(  "Either specify a location to look for the test output")
            print(  "files with the -o/--output-dir option, or the test output")
            print(  "file with the -f/--output-file option.\n")
            sys.exit(1)

        args.output_file = os.path.basename(newest)
    
    filein = NetCDFFile(os.path.join(args.output_dir, args.output_file),'r')
     
    return filein


# cmap_discretize returns a discretized version of a colormap
# From http://www.scipy.org/Cookbook/Matplotlib/ColormapTransformations 
def cmap_discretize(cmap, N):
     from scipy import interpolate
     cdict = cmap._segmentdata.copy()
     colors_i = numpy.linspace(0,1,N)
     indices  = numpy.linspace(0,1,N+1)
     for key in ('red','green','blue'):
         D = numpy.array(cdict[key])
         I = interpolate.interp1d(D[:,0], D[:,1])
         clrs = I(colors_i)
         A = numpy.zeros((N+1,3), float)
         A[:,0] = indices
         A[1:,1] = clrs
         A[:-1,2] = clrs
         L = []
         for l in A:
             L.append(tuple(l))
         cdict[key] = tuple(L)
     return colors.LinearSegmentedColormap('colormap',cdict,1024)


# =========================
# Actual script starts here
# =========================
def main():
    """
    Plot the Ross test results.
    """

    filein = get_in_file()

    velnorm = numpy.array(filein.variables['velnorm'][0,0,:,:])
    if netCDF_module == 'Scientific.IO.NetCDF':
        velnorm = velnorm * filein.variables['velnorm'].scale_factor
    filein.close()

    if args.use_mask:
        mask_file = NetCDFFile(os.path.join(args.output_dir, args.output_file.replace('out.nc','.nc')), 'r')
        mask = numpy.array(mask_file.variables['kinbcmask'][0,:,:])
        mask_file.close()

    # Read the grid coordinates 
    gird_file = open(os.path.join('data','111by147Grid.dat'))
    nx,ny = 147,111
    x,y = list(),list()
    for i in range(4): 
        gird_file.readline()
    for i in range(ny):
        y.append(float(gird_file.readline()))
    for i in range(3):
        gird_file.readline()
    for i in range(nx):
        x.append(float(gird_file.readline()))
    gird_file.close()

    # Read the RIGGS data
    riggs_file = open(os.path.join('data','riggs_clean.dat'))
    latitude,longitude,riggs_velocity,glimmer_velocity = list(),list(),list(),list()
    lat0,lon0,vel0 = list(),list(),list()
    for line in riggs_file:
        tokens = [float(data) for data in line.split()]
        # Convert from degrees, minutes and seconds
        lat = -(tokens[3] + tokens[4]/60 + tokens[5]/60**2)
        lon = -(tokens[6] + tokens[7]/60 + tokens[8]/60**2)*tokens[9]
        if y[1] < lat < y[-2] and x[1] < lon < x[-2]:
            # Interpolate the Glimmer data onto this point
            for i in range(nx):
                if lon < x[i+1]: break
            for j in range(ny):
                if lat < y[j+1]: break
            alpha = (lon-x[i])/(x[i+1]-x[i])
            beta  = (lat-y[j])/(y[j+1]-y[j])
            v = (1-alpha)*((1-beta)*velnorm[j,i]  +beta*velnorm[j+1,i]) \
               +   alpha *((1-beta)*velnorm[j,i+1]+beta*velnorm[j+1,i+1])
            #print(alpha)
            #print(beta)
            #print(v)
            if v != 0:
                latitude.append(lat)
                longitude.append(lon)
                riggs_velocity.append(tokens[10])
                glimmer_velocity.append(v)
            else:
                print 'Glimmer velocity is zero at (%g, %g); RIGGS velocity is %g' % (lon,lat,tokens[10])
                lat0.append(lat)
                lon0.append(lon)
                vel0.append(tokens[10])
        else:
            print 'RIGGS data is off the grid at (%g, %g); RIGGS velocity is %g' % (lon,lat,tokens[10])
    print 'These points are not included in Chi-squared.'
    riggs_file.close()

    # Present the output
    # ------------------
    sigma = 30
    X2 = sum([((v1-v2)/sigma)**2 for v1,v2 in zip(riggs_velocity,glimmer_velocity)])
    print
    print 'Chi-squared for',len(riggs_velocity),'points is', X2
    print 'The maximum velocity from Glimmer/CISM is', numpy.max(velnorm)
    print 'The maximum velocity from the RIGGS data is', max(riggs_velocity)

    # Create a scatter plot
    pyplot.figure(1)
    pyplot.clf()
    pyplot.plot([0,1800],[0,1800],color='red',linewidth=2,alpha=0.5)
    pyplot.scatter(glimmer_velocity,riggs_velocity,color='blue')
    pyplot.axis('equal')
    pyplot.axis([0,1800,0,1800])
    pyplot.xticks(numpy.linspace(0.0, 1800.0, num=10, endpoint=True))
    pyplot.yticks(numpy.linspace(0.0, 1800.0, num=10, endpoint=True))
    pyplot.ylabel('Measured velocity from RIGGS (meters/year)')
    pyplot.xlabel('Model velocity from Glimmer/CISM (meters/year)')
    pyplot.title('Ross Ice Shelf Experiment')

    # Create a color plot of Glimmer and RIGGS velocities
    pyplot.figure(2)
    pyplot.clf()
    plot_discarded_points = False
    # The dimensions of X and Y for pcolor should be one greater than those of the data
    midx = sum(map(numpy.array,(x[:-1],x[1:])))/2
    midy = sum(map(numpy.array,(y[:-1],y[1:])))/2
    x = x[:1] + list(midx) + x[-1:]
    y = y[:1] + list(midy) + y[-1:]
    x,y = numpy.meshgrid(x,y)
    # Set plot parameters
    vmax = numpy.max(velnorm)
    if args.vmax > 0:
        vmax = args.vmax
    norm  = colors.Normalize(0,vmax)
    cmap  = eval('pyplot.cm.' + args.colormap)
    ticks = None
    if args.ncontours > 1:
        cmap = cmap_discretize(cmap, args.ncontours)
        ticks = numpy.arange(args.ncontours+1) * vmax/args.ncontours
    if args.use_mask:
        velnorm = numpy.ma.masked_array(velnorm,mask==1)
    # Make the actual plot
    pyplot.pcolor(x,y,velnorm,norm=norm,cmap=cmap)
    if plot_discarded_points:
        pyplot.scatter(lon0,lat0,s=50,c='white')
    pyplot.scatter(longitude,latitude,s=50,c=riggs_velocity,norm=norm,cmap=cmap)
    pyplot.axis([x[0,0],x[-1,-1],y[0,0],y[-1,-1]])
    pyplot.axes().set_aspect('equal')
    pyplot.colorbar(orientation='vertical',fraction=0.05,pad=0.01,shrink=0.93,ticks=ticks)
    pyplot.title('Ross Ice Shelf Experiment - Velocity (meters/year)')
    print
    print 'FIGURE 2:'
    print 'If Glimmer/CISM perfectly predicted the velocities measured by the Ross'
    print 'Ice Shelf Geophysical and Glaciological Survey (RIGGS), the colors in'
    print 'the circles (which represent the RIGGS velocities) would be the same as'
    print 'the color of the background (which represents the velocities calculated'
    print 'by Glimmer/CISM).'
    print


    # Create a contour plot of Glimmer to compare to paper Figure 3.
    pyplot.figure(3, facecolor='w', figsize=(12, 4), dpi=72)
    pyplot.clf()
    pyplot.subplot(1,2,1)
    pyplot.contour(velnorm, numpy.linspace(0.0, 5000.0, num=51, endpoint=True) ) 
    pyplot.axis('equal')

    pyplot.subplot(1,4,3)
    pyplot.contour(velnorm, numpy.concatenate((numpy.linspace(0.0, 175.0, num=8, endpoint=True), numpy.linspace(200.0, 2000.0, num=10, endpoint=True))) ) 
    pyplot.axis('equal')
    pyplot.xlim((110,147))
    pyplot.ylim((1,58))

    print
    print 'FIGURE 3:'
    print 'Compare this plot to Figure 3 in the EISMINT paper: '
    print 'MacAyeal, et al.: An ice-shelf model test based on the Ross Ice Shelf, Antarctica, Ann. Glaciol., 23, 46-51, 1996'
    print 'Try URL: http://www.igsoc.org/annals.old/23/igs_annals_vol23_year1996_pg46-51.pdf'
    print

    pyplot.show()

# Run only if this is being run as a script.
if __name__=='__main__':
    
    # get the command line arguments
    args = parser.parse_args()
    
    # run the script
    sys.exit(main())

