#!/usr/bin/env python
#
# Magnus Hagdorn
#
# plot difference between numeric and exact solution

import os.path
import pylab, matplotlib
import Scientific.IO.NetCDF
import numpy.ma
from optparse import OptionParser


if __name__ == '__main__':

    usage = """usage: %prog [options] file.nc

plot difference between exact and simulated solution at specified time"""

    parser = OptionParser(usage=usage)
    parser.add_option("-o","--output",metavar="FILE",help="write image to file. image type is determined by file suffix")
    parser.add_option("-T","--time-slice",metavar="T",default=-1,type='int',help="extract data for time slice T")
    (options, args) = parser.parse_args()

    if len(args)!=1:
        parser.error('Expecting one input file')

    infile = Scientific.IO.NetCDF.NetCDFFile(args[0],'r')
    diff = infile.variables['thk'][options.time_slice,:,:] - infile.variables['thke'][options.time_slice,:,:]
    mask = numpy.where(infile.variables['thk'][options.time_slice,:,:] + infile.variables['thke'][options.time_slice,:,:] > 0, False, True)
    diff = numpy.ma.array(diff,mask=mask)
    extent = [infile.variables['x1'][0]/1000.,infile.variables['x1'][-1]/1000.,
                           infile.variables['y1'][0]/1000.,infile.variables['y1'][-1]/1000.]
    #()

    pylab.title("Experiment %s"%infile.title)
    pylab.imshow(diff,origin='lower',
                 norm=matplotlib.colors.Normalize(vmin=-200,vmax=200),
                 cmap=matplotlib.cm.RdBu_r,
                 extent=extent)
    pylab.colorbar()
    pylab.contour(diff,extent=extent,colors='k')

    pylab.contour(infile.variables['thk'][options.time_slice,:,:],[0],colors='grey',extent=extent)

    if options.output!=None:
        pylab.savefig(options.output)
    else:
        pylab.show()
    
