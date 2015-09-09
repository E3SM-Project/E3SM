#!/usr/bin/env python
#
# Magnus Hagdorn
#
# plot absolute dome and max errors


import numpy, Scientific.IO.NetCDF
import pylab
from optparse import OptionParser

def parse_title(title):
    """Parse title string."""

    t = title.split(',')

    exp_name = t[0][-1]
    solver = t[1].strip()
    dx = t[2].strip()[:-2]
    dt = t[3].strip()[:-1]

    return (exp_name,solver,dx,dt)

if __name__ == '__main__':

    usage = """usage: %prog [options] file1.nc [filen.nc]

plot model errors as function of grid spacing"""

    parser = OptionParser(usage=usage)
    parser.add_option("-o","--output",metavar="FILE",help="write image to file. image type is determined by file suffix")
    (options, args) = parser.parse_args()

    if len(args)<1:
        parser.error('Expecting at least one file')

    runs = {}

    # load data
    for f in args:
        cffile = Scientific.IO.NetCDF.NetCDFFile(f,'r')
        (exp_name,solver,dx,dt) = parse_title(cffile.title)
        if solver not in runs:
            runs[solver] = {}
        if exp_name not in runs[solver]:
            runs[solver][exp_name]={}
        diff = cffile.variables['thke'][-1,:,:] - cffile.variables['thk'][-1,:,:]
        centre = (numpy.shape(diff)[0]-1)/2
        dome_e = diff[centre,centre]
        diff = numpy.ravel(diff)
        max_e = max(abs(diff))
        runs[solver][exp_name][int(dx)] = [dome_e,max_e]
        cffile.close()

    # start plotting
    pylab.figure(1)

    area_dome = pylab.subplot(211)
    area_dome.set_xlabel("dx [km]")
    area_dome.set_ylabel("absolute dome error [m]")

    area_max = pylab.subplot(212)
    area_max.set_xlabel("dx [km]")
    area_max.set_ylabel("absolute maximum error [m]")

    styles = {'non-lin':'-','lin':':','ADI':'--'}
    colours = ['red','green','blue','cyan','yellow','orange','magenta','pink']

    done_grid = {}
    i = 0
    for s in runs:
        for e in runs[s]:
            if e not in done_grid:
                done_grid[e] = i
                i = i + 1
                title = e
            else:
                title = None

            dx = runs[s][e].keys()
            dx.sort()

            # get errors
            error_dome = []
            error_max = []
            for x in dx:
                error_dome.append(runs[s][e][x][0])
                error_max.append(runs[s][e][x][1])

            # plot line
            area_dome.plot(dx,error_dome,ls=styles[s],color=colours[done_grid[e]],label=title)
            area_max.plot(dx,error_max,color=colours[done_grid[e]])


    area_dome.legend()

    if options.output!=None:
        pylab.savefig(options.output)
    else:
        pylab.show()
