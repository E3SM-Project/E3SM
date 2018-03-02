#!/usr/bin/env python
#
# Magnus Hagdorn
#
# plot relative volume errors


import Scientific.IO.NetCDF
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

def calc_verror(cffile):

    verror = (cffile.variables['ivol'][:]-cffile.variables['ivole'][:])/cffile.variables['ivole'][:]
    verror[0] = 0.
    return verror

if __name__ == '__main__':

    usage = """usage: %prog [options] file1.nc [filen.nc]

plot relative volume error as a function of time"""

    parser = OptionParser(usage=usage)
    parser.add_option("-o","--output",metavar="FILE",help="write image to file. image type is determined by file suffix")
    (options, args) = parser.parse_args()

    if len(args)<1:
        parser.error('Expecting at least one file')


    # start plotting
    pylab.figure(1)
    vol_err = pylab.subplot(111)
    vol_err.set_xlabel = 'time [ka]'
    vol_err.set_ylabel = 'relative volume error'

    styles = {'non-lin':'-','lin':':','ADI':'--'}
    colours = ['red','green','blue','cyan','yellow','orange','magenta','pink']

    plotted = {}
    i = 0
    for f in args:
        cffile = Scientific.IO.NetCDF.NetCDFFile(f,'r')
        (exp_name,solver,dx,dt) = parse_title(cffile.title)
        title = 'dx=%skm, dt=%sa'%(dx,dt)
        if title not in plotted:
             plotted[title] = i
             i = i+1
             label = title
        else:
            label = None
        vol_err.plot(cffile.variables['time'][:],calc_verror(cffile),ls=styles[solver],color=colours[plotted[title]],label=label)
        cffile.close()
    
    vol_err.legend()

    if options.output!=None:
        pylab.savefig(options.output)
    else:
        pylab.show()
