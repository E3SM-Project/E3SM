#!/usr/bin/env python
#
# Magnus Hagdorn
#
# plot time steps and run time


import os.path
import pylab
from optparse import OptionParser

def parse_title(title):
    """Parse title string."""

    config = title.split('_')
    exp = config[1]
    solver = config[2]
    dx = int(config[3][:-2])
    dt = float(config[4][:-1])

    return (exp,solver,dx,dt)

if __name__ == '__main__':

    usage = """usage: %prog [options] [RESULTS]

parse RESULTS file and plot model runtimes"""

    parser = OptionParser(usage=usage)
    parser.add_option("-o","--output",metavar="FILE",help="write image to file. image type is determined by file suffix")
    (options, args) = parser.parse_args()

    if len(args)>1:
        parser.error('Expecting at most one results file')
    elif len(args)==1:        
        inname = args[0]
    else:
        inname = 'results'


    # parse results file
    runs = {}
    for line in file(inname).readlines():
        # ignore comments and empty lines
        line = line.strip()
        pos = line.find('#')
        if pos>-1:
            line = line[:pos]
        if len(line)==0:
            continue
        line = line.split()
        fname = line[0][:-len('.config')]
        (exp_name,solver,dx,dt) = parse_title(fname)
        time = float(line[1])
        if solver not in runs:
            runs[solver] = {}
        if exp_name not in runs[solver]:
            runs[solver][exp_name]={}
        runs[solver][exp_name][dx] = [time,dt]


    # start plotting
    pylab.figure(1)
    area_dt = pylab.subplot(211)
    area_dt.set_xlabel("dx [km]")
    area_dt.set_ylabel("dt [a]")
    
    area_rt = pylab.subplot(212)
    area_rt.set_xlabel("dx [km]")
    area_rt.set_ylabel("run time [s]")

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

            rt = []
            dt = []
            for x in dx:
                rt.append(runs[s][e][x][0])
                dt.append(runs[s][e][x][1])

            # plot line
            area_rt.semilogy(dx,rt,ls=styles[s],color=colours[done_grid[e]],label=title)
            area_dt.plot(dx,dt,ls=styles[s],color=colours[done_grid[e]])

    area_rt.legend()

    if options.output!=None:
        pylab.savefig(options.output)
    else:
        pylab.show()
    
