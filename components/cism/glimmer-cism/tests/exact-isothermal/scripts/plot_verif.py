#!/usr/bin/env python
#
# Magnus Hagdorn
#
# fill crazy dictionary and plot it

import Numeric, Scientific.IO.NetCDF, pygsl, sys
import PyGMT
from pygsl import histogram

def parse_title(title):
    """Parse title string."""

    t = title.split(',')

    exp_name = t[0][-1]
    solver = t[1].strip()
    dx = t[2].strip()[:-2]
    dt = t[3].strip()[:-1]

    return (exp_name,solver,dx,dt)

def process(fname):
    """Process verif output file.

    fname: name of file

    return: tuple of (exp_name, solver, dx, dt, dome_e, max_e, min_e, mean_e, sd_e)"""

    # extract errors
    ncf = Scientific.IO.NetCDF.NetCDFFile(fname)
    diff = ncf.variables['thke'][-1,:,:] - ncf.variables['thk'][-1,:,:]
    centre = (Numeric.shape(diff)[0]-1)/2
    dome_e = diff[centre,centre]
    diff = Numeric.ravel(diff)
    max_e = max(abs(diff))
    min_e = min(abs(diff))
    hist = histogram.histogram(100)
    hist.set_ranges_uniform(PyGMT.round_down(min_e),PyGMT.round_up(max_e))
    for e in diff.tolist():
        hist.increment(e)
    mean_e = hist.mean()
    sd_e = hist.sigma()

    (exp_name,solver,dx,dt) = parse_title(ncf.title)
    
    ncf.close()

    return (exp_name, solver, dx, dt, dome_e, max_e, min_e, mean_e, sd_e)

def munch_verif(fname, verif):
    """Read in verif result file and fill dictionary."""

    (exp_name, solver, dx, dt, dome_e, max_e, min_e, mean_e, sd_e) = process(fname)

    if exp_name not in verif:
        verif[exp_name] = {}
    data_exp = verif[exp_name]
    
    if solver not in data_exp:
        data_exp[solver] = {}
    data_solver = data_exp[solver]

    if dx not in data_solver:
        data_solver[dx] = {}
    data_dx = data_solver[dx]

    data_dx[float(dt)] = (dome_e, max_e, min_e, mean_e, sd_e)

    return verif

if __name__ == '__main__':

    # read data
    verif = {}
    for f in sys.argv[1:]:
        munch_verif(f,verif)

    exp = verif.keys()[0]

    grids = {}
    gn = 0
    for s in verif[exp]:
        for g in verif[exp][s]:
            if g not in grids:
                grids[g] = gn
                gn = gn+1
    print grids
    key_y=2.5
    ysize = 7.
    dy = 1.

    # create plot
    plot = PyGMT.Canvas('blub.ps',size='A4')
    plot.defaults['LABEL_FONT_SIZE']='12p'
    plot.defaults['ANOT_FONT_SIZE']='10p'

    bigarea = PyGMT.AreaXY(plot,size=[30,30])

    bigarea.text([7.5,2*ysize+3*dy+key_y],"Experiment %s"%exp,textargs='20 0 0 CM')

    area_dome = PyGMT.AutoXY(bigarea,pos=[0.,key_y+dy],size=[15,ysize],logx=True)
    area_dome.xlabel = 'time step [a]'
    area_dome.ylabel = 'dome error [m]'
    area_dome.axis = 'WeSn'

    area_max = PyGMT.AutoXY(bigarea,pos=[0.,ysize+2*dy+key_y],size=[15,ysize],logx=True)
    area_max.xlabel = 'time step [a]'
    area_max.ylabel = 'max error [m]'
    area_max.axis = 'Wesn'

    s_keyarea = PyGMT.KeyArea(bigarea,pos=[0.,-0.8],size=[5,key_y])
    s_keyarea.num=[1,7]

    e_keyarea = PyGMT.KeyArea(bigarea,pos=[5.,-0.8],size=[10.,key_y])
    e_keyarea.num=[2,7]

    styles = {'non-lin':'','lin':'to','ADI':'ta'}
    colours = ['255/0/0','0/255/0','0/0/255','0/255/255','255/0/255','255/255/0','127/0/0','0/127/0','0/0/127','0/127/127','127/0/127','127/127/0']

    done_grid = []
    # loop over solvers
    sn = 0
    for s in verif[exp]:
        # loop over grid sizes
        s_keyarea.plot_line(s,'3/0/0/0%s'%styles[s])
        for g in verif[exp][s]:
            # get time steps
            dt = verif[exp][s][g].keys()
            dt.sort()

            # get errors
            error_dome = []
            error_max = []
            for t in dt:
                error_dome.append(verif[exp][s][g][t][0])
                error_max.append(verif[exp][s][g][t][1])
            # plot line
            area_dome.line('-W3/%s%s'%(colours[grids[g]],styles[s]),dt,error_dome)
            area_dome.plotsymbol(dt,error_dome,size=0.1,args='-W1/%s'%(colours[grids[g]]))

            area_max.line('-W3/%s%s'%(colours[grids[g]],styles[s]),dt,error_max)
            area_max.plotsymbol(dt,error_max,size=0.1,args='-W1/%s'%(colours[grids[g]]))

            if g not in done_grid:
                e_keyarea.plot_line('@~D@~x=%skm'%g,'3/%s'%(colours[grids[g]]))
                done_grid.append(g)

            gn = gn +1
        sn = sn+1
    

    area_dome.finalise()
    area_dome.coordsystem()

    area_max.finalise()
    area_max.coordsystem()    

    plot.close()
