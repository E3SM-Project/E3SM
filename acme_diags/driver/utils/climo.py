import copy
import numpy as np
import numpy.ma as ma
import cdms2


def climo(var, season):
    """
    Compute the climatology for var for the given season.
    The returned variable must be 2 dimensional.
    """
    season_idx = {
        '01': [1,0,0,0,0,0,0,0,0,0,0,0],
        '02': [0,1,0,0,0,0,0,0,0,0,0,0],
        '03': [0,0,1,0,0,0,0,0,0,0,0,0],
        '04': [0,0,0,1,0,0,0,0,0,0,0,0],
        '05': [0,0,0,0,1,0,0,0,0,0,0,0],
        '06': [0,0,0,0,0,1,0,0,0,0,0,0],
        '07': [0,0,0,0,0,0,1,0,0,0,0,0],
        '08': [0,0,0,0,0,0,0,1,0,0,0,0],
        '09': [0,0,0,0,0,0,0,0,1,0,0,0],
        '10': [0,0,0,0,0,0,0,0,0,1,0,0],
        '11': [0,0,0,0,0,0,0,0,0,0,1,0],
        '12': [0,0,0,0,0,0,0,0,0,0,0,1],
        'DJF':[1,1,0,0,0,0,0,0,0,0,0,1],
        'MAM':[0,0,1,1,1,0,0,0,0,0,0,0],
        'JJA':[0,0,0,0,0,1,1,1,0,0,0,0],
        'SON':[0,0,0,0,0,0,0,0,1,1,1,0],
        'ANN':[1,1,1,1,1,1,1,1,1,1,1,1],
    }

    # Redefine time to be in the middle of the time interval
    var_time = var.getTime()
    tbounds = var_time.getBounds()
    var_time[:] = 0.5*(tbounds[:,0]+tbounds[:,1])
    var_time_absolute = var_time.asComponentTime()

    # Compute time length
    dt = tbounds[:,1] - tbounds[:,0]

    # Convert to masked array
    v = var.asma()

    # Compute climatology
    if season == 'ANNUALCYCLE':
        cycle = ['01','02','03','04','05','06','07','08','09','10','11','12']
    elif season == 'SEASONALCYCLE':
        cycle = ['DJF','MAM','JJA', 'SON']
    else:
        cycle = [ season ]

    ncycle = len(cycle)
    climo = ma.zeros([ncycle]+list(np.shape(v))[1:])
    for n in range(ncycle):
        idx = np.array( [ season_idx[cycle[n]][var_time_absolute[i].month-1]
                          for i in range(len(var_time_absolute)) ], dtype=np.int).nonzero()
        climo[n] = ma.average(v[idx], axis=0, weights=dt[idx])

    print('in climo.py')

    temp = cdms2.asVariable(climo)(squeeze=1)

    temp.attributes = var.attributes
    temp.units = var.units
    temp.setGrid(var.getGrid())
    temp.setAxis(0, var.getAxis(1))
    temp.setAxis(1, var.getAxis(2))
    return temp

    """
    t_var = copy.deepcopy(var)
    t_var = t_var(squeeze=1)  # Doesn't do anything b/c shape is (6000, _, _)
    #climo = climo(squeeze=1)
    #t_var.data = climo.data
    stuff = cdms2.asVariable(climo)(squeeze=1)
    print(stuff.shape)
    print(t_var._data.shape)
    print(stuff._data.shape)
    t_var._data = stuff._data
    return t_var

    #print(type(var))
    #print('dir(climo)', dir(climo))
    print(var.data.shape)
    print(climo.data.shape)
    var[:,] = climo
    var = var[:1]
    #var.data = climo.data
    #quit()
    return var
    """
    '''
    t_var = cdms2.asVariable(climo)
    print(dir(t_var))
    quit()
    t_var = t_var(squeeze=1)
    #print(type(t_var))
    #print(dir(t_var))
    t_var.attributes = var.attributes
    t_var.units = var.units
    print('t_var.getGrid()', t_var.getGrid())
    t_var.setGrid(var.getGrid())
    print('t_var.getGrid()', t_var.getGrid())
    #print(type(t_var))
    #print(type(t_var(squeeze=1)))
    # return t_var(squeeze=1)
    return t_var
    '''
