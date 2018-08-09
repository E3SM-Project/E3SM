import time
import numpy as np
import numpy.ma as ma


def climo(var, season):
    """
    Compute the climatology for var for the given season.
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

    return climo
