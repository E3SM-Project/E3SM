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
    if var_time is None:
        # Climo cannot be ran on this variable.
        return var

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
        cycle = ['DJF','MAM','JJA','SON']
    else:
        cycle = [season]

    ncycle = len(cycle)
    climo = ma.zeros([ncycle]+list(np.shape(v))[1:])
    for n in range(ncycle):
        idx = np.array( [ season_idx[cycle[n]][var_time_absolute[i].month-1]
                          for i in range(len(var_time_absolute)) ], dtype=np.int).nonzero()
        climo[n] = ma.average(v[idx], axis=0, weights=dt[idx])


    trans_var = cdms2.createVariable(climo)(squeeze=1)
    # Losing the grid after a squeeze is normal, we need to set it again.
    trans_var.setGrid(var.getGrid())
    # Set the correct axis as well.
    trans_var.setAxis(0, var.getAxis(1))
    trans_var.setAxis(1, var.getAxis(2))
    if var.getLevel():
        # If it's a 3D variable, set the last axis.
        trans_var.setAxis(2, var.getAxis(3))
    
    # Copy any missing attributes from var to trans_var.
    # The below doesn't work, since MaskedArrays apparently
    # can't be set as attributes from another object. 
    # for attr in dir(var):
    #     if not attr.startswith('__'):
    #         old_attr = getattr(var, attr)
    #         setattr(trans_var, attr, old_attr)

    trans_var.attributes = var.attributes
    trans_var.units = var.units
    trans_var.id = var.id
    trans_var.long_name = var.long_name

    return trans_var

