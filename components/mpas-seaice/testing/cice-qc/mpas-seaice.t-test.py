#!/usr/bin/env python
'''
This script performs the t-test validation for non-bit-for-bit results for the
CICE model.

Written by: Matthew Turner
Date: October, 2017

Adapted for MPAS-Seaice
Darin Comeau
September 2022
'''

import os
import sys
import logging
import glob
import numpy as np
import numpy.ma as ma
import netCDF4 as nc

import re
import fnmatch

def maenumerate(marr):
    '''
    This function provides the enumerate functionality for masked arrays
    '''
    mask = ~marr.mask.ravel()
    try:   # Python 2
        import itertools
        for i, m in itertools.izip(np.ndindex(marr.shape[-2:]), mask):
            if m: yield i
    except:  # Python 3
        for i, m in zip(np.ndindex(marr.shape[-2:]), mask):
            if m: yield i

def gen_filenames(base_dir, test_dir):
    '''
    This function is passed the directories of the baseline and test history
    files and generates a list of filenames for each.
    '''
    # The path to output files for simulation 'a' (the '-bc' simulation)
    path_a = base_dir

    # The path to output files for simulation 'b' (the test simulation)
    path_b = test_dir

    # Find the number of output files to be read in
    files_a = fnmatch.filter(os.listdir(path_a+'/'), '*mpassi.hist.000*')
    files_b = fnmatch.filter(os.listdir(path_b+'/'), '*mpassi.hist.000*')


    if not len(files_a) == len(files_b):
        logger.error("Number of output files for baseline simulation does not match the number" + \
              " of files for the test simulation.  Exiting...\n" + \
              "Baseline directory: {}\n".format(path_a) + \
              "   # of files: {}\n".format(len(files_a)) + \
              "Test directory: {}\n".format(path_b) + \
              "   # of files: {}".format(len(files_b)))
        sys.exit(-1)

    if len(files_a) < 60:
        logger.error("Number of output files too small, expecting at least 60." + \
              " Exiting...\n" + \
              "Baseline directory: {}\n".format(path_a) + \
              "   # of files: {}\n".format(len(files_a)) + \
              "Test directory: {}\n".format(path_b) + \
              "   # of files: {}".format(len(files_b)))
        sys.exit(-1)

    logger.info("Number of files: %d", len(files_a))

    return path_a, path_b, files_a, files_b

def get_geom(path, file):
    '''
    This function reads the latCell, lonCell variables from a netcdf file
    '''
    fid = nc.Dataset("{}/{}".format(path, file), 'r')
    latCell = fid.variables['latCell'][:]
    lonCell = fid.variables['lonCell'][:]
    fid.close()

    return latCell, lonCell

def read_data(path_a, path_b, files_a, files_b):
    '''
    Read the baseline and test data for sea ice thickness.  The calculate
    the difference for all locations where sea ice thickness is greater
    than 0.01 meters.
    '''
    def fill_data_array(path, files):
        '''Function to fill the data arrays'''
        # # Initialize the data array
        data = []
        cnt = 0
        logging.info("Reading data")
        for fname in sorted(files):
            nfid = nc.Dataset("{}/{}".format(path, fname), 'r')
            tmp = nfid.variables[var][:][:]
            if cnt == 0:
                data = tmp
            else:
                data = np.append(data, tmp, axis=0)
            cnt += 1
            nfid.close()

        return data # (nDays, nCells)

    def calc_diff(data_a, data_b):
        '''
        Calculate the difference and mask the points where the the difference at
        every timestep is 0, or the sea ice thickness for every timestep is < 0.01 meters
        for data_a or data_b
        '''
        logging.info("Diffing data")
        data_d = data_a - data_b
        mask_d = np.logical_or(\
                      np.logical_or(\
                           np.all(np.equal(data_d, 0.), axis=0), np.all(data_a < 0.01, axis=0))\
                      , np.all(data_b < 0.01, axis=0))
        mask_array_a = np.zeros_like(data_d)

        for x, value in np.ndenumerate(mask_d):
            iCell = x
            mask_array_a[:, iCell] = value
        del mask_d

        data_a = ma.masked_array(data_a, mask=mask_array_a)
        data_b = ma.masked_array(data_b, mask=mask_array_a)
        data_d = ma.masked_array(data_d, mask=mask_array_a)

        del mask_array_a

        return data_a, data_b, data_d

    var = 'iceThicknessCell'

    data_a = fill_data_array(path_a, files_a)
    data_b = fill_data_array(path_b, files_b)

    data_a, data_b, data_d = calc_diff(data_a, data_b)

    return data_a, data_b, data_d

def two_stage_test(data_a, num_samp, data_d, fname, path):
    '''
    This function performs the Two-Stage Paired Thickness Test
    '''
    def stage_one(data_d, num_samp, mean_d, variance_d):
        logging.info("Testing data")
        logger.debug('Running step 1 of 2-stage test')

        # Calculate the mean from 1:end-1 and 2:end
        mean_nm1_d = np.mean(data_d[:-1, :], axis=0)
        mean_2n_d = np.mean(data_d[1:, :], axis=0)

        # Calculate equation (5) for both simulations
        r1_num = np.zeros_like(mean_d)
        r1_den1 = np.zeros_like(mean_d)
        r1_den2 = np.zeros_like(mean_d)
        for i in np.arange(np.size(data_a, axis=0)-1):
            r1_num = r1_num + (data_d[i, :]-mean_nm1_d[:])*(data_d[i+1, :]-mean_2n_d[:])
            r1_den1 = r1_den1 + np.square(data_d[i, :]-mean_nm1_d[:])

        for i in np.arange(1, np.size(data_a, axis=0)):
            r1_den2 = r1_den2 + np.square(data_d[i, :] - mean_2n_d[:])

        r1 = r1_num / np.sqrt(r1_den1*r1_den2)

        # Calculate the effective sample size
        n_eff = num_samp * ((1.-r1) / (1.+r1))
        n_eff[n_eff < 2] = 2
        n_eff[n_eff > num_samp] = num_samp

        # Calculate the t-statistic with n_eff
        t_val = mean_d / np.sqrt(variance_d / n_eff)

        # Effective degrees of freedom
        df = n_eff - 1

        # Read in t_crit table
        if os.path.exists('/lcrc/group/e3sm/public_html/inputdata/ice/mpas-seaice/general/CICE_t_critical_p0.8.nc'):
          nfid = nc.Dataset("/lcrc/group/e3sm/public_html/inputdata/ice/mpas-seaice/general/CICE_t_critical_p0.8.nc", 'r')
        else:
          logger.error("Cannot locate file CICE_t_critical_p0.8.nc")
        df_table = nfid.variables['df'][:]
        t_crit_table = nfid.variables['tcrit'][:]
        nfid.close()
        t_crit = np.zeros_like(t_val)

        # Calculate critical t-value for each grid cell, based on the t_crit table
        data_d = np.squeeze(data_d[:,0])
        for x in maenumerate(data_d):
            min_val = np.min(np.abs(df[x]-df_table))
            idx = np.where(np.abs(df[x]-df_table) == min_val)
            # Handle the cases where the data point falls exactly half way between
            # 2 critical T-values (i.e., idx has more than 1 value in it)
            while True:
                try:
                    idx = idx[0]
                except:
                    break
            t_crit[x] = t_crit_table[idx]

        # Create an array of Pass / Fail values for each grid cell
        H1 = np.abs(t_val) > t_crit

        return n_eff, H1, r1, t_crit

    # Calculate the mean of the difference
    mean_d = np.mean(data_d, axis=0)

    # Loop through each timestep and calculate the square of the difference.
    # This is required (instead of just np.square(data_d - mean_d) to reduce
    # the memory footprint of the script.
    tmp1 = np.zeros_like(data_d)
    for i in np.arange(np.shape(data_d)[0]):
        tmp1[i, :] = np.square(data_d[i, :] - mean_d)
    variance_d = np.sum(tmp1) / float(num_samp - 1)

    n_eff, H1, r1, t_crit = stage_one(data_d, num_samp, mean_d, variance_d)

    if np.all(H1 == False) and np.all(n_eff >= 30):
        # H0 confirmed in all cells, and all effective sample size >= 30
        logger.debug('H0 confirmed in all cells, and all effective sample size >= 30')
        logger.info('2 stage test passed')
        return True, H1

    elif np.all(H1):
        # H1 in every grid cell
        logger.debug('H1 in all cells')
        logger.info('2 stage test failed')
        return False, H1

########### H0 confirmed for some grid cells with n_eff < 30 ############
    logger.debug('Number of H1 grid cells after stage 1 = %d', np.sum(H1))

    logger.debug('Running step 2 of 2-stage test')

    # Find the indices where n_eff is less than 30, and H0 is confirmed
    tmp_idx = np.where(n_eff < 30) and np.where(H1 == False)

    # Calculate the T-statistic using actual sample size
    t_val = mean_d / np.sqrt(variance_d / num_samp)

    # Find t_crit from the nearest value on the Lookup Table Test
    if os.path.exists('/lcrc/group/e3sm/public_html/inputdata/ice/mpas-seaice/general/CICE_Lookup_Table_p0.8_n1825.nc'):
      nfid = nc.Dataset("/lcrc/group/e3sm/public_html/inputdata/ice/mpas-seaice/general/CICE_Lookup_Table_p0.8_n1825.nc", 'r')
    else:
       logger.error("Cannot locate file CICE_Lookup_Table_p0.8_n1825.nc")
    r1_table = nfid.variables['r1'][:]
    t_crit_table = nfid.variables['tcrit'][:]
    nfid.close()

    # Fill t_crit based on lookup table
    data_d = np.squeeze(data_d[:,0])
    for x in maenumerate(data_d):
        min_val = np.min(np.abs(r1[x]-r1_table))
        idx = np.where(np.abs(r1[x]-r1_table) == min_val)
        # Handle the cases where the data point falls exactly half way between
        # 2 critical T-values (i.e., idx has more than 1 value in it)
        while True:
            try:
                idx = idx[0]
            except:
                break
        t_crit[x] = t_crit_table[idx]

    # Create an array showing locations of Pass / Fail grid cells
    H1[tmp_idx] = abs(t_val[tmp_idx]) > t_crit[tmp_idx]

    logger.debug('Number of H1 grid cells after stage 2 = %d', np.sum(H1))

    if np.all(H1):
        # H1 in all grid cells
        logger.debug('H1 in all cells, stage 2')
        logger.info('2 Stage Test Failed')
        return False, H1

    elif np.all(H1 == False):
        # H0 confirmed in all grid cells
        logger.debug('H0 confirmed in all cells with n_eff < 30')
        logger.info('2 Stage Test Passed')
        return True, H1

####### Some grid cells have H0 confirmed, and some do not ######

    # Calculate the area-weighted fraction of the test region that failed (f_val).
    #   If f_val is greater than or equal to the critical fraction, the test fails"
    f_val = critical_fraction(data_a, H1, fname, path)
    f_crit = 0.5
    if f_val >= f_crit:
        logger.info('2 Stage Test Failed')
        logger.debug('Area-weighted fraction of failures is greater than ' + \
                    'critical fraction.  Test failed.')
        logger.debug('Area-weighted fraction of failures = %f', f_val)
        return False, H1
    else:
        logger.info('2 Stage Test Passed')
        logger.debug('Area-weighted fraction of failures = %f', f_val)
        return True, H1

def critical_fraction(data_a, failures, fname, path_a):
    '''
    This function calculates the area-weighted average of cells where H1 is true.
    '''
    logger.debug('Calculating area-weighted average of H1 cells')
    # First calculate the weight attributed to each grid point (based on Area)
    nfid = nc.Dataset("{}/{}".format(path_a, fname), 'r')
    tarea = nfid.variables['areaCell'][:]
    nfid.close()
    tarea = ma.masked_array(tarea, mask=data_a[0, :].mask)
    area_weight = tarea / np.sum(tarea)

    # Calculate the area weight of the failing grid cells
    weight_tot = 0
    weight_fail = 0
    data_a = np.squeeze(data_a[:,0])
    for x in maenumerate(data_a):
        weight_tot += area_weight[x]
        if failures[x]:
            weight_fail += area_weight[x]

    return weight_fail/weight_tot

def skill_test(path_a, fname, data_a, data_b, num_samp, hemisphere):
    '''Calculate Taylor Skill Score'''
    # First calculate the weight attributed to each grid point (based on Area)
    nfid = nc.Dataset("{}/{}".format(path_a, fname), 'r')
    tarea = nfid.variables['areaCell'][:]
    nfid.close()
    tarea = ma.masked_array(tarea, mask=data_a[0, :].mask)
    area_weight = tarea / np.sum(tarea)

    weighted_mean_a = 0
    weighted_mean_b = 0
    for i in np.arange(num_samp):
        weighted_mean_a = weighted_mean_a + np.sum(area_weight*data_a[i, :])
        weighted_mean_b = weighted_mean_b + np.sum(area_weight*data_b[i, :])

    weighted_mean_a = weighted_mean_a / num_samp
    weighted_mean_b = weighted_mean_b / num_samp

    nonzero_weights = np.count_nonzero(area_weight)
    area_var_a = 0
    area_var_b = 0
    for t in np.arange(num_samp):
        area_var_a = area_var_a + np.sum(area_weight*np.square(data_a[t, :]-weighted_mean_a))
        area_var_b = area_var_b + np.sum(area_weight*np.square(data_b[t, :]-weighted_mean_b))

    area_var_a = nonzero_weights / (num_samp * nonzero_weights - 1.) * area_var_a
    area_var_b = nonzero_weights / (num_samp * nonzero_weights - 1.) * area_var_b
    std_a = np.sqrt(area_var_a)
    std_b = np.sqrt(area_var_b)

    combined_cov = 0
    for i in np.arange(num_samp):
        combined_cov = combined_cov + np.sum(area_weight*(data_a[i, :]-weighted_mean_a)*\
                                            (data_b[i, :]-weighted_mean_b))

    combined_cov = nonzero_weights / (num_samp * nonzero_weights - 1.) * combined_cov

    weighted_r = combined_cov / (std_a*std_b)

    s = np.square((1+weighted_r)*(std_a*std_b)/\
                 (area_var_a + area_var_b))

    logger.debug('%s Hemisphere skill score = %f', hemisphere, s)

    s_crit = 0.99
    if s < 0 or s > 1:
        logger.error('Skill score out of range for %s Hemisphere', hemisphere)
        return False
    elif s > s_crit:
        logger.info('Quadratic Skill Test Passed for %s Hemisphere', hemisphere)
        return True
    else:
        logger.info('Quadratic Skill Test Failed for %s Hemisphere', hemisphere)
        return False

### DC plot_data not yet working, produces empty plots, needs updating

def plot_data(data, lat, lon, units, case, plot_type):
    '''This function plots MPAS-Seaice data and creates a .png file.'''

    try:
        # Load the necessary plotting libraries
        import matplotlib.pyplot as plt
        import cartopy.crs as ccrs
        import cartopy.feature as cfeature
    except ImportError:
        logger.warning('Error loading necessary Python modules in plot_data function')
        return

    # Suppress Matplotlib deprecation warnings
    import warnings
    warnings.filterwarnings("ignore", category=UserWarning)

    # define north and south polar stereographic coord ref system
    npstereo = ccrs.NorthPolarStereo(central_longitude=-90.0) # define projection
    spstereo = ccrs.SouthPolarStereo(central_longitude= 90.0) # define projection
   
    # define figure
    fig = plt.figure(figsize=[14,7])  

    # add axis for each hemishpere
    ax1 = fig.add_subplot(121,projection=npstereo)
    ax2 = fig.add_subplot(122,projection=spstereo)   

    # set plot extents
    ax1.set_extent([-180.,180.,35.,90.],ccrs.PlateCarree())
    ax2.set_extent([-180.,180.,-90.,-35.],ccrs.PlateCarree())
   
    # add land features NH plot
    ax1.add_feature(cfeature.LAND, color='lightgray')
    ax1.add_feature(cfeature.BORDERS)
    ax1.add_feature(cfeature.COASTLINE)

    # add land features SH plot
    ax2.add_feature(cfeature.LAND, color='lightgray')
    ax2.add_feature(cfeature.BORDERS)
    ax2.add_feature(cfeature.COASTLINE)

    # add grid lines
    dlon = 30.0
    dlat = 15.0
    mpLons = np.arange(-180.  ,180.0+dlon,dlon)
    mpLats = np.arange(-90.,90.0+dlat ,dlat)

    g1 = ax1.gridlines(xlocs=mpLons,ylocs=mpLats,
                       draw_labels=True,
                       x_inline=False,y_inline=False)

    g2 = ax2.gridlines(xlocs=mpLons,ylocs=mpLats,
                       draw_labels=True,
                       x_inline=False,y_inline=False)


    # Specify Min/max colors for each hemisphere
    # check for minus to see if it is a difference plot
    if '\n- ' in case: # this is a difference plot
        # specify colormap
        mycmap = 'seismic' # blue,white,red with white centered colormap

        # determine max absolute value to use for color range
        # intent is use same min/max with center zero
        dmin = np.abs(data.min())
        dmax = np.abs(data.max())
        clim = np.max([dmin,dmax])

        # this specifies both hemishperes the same range.
        cminNH = -clim
        cmaxNH =  clim
        cminSH = -clim
        cmaxSH =  clim

    else:  # not a difference plot
        # specify colormap
        mycmap = 'RdBu'

        # arbitrary limits for each Hemishpere
        cminNH = 0.0
        cmaxNH = 5.0
        cminSH = 0.0
        cmaxSH = 2.0

    if plot_type == 'scatter':
        # plot NH
        scNH = ax1.scatter(lon,lat,c=data,cmap=mycmap,s=4,edgecolors='none',
                           vmin=cminNH, vmax=cmaxNH,
                           transform=ccrs.PlateCarree())
        
        # plot SH
        scSH = ax2.scatter(lon,lat,c=data,cmap=mycmap,s=4,edgecolors='none',
                           vmin=cminSH, vmax=cmaxSH,
                           transform=ccrs.PlateCarree())

    else: 
        if plot_type == 'contour':
            print("contour plot depreciated. using pcolor.")
            
        scNH = ax1.pcolormesh(lon,lat,data,cmap=mycmap,
                              vmin=cminNH, vmax=cmaxNH,
                              transform=ccrs.PlateCarree())
        
        scSH = ax2.pcolormesh(lon,lat,data,cmap=mycmap,
                              vmin=cminSH, vmax=cmaxSH,
                              transform=ccrs.PlateCarree())

    plt.suptitle(f'MPAS-Seaice Mean Ice Thickness\n{case:s}')

    # add more whitespace between plots for colorbar. 
    plt.subplots_adjust(wspace=0.4)

    # add separate axes for colorbars
    # first get position/size of current axes
    pos1 = ax1.get_position()
    pos2 = ax2.get_position()

    # now add new colormap axes using the position ax1, ax2 as reference
    cax1 = fig.add_axes([pos1.x0+pos1.width+0.03,
                         pos1.y0,
                         0.02,
                         pos1.height])
 
    cax2 = fig.add_axes([pos2.x0+pos2.width+0.03,
                         pos2.y0,
                         0.02,
                         pos2.height])


    if '\n- ' in case:
        # If making a difference plot, use scientific notation for colorbar
        cbNH = plt.colorbar(scNH, cax=cax1, orientation="vertical", 
                            pad=0.1, format="%.1e")
        cbSH = plt.colorbar(scSH, cax=cax2, orientation="vertical", 
                            pad=0.1, format="%.1e")

    else:
        #pass
        # If plotting non-difference data, do not use scientific notation for colorbar
        cbNH = plt.colorbar(scNH, cax=cax1, orientation="vertical", 
                            pad=0.1, format="%.2f")
        cbSH = plt.colorbar(scSH, cax=cax2, orientation="vertical", 
                            pad=0.1, format="%.2f")

    cbNH.set_label(units, loc='center')
    cbSH.set_label(units, loc='center')

    outfile = 'ice_thickness_{}.png'.format(case.replace('\n- ','_minus_'))
    logger.info('Creating map of the data ({})'.format(outfile))
    plt.savefig(outfile, dpi=300, bbox_inches='tight')

def plot_two_stage_failures(data, lat, lon):
    '''This function plots each grid cell and whether or not it Passed or Failed
       the two-stage test.  It then either creates a .png file
       (two_stage_test_failure_map.png), or saves the failure locations to a
       text file.
    '''

    # Convert the boolean array (data) to an integer array
    int_data = data.astype(int)

    try:
        logger.info('Creating map of the failures (two_stage_test_failure_map.png)')
        # Load the necessary plotting libraries
        import matplotlib.pyplot as plt
        import cartopy.crs as ccrs
        import cartopy.feature as cfeature
        from mpl_toolkits.axes_grid1 import make_axes_locatable
        from matplotlib.colors import LinearSegmentedColormap

        # Suppress Matplotlib deprecation warnings
        import warnings
        warnings.filterwarnings("ignore", category=UserWarning)

        # Create the figure
        fig = plt.figure(figsize=(12, 8))
        
        # define plot projection and create axis
        pltprj = ccrs.Mollweide(central_longitude=0.0)
        ax = fig.add_subplot(111,projection=pltprj)
        
        # add land
        ax.add_feature(cfeature.LAND, color='lightgray')
        ax.add_feature(cfeature.BORDERS)
        ax.add_feature(cfeature.COASTLINE)
        #gshhs = cfeature.GSHHSFeature(scale='auto',facecolor='lightgray',edgecolor='none')
        #ax.add_feature(gshhs)

        # Create the custom colormap
        colors = [(0, 0, 1), (1, 0, 0)]  # Blue, Red
        cmap_name = 'RB_2bins'
        cm = LinearSegmentedColormap.from_list(cmap_name, colors, N=2)

        # Plot the data as a scatter plot
        sc = ax.scatter(lon,lat,c=int_data,cmap=cm,s=4,lw=0,
                        vmin=0.,vmax=1.,
                        transform=ccrs.PlateCarree())

        # add grid lines
        dlon = 60.0
        dlat = 30.0
        mpLons = np.arange(-180.  ,180.0+dlon,dlon)
        mpLats = np.arange(-90.,90.0+dlat ,dlat)
        mpLabels = {"left":   "y",
                    "bottom": "x"}
    
        ax.gridlines(xlocs=mpLons,ylocs=mpLats,
                     draw_labels=mpLabels)

        plt.title('MPAS-Seaice Two-Stage Test Failures')

        # Create the colorbar and add Pass / Fail labels
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("bottom", size="5%", pad=0.5)
        cb = plt.colorbar(sc, cax=cax, orientation="horizontal", format="%.0f")
        cb.set_ticks([])
        cb.ax.text(-0.01, -0.5, 'PASS')
        cb.ax.text(0.99, -0.5, 'FAIL')

        plt.savefig('two_stage_test_failure_map.png', dpi=300)
    except:
        logger.warning('')
        logger.warning('Unable to plot the data.  Saving latitude and longitude')
        logger.warning('for ONLY failures to two_stage_failure_locations.txt')

        # Create a file and write the failures only to the file
        f = open('two_stage_failure_locations.txt', 'w')
        f.write('# CICE Two-stage test failures\n')
        f.write('# Longitude,Latitude\n')
        for i in range(data.shape[0]):
            for j in range(data.shape[1]):
                if (not data.mask[i, j]) and data[i, j]:
                    f.write('{},{}\n'.format(lon[i, j], lat[i, j]))

        f.close()

def main():
    import argparse
    parser = argparse.ArgumentParser(description='This script performs the T-test for \
                           CICE simulations that should be bit-for-bit, but are not.')
    parser.add_argument('base_dir', \
                help='Path to the baseline history (iceh_inst*) files.  REQUIRED')
    parser.add_argument('test_dir', \
                help='Path to the test history (iceh_inst*) files.  REQUIRED')
    parser.add_argument('-v', '--verbose', dest='verbose', help='Print debug output?', \
                        action='store_true')
    parser.add_argument('-pt','--plot_type', dest='plot_type', help='Specify type of plot \
                        to create', choices=['scatter','contour','pcolor'])

    parser.set_defaults(verbose=False)
    parser.set_defaults(plot_type='pcolor')

    # If no arguments are provided, print the help message
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args()

    # Set up the logger
    global logger
    if args.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)
    # Log to log file as well as stdout
    fh = logging.FileHandler(r'qc_log.txt', 'w')
    logger = logging.getLogger(__name__)
    logger.addHandler(fh)

    logger.info('Running QC test on the following directories:')
    logger.info('  {}'.format(args.base_dir))
    logger.info('  {}'.format(args.test_dir))

    dir_a, dir_b, files_base, files_test = gen_filenames(args.base_dir, args.test_dir)

    t_lat, t_lon = get_geom(dir_a, files_base[0])

    data_base, data_test, data_diff = read_data(dir_a, dir_b, files_base, files_test)

    num_samp, nCells = np.shape(data_base)

    if np.ma.all(data_diff.mask):
        logger.info("Data is bit-for-bit.  No need to run QC test")
        sys.exit(0)

    # Run the two-stage test
    PASSED, H1_array = two_stage_test(data_base, num_samp, data_diff, files_base[0], dir_a)

    # Delete arrays that are no longer necessary
    del data_diff

    # If test failed, attempt to create a plot of the failure locations
    if not PASSED:
        plot_two_stage_failures(H1_array, t_lat, t_lon)
        
        # # Create plots of mean ice thickness
        baseDir = os.path.abspath(args.base_dir).rstrip('run/').rstrip(\
                                                        'run').split('/')[-1]
        testDir = os.path.abspath(args.test_dir).rstrip('run/').rstrip( \
                                                        'run').split('/')[-1]
        plot_data(np.mean(data_base,axis=0), t_lat, t_lon, 'm', baseDir, args.plot_type)
        plot_data(np.mean(data_test,axis=0), t_lat, t_lon, 'm', testDir, args.plot_type)
        plot_data(np.mean(data_base-data_test,axis=0), t_lat, t_lon, 'm', '{}\n- {}'.\
                  format(baseDir,testDir), args.plot_type)

        logger.error('Quality Control Test FAILED')
        sys.exit(-1)

    # Create a northern hemisphere and southern hemisphere mask
    mask_tlat = t_lat < 0
    mask_nh = np.zeros_like(data_base)
    mask_sh = np.zeros_like(data_base)
    for x, val in np.ndenumerate(mask_tlat):
        mask_nh[:, x] = val
        mask_sh[:, x] = not val    

    # Run skill test on northern hemisphere
    data_nh_a = ma.masked_array(data_base, mask=mask_nh)
    data_nh_b = ma.masked_array(data_test, mask=mask_nh)
    if np.ma.all(data_nh_a.mask) and np.ma.all(data_nh_b.mask):
        logger.info("Northern Hemisphere data is bit-for-bit")
        PASSED_NH = True
    else:
        PASSED_NH = skill_test(dir_a, files_base[0], data_nh_a, data_nh_b, num_samp, 'Northern')

    # Run skill test on southern hemisphere
    data_sh_a = ma.masked_array(data_base, mask=mask_sh)
    data_sh_b = ma.masked_array(data_test, mask=mask_sh)
    if np.ma.all(data_sh_a.mask) and np.ma.all(data_sh_b.mask):
        logger.info("Southern Hemisphere data is bit-for-bit")
        PASSED_SH = True
    else:
        PASSED_SH = skill_test(dir_a, files_base[0], data_sh_a, data_sh_b, num_samp, 'Southern')

    PASSED_SKILL = PASSED_NH and PASSED_SH

    # Plot the ice thickness data for the base and test cases
    baseDir = os.path.abspath(args.base_dir).rstrip('run/').rstrip( \
                                                    'run').split('/')[-1]
    testDir = os.path.abspath(args.test_dir).rstrip('run/').rstrip( \
                                                    'run').split('/')[-1]
    plot_data(np.mean(data_base,axis=0), t_lat, t_lon, 'm', baseDir, args.plot_type)
    plot_data(np.mean(data_test,axis=0), t_lat, t_lon, 'm', testDir, args.plot_type)
    plot_data(np.mean(data_base-data_test,axis=0), t_lat, t_lon, 'm', '{}\n- {}'.\
              format(baseDir,testDir), args.plot_type)

    logger.info('')
    if not PASSED_SKILL:
        logger.error('Quality Control Test FAILED')
        sys.exit(1)  # exit with an error return code
    else:
        logger.info('Quality Control Test PASSED')
        sys.exit(0)  # exit with successfull return code

if __name__ == "__main__":
    main()
