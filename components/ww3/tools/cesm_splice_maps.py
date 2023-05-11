#!/usr/bin/python

# netcdf4 reader for cesm remapping file

###
# creates land/ocean (red/blue) maps of src grid and dst grid
# creates map of what points are being interpolated from (red/blue)
# creates map of what points are being interpolated to (red/blue)
# needs map of what points are missing in destination (compare blin_dst_mask and blin_row_mask)?
# needs map of what points are being included that shouldn't be (compare blin_src_mask and blin_col_mask)?

import netCDF4  # for netCDF
from mpl_toolkits.basemap import Basemap  # mapping toolkit
import matplotlib.pyplot as plt  # shorthand mapping
import matplotlib.gridspec as gridspec # n x m figures
import numpy as np # num num numpy
import math
import glob
import os
import shutil # file operations

os.chdir('/Users/arbetter/Desktop/python/cesm_maps/RemappingFiles/')
os.getcwd()


cesm_remap_blin_file = ['map_gx1v6_TO_ww3a_blin.150511.nc','map_ww3a_TO_gx1v6_blin.150511.nc']
cesm_remap_stod_file = ['map_gx1v6_TO_ww3a_neareststod.150511.nc','map_ww3a_TO_gx1v6_neareststod.150511.nc']
cesm_remap_splice_file = ['map_gx1v6_TO_ww3a_splice.nc','map_ww3a_TO_gx1v6_splice.nc']

#cesm_remap_blin_file = ['map_gx1v6_TO_ww3a_blin.150430.nc']
#cesm_remap_stod_file = ['map_gx1v6_TO_ww3a_neareststod.150430.nc']
#cesm_remap_splice_file = ['map_gx1v6_TO_ww3a_splice.nc']

#cesm_remap_blin_file = ['map_ww3a_TO_gx3v7_blin.150311.nc']
#cesm_remap_stod_file = ['map_ww3a_TO_gx3v7_neareststod.150311.nc']
#cesm_remap_splice_file = ['map_ww3a_TO_gx3v7_splice.nc']


nfiles = len(cesm_remap_blin_file)

showplots = 0  # set to 0 to skip plotting to the screen

if nfiles == 1:
    cleanup = False
else:
    cleanup = True


print '#files =',nfiles
print 'cleanup =',cleanup


for k in range(0,nfiles):

### begin load blin block

#for k in range(0,-1):
    print k, cesm_remap_blin_file[k]

# get the infrom from grid map
    blingroup = netCDF4.Dataset(cesm_remap_blin_file[k],'r')
    print blingroup.data_model

# Source variables for blin
    blin_src_griddims = blingroup.variables['src_grid_dims'][:]
    nx_src = blin_src_griddims[0]
    ny_src = blin_src_griddims[1]
    blin_src_mask = blingroup.variables['mask_a'][:]
    blin_src_lon = blingroup.variables['xc_a'][:]
    blin_src_lat = blingroup.variables['yc_a'][:]

# DST varaiables: blin

    blin_dst_griddims = blingroup.variables['dst_grid_dims'][:]
    nx_dst = blin_dst_griddims[0]
    ny_dst = blin_dst_griddims[1]
    blin_dst_mask = blingroup.variables['mask_b'][:]
    blin_dst_lon  = blingroup.variables['xc_b'][:]
    blin_dst_lat = blingroup.variables['yc_b'][:]

# the blinear weighting function
    s_n_blin = len(blingroup.dimensions['n_s'])
    src_n_blin = len(blingroup.dimensions['n_a'])
    dst_n_blin = len(blingroup.dimensions['n_b'])

# subtract 1 from row, col so index goes from 0 to nx*ny-1
    s_row_blin = blingroup.variables['row'][:] -1 # row is the 1D index of the dst array
    s_col_blin = blingroup.variables['col'][:] -1 # col is the 1D index of the src array
    S_weight_blin = blingroup.variables['S'][:]

# close the grid file to avoid over-writing
    blingroup.close()


    if ('map_gx' in cesm_remap_blin_file[k]):  # gx3v7 in radians
        print 'gxin radians'
        radlat = 180.*blin_src_lat[:]/math.pi
        blin_src_lat = radlat
        radlon = 180*blin_src_lon[:]/math.pi
        blin_src_lon = radlon

    if ('TO_gx' in cesm_remap_blin_file[k]):  # gx3v7 in radians
        print 'gx in radians'
        radlat = 180.*blin_dst_lat[:]/math.pi
        blin_dst_lat = radlat
        radlon = 180*blin_dst_lon[:]/math.pi
        blin_dst_lon = radlon

    print cesm_remap_blin_file[k]
    print 'blin_src_grid_dims', blin_src_griddims, src_n_blin
    print 'blin_dst_grid_dims', blin_dst_griddims, dst_n_blin
    print '#points used:', s_n_blin

# col is the source grid
    blin_col_lon = blin_src_lon[s_col_blin[:]]
    blin_col_lat = blin_src_lat[s_col_blin[:]]
# row is the destination grid
    blin_row_lon = blin_dst_lon[s_row_blin[:]]
    blin_row_lat = blin_dst_lat[s_row_blin[:]]

# end load blin block



# begin load stod block


#for k in range(0,-1):
    print k, cesm_remap_stod_file[k]

# get the infrom from grid map
    stodgroup = netCDF4.Dataset(cesm_remap_stod_file[k],'r')
    print stodgroup.data_model

# Source variables for stod
    stod_src_griddims = stodgroup.variables['src_grid_dims'][:]
    nx_src = stod_src_griddims[0]
    ny_src = stod_src_griddims[1]
    stod_src_mask = stodgroup.variables['mask_a'][:]
    stod_src_lon = stodgroup.variables['xc_a'][:]
    stod_src_lat = stodgroup.variables['yc_a'][:]

# DST varaiables: stod

    stod_dst_griddims = stodgroup.variables['dst_grid_dims'][:]
    nx_dst = stod_dst_griddims[0]
    ny_dst = stod_dst_griddims[1]
    stod_dst_mask = stodgroup.variables['mask_b'][:]
    stod_dst_lon  = stodgroup.variables['xc_b'][:]
    stod_dst_lat = stodgroup.variables['yc_b'][:]

# the stodear weighting function
    s_n_stod = len(stodgroup.dimensions['n_s'])
    src_n_stod = len(stodgroup.dimensions['n_a'])
    dst_n_stod = len(stodgroup.dimensions['n_b'])

# subtract 1 from row, col so index goes from 0 to nx*ny-1
    s_row_stod = stodgroup.variables['row'][:] -1 # row is the 1D index of the dst array
    s_col_stod = stodgroup.variables['col'][:] -1 # col is the 1D index of the src array
    S_weight_stod = stodgroup.variables['S'][:]

# close the grid file to avoid over-writing
    stodgroup.close()


    if ('map_gx' in cesm_remap_stod_file[k]):  # gx3v7 in radians
        print 'gx in radians'
        radlat = 180.*stod_src_lat[:]/math.pi
        stod_src_lat = radlat
        radlon = 180*stod_src_lon[:]/math.pi
        stod_src_lon = radlon

    if ('TO_gx' in cesm_remap_stod_file[k]):  # gx3v7 in radians
        print 'gxin radians'
        radlat = 180.*stod_dst_lat[:]/math.pi
        stod_dst_lat = radlat
        radlon = 180*stod_dst_lon[:]/math.pi
        stod_dst_lon = radlon


    print cesm_remap_stod_file[k]
    print 'stod_src_grid_dims', stod_src_griddims, src_n_stod
    print 'stod_dst_grid_dims', stod_dst_griddims, dst_n_stod
    print '#points used:', s_n_stod

# col is the source grid
    stod_col_lon = stod_src_lon[s_col_stod[:]]
    stod_col_lat = stod_src_lat[s_col_stod[:]]
# row is the destination grid
    stod_row_lon = stod_dst_lon[s_row_stod[:]]
    stod_row_lat = stod_dst_lat[s_row_stod[:]]

#end load stod block

### Blin mapping block
# size the plot
    blinfig = plt.figure(figsize=(12,7))

### Create a 2x2 plot grid
    blinfig = gridspec.GridSpec(2,3)

### Create a map
    blinmap = Basemap(llcrnrlon=0,llcrnrlat=-90,urcrnrlon=360,urcrnrlat=90,projection='mill')


# using blin_src_mask as a mask, create a masked array where 0 is water (not 1)
    blin_src_omask = np.asarray(blin_src_mask)   ### 0 = land, 1 = ocn
    blin_src_ocnmask = np.flatnonzero(blin_src_omask) ### indices where ocn
    blin_src_lmask = np.asarray(1 - blin_src_mask) ### 0 = ocn, 1 = land
    blin_src_landmask = np.flatnonzero(blin_src_lmask) ### indices where land
    blin_src_allmask = np.arange(len(blin_src_mask)) # an array with all indices

    blin_src_lonocn = blin_src_lon[blin_src_ocnmask]
    blin_src_latocn = blin_src_lat[blin_src_ocnmask]
    blin_src_ocn = blin_src_omask[blin_src_ocnmask] # shouldbe an array of 1s
    blin_src_lonland = blin_src_lon[blin_src_landmask]
    blin_src_latland = blin_src_lat[blin_src_landmask]
    blin_src_land = blin_src_lmask[blin_src_landmask] # should be an array of 1s

# also do for blin_dst_mask
    blin_dst_omask = np.asarray(blin_dst_mask)
    blin_dst_ocnmask = np.flatnonzero(blin_dst_omask)
    blin_dst_lmask = np.asarray(1-blin_dst_mask)
    blin_dst_landmask = np.flatnonzero(blin_dst_lmask)
    blin_dst_allmask = np.arange(len(blin_dst_mask))

    blin_dst_lonocn = blin_dst_lon[blin_dst_ocnmask]
    blin_dst_latocn = blin_dst_lat[blin_dst_ocnmask]
    blin_dst_ocnpts = blin_dst_omask[blin_dst_ocnmask] # should be an array of 1s
    blin_dst_lonland = blin_dst_lon[blin_dst_landmask]
    blin_dst_latland = blin_dst_lat[blin_dst_landmask]
    blin_dst_landpts = blin_dst_lmask[blin_dst_landmask] #should be an array of 1s

#### now find which points are used in the mapping

#  s_col_blin specifies which src points are calcuated (col/src)
    blin_col_omask = np.asarray(blin_src_mask[s_col_blin[:]]) # 0/1 array
    blin_col_lonocn = blin_src_lon[s_col_blin[:]]*blin_src_mask[s_col_blin[:]]
    blin_col_latocn = blin_src_lat[s_col_blin[:]]*blin_src_mask[s_col_blin[:]]
    blin_col_lonland = blin_src_lon[s_col_blin[:]]*(1-blin_src_mask[s_col_blin[:]])
    blin_col_latland = blin_src_lat[s_col_blin[:]]*(1-blin_src_mask[s_col_blin[:]])

# blin_row_ocn are points calculated in the destination (row/dst)
    blin_row_omask = np.asarray(blin_dst_mask[s_row_blin[:]]) # 0/1 array
    blin_row_lonocn = blin_dst_lon[s_row_blin[:]]*blin_dst_mask[s_row_blin[:]]
    blin_row_latocn = blin_dst_lat[s_row_blin[:]]*blin_dst_mask[s_row_blin[:]]
    blin_row_lonland = blin_dst_lon[s_row_blin[:]]*(1-blin_dst_mask[s_row_blin[:]])
    blin_row_latland = blin_dst_lat[s_row_blin[:]]*(1-blin_dst_mask[s_row_blin[:]])

# compare blin_dst_ocnmask and blin_row_ocnmask (row/dst)
    blin_dst_uncalc = np.setdiff1d(blin_dst_allmask,s_row_blin)
#    blin_dst_lonuncalc = blin_dst_lon[blin_dst_uncalc[:]]*blin_dst_unmask[blin_dst_uncalc[:]]
#    blin_dst_latuncalc = blin_dst_lat[blin_dst_uncalc[:]]*blin_dst_unmask[blin_dst_uncalc[:]]
#    blin_dst_lonuncalcland = blin_dst_lon[blin_dst_uncalc[:]]*(1-blin_dst_unmask[blin_dst_uncalc[:]])
#    blin_dst_latuncalcland = blin_dst_lat[blin_dst_uncalc[:]]*(1-blin_dst_unmask[blin_dst_uncalc[:]])
    blin_dst_lonuncalc = blin_dst_lon[blin_dst_uncalc[:]]*blin_dst_mask[blin_dst_uncalc[:]]
    blin_dst_latuncalc = blin_dst_lat[blin_dst_uncalc[:]]*blin_dst_mask[blin_dst_uncalc[:]]
    blin_dst_lonuncalcland = blin_dst_lon[blin_dst_uncalc[:]]*(1-blin_dst_mask[blin_dst_uncalc[:]])
    blin_dst_latuncalcland = blin_dst_lat[blin_dst_uncalc[:]]*(1-blin_dst_mask[blin_dst_uncalc[:]])



# blin_src_unused is the points not used in blin_row_ocn (not mapped) (col/src)
    blin_src_unused = np.setdiff1d(blin_src_allmask,s_col_blin)
    blin_src_lonunused = blin_src_lon[blin_src_unused[:]]*blin_src_mask[blin_src_unused[:]]
    blin_src_latunused = blin_src_lat[blin_src_unused[:]]*blin_src_mask[blin_src_unused[:]]
    blin_src_lonunusedland = blin_src_lon[blin_src_unused[:]]*(1-blin_src_mask[blin_src_unused[:]])
    blin_src_latunusedland = blin_src_lat[blin_src_unused[:]]*(1-blin_src_mask[blin_src_unused[:]])


    blin_src_x, blin_src_y = blinmap(blin_src_lon,blin_src_lat)
    blin_dst_x, blin_dst_y = blinmap(blin_dst_lon,blin_dst_lat)
    blin_src_xocn, blin_src_yocn = blinmap(blin_src_lonocn,blin_src_latocn)
    blin_dst_xocn, blin_dst_yocn = blinmap(blin_dst_lonocn,blin_dst_latocn)
    blin_src_xland, blin_src_yland = blinmap(blin_src_lonland,blin_src_latland)
    blin_dst_xland, blin_dst_yland = blinmap(blin_dst_lonland,blin_dst_latland)
    blin_col_x, blin_col_y = blinmap(blin_col_lon,blin_col_lat)
    blin_row_x, blin_row_y = blinmap(blin_row_lon,blin_row_lat)
    blin_col_xocn, blin_col_yocn = blinmap(blin_col_lonocn,blin_col_latocn)
    blin_row_xocn, blin_row_yocn = blinmap(blin_row_lonocn,blin_row_latocn)
    blin_col_xland,blin_col_yland = blinmap(blin_col_lonland,blin_col_latland)
    blin_row_xland,blin_row_yland = blinmap(blin_row_lonland,blin_row_latland)
    blin_src_xunused,blin_src_yunused = blinmap(blin_src_lonunused,blin_src_latunused)
    blin_dst_xuncalc,blin_dst_yuncalc = blinmap(blin_dst_lonuncalc,blin_dst_latuncalc)
    blin_src_xunusedland,blin_src_yunusedland = blinmap(blin_src_lonunusedland,blin_src_latunusedland)
    blin_dst_xuncalcland,blin_dst_yuncalcland = blinmap(blin_dst_lonuncalcland,blin_dst_latuncalcland)


    print 'read stod & blin'

# here want to mark the points that are being used to create
#    blin_col_field = list()
#    blin_row_field = list()
#    blin_col_field[:] = blin_src_landmask[s_col_blin[:]]
#blin_row_field[:] = blin_dst_landmask[s_row_blin[:]]



# and a subplot for src
    srcmap = plt.subplot2grid((2,3),(0,0))

    blinmap.drawcoastlines()

# scatter markers
    blinmap.scatter(blin_src_xland,blin_src_yland,5,color='red',marker='+')
    blinmap.scatter(blin_src_xocn,blin_src_yocn,3,color='blue')

#    plt.title(cesm_remap_file[i])
    plt.title('SRC GRID')

# subplot for points used in interp
    srcfig = plt.subplot2grid((2,3),(0,1))
    blinmap.scatter(blin_col_xland,blin_col_yland,color='red',marker='+')
    blinmap.scatter(blin_col_xocn,blin_col_yocn,3,color='blue')

#   blinmap.drawcoastlines()
    plt.title('SRC IN')

# subplot for ocean points not used in interp
    unusedfig = plt.subplot2grid((2,3),(0,2))
    blinmap.scatter(blin_src_xunusedland,blin_src_yunusedland,5,color='red',marker='+')
    blinmap.scatter(blin_src_xunused,blin_src_yunused,3,color='blue')

    plt.title('SRC UNUSED')

# subplot for ocean points not calculated

    uncalcfig = plt.subplot2grid((2,3),(1,2))
    blinmap.scatter(blin_dst_xuncalcland,blin_dst_yuncalcland,5,color='red',marker='+')
    blinmap.scatter(blin_dst_xuncalc,blin_dst_yuncalc,3,color='blue')

    plt.title('DST UNCALC')

# subplot for points calculated in interp
    dstfig=plt.subplot2grid((2,3),(1,1))
#    blinmap.drawcoastlines()

    blinmap.scatter(blin_row_xland,blin_row_yland,color='red',marker='+')
    blinmap.scatter(blin_row_xocn,blin_row_yocn,3,color='blue')

#    for i in range(0,s_n):
#        if blin_row_field[i] == 1:
#            blinmap.scatter(blin_row_x[i],blin_row_y[i],5,color='red',marker='+')
#        else:
#            blinmap.scatter(blin_row_x[i],blin_row_y[i],3,color='blue')

    plt.title('DST OUT')


# and a subplot for dst
    dstmap = plt.subplot2grid((2,3),(1,0))

    blinmap.drawcoastlines()

# scatter markers
    blinmap.scatter(blin_dst_xland,blin_dst_yland,5,color='red',marker='+')
    blinmap.scatter(blin_dst_xocn,blin_dst_yocn,3,color='blue')

    plt.title('DST GRID')

    plt.suptitle(cesm_remap_blin_file[k])

    if showplots == 1:
        plt.show()
    plt.savefig(cesm_remap_blin_file[k]+'.png')

    if showplots == 0:
        plt.close()


### End Blin Map block

### StoDmapping block
# size the plot
    stodfig = plt.figure(figsize=(12,7))

### Create a 2x2 plot grid
    stodfig = gridspec.GridSpec(2,3)

### Create a map
    stodmap = Basemap(llcrnrlon=0,llcrnrlat=-90,urcrnrlon=360,urcrnrlat=90,projection='mill')


# using stod_src_mask as a mask, create a masked array where 0 is water (not 1)
    stod_src_omask = np.asarray(stod_src_mask)   ### 0 = land, 1 = ocn
    stod_src_ocnmask = np.flatnonzero(stod_src_omask) ### indices where ocn
    stod_src_lmask = np.asarray(1 - stod_src_mask) ### 0 = ocn, 1 = land
    stod_src_landmask = np.flatnonzero(stod_src_lmask) ### indices where land
    stod_src_allmask = np.arange(len(stod_src_mask)) # an array with all indices

    stod_src_lonocn = stod_src_lon[stod_src_ocnmask]
    stod_src_latocn = stod_src_lat[stod_src_ocnmask]
    stod_src_ocn = stod_src_omask[stod_src_ocnmask] # shouldbe an array of 1s
    stod_src_lonland = stod_src_lon[stod_src_landmask]
    stod_src_latland = stod_src_lat[stod_src_landmask]
    stod_src_land = stod_src_lmask[stod_src_landmask] # should be an array of 1s

# also do for stod_dst_mask
    stod_dst_omask = np.asarray(stod_dst_mask)
    stod_dst_ocnmask = np.flatnonzero(stod_dst_omask)
    stod_dst_lmask = np.asarray(1-stod_dst_mask)
    stod_dst_landmask = np.flatnonzero(stod_dst_lmask)
    stod_dst_allmask = np.arange(len(stod_dst_mask))

    stod_dst_lonocn = stod_dst_lon[stod_dst_ocnmask]
    stod_dst_latocn = stod_dst_lat[stod_dst_ocnmask]
    stod_dst_ocnpts = stod_dst_omask[stod_dst_ocnmask] # should be an array of 1s
    stod_dst_lonland = stod_dst_lon[stod_dst_landmask]
    stod_dst_latland = stod_dst_lat[stod_dst_landmask]
    stod_dst_landpts = stod_dst_lmask[stod_dst_landmask] #should be an array of 1s

#### now find which points are used in the mapping

#  s_col_stod specifies which src points are calcuated (col/src)
    stod_col_omask = np.asarray(stod_src_mask[s_col_stod[:]]) # 0/1 array
    stod_col_lonocn = stod_src_lon[s_col_stod[:]]*stod_src_mask[s_col_stod[:]]
    stod_col_latocn = stod_src_lat[s_col_stod[:]]*stod_src_mask[s_col_stod[:]]
    stod_col_lonland = stod_src_lon[s_col_stod[:]]*(1-stod_src_mask[s_col_stod[:]])
    stod_col_latland = stod_src_lat[s_col_stod[:]]*(1-stod_src_mask[s_col_stod[:]])

# stod_row_ocn are points calculated in the destination (row/dst)
    stod_row_omask = np.asarray(stod_dst_mask[s_row_stod[:]]) # 0/1 array
    stod_row_lonocn = stod_dst_lon[s_row_stod[:]]*stod_dst_mask[s_row_stod[:]]
    stod_row_latocn = stod_dst_lat[s_row_stod[:]]*stod_dst_mask[s_row_stod[:]]
    stod_row_lonland = stod_dst_lon[s_row_stod[:]]*(1-stod_dst_mask[s_row_stod[:]])
    stod_row_latland = stod_dst_lat[s_row_stod[:]]*(1-stod_dst_mask[s_row_stod[:]])

# compare stod_dst_ocnmask and stod_row_ocnmask (row/dst)
    stod_dst_uncalc = np.setdiff1d(stod_dst_allmask,s_row_stod)
#    stod_dst_lonuncalc = stod_dst_lon[stod_dst_uncalc[:]]*stod_dst_unmask[stod_dst_uncalc[:]]
#    stod_dst_latuncalc = stod_dst_lat[stod_dst_uncalc[:]]*stod_dst_unmask[stod_dst_uncalc[:]]
#    stod_dst_lonuncalcland = stod_dst_lon[stod_dst_uncalc[:]]*(1-stod_dst_unmask[stod_dst_uncalc[:]])
#    stod_dst_latuncalcland = stod_dst_lat[stod_dst_uncalc[:]]*(1-stod_dst_unmask[stod_dst_uncalc[:]])
    stod_dst_lonuncalc = stod_dst_lon[stod_dst_uncalc[:]]*stod_dst_mask[stod_dst_uncalc[:]]
    stod_dst_latuncalc = stod_dst_lat[stod_dst_uncalc[:]]*stod_dst_mask[stod_dst_uncalc[:]]
    stod_dst_lonuncalcland = stod_dst_lon[stod_dst_uncalc[:]]*(1-stod_dst_mask[stod_dst_uncalc[:]])
    stod_dst_latuncalcland = stod_dst_lat[stod_dst_uncalc[:]]*(1-stod_dst_mask[stod_dst_uncalc[:]])



# stod_src_unused is the points not used in stod_row_ocn (not mapped) (col/src)
    stod_src_unused = np.setdiff1d(stod_src_allmask,s_col_stod)
    stod_src_lonunused = stod_src_lon[stod_src_unused[:]]*stod_src_mask[stod_src_unused[:]]
    stod_src_latunused = stod_src_lat[stod_src_unused[:]]*stod_src_mask[stod_src_unused[:]]
    stod_src_lonunusedland = stod_src_lon[stod_src_unused[:]]*(1-stod_src_mask[stod_src_unused[:]])
    stod_src_latunusedland = stod_src_lat[stod_src_unused[:]]*(1-stod_src_mask[stod_src_unused[:]])


    stod_src_x, stod_src_y = stodmap(stod_src_lon,stod_src_lat)
    stod_dst_x, stod_dst_y = stodmap(stod_dst_lon,stod_dst_lat)
    stod_src_xocn, stod_src_yocn = stodmap(stod_src_lonocn,stod_src_latocn)
    stod_dst_xocn, stod_dst_yocn = stodmap(stod_dst_lonocn,stod_dst_latocn)
    stod_src_xland, stod_src_yland = stodmap(stod_src_lonland,stod_src_latland)
    stod_dst_xland, stod_dst_yland = stodmap(stod_dst_lonland,stod_dst_latland)
    stod_col_x, stod_col_y = stodmap(stod_col_lon,stod_col_lat)
    stod_row_x, stod_row_y = stodmap(stod_row_lon,stod_row_lat)
    stod_col_xocn, stod_col_yocn = stodmap(stod_col_lonocn,stod_col_latocn)
    stod_row_xocn, stod_row_yocn = stodmap(stod_row_lonocn,stod_row_latocn)
    stod_col_xland,stod_col_yland = stodmap(stod_col_lonland,stod_col_latland)
    stod_row_xland,stod_row_yland = stodmap(stod_row_lonland,stod_row_latland)
    stod_src_xunused,stod_src_yunused = stodmap(stod_src_lonunused,stod_src_latunused)
    stod_dst_xuncalc,stod_dst_yuncalc = stodmap(stod_dst_lonuncalc,stod_dst_latuncalc)
    stod_src_xunusedland,stod_src_yunusedland = stodmap(stod_src_lonunusedland,stod_src_latunusedland)
    stod_dst_xuncalcland,stod_dst_yuncalcland = stodmap(stod_dst_lonuncalcland,stod_dst_latuncalcland)


    print 'here'

# here want to mark the points that are being used to create
#    stod_col_field = list()
#    stod_row_field = list()
#    stod_col_field[:] = stod_src_landmask[s_col_stod[:]]
#stod_row_field[:] = stod_dst_landmask[s_row_stod[:]]



# and a subplot for src
    srcmap = plt.subplot2grid((2,3),(0,0))

    stodmap.drawcoastlines()

# scatter markers
    stodmap.scatter(stod_src_xland,stod_src_yland,5,color='red',marker='+')
    stodmap.scatter(stod_src_xocn,stod_src_yocn,3,color='blue')

#    plt.title(cesm_remap_file[i])
    plt.title('SRC GRID')

# subplot for points used in interp
    srcfig = plt.subplot2grid((2,3),(0,1))
    stodmap.scatter(stod_col_xland,stod_col_yland,color='red',marker='+')
    stodmap.scatter(stod_col_xocn,stod_col_yocn,3,color='blue')

#   stodmap.drawcoastlines()
    plt.title('SRC IN')

# subplot for ocean points not used in interp
    unusedfig = plt.subplot2grid((2,3),(0,2))
    stodmap.scatter(stod_src_xunusedland,stod_src_yunusedland,5,color='red',marker='+')
    stodmap.scatter(stod_src_xunused,stod_src_yunused,3,color='blue')

    plt.title('SRC UNUSED')

# subplot for ocean points not calculated

    uncalcfig = plt.subplot2grid((2,3),(1,2))
    stodmap.scatter(stod_dst_xuncalcland,stod_dst_yuncalcland,5,color='red',marker='+')
    stodmap.scatter(stod_dst_xuncalc,stod_dst_yuncalc,3,color='blue')

    plt.title('DST UNCALC')

# subplot for points calculated in interp
    dstfig=plt.subplot2grid((2,3),(1,1))
#    stodmap.drawcoastlines()

    stodmap.scatter(stod_row_xland,stod_row_yland,color='red',marker='+')
    stodmap.scatter(stod_row_xocn,stod_row_yocn,3,color='blue')

#    for i in range(0,s_n):
#        if stod_row_field[i] == 1:
#            stodmap.scatter(stod_row_x[i],stod_row_y[i],5,color='red',marker='+')
#        else:
#            stodmap.scatter(stod_row_x[i],stod_row_y[i],3,color='blue')

    plt.title('DST OUT')


# and a subplot for dst
    dstmap = plt.subplot2grid((2,3),(1,0))

    stodmap.drawcoastlines()

# scatter markers
    stodmap.scatter(stod_dst_xland,stod_dst_yland,5,color='red',marker='+')
    stodmap.scatter(stod_dst_xocn,stod_dst_yocn,3,color='blue')

    plt.title('DST GRID')

    plt.suptitle(cesm_remap_stod_file[k])

    if showplots == 1:
        plt.show()

    plt.savefig(cesm_remap_stod_file[k]+'.png')

    if showplots == 0:
        plt.close()


### End stod Map block

# Now compare arrays s_row_stod and s_row_blin
#  to get the values of the SRC points in DST not mapped in ROW

# this checks for elements of s_row_stod not in s_row_blin
    blin_not_mapped = np.setdiff1d(s_row_stod,s_row_blin)
# this returns points that are unique to one or the other grid (so not what we want)
#    blin_not_mapped2 = np.setxor1d(s_row_stod,s_row_blin)

# copy to a new array
    s_row_xtra = np.copy(blin_not_mapped)

# for S to D, there will be exactly one element in COL corresponding to ROW
# so get the indices of s_row_stod that match the values of s_row_xtra

    xtra_row_list = []
    for ii in s_row_xtra:
        for j in range(len(s_row_stod)):
            if (ii == s_row_stod[j]):
                print 'Appending: ', ii, j, len(s_row_xtra), len(s_row_stod)
                xtra_row_list.append(j)

    xtra_index = np.asarray(xtra_row_list)

    s_col_xtra = s_col_stod[xtra_index]

    S_weight_xtra = S_weight_stod[xtra_index]


# now get the corresponding col points

    # separate into ocn and land points
    xtra_dst_ocnpts = s_row_xtra[np.nonzero(blin_dst_omask[s_row_xtra])]

    xtra_dst_landpts = s_row_xtra[np.nonzero(blin_dst_lmask[s_row_xtra])]

    xtra_dst_lonocn = blin_dst_lon[xtra_dst_ocnpts]
    xtra_dst_latocn = blin_dst_lat[xtra_dst_ocnpts]

    xtra_dst_lonland = blin_dst_lon[xtra_dst_landpts]
    xtra_dst_latland = blin_dst_lat[xtra_dst_landpts]

    xtra_src_ocnpts = s_col_xtra[np.nonzero(blin_src_omask[s_col_xtra])]
    xtra_src_landpts = s_col_xtra[np.nonzero(blin_src_lmask[s_col_xtra])]

    xtra_src_lonocn = blin_src_lon[xtra_src_ocnpts]
    xtra_src_latocn = blin_src_lat[xtra_src_ocnpts]

    xtra_src_lonland = blin_src_lon[xtra_src_landpts]
    xtra_src_latland = blin_src_lat[xtra_src_landpts]

# size the plot
    xtrafig = plt.figure(figsize=(12,7))

### Create a 2x2 plot grid
    xtrafig = gridspec.GridSpec(2,3)

### Create a map
    xtramap = Basemap(llcrnrlon=0,llcrnrlat=-90,urcrnrlon=360,urcrnrlat=90,projection='mill')

    xtra_dst_xocn, xtra_dst_yocn = xtramap(xtra_dst_lonocn,xtra_dst_latocn)
    xtra_dst_xland, xtra_dst_yland = xtramap(xtra_dst_lonland,xtra_dst_latland)

    xtra_src_xocn, xtra_src_yocn = xtramap(xtra_src_lonocn,xtra_src_latocn)
    xtra_src_xland, xtra_src_yland = xtramap(xtra_src_lonland,xtra_src_latland)

    # subplot for points used in interp
    srcfig = plt.subplot2grid((2,3),(0,2))
#    xtramap.scatter(xtra_not_x,xtra_not_y,3,color='purple')
    xtramap.scatter(xtra_dst_xland,xtra_dst_yland,color='red',marker='+')
    xtramap.scatter(xtra_dst_xocn,xtra_dst_yocn,3,color='blue')

    xtramap.drawcoastlines()
    plt.title('xtra DST OUT')

# subplot for points calculated in interp
    dstfig=plt.subplot2grid((2,3),(0,1))
#    stodmap.drawcoastlines()

    stodmap.scatter(stod_row_xland,stod_row_yland,color='red',marker='+')
    stodmap.scatter(stod_row_xocn,stod_row_yocn,3,color='blue')

    plt.title('STOD DST OUT')

    # subplot for points calculated in interp
    dstfig=plt.subplot2grid((2,3),(0,0))
#    blinmap.drawcoastlines()

    blinmap.scatter(blin_row_xland,blin_row_yland,color='red',marker='+')
    blinmap.scatter(blin_row_xocn,blin_row_yocn,3,color='blue')

    plt.title('BLIN DST OUT')

        # subplot for points used in interp
    srcfig = plt.subplot2grid((2,3),(1,2))
#    xtramap.scatter(xtra_not_x,xtra_not_y,3,color='purple')
    xtramap.scatter(xtra_src_xland,xtra_src_yland,color='red',marker='+')
    xtramap.scatter(xtra_src_xocn,xtra_src_yocn,3,color='blue')

    xtramap.drawcoastlines()
    plt.title('xtra SRC IN')

# subplot for points used in interp
    srcfig = plt.subplot2grid((2,3),(1,1))
    stodmap.scatter(stod_col_xland,stod_col_yland,color='red',marker='+')
    stodmap.scatter(stod_col_xocn,stod_col_yocn,3,color='blue')

#   stodmap.drawcoastlines()
    plt.title('STOD SRC IN')

    srcfig = plt.subplot2grid((2,3),(1,0))
    blinmap.scatter(blin_col_xland,blin_col_yland,color='red',marker='+')
    blinmap.scatter(blin_col_xocn,blin_col_yocn,3,color='blue')

#   blinmap.drawcoastlines()
    plt.title('BLIN SRC IN')


    if showplots == 1:
        plt.show()

    plt.savefig(cesm_remap_blin_file[k]+'.png')

    if showplots == 0:
        plt.close()


#####

# Now put it all together

# To copy the variables of the netCDF file
    S_weight_splice = np.append(S_weight_blin,S_weight_xtra)
    s_row_splice = np.append(s_row_blin,s_row_xtra)
    s_col_splice = np.append(s_col_blin,s_col_xtra)
    s_n_splice = len(S_weight_splice)

# delete the splice file if it exists
    try:
        os.remove(cesm_remap_splice_file[k])
    except OSError:
        pass

# copy the blin file to the splice file (overwrites existing file)
#    shutil.copyfile(cesm_remap_blin_file[k],cesm_remap_splice_file[k])

# To copy the global attributes of the netCDF file
# get the infrom from grid map
    blingroup = netCDF4.Dataset(cesm_remap_blin_file[k],'r')
    print blingroup.data_model

    # create the new spicegroup model
    splicegroup = netCDF4.Dataset(cesm_remap_splice_file[k],'w')
    print splicegroup.data_model

    print 'Attributes'
    for attname in blingroup.ncattrs():
        print attname
        setattr(splicegroup,attname,getattr(blingroup,attname))


# To copy the dimension of the netCDF file
    print 'Dimensions'
    for dimname,dim in blingroup.dimensions.iteritems():

# if you want to make changes in the dimensions of the new file
# you should add your own conditions here before the creation of the dimension.
#   print dimname, dim
        if (dimname == 'n_s'):
            print 'change n_s: old',dim, dimname, len(dim)
            dim = splicegroup.createDimension(dimname,size=s_n_splice)
            print 'change n_s: new',dim, dimname, len(dim)
        else:
            splicegroup.createDimension(dimname,len(dim))



    print 'Variables'
    for varname,ncvar in blingroup.variables.iteritems():
        print varname, ncvar
        print 'Type:',ncvar, type(ncvar)
# if you want to make changes in the variables of the new file
# you should add your own conditions here before the creation of the variable.
        if (varname == 'row'):
            print 'Change row here', varname, ncvar
            var = splicegroup.createVariable(varname,ncvar.dtype,dimensions=([u'n_s',]),fill_value=False)
#Proceed to copy the variable attributes
            for attname in ncvar.ncattrs():
                setattr(var,attname,getattr(ncvar,attname))
                print 'Attribute: ', attname
            var[:] = s_row_splice[:]+1
        elif (varname == 'col'):
            var = splicegroup.createVariable(varname,ncvar.dtype,dimensions=([u'n_s',]),fill_value=False)
            print 'Change col here', ncvar
#Proceed to copy the variable attributes
            for attname in ncvar.ncattrs():
                setattr(var,attname,getattr(ncvar,attname))
                print 'Attribute: ', attname
            var[:] = s_col_splice[:]+1
        elif (varname == 'S'):
            print 'Change S here',ncvar
#Proceed to copy the variable attributes
            for attname in ncvar.ncattrs():
                setattr(var,attname,getattr(ncvar,attname))
                print 'Attribute: ', attname
            var = splicegroup.createVariable(varname,ncvar.dtype,dimensions=([u'n_s',]),fill_value=False)
            var[:] = S_weight_splice
        else:
            var = splicegroup.createVariable(varname,ncvar.dtype,ncvar.dimensions,fill_value=False)
#Proceed to copy the variable attributes
            for attname in ncvar.ncattrs():
                setattr(var,attname,getattr(ncvar,attname))
                print 'Attribute: ', attname
#Finally copy the variable data to the new created variable
            print 'filling var with...', var, ncvar
            var[:] = ncvar[:]

    blingroup.close()
    splicegroup.close()


    # create the new spicegroup model
    splicegroup = netCDF4.Dataset(cesm_remap_splice_file[k],'r')
    print splicegroup.data_model

# Source variables for splice ... unchanged
    splice_src_griddims = splicegroup.variables['src_grid_dims'][:]
    nx_src = splice_src_griddims[0]
    ny_src = splice_src_griddims[1]
    splice_src_mask = splicegroup.variables['mask_a'][:]
    splice_src_lon = splicegroup.variables['xc_a'][:]
    splice_src_lat = splicegroup.variables['yc_a'][:]

# DST varaiables: splice (should be unchanged)

    splice_dst_griddims = splicegroup.variables['dst_grid_dims'][:]
    nx_dst = splice_dst_griddims[0]
    ny_dst = splice_dst_griddims[1]
    splice_dst_mask = splicegroup.variables['mask_b'][:]
    splice_dst_lon  = splicegroup.variables['xc_b'][:]
    splice_dst_lat = splicegroup.variables['yc_b'][:]

# the splice weighting function ... to be created from blin + xtra
    s_n_splice = len(splicegroup.dimensions['n_s'])
    src_n_splice = len(splicegroup.dimensions['n_a'])
    dst_n_splice = len(splicegroup.dimensions['n_b'])

# subtract 1 from row, col so index goes from 0 to nx*ny-1
    s_row_splice = splicegroup.variables['row'][:] -1 # row is the 1D index of the dst array
    s_col_splice = splicegroup.variables['col'][:] -1 # col is the 1D index of the src array
    S_weight_splice = splicegroup.variables['S'][:]

#    splicegroup.sync()

# close the grid file to avoid over-writing
    splicegroup.close()


    if ('map_gx' in cesm_remap_splice_file[k]):  # gx3v7 in radians
        print 'gx in radians'
        radlat = 180.*splice_src_lat[:]/math.pi
        splice_src_lat = radlat
        radlon = 180*splice_src_lon[:]/math.pi
        splice_src_lon = radlon

    if ('TO_gx' in cesm_remap_splice_file[k]):  # gx3v7 in radians
        print 'gx in radians'
        radlat = 180.*splice_dst_lat[:]/math.pi
        splice_dst_lat = radlat
        radlon = 180*splice_dst_lon[:]/math.pi
        splice_dst_lon = radlon


    print cesm_remap_stod_file[k]
    print 'stod_src_grid_dims', splice_src_griddims, src_n_splice
    print 'stod_dst_grid_dims', splice_dst_griddims, dst_n_splice
    print '#points used:', s_n_splice

# col is the source grid
    splice_col_lon = splice_src_lon[s_col_splice[:]]
    splice_col_lat = splice_src_lat[s_col_splice[:]]
# row is the destination grid
    splice_row_lon = splice_dst_lon[s_row_splice[:]]
    splice_row_lat = splice_dst_lat[s_row_splice[:]]

#end load stod block

### Spice mapping block
# size the plot
    splicefig = plt.figure(figsize=(12,7))

### Create a 2x2 plot grid
    splicefig = gridspec.GridSpec(2,3)

### Create a map
    splicemap = Basemap(llcrnrlon=0,llcrnrlat=-90,urcrnrlon=360,urcrnrlat=90,projection='mill')


# using splice_src_mask as a mask, create a masked array where 0 is water (not 1)
    splice_src_omask = np.asarray(splice_src_mask)   ### 0 = land, 1 = ocn
    splice_src_ocnmask = np.flatnonzero(splice_src_omask) ### indices where ocn
    splice_src_lmask = np.asarray(1 - splice_src_mask) ### 0 = ocn, 1 = land
    splice_src_landmask = np.flatnonzero(splice_src_lmask) ### indices where land
    splice_src_allmask = np.arange(len(splice_src_mask)) # an array with all indices

    splice_src_lonocn = splice_src_lon[splice_src_ocnmask]
    splice_src_latocn = splice_src_lat[splice_src_ocnmask]
    splice_src_ocn = splice_src_omask[splice_src_ocnmask] # shouldbe an array of 1s
    splice_src_lonland = splice_src_lon[splice_src_landmask]
    splice_src_latland = splice_src_lat[splice_src_landmask]
    splice_src_land = splice_src_lmask[splice_src_landmask] # should be an array of 1s

# also do for splice_dst_mask
    splice_dst_omask = np.asarray(splice_dst_mask)
    splice_dst_ocnmask = np.flatnonzero(splice_dst_omask)
    splice_dst_lmask = np.asarray(1-splice_dst_mask)
    splice_dst_landmask = np.flatnonzero(splice_dst_lmask)
    splice_dst_allmask = np.arange(len(splice_dst_mask))

    splice_dst_lonocn = splice_dst_lon[splice_dst_ocnmask]
    splice_dst_latocn = splice_dst_lat[splice_dst_ocnmask]
    splice_dst_ocnpts = splice_dst_omask[splice_dst_ocnmask] # should be an array of 1s
    splice_dst_lonland = splice_dst_lon[splice_dst_landmask]
    splice_dst_latland = splice_dst_lat[splice_dst_landmask]
    splice_dst_landpts = splice_dst_lmask[splice_dst_landmask] #should be an array of 1s

#### now find which points are used in the mapping

#  s_col_splice specifies which src points are calcuated (col/src)
    splice_col_omask = np.asarray(splice_src_mask[s_col_splice[:]]) # 0/1 array
    splice_col_lonocn = splice_src_lon[s_col_splice[:]]*splice_src_mask[s_col_splice[:]]
    splice_col_latocn = splice_src_lat[s_col_splice[:]]*splice_src_mask[s_col_splice[:]]
    splice_col_lonland = splice_src_lon[s_col_splice[:]]*(1-splice_src_mask[s_col_splice[:]])
    splice_col_latland = splice_src_lat[s_col_splice[:]]*(1-splice_src_mask[s_col_splice[:]])

# splice_row_ocn are points calculated in the destination (row/dst)
    splice_row_omask = np.asarray(splice_dst_mask[s_row_splice[:]]) # 0/1 array
    splice_row_lonocn = splice_dst_lon[s_row_splice[:]]*splice_dst_mask[s_row_splice[:]]
    splice_row_latocn = splice_dst_lat[s_row_splice[:]]*splice_dst_mask[s_row_splice[:]]
    splice_row_lonland = splice_dst_lon[s_row_splice[:]]*(1-splice_dst_mask[s_row_splice[:]])
    splice_row_latland = splice_dst_lat[s_row_splice[:]]*(1-splice_dst_mask[s_row_splice[:]])

# compare splice_dst_ocnmask and splice_row_ocnmask (row/dst)
    splice_dst_uncalc = np.setdiff1d(splice_dst_allmask,s_row_splice)
#    splice_dst_lonuncalc = splice_dst_lon[splice_dst_uncalc[:]]*splice_dst_unmask[splice_dst_uncalc[:]]
#    splice_dst_latuncalc = splice_dst_lat[splice_dst_uncalc[:]]*splice_dst_unmask[splice_dst_uncalc[:]]
#    splice_dst_lonuncalcland = splice_dst_lon[splice_dst_uncalc[:]]*(1-splice_dst_unmask[splice_dst_uncalc[:]])
#    splice_dst_latuncalcland = splice_dst_lat[splice_dst_uncalc[:]]*(1-splice_dst_unmask[splice_dst_uncalc[:]])
    splice_dst_lonuncalc = splice_dst_lon[splice_dst_uncalc[:]]*splice_dst_mask[splice_dst_uncalc[:]]
    splice_dst_latuncalc = splice_dst_lat[splice_dst_uncalc[:]]*splice_dst_mask[splice_dst_uncalc[:]]
    splice_dst_lonuncalcland = splice_dst_lon[splice_dst_uncalc[:]]*(1-splice_dst_mask[splice_dst_uncalc[:]])
    splice_dst_latuncalcland = splice_dst_lat[splice_dst_uncalc[:]]*(1-splice_dst_mask[splice_dst_uncalc[:]])



# splice_src_unused is the points not used in splice_row_ocn (not mapped) (col/src)
    splice_src_unused = np.setdiff1d(splice_src_allmask,s_col_splice)
    splice_src_lonunused = splice_src_lon[splice_src_unused[:]]*splice_src_mask[splice_src_unused[:]]
    splice_src_latunused = splice_src_lat[splice_src_unused[:]]*splice_src_mask[splice_src_unused[:]]
    splice_src_lonunusedland = splice_src_lon[splice_src_unused[:]]*(1-splice_src_mask[splice_src_unused[:]])
    splice_src_latunusedland = splice_src_lat[splice_src_unused[:]]*(1-splice_src_mask[splice_src_unused[:]])


    splice_src_x, splice_src_y = splicemap(splice_src_lon,splice_src_lat)
    splice_dst_x, splice_dst_y = splicemap(splice_dst_lon,splice_dst_lat)
    splice_src_xocn, splice_src_yocn = splicemap(splice_src_lonocn,splice_src_latocn)
    splice_dst_xocn, splice_dst_yocn = splicemap(splice_dst_lonocn,splice_dst_latocn)
    splice_src_xland, splice_src_yland = splicemap(splice_src_lonland,splice_src_latland)
    splice_dst_xland, splice_dst_yland = splicemap(splice_dst_lonland,splice_dst_latland)
    splice_col_x, splice_col_y = splicemap(splice_col_lon,splice_col_lat)
    splice_row_x, splice_row_y = splicemap(splice_row_lon,splice_row_lat)
    splice_col_xocn, splice_col_yocn = splicemap(splice_col_lonocn,splice_col_latocn)
    splice_row_xocn, splice_row_yocn = splicemap(splice_row_lonocn,splice_row_latocn)
    splice_col_xland,splice_col_yland = splicemap(splice_col_lonland,splice_col_latland)
    splice_row_xland,splice_row_yland = splicemap(splice_row_lonland,splice_row_latland)
    splice_src_xunused,splice_src_yunused = splicemap(splice_src_lonunused,splice_src_latunused)
    splice_dst_xuncalc,splice_dst_yuncalc = splicemap(splice_dst_lonuncalc,splice_dst_latuncalc)
    splice_src_xunusedland,splice_src_yunusedland = splicemap(splice_src_lonunusedland,splice_src_latunusedland)
    splice_dst_xuncalcland,splice_dst_yuncalcland = splicemap(splice_dst_lonuncalcland,splice_dst_latuncalcland)


    print 'splicing done'

# here want to mark the points that are being used to create
#    splice_col_field = list()
#    splice_row_field = list()
#    splice_col_field[:] = splice_src_landmask[s_col_splice[:]]
#splice_row_field[:] = splice_dst_landmask[s_row_splice[:]]



# and a subplot for src
    srcmap = plt.subplot2grid((2,3),(0,0))

    splicemap.drawcoastlines()

# scatter markers
    splicemap.scatter(splice_src_xland,splice_src_yland,5,color='red',marker='+')
    splicemap.scatter(splice_src_xocn,splice_src_yocn,3,color='blue')

#    plt.title(cesm_remap_file[i])
    plt.title('SRC GRID')

# subplot for points used in interp
    srcfig = plt.subplot2grid((2,3),(0,1))
    splicemap.scatter(splice_col_xland,splice_col_yland,color='red',marker='+')
    splicemap.scatter(splice_col_xocn,splice_col_yocn,3,color='blue')

#   splicemap.drawcoastlines()
    plt.title('SRC IN')

# subplot for ocean points not used in interp
    unusedfig = plt.subplot2grid((2,3),(0,2))
    splicemap.scatter(splice_src_xunusedland,splice_src_yunusedland,5,color='red',marker='+')
    splicemap.scatter(splice_src_xunused,splice_src_yunused,3,color='blue')

    plt.title('SRC UNUSED')

# subplot for ocean points not calculated

    uncalcfig = plt.subplot2grid((2,3),(1,2))
    splicemap.scatter(splice_dst_xuncalcland,splice_dst_yuncalcland,5,color='red',marker='+')
    splicemap.scatter(splice_dst_xuncalc,splice_dst_yuncalc,3,color='blue')

    plt.title('DST UNCALC')

# subplot for points calculated in interp
    dstfig=plt.subplot2grid((2,3),(1,1))
#    splicemap.drawcoastlines()

    splicemap.scatter(splice_row_xland,splice_row_yland,color='red',marker='+')
    splicemap.scatter(splice_row_xocn,splice_row_yocn,3,color='blue')

#    for i in range(0,s_n):
#        if splice_row_field[i] == 1:
#            splicemap.scatter(splice_row_x[i],splice_row_y[i],5,color='red',marker='+')
#        else:
#            splicemap.scatter(splice_row_x[i],splice_row_y[i],3,color='blue')

    plt.title('DST OUT')


# and a subplot for dst
    dstmap = plt.subplot2grid((2,3),(1,0))

    splicemap.drawcoastlines()

# scatter markers
    splicemap.scatter(splice_dst_xland,splice_dst_yland,5,color='red',marker='+')
    splicemap.scatter(splice_dst_xocn,splice_dst_yocn,3,color='blue')

    plt.title('DST GRID')

    plt.suptitle(cesm_remap_splice_file[k])

    plt.show()

    if showplots == 1:
        plt.show()
    plt.savefig(cesm_remap_splice_file[k]+'.png')

    if showplots == 0:
       plt.close()


### End splice Map block



    if (cleanup):

        del blin_col_x
        del blin_col_y
        del blin_row_x
        del blin_row_y
        del blin_src_x
        del blin_src_y
        del blin_dst_x
        del blin_dst_y
        del blin_dst_ocnpts
        del blin_dst_landpts
        del blin_src_lon
        del blin_src_lat
        del blin_dst_lon
        del blin_dst_lat
        del blin_col_lon
        del blin_col_lat
        del blin_row_lon
        del blin_row_lat

        del blinmap

        del stod_col_x
        del stod_col_y
        del stod_row_x
        del stod_row_y
        del stod_src_x
        del stod_src_y
        del stod_dst_x
        del stod_dst_y
        del stod_dst_ocnpts
        del stod_dst_landpts
        del stod_src_lon
        del stod_src_lat
        del stod_dst_lon
        del stod_dst_lat
        del stod_col_lon
        del stod_col_lat
        del stod_row_lon
        del stod_row_lat

        del stodmap

        del splice_col_x
        del splice_col_y
        del splice_row_x
        del splice_row_y
        del splice_src_x
        del splice_src_y
        del splice_dst_x
        del splice_dst_y
        del splice_dst_ocnpts
        del splice_dst_landpts
        del splice_src_lon
        del splice_src_lat
        del splice_dst_lon
        del splice_dst_lat
        del splice_col_lon
        del splice_col_lat
        del splice_row_lon
        del splice_row_lat

        del splicemap

        del s_row_xtra
        del xtra_row_list
        del xtra_index
        del s_col_xtra
        del S_weight_xtra
