#!/usr/bin/env python3
import os, sys, csv, time, math
from optparse import OptionParser
import numpy
import netcdf4_functions as nffun

parser = OptionParser()

parser.add_option("--site", dest="site", default='', \
                  help = '6-character FLUXNET code to run (required)')
parser.add_option("--sitegroup", dest="sitegroup", default="AmeriFlux", \
                  help = "site group to use (default AmeriFlux)")
parser.add_option("--lat_bounds", dest="lat_bounds", default='-999,-999', \
                  help = 'latitude range for regional run')
parser.add_option("--lon_bounds", dest="lon_bounds", default='-999,-999', \
                  help = 'longitude range for regional run')
parser.add_option("--lai", dest="lai", default=-999, \
                  help = 'Set constant LAI (SP mode only)')
parser.add_option("--model", dest="mymodel", default='', \
                 help = 'Model to use (CLM5 or ELM)')
parser.add_option("--mask", dest="mymask", default='', \
                  help = 'Mask file to use (regional only)')
parser.add_option("--pft", dest="mypft", default=-1, \
                  help = 'Replace all gridcell PFT with this value')
parser.add_option("--point_list", dest="point_list", default='', \
                  help = 'File containing list of points to run (unstructured)')
parser.add_option("--keep_duplicates", dest="keep_duplicates", default=False, \
                  help = 'Keep duplicate points', action='store_true')
parser.add_option("--ccsm_input", dest="ccsm_input", \
                  default='../../../../ccsm_inputdata', \
                  help = "input data directory for CESM (required)")
parser.add_option("--metdir", dest="metdir", default="none", \
                  help = 'subdirectory for met data forcing')
parser.add_option("--makemetdata", dest="makemet", default=False, \
		  help = 'Generate meteorology', action="store_true")
parser.add_option("--surfdata_grid", dest="surfdata_grid", default=False, \
                  help = 'Use gridded soil data', action="store_true")
#parser.add_option("--include_nonveg", dest="include_nonveg", default=False, \
#                      help = "Include non-vegetated fractions from surface data file")
parser.add_option("--res", dest="res", default="hcru_hcru", \
                     help = 'Resolution of global files')
parser.add_option("--nopftdyn", dest="nopftdyn", default=False, \
                     action='store_true', help='Do not make transient PFT file')
parser.add_option("--mysimyr", dest="mysimyr", default=1850, \
                     help = 'Simulation year (1850 or 2000)')
(options, args) = parser.parse_args()


ccsm_input = os.path.abspath(options.ccsm_input)

#------------------- get site information ----------------------------------

#Remove existing temp files
os.system('rm temp/*.nc')

lat_bounds = options.lat_bounds.split(',')
lon_bounds = options.lon_bounds.split(',')
lat_bounds = [float(l) for l in lat_bounds]
lon_bounds = [float(l) for l in lon_bounds]

mysimyr=int(options.mysimyr)

if ('hcru' in options.res):
    resx = 0.5
    resy = 0.5
    domainfile_orig = ccsm_input+'/share/domains/domain.clm/domain.lnd.360x720_cruncep.100429.nc'

    if (options.mymodel == 'CLM5'):
        surffile_orig = ccsm_input+'/lnd/clm2/surfdata_map/surfdata_360x720cru_16pfts_Irrig_CMIP6_simyr1850_c170824.nc'
    else:
        if (mysimyr == 2000):
            surffile_orig =  ccsm_input+'/lnd/clm2/surfdata_map/surfdata_360x720cru_simyr2000_c160307.nc'
        else:
            #CMIP6 stype (Hurtt v2)
            surffile_orig = ccsm_input+'/lnd/clm2/surfdata_map/surfdata_360x720cru_simyr1850_c180216.nc'

    pftdyn_orig = ccsm_input+'/lnd/clm2/surfdata_map/landuse.timeseries_360x720cru_hist_simyr1850-2015_c180220.nc'
    nyears_landuse=166
elif ('f19' in options.res):
    domainfile_orig = ccsm_input+'/share/domains/domain.lnd.fv1.9x2.5_gx1v6.090206.nc'
    surffile_orig =  ccsm_input+'/lnd/clm2/surfdata_map/surfdata_1.9x2.5_simyr1850_c171002.nc'
    pftdyn_orig = ccsm_input+'/lnd/clm2/surfdata_map/landuse.timeseries_1.9x2.5_rcp8.5_simyr1850-2100_c141219.nc'
    nyears_landuse=251
    resx = 2.5
    resy = 1.9
elif ('f09' in options.res):
    domainfile_orig = ccsm_input+'/share/domains/domain.lnd.fv0.9x1.25_gx1v6.090309.nc'
    surffile_orig =  ccsm_input+'/lnd/clm2/surfdata_map/surfdata_0.9x1.25_simyr1850_c171002.nc'
    pftdyn_orig = ccsm_input+'/lnd/clm2/surfdata_map/landuse.timeseries_0.9x1.25_rcp8.5_simyr1850-2100_c141219.nc'
    nyears_landuse=251
    resx = 1.25
    resy = 0.9
elif ('ne30' in options.res):
    #domainfile_orig = ccsm_input+'/share/domains/domain.lnd.ne30np4_oEC60to30.20151214.nc'
    #water cycle experiment
    domainfile_orig = ccsm_input+'/share/domains/domain.lnd.ne30np4_oEC60to30v3.161222.nc'  
    #surffile_orig   = ccsm_input+'/lnd/clm2/surfdata_map/surfdata_ne30np4_simyr1850_2015_c171018.nc'
    surffile_orig   = ccsm_input+'/lnd/clm2/surfdata_map/surfdata_ne30np4_simyr1850_c180306.nc'
    pftdyn_orig     = ccsm_input+'/lnd/clm2/surfdata_map/landuse.timeseries_ne30np4_hist_simyr1850_2015_c20171018.nc'
    nyears_landuse  = 166

n_grids=1
issite = False
isglobal = False
lat=[]
lon=[]
if (lat_bounds[0] > -90 and lon_bounds[0] > -180):
    print( '\nCreating regional datasets using '+options.res+ 'resolution')
    if (lon_bounds[0] < 0):
        lon_bounds[0] = lon_bounds[0]+360.
    if (lon_bounds[1] < 0):
        lon_bounds[1] = lon_bounds[1]+360.
elif (options.point_list != ''):
    issite=True
    input_file = open(options.point_list,'r')
    n_grids=0
    point_pfts=[]
    for s in input_file:
	if (n_grids == 0):
	     header = s.split()
        else:
             data = s.split()
             dnum=0
             point_pfts.append(-1)
             for d in data:
                 if ('lon' in header[dnum]): 
                      mylon = float(d)
                      if (mylon < 0):
                          mylon = mylon+360
                      lon.append(mylon)
                 elif ('lat' in header[dnum]):
                      lat.append(float(d))
                 elif ('pft' in header[dnum]):
                      point_pfts[n_grids-1] = int(d)
                 if (int(options.mypft) >= 0):    #overrides info in file
                     point_pfts[n_grids-1] = options.mypft
                 dnum=dnum+1
        n_grids=n_grids+1
    input_file.close()
    n_grids = n_grids-1
elif (options.site != ''):
    print('\nCreating datasets for '+options.site+' using '+options.res+' resolution')
    issite = True
    AFdatareader = csv.reader(open(ccsm_input+'/lnd/clm2/PTCLM/'+options.sitegroup+'_sitedata.txt',"rb"))
    for row in AFdatareader:
        if row[0] == options.site:
            mylon=float(row[3])
            if (mylon < 0):
                mylon=360.0+float(row[3])
            lon.append(mylon)
            lat.append(float(row[4]))
            if ('US-SPR' in options.site):
                lon.append(mylon)
                lat.append(float(row[4]))
                n_grids = 2
            startyear=int(row[6])
            endyear=int(row[7])
            alignyear = int(row[8])
else:
    isglobal=True

#get corresponding 0.5x0.5 and 1.9x2.5 degree grid cells
if (options.res == 'hcru_hcru'):
     longxy = (numpy.cumsum(numpy.ones([721]))-1)*0.5
     latixy = (numpy.cumsum(numpy.ones([361]))-1)*0.5 -90.0
elif (options.res == 'f19_f19'):
    longxy = (numpy.cumsum(numpy.ones([145]))-1)*2.5-1.25
    latixy_centers = (numpy.cumsum(numpy.ones([96]))-1)*(180.0/95) - 90.0
    latixy = numpy.zeros([97], numpy.float)
    longxy[0]   = 0
    latixy[0]   =  -90
    latixy[96]  =  90
    for i in range(1,96):
        latixy[i] = (latixy_centers[i-1]+latixy_centers[i])/2.0
else:
    longxy = nffun.getvar(surffile_orig, 'LONGXY')
    latixy = nffun.getvar(surffile_orig, 'LATIXY')

xgrid_min=[]
xgrid_max=[]
ygrid_min=[]
ygrid_max=[]
for n in range(0,n_grids):
    if (issite):
        lon_bounds = [lon[n],lon[n]]
        lat_bounds = [lat[n],lat[n]]
    xgrid_min.append(-1)
    xgrid_max.append(-1)
    ygrid_min.append(-1)
    ygrid_max.append(-1)
    if ('ne' in options.res):
      if (lon_bounds[0] != lon_bounds[1] or lat_bounds[0] != lat_bounds[1]):
        print('Regional subsets not allowed for ne resolutions.  Use point lists instead')
        sys.exit()
      ygrid_min[n] = 0
      ygrid_max[n] = 0
      mindist=99999
      for i in range(0,longxy.shape[0]-1):
        thisdist = (lon_bounds[0]-longxy[i])**2 + (lat_bounds[0]-latixy[i])**2
        if (thisdist < mindist):
          xgrid_min[n] = i
          xgrid_max[n] = i
          mindist=thisdist
    else:
      for i in range(0,longxy.shape[0]-1):
        if (lon_bounds[0] >= longxy[i]):
            xgrid_min[n] = i
            xgrid_max[n] = i
        elif (lon_bounds[1] >= longxy[i+1]):
            xgrid_max[n] = i
      if (lon_bounds[0] == 180 and lon_bounds[1] == 180):  #global
        xgrid_min[n] = 0
        xgrid_max[n] = longxy.shape[0]-2
      for i in range(0,latixy.shape[0]-1):
        if (lat_bounds[0] >= latixy[i]):
            ygrid_min[n] = i
            ygrid_max[n] = i
        elif (lat_bounds[1] >= latixy[i+1]):
            ygrid_max[n] = i
    #print n, lat[n], lon[n], xgrid_max[n], ygrid_max[n]
if (n_grids > 1 and options.site == ''):       #remove duplicate points
  n_grids_uniq = 1
  n_dups = 0
  xgrid_min_uniq = [xgrid_min[0]]
  ygrid_min_uniq = [ygrid_min[0]]
  lon_uniq = [lon[0]]
  lat_uniq = [lat[0]]
  point_pfts_uniq = [point_pfts[0]]
  point_index = [1]
  myoutput = open('point_list_output.txt','w')
  myoutput.write(str(lon[0])+','+str(lat[0])+','+str(point_index[0])+'\n')
  for n in range (1,n_grids):
      is_unique = True
      for m in range(0,n_grids_uniq):
          if (xgrid_min[n] == xgrid_min_uniq[m] and ygrid_min[n] == ygrid_min_uniq[m] \
              and point_pfts[n] == point_pfts_uniq[m]):
               n_dups = n_dups+1
               is_unique = False
               #point_index.append(m+1)
      if (is_unique or options.keep_duplicates):
          xgrid_min_uniq.append(xgrid_min[n])
          ygrid_min_uniq.append(ygrid_min[n])
          point_pfts_uniq.append(point_pfts[n])      
          lon_uniq.append(lon[n])
          lat_uniq.append(lat[n])
          n_grids_uniq = n_grids_uniq+1
          point_index.append(n_grids_uniq)
      myoutput.write(str(lon[n])+','+str(lat[n])+','+str(point_index[n])+'\n')
  myoutput.close()
  xgrid_min = xgrid_min_uniq
  xgrid_max = xgrid_min_uniq
  ygrid_min = ygrid_min_uniq
  ygrid_max = ygrid_min_uniq
  lon = lon_uniq
  lat = lat_uniq
  point_pfts = point_pfts_uniq
  n_grids = n_grids_uniq
  print n_grids, ' Unique points'
  print n_dups, ' duplicate points removed'
  print len(point_index)
  print point_index
#---------------------Create domain data --------------------------------------------------

print('Creating domain data')
os.system('mkdir -p ./temp')

domainfile_list=''
for n in range(0,n_grids):
    nst = str(100000+n)[1:]
    domainfile_new = './temp/domain'+nst+'.nc'
    if (not os.path.exists(domainfile_orig)):
        print('Error:  '+domainfile_orig+' does not exist.  Aborting')
        sys.exit(1)

    if (isglobal):
        os.system('cp '+domainfile_orig+' '+domainfile_new)
    else:
        os.system('ncks -d ni,'+str(xgrid_min[n])+','+str(xgrid_max[n])+' -d nj,'+str(ygrid_min[n])+ \
              ','+str(ygrid_max[n])+' '+domainfile_orig+' '+domainfile_new)

    if (issite):
        frac = nffun.getvar(domainfile_new, 'frac')
        mask = nffun.getvar(domainfile_new, 'mask')
        xc = nffun.getvar(domainfile_new, 'xc')
        yc = nffun.getvar(domainfile_new, 'yc')
        xv = nffun.getvar(domainfile_new, 'xv')
        yv = nffun.getvar(domainfile_new, 'yv')
        area = nffun.getvar(domainfile_new, 'area')
        frac[0] = 1.0
        mask[0] = 1
        if (options.site != ''):
            xc[0] = lon[n]
            yc[0] = lat[n]
            xv[0][0][0] = lon[n]-resx/2
            xv[0][0][1] = lon[n]+resx/2
            xv[0][0][2] = lon[n]-resx/2
            xv[0][0][3] = lon[n]+resx/2
            yv[0][0][0] = lat[n]-resy/2
            yv[0][0][1] = lat[n]-resy/2
            yv[0][0][2] = lat[n]+resy/2
            yv[0][0][3] = lat[n]+resy/2
            area[0] = resx*resy*math.pi/180*math.pi/180
            ierr = nffun.putvar(domainfile_new, 'xc', xc)
            ierr = nffun.putvar(domainfile_new, 'yc', yc)
            ierr = nffun.putvar(domainfile_new, 'xv', xv)
            ierr = nffun.putvar(domainfile_new, 'yv', yv)
            ierr = nffun.putvar(domainfile_new, 'area', area)
       
        ierr = nffun.putvar(domainfile_new, 'frac', frac)
        ierr = nffun.putvar(domainfile_new, 'mask', mask)
        os.system('ncks -O --mk_rec_dim nj '+domainfile_new+' '+domainfile_new)
    elif (options.mymask != ''):
       print 'Applying mask from '+options.mymask
       os.system('ncks -d lon,'+str(xgrid_min[n])+','+str(xgrid_max[n])+' -d lat,'+str(ygrid_min[n])+ \
              ','+str(ygrid_max[n])+' '+options.mymask+' mask_temp.nc')
       newmask = nffun.getvar('mask_temp.nc', 'PNW_mask')
       ierr = nffun.putvar(domainfile_new, 'mask', newmask)
       os.system('rm mask_temp.nc')

    domainfile_list = domainfile_list+' '+domainfile_new

domainfile_new = './temp/domain.nc'
if (n_grids > 1):
    os.system('ncrcat '+domainfile_list+' '+domainfile_new)
    os.system('nccopy -u  '+domainfile_new+' '+domainfile_new+'.tmp')
    os.system('ncpdq -O -a ni,nj '+domainfile_new+'.tmp '+domainfile_new)
    #os.system('ncwa -O -a ni -d ni,0,0 '+domainfile_new+'.tmp1 '+domainfile_new+'.tmp2')
    os.system('ncrename -h -O -d ni,ni_temp '+domainfile_new+' '+domainfile_new+' ')
    os.system('ncrename -h -O -d nj,ni '+domainfile_new+' '+domainfile_new+' ')
    os.system('ncrename -h -O -d ni_temp,nj '+domainfile_new+' '+domainfile_new+' ')
    os.system('rm ./temp/domain?????.nc*')
    #os.system('mv '+domainfile_new+'.tmp3 '+domainfile_new)
    #os.system('rm '+domainfile_new+'.tmp*')
else:
    os.system('mv '+domainfile_list+' '+domainfile_new)

#-------------------- create surface data ----------------------------------
print('Creating surface data')

surffile_list = ''
for n in range(0,n_grids):
    nst = str(100000+n)[1:]
    surffile_new =  './temp/surfdata'+nst+'.nc'
    if (not os.path.exists(surffile_orig)):
        print 'Error:  '+surffile_orig+' does not exist.  Aborting'
        sys.exit(1)
    if (isglobal):
        os.system('cp '+surffile_orig+' '+surffile_new)
    else:
        if ('ne' in options.res):
          os.system('ncks --fix_rec_dmn time -d gridcell,'+str(xgrid_min[n])+','+str(xgrid_max[n])+ \
            ' '+surffile_orig+' '+surffile_new)
        else:
          os.system('ncks --fix_rec_dmn time -d lsmlon,'+str(xgrid_min[n])+','+str(xgrid_max[n])+ \
             ' -d lsmlat,'+str(ygrid_min[n])+','+str(ygrid_max[n])+' '+surffile_orig+' '+surffile_new)

    if (issite):
        landfrac_pft = nffun.getvar(surffile_new, 'LANDFRAC_PFT')
        pftdata_mask = nffun.getvar(surffile_new, 'PFTDATA_MASK')
        longxy       = nffun.getvar(surffile_new, 'LONGXY')
        latixy       = nffun.getvar(surffile_new, 'LATIXY')
        area         = nffun.getvar(surffile_new, 'AREA')
        pct_wetland  = nffun.getvar(surffile_new, 'PCT_WETLAND')
        pct_lake     = nffun.getvar(surffile_new, 'PCT_LAKE')
        pct_glacier  = nffun.getvar(surffile_new, 'PCT_GLACIER')
        pct_urban    = nffun.getvar(surffile_new, 'PCT_URBAN')
        if (options.mymodel == 'CLM5'):
            pct_crop = nffun.getvar(surffile_new, 'PCT_CROP')
        else:
          soil_order   = nffun.getvar(surffile_new, 'SOIL_ORDER')
          labilep      = nffun.getvar(surffile_new, 'LABILE_P')
          primp        = nffun.getvar(surffile_new, 'APATITE_P')
          secondp      = nffun.getvar(surffile_new, 'SECONDARY_P')
          occlp        = nffun.getvar(surffile_new, 'OCCLUDED_P')

        #input from site-specific information
        soil_color   = nffun.getvar(surffile_new, 'SOIL_COLOR')
        pct_sand     = nffun.getvar(surffile_new, 'PCT_SAND')
        pct_clay     = nffun.getvar(surffile_new, 'PCT_CLAY')
        organic      = nffun.getvar(surffile_new, 'ORGANIC')
        fmax         = nffun.getvar(surffile_new, 'FMAX')
        pct_nat_veg  = nffun.getvar(surffile_new, 'PCT_NATVEG')
        pct_pft      = nffun.getvar(surffile_new, 'PCT_NAT_PFT') 
        monthly_lai  = nffun.getvar(surffile_new, 'MONTHLY_LAI')
        monthly_sai  = nffun.getvar(surffile_new, 'MONTHLY_SAI')
        monthly_height_top = nffun.getvar(surffile_new, 'MONTHLY_HEIGHT_TOP')
        monthly_height_bot = nffun.getvar(surffile_new, 'MONTHLY_HEIGHT_BOT')

        npft = 17
        if (options.mymodel == 'CLM5'):
            npft = 15

        #read file for site-specific PFT information
        mypft_frac=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
        mypct_sand = 0.0 
        mypct_clay = 0.0
 
        if (options.surfdata_grid == False and options.site != ''):
            AFdatareader = csv.reader(open(ccsm_input+'/lnd/clm2/PTCLM/'+options.sitegroup+'_pftdata.txt','rb'))
            for row in AFdatareader:
                if row[0] == options.site:
                    for thispft in range(0,5):
                        mypft_frac[int(row[2+2*thispft])]=float(row[1+2*thispft])
            if (sum(mypft_frac[0:npft]) == 0.0):
                print('*** Warning:  PFT data NOT found.  Using gridded data ***')
        #read file for site-specific soil information
            AFdatareader = csv.reader(open(ccsm_input+'/lnd/clm2/PTCLM/'+options.sitegroup+'_soildata.txt','rb'))
            for row in AFdatareader:
                if row[0] == options.site:
                    mypct_sand = row[4]
                    mypct_clay = row[5]
            if (mypct_sand == 0.0 and mypct_clay == 0.0):
                print('*** Warning:  Soil data NOT found.  Using gridded data ***')
        else:
          try:
            mypft_frac[point_pfts[n]] = 100.0
          except NameError:
            print('using PFT information from surface data')

        #landfrac_pft[0][0] = 1.0
        #pftdata_mask[0][0] = 1

        if (options.site != ''):
            longxy[0][0] = lon[n]
            latixy[0][0] = lat[n]
            area[0] = 111.2*resy*111.321*math.cos((lon[n]*resx)*math.pi/180)*resx

        if (not options.surfdata_grid or sum(mypft_frac[0:npft]) > 0.0):
            pct_wetland[0][0] = 0.0
            pct_lake[0][0]    = 0.0
            pct_glacier[0][0] = 0.0
            if (options.mymodel == 'CLM5'):
                pct_crop[0][0] = 0.0
            if ('US-SPR' in options.site and mysimyr !=2000):
                #SPRUCE P initial data
                soil_order[0][0] = 3
                labilep[0][0]    = 4.0
                primp[0][0]      = 1.0
                secondp[0][0]    = 10.0
                occlp[0][0]      = 5.0

            pct_nat_veg[0][0] = 100.0
            for k in range(0,3):
                pct_urban[k][0][0] = 0.0
            for k in range(0,10):
                if (mypct_sand > 0.0 or mypct_clay > 0.0):
                    if (k == 0):
                       print 'Setting %sand to '+str(mypct_sand)
                       print 'Setting %clay to '+str(mypct_clay)
                    pct_sand[k][0][0]   = mypct_sand
                    pct_clay[k][0][0]   = mypct_clay
                if ('US-SPR' in options.site):
                    if (k < 8):
                        organic[k][0][0] = 130.0
                    elif (k == 8):
                        organic[k][0][0] = 65.0
            pft_names=['Bare ground','ENF Temperate','ENF Boreal','DNF Boreal','EBF Tropical', \
                       'EBF Temperate', 'DBF Tropical', 'DBF Temperate', 'DBF Boreal', 'EB Shrub' \
                       , 'DB Shrub Temperate', 'BD Shrub Boreal', 'C3 arctic grass', \
                       'C3 non-arctic grass', 'C4 grass', 'Crop','xxx','xxx']
            if (options.mypft >= 0):
              print 'Setting PFT '+str(options.mypft)+'('+pft_names[int(options.mypft)]+') to 100%'
              pct_pft[:,0,0] = 0.0
              pct_pft[int(options.mypft),0,0] = 100.0
            else:
              for p in range(0,npft):
                if (sum(mypft_frac[0:npft]) > 0.0):
                    if (mypft_frac[p] > 0.0):
                        if (p < 16):
                           print 'Setting PFT '+str(p)+'('+pft_names[p]+') to '+ \
                           str(mypft_frac[p])+'%'
                        else:
                           print 'Setting PFT '+str(p)+' to '+str(mypft_frac[p])+'%'
                    pct_pft[p][0][0] = mypft_frac[p]
                #maxlai = (monthly_lai).max(axis=0)
                for t in range(0,12):
                    if (float(options.lai) > 0):
                      monthly_lai[t][p][0][0] = float(options.lai)
                    #monthly_lai[t][p][j][i] = monthly_lai[t][p][0][0] 
                    #monthly_sai[t][p][j][i] = monthly_sai[t][p][0][0]
                    #monthly_height_top[t][p][j][i] = monthly_height_top[t][p][0][0]
                    #monthly_height_bot[t][p][j][i] = monthly_height_bot[t][p][0][0]

        ierr = nffun.putvar(surffile_new, 'LANDFRAC_PFT', landfrac_pft)
        ierr = nffun.putvar(surffile_new, 'PFTDATA_MASK', pftdata_mask)
        ierr = nffun.putvar(surffile_new, 'LONGXY', longxy)
        ierr = nffun.putvar(surffile_new, 'LATIXY', latixy)
        ierr = nffun.putvar(surffile_new, 'AREA', area)
        ierr = nffun.putvar(surffile_new, 'PCT_WETLAND', pct_wetland)
        ierr = nffun.putvar(surffile_new, 'PCT_LAKE', pct_lake)
        ierr = nffun.putvar(surffile_new, 'PCT_GLACIER',pct_glacier)
        ierr = nffun.putvar(surffile_new, 'PCT_URBAN', pct_urban)
        if (options.mymodel == 'CLM5'):
            ierr = nffun.putvar(surffile_new, 'PCT_CROP', pct_crop)
        else:
            ierr = nffun.putvar(surffile_new, 'SOIL_ORDER', soil_order)
            ierr = nffun.putvar(surffile_new, 'LABILE_P', labilep)
            ierr = nffun.putvar(surffile_new, 'APATITE_P', primp)
            ierr = nffun.putvar(surffile_new, 'SECONDARY_P', secondp)
            ierr = nffun.putvar(surffile_new, 'OCCLUDED_P', occlp)
        ierr = nffun.putvar(surffile_new, 'SOIL_COLOR', soil_color)
        ierr = nffun.putvar(surffile_new, 'FMAX', fmax)
        ierr = nffun.putvar(surffile_new, 'ORGANIC', organic)
        ierr = nffun.putvar(surffile_new, 'PCT_SAND', pct_sand)
        ierr = nffun.putvar(surffile_new, 'PCT_CLAY', pct_clay)
        ierr = nffun.putvar(surffile_new, 'PCT_NATVEG', pct_nat_veg)
        ierr = nffun.putvar(surffile_new, 'PCT_NAT_PFT', pct_pft)
        ierr = nffun.putvar(surffile_new, 'MONTHLY_HEIGHT_TOP', monthly_height_top)
        ierr = nffun.putvar(surffile_new, 'MONTHLY_HEIGHT_BOT', monthly_height_bot)
        ierr = nffun.putvar(surffile_new, 'MONTHLY_LAI', monthly_lai)
    else:
        if (int(options.mypft) >= 0):
          pct_pft      = nffun.getvar(surffile_new, 'PCT_NAT_PFT')
          pct_pft[:,:,:] = 0.0
          pct_pft[int(options.mypft),:,:] = 100.0
          ierr = nffun.putvar(surffile_new, 'PCT_NAT_PFT', pct_pft)

    surffile_list = surffile_list+' '+surffile_new

surffile_new = './temp/surfdata.nc'

if (n_grids > 1):
  os.system('ncecat '+surffile_list+' '+surffile_new)
  os.system('rm ./temp/surfdata?????.nc*')
  #remove ni dimension
  os.system('ncwa -O -a lsmlat -d lsmlat,0,0 '+surffile_new+' '+surffile_new+'.tmp')
  os.system('nccopy -3 -u '+surffile_new+'.tmp'+' '+surffile_new+'.tmp2')
  os.system('ncpdq -a lsmlon,record '+surffile_new+'.tmp2 '+surffile_new+'.tmp3')
  os.system('ncwa -O -a lsmlon -d lsmlon,0,0 '+surffile_new+'.tmp3 '+surffile_new+'.tmp4')
  os.system('ncrename -h -O -d record,gridcell '+surffile_new+'.tmp4 '+surffile_new+'.tmp5')

  os.system('mv '+surffile_new+'.tmp5 '+surffile_new)
  os.system('rm '+surffile_new+'.tmp*')
else:
  os.system('mv '+surffile_list+' '+surffile_new)


#-------------------- create pftdyn surface data ----------------------------------

if (options.nopftdyn == False):

  print('Creating dynpft data')
  pftdyn_list = ''
  for n in range(0,n_grids):
    nst = str(100000+n)[1:]
    pftdyn_new = './temp/surfdata.pftdyn'+nst+'.nc'
    
    if (not os.path.exists(pftdyn_orig)):
        print 'Error: '+pftdyn_orig+' does not exist.  Aborting'
        sys.exit(1)
    if (isglobal):
        os.system('cp '+pftdyn_orig+' '+pftdyn_new)
    else:
        if ('ne' in options.res):
          os.system('ncks --fix_rec_dmn time -d gridcell,'+str(xgrid_min[n])+','+str(xgrid_max[n])+ \
                  ' '+pftdyn_orig+' '+pftdyn_new)
        else:
          os.system('ncks --fix_rec_dmn time -d lsmlon,'+str(xgrid_min[n])+','+str(xgrid_max[n])+ \
                  ' -d lsmlat,'+str(ygrid_min[n])+','+str(ygrid_max[n])+' '+pftdyn_orig+' '+pftdyn_new)
    if (issite):
        landfrac     = nffun.getvar(pftdyn_new, 'LANDFRAC_PFT')
        pftdata_mask = nffun.getvar(pftdyn_new, 'PFTDATA_MASK')
        longxy       = nffun.getvar(pftdyn_new, 'LONGXY')
        latixy       = nffun.getvar(pftdyn_new, 'LATIXY')
        area         = nffun.getvar(pftdyn_new, 'AREA')
        pct_pft      = nffun.getvar(pftdyn_new, 'PCT_NAT_PFT')
        pct_lake_1850    = nffun.getvar(surffile_new, 'PCT_LAKE')
        pct_glacier_1850 = nffun.getvar(surffile_new, 'PCT_GLACIER')
        pct_wetland_1850 = nffun.getvar(surffile_new, 'PCT_WETLAND')
        pct_urban_1850   = nffun.getvar(surffile_new, 'PCT_URBAN')
        pct_pft_1850     = nffun.getvar(surffile_new, 'PCT_NAT_PFT')
        if (options.mymodel == 'CLM5'):
            pct_crop_1850    = nffun.getvar(surffile_new, 'PCT_CROP')
        grazing      = nffun.getvar(pftdyn_new, 'GRAZING')
        harvest_sh1  = nffun.getvar(pftdyn_new, 'HARVEST_SH1')
        harvest_sh2  = nffun.getvar(pftdyn_new, 'HARVEST_SH2')
        harvest_sh3  = nffun.getvar(pftdyn_new, 'HARVEST_SH3')
        harvest_vh1  = nffun.getvar(pftdyn_new, 'HARVEST_VH1')
        harvest_vh2  = nffun.getvar(pftdyn_new, 'HARVEST_VH2')
        
        #read file for site-specific PFT information
        dynexist = False
        mypft_frac=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
        if (options.surfdata_grid == False and options.site != ''):
            AFdatareader = csv.reader(open(ccsm_input+'/lnd/clm2/PTCLM/'+options.sitegroup+'_pftdata.txt','rb'))
            for row in AFdatareader:
                #print(row[0], row[1], options.site)
                if row[0] == options.site:
                    for thispft in range(0,5):
                        mypft_frac[int(row[2+2*thispft])]=float(row[1+2*thispft])

            if (os.path.exists(ccsm_input+'/lnd/clm2/PTCLM/'+options.site+'_dynpftdata.txt')):
                dynexist = True
                DYdatareader = csv.reader(open(ccsm_input+'/lnd/clm2/PTCLM/'+options.site+'_dynpftdata.txt','rb'))
                dim = (19,200)
                pftdata = numpy.zeros(dim)
                for row in DYdatareader:
                    if row[0] == '1850':
                        nrows=1
                        for i in range(0,19):
                            pftdata[i][0] = float(row[i])
                    elif row[0] != 'trans_year':
                        nrows += 1
                        for i in range(0,19):
                            pftdata[i][nrows-1] = float(row[i])
            else:
                print('Warning:  Dynamic pft file for site '+options.site+' does not exist')
                print('Using constant 1850 values')

        landfrac_pft[0][0] = 1.0
        pftdata_mask[0][0] = 1
        if (options.site != ''):
            longxy[0][0] = lon[n]
            latixy[0][0] = lat[n]
            area[0][0] = 111.2*resy*111.321*math.cos((lon[n]*resx)*math.pi/180)*resx
        thisrow = 0
        for t in range(0,nyears_landuse):     
            if (options.surfdata_grid == False):
                if (dynexist):
                    for p in range(0,npft):
                        pct_pft[t][p][0][0] = 0.
                    harvest_thisyear = False
                    if pftdata[0][thisrow+1] == 1850+t:
                        thisrow = thisrow+1
                        harvest_thisyear = True
                    if (t == 0 or pftdata[16][thisrow] == 1):
                        harvest_thisyear = True
                    for k in range(0,5):
                        pct_pft[t][int(pftdata[k*2+2][thisrow])][0][0] = \
                            pftdata[k*2+1][thisrow]
                        grazing[t][0][0] = pftdata[17][thisrow]
                        if (harvest_thisyear):
                            harvest_sh1[t][0][0] = pftdata[13][thisrow]
                            harvest_sh2[t][0][0] = pftdata[14][thisrow]
                            harvest_sh3[t][0][0] = pftdata[15][thisrow]
                            harvest_vh1[t][0][0] = pftdata[11][thisrow]
                            harvest_vh2[t][0][0] = pftdata[12][thisrow]
                        else:
                            harvest_sh1[t][0][0] = 0.
                            harvest_sh2[t][0][0] = 0.
                            harvest_sh3[t][0][0] = 0.
                            harvest_vh1[t][0][0] = 0.
                            harvest_vh2[t][0][0] = 0.
                else:
                    for p in range(0,npft):
                        if (sum(mypft_frac[0:16]) == 0.0):
                            #No dyn file - use 1850 values from gridded file
                            pct_pft[t][p][0][0] = pct_pft_1850[p][n]
                        else:
                            #Use specified 1850 values
                            pct_pft[t][p][0][0] = mypft_frac[p]
                    grazing[t][0][0] = 0.
                    harvest_sh1[t][0][0] = 0.
                    harvest_sh2[t][0][0] = 0.
                    harvest_sh3[t][0][0] = 0.
                    harvest_vh1[t][0][0] = 0.
                    harvest_vh2[t][0][0] = 0.
            else:
                #use time-varying files from gridded file
                print 'using '+surffile_new+' for 1850 information'
                nonpft = float(pct_lake_1850[n]+pct_glacier_1850[n]+ \
                               pct_wetland_1850[n]+pct_urban_1850[n])
                if (options.mymodel == 'CLM5'):
                    nonpft = nonpft+float(pct_crop_1850[n])
                sumpft = 0.0
                pct_pft_temp = pct_pft
                for p in range(0,npft):
                    sumpft = sumpft + pct_pft_temp[t][p][0][0]
                for p in range(0,npft):
                    if (t == 0):
                        #Force 1850 values to surface data file
                        pct_pft[t][p][0][0] = pct_pft_1850[p][n]
                    else:
                        #Scale time-varying values to non-pft fraction
                        #which might not agree with 1850 values
                        #WARNING: - large errors may result if files are inconsistent
                        pct_pft[t][p][0][0] = pct_pft[t][p][0][0]/sumpft*(100.0) #-nonpft)

        ierr = nffun.putvar(pftdyn_new, 'LANDFRAC_PFT', landfrac)
        ierr = nffun.putvar(pftdyn_new, 'PFTDATA_MASK', pftdata_mask)
        ierr = nffun.putvar(pftdyn_new, 'LONGXY', longxy)
        ierr = nffun.putvar(pftdyn_new, 'LATIXY', latixy)
        ierr = nffun.putvar(pftdyn_new, 'AREA', area)
        ierr = nffun.putvar(pftdyn_new, 'PCT_NAT_PFT', pct_pft)
        ierr = nffun.putvar(pftdyn_new, 'GRAZING', grazing)
        ierr = nffun.putvar(pftdyn_new, 'HARVEST_SH1', harvest_sh1)
        ierr = nffun.putvar(pftdyn_new, 'HARVEST_SH2', harvest_sh2)
        ierr = nffun.putvar(pftdyn_new, 'HARVEST_SH3', harvest_sh3)
        ierr = nffun.putvar(pftdyn_new, 'HARVEST_VH1', harvest_vh1)
        ierr = nffun.putvar(pftdyn_new, 'HARVEST_VH2', harvest_vh2)
    pftdyn_list = pftdyn_list+' '+pftdyn_new

  pftdyn_new = './temp/surfdata.pftdyn.nc'
  if (os.path.isfile(pftdyn_new)):
      print('Warning:  Removing existing pftdyn data file')
      os.system('rm -rf '+pftdyn_new)

  if (n_grids > 1):
      os.system('ncecat '+pftdyn_list+' '+pftdyn_new)

      os.system('rm ./temp/surfdata.pftdyn?????.nc*')
      #remove ni dimension
      os.system('ncwa -O -a lsmlat -d lsmlat,0,0 '+pftdyn_new+' '+pftdyn_new+'.tmp')
      os.system('nccopy -3 -u '+pftdyn_new+'.tmp'+' '+pftdyn_new+'.tmp2')
      os.system('ncpdq -a lsmlon,record '+pftdyn_new+'.tmp2 '+pftdyn_new+'.tmp3')
      os.system('ncwa -O -a lsmlon -d lsmlon,0,0 '+pftdyn_new+'.tmp3 '+pftdyn_new+'.tmp4')
      os.system('ncrename -h -O -d record,gridcell '+pftdyn_new+'.tmp4 '+pftdyn_new+'.tmp5')

      os.system('mv '+pftdyn_new+'.tmp5 '+pftdyn_new)
      os.system('rm '+pftdyn_new+'.tmp*')
  else:
      os.system('mv '+pftdyn_list+' '+pftdyn_new)
