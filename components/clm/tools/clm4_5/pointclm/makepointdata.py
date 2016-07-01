#!/usr/bin/env python
import os, sys, csv, time, math
from optparse import OptionParser
import numpy
import netcdf_functions as nffun

print('\n')
print('Makepointdata.py version 0.2')
print('Utillity to create point-level data from 0.5x0.5 gridded datasets')
print('Contact: ricciutodm@ornl.gov')

parser = OptionParser()
parser.add_option("--compset", dest="compset", default='I1850CN', \
                  help = "component set to use (required)")
parser.add_option("--casename", dest="mycase", default="", \
                  help="Full name of case")
parser.add_option("--caseroot", dest="caseroot", default="../..")
parser.add_option("--site", dest="site", default='', \
                  help = '6-character FLUXNET code to run (required)')
parser.add_option("--sitegroup", dest="sitegroup", default="AmeriFlux", \
                  help = "site group to use (default AmeriFlux)")
parser.add_option("--regional", action="store_true", \
                   dest="regional", default=False, \
                   help="Flag for regional run (2x2 or greater)")
parser.add_option("--csmdir", dest="csmdir", default='../../..', \
                  help = "base CESM directory (default = ../../..)")
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
parser.add_option("--clm40", dest="clm40", default=False, \
                      help = "Flag for CLM 4.0", action="store_true")
parser.add_option("--xpts", dest="xpts", default=1, \
                      help = "number of x points (regional only)")
parser.add_option("--ypts", dest="ypts", default=1, \
                      help = "number of y points (regional only)")


(options, args) = parser.parse_args()


csmdir=os.path.abspath(options.csmdir)
ccsm_input = os.path.abspath(options.ccsm_input)

#------------------- get site information ----------------------------------

os.chdir(ccsm_input+'/lnd/clm2/PTCLM/')
AFdatareader = csv.reader(open(options.sitegroup+'_sitedata.txt',"rb"))
for row in AFdatareader:
    if row[0] == options.site:
        lon=float(row[3])
        if (lon < 0):
            lon=360.0+float(row[3]) 
        lat=float(row[4])
        startyear=int(row[6])
        endyear=int(row[7])
        alignyear = int(row[8])
        if options.regional == True:
           numxpts=int(options.xpts)
           numypts=int(options.ypts)
           resx=float(row[11])
           resy=float(row[12])
        else:
           numxpts=1
           numypts=1
           resx=0.1      #longitudinal resolution (degrees) 
           resy=0.1      #latitudinal resolution (degrees)
           #lat=lat-resy/2
           #lon=lon-resx/2
        if (options.makemet):
            print(" Making meteorological data for site")
            metcmd = 'python '+csmdir+'/components/clm/tools/clm4_5/pointclm/makemetdata.py' \
                          +' --site '+options.site+' --lat '+row[4]+' --lon '+ \
                          row[3]+' --ccsm_input '+ccsm_input+ \
                          ' --startyear '+row[6]+' --endyear '+row[7]+' --numxpts '+ \
                          str(numxpts)+' --numypts '+str(numypts)
            if (options.metdir != 'none'):
                metcmd = metcmd + ' --metdir '+options.metdir
            os.system(metcmd)
        else:
            print('Met data not requested.  Model will not run if data do not exist')

#get corresponding 0.5x0.5 and 1.9x2.5 degree grid cells
if (lon < 0):
    xgrid        = 720+int(lon/0.5)
    xgrid_lores  = 144+int(lon/2.5)
    xgrid_T62    = 192+int(lon/1.875)
else:
    xgrid        = int(lon*2)
    xgrid_lores  = int(lon/2.5)
    xgrid_T62    = int(lon/1.875)
ygrid         = int((lat+90)/0.5)
ygrid_lores   = int((lat+90)/1.9)
ygrid_T62     = int((lat+90)/1.9)

#---------------------Create domain data --------------------------------------------------

print('Creating domain data')
os.system('mkdir -p '+csmdir+'/components/clm/tools/clm4_5/pointclm/temp')
domainfile_orig = ccsm_input+'/atm/datm7/domain.clm/' \
    +'domain.360x720_ORCHIDEE0to360.100409.nc'
domainfile_new = csmdir+'/components/clm/tools/clm4_5/pointclm/temp/' \
    +'domain.lnd.'+str(numxpts)+'x'+str(numypts)+'pt_'+options.site+'_navy.nc'
if (os.path.isfile(domainfile_new)):
    print('Warning:  Removing existing domain file')
    os.system('rm -rf '+domainfile_new)
os.system('ncks -d ni,'+str(xgrid)+','+str(xgrid+numxpts-1)+' -d nj,'+str(ygrid)+ \
          ','+str(ygrid+numypts-1)+' '+domainfile_orig+' '+domainfile_new)
#domainfile_new_nc = NetCDF.NetCDFFile(domainfile_new, "a")

frac = nffun.getvar(domainfile_new, 'frac')
mask = nffun.getvar(domainfile_new, 'mask')
xc = nffun.getvar(domainfile_new, 'xc')
yc = nffun.getvar(domainfile_new, 'yc')
xv = nffun.getvar(domainfile_new, 'xv')
yv = nffun.getvar(domainfile_new, 'yv')
area = nffun.getvar(domainfile_new, 'area')

for i in range(0,numxpts):
    for j in range(0,numypts):
        frac[j][i] = 1.0
        mask[j][i] = 1
        xc[j][i] = lon+i*resx
        yc[j][i] = lat+j*resy
        #print(i,j,lon,lat)
        xv[j][i][0] = lon-resx/2+i*resx
        xv[j][i][1] = lon+resx/2+i*resx
        xv[j][i][2] = lon-resx/2+i*resx
        xv[j][i][3] = lon+resx/2+i*resx
        yv[j][i][0] = lat-resy/2+j*resy
        yv[j][i][1] = lat-resy/2+j*resy
        yv[j][i][2] = lat+resy/2+j*resy
        yv[j][i][3] = lat+resy/2+j*resy
        area[j][i] = resx*resy*math.pi/180*math.pi/180

ierr = nffun.putvar(domainfile_new, 'frac', frac)
ierr = nffun.putvar(domainfile_new, 'mask', mask)
ierr = nffun.putvar(domainfile_new, 'xc', xc)
ierr = nffun.putvar(domainfile_new, 'yc', yc)
ierr = nffun.putvar(domainfile_new, 'xv', xv)
ierr = nffun.putvar(domainfile_new, 'yv', yv)
ierr = nffun.putvar(domainfile_new, 'area', area)

#os.system('cp 'ccsm_input+'/share/domains/domain.clm/domain.lnd.'+str(numxpts)+'x'+str(numypts)+ \
#              'pt_'+options.site+'_navy.nc '+ccsm_input+'/atm/datm7/domain.clm')

#-------------------- create surface data ----------------------------------
print('Creating surface data')
mysimyr=1850
if (options.compset == 'I2000CN'):
    mysimyr=2000
#surffile_orig = ccsm_input+'/lnd/clm2/surfdata/' \
#    +'surfdata_0.5x0.5_simyr1850.nc'
if (options.clm40):
    surffile_orig = ccsm_input+'/lnd/clm2/surfdata/surfdata_360x720_nourb_simyr1850_c120717.nc'
    surffile_new  = csmdir+'/components/clm/tools/clm4_5/pointclm/temp/surfdata_'+str(numxpts)+'x' \
                    +str(numypts)+'pt_'+options.mycase+'_simyr'+str(mysimyr)+'.nc'
else:
    #surffile_orig = ccsm_input+'/lnd/clm2/surfdata_map/surfdata_360x720cru_simyr1850_c130415.nc'
    surffile_orig = ccsm_input+'/lnd/clm2/surfdata_map/surfdata_360x720cru_simyr1850_c150626.nc'
    surffile_new =  csmdir+'/components/clm/tools/clm4_5/pointclm/temp/surfdata_'+str(numxpts)+'x'+str(numypts)+'pt_'+ \
                    options.mycase+'_simyr'+str(mysimyr)+'.nc'

print(surffile_orig)
if (os.path.isfile(surffile_new)):
    print('Warning:  Removing existing surface file')
    os.system('rm -rf '+surffile_new)
print('ncks --fix_rec_dmn time -d lsmlon,'+str(xgrid)+','+str(xgrid+numxpts-1)+' -d lsmlat,'+str(ygrid)+ \
          ','+str(ygrid+numypts-1)+' '+surffile_orig+' '+surffile_new)
os.system('ncks --fix_rec_dmn time -d lsmlon,'+str(xgrid)+','+str(xgrid+numxpts-1)+' -d lsmlat,'+str(ygrid)+ \
          ','+str(ygrid+numypts-1)+' '+surffile_orig+' '+surffile_new)
    
#surffile_new_nc = NetCDF.NetCDFFile(surffile_new, "a")
#if (options.clm40 == False):

landfrac_pft = nffun.getvar(surffile_new, 'LANDFRAC_PFT')
pftdata_mask = nffun.getvar(surffile_new, 'PFTDATA_MASK')
longxy       = nffun.getvar(surffile_new, 'LONGXY')
latixy       = nffun.getvar(surffile_new, 'LATIXY')
area         = nffun.getvar(surffile_new, 'AREA')
pct_wetland  = nffun.getvar(surffile_new, 'PCT_WETLAND')
pct_lake     = nffun.getvar(surffile_new, 'PCT_LAKE')
pct_glacier  = nffun.getvar(surffile_new, 'PCT_GLACIER')
pct_urban    = nffun.getvar(surffile_new, 'PCT_URBAN')

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

#read file for site-specific PFT information
if (options.surfdata_grid == False):
    AFdatareader = csv.reader(open(options.sitegroup+'_pftdata.txt','rb'))
    mypft_frac=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
    for row in AFdatareader:
        if row[0] == options.site:
            for thispft in range(0,5):
                mypft_frac[int(row[2+2*thispft])]=float(row[1+2*thispft])
    if (sum(mypft_frac[0:17]) == 0.0):
        print('*** Warning:  PFT data NOT found.  Using gridded data ***')


#read file for site-specific soil information
    mypct_sand = 0.0
    mypct_clay = 0.0
    AFdatareader = csv.reader(open(options.sitegroup+'_soildata.txt','rb'))
    for row in AFdatareader:
    	if row[0] == options.site:
            mypct_sand = row[4]
            mypct_clay = row[5]
    if (mypct_sand == 0.0 and mypct_clay == 0.0):
        print('*** Warning:  Soil data NOT found.  Using gridded data ***')

for i in range(0,numxpts):
    for j in range(0,numypts):
#        if (options.clm40 == False):
        landfrac_pft[j][i] = 1.0
        pftdata_mask[j][i] = 1
        longxy[j][i] = lon+i*resx
        latixy[j][i] = lat+j*resy
        area[j][i] = 111.2*resy*111.321*math.cos((lon+i*resx)*math.pi/180)*resx
        if ((not options.surfdata_grid)):
            pct_wetland[j][i] = 0.0
            pct_lake[j][i]    = 0.0
            pct_glacier[j][i] = 0.0
            pct_nat_veg[j][i] = 100.0
            if (options.clm40):
                pct_urban[j][i] = 0.0
            else:
                for k in range(0,3):
                    pct_urban[k][j][i] = 0.0
        else:
            pct_wetland[j][i] = pct_wetland[0][0]
            pct_lake[j][i]    = pct_lake[0][0]
            pct_glacier[j][i] = pct_glacier[0][0]      
            pct_nat_veg[j][i] = pct_nat_veg[0][0]      
            if (options.clm40):
                pct_urban[j][i] = pct_urban[0][0]
            else:
                for k in range(0,3):
                    pct_urban[k][j][i]   = pct_urban[k][0][0]
        soil_color[j][i] = soil_color[0][0]
        fmax[j][i] = fmax[0][0]
        for k in range(0,10):
            organic[k][j][i] = organic[k][0][0]
            if (options.site == 'US-SPR'):
                if (k < 8):
                    organic[k][j][i] = 130.
                elif (k == 8):
                    organic[k][j][i] = 65.
                elif (k == 9):
                    organic[k][j][i] = 0.
            if (options.surfdata_grid == False and (mypct_sand > 0.0 or mypct_clay > 0.0)):
                pct_sand[k][j][i]   = mypct_sand
                pct_clay[k][j][i]   = mypct_clay
            else:
                pct_sand[k][j][i]   = pct_sand[k][0][0]
                pct_clay[k][j][i]   = pct_clay[k][0][0]
        for p in range(0,npft):
            if (options.surfdata_grid == False and sum(mypft_frac[0:17]) > 0.0):
                pct_pft[p][j][i] = mypft_frac[p]
            else:
                pct_pft[p][j][i] = pct_pft[p][0][0]
            
            maxlai = (monthly_lai).max(axis=0)
            for t in range(0,12):
                monthly_lai[t][p][j][i] = monthly_lai[t][p][0][0] 
                if (options.site == 'US-Ne3'):
                    if (monthly_lai[t][p][j][i] > 0):
                        monthly_lai[t][p][j][i] = monthly_lai[t][p][0][0]/maxlai[p]
                        if (t <= 3 or t >= 9):
                            monthly_lai[t][p][j][i] = 0.1 #0.2
                        elif (t > 3 and t <= 6):
                            monthly_lai[t][p][j][i] = (t-3)*2.5 #0.8
                        elif (t == 7):
                            monthly_lai[t][p][j][i] = 7.5 #0.75
                        elif (t == 8):
                            monthly_lai[t][p][j][i] = 3.75 #0.3
                        #print p, t, monthly_lai[t][p][j][i], maxlai[p]
                monthly_sai[t][p][j][i] = monthly_sai[t][p][0][0]
                monthly_height_top[t][p][j][i] = monthly_height_top[t][p][0][0]
                monthly_height_bot[t][p][j][i] = monthly_height_bot[t][p][0][0]

#if (options.clm40 == False):
ierr = nffun.putvar(surffile_new, 'LANDFRAC_PFT', landfrac_pft)
ierr = nffun.putvar(surffile_new, 'PFTDATA_MASK', pftdata_mask)
ierr = nffun.putvar(surffile_new, 'LONGXY', longxy)
ierr = nffun.putvar(surffile_new, 'LATIXY', latixy)
ierr = nffun.putvar(surffile_new, 'AREA', area)
ierr = nffun.putvar(surffile_new, 'PCT_WETLAND', pct_wetland)
ierr = nffun.putvar(surffile_new, 'PCT_LAKE', pct_lake)
ierr = nffun.putvar(surffile_new, 'PCT_GLACIER',pct_glacier)
ierr = nffun.putvar(surffile_new, 'PCT_URBAN', pct_urban)
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

#os.system('cp '+ccsm_input+'/lnd/clm2/surfdata_map/surfdata_'+str(numxpts)+'x'+ \
#              str(numypts)+'pt_'+options.mycase+'_simyr'+str(mysimyr)+'.nc '+ccsm_input+ \
#              '/lnd/clm2/surfdata')

#-------------------- create pftdyn surface data ----------------------------------

if ('20TR' in options.compset):

    print('Creating dynpft data')


    #pftdyn_orig = ccsm_input+'/lnd/clm2/surfdata/' \
    #  +'surfdata.pftdyn_0.5x0.5_simyr1850-2010_ACME.nc'
    pftdyn_orig = ccsm_input+'/lnd/clm2/surfdata_map/' \
       +'landuse.timeseries_360x720cru_rcp4.5_simyr1850-2100_c141219.nc'

    print 'using '+surffile_new+' for 1850 information'
    if (options.clm40):
        pftdyn_new = csmdir+'/components/clm/tools/clm4_5/pointclm/temp/' \
          +'surfdata.pftdyn_'+str(numxpts)+'x'+str(numypts)+'pt_'+options.mycase+'.nc'
    else:
        pftdyn_new = csmdir+'/components/clm/tools/clm4_5/pointclm/temp/' \
          +'surfdata.pftdyn_'+str(numxpts)+'x'+str(numypts)+'pt_'+options.mycase+'.nc'
        
    print pftdyn_new
    if (os.path.isfile(pftdyn_new)):
        print('Warning:  Removing existing pftdyn file')
        os.system('rm -rf '+pftdyn_new)
    os.system('ncks --fix_rec_dmn time -d lsmlon,'+str(xgrid)+','+str(xgrid+numxpts-1)+' -d lsmlat,'+str(ygrid)+ \
                  ','+str(ygrid+numypts-1)+' '+pftdyn_orig+' '+pftdyn_new)
    print('ncks --fix_rec_dmn time -d lsmlon,'+str(xgrid)+','+str(xgrid+numxpts-1)+' -d lsmlat,'+str(ygrid)+ \
                  ','+str(ygrid+numypts-1)+' '+pftdyn_orig+' '+pftdyn_new)
    #pftdyn_new_nc = NetCDF.NetCDFFile(pftdyn_new, "a")
    landfrac     = nffun.getvar(pftdyn_new, 'LANDFRAC_PFT')
    pftdata_mask = nffun.getvar(pftdyn_new, 'PFTDATA_MASK')
    longxy       = nffun.getvar(pftdyn_new, 'LONGXY')
    latixy       = nffun.getvar(pftdyn_new, 'LATIXY')
    area         = nffun.getvar(pftdyn_new, 'AREA')
    #pct_wetland  = nffun.getvar(pftdyn_new, 'PCT_WETLAND')
    #pct_lake     = nffun.getvar(pftdyn_new, 'PCT_LAKE')
    #pct_glacier  = nffun.getvar(pftdyn_new, 'PCT_GLACIER')
    #pct_urban    = nffun.getvar(pftdyn_new, 'PCT_URBAN')
    pct_pft      = nffun.getvar(pftdyn_new, 'PCT_NAT_PFT')
    #pct_lake_1850    = nffun.getvar(surffile_new, 'PCT_LAKE')
    #pct_glacier_1850 = nffun.getvar(surffile_new, 'PCT_GLACIER')
    #pct_wetland_1850 = nffun.getvar(surffile_new, 'PCT_WETLAND')
    #pct_urban_1850   = nffun.getvar(surffile_new, 'PCT_URBAN')
    pct_pft_1850     = nffun.getvar(surffile_new, 'PCT_NAT_PFT')
    grazing      = nffun.getvar(pftdyn_new, 'GRAZING')
    harvest_sh1  = nffun.getvar(pftdyn_new, 'HARVEST_SH1')
    harvest_sh2  = nffun.getvar(pftdyn_new, 'HARVEST_SH2')
    harvest_sh3  = nffun.getvar(pftdyn_new, 'HARVEST_SH3')
    harvest_vh1  = nffun.getvar(pftdyn_new, 'HARVEST_VH1')
    harvest_vh2  = nffun.getvar(pftdyn_new, 'HARVEST_VH2')

    npft = 17

    #read file for site-specific PFT information
    AFdatareader = csv.reader(open(options.sitegroup+'_pftdata.txt','rb'))
    mypft_frac=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
    for row in AFdatareader:
        #print(row[0], row[1], options.site)
        if row[0] == options.site:
            for thispft in range(0,5):
                mypft_frac[int(row[2+2*thispft])]=float(row[1+2*thispft])

    if (os.path.exists(options.site+'_dynpftdata.txt')):
        dynexist = True
        DYdatareader = csv.reader(open(options.site+'_dynpftdata.txt','rb'))
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
        dynexist = False
        print('Warning:  Dynamic pft file for site '+options.site+' does not exist')
        print('Using constant 1850 values')

    for i in range(0,numxpts):
        for j in range(0,numypts):
            landfrac_pft[j][i] = 1.0
            pftdata_mask[j][i] = 1
            longxy[j][i] = lon+i*resx
            latixy[j][i] = lat+j*resy
            area[j][i] = 111.2*resy*111.321*math.cos((lon+i*resx)*math.pi/180)*resx
            #if (not options.surfdata_grid):
            #    pct_wetland[j][i] = 0.0
            #    pct_lake[j][i]    = 0.0
            #    pct_glacier[j][i] = 0.0
            #    pct_urban[j][i]   = 0.0
            #else:
            #    pct_wetland[j][i] = pct_wetland_1850[0][0]
            #    pct_lake[j][i] = pct_lake_1850[0][0]
            #    pct_glacier[j][j] = pct_glacier_1850[0][0]
            #    pct_urban[j][i] = pct_urban_1850[0][0]
            thisrow = 0
            for t in range(0,161):     
                if (options.surfdata_grid == False):
                    if (dynexist):
                        for p in range(0,npft):
                            pct_pft[t][p][j][i] = 0.
                        harvest_thisyear = False
                        if pftdata[0][thisrow+1] == 1850+t:
                            thisrow = thisrow+1
                            harvest_thisyear = True
                        if (t == 0 or pftdata[16][thisrow] == 1):
                            harvest_thisyear = True
                        for k in range(0,5):
                            pct_pft[t][int(pftdata[k*2+2][thisrow])][j][i] = \
                                pftdata[k*2+1][thisrow]
                            grazing[t][j][i] = pftdata[17][thisrow]
                            if (harvest_thisyear):
                                harvest_sh1[t][j][i] = pftdata[13][thisrow]
                                harvest_sh2[t][j][i] = pftdata[14][thisrow]
                                harvest_sh3[t][j][i] = pftdata[15][thisrow]
                                harvest_vh1[t][j][i] = pftdata[11][thisrow]
                                harvest_vh2[t][j][i] = pftdata[12][thisrow]
                            else:
                                harvest_sh1[t][j][i] = 0.
                                harvest_sh2[t][j][i] = 0.
                                harvest_sh3[t][j][i] = 0.
                                harvest_vh1[t][j][i] = 0.
                                harvest_vh2[t][j][i] = 0.
                    else:
                        for p in range(0,npft):
                            if (sum(mypft_frac[0:16]) == 0.0):
                                print p, pct_pft[0][p][0][0]
                                #No dyn file - use 1850 values from gridded file
                                pct_pft[t][p][j][i] = pct_pft_1850[p][0][0]
                                
                            else:
                                #Use specified 1850 values
                                pct_pft[t][p][j][i] = mypft_frac[p]
                        grazing[t][j][i] = 0.
                        harvest_sh1[t][j][i] = 0.
                        harvest_sh2[t][j][i] = 0.
                        harvest_sh3[t][j][i] = 0.
                        harvest_vh1[t][j][i] = 0.
                        harvest_vh2[t][j][i] = 0.
                else:
                    #use time-varying files from gridded file
                    nonpft = float(pct_lake_1850[0][0]+pct_glacier_1850[0][0]+ \
                                       pct_wetland_1850[0][0]+pct_urban_1850[0][0])
                    sumpft = 0.0
                    pct_pft_temp = pct_pft
                    for p in range(0,npft):
                        sumpft = sumpft + pct_pft_temp[t][p][0][0]
                    for p in range(0,npft):
                        if (t == 0):
                            #Force 1850 values to surface data file
                            pct_pft[t][p][j][i] = pct_pft_1850[p][0][0]
                        else:
                            #Scale time-varying values to non-pft fraction
                            #which might not agree with 1850 values
                            #WARNING: - large errors may result if files are inconsistent
                            nonpft = float(pct_lake_1850[0][0]+pct_glacier_1850[0][0]+ \
                                pct_wetland_1850[0][0]+pct_urban_1850[0][0])
                            pct_pft[t][p][j][i] = pct_pft[t][p][0][0]/sumpft*(100.0) #-nonpft)
                    harvest_sh1[t][j][i] = harvest_sh1[t][0][0]
                    harvest_sh2[t][j][i] = harvest_sh2[t][0][0]
                    harvest_sh3[t][j][i] = harvest_sh3[t][0][0]
                    harvest_vh1[t][j][i] = harvest_vh1[t][0][0]
                    harvest_vh2[t][j][i] = harvest_vh2[t][0][0]
                    grazing[t][j][i] = grazing[t][0][0]


    ierr = nffun.putvar(pftdyn_new, 'LANDFRAC_PFT', landfrac)
    ierr = nffun.putvar(pftdyn_new, 'PFTDATA_MASK', pftdata_mask)
    ierr = nffun.putvar(pftdyn_new, 'LONGXY', longxy)
    ierr = nffun.putvar(pftdyn_new, 'LATIXY', latixy)
    ierr = nffun.putvar(pftdyn_new, 'AREA', area)
    #ierr = nffun.putvar(pftdyn_new, 'PCT_WETLAND', pct_wetland)
    #ierr = nffun.putvar(pftdyn_new, 'PCT_LAKE', pct_lake)
    #ierr = nffun.putvar(pftdyn_new, 'PCT_GLACIER', pct_glacier)
    #ierr = nffun.putvar(pftdyn_new, 'PCT_URBAN', pct_urban)
    ierr = nffun.putvar(pftdyn_new, 'PCT_NAT_PFT', pct_pft)
    ierr = nffun.putvar(pftdyn_new, 'GRAZING', grazing)
    ierr = nffun.putvar(pftdyn_new, 'HARVEST_SH1', harvest_sh1)
    ierr = nffun.putvar(pftdyn_new, 'HARVEST_SH2', harvest_sh2)
    ierr = nffun.putvar(pftdyn_new, 'HARVEST_SH3', harvest_sh3)
    ierr = nffun.putvar(pftdyn_new, 'HARVEST_VH1', harvest_vh1)
    ierr = nffun.putvar(pftdyn_new, 'HARVEST_VH2', harvest_vh2)

