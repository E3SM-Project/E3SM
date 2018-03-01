# Copyright (C) 2004, 2005, 2009, 2010
# Glimmer-CISM contributors - see AUTHORS file for list of contributors
#
# This file is part of Glimmer-CISM.
#
# Glimmer-CISM is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or (at
# your option) any later version.
#
# Glimmer-CISM is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Glimmer-CISM.  If not, see <http://www.gnu.org/licenses/>.
#
# Glimmer-CISM is hosted on BerliOS.de:
# https://developer.berlios.de/projects/glimmer-cism/

"""Loading CF files."""

__all__=['GCloadfile','GCvariable','GCchecklist']

from GC_netcdf import GCNetCDFFile
import numpy, os
from GC_proj import *
#MH#from GC_colourmap import *
from GC_file import *
from GC_createfile import *
#MH#from TwoDspline import TwoDspline
import scipy.ndimage

temperatures = ['btemp','temp']

def GCchecklist(section,variable):
    """Check if section is a list.

    if variable is a None,     return (True,[0,len(variable)-1])
    if variable is a list/etc, return (True,[section[0],section[1]])
    if variable is a value,    return (False,val) """

    if section is None:
        return (True, [0,len(variable)-1])
    elif type(section) == list or type(section) == tuple or type(section) == numpy.ndarray:
        return (True, [section[0],section[1]])
    else:
        return (False,section)

class GCloadfile(GCfile):
    """Loading a CF netCDF file."""

    def __init__(self,fname):
        """Initialise.

        fname: name of CF file."""

        GCfile.__init__(self,fname)

        self.file = GCNetCDFFile(self.fname,'r')

        self.timescale = 0.001
        # get mapping variable name
        for var in self.file.variables.keys():
            if hasattr(self.file.variables[var],'grid_mapping_name'):
                self.mapvarname = var
                break
        self.reset_bb()
        # initialising variable dictionary
        self.__vars = {}
        # RSL residuals
        self.__rslres = {}

    def time(self,t):
        """Return selected time value."""

        (isar,sel) = GCchecklist(t,self.file.variables['time'])

        if isar:
            return self.file.variables['time'][sel[0]:sel[1]+1]*self.timescale
        else:
            return self.file.variables['time'][sel]*self.timescale

    def timeslice(self,time,round='n'):
        """Get the time slice.

        time: time to look up in ISM file
        round: 'n' round to nearest
               'u' round up
               'd' round down"""

        if round not in ['n','u','d']:
            raise ValueError, "Expected one of 'n', 'u', 'd'"

        t0 = 0
        t1 = len(self.file.variables['time'][:])-1
        if time < self.time(t0) or time > self.time(t1):
            raise ValueError, 'Selected time slice [%f] is outside file %s: [%f, %f]'%(time,self.fname,self.time(t0),self.time(t1))
        if time == self.time(t0): return t0
        if time == self.time(t1): return t1
        # use Newton bisection
        tmid = int((t1-t0)/2)
        while tmid > 0:
            if time < self.time(t0+tmid):
                t1 = t0+tmid
            elif time > self.time(t0+tmid):
                t0 = t0+tmid
            else:
                return t0+tmid
            tmid = int((t1-t0)/2)
        if round == 'u':
            return t1
        elif round == 'd':
            return t0
        else:
            if (time-self.time(t0)) < (self.time(t1) - time):
                return t0
            else:
                return t1
        raise AssertionError, 'Why are we here?'

    def getvar(self,var):
        """Get a variable from file.

        var: name of variables

        this method caches the return variable structure."""

        if var not in self.__vars:
            self.__vars[var] = GCvariable(self,var)
        return self.__vars[var]

    def getIceArea(self,time=None,scale=1.):
        """Get area covered by ice.
        
        time: if None, return data for all time slices
              if list/etc of size two, interpret as array selection
              if single value, get only this time slice"""

        (tarray,t) = GCchecklist(time,self.file.variables['time'])
        values = []
        fact = self.deltax*self.deltay*scale
        if tarray:
            for i in range(t[0],t[1]+1):
                ih = numpy.where(self.file.variables['thk'][i,:,:]>0.,1,0).flat
                values.append(sum(ih)*fact)
            return values
        ih = numpy.where(self.file.variables['thk'][t,:,:]>0.,1,0).flat
        return sum(ih)*fact

    def getIceVolume(self,time=None,scale=1.):
        """Get ice volume
        
        time: if None, return data for all time slices
              if list/etc of size two, interpret as array selection
              if single value, get only this time slice"""

        (tarray,t) = GCchecklist(time,self.file.variables['time'])
        values = []
        fact = self.deltax*self.deltay*scale
        if tarray:
            for i in range(t[0],t[1]+1):
                ih = numpy.where(self.file.variables['thk'][i,:,:]>0.,self.file.variables['thk'][i,:,:],0.).flat
                values.append(sum(ih)*fact)
            return values
        ih = self.file.variables['thk'][t,:,:].flat
        return sum(ih)*fact

    def getFracMelt(self,time=None,scale=1.):
        """Get fractional area where basal melting occurs.

        time: if None, return data for all time slices
              if list/etc of size two, interpret as array selection
              if single value, get only this time slice"""

        (tarray,t) = GCchecklist(time,self.file.variables['time'])
        values = []
        fact = self.deltax*self.deltay*scale
        if tarray:
            for i in range(t[0],t[1]+1):
                ih = self.getIceArea(time=i,scale=scale)
                if ih>0:
                    mlt = numpy.where(self.file.variables['bmlt'][i,:,:]>0.,1,0).flat
                    values.append(sum(mlt)*fact/ih)
                else:
                    values.append(0.)
            return values
        ih = self.getIceArea(time=t,scale=scale)
        if ih>0:
            mlt = numpy.where(self.file.variables['bmlt'][t,:,:]>0.,1,0).flat
            return sum(mlt)*fact/ih
        else:
            return 0.

##     def getRSL(self,loc,time,clip=True):
##         """Get RSL data.

##         loc: array,list,tuple containing longitude and latitude of RSL location
##         time: if None, return data for all time slices
##               if list/etc of size two, interpret as array selection
##               if single value, get only this time slice
##         clip: if set to true only extract RSL for ice free locations"""

##         # get times
##         (tarray,t) = GCchecklist(time,self.file.variables['time'])
##         # get location
##         xyloc = self.project(loc)
##         if not self.inside(xyloc):
##             raise RuntimeError, 'Point outside grid'
##         data = self.getvar('isobase')
##         # extract data
##         values = []
##         if tarray:
##             if clip:
##                 ih_data = self.getvar('thk')
##                 for i in range(t[0],t[1]+1):
##                     ih = ih_data.spline(xyloc,i)
##                     if ih>0.:
##                         values.append('nan')
##                     else:
##                         values.append(data.spline(xyloc,i))
##             else:
##                 for i in range(t[0],t[1]+1):
##                     values.append(data.spline(xyloc,i))
##             return values

##         return data.spline(xyloc,t)

##     def getRSLresiduals(self,rsldb,time=None):
##         """Get RSL residuals.

##         rsldb: RSL data base
##         time: time interval to be processed"""

##         hnx = 50
##         hny = 50
        
##         # get times
##         if time==None:
##             t = [self.timeslice(rsldb.mint*self.timescale,'d'),self.timeslice(0.)]
##         else:
##             t = [self.timeslice(time[0],'d'),self.timeslice(time[1],'u')]
##         times = self.time(t)
        
##         # loop over locations
##         res_times = []
##         residuals = []
##         for loc in rsldb.getLocationRange(self.minmax_long,self.minmax_lat):
##             try:
##                 res = self.get_rslres(rsldb,loc[0])
##             except:
##                 continue
##             for i in range(0,len(res[0])):
##                 res_times.append(res[0][i])
##                 residuals.append(res[1][i])
##         # create histogram
##         hist = histogram.histogram2d(hnx,hny)
##         hist.set_ranges_uniform(times[0],times[-1],PyGMT.round_down(min(residuals)),PyGMT.round_up(max(residuals)))

##         for i in range(0,len(residuals)):
##             hist.increment(res_times[i],residuals[i])
##         # turn into a grid
##         grid = PyGMT.Grid()
##         grid.x_minmax = [times[0],times[-1]]
##         grid.y_minmax = [PyGMT.round_down(min(residuals)),PyGMT.round_up(max(residuals))]
##         grid.data=numpy.zeros([hnx,hny],'f')
##         for j in range(0,hny):
##             for i in range(0,hnx):
##                 grid.data[i,j] = hist[i,j]
##         return grid

##     def get_rslres(self,rsldb,lid,avg=False):
##         """Get RSL residual.

##         rsldb: RSL database
##         lid: location id.
##         avg: set to True to get average"""

##         # check if residuals are cached
##         if lid not in self.__rslres:
##             # get coordinates
##             cu = rsldb.db.cursor()
##             cu.execute('SELECT longitude,latitude FROM location WHERE location_id == %i',(lid))
##             loc = cu.fetchone()

##             # get data
##             times = []
##             obs = []
##             cu.execute('SELECT time,rsl FROM measurement WHERE location_id == %i',(lid))
##             for o in cu.fetchall():
##                 times.append(o[0]*self.timescale)
##                 obs.append(o[1])
##             ti = [self.timeslice(min(times)-2.,'d'), self.timeslice(max(times),'u')]
##             ts = spline.cspline(ti[1]-ti[0]+1)
##             ts.init(self.time(ti),self.getRSL(list(loc),ti,clip=False))
##             residuals = []
##             for i in range(0,len(times)):
##                 residuals.append(obs[i]-ts.eval(times[i]))
##             self.__rslres[lid] = (times,residuals)
##         if avg:
##             r = 0.
##             if len(self.__rslres[lid][1])>0:
##                 r = sum(self.__rslres[lid][1])/float(len(self.__rslres[lid][1]))
##             return r
##         else:
##             return self.__rslres[lid]

    def clone(self,fname):
        """Clone self.

        create a new CF file with name fname and
        copy dimensions, mapping and global metadata."""

        newcf = GCcreatefile(fname)
        # copy global attributes
        for attrib in ['title','institution','source','references','comment','history']:
            if hasattr(self,attrib):
                setattr(newcf,attrib,getattr(self,attrib))
        # create dimensions
        for dim in self.file.dimensions.keys():
            newcf.createDimension(dim,self.file.dimensions[dim])
            # create dim variables
            var = newcf.createVariable(dim)
            if dim != 'time':
                var[:] = self.file.variables[dim][:]
        # copy mapping
        if self.mapvarname in self.file.variables.keys():
            varmap=newcf.file.createVariable(self.mapvarname,'c',())
            copyGCMap(self.file.variables[self.mapvarname],varmap)

        return newcf
              
class GCvariable(object):
    """Handling CF variables."""

    def __init__(self,cffile,var):
        """Initialise.

        CFFile: CF file
        var: name of variable"""

        self.cffile = cffile
        self.file = cffile.file
        if var[-4:] == '_avg':
            self.name = var[:-4]
            self.average = True
        else:
            self.name = var
            self.average = False
        if self.name=='is':
            if 'topg' not in self.file.variables.keys() or 'thk' not in self.file.variables.keys():
                raise KeyError, 'Variable not in file'
        elif self.name=='isobase':
            if 'slc' not in self.file.variables.keys():
                raise KeyError, 'Variable not in file'
        elif self.name=='pmp':
            if 'thk' not in self.file.variables.keys():
                raise KeyError, 'Variable not in file'
        elif self.name=='vel':
            if 'uvel' not in self.file.variables.keys() or 'vvel' not in self.file.variables.keys():
                raise KeyError, 'Variable not in file'
        elif self.name=='bvel':
            if 'ubas' not in self.file.variables.keys() or 'vbas' not in self.file.variables.keys():
                raise KeyError, 'Variable not in file'
        elif self.name=='bvel_tavg':
            if 'ubas_tavg' not in self.file.variables.keys() or 'vbas_tavg' not in self.file.variables.keys():
                raise KeyError, 'Variable not in file'
        elif self.name=='tau':
            if 'taux' not in self.file.variables.keys() or 'tauy' not in self.file.variables.keys():
                raise KeyError, 'Variable not in file'
        elif self.name not in self.file.variables.keys():
            raise KeyError, 'Variable not in file: %s'%self.name
#MH#        self.__colourmap = GCcolourmap(self)
        self.pmt = False
        self.__varcache = None

    def __get_units(self):
        try:
            if self.name == 'is':
                return self.file.variables['topg'].units
            elif self.name=='isobase':
                return self.file.variables['slc'].units
            elif self.name == 'pmp':
                return 'degree_Celsius'
            elif self.name == 'vel':
                return self.file.variables['uvel'].units
            elif self.name == 'bvel':
                return self.file.variables['ubas'].units
            elif self.name == 'bvel_tavg':
                return self.file.variables['ubas_tavg'].units
            elif self.name == 'tau':
                return self.file.variables['taux'].units
            else:
                return self.file.variables[self.name].units
        except:
            return ''
    units = property(__get_units)

    def __get_long_name(self):
        try:
            if self.name in temperatures and self.pmt and 'thk' in self.file.variables.keys():
                name = 'homologous %s'%self.file.variables[self.name].long_name
            elif self.name == 'is':
                name =  'ice surface elevation'
            elif self.name=='isobase':
                name = 'isobase'
            elif self.name == 'pmp':
                name = 'pressure melting point of ice'
            elif self.name == 'vel':
                name = 'horizontal velocity'
            elif self.name == 'bvel':
                name = 'horizontal basal velocity'
            elif self.name == 'bvel_tavg':
                name = 'horizontal basal velocity (time average)'
            elif self.name == 'tau':
                name = 'basal shear stress'
            else:
                name = self.file.variables[self.name].long_name
        except:
            name = ''
        if self.average:
            name = 'vertically averaged %s'%name
        return name
    long_name = property(__get_long_name)

    def __get_standard_name(self):
        try:
            return self.file.variables[self.name].standard_name
        except:
            return ''
    standard_name = property(__get_standard_name)

    def __get_xdimension(self):
        if self.name=='is':
            return self.file.variables['topg'].dimensions[-1]
        elif self.name=='isobase':
            return self.file.variables['slc'].dimensions[-1]
        elif self.name == 'pmp':
            return self.file.variables['thk'].dimensions[-1]
        elif self.name == 'vel':
            return self.file.variables['uvel'].dimensions[-1]
        elif self.name == 'bvel':
            return self.file.variables['ubas'].dimensions[-1]
        elif self.name == 'bvel_tavg':
            return self.file.variables['ubas_tavg'].dimensions[-1]
        elif self.name == 'tau':
            return self.file.variables['taux'].dimensions[-1]
        else:
            return self.file.variables[self.name].dimensions[-1]
    xdimension = property(__get_xdimension)

    def __get_xdim(self):
        return self.file.variables[self.xdimension]
    xdim = property(__get_xdim)

    def __get_ydim(self):
        if self.name=='is':
            return self.file.variables[self.file.variables['topg'].dimensions[-2]]
        elif self.name=='isobase':
            return self.file.variables[self.file.variables['slc'].dimensions[-2]]
        elif self.name == 'pmp':
            return self.file.variables[self.file.variables['thk'].dimensions[-2]]
        elif self.name == 'vel':
            return self.file.variables[self.file.variables['uvel'].dimensions[-2]]
        elif self.name == 'bvel':
            return self.file.variables[self.file.variables['ubas'].dimensions[-2]]
        elif self.name == 'bvel_tavg':
            return self.file.variables[self.file.variables['ubas_tavg'].dimensions[-2]]
        elif self.name == 'tau':
            return self.file.variables[self.file.variables['taux'].dimensions[-2]]
        else:
            return self.file.variables[self.file.variables[self.name].dimensions[-2]]
    ydim = property(__get_ydim)

    def __is3d(self):
        is3d = False
        if self.name not in ['is', 'isobase', 'bvel', 'bvel_tavg','pmp', 'tau']:
            if self.name == 'vel':
                is3d = True
            elif 'level' in self.file.variables[self.name].dimensions :
                is3d = True
            elif 'lithoz' in self.file.variables[self.name].dimensions :
                is3d = True
        return is3d
    is3d = property(__is3d)

    def __get_var(self):
        if self.name=='is':
            if self.__varcache == None:
                self.__varcache = self.file.variables['topg'][:,:,:]+self.file.variables['thk'][:,:,:]
            return self.__varcache

        elif self.name == 'pmp':
            if self.__varcache == None:
                ih = self.file.variables['thk'][time,:,:]
                self.__varcache = calc_pmp(ih)
            return self.__varcache

        elif self.name == 'vel':
            if self.__varcache == None:
                self.__varcache = numpy.sqrt(self.file.variables['uvel'][:,:,:,:]*self.file.variables['uvel'][:,:,:,:] +
                                               self.file.variables['vvel'][:,:,:,:]*self.file.variables['vvel'][:,:,:,:])
            return self.__varcache

        elif self.name == 'bvel':
            if self.__varcache == None:
                self.__varcache = numpy.sqrt(self.file.variables['ubas'][:,:,:]*self.file.variables['ubas'][:,:,:]+
                                               self.file.variables['vbas'][:,:,:]*self.file.variables['vbas'][:,:,:])
            return self.__varcache

        elif self.name == 'bvel_tavg':
            if self.__varcache == None:
                self.__varcache = numpy.sqrt(self.file.variables['ubas_tavg'][:,:,:]*self.file.variables['ubas_tavg'][:,:,:]+
                                               self.file.variables['vbas_tavg'][:,:,:]*self.file.variables['vbas_tavg'][:,:,:])
            return self.__varcache

        elif self.name == 'tau':
            if self.__varcache == None:
                self.__varcache = numpy.sqrt(self.file.variables['taux'][:,:,:]*self.file.variables['taux'][:,:,:]+
                                               self.file.variables['tauy'][:,:,:]*self.file.variables['tauy'][:,:,:])
            return self.__varcache
        
        else:
            return self.file.variables[self.name]
    var = property(__get_var)

    def __get_isvelogrid(self):
        return self.xdimension=='x0'
    isvelogrid = property(__get_isvelogrid)

    def get2Dfield(self,time,level=0,velogrid=False,clip=None):
        """Get a 2D field.

        time: time slice
        level: horizontal slice
        velogrid: set to true to interpolate onto velocity grid."""

        if self.average:
            if not self.is3d:
                raise RuntimeError, 'Variable %s is not 3D.'%self.name
            # integrate
            grid = numpy.zeros((self.xdim.shape[0],self.ydim.shape[0]),'f')

            sigma = self.file.variables['level']
            sliceup = self.__get2Dfield(time,level=-1,velogrid=velogrid,clip=clip)
            for k in range(sigma.shape[0]-2,-1,-1):
                g_slice = self.__get2Dfield(time,level=k,velogrid=velogrid,clip=clip)
                grid = grid+(sliceup+g_slice)*(sigma[k+1]-sigma[k])
                sliceup = self.__get2Dfield(time,level=k,velogrid=velogrid,clip=clip)
            grid = 0.5*grid
        else:
            grid = self.__get2Dfield(time,level=level,velogrid=velogrid,clip=clip)

        return grid
    
    def __get2Dfield(self,time,level=0,velogrid=False,clip=None):
        """Get a 2D field.

        time: time slice
        level: horizontal slice
        velogrid: set to true to interpolate onto velocity grid."""

        if self.is3d:
            if self.name == 'vel':
                grid = numpy.sqrt(
                    self.file.variables['uvel'][time,level,:,:]*self.file.variables['uvel'][time,level,:,:]+
                    self.file.variables['vvel'][time,level,:,:]*self.file.variables['vvel'][time,level,:,:])
            else:
                grid = self.file.variables[self.name][time,level,:,:]
        else:
            if self.name == 'is':
                grid = self.file.variables['topg'][time,:,:] + self.file.variables['thk'][time,:,:]
            elif self.name=='isobase':
                grid = self.file.variables['slc'][time,:,:]
            elif self.name == 'pmp':
                ih = self.file.variables['thk'][time,:,:]
                grid = calc_pmp(ih)
            elif self.name == 'bvel':
                grid = numpy.sqrt(
                    self.file.variables['ubas'][time,:,:]*self.file.variables['ubas'][time,:,:]+
                    self.file.variables['vbas'][time,:,:]*self.file.variables['vbas'][time,:,:])
            elif self.name == 'bvel_tavg':
                grid = numpy.sqrt(
                    self.file.variables['ubas_tavg'][time,:,:]*self.file.variables['ubas_tavg'][time,:,:]+
                    self.file.variables['vbas_tavg'][time,:,:]*self.file.variables['vbas_tavg'][time,:,:])
            elif self.name == 'tau':
                grid = numpy.sqrt(
                    self.file.variables['taux'][time,:,:]*self.file.variables['taux'][time,:,:]+
                    self.file.variables['tauy'][time,:,:]*self.file.variables['tauy'][time,:,:])
            else:
                grid = self.file.variables[self.name][time,:,:]
        if self.name in ['topg','is']:
            if 'eus' in self.file.variables.keys():
                grid = grid - self.file.variables['eus'][time]
        if self.name=='isobase':
            if 'eus' in self.file.variables.keys():
                grid = grid + self.file.variables['eus'][time]
        # correct temperature
        if self.name in temperatures:
            if self.pmt:
                if 'thk' not in self.file.variables.keys():
                    print 'Warning, cannot correct for pmt because ice thicknesses are not in file'
                else:
                    ih = self.file.variables['thk'][time,:,:]
                    if self.name == 'btemp':
                        fact = 1.
                    else:
                        fact = self.file.variables['level'][level]
                    grid = grid - calc_pmp(ih,fact)

        if velogrid:
            if not self.isvelogrid:
                grid = 0.25*(grid[:-1,:-1]+grid[1:,1:]+grid[:-1,1:]+grid[1:,:-1])
        if clip!=None:
            m = GCvariable(self.cffile,clip).get2Dfield(time,level=level)
            if grid.shape[0]!=m.shape[0]:
                m = 0.25*(m[:-1,:-1]+m[1:,1:]+m[:-1,1:]+m[1:,:-1])
            maskArray = numpy.where(m>0.,False,True)
            
            grid = numpy.ma.array(grid,mask=maskArray)            
        return grid

    def interpolate(self,x,y,time,level=0):
        """interpolate 2D field
        x: list of x coordinates
        y: list of y coordinates
        time: the time slice
        level : the vertical slice
        kind: ['linear', 'cubic', 'quintic'] - the kind of interpolation to use"""
        
        data = self.get2Dfield(time,level=level)
        x0 = self.xdim[0]
        y0 = self.ydim[0]
        dx = self.xdim[1]-x0
        dy = self.ydim[1]-y0
        ivals = (numpy.array(x)-x0)/dx
        jvals = (numpy.array(y)-y0)/dy
        coords = numpy.array([ivals, jvals])
        p =  scipy.ndimage.map_coordinates(data,coords)
        return p

    def getSpotIJ(self,node,time=None,level=0):
        """Get data at a grid node.

        node: list/tuple/array of size 2 selecting node
        time: if None, return data for all time slices
              if list/etc of size two, interpret as array selection
              if single value, get only this time slice
        level: if None get data for all levels (time must be a single value)
               otherwise get a specific level"""

        if node[0] < 0 or node[0] >= self.xdim.shape[0] or node[1] < 0 or node[1] >= self.ydim.shape[0]:
            raise RuntimeError, 'node is outside bounds'

        (tarray,t) = GCchecklist(time,self.file.variables['time'])
        (larray,l) = GCchecklist(level,self.file.variables['level'])

        if 'level' in self.file.variables[self.name].dimensions:
            (larray,l) = GCchecklist(level,self.file.variables['level'])
        elif 'lithoz' in self.file.variables[self.name].dimensions:
            (larray,l) = GCchecklist(level,self.file.variables['lithoz'])
        else:
            larray = False
            l = 0

        if larray and tarray:
            raise RuntimeError, 'Cannot select both multiple times and vertical slices'

        values = []
        if tarray:
            for i in range(t[0],t[1]+1):
                values.append(self.get2Dfield(i,l)[node[0],node[1]])
            return values
        if larray:
            for i in range(l[0],l[1]+1):
                values.append(self.get2Dfield(t,i)[node[0],node[1]])
            return values

        return self.get2Dfield(t,l)[node[0],node[1]]
        
def calc_pmp(ice_thickness, sigma = 1.):
    """Calculate pressure melting point of ice.

    ice_thickness: thickness of ice
    sigma: sigma level at which to calculate pmp"""

    return -8.7e-4*ice_thickness*sigma

if __name__ == '__main__':
    
    import sys

    cffile = GCloadfile(sys.argv[1])

    print cffile.title
    print cffile.institution
    print cffile.source
    print cffile.references
    print cffile.comment
    print cffile.history
    print cffile.ll_xy,cffile.ur_xy
    print cffile.projection.proj4_params()

    cffile.close()
