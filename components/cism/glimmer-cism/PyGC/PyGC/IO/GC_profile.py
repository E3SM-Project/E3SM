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

"""Loading GC profiles files."""

__all__=['GCprofile','GCloadprofile','GCprofvar']

import numpy
from GC_loadfile import *
from GC_readlines import GCreadlines

class GCprofile(object):
    """Handling Profile lines."""

    def __init__(self,cffile,interval=10000.,xrange=[None, None],xscale=0.001):
        """Initialise.
        
        cffile: CF file object
        interval: linearly interpolate profile
        xrange: clip profile (Default [None,None])
        xscale: scale x values
        """

        self.cffile = cffile
        self.interval = interval
        self.xscale = xscale
        self.xrange = xrange

    def set_control_points(self,xloc,yloc,projected=True):
        """Set control points.

        xloc/yloc: control points along profile.
        projected: set to True (Default) if xloc/yloc are in projected coord system
        """

        if projected:
            self.xloc = xloc
            self.yloc = yloc
        else:
            self.xloc = []
            self.yloc = []
            for i in range(0,len(xloc)):
                point = self.cffile.project([xloc[i],yloc[i]])
                self.xloc.append(point[0])
                self.yloc.append(point[1])
                
        self.interpolated = GCinterpolate_xy([self.xloc,self.yloc],self.interval)
        
        self.interval = self.interval*self.xscale
        if self.xrange[0] == None:
            start = 0
        else:
            start = int(self.xrange[0]/self.interval)
        if self.xrange[1] == None:
            end = -1
        else:
            end = int(self.xrange[1]/self.interval+0.9999)
        self.interpolated = self.interpolated[:,start:end]
        self.xrange = [start*self.interval, end*self.interval]
        self.xvalues = []
        for i in range(0,len(self.interpolated[1,:])):
            self.xvalues.append(self.xrange[0]+i*self.interval)

        # clip to region of file
        interpx = []
        interpy = []
        xvalues = []
        for i in range(0,len(self.interpolated[1,:])):
            if self.cffile.inside(self.interpolated[:,i]):
                interpx.append(self.interpolated[0,i])
                interpy.append(self.interpolated[1,i])
                xvalues.append(self.xvalues[i])
        self.xvalues = xvalues
        self.interpolated = numpy.array([interpx,interpy],'f')

    def coords_file(self,fname,projected=True):
        """Read control points from file.

        fname: name of file to be read
        projected: set to True (Default) if xloc/yloc are in projected coord system
        """

        xdata = []
        ydata = []
        infile = file(fname)
        for line in GCreadlines(infile):
            l = line.split()
            xdata.append(float(l[0]))
            ydata.append(float(l[1]))
        infile.close()
        self.set_control_points(xdata,ydata,projected=projected)
                                          
                                          
class GCloadprofile(GCloadfile):
    """Loading a profile line from a CF netCDF file."""

    def __init__(self,fname,pdata):
        """Initialise.

        fname: name of CF file.
        pdata: GCprofile instance
        xrange: clip profile (Default [None,None])
        xscale: scale x values
        """
        
        GCloadfile.__init__(self,fname)

        self.profiledata = pdata
        self.xscale = self.profiledata.xscale
        self.xloc = self.profiledata.xloc
        self.yloc = self.profiledata.yloc
        self.interval = self.profiledata.interval
        self.interpolated = self.profiledata.interpolated
        self.xvalues = self.profiledata.xvalues
        
        # caching profiles
        self.__profiles = {}

    def getprofile(self,var):
        """Get a profile variable from file.

        var: name of variables

        this method caches the return profile variable structure."""

        if var not in self.__profiles:
            self.__profiles[var] = GCprofvar(self,var)
        return self.__profiles[var]

    def getExtent(self,time=None,interval=1):
        """Get ice extent along profile.

        time: if None, return data for all time slices
              if list/etc of size two, interpret as array selection
              if single value, get only this time slice

        interval: extract every interval timeslice      

        this is a hack, we start looking from the end of the profile and stop when we found
        a change from no ice to ice."""

        (tarray,t) = GCchecklist(time,self.file.variables['time'])
        values = []
        if tarray:
            for i in range(t[0],t[1]+1,interval):
                values.append(self.__getExtent(i))
            return values
        return self.__getExtent(t)

    def __getExtent(self,time):
        data = self.getprofile('thk').getProfile(time)
        #right BC
        if data[-1] > 0. and data[-1]<1.e10:
            return self.xvalues[-1]
        for i in range(len(data[:])-2,0,-1):
            if data[i] > 0. and data[i]<1.e10:
                return self.xvalues[i]
        # if we are at the end we just return position of last value to close data set
        return self.xvalues[0]

class GCprofvar(GCvariable):
    """Handling CF Profiles."""

    def __init__(self,cfprofile,var):
        """Initialise.

        cfprofile: CF Profile file
        var: name of variable"""

        self.yres = 10.

        if not isinstance(cfprofile,GCloadprofile):
            raise ValueError, 'Not a profile file'

        GCvariable.__init__(self,cfprofile,var)
        self.showpmp = False

        #cache profile data
        self.__data = {}
        self.__data2d = {}
        
    def getProfile(self,time,level=0):
        """Get a profile.

        time: time slice
        level: horizontal slice."""

        #if time not in self.__data:
        #    self.__data[time] = {}
        #if level not in self.__data[time]:
        #    var = self.getGMTgrid(time,level=level)
        #    self.__data[time][level] = var.grdtrack(self.cffile.interpolated[0,:],self.cffile.interpolated[1,:])

        #return self.__data[time][level]

        return self.interpolate(self.cffile.interpolated[0,:],self.cffile.interpolated[1,:],time,level=level)


##     def getProfile2D_litho(self,time):
##         """Get a 2D profile which is not in sigma coordsystem.

##         i.e. litho temperature.

##         time: time slice

##         returns a GMT grid"""

##         if 'lithoz' not in self.var.dimensions:
##             raise ValueError, 'Not a 3D variable'


##         ymin = self.file.variables['lithoz'][-1]
##         ymax = self.file.variables['lithoz'][0]
##         yres = 50.
##         y = self.file.variables['lithoz'][:].tolist()
##         y.reverse()

##         # load data
##         data = numpy.zeros([len(self.file.variables['lithoz']), len(self.cffile.xvalues)], 'f')
##         for i in range(0,len(self.file.variables['lithoz'])):
##             prof = self.getProfile(time,level=len(self.file.variables['lithoz'])-i-1)
##             data[i,:] = prof

##         # setup output grid
##         grid = Grid()
##         grid.x_minmax = [0,self.cffile.xvalues[-1]]
##         grid.y_minmax = [ymin,ymax]
##         numy=int((ymax-ymin)/yres)+1
##         grid.data = numpy.zeros([len(self.cffile.xvalues),numy], 'f')

##         # interpolate
##         pos = (numpy.arange(numy)-numy+1)*yres
##         for j in range(0,len(self.cffile.xvalues)):
##             interpolated = CFinterpolate_linear(y,data[:,j],pos)
##             grid.data[j,:] = interpolated

##         return grid

    def getProfile2D(self,time):
        """Get a 2D profile.

        time: time slice

        returns a GMT grid"""

        if 'level' not in self.var.dimensions:
            raise ValueError, 'Not a 3D variable'

        if time not in self.__data2d:
            # load ice thickness and bedrock profiles
            ihprof = numpy.array(GCprofvar(self.cffile,'thk').getProfile(time))
            try:
                rhprof = numpy.array(GCprofvar(self.cffile,'topg').getProfile(time))
            except:
                rhprof = numpy.zeros(len(ihprof))
            ymin=min(rhprof)
            ymax=max(rhprof+ihprof)
            numy=int((ymax-ymin)/self.yres)+1

            # load data
            data = numpy.zeros([len(self.file.variables['level']), len(self.cffile.xvalues)], 'f')
            for i in range(0,len(self.file.variables['level'])):
                prof = self.getProfile(time,level=len(self.file.variables['level'])-i-1)
                data[i,:] = prof

            grid = numpy.zeros([len(self.cffile.xvalues),numy], 'f')
            mask = numpy.ones([len(self.cffile.xvalues),numy], bool)
            # interpolate
            rhprof = rhprof-ymin
            for j in range(0,len(self.cffile.xvalues)):
                if ihprof[j]>0.:
                    start = int(rhprof[j]/self.yres)
                    end   = int((rhprof[j]+ihprof[j])/self.yres)+1
                    pos = numpy.arange(start*self.yres,end*self.yres,self.yres)
                    grid[j,start:end] = GCinterpolate_linear(rhprof[j]+ihprof[j]*self.file.variables['level'][:],
                                                             data[:,j],
                                                             pos)
                    mask[j,start:end] = False

            #self.__data2d[time] = grid
            return [ymin,ymax],numpy.transpose(numpy.ma.array(grid,mask=mask))
        #return self.__data2d[time]

    def getProfileTS(self,time=None,level=0):
        """Get a time-distance data.

        
        time: if None, return data for all time slices
              if list/etc of size two, interpret as array selection
              if single value, get only this time slice
        level: horizontal slice."""

        (tarray,t) = GCchecklist(time,self.file.variables['time'])

        if tarray:
            data = []
            for i in range(t[0],t[1]+1):
                data.append(self.getProfile(i,level=level))
            data = numpy.transpose(numpy.array(data))
            return data
        else:
            return self.getProfile(t,level=level)

def GCinterpolate_xy(profile,interval):
    "linearly interpolate profile, return interpolated array and number of points outside region"

    data = []
    p = numpy.array(profile)
    for i in range(0,len(p[0,:])):
        data.append([float(p[0,i]),float(p[1,i])])
    remainder = lr = 0.
    d0 = data[0]
    ix = []
    iy = []
    ix.append(d0[0])
    iy.append(d0[1])

    for d in data[1:]:
        x = d[0] - d0[0]
        y = d[1] - d0[1]
        dist = numpy.sqrt(x*x + y*y)
        remainder = remainder + dist
        cm = x/dist
        sm = y/dist
        while (remainder - interval) >= 0.:
            d0[0] = d0[0] + (interval-lr)*cm
            d0[1] = d0[1] + (interval-lr)*sm
            lr = 0.
            ix.append(d0[0])
            iy.append(d0[1])
            remainder = remainder - interval
        lr = remainder
        d0 = d    
    
    return numpy.array([ix,iy],'f')

def GCinterpolate_linear(x,y,pos):
    """Linear interpolation.

    x,y: data points
    pos: list of new x positions onto which y should be interpolated.

    we assume monotonic sequences in x and pos"""

    if len(x)!=len(y):
        raise ValueError, 'x and y are not of the same length.'

    res = []
    j = 0
    for i in range(0,len(pos)):
        # boundary conditions
        if pos[i] <= x[0]:
            res.append(y[0])
            continue
        if pos[i] >= x[-1]:
            res.append(y[-1])
            continue
        while (pos[i]>x[j]):
            j = j + 1
        res.append(y[j-1]+(pos[i]-x[j-1])*(y[j]-y[j-1])/(x[j]-x[j-1]))
    return res
