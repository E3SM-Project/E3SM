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

"""Handle command line options in a standard way."""

__all__ = ['GCOptParser','GCOptions']

import optparse, sys, os.path
from IO.GC_loadfile import *
from GC_colourmap import *
from IO.GC_profile import *

## from CF_rsl import CFRSLlocs
## from CF_IOmisc import CFreadlines

class GCOptParser(optparse.OptionParser):
    """Handle options."""

    def __init__(self,usage = "usage: %prog [options] infile"):
        """Initialise.

        usage: usage string.
        """
        optparse.OptionParser.__init__(self,usage)

        self.rsldb = None
        
    def plot(self):
        """Plot options."""

        group = optparse.OptionGroup(self,"Plot Options","These options are used to control the appearance of the plot")
        group.add_option("--size",type="float",dest="size",nargs=2,metavar="WIDTH HEIGHT", help="size of plot in cm")
        #group.add_option("--verbose",action="store_true", dest="verbose",default=False,help="Be verbose")
        group.add_option("-o","--output",metavar="NAME",help="write image to FILE")
        group.add_option("--title",action="store_true",default=False,help="display a title")
        group.add_option("--timestamp",action="store_true",default=False,help="display a timestamp")
        self.add_option_group(group)

    def region(self):
        """Specifying region of interest."""

        group = optparse.OptionGroup(self,"Region Options","These options are used to control the region of interest.")
        group.add_option("--llx",dest='llx',metavar="X Y",type="float",nargs=2,help="lower left corner in projected coordinate system")
        group.add_option("--urx",dest='urx',metavar="X Y",type="float",nargs=2,help="upper right corner in projected coordinate system")
        group.add_option("--llg",dest='llg',metavar="X Y",type="float",nargs=2,help="lower left corner in geographic coordinate system")
        group.add_option("--urg",dest='urg',metavar="X Y",type="float",nargs=2,help="upper right corner in geographic coordinate system")
        self.add_option_group(group)

    def region1d(self,onlyx=False,onlyy=False):
        """Specifying axis ranges."""

        group = optparse.OptionGroup(self,"Axis Options","These options are used to control the x and y axis.")
        if not onlyy:
            group.add_option("--xrange",type="float",nargs=2,metavar="X1 X2",help="set x-axis range to X1:X2")
        if not onlyx:
            group.add_option("--yrange",type="float",nargs=2,metavar="Y1 Y2",help="set y-axis range to Y1:Y2")
        self.add_option_group(group)


    def eisforcing(self):
        """Options for handling EIS forcing time series."""

        group = optparse.OptionGroup(self,"EIS forcing","Files containing time series used for forcing EIS.")
        group.add_option("--ela",dest='elafile',metavar="FILE",type="string",help="Name of file containing ELA forcing")
        group.add_option("--temp",dest='tempfile',metavar="FILE",type="string",help="Name of file containing temperature forcing")
        group.add_option("--type_temp",type="choice",metavar="TYPE",choices=['poly','exp'],default="poly",help="Select temperature calculations (default: poly)")
        group.add_option("--lat0_temp",type="float",metavar="LAT",default=44.95,help="Origin latitude for temperature calculations using exponential type (default: 44.95)")
        group.add_option("--slc",dest='slcfile',metavar="FILE",type="string",help="Name of file containing SLC forcing")
        self.add_option_group(group)

    def __var(self):
        # variable options
        self.add_option("-v","--variable",metavar='NAME',action='append',type="string",dest='vars',help="variable to be processed (this option can be used more than once), append _avg to get the vertically integrated average")
        self.add_option("-l","--level",metavar="LEV",type='int',dest='level',help='level to be plotted')
        self.add_option("--pmt",action="store_true", dest="pmt",default=False,help='Correct temperature for temperature dependance on pressure')

    def var_options(self):
        """Extra variable stuff"""
        
        self.add_option("-c","--clip",metavar='VAR',type="choice",dest='clip',choices=['thk','topg','usurf','is'],help="display variable only where ['thk','topg','usurf','is']>0.")
        self.add_option("-i","--illuminate",metavar='VAR',type="choice",dest='illuminate',choices=['thk','topg','usurf','is'],help="illuminate surface using gradient of ['thk','topg','usurf','is']")
        self.add_option("--vectors",metavar='VAR',type="choice",choices=['vel','bvel'],help="plot velocity vectors of ['vel','bvel']")
        self.add_option("--land",action="store_true", dest="land",default=False,help="Indicate area above SL")
        self.add_option("--mask",action="store_true",default=False,help="Indicate ice thickness mask")
        self.add_option("--no-geo-coords",action="store_false",dest="geo_coord",default=True,help="do not plot geographic coordinate system")
        try:
            self.add_option("--colourmap",type="string",dest="colourmap",help="name of GMT cpt file to be used (autogenerate one when set to None)")
        except:
            pass
        try:
            self.add_option("--legend",metavar='L',type="choice",choices=['h','v'],default=None,help="Plot a colour legend, specify 'h' for a horizontal or 'v' for a vertical legend")
        except:
            pass
            
    def variable(self):
        """Variable option."""

        self.__var()
        self.var_options()
        
    def spot(self):
        """Spot options."""

        self.__var()
        self.add_option("--ij",dest='ij',metavar="I J",type="int",nargs=2,action='append',help="node to be plotted (this option can be used more than once)")

    def profile_file(self,plist=False):
        """Options for profile files.

        plist: set to True if a number of profiles can be specified"""

        if plist:
            self.add_option("-p","--profile",action="append",metavar='PROFILE',type='string',dest='profname',help="name of file containing profile control points (this option can be used more than once)")
        else:
            self.add_option("-p","--profile",metavar='PROFILE',type='string',dest='profname',help="name of file containing profile control points")
        self.add_option("--not_projected",action="store_false",default=True,dest="prof_is_projected",help="Set this flag if the profile data is not projected.")        
        self.add_option("--interval",type="float",metavar='INTERVAL',default=10000.,help="set interval to INTERVAL (default = 10000.m)")

    def profile(self,vars=True):
        """Profile options.

        vars: set to False if only profile is needed"""

        if vars:
            self.__var()
        self.profile_file()
        self.add_option("--showpmp",action="store_true", dest="showpmp",default=False,help='Indicate pressure melting point of ice (only used for temperatures)')
        try:
            self.add_option("--colourmap",type="string",dest="colourmap",help="name of GMT cpt file to be used (autogenerate one when set to None)")
        except:
            pass
        try:
            self.add_option("--legend",action="store_true", dest="dolegend",default=False,help="Plot a colour legend")
        except:
            pass
        
    def time(self):
        """Time option."""
        self.add_option("-t","--time",metavar='TIME',action='append',type="float",dest='times',help="time to be processed (this option can be used more than once)")
        self.add_option("-T","--timeslice",metavar='N',action='append',type="int",help="time slice to be processed (this option can be used more than once)")

    def timeint(self):
        """Time interval options."""
        self.add_option("-t","--time",metavar='T0 T1',type="float",nargs=2,dest='times',help="specify time interval T0 T1, if none process entire file.")

    def epoch(self):
        """Glacial Stages."""
        self.add_option("-e","--epoch",metavar='NAME',type="string",help='load glacial stages from file and plot them on time axis')

##     def rsl(self):
##         """RSL options."""
##         self.add_option("-r","--rsldb",metavar='DB',type="string",default=self.rsldb,help="name of RSL database file [%s]"%self.rsldb)
##         self.add_option("--rsl_selection",type="choice",choices=CFRSLlocs.keys(),default='fenscan',help="Change selection of RSL locations to be plotted, can be one of %s (default: fenscan)"%str(CFRSLlocs.keys()))

class GCOptions(object):
    """Do some option/argument massaging."""

    def __init__(self,parser,numargs=None):
        """Initialise.

        parser: Option parser.
        numargs: the number of arguments expected. A negative numargs implies the minimum number of arguments."""

        self.parser = parser

        (self.options, self.args) = self.parser.parse_args()

        if numargs != None:
            if numargs>=0:
                if len(self.args)!=numargs:
                    self.parser.error('Error, expected %d arguments and got %d arguments\n'%(numargs,len(self.args)))
            else:
                if len(self.args)<-numargs:
                    self.parser.error('Error, expected at least %d arguments and got %d arguments\n'%(-numargs,len(self.args)))
                    
    def __get_nfiles(self):
        return len(self.args)
    nfiles = property(__get_nfiles)
    
    def __get_nvars(self):
        try:
            return len(self.options.vars)
        except:
            return 1
    nvars = property(__get_nvars)

    def __get_ntimes(self):
        try:
            return len(self.options.times)
        except:
            return 1
    ntimes = property(__get_ntimes)

    def __get_plotsize(self):
        if self.options.size != None:
            inch = 2.54
            return [self.options.size[0]*inch,self.options.size[1]*inch]
    plotsize = property(__get_plotsize)

    def gcfile(self,argn=0):
        """Load CF file.

        argn: number of argument holding CF file name."""

        infile = GCloadfile(self.args[argn])

        if hasattr(self.options,'llx'):
            if self.options.llx != None:
                infile.ll_xy = list(self.options.llx)
        if hasattr(self.options,'urx'):
            if self.options.urx != None:
                infile.ur_xy = list(self.options.urx)
        if hasattr(self.options,'llg'):
            if self.options.llg != None:
                infile.ll_geo = list(self.options.llg)
        if hasattr(self.options,'urg'):
            if self.options.urg != None:
                infile.ur_geo = list(self.options.urg)
        
        return infile

    def cfprofile(self,argn=0):
        """Load CF profile.

        argn: number of argument holding CF file name."""

        try:
            xrange=self.options.xrange
        except:
            xrange=None
        if xrange==None:
            xrange=[None,None]
        # load profile data
        profileData = GCprofile(GCloadfile(self.args[argn]),interval=self.options.interval,xrange=xrange)
        profileData.coords_file(self.options.profname,self.options.prof_is_projected)
        profile = GCloadprofile(self.args[argn],profileData)

        return profile

    def vars(self,gcfile,varn=0):
        """Get variable.

        gcfile: CF netCDF file
        varn: variable number
        """

        try:
            var = gcfile.getvar(self.options.vars[varn])
        except KeyError:
            self.parser.error("Cannot find variable %s in file %s"%(self.options.vars[varn],gcfile.fname))
        try:
            if self.options.colourmap == 'None':
                var.colourmap = 'auto'
            elif self.options.colourmap != None:
                var.colourmap = self.options.colourmap
        except:
            var.colourmap = 'auto'
        var.pmt = self.options.pmt
        return var

    def colourmap(self,gcfile,varn=0):
        """Get colourmap.

        gcfile: CF netCDF file
        varn: variable number
        """

        try:
            var = gcfile.getvar(self.options.vars[varn])
        except KeyError:
            self.parser.error("Cannot find variable %s in file %s"%(self.options.vars[varn],gcfile.fname))

        cmap = GCcolourmap(var)

        if self.options.colourmap == 'None':
            return (None,None,cmap.title)
        else:
            return (cmap.norm,cmap.colourmap,cmap.title)

    def profs(self,cffile,varn=0):
        """Get profiles.

        cffile: CF netCDF profile file
        varn: variable number
        """

        prof = cffile.getprofile(self.options.vars[varn])
        prof.pmt = self.options.pmt
        try:
            prof.showpmp = self.options.showpmp
        except:
            pass
        prof.colourmap = self.colourmap(cffile,varn=varn)
        return prof

    def times(self,gcfile,timen=0):
        """Get time slice.

        timen: time number."""

        if self.options.times != None:
            return gcfile.timeslice(self.options.times[timen])
        elif self.options.timeslice !=None:
            return self.options.timeslice[timen]
        else:
            return 0
