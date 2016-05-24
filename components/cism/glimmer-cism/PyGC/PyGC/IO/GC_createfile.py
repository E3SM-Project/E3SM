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

"""Creating GC files."""

__all__=['GCVariableDef','GCcreatefile']

import numpy,ConfigParser,os,re,string, glob
from GC_netcdf import GCNetCDFFile
from GC_file import *

NOATTRIB = ['name','dimensions','data','factor','load','f90file','hot','type','dimlen']

class GCVariableDef(dict):
    """Dictionary containing variable definitions."""

    def __init__(self,filename):
        """Initialise Variable class.

        filename: name or list of names of file(s) containing variable definitions."""

        dict.__init__(self)

        # reading variable configuration file
        vars = ConfigParser.ConfigParser()
        vars.read(filename)

        for v in vars.sections():
            vardef = {}
            vardef['name'] = v
            for (name, value) in vars.items(v):
                vardef[name] = value
            self.__setitem__(v,vardef)

    def keys(self):
        """Reorder standard keys alphabetically."""
        dk = []
        vk = []
        for v in dict.keys(self):
            if is_dimvar(self.__getitem__(v)):
                dk.append(v)
            else:
                vk.append(v)
        dk.sort()
        vk.sort()
        return dk+vk

class GCcreatefile(GCfile):
    """Creating a GC netCDF file."""

    def __init__(self,fname,append=False):
        """Initialise.

        fname: name of GC file.
        append: set to true if file should be open rw"""

        GCfile.__init__(self,fname)
        self.mapvarname = 'mapping'
        # get variable definitions
        try:
            vname=os.environ['GLIMMER_PREFIX']
        except KeyError:
            vname = os.path.expanduser(os.path.join('~','glimmer'))
        vname = os.path.join(vname,'share','glimmer')
        if not os.path.exists(vname):
            raise RuntimeError, 'Cannot find ncdf_vars.def\nPlease set GLIMMER_HOME to where glimmer is installed'
        self.vars = GCVariableDef(glob.glob(vname+'/*.def'))

        if append:
            self.file = GCNetCDFFile(self.fname,'a')
        else:
            self.file = GCNetCDFFile(self.fname,'w')
        self.file.Conventions = "GC-1.0"
        
    def createDimension(self,name, length):
        """Create a dimension.

        Creates a new dimension with the given name and length.
        length must be a positive integer or None, which stands for
        the unlimited dimension. Note that there can be only one
        unlimited dimension in a file."""
        self.file.createDimension(name,length)

    def createVariable(self,name):
        """Create a GC variable.

        name: name of variable."""

        if name not in self.vars:
            raise KeyError, 'Cannot find definition for variable %s'%name
        v = self.vars[name]
        var = self.file.createVariable(name,'f',tuple(string.replace(v['dimensions'],' ','').split(',')))
        for a in v:
            if a not in NOATTRIB:
                setattr(var,a,v[a])
        if self.mapvarname != '' and 'x' in v['dimensions'] and 'y' in v['dimensions']:
            var.grid_mapping = self.mapvarname
        return var

if __name__ == '__main__':
    # creating a test netCDF file

    import GC_proj

    filename="test.nc"
    numx=100
    numy=150

    proj = GC_proj.DummyProj()
    proj.grid_mapping_name='albers_conical_equal_area'
    proj.false_easting = [1903971.]
    proj.false_northing = [898179.3]
    proj.longitude_of_central_meridian = [33.5]
    proj.latitude_of_projection_origin = [60.5]
    proj.standard_parallel = [52.83333, 68.16666]

    cffile = GCcreatefile(filename)
    cffile.title = "Test GC file"
    cffile.institution = "University of Edinburgh"
    cffile.source = "None"
    cffile.comment = "Testing if our netCDF files conform with GC standard"

    # creating dimensions
    cffile.createDimension('x0',numx-1)
    cffile.createDimension('x1',numx)
    cffile.createDimension('y0',numy-1)
    cffile.createDimension('y1',numy)
    cffile.createDimension('level',1)
    cffile.createDimension('staglevel',1)
    cffile.createDimension('lithoz',1)
    cffile.createDimension('time',None)

    cffile.projection=proj
    
    # creating variables
    var=cffile.createVariable('x0')
    var[:]=numpy.arange(numx-1).astype('f')
    var=cffile.createVariable('x1')
    var[:]=numpy.arange(numx).astype('f')
    var=cffile.createVariable('y0')
    var[:]=numpy.arange(numy-1).astype('f')
    var=cffile.createVariable('y1')
    var[:]=numpy.arange(numy).astype('f')

   
    for v in cffile.vars:
         if 'spot' not in v and v not in ['VARSET','level','staglevel','x0','y0','x1','y1']:
             print 'Creating variable %s'%v
             var = cffile.createVariable(v)

    cffile.close()
