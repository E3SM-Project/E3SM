# Copyright (C) 2004, 2005, 2006, 2009, 2010
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

"""Handling GC grid mappings"""

__all__=['DummyProj','GCProj','GCProj_aea','GCProj_lcc','getGCProj','copyGCMap',
         'GCProj_parse_GMTproj','GCProj_parse_ESRIprj','GCProj_printGMThelp']

import numpy, string
from PyGC import Proj

class DummyProj:
    """Emptiy class for storing GC projection."""
    pass

class GCProj(object):
    def __init__(self,var):
        """Initialise.

        var: GC grid mapping variable."""
        self.params = {}
        self.params['proj'] = None
        self.params['ellps'] = 'WGS84' # set default ellipsoid
        try:
            self.params['x_0'] = numpy.squeeze(var.false_easting)
        except:
            self.params['x_0'] = 0.
        try:
            self.params['y_0'] = numpy.squeeze(var.false_northing)
        except:
            self.params['y_0'] = 0.
        self.gmt_type = ''
        self.__Proj4 = None 

    def __get_Proj4(self):
        if self.__Proj4 == None:
            self.__Proj4 = Proj(self.proj4_params())
        return self.__Proj4
    Proj4 = property(__get_Proj4)

    def proj4_params(self):
        params = []
        for p in self.params.keys():
            params.append('%s=%s'%(p,self.params[p]))
        return params

    def proj4(self,point,inv=False):
        """Do projection.

        point: 2D point
        inv:   if True do the inverse projection."""

        if inv:
            return self.Proj4.inv(point)
        else:
            return self.Proj4.fwd(point)

    def setOrigin(self,lon0,lat0):
        """Set origin of projected grid.

        lon0: Logitude of origin
        lat0: Latitude of origin."""

        orig = self.proj4([lon0,lat0])
        self.params['x_0'] = -orig[0]
        self.params['y_0'] = -orig[1]
        self.__Proj4 = Proj(self.proj4_params())

    def getGMTregion(self,ll,ur):
        """Get GMT region string.

        ll: coordinates of lower left corner in projected grid
        ur: coordinates of upper right corner in projected grid."""

        ws = self.Proj4.inv(ll)
        en = self.Proj4.inv(ur)
        return '%f/%f/%f/%fr'%(ws[0],ws[1],en[0],en[1])

class GCProj_stere(GCProj):
    """Stereographic Projections."""

    def __init__(self,var):
        """Initialise.
        
        var: GC grid mapping variable."""
        GCProj.__init__(self,var)
        self.params['proj'] = 'stere'
        # polar variation
        if var.grid_mapping_name == 'polar_stereographic':
            self.params['lon_0'] = numpy.squeeze(var.straight_vertical_longitude_from_pole)
            self.params['lat_0'] = numpy.squeeze(var.latitude_of_projection_origin)
            if hasattr(var,'standard_parallel'):
                self.params['lat_ts'] = numpy.squeeze(var.standard_parallel)
            elif hasattr(var,'scale_factor_at_projection_origin'):
                self.params['k_0'] = numpy.squeeze(var.scale_factor_at_projection_origin)
        else: # others
            self.params['lat_0'] = numpy.squeeze(var.latitude_of_projection_origin)
            self.params['lon_0'] = numpy.squeeze(var.longitude_of_projection_origin)
            self.params['k_0'] = numpy.squeeze(var.scale_factor_at_projection_origin)
        self.gmt_type = 's'

    def getGMTprojection(self,mapwidth=None):
        """Get GMT projection string."""

        if mapwidth==None:
            if 'k_0' in self.params:
                gmt = '%s%f/%f/%f'%(self.gmt_type,
                                    self.params['lon_0'],self.params['lat_0'],self.params['k_0'])
            else:
                gmt = '%s%f/%f/%f'%(self.gmt_type,
                                    self.params['lon_0'],self.params['lat_0'],self.params['lat_ts'])
        else:
            gmt = '%s%f/%f'%(self.gmt_type,
                             self.params['lon_0'],self.params['lat_0'])
        return gmt
        
class GCProj_laea(GCProj):
    """Lambert Azimuthal Equal Area"""
    def __init__(self,var):
        """Initialise.
        
        var: GC grid mapping variable."""
        GCProj.__init__(self,var)
        self.params['proj'] = 'laea'
        self.params['lat_0'] = numpy.squeeze(var.latitude_of_projection_origin)
        self.params['lon_0'] = numpy.squeeze(var.longitude_of_central_meridian)
        self.gmt_type = 'a'

    def getGMTprojection(self,mapwidth=None):
        """Get GMT projection string."""

        gmt = '%s%f/%f'%(self.gmt_type,
                         self.params['lon_0'],self.params['lat_0'])
        return gmt
    
class GCProj_aea(GCProj_laea):
    """Albers Equal-Area Conic."""

    def __init__(self,var):
        """Initialise.
        
        var: GC grid mapping variable."""
        GCProj_laea.__init__(self,var)
        if len(var.standard_parallel) == 2:
            self.params['lat_1'] = var.standard_parallel[0]
            self.params['lat_2'] = var.standard_parallel[1]
        elif len(var.standard_parallel) < 1 or len(var.standard_parallel) > 2:
            raise RuntimeError, 'Wrong size of standard_parallel attribute'
        else:
            self.params['lat_1'] = numpy.squeeze(var.standard_parallel)
        self.params['proj'] = 'aea'
        self.gmt_type = 'b'

    def getGMTprojection(self,mapwidth=None):
        """Get GMT projection string."""

        gmt = '%s%f/%f/%f/%f'%(self.gmt_type,
                               self.params['lon_0'],self.params['lat_0'],
                               self.params['lat_1'],self.params['lat_2'])
        return gmt

class GCProj_lcc(GCProj_aea):
    """Lambert Conic Conformal."""

    def __init__(self,var):
        """Initialise.
        
        var: GC grid mapping variable."""
        GCProj_aea.__init__(self,var)
        self.params['proj'] = 'lcc'
        self.gmt_type = 'l'

def getGCProj(var):
    """Get projection from GC grid mapping variable.
    
    var: GC grid mapping variable."""
    
    GCProj_MAP = {'lambert_azimuthal_equal_area' : GCProj_laea,
                  'albers_conical_equal_area' : GCProj_aea,
                  'lambert_conformal_conic' : GCProj_lcc,
                  'polar_stereographic' : GCProj_stere,
                  'stereographic' : GCProj_stere}

    if var.grid_mapping_name not in GCProj_MAP:
        raise KeyError, 'Error, no idea how to handle projection: %s'%var.grid_mapping_name

    return GCProj_MAP[var.grid_mapping_name](var)

def copyGCMap(orig,copy):
    """Copy GC mapping variable."""

    if 'grid_mapping_name' in dir(orig):
        copy.grid_mapping_name = orig.grid_mapping_name
    if 'standard_parallel' in dir(orig):
        copy.standard_parallel = numpy.array(orig.standard_parallel)
    if 'longitude_of_central_meridian' in dir(orig):
        copy.longitude_of_central_meridian = numpy.array(orig.longitude_of_central_meridian)
    if 'latitude_of_projection_origin' in dir(orig):
        copy.latitude_of_projection_origin = numpy.array(orig.latitude_of_projection_origin)
    if 'longitude_of_projection_origin' in dir(orig):
        copy.longitude_of_projection_origin = numpy.array(orig.longitude_of_projection_origin)
    if 'false_easting' in dir(orig):
        copy.false_easting = numpy.array(orig.false_easting)
    if 'false_northing' in dir(orig):
        copy.false_northing = numpy.array(orig.false_northing)
    if 'straight_vertical_longitude_from_pole' in dir(orig):
        copy.straight_vertical_longitude_from_pole = numpy.array(orig.straight_vertical_longitude_from_pole)
    if 'scale_factor_at_projection_origin' in dir(orig):
        copy.scale_factor_at_projection_origin = numpy.array(orig.scale_factor_at_projection_origin)

def GCProj_parse_GMTproj(projstring):
    """Parse GMT projection string."""

    proj = DummyProj()
    
    ps = projstring[1:].split('/')
    if projstring[0] in ['b','B']:
        if len(ps) != 4:
            print 'Error, wrong number of projection arguments'
            usage()
            sys.exit(1)
        proj.grid_mapping_name='albers_conical_equal_area'
        proj.longitude_of_central_meridian = float(ps[0])
        proj.latitude_of_projection_origin = float(ps[1])
        proj.standard_parallel = [float(ps[2]),float(ps[3])]
    elif projstring[0] in ['l','L']:
        if len(ps) != 4:
            print 'Error, wrong number of projection arguments'
            usage()
            sys.exit(1)
        proj.grid_mapping_name='lambert_conformal_conic'
        proj.longitude_of_central_meridian = float(ps[0])
        proj.latitude_of_projection_origin = float(ps[1])
        proj.standard_parallel = [float(ps[2]),float(ps[3])]
    elif projstring[0] in ['a','A']:
        if len(ps) != 2:
            print 'Error, wrong number of projection arguments'
            usage()
            sys.exit(1)
        proj.grid_mapping_name='lambert_azimuthal_equal_area'
        proj.longitude_of_central_meridian = float(ps[0])
        proj.latitude_of_projection_origin = float(ps[1])
    elif projstring[0] in ['s','S']:
        if len(ps) == 2 or len(ps) == 3:
            proj.latitude_of_projection_origin = float(ps[1])
            if len(ps) == 3:
                proj.grid_mapping_name='polar_stereographic'
                proj.straight_vertical_longitude_from_pole = float(ps[0])
                proj.standard_parallel = float(ps[2])
            else:
                proj.grid_mapping_name='stereographic'
                proj.longitude_of_projection_origin = float(ps[0])
                proj.scale_factor_at_projection_origin = 1.
        else:
            print 'Error, wrong number of projection arguments'
            usage()
            sys.exit(1)
    else:
        print 'Error, no idea about projection: ',projstring
        usage()
        sys.exit(1)
    return proj

def GCProj_parse_ESRIprj(pFile):
    """Parse an ESRI projection file.

    pFile: file object."""

    proj = DummyProj()

    pName = pFile.readline()
    pName=pName[14:]
  
    pStrip = string.strip(pName)
      
    if pStrip == 'ALBERS':
        proj.grid_mapping_name='albers_conical_equal_area'
        line = pFile.readline()
        while line != '':
            line = pFile.readline()
            if string.strip(line) == 'Parameters':
                standard1=DmsParse(pFile)
                standard2=DmsParse(pFile)
                proj.standard_parallel = [standard1,standard2]
                long_central=DmsParse(pFile)
                if long_central < 0.:
                    long_central = 360.+long_central
                    proj.longitude_of_central_meridian = long_central
                else:
                    proj.longitude_of_central_meridian = long_central
                lat_proj_origin=DmsParse(pFile)
                proj.latitude_of_projection_origin = lat_proj_origin
                proj.false_easting=0
                proj.false_northing=0

  
    elif pStrip == 'LAMBERT':
        proj.grid_mapping_name='lambert_conformal_conic'
        line = pFile.readline()
        while line != '':
            line = pFile.readline()
            if string.strip(line) == 'Parameters':
                standard1=DmsParse(pFile)
                standard2=DmsParse(pFile)
                proj.standard_parallel = [standard1,standard2]
                long_central=DmsParse(pFile)
                if long_central < 0.:
                    long_central = 360.+long_central
                    proj.longitude_of_central_meridian = long_central
                else:
                    proj.longitude_of_central_meridian = long_central
                lat_proj_origin=DmsParse(pFile)
                proj.latitude_of_projection_origin = lat_proj_origin
                proj.false_easting=0
                proj.false_northing=0
                  
    elif pStrip == 'LAMBERT_AZIMUTHAL':
        proj.grid_mapping_name='lambert_azimuthal_equal_area'
        line = pFile.readline()
        while line != '':
            line = pFile.readline()
            if line != '':
          	first = string.split(line)
          	if first[0] == 'Parameters':
                    skip = pFile.readline() #skip a line as data not required
                    long_central=DmsParse(pFile)
                    if long_central < 0.:
                        long_central = 360.+long_central
                        proj.longitude_of_central_meridian = long_central
                    else:
                        proj.longitude_of_central_meridian = long_central
                    lat_proj_origin=DmsParse(pFile)
                    proj.latitude_of_projection_origin = lat_proj_origin
                    proj.false_easting=0
                    proj.false_northing=0



    # Due to inconsistencies between ESRI projection files and the parameters required by GCProj, the stereographic
    # projections may produce incorrect results so CHECK them carefully!. 
    elif pStrip == 'STEREOGRAPHIC':
      
        line = pFile.readline()
        while line != '':
            line = pFile.readline()
            if line != '':
                first = string.split(line)
                if first[0] == 'Parameters':

                    l = string.split(pFile.readline())
                    ptype = l[0]
                    if ptype == '2':
                        flag = 0
                        if flag == 0:
                            long_central=ProjParser.DmsParse(pFile)
                            lat_proj_origin=DmsParse(pFile)
                            proj.latitude_of_projection_origin = lat_proj_origin
                            pol_v_eq = string.split(pFile.readline()) # get whether a polar or equatorial view.
                            pol_v_eq = pol_v_eq[0]
                            proj.false_easting=0
                            proj.false_northing=0
                            flag = 1
                  	
                        if flag == 1:
                            if pol_v_eq == 'EQUATORIAL':
                      		proj.grid_mapping_name='stereographic' # type 2 stereographic projection with equatorial
                      		proj.longitude_of_central_meridian = long_central # long of cent projection
                      		l = string.split(pFile.readline())
                      		scale_factor=l[0]
                      		proj.scale_factor_at_projection_origin = scale_factor # scale factor
                            else:
                      		proj.grid_mapping_name='polar_stereographic' # type 2 stereographic projection with north or southpole
                      		proj.straight_vertical_longitude_from_pole = long_central # geotiff=central meridian arc=long of cent meridian. Not sure that this is the right parameter.
                      		#l = string.split(pFile.readline())
                      		#lat_std_par = int(float(l[0]))
                      		lat_std_par = DmsParse(pFile)
                      		proj.standard_parallel = lat_std_par # lat of standard parallel
                      		
              
                    elif ptype == '1':
                        proj.grid_mapping_name='stereographic'
                        line=pFile.readline()
                        long_central=DmsParse(pFile)
                        proj.longitude_of_projection_origin = long_central
                        lat_proj_origin=DmsParse(pFile)
                        proj.latitude_of_projection_origin = lat_proj_origin
                        scale_factor=1. #defaults to 1 as Arc projection files do not include scale factors
                        proj.scale_factor_at_projection_origin = scale_factor
                        #print proj.scale_factor_at_projection_origin
                        false_easting=string.split(string.strip(pFile.readline()))
                        false_easting=false_easting[0]
                        proj.false_easting=0
                        false_northing=string.split(string.strip(pFile.readline()))
                        false_northing=false_northing[0]
                        proj.false_northing=0
    
    else:
        print 'Projection ',pStrip,' not recognized'
        print 'This program will only recognise albers equal area, lambert, lambert azimuthal and stereographic projections.'
        sys.exit(1)    

    return proj

def DmsParse(pFile): 
    # This parses a string from the projection file
    l = string.split(pFile.readline())
    deg = float(l[0])
    minute = float(l[1])/60.
    sec = float(l[2])/3600.
    dms=(deg)+(minute)+(sec)
    return dms

def GCProj_printGMThelp():
    print '  -Jspec\n\tGMT like projection specification'
    print '\t  -Jalon0/lat0\n\t      Lambert Azimuthal Equal Area. Give projection center'
    print '\t  -Jblon0/lat0/lat1/lat2\n\t      Albers Equal-Area Conic. Give projection center and\n\t      two  standard  parallels.'
    print '\t  -Jllon0/lat0/lat1/lat2\n\t      Lambert Conic Conformal. Give projection center and\n\t      two  standard  parallels.'
    print '\t  -Jslon0/lat0[/slat]\n\t      Stereographic projection. Give projection center and\n\t      optionally standard parallel.'
