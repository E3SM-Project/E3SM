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

"""Handling colourmaps."""

__all__=['GCcolourmap']

import os,numpy
import colorsys
import sys,os.path

GC_sharedir = os.path.join(os.path.dirname(sys.argv[0]),"..","data")
try:
    from GC_autoconf import *
except:
    pass

import matplotlib.colors


NC_FILL_VALUE = 1.9938419936773738e+37

class GCcolourmap(object):
    """Colourmaps."""

    STDN_MAP = { 'bedrock_altitude' : 'topo.cpt',
                 'land_ice_thickness' : 'ice.cpt',
                 'land_ice_surface_mass_balance' : 'mb.cpt',
                 'land_ice_basal_melt_rate' : 'mb.cpt',
                 'land_ice_basal_x_velocity' : 'velo.cpt',
                 'land_ice_basal_y_velocity' : 'velo.cpt',
                 'surface_temperature' : 'surf_temp.cpt',
                 'lwe_precipitation_rate' : 'mb.cpt',
                 'land_ice_x_velocity' : 'velo.cpt',
                 'land_ice_y_velocity' : 'velo.cpt',
                 'land_ice_temperature' : 'temp.cpt'}
    VARN_MAP = { 'relx' : 'topo.cpt',
                 'presprcp' : 'mb.cpt',
                 'presusrf' : 'ice.cpt',
                 'ubas' : 'velo.cpt',
                 'vbas' : 'velo.cpt',
                 'thk' : 'ice.cpt',
                 'is' : 'ice.cpt',
                 'usurf' : 'ice.cpt',
                 'lsurf' : 'topo.cpt',
                 'topg' : 'topo.cpt',
                 'acab' : 'mb.cpt',
                 'bmlt' : 'mb.cpt',
                 'artm' : 'surf_temp.cpt',
                 'btemp' : 'temp.cpt',
                 'prcp' : 'mb.cpt',
                 'ablt' : 'mb.cpt',
                 'uvel' : 'velo.cpt',
                 'vvel' : 'velo.cpt',
                 'wvel' : 'velo.cpt',
                 'vel'  : 'velo.cpt',
                 'bvel'  : 'velo.cpt',
                 'bvel_tavg'  : 'velo.cpt',
                 'bheatflx': 'gthf.cpt',
                 'litho_temp' : 'litho_temp.cpt',
                 'temp' : 'temp.cpt'}

    def __init__(self,var,filename=None):
        """Initialise.

        var: CF variable
        filename : name of GMT CPT file, if None, do all the automagic."""

        self.var = var
        self.name = var.name
        self.long_name = var.long_name
        self.units = var.units
        self.__cpt = None
        self.__norm = None

        self.__cpt=None
        if filename!=None:
            if os.path.exists(filename):
                fn = filename
            else:
                fn = os.path.join(GC_sharedir,filename)
        else:
            if self.name in self.VARN_MAP:
                fn = os.path.join(GC_sharedir,self.VARN_MAP[self.name])
            elif var.standard_name in self.STDN_MAP:
                fn = os.path.join(GC_sharedir,self.STDN_MAP[var.standard_name])
            else:
                fn = None
        if fn!=None:
            (self.__norm,self.__cpt) = gmtColormap(fn)
                
    def __get_cpt(self):
        return self.__cpt
    colourmap = property(__get_cpt)
    def __get_norm(self):
        return self.__norm
    norm = property(__get_norm)
            
    def __get_title(self):
        title = self.long_name
        if len(self.units)>0 and len(self.units)>0:
            title = title + ' '
        if len(self.units)>0:
            title = title + '[%s]'%self.units
        return title
    title = property(__get_title)

def gmtColormap(filePath,numColours=256):
      try:
          f = open(filePath)
      except:
          print "file ",filePath, "not found"
          return None

      lines = f.readlines()
      f.close()

      x = []
      colours = []
      under = None
      over = None
      bad = None
      colorModel = "RGB"
      for l in lines:
          ls = l.split()
          if l[0] == "#":
             if ls[-1] == "HSV":
                 colorModel = "HSV"
                 continue
             else:
                 continue
          c = convertColours(numpy.array([float(ls[1]),
                                          float(ls[2]),
                                          float(ls[3])],numpy.float),colorModel)
          if ls[0] == "B":
              under = c
          elif ls[0] == "F":
              over = c
          elif ls[0] == "N":
              bad = c
          else:
              x.append(float(ls[0]))
              xtemp =float(ls[4])
              colours.append(c)

      x.append(xtemp)
      x = numpy.array(x,numpy.float)
      colours = numpy.array(colours,float)

      xNorm = (x - x[0])/(x[-1] - x[0])

      norm = matplotlib.colors.BoundaryNorm(x,len(x)-1)
      cmap = matplotlib.colors.ListedColormap(colours,os.path.basename(filePath))

      if over != None:
          cmap.set_over = over
      if under != None:
          cmap.set_under = under
      if bad != None:
          cmap.set_bad = bad
      
      return (norm,cmap)

def convertColours(colour,model):
    if model=='HSV':
        rr,gg,bb = colorsys.hsv_to_rgb(colours[0]/360.,colours[1],colours[2])
        return numpy.array([rr,gg,bb])
    if model=='RGB':
        return colour/255.
    raise RuntimeError

if __name__ == '__main__':
    import sys

    norm,cmap = gmtColormap(sys.argv[1])
    print  norm([-15000+1,0,4000])
    print dir(norm)
    print norm.boundaries
    print dir(cmap)
