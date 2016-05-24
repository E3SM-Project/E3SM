# Copyright (C) 2010
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

__all__ = ['GCNetCDFFile','ScaledMaskedVariable']

import numpy

HAVE_NCFILE=False

if not HAVE_NCFILE:
    HAVE_NCFILE=True
    try:
        from netCDF4 import Dataset as NetCDFFile
    except:
        HAVE_NCFILE=False

if not HAVE_NCFILE:
    HAVE_NCFILE=True
    try:
        from Scientific.IO.NetCDF import NetCDFFile
    except:
        HAVE_NCFILE=False

if not HAVE_NCFILE:
    raise ImportError, "Require either netCDF4 or Scientific.IO.NetCDF module"


defaultFillVals = {'byte'      : -127,
                   'short'     : -32767,
                   'intc'      : -2147483647L,
                   'int_'      : -2147483647L,
                   'longlong'  : -9223372036854775806L,
                   'int8'      : -127,
                   'int16'     : -32767,
                   'int32'     : -2147483647L,
                   'int64'     : -9223372036854775806L,
                   'ubyte'     : 255,
                   'ushort'    : 65535,
                   'uintc'     : 4294967295L,
                   'uint'      : 4294967295L,
                   'ulonglong' : 18446744073709551614L,
                   'uint8'     : 255,
                   'uint16'    : 65535,
                   'uint32'    : 4294967295L,
                   'uint64'    : 18446744073709551614L,
                   'single'    : 9.9692099683868690e+36,
                   'double'    : 9.9692099683868690e+36,
                   'float_'    : 9.9692099683868690e+36,
                   'float32'   : 9.9692099683868690e+36,
                   'float64'   : 9.9692099683868690e+36}


class GCNetCDFFile(object):

    _OWN_ATTRIBS = ['_ncfile','variables']
    
    def __init__(self,filename,mode):
        self._ncfile = NetCDFFile(filename,mode)
        self.variables = {}
        for v in self._ncfile.variables:
            self.variables[v] = MaskedAndScaledVariable(self._ncfile.variables[v])

    def __getattr__(self,name):
        if hasattr(self._ncfile,name):
            return getattr(self._ncfile,name)
        else:
            raise AttributeError, 'no such attribute %s'%name

    def __setattr__(self,name,value):
        if name in self._OWN_ATTRIBS:
            object.__setattr__(self, name, value)
            return
        setattr(self._ncfile,name,value)

    def __delattr__(self,name):
        if name in self.__dict__:
            object.__delattr__(self,name)
        else:
            delattr(self._ncfile,name)

class MaskedAndScaledVariable(object):

    _OWN_ATTRIBS = ['_var','_maskandscale','maskandscale','_sf','_off','dtype']
    
    def __init__(self,var,maskandscale=True):
        self._var = var
        self._maskandscale = True
        self.maskandscale = maskandscale

        self._sf = 1.
        if hasattr(self._var,'scale_factor'):
            self._sf = numpy.squeeze(self._var.scale_factor)
        self._off = 0.
        if hasattr(self._var,'add_offset'):
            self._off = numpy.squeeze(self._var.add_offset)

        if hasattr(self._var,'dtype'):
            self.dtype = self._var.dtype
        else:
            self.dtype = numpy.dtype(self._var.typecode())

    def __len__(self):
        return len(self._var)

    def __getMS(self):
        return self._maskandscale
    def __setMS(self,val):
        assert(isinstance(val,bool))
        self._maskandscale = val
    maskandscale = property(__getMS,__setMS)
        
    def __getitem__(self, index):
        datout = numpy.squeeze(self._var[index])

        if self.maskandscale:
            # mask data
            totalmask = numpy.zeros(datout.shape, bool)
            fill_value = None
            if hasattr(self._var,'missing_value') and (datout == self._var.missing_value).any():
                mask=(data==self._var.missing_value)
                fill_value = self._var.missing_value
                totalmask += mask
            if hasattr(self._var, '_FillValue') and (datout == self._var._FillValue).any():
                mask=data==self._var._FillValue
                if fill_value is None:
                    fill_value = self._var._FillValue
                totalmask += mask
            else:
                fillval = defaultFillVals[str(datout.dtype)]
                if (datout == fillval).any():
                    mask=datout==fillval
                    if fill_value is None:
                        fill_value = fillval
                    totalmask += mask
            if fill_value is not None:
                datout = numpy.ma.masked_array(datout,mask=totalmask,fill_value=fill_value)

            datout = self._off + self._sf * datout
        return datout

    def __setitem__(self,index,value):
        if self.maskandscale:
            v = (value-self._off)/self._sf

            if hasattr(v,'mask'):
                v = v.filled(v.fill_value)
        else:
            v = value
        self._var[index] = v

    def __getattr__(self,name):
        if hasattr(self._var,name):
            return getattr(self._var,name)
        else:
            raise AttributeError, 'no such attribute %s'%name

    def __setattr__(self,name,value):
        if name in self._OWN_ATTRIBS:
            object.__setattr__(self, name, value)
            return
        setattr(self._var,name,value)

    def __delattr__(self,name):
        if name in self.__dict__:
            object.__delattr__(self,name)
        else:
            delattr(self._var,name)

if __name__ == '__main__':
    import sys

    nc = GCNetCDFFile(sys.argv[1],'r+')

    #v = MaskedAndScaledVariable(nc.variables['thk'])
    v = nc.variables['thk']
    v.blub='hello'
    print v[:,30,30]
    print v.long_name

    print v.blub

    del v.blub

    v[:,:,:] = 0
    v[0,0,0] = -2000

    print v[0,0,0]
    print hasattr(v,'blub')

    nc.close()
