# This script allows use of any of three python netCDF modules:
# Scientific.IO.NetCDF, netCDF4, or pycdf
# To use whichever netCDF module you might have installed put
# from netCDF import *
# in your script.
# Programs should use the Scientific.IO.NetCDF syntax;
# Generally, netCDF4 matches the Scientific.IO.NetCDF syntax and functionality.  However there are some differences which require knowing which module is in use to properly call methods.  The variable netCDF_module is provided to accomplish this.
# If the pycdf module is to be used, an appropriate "translation" of the method calls is provided.
# Written March 16, 2010 by Glen Granzow

try:
  from Scientific.IO.NetCDF import NetCDFFile
  netCDF_module = 'Scientific.IO.NetCDF'
except ImportError:
  try:
    from netCDF4 import Dataset as NetCDFFile
    netCDF_module = 'netCDF4'
  except ImportError:
    try:
      import pycdf
      netCDF_module = 'pycdf'
    except ImportError:
      print 'Unable to import any of the following python modules:'
      print '  Scientific.IO.NetCDF \n  netcdf4 \n  pycdf'
      print 'One of them must be installed.'
      raise ImportError('No netCDF module found')

    def NCtype(value):
      if isinstance(value,int):   return pycdf.NC.INT
      if isinstance(value,float): return pycdf.NC.FLOAT
      if isinstance(value,str):   return pycdf.NC.CHAR   

    class NetCDFFile(object):
      def __init__(self,filename,mode):
        if mode == 'w':
          self.FILE = pycdf.CDF(filename, pycdf.NC.WRITE | pycdf.NC.CREATE | pycdf.NC.TRUNC)
        if mode == 'r':
          self.FILE = pycdf.CDF(filename, pycdf.NC.NOWRITE)
        if mode == 'a':
          self.FILE = pycdf.CDF(filename, pycdf.NC.WRITE)
        self.FILE.automode()
      
      def __setattr__(self,name,value):
        if name == 'FILE':
          object.__setattr__(self,name,value)
        else: # Used to assign global attributes
          self.FILE.attr(name).put(NCtype(value),value)
      
      def __getattr__(self,name):
        if name == 'dimensions':
          return self.FILE.dimensions()
        if name == 'variables':
          dictionary = dict()
          for variable in self.FILE.variables().keys():
            dictionary[variable] = NetCDFvariable(self.FILE.var(variable),None,None,None)
          return dictionary
        global_attributes = self.FILE.attributes()
        if name in global_attributes:
          return global_attributes[name]
        return object.__getattribute__(self,name)
      
      def __dir__(self):
        return self.FILE.attributes().keys()

      def hasattr(self,name):
        return name in dir()

      def createDimension(self,name,size):
        self.FILE.def_dim(name,size)

      def createVariable(self,name,datatype,dimensions): 
        dictNC = {'f':pycdf.NC.FLOAT, 'd':pycdf.NC.DOUBLE, 'i':pycdf.NC.INT, 'c':pycdf.NC.CHAR, 'b':pycdf.NC.BYTE}
        return NetCDFvariable(self.FILE,name,dictNC[datatype],dimensions)

      def sync(self):
        self.FILE.sync()

      def close(self):
        self.FILE.close()

    class NetCDFvariable(object):
      def __init__(self,FILE,name,datatype,dimensions):
        if isinstance(FILE,pycdf.pycdf.CDFVar): 
        # FILE is an already defined netCDF variable (not a file)
          self.VARIABLE = FILE
        else:
        # Create a new variable in the netCDF file
          self.VARIABLE = FILE.def_var(name,datatype,dimensions)
        self.shape = self.VARIABLE.shape()

      def __setitem__(self,elem,data):
        self.VARIABLE[elem] = data

      def __getitem__(self,elem):
        return self.VARIABLE[elem]

      def __setattr__(self,name,value):
        if name in ('VARIABLE','shape'):
          object.__setattr__(self,name,value)
        else: # Used to assign variable attributes
          self.VARIABLE.attr(name).put(NCtype(value),value)
      
      def __getattr__(self,name):
        if name == 'dimensions':
          return self.VARIABLE.dimensions()
        variable_attributes = self.VARIABLE.attributes()
        if name in variable_attributes:
          return variable_attributes[name]
        return object.__getattribute__(self,name)
      
      def assignValue(self,value):
        self.VARIABLE.put(value)

      def getValue(self):
        return self.VARIABLE.get()

      def __dir__(self):
        return self.VARIABLE.attributes().keys()

      def typecode(self):
        NCdict = {pycdf.NC.FLOAT:'f', pycdf.NC.DOUBLE:'d', pycdf.NC.INT:'i', pycdf.NC.CHAR:'c', pycdf.NC.BYTE:'b'}
        return NCdict[self.VARIABLE.inq_type()]

print 'Using',netCDF_module,'for netCDF file I/O'
