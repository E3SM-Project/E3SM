from utils import expect
#  import sys
#  sys.path.insert(0,'/home/luca/temp/netcdf4-python/install/usr/local/lib64/python3.6/site-packages')

from netCDF4 import Dataset
from netCDF4 import Variable as NCVar

import re, pathlib

###############################################################################
class NcFileInit(object):
###############################################################################

    ###########################################################################
    def __init__(self,filename,create=False,ne=0,np=0,ncol=0,nlev=0,phys_grid="gll",overwrite=False,
                 add_mid_scalars_1d=None,add_int_scalars_1d=None,
                 add_scalars_2d=None,add_vectors_2d=None,
                 add_mid_scalars_3d=None,add_int_scalars_3d=None,
                 add_mid_vectors_3d=None,add_int_vectors_3d=None,
                 input_file=None,import_variables=None,
                 compute_variables=None,remove_variables=None):
    ###########################################################################

        #### TODO: change var-names to scalars2d, scalars3d_mid, scalars3d_int, vectors_2d...

        import netCDF4
        self._overwrite = overwrite
        self._1d_mid_s  = add_mid_scalars_1d
        self._1d_int_s  = add_int_scalars_1d
        self._2d_s      = add_scalars_2d
        self._2d_v      = add_vectors_2d
        self._3d_mid_s  = add_mid_scalars_3d
        self._3d_int_s  = add_int_scalars_3d
        self._3d_mid_v  = add_mid_vectors_3d
        self._3d_int_v  = add_int_vectors_3d
        self._ifile     = input_file
        self._ivars     = import_variables
        self._cvars     = compute_variables
        self._fname     = filename
        self._rvars     = remove_variables

        if create:
            self.create_database(filename,ncol,ne,np,nlev,phys_grid)
        self._ds = self.get_database(filename,'a')

    ###########################################################################
    def create_database(self,fname,ncol,ne,np,nlev,phys_grid):
    ###########################################################################

        expect (phys_grid=="gll", "Error! So far, only GLL physics grid supported.")

        f = pathlib.Path(fname).resolve().absolute()
        expect (not f.exists(),
                "Error! Output nc file alrady exists. Please, move/delete existing file first.")

        ds = Dataset(f,"w",format="NETCDF4",persist=True)

        expect (ncol >= 1 or (ne>1 and np>1),
                "Error! In order to create a database, do one of the following:\n"
                "   - specify ncol via -ncol flag (ncol>=1)\n"
                "   - specify ne via -ne flag, and possibly np via -np flag (ne>=2, np>=2)\n")
        expect (nlev > 1, "Error! Invalid value for nlev (must be >=2).")
        expect (ncol==0 or ne==0,
                "Error! You can specify either -ncol or -ne, but not both (or I don't know which one to use).")

        # Allow user to specify arbitrary ncols, which can be handy for small unit tests
        if ncol>0:
            ncols = ncol
        else:
            npm1 = np - 1
            nelems = ne*ne*6
            ncols = nelems*npm1*npm1 + 2
            ds.ne   = ne
            ds.np   = np

        ds.createDimension("COL",ncols)
        ds.createDimension("LEV",nlev)
        ds.createDimension("ILEV",nlev+1)

        ds.ncol = ncols
        ds.nlev = nlev

        ds.close()

    ###########################################################################
    def get_database(self,fname,mode):
    ###########################################################################

        f = pathlib.Path(fname).resolve().absolute()
        expect (f.exists(), "Error! Nc file {} does not exists.".format(f))

        ds = Dataset(f,mode,persist=True)

        return ds

    ###########################################################################
    def get_scalar_name_and_value(self,arg_str):
    ###########################################################################
        if arg_str.find('='):
            name_val = arg_str.split('=')
            expect (len(name_val)==2,
                    "Error! Scalar variables must be specified with 'var_name' or 'var_name=value'."
                    "       You cannot have '=' in the variable name.")
            name = name_val[0]
            val = float(name_val[1])
        else:
            name = arg_str
            val = 0.0

        return name, val

    ###########################################################################
    def get_vector_name_and_values(self,arg_str):
    ###########################################################################
        if ':' in arg_str:
            var_dim = arg_str.split(':')
            expect (len(var_dim)==2,
                    "Error! Vector variables must be specified with 'var_name:vec_dim' or 'var_name=[val0,...,valN]'."
                    "       You cannot have ':' in the variable name.")
            name = var_dim[0]
            dim = var_dim[1];
            expect (dim.isdigit() and int(dim)>0,
                    "Error! Invalid vector dimension '{}' for variable {}".format(dim,name))
            values = [0.0] * int(dim)
        else:
            expect ('=' in arg_str,
                    "Error! Vector variables must be specified with 'name:dim' or 'name=[val0,...,valN]'.")

            var_vals = arg_str.split('=')
            name = var_vals[0]
            vals = var_vals[1]
            expect(vals[0]=='[' and vals[-1]==']',
                    "Error! Vector variable must be specified with 'name:dim' or 'name=[val0,...,valN]'.")
            values = [float(d) for d in vals[1:-1].split(',')]
            dim = len(values)

        return name, values

    ###########################################################################
    def correct_dims(self,var,dim_names):
    ###########################################################################
        var_dims = var.get_dims()
        l = len(dim_names)
        if len(var_dims)!=l:
            return False

        for i in range(1,l):
            if var_dims[i].name != dim_names[i]:
                return False

        return True

    ###########################################################################
    def check_var_name(self,var_name):
    ###########################################################################
        # re.match(r'^\w+$', string) is a very compact regexp check, to ensure
        # the string only contains alphanumeric chars and underscores
        expect(re.match(r'^\w+$', var_name),
                "Error! Variable names must contain only alphanumeric characters or underscores.\n")

    ###########################################################################
    def check_overwrite_var(self,var_name):
    ###########################################################################
        expect (not var_name in self._ds.variables.keys() or self._overwrite,
                "Error! Variable '{}' already exists. To overwrite values, use -o flag.".format(var_name))

    ###########################################################################
    def add_scalar_variable(self,ds,name,dims,value):
    ###########################################################################

        self.check_var_name(name)

        self.check_overwrite_var(name)
        if name in ds.variables.keys():
            var = ds.variables[name]
            var_dims = var.get_dims()

            expect (self.correct_dims(var,dims),
                    "Error! Trying to overwrite variable {} using wrong dimensions: ({}) insted of ({}).".
                        format (name,
                                ",".join(dims),
                                ",".join([dim.name for dim in var_dims])))
        else:
            var = ds.createVariable(name,"f8",dims)
        var[:] = value

    ###########################################################################
    def add_vector_variable(self,ds,base_name,dims,values):
    ###########################################################################

        for i in range(len(values)):
            self.add_scalar_variable(ds,"{}_{}".format(base_name,i),dims,values[i])

    ###########################################################################
    def add_variables(self,ds):
    ###########################################################################

        # 1d variables (midpoints and interfaces)
        for item in self._1d_mid_s:
            name, val = self.get_scalar_name_and_value(item)
            self.add_scalar_variable (ds,name,("LEV"),val)

        for item in self._1d_int_s:
            name, val = self.get_scalar_name_and_value(item)
            self.add_scalar_variable (ds,name,("ILEV"),val)

        # 2d variables (scalars and vectors)
        for item in self._2d_s:
            name, val = self.get_scalar_name_and_value(item)
            self.add_scalar_variable (ds,name,("COL"),val)

        for item in self._2d_v:
            name, values = self.get_vector_name_and_values(item)
            self.add_vector_variable (ds,name,("COL"),values)

        # 3d variables (scalars and vectors, at midpoints and interfaces)
        for item in self._3d_mid_s:
            name, val = self.get_scalar_name_and_value(item)
            self.add_scalar_variable (ds,name,("COL","LEV"),val)

        for item in self._3d_int_s:
            name, val = self.get_scalar_name_and_value(item)
            self.add_scalar_variable (ds,name,("COL","ILEV"),val)

        for item in self._3d_mid_v:
            name, values = self.get_vector_name_and_values(item)
            self.add_vector_variable (ds,name,("COL","LEV"),values)

        for item in self._3d_int_v:
            name, values = self.get_vector_name_and_values(item)
            self.add_vector_variable (ds,name,("COL","ILEV"),values)

    ###########################################################################
    def check_dims(self,ds,dims):
    ###########################################################################
        for dim in dims:
            expect (dim in self._ds.dimensions,
                    "Error! Dimension {} not found in the output nc file.".format(dim))
            expect (ds.dimensions[dim].size==self._ds.dimensions[dim].size,
                    "Error! Dimension {} in input file {} has a different extent than in {}.\n"
                    "   - {}: {}\n"
                    "   - {}: {}".format(dim,self._ifile,self._fname,
                                         self._ifile,ds.dimensions[dim].size,
                                         self._fname,self._ds.dimensions[dim].size))

    ###########################################################################
    def import_variables(self):
    ###########################################################################

        if len(self._ivars)>0:
            ds = self.get_database(self._ifile,'r')

            for item in self._ivars:
                if '=' in item:
                    tokens = item.split('=')
                    expect(len(tokens)==2,
                           "Error! Import variable with either 'name' or 'name=name_in', where name_in is\n"
                           "       is the var name in the input file, and name is the var name in the output file.")
                    name_out = tokens[0]
                    name_in  = tokens[1]
                else:
                    name_out = item
                    name_in  = item

                var = ds.variables[name_in]

                # Make sure this var's dims are in our output nc file
                self.check_dims(ds,var.dimensions)
                self.check_overwrite_var(name_out)
                if name_out not in self._ds.variables.keys():
                    self.check_var_name(name_out)
                    self._ds.createVariable(name_out,var.dtype,var.dimensions)

                self._ds.variables[name_out][:] = var[:]
                
            ds.close()

    ###########################################################################
    def compute_variables(self):
    ###########################################################################

        if len(self._cvars)>0:
            from nco import Nco
            n = Nco()

            for expr in self._cvars:
                # Split the expression, to get the output var name
                tokens = expr.split('=')
                expect(len(tokens)==2,"Error! Compute variables with 'var=expr' syntax.")

                var_name = tokens[0]

                self.check_var_name(var_name)
                self.check_overwrite_var(var_name)

                n.ncap2(self._fname,output=self._fname,spt=expr)

    ###########################################################################
    def remove_variables(self):
    ###########################################################################
        if len(self._rvars)>0:
            from nco import Nco
            n = Nco()
            n.ncks(self._fname,output=self._fname,exclude=True,force=True,variable=",".join(self._rvars))

    ###########################################################################
    def process_nc_file(self):
    ###########################################################################

        # Add vars, initing to constant value
        self.add_variables(self._ds)

        # Import vars from another file
        self.import_variables()

        # Compute variables from math expressions (possibly involving other stored vars)
        self.compute_variables()

        # Remove variables
        self.remove_variables()

        # Sync to file and close
        self._ds.sync()
        self._ds.close()

        return True
