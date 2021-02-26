from utils import expect

from netCDF4 import Dataset

import re, pathlib

###############################################################################
class NcFileInit(object):
###############################################################################

    ###########################################################################
    def __init__(self,filename,create=False,ne=0,np=0,nlev=0,phys_grid="gll",
                 mid_scalars_1d=None,int_scalars_1d=None,
                 scalars_2d=None,vectors_2d=None,
                 mid_scalars_3d=None,int_scalars_3d=None,
                 mid_vectors_3d=None,int_vectors_3d=None):
    ###########################################################################

        #### TODO: change var-names to scalars2d, scalars3d_mid, scalars3d_int, vectors_2d...

        self._create   = create
        self._fname    = filename
        self._ne       = ne
        self._np       = np
        self._nlev     = nlev
        self._pg       = phys_grid
        self._1d_mid_s = mid_scalars_1d
        self._1d_int_s = int_scalars_1d
        self._2d_s     = scalars_2d
        self._2d_v     = vectors_2d
        self._3d_mid_s = mid_scalars_3d
        self._3d_int_s = int_scalars_3d
        self._3d_mid_v = mid_vectors_3d
        self._3d_int_v = int_vectors_3d

        expect (filename is not None, "Error! Missing output file name.")
        expect (create is False or ne > 1, "Error! Invalid value for ne.")
        expect (create is False or nlev > 1, "Error! Invalid value for nlev.")
        expect (phys_grid=="gll", "Error! So far, only GLL physics grid supported.")

    ###########################################################################
    def get_database(self):
    ###########################################################################

        if self._create:
            npm1 = self._np - 1
            nelems = self._ne*self._ne*6
            ncols = nelems*npm1*npm1 + 2

            f = pathlib.Path(self._fname).resolve().absolute()
            expect (not f.exists(),
                    "Error! Output nc file alrady exists. Please, move/delete existing file first.")
            ds = Dataset(f,"w",format="NETCDF4",persist=True)

            ds.createDimension("COL",ncols)
            ds.createDimension("LEV",self._nlev)
            ds.createDimension("ILEV",self._nlev+1)
            ds.ne   = self._ne
            ds.np   = self._np
            ds.nlev = self._nlev
        else:
            ds = Dataset(self._fname,"a",persist=True)

        return ds

    ###########################################################################
    def add_variables(self,ds):
    ###########################################################################
        # re.match(r'^\w+$', string) is a very compact regexp check, to ensure
        # the string only contains alphanumeric chars and underscores

        # 1d variables
        for item in self._1d_mid_s:
            if item.find('='):
                name_val = item.split('=')
                expect (len(name_val)==2,
                        "Error! Scalar variables must be specified with 'var_name' or 'var_name=value'."
                        "       You cannot have '=' in the variable name.")
                name = name_val[0]
                val = float(name_val[1])
            else:
                name = item
                val = 0.0

            expect(re.match(r'^\w+$', name),
                    "Error! Variable names must contain only alphanumeric characters or underscores.\n")

            var = ds.createVariable(name,"f8",("LEV"))
            var[:] = val
        for item in self._1d_int_s:
            if item.find('='):
                name_val = item.split('=')
                expect (len(name_val)==2,
                        "Error! Scalar variables must be specified with 'var_name' or 'var_name=value'."
                        "       You cannot have '=' in the variable name.")
                name = name_val[0]
                val = float(name_val[1])
            else:
                name = item
                val = 0.0

            expect(re.match(r'^\w+$', name),
                    "Error! Variable names must contain only alphanumeric characters or underscores.\n")

            var = ds.createVariable(name,"f8",("ILEV"))
            var[:] = val

        # 2d variables
        for item in self._2d_s:
            if item.find('='):
                name_val = item.split('=')
                expect (len(name_val)==2,
                        "Error! Scalar variables must be specified with 'var_name' or 'var_name=value'."
                        "       You cannot have '=' in the variable name.")
                name = name_val[0]
                val = float(name_val[1])
            else:
                name = item
                val = 0.0

            expect(re.match(r'^\w+$', name),
                   "Error! Variable names must contain only alphanumeric characters or underscores.\n")

            var = ds.createVariable(name,"f8",("COL"))
            var[:] = val
        for item in self._2d_v:
            values = []
            if ':' in item:
                var_dim = item.split(':')
                expect (len(var_dim)==2,
                        "Error! Vector variables must be specified with 'var_name:vec_dim' or 'var_name=[val0,...,valN]'."
                        "       You cannot have ':' in the variable name.")
                name = var_dim[0]
                dim = var_dim[1];
                expect (dim.isdigit() and int(dim)>0,
                        "Error! Invalid vector dimension '{}' for variable {}".format(dim,name))
                values = [0.0] * int(dim)
            else:
                expect ('=' in item,
                        "Error! Vector variable must be specified with 'name:dim' or 'name=[val0,...,valN]'.")

                var_vals = item.split('=')
                name = var_vals[0]
                vals = var_vals[1]
                expect(vals[0]=='[' and vals[-1]==']',
                        "Error! Vector variable must be specified with 'name:dim' or 'name=[val0,...,valN]'.")
                values = [float(d) for d in vals[1:-1].split(',')]
                dim = len(values)

            for i in range(dim):
                var = ds.createVariable("{}_{}".format(name,i),"f8",("COL"))
                var[:] = values[i]

        # 3d variables
        for item in self._3d_mid_s:
            if item.find('='):
                name_val = item.split('=')
                expect (len(name_val)==2,
                        "Error! Scalar variables must be specified with 'var_name' or 'var_name=value'."
                        "       You cannot have '=' in the variable name.")
                name = name_val[0]
                val = float(name_val[1])
            else:
                name = item
                val = 0.0

            expect(re.match(r'^\w+$', name),
                   "Error! Variable names must contain only alphanumeric characters or underscores.\n")

            var = ds.createVariable(name,"f8",("COL","LEV"))
            var[:] = val
        for item in self._3d_int_s:
            if item.find('='):
                name_val = item.split('=')
                expect (len(name_val)==2,
                        "Error! Scalar variables must be specified with 'var_name' or 'var_name=value'."
                        "       You cannot have '=' in the variable name.")
                name = name_val[0]
                val = float(name_val[1])
            else:
                name = item
                val = 0.0

            expect(re.match(r'^\w+$', name),
                   "Error! Variable names must contain only alphanumeric characters or underscores.\n")

            var = ds.createVariable(name,"f8",("COL","ILEV"))
            var[:] = val

        for item in self._3d_mid_v:
            values = []
            if ':' in item:
                var_dim = item.split(':')
                expect (len(var_dim)==2,
                        "Error! Vector variables must be specified with 'var_name:vec_dim' or 'var_name=[val0,...,valN]'."
                        "       You cannot have ':' in the variable name.")
                name = var_dim[0]
                dim = var_dim[1];
                expect (dim.isdigit() and int(dim)>0,
                        "Error! Invalid vector dimension '{}' for variable {}".format(dim,name))
                values = [0.0] * int(dim)
            else:
                expect ('=' in item,
                        "Error! Vector variable must be specified with 'name:dim' or 'name=[val0,...,valN]'.")

                var_vals = item.split('=')
                name = var_vals[0]
                vals = var_vals[1]
                expect(vals[0]=='[' and vals[-1]==']',
                        "Error! Vector variable must be specified with 'name:dim' or 'name=[val0,...,valN]'.")
                values = [float(d) for d in vals[1:-1].split(',')]
                dim = len(values)

            for i in range(dim):
                var = ds.createVariable("{}_{}".format(name,i),"f8",("COL","LEV"))
                var[:] = values[i]
        for item in self._3d_int_v:
            values = []
            if ':' in item:
                var_dim = item.split(':')
                expect (len(var_dim)==2,
                        "Error! Vector variables must be specified with 'var_name:vec_dim' or 'var_name=[val0,...,valN]'."
                        "       You cannot have ':' in the variable name.")
                name = var_dim[0]
                dim = var_dim[1];
                expect (dim.isdigit() and int(dim)>0,
                        "Error! Invalid vector dimension '{}' for variable {}".format(dim,name))
                values = [0.0] * int(dim)
            else:
                expect ('=' in item,
                        "Error! Vector variable must be specified with 'name:dim' or 'name=[val0,...,valN]'.")

                var_vals = item.split('=')
                name = var_vals[0]
                vals = var_vals[1]
                expect(vals[0]=='[' and vals[-1]==']',
                        "Error! Vector variable must be specified with 'name:dim' or 'name=[val0,...,valN]'.")
                values = [float(d) for d in vals[1:-1].split(',')]
                dim = len(values)

            for i in range(dim):
                var = ds.createVariable("{}_{}".format(name,i),"f8",("COL","ILEV"))
                var[:] = values[i]

    ###########################################################################
    def process_nc_file(self):
    ###########################################################################
        ds = self.get_database()

        self.add_variables(ds)

        ds.sync()
        ds.close()

        return True
