from utils import expect

from netCDF4 import Dataset

import re, pathlib

###############################################################################
class NcFileInit(object):
###############################################################################

    ###########################################################################
    def __init__(self,filename,create=False,ne=0,np=0,ncol=0,nlev=0,phys_grid="gll",
                 overwrite=False,import_file=None,
                 add_variables=None,import_variables=None,
                 compute_variables=None,remove_variables=None):
    ###########################################################################

        self._overwrite = overwrite

        self._avars     = add_variables
        self._ifile     = import_file
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

        ds.createDimension("time",1)
        ds.createDimension("ncol",ncols)
        ds.createDimension("lev",nlev)
        ds.createDimension("ilev",nlev+1)

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
    def same_dims(self,dims1,dims2):
    ###########################################################################
        l = len(dims1)
        if len(dims2)!=l:
            return False

        for i in range(1,l):
            if dims1[i] != dims2[i]:
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
    def add_variable(self,name,dims,value):
    ###########################################################################

        # Check the var name is good (does not contain bad chars)
        self.check_var_name(name)

        # Check we are not overwriting (unless allowed)
        self.check_overwrite_var(name)

        if name in self._ds.variables.keys():
            # If overwriting, make sure the dimensions are the same
            var = self._ds.variables[name]
            var_dims = var.get_dims()

            expect (self.same_dims(var.get_dims(),dims),
                    "Error! Trying to overwrite variable {} using wrong dimensions: ({}) insted of ({}).".
                        format (name,
                                ",".join(dims),
                                ",".join([dim.name for dim in var_dims])))
        else:
            var = self._ds.createVariable(name,"f8",dims)
        var[:] = value

    ###########################################################################
    def get_name(self,name_dims):
    ###########################################################################
        open = name_dims.find('(')
        close = name_dims.find(')')

        # Check format
        expect (open!=-1,"Error! Var declaration should be 'name(dim1,...,dimN)'.")
        expect (close!=-1,"Error! Var declaration should be 'name(dim1,...,dimN)'.")
        expect (close>open,"Error! Var declaration should be 'name(dim1,...,dimN)'.")

        name = name_dims[0:open]
        return name

    ###########################################################################
    def get_dims(self,name_dims):
    ###########################################################################
        opn = name_dims.find('(')
        cls = name_dims.find(')')

        # Check format
        expect (opn!=-1,"Error! Var declaration should be 'name(dim1,...,dimN)'.")
        expect (cls!=-1,"Error! Var declaration should be 'name(dim1,...,dimN)'.")
        expect (cls>opn,"Error! Var declaration should be 'name(dim1,...,dimN)'.")
        expect (cls==len(name_dims)-1,"Error! Var declaration should be 'name(dim1,...,dimN)'.")

        dims = name_dims[opn+1:cls].split(',')
        expect (len(dims)>0,"Error! Var declaration should be 'name(dim1,...,dimN)'.")

        return dims

    ###########################################################################
    def is_vector_layout(self,dims):
    ###########################################################################
        valid = ["COL", "LEV", "ILEV"]
        for dim in dims:
            if dim not in valid:
                expect (dim.isdigit(), "Error! Unexpected dimension '{}'".format(dim))
                return True
        return False

    ###########################################################################
    def get_scalar_dims(self,dims):
    ###########################################################################
        valid = ["COL", "LEV", "ILEV"]
        s_dims = []
        vec_dim_id = -1 
        for i in range(0,len(dims)):
            if dims[i] in valid:
                s_dims.append(dims[i])
            else:
                expect (vec_dim_id==-1,
                        "Error! Multiple integer extents found in dims specification '{}'.\n"
                        "       Only vectors are supported, for non-scalar layouts.".format(dims))
                vec_dim_id = i

        expect(vec_dim_id>0, "Error! Something went wrong while detecting vector dim id from '{}'.".format(dims))

        return vec_dim_id, s_dims

    ###########################################################################
    def get_values(self,vals_str,vec_dim):
    ###########################################################################
        if vals_str=="":
            # User did not specify values. Use all zeros.
            return [0]*vec_dim
        else:
            try:
                # Try to convert vals_str to a single float. If successful,
                # the user passed a single value for all vector components
                value = float(vals_str)
                return [value]*vec_dim
            except ValueError:
                # Didn't work. Then we must have a [v1,...,vN] format
                opn = vals_str.find('[')
                cls = vals_str.find(']')
                expect (opn!=-1,"Error! Var values specification should be '..=val' or '..=[val1,...,valN]'.")
                expect (cls!=-1,"Error! Var values specification should be '..=val' or '..=[val1,...,valN]'.")
                expect (cls>opn,"Error! Var values specification should be '..=val' or '..=[val1,...,valN]'.")
                expect (cls==len(vals_str)-1,"Error! Var values specification should be '..=val' or '..=[val1,...,valN]'.")
                vals = vals_str[opn+1:cls].split(',')
                expect (len(vals)==vec_dim,
                        "Error! Var values specification has the wrong length: {} instead of {}."
                        .format(len(vals),vec_dim))
                try:
                    values = [float(v) for v in vals]
                    return values
                except ValueError:
                    expect(False, "Error! Something went wrong converting strings '{}' to floats.".format(vals))


    ###########################################################################
    def add_variables(self):
    ###########################################################################

        for item in self._avars:
            if '=' in item:
                # User provided initial values
                name_dims_vals = item.split('=')
                expect (len(name_dims_vals)==2, "Error! Invalid variable declaration: {}".format(item))
                name_dims = name_dims_vals[0]
                vals_str = name_dims_vals[1]
            else:
                name_dims = item
                vals_str = ""

            # From the string name(dim1,...,dimN) extract name and dimensions
            name = self.get_name(name_dims)
            dims = self.get_dims(name_dims)

            is_vector = self.is_vector_layout(dims)

            if is_vector:
                # From the list (dim1,...,dimN), check if it is a vector field,
                # and if so, get the idx of the vector dimension, the extent
                # along that dimension, and the dims list without the vector dim.
                vec_dim_id, scalar_dims = self.get_scalar_dims(dims)

                vec_dim = 1 if vec_dim_id==-1 else int(dims[vec_dim_id])

                # From the string after the = (if any), get the initialization
                # values. The string can be a single value (for scalar or vector
                # fields) or a list of values [v1,...,vn] (for vector field)
                values = self.get_values(vals_str,vec_dim)

                for i in range(0,len(values)):
                    self.add_variable("{}_{}".format(name,i),scalar_dims,values[i])
            else:
                value = 0.0 if vals_str=="" else float(vals_str)
                self.add_variable(name,dims,value)

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
        self.add_variables()

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
