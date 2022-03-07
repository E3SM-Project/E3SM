from utils import expect, run_cmd_no_fail

from netCDF4 import Dataset
from nco import Nco
from nco.custom import Atted
import numpy as np

import re, pathlib

###############################################################################
class PopulateNcFile(object):
###############################################################################

    ###########################################################################
    def __init__(self,nc_file,import_file=None,map_file=None,overwrite=False,
                 add_dimensions=None,add_variables=None,import_variables=None,
                 compute_variables=None,remove_variables=None,slice_variables=None,
                 vector_variables=None, prune_history=False):
    ###########################################################################

        self._overwrite = overwrite
        self._prun_hist = prune_history

        self._ofile = pathlib.Path(nc_file).resolve().absolute()
        self._ifile = pathlib.Path(import_file).resolve().absolute()
        self._mfile = pathlib.Path(map_file).resolve().absolute()

        self._adims     = add_dimensions

        self._avars     = add_variables
        self._ivars     = import_variables
        self._cvars     = compute_variables
        self._rvars     = remove_variables
        self._svars     = slice_variables
        self._vvars     = vector_variables

        expect (self._ofile.exists(), "Error! File '{}' does not exist.".format(self._ofile))

    ###########################################################################
    def get_database(self,nc_file,mode,check_dims=True):
    ###########################################################################

        expect (nc_file.exists(), "Error! Nc file {} does not exists.".format(nc_file))

        ds = Dataset(nc_file,mode,persist=True)
        if check_dims:
            expect ('ncol' in ds.dimensions,
                    "Error! NetCDF file '{}' does not contain a 'ncol' dimension."
                    .format(nc_file))
            expect ('ncol' in ds.dimensions and ds.dimensions['ncol'].size>0,
                    "Error! NetCDF file '{}' does not contain a valid 'ncol' dimension."
                    .format(nc_file))
            expect ('lev' in ds.dimensions and ds.dimensions['lev'].size>0,
                    "Error! NetCDF file '{}' does not contain a valid 'lev' dimension."
                    .format(nc_file))

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
        ds = self.get_database(self._ofile,'r')
        expect (not var_name in ds.variables.keys() or self._overwrite,
                "Error! Variable '{}' already exists. To overwrite values, use -o flag.".format(var_name))
        ds.close()

    ###########################################################################
    def add_variable(self,ds,name,dims,value):
    ###########################################################################

        # Check the var name is good (does not contain bad chars)
        self.check_var_name(name)

        # Check we are not overwriting (unless allowed)
        self.check_overwrite_var(name)

        if name in ds.variables.keys():
            # If overwriting, make sure the dimensions are the same
            var = ds.variables[name]
            var_dims = var.get_dims()

            expect (self.same_dims(var.get_dims(),dims),
                    "Error! Trying to overwrite variable {} using wrong dimensions: ({}) insted of ({}).".
                        format (name,
                                ",".join(dims),
                                ",".join([dim.name for dim in var_dims])))
        else:
            var = ds.createVariable(name,"f8",dims)
        var[:] = value

    ###########################################################################
    def get_name(self,name_dims):
    ###########################################################################
        opn = name_dims.find('(')
        cls = name_dims.find(')')

        # Check format
        expect (opn!=-1,"Error! Var declaration should be 'name(dim1,...,dimN)'.")
        expect (cls!=-1,"Error! Var declaration should be 'name(dim1,...,dimN)'.")
        expect (cls>opn,"Error! Var declaration should be 'name(dim1,...,dimN)'.")

        name = name_dims[0:opn]
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
        valid = ["time", "ncol", "lev", "ilev"]
        for dim in dims:
            if dim not in valid:
                expect (dim.isdigit(), "Error! Unexpected dimension '{}'".format(dim))
                return True
        return False

    ###########################################################################
    def get_scalar_dims(self,dims):
    ###########################################################################
        valid = ["ncol", "lev", "ilev"]
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
    def split_braket_list(self,string):
    ###########################################################################
        # Parse a string of the form "[a0,...,aN]", and return the list 'a0,...,aN'
        import re

        valid_re = re.compile(r'[[][a-zA-Z_,]+[]]')
        expect (valid_re.match(string),
                "Error! Braket list should be of the form '[a0,...,aN]'\n"
                "       Input string: {}".format(string))

        return string[1:-1].split(',')

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
                vals = self.split_braket_list(vals_str)
                expect (len(vals)==vec_dim,
                        "Error! Var values specification has the wrong length: {} instead of {}."
                        .format(len(vals),vec_dim))
                try:
                    values = [float(v) for v in vals]
                    return values
                except ValueError:
                    expect(False, "Error! Something went wrong converting strings '{}' to floats.".format(vals))

    ###########################################################################
    def add_dimensions(self):
    ###########################################################################

        ds = self.get_database(self._ofile,'a',check_dims=False)
        for item in self._adims:
            expect ('=' in item,
                    "Error! Add dimension using 'name=length' format.\n")
            tokens = item.split('=')
            expect (len(tokens)==2,
                    "Error! Add dimension using 'name=length' format.\n")

            name = tokens[0]
            size = int(tokens[1])

            expect (not name in ds.dimensions,
                    "Error! Dimension {} already exists in file {}/\n"
                    .format(name,self._ofile))
            expect (size>0,
                    "Error! Invalid extent for dimension {}".format(name))

            ds.createDimension(name,size)

    ###########################################################################
    def add_variables(self):
    ###########################################################################

        ds = self.get_database(self._ofile,'a')
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
                    self.add_variable(ds,"{}_{}".format(name,i),scalar_dims,values[i])
            else:
                value = 0.0 if vals_str=="" else float(vals_str)
                self.add_variable(ds,name,dims,value)

        ds.sync()
        ds.close()

    ###########################################################################
    def check_dims(self,ds_in,ds_out,dims):
    ###########################################################################
        for dim in dims:
            expect (dim in ds_out.dimensions,
                    "Error! Dimension {} not found in the output file '{}'.".format(dim,self._ofile))
            expect (ds_in.dimensions[dim].size==ds_out.dimensions[dim].size,
                    "Error! Dimension {} in input file '{}' has a different extent than in output file '{}'.\n"
                    "   - {}: {}\n"
                    "   - {}: {}".format(dim,self._ifile,self._ofile,
                                         self._ifile,ds_in.dimensions[dim].size,
                                         self._ofile,ds_out.dimensions[dim].size))

    ###########################################################################
    def import_variables(self):
    ###########################################################################

        if len(self._ivars)>0:
            expect (self._ifile.exists(),
                    "Error! Import file '{}' does not exist.".format(self._ifile))

            ds_out = self.get_database(self._ofile,'a')
            ds_in  = self.get_database(self._ifile,'r')

            expect ('ncol' in ds_in.dimensions,
                    "Error! 'ncol' not found in input file dimensions'")
            expect ('lev' in ds_in.dimensions,
                    "Error! 'lev' not found in input file dimensions'")

            ncol_out = ds_out.dimensions['ncol'].size
            nlev_out = ds_out.dimensions['lev'].size
            ncol_in  = ds_in.dimensions['ncol'].size
            nlev_in  = ds_in.dimensions['lev'].size

            ds_in.close()
            ds_out.close()

            expect (nlev_in==nlev_out,
                    "Error! Vertical remapping unavailable, due to ncremap assumption that level idx strides slower than column idx.")

            if ncol_in==ncol_out:
                self.import_variables_no_remap(self._ifile)
            else:
                self.import_variables_horiz_remap()

            # To protect against the possiblity that the input file stored vars with
            # a layout different from scream (e.g., T(time,lev,ncol) instead of
            # T(time,ncol,lev)), we run ncpdq to rearrange (if need be) the dimensions
            run_cmd_no_fail ("ncpdq -a ncol,lev -O {} {}".format(self._ofile,self._ofile))
            run_cmd_no_fail ("ncpdq -a ncol,ilev -O {} {}".format(self._ofile,self._ofile))

    ###########################################################################
    def import_variables_no_remap(self,ifile):
    ###########################################################################

        ds_out = self.get_database(self._ofile,'a')
        ds_in  = self.get_database(ifile,'r')
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

            var = ds_in.variables[name_in]

            # Make sure this var's dims are in our output nc file
            self.check_dims(ds_in,ds_out,var.dimensions)
            self.check_overwrite_var(name_out)
            if name_out not in ds_out.variables.keys():
                self.check_var_name(name_out)
                ds_out.createVariable(name_out,var.dtype,var.dimensions)

            ds_out.variables[name_out][:] = var[:]

        ds_in.close()
        ds_out.sync()
        ds_out.close()

    ###########################################################################
    def import_variables_horiz_remap(self):
    ###########################################################################
        
        # Use a temp file, cause ncremap adds a lot of auxiliary fields
        temp_file = pathlib.Path("./pncf_tmp.nc").resolve().absolute()
        cmd = " ncremap -m {} -i {} -o {} -v {}".\
                format(self._mfile,self._ifile,temp_file,",".join(self._ivars))
        run_cmd_no_fail(cmd)

        # Import only desired vars from temp file to the actual nc file
        self.import_variables_no_remap(temp_file)

        # Clean up temp file
        run_cmd_no_fail("rm {}".format(temp_file))

    ###########################################################################
    def compute_variables(self):
    ###########################################################################

        for expr in self._cvars:
            # Split the expression, to get the output var name
            tokens = expr.split('=')
            expect(len(tokens)==2,"Error! Compute variables with 'var=expr' syntax.")

            var_name = tokens[0]

            self.check_var_name(var_name)
            self.check_overwrite_var(var_name)

            Nco().ncap2(str(self._ofile),output=str(self._ofile),spt=expr)

    ###########################################################################
    def remove_variables(self):
    ###########################################################################
        
        if len(self._rvars)>0:
            Nco().ncks(str(self._ofile),output=str(self._ofile),exclude=True,force=True,variable=",".join(self._rvars))

    ###########################################################################
    def slice_variables(self):
    ###########################################################################

        ds = self.get_database(self._ofile,'a')
        for item in self._svars:
            # Syntax: new_var=old_var(dim_name=N) to extract N-th slice along dim $dim_name
            tokens = item.split('=')

            expect (len(tokens)==3,
                    "Invalid variable declaration: {}\n"
                    "       Syntax: new_var=old_var(dim_name=slice_idx).".format(item))

            # First token is 'new_var'
            new_var_name = tokens[0]
            expect (not new_var_name in ds.variables.keys(),
                    "Variable '{0}' already exists in the database.\n"
                    "       Please, remove it first, using '-rvars {0}'".format(new_var_name))

            # Second token is 'old_var(dim_name'
            expect(len(tokens[1].split('('))==2,
                  "Invalid variable declaration: {}\n"
                  "       Syntax: new_var=old_var(dim_name=slice_idx).".format(item))

            old_var_name = tokens[1].split('(')[0]
            dim_name     = tokens[1].split('(')[1]
            expect (old_var_name in ds.variables,
                    "RHS variable '{}' not found in the database.".format(old_var_name))

            old_var = ds.variables[old_var_name];
            expect (dim_name in old_var.dimensions,
                    "RHS variable '{0}' does not have dimension '{1}'\n"
                    "       {0} layout: {2}".format(old_var_name,dim_name,old_var.dimensions))
            
            # Third token is 'slice_idx)'
            expect(tokens[2][-1]==')',
                   "Invalid variable declaration: {}\n"
                   "       Syntax: new_var=old_var(dim_name=slice_idx).".format(item))

            # Check input slice idx
            try:
                slice_idx = int(tokens[2][:-1])
            except ValueError:
                print("ERROR! Invalid variable declaration: {}\n"
                      "       Syntax: new_var=old_var(dim_name=slice_idx).".format(item))
                raise

            expect (slice_idx>=0 and slice_idx<ds.dimensions[dim_name].size,
                    "Slice index for dim {} out of bounds.\n"
                    "     - slice index: {}\n"
                    "     - dimension size: {}".format(dim_name,slice_idx,ds.dimensions[dim_name].size))

            new_dims = list(old_var.dimensions)
            dim_idx = new_dims.index(dim_name)
            new_dims.remove(dim_name)

            new_var = ds.createVariable(new_var_name,"f8",new_dims)
            new_var[:] = old_var[:].take(slice_idx,axis=dim_idx)

        ds.sync()
        ds.close()


    ###########################################################################
    def vector_variables(self):
    ###########################################################################

        import copy

        ds = self.get_database(self._ofile,'a')
        for item in self._vvars:
            expect('=' in item,
                    "Error! Define vector variables with 'var=[var1,...,varN]'")

            name_and_components = item.split('=')
            expect (len(name_and_components)==2, "Error! Invalid variable declaration: {}".format(item))

            name = name_and_components[0]
            components = self.split_braket_list(name_and_components[1])

            vec_len = len(components)
            dim_tag = "dim{}".format(vec_len)

            if not dim_tag in ds.dimensions.keys():
                ds.createDimension(dim_tag,vec_len)
            vec_dim = ds.dimensions[dim_tag]

            dims = []
            vars = []
            for n in components:
                expect (n in ds.variables.keys(),
                        "Error! Vector variables must be declared in terms of existing variables.\n"
                        "       Variable '{}' not found in the database".format(n))
                var = ds.variables[n]
                if dims == []:
                    dims = list(var.dimensions)
                else:
                    expect (dims==list(var.dimensions),
                            "Error! Vector variables must be declared in terms of variables *of the same dimensions*")
                vars.append(var)

            scalar_dims = copy.deepcopy(dims)
            num_scalar_dims = len(dims)

            col_dim = ds.dimensions["ncol"]

            if col_dim.name in dims:
                # Insert right after the 'ncol' dim
                dims.insert(dims.index(col_dim.name)+1,dim_tag)
                vec_dim_id = dims.index(dim_tag)
            else:
                # Insert right after 'time'
                dims.insert(1,dim_tag)

            vec_var = ds.createVariable(name,"f8",dims)

            for i in range(0,vec_len):
                if num_scalar_dims==2:
                    if scalar_dims[1]=="ncol":
                        vec_var[:,:,i] = vars[i][:,:]
                    elif scalar_dims[1]in ["lev","ilev"]:
                        vec_var[:,i,:] = vars[i][:,:]
                    else:
                        raise NotImplementedError("Scalar variables dimensions not supported.")

                elif num_scalar_dims==3:
                    expect (scalar_dims[1]=="ncol",
                            "Error! Scalar variable dimensions {} not supported.".format(scalar_dims))
                    vec_var[:,:,i,:] = vars[i][:,:,:]
                else:
                    raise NotImplementedError("Scalar variable dimensions not supported.")

        ds.sync()
        ds.close()


    ###########################################################################
    def run(self):
    ###########################################################################

        # Add dimensions
        self.add_dimensions()

        # Add vars, initing to constant value
        self.add_variables()

        # Import vars from another file
        self.import_variables()

        # Create vector variables
        self.vector_variables()

        # Compute variables from math expressions (possibly involving other stored vars)
        self.compute_variables()

        # Extract slices from existing variables
        self.slice_variables()

        # Remove variables
        self.remove_variables()

        if self._prun_hist:
            opt = [
                    Atted(mode='delete',att_name='history',var_name='global'),
                    "-h"
            ]

            Nco().ncatted(input=str(self._ofile),output=str(self._ofile),options=opt,use_shell=True)

        return True
