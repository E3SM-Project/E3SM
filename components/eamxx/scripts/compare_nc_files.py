from utils import expect, ensure_netcdf4

ensure_netcdf4()

from utils import _ensure_pylib_impl
_ensure_pylib_impl("xarray")

import xarray as xr
import numpy as np

import pathlib

###############################################################################
class CompareNcFiles(object):
###############################################################################

    ###########################################################################
    def __init__(self,src_file,tgt_file=None,tolerance=0,compare=None,allow_transpose=False):
    ###########################################################################

        self._src_file = pathlib.Path(src_file).resolve().absolute()
        expect (self._src_file.exists(),
                "Error! File '{}' does not exist.".format(self._src_file))

        self._compare  = compare
        self._tol = tolerance
        self._allow_transpose = allow_transpose

        if tgt_file is None:
            self._tgt_file = self._src_file
        else:
            self._tgt_file = pathlib.Path(tgt_file).resolve().absolute()
            expect (self._tgt_file.exists(),
                    "Error! File '{}' does not exist.".format(self._tgt_file))

    ###########################################################################
    def get_name_and_dims(self,name_dims):
    ###########################################################################
        opn = name_dims.find('(')
        cls = name_dims.find(')')
        nopn = name_dims.count('(')
        ncls = name_dims.count(')')
        last = len(name_dims)-1

        # Check format
        expect (nopn==ncls and nopn<=1,
                f"Format error in var specification string '{name_dims}'. Unmatched parentheses.")
        expect (cls==-1 or cls==last,
                f"Format error in var specification string '{name_dims}'. Traling characters.")

        if nopn>0:
            name = name_dims[0:opn]
            dims = name_dims[opn+1:cls].split(',')

            # Sanity check
            for d in dims:
                expect (d==":" or d.isdigit(),
                        f"Format error in var specification string '{name_dims}'.\n"
                        f"Each dimension must be ':' or a digit.\n")
        else:
            name = name_dims
            dims = []

        return name,dims

    ###########################################################################
    def compare_variables(self):
    ###########################################################################

        ds_src = xr.open_dataset(self._src_file)
        ds_tgt = xr.open_dataset(self._tgt_file)

        success = True

        # If no comparison is passed, compare all variables.
        if self._compare is None or self._compare==[]:
            self._compare = []
            for var in ds_src.variables:
                if var not in ds_tgt.variables:
                    print (f" Comparison failed! Variable not found.\n"
                           f"   - var name: {var}\n"
                           f"   - file name: {self._tgt_file}")
                    success = False
                    continue
                self._compare.append(var+"="+var)

        for expr in self._compare:
            # Split the expression, to get the output var name
            tokens = expr.split('=')
            expect(len(tokens)==2,"Error! Compare variables with 'lhs=rhs' syntax.")

            lhs = tokens[0]
            rhs = tokens[1]

            lname, ldims = self.get_name_and_dims(lhs)
            rname, rdims = self.get_name_and_dims(rhs)

            if lname not in ds_src.variables:
                print (f" Comparison failed! Variable not found.\n"
                       f"   - var name: {lname}\n"
                       f"   - file name: {self._src_file}")
                success = False
                continue
            if rname not in ds_tgt.variables:
                print (f" Comparison failed! Variable not found.\n"
                       f"   - var name: {rname}\n"
                       f"   - file name: {self._tgt_file}")
                success = False
                continue
            lvar = ds_src[lname];
            rvar = ds_tgt[rname];

            lvar_rank = len(lvar.dims)
            rvar_rank = len(rvar.dims)

            expect (len(ldims)==0 or len(ldims)==lvar_rank,
                    f"Invalid slice specification for {lname}.\n"
                    f"  input request: ({','.join(ldims)})\n"
                    f"  variable rank: {lvar_rank}")
            expect (len(rdims)==0 or len(rdims)==rvar_rank,
                    f"Invalid slice specification for {rname}.\n"
                    f"  input request: ({','.join(rdims)})\n"
                    f"  variable rank: {rvar_rank}")

            lslices = {lvar.dims[idim]:int(slice)-1 for idim,slice in enumerate(ldims) if slice!=":"}
            rslices = {rvar.dims[idim]:int(slice)-1 for idim,slice in enumerate(rdims) if slice!=":"}
            lvar_sliced = lvar.sel(lslices)
            rvar_sliced = rvar.sel(rslices)
            expect (set(lvar_sliced.dims) == set(rvar_sliced.dims),
                    f"Error, even when sliced these two elements do not share the same dimensionsn\n"
                    f"   - left var name   : {lname}\n"
                    f"   - right var name  : {rname}\n"
                    f"   - left dimensions : {lvar_sliced.dims}\n"
                    f"   - right dimensions: {rvar_sliced.dims}\n")

            if self._allow_transpose:
              rvar_sliced = rvar_sliced.transpose(*lvar_sliced.dims)

            equal = (lvar_sliced.data==rvar_sliced.data).all()
            if not equal:
                rse = np.sqrt((lvar_sliced.data-rvar_sliced.data)**2)
                nonmatch_count = np.count_nonzero(rse)
                print (f" Comparison failed. Values differ at {nonmatch_count} out of {rse.size} locations.\n"
                      f"  - input comparison: {expr}\n"
                      f'  - max L2 error, {rse.max()}\n'
                      f'  - max L2 location, [{",".join(map(str,(np.array(np.unravel_index(rse.argmax(),rse.shape))+1).tolist()))}]\n'
                      f'  - dimensions, {lvar_sliced.dims}')
                success = False

        return success

    ###########################################################################
    def run(self):
    ###########################################################################

        # Compare variables
        return self.compare_variables()
