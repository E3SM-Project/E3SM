from utils import expect, ensure_netcdf4

ensure_netcdf4()

from netCDF4 import Dataset
import numpy as np

import pathlib

###############################################################################
class CompareNcFiles(object):
###############################################################################

    ###########################################################################
    def __init__(self,src_file,tgt_file=None,compare=None):
    ###########################################################################

        self._src_file = pathlib.Path(src_file).resolve().absolute()
        expect (self._src_file.exists(),
                "Error! File '{}' does not exist.".format(self._src_file))

        self._compare  = compare

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

        ds_src = Dataset(self._src_file,'r')
        ds_tgt = Dataset(self._tgt_file,'r')

        success = True
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
            lvar = ds_src.variables[lname];
            rvar = ds_tgt.variables[rname];

            lvar_rank = len(lvar.dimensions)
            rvar_rank = len(rvar.dimensions)

            expect (len(ldims)==0 or len(ldims)==lvar_rank,
                    f"Invalid slice specification for {lname}.\n"
                    f"  input request: ({','.join(ldims)})\n"
                    f"  variable rank: {lvar_rank}")
            expect (len(rdims)==0 or len(rdims)==rvar_rank,
                    f"Invalid slice specification for {rname}.\n"
                    f"  input request: ({','.join(rdims)})\n"
                    f"  variable rank: {rvar_rank}")


            lslices = [[idim,slice] for idim,slice in enumerate(ldims) if slice!=":"]
            rslices = [[idim,slice] for idim,slice in enumerate(rdims) if slice!=":"]

            lrank = lvar_rank - len(lslices)
            rrank = rvar_rank - len(rslices)

            if lrank!=rrank:
                print (f" Comparison failed. Rank mismatch.\n"
                       f"  - input comparison: {expr}\n"
                       f"  - upon slicing, rank({lname}) = {lrank}\n"
                       f"  - upon slicing, rank({rname}) = {rrank}")
                success = False
                continue

            lvals = self.slice_variable(lvar,lvar[:],lslices)
            rvals = self.slice_variable(rvar,rvar[:],rslices)

            if not np.array_equal(lvals,rvals):
                #  print (f"lvals: {lvals}")
                #  print (f"rvals: {rvals}")
                item = np.argwhere(lvals!=rvals)[0]
                rval = self.slice_variable(rvar,rvals,
                                           [[idim,slice] for idim,slice in enumerate(item)])
                lval = self.slice_variable(lvar,lvals,
                                           [[idim,slice] for idim,slice in enumerate(item)])
                loc = ",".join([str(i+1) for i in item])
                print (f" Comparison failed. Values differ.\n"
                       f"  - input comparison: {expr}\n"
                       f'  - upon slicing, {lname}({loc}) = {lval}\n'
                       f'  - upon slicing, {rname}({loc}) = {rval}')
                success = False

        return success

    ###########################################################################
    def slice_variable(self,var,vals,slices):
    ###########################################################################

        if len(slices)==0:
            return vals

        idim, slice_idx = slices.pop(-1)
        vals = vals.take(int(slice_idx)-1,axis=int(idim))

        return self.slice_variable(var,vals,slices)

    ###########################################################################
    def run(self):
    ###########################################################################

        # Compare variables
        return self.compare_variables()
