from utils import expect, ensure_netcdf4

ensure_netcdf4()

from netCDF4 import Dataset
import numpy as np

import pathlib

###############################################################################
class CheckTendencies(object):
###############################################################################

    ###########################################################################
    def __init__(self,file,variables=None,tendencies=None):
    ###########################################################################

        self._file = pathlib.Path(file).resolve().absolute()
        self._vars  = variables
        self._tends = tendencies

        expect (len(self._vars)==len(self._tends),
                "Error! Input variables and tendencies names lists have different lengths.\n"
                f" - variables : {self._vars}\n"
                f" - tendencies: {self._tends}\n")

        expect (self._file.exists(),
                "Error! File '{}' does not exist.".format(self._file))

    ###########################################################################
    def run(self):
    ###########################################################################

        success = True

        ds = Dataset(self._file,'r')

        # Get time dim extent
        expect ("time" in ds.dimensions.keys() and "time" in ds.variables.keys(),
                f"Error! Missing time dimension and/or variable in file {self._file}")
        nt = ds.dimensions["time"].size
        expect (nt>=2,
                "Error! Not enough time slices in netcdf file (must be at least 2).\n"
                f" - file name: {self._file}\n"
                f" - time size: {nt}")
        time = ds.variables["time"]

        # Loop over var-tend pairs
        for (vname,tname) in zip(self._vars,self._tends):
            # Get var/tend
            expect (vname in ds.variables,
                    "Error! Variable not found in netcdf file.\n"
                    f" - file name: {self._file}\n"
                    f" - var name : {vname}")
            expect (tname in ds.variables,
                    "Error! Variable tendency not found in netcdf file.\n"
                    f" - file name: {self._file}\n"
                    f" - tend name : {tname}")

            var  = ds.variables[vname]
            tend = ds.variables[tname]
            expect (var.dtype==tend.dtype,
                    "Error! Cannot compare variables with different data types"
                    f" - file name : {self._file}\n"
                    f" - var name  : {vname}\n"
                    f" - tend name : {tname}\n"
                    f" - var dtype : {var.dtype}\n"
                    f" - tend dtype: {tend.dtype}\n")
            expect (var.dtype in [np.float32, np.float64],
                    "Error! Only single and double precision supported.\n"
                    f" - file name : {self._file}\n"
                    f" - var name  : {vname}\n"
                    f" - var dtype : {var.dtype}\n")

            tol = np.finfo(var.dtype).eps * 10

            # Sanity checks on var/tend
            expect (var.dimensions==tend.dimensions,
                    "Error! Variable and tendency have different dimensions\n"
                    f" - var name : {vname}\n"
                    f" - tend name: {tname}\n"
                    f" - var dims : {var.dimensions}\n"
                    f" - tend dims: {tend.dimensions}")
            expect (var.dimensions[0]=="time",
                    "Error! Variable first dimension should be time.\n"
                    f" - file name: {self._file}\n"
                    f" - var name : {vname}\n"
                    f" - var dims : {var.dimensions}")
            expect (tend.dimensions[0]=="time",
                    "Error! Tendency first dimension should be time.\n"
                    f" - file name: {self._file}\n"
                    f" - tend name: {tname}\n"
                    f" - tend dims: {tend.dimensions}")

            # Loop over all time intervals [t^k, t^k+1]
            for k in range(1,nt):
                t_k   = time[:][k]
                t_km1 = time[:][k-1]
                dt_days = t_k - t_km1
                dt = dt_days*86400

                var_k   = var[:].take(k,axis=0)
                var_km1 = var[:].take(k-1,axis=0)
                tend_k  = tend[:].take(k,axis=0)

                computed = (var_k - var_km1) / dt
                if not np.isclose(computed,tend_k,tol,tol).all():
                    diff = np.abs(computed - tend_k);
                    max_diff = diff.max()
                    ind_max_1d = np.argmax(diff)
                    shape = np.shape(diff)
                    inds_max = np.unravel_index(ind_max_1d,shape)
                    inds_max_km1 = (k-1,) + np.unravel_index(ind_max_1d,shape)
                    print (f" Check failed for var {vname}.\n"
                           f"  - (t_n,t_n+1) = ({t_km1}, {t_k})\n"
                           f"  - dt = {dt}\n"
                           f"  - rel tol = abs_tol = {tol}\n"
                           f"  - max(|d({vname})/dt - {tname}|) = {max_diff}|)\n"
                           f"  - argmax(|d{vname}/dt - {tname}| = {inds_max}|)\n"
                           f"  - {vname}{(k-1,)+inds_max} = {var_km1[inds_max]}\n"
                           f"  - {vname}{(k,)+inds_max} = {var_k[inds_max]}\n"
                           f"  - {tname}{(k,)+inds_max} = {tend_k[inds_max]}\n"
                           f"  - d{vname}/dt{(k,)+inds_max} = {computed[inds_max]}\n"
                           )
                    success = False

        return success
