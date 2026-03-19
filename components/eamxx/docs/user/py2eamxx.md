# PyEAMxx

PyEAMxx is currently under heavy development and may contain some
rough edges. If you encounter any issues, please report them on the
team on
[github discussions](https://github.com/E3SM-Project/E3SM/labels/eamxx).
Likewise, if you have questions or would like to request features,
please post them on the
[github discussions](https://github.com/E3SM-Project/E3SM/labels/eamxx).

## Quick Start

For now, the only way to use PyEAMxx is to either build it on your own
or use our prebuilt conda binaries. We prefer for you to use the latter.
In a conda environment, please use the following command to install it:

```bash
conda install -c mahf708 pyeamxx=0.0.2
```

It is recommended to use the latest version of PyEAMxx, which is
currently 0.0.2. As you can see, it is a young package with a lot of
potential. We do not guarantee that the API will remain stable, but we
will try to document any changes as frequently as we could.

## Examples

We provide an example to demo calling the radiation process (RRTMGP).
More examples are on the way. If you'd like to add your example,
please feel free to submit a PR.

### RRTMGP

```python
from mpi4py import MPI
import pyeamxx

pyeamxx.init()

dt = 1800
t0_str = "2020-10-10-00000"

ic_file = "/lcrc/group/e3sm/public_html/inputdata/atm/scream/init/screami_unit_tests_ne2np4L72_20220822.nc"
ncols = 218
nlevs = 72
pyeamxx.create_grids_manager(ncols,nlevs, ic_file)

rad_dict = {
    "column_chunk_size": 123,
    "active_gases": ["h2o", "co2", "o3", "n2o", "co" , "ch4", "o2", "n2"],
    "orbital_year": 1990,
    "log_level": "info",
    "do_aerosol_rad": False,
    "rrtmgp_coefficients_file_sw": "/lcrc/group/e3sm/data/inputdata/atm/scream/init/rrtmgp-data-sw-g112-210809.nc",
    "rrtmgp_coefficients_file_lw": "/lcrc/group/e3sm/data/inputdata/atm/scream/init/rrtmgp-data-lw-g128-210809.nc",
    "rrtmgp_cloud_optics_file_sw": "/lcrc/group/e3sm/data/inputdata/atm/scream/init/rrtmgp-cloud-optics-coeffs-sw.nc",
    "rrtmgp_cloud_optics_file_lw": "/lcrc/group/e3sm/data/inputdata/atm/scream/init/rrtmgp-cloud-optics-coeffs-lw.nc",
}

rad = pyeamxx.AtmProc(rad_dict, 'RRTMGP')
rad.read_ic(ic_file)
rad.initialize(t0_str)

t = rad.get_field("T_mid")
tm = t.get()

print(tm[5,5], flush=True)

rad.run(dt)
rad.run(dt)

print(tm[5,5], flush=True)
```
