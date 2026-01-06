On Compy, one can compile `da.pyx` using `cython` via:

```
module purge
module load cmake/3.19.6 gcc/8.1.0 intel/20.0.0 intelmpi/2020 python/anaconda3-2020.02 netcdf/4.6.3 pnetcdf/1.9.0 mkl/2020
cython da.pyx
```

