# ORNL Baseline 

## Environments

```
module load cmake
module load netcdf-c
module load netcdf-fortran
module load hdf5
module load parallel-netcdf

ulimit -n 102400
ulimit -s 819200
ulimit -c unlimited
```

## uELM Configuration 

### Repo and Branch 
```
git clone git@github.com:daliwang/uELM.git
cd uELM
export uELM_home=$PWD
git checkout elm_datmmode_uELM
git submodule update --init --recursive
```

### Case generation, Compile and Run 
```
cd $uELM_home/case_gene/ORNL_baseline
sh uELM_caseGEN_AKSP.sh
```

Generated case directory:
```
$uELM_home/e3sm_cases
```

Run directory:
```
/gpfs/wolf2/cades/cli185/scratch/$USER 
```
