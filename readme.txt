## How to work with E3SM during Hackathon 

# Getting the source code 
git clone -b zz217/cam/openacc_openmp --recursive https://github.com/zhengjizhao/E3SM.git 

# Building a test case and run

## For Cori GPUs
module purge
module load esslurm
salloc -N 1 -c 2 -C gpu -G 1 -A m1759 -t 4:00:00

cd ~/E3SM/cime/scripts
./create_test SMS_PS_Ld5.ne4_ne4.FC5AV1C-L -t report --compiler pgiacc --walltime 00:30 --no-run -m cori-gpu --mpilib mvapich2
cd $SCRATCH/e3sm_scratch/cori-gpu/SMS_PS_Ld5.ne4_ne4.FC5AV1C-L.cori-gpu_pgiacc.report
./case.build   
./case.submit

##For Summit
cd ~/E3SM/cime/scripts
./create_test SMS_PS_Ld5.ne4_ne4.FC5AV1C-L -t report --compiler pgiacc --walltime 00:30 --no-run 
cd /gpfs/alpine/csc399/proj-shared/zz217/e3sm_scratch/SMS_PS_Ld5.ne4_ne4.FC5AV1C-L.summit_pgiacc.report
./case.build
./case.submit

##For OpenMP offload with IBM compilers: 
./create_test SMS_PS_Ld5.ne4_ne4.FC5AV1C-L -t report --compiler ibmgpu --walltime 00:30 --no-run 
Add source codes that have the OpenMP target directives in the `~/E3SM/cime/config/e3sm/machines/Depends.Summit.ibmgpu` file.

### Two larger teset cases are
`SMS_PS_Ld5.ne16_ne16.FC5AV1C-L` 
`SMS_PS_Ld5.ne30_ne30.FC5AV1C-L`

# Configuration files 
~/E3SM/cime/config/e3sm/machines/config_machines.xml
~/E3SM/cime/config/e3sm/machines/config_compilers.xml
~/E3SM/cime/config/e3sm/machines/config_batch.xml
~/E3SM/cime/config/e3sm/machines/Depend.<machine>.<compiler> 
~/E3SM/cime/config/e3sm/allactive/config_pesall.xml

# Checking correctness after code modifications
zgrep "nstep, te" atm.log.*.gz #or
cat  atm.log.*.gz| head -n 10 | md5sum

