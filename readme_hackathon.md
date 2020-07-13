## How to work with E3SM during Hackathon 

### Getting the source code 

We will use Shan's [OpenACC port of MAM/E3SM](https://github.com/E3SM-Project/E3SM/tree/shz0116/cam/cam_openacc), 
but have rebased it to the current master branch on July 8, 2020 (commit f37bf3e054b3d12ec77191c0463b7bcd0b75db04). 

 
``` shell
git clone -b zz217/cam/openacc_openmp --recursive https://github.com/zhengjizhao/E3SM.git 
```

### Building a test case and run

**For Cori GPUs**

1. Get on a compute Cori GPU node:
```shell
module purge
module load esslurm
salloc -N 1 -c 2 -C gpu -G 1 -A m1759 -t 4:00:00
```
2. Build a case and run: 
```shell
cd ~/E3SM/cime/scripts
./create_test SMS_PS_Ld5.ne4_ne4.FC5AV1C-L -t report --compiler pgiacc --walltime 00:30 --no-run -m cori-gpu --mpilib mvapich2
cd $SCRATCH/e3sm_scratch/cori-gpu/SMS_PS_Ld5.ne4_ne4.FC5AV1C-L.cori-gpu_pgiacc.report
./case.build   
./case.submit
```

The above `create_test` command generates a test case at the directory 
`$SCRATCH/e3sm_scratch/cori-gpu/SMS_PS_Ld5.ne4_ne4.FC5AV1C-L.cori-gpu_pgiacc.report` (`$case_dir`). 
Where `SMS_PS_Ld5.ne4_ne4.FC5AV1C-L` is a test case name. Two larger test cases, `SMS_PS_Ld5.ne16_ne16.FC5AV1C-L` and
`SMS_PS_Ld5.ne30_ne30.FC5AV1C-L`, were used in [the WACCPD paper](link.springer.com/chapter/10.1007/978-3-030-49943-3_3).
The compilation will be done in the `$case_dir/bld` directory, and the resulting executable is `$case_dir/bld/e3sm.exe`. 
The job will run in the `$case_dir/run` directory. 
The timing is in the `$case_dir/run/timing/model_timing.0` file.  

A sample `$case_dir` is available [here](https://portal.nersc.gov/project/m1759/e3sm/SMS_PS_Ld5.ne4_ne4.FC5AV1C-L.cori-gpu_pgiacc.report). 
The standard output from executing the above three commands is available [here](https://portal.nersc.gov/project/m1759/e3sm/screendump.txt) 

!!! note 
    * Remove the `--no-run` option for `create_test` to buid and run in batch (skip the manual executions of case.build and case.submit). 
    * Use the `-h` option to see the available command line options of the `create_test`, `case.build`, and `case.submit` commands.  
    * Use the -u (--use-existing) option of `create_test` to resume after fixing compilation errors instead of starting from scratch.
    * In the `case_dir`, type `./case.build` to update the build after modifying source codes. There are multiple --clean options available if needed.  
    * Use the `preview_run` command to check your job before submission. 
    * The `.case.test` file in the `$case_dir` directory contains the batch directives. 

**For Summit**

```shell
cd ~/E3SM/cime/scripts
./create_test SMS_PS_Ld5.ne4_ne4.FC5AV1C-L -t report --compiler pgiacc --walltime 00:30 --no-run 
cd /gpfs/alpine/csc399/proj-shared/zz217/e3sm_scratch/SMS_PS_Ld5.ne4_ne4.FC5AV1C-L.summit_pgiacc.report
./case.build
./case.submit
```

To use the IBM compilers for OpenMP offloading, use `--compiler ibmgpu` in the above `create_test` command. 
And put the source codes that have the OpenMP target directives in the `~/E3SM/cime/config/e3sm/machines/Depends.Summit.ibmgpu` file.


**Ascent@ORNL**

Similar steps to follow as for Summit above, with addition of the following:

```
./create_test SMS_PS_Ld5.ne4_ne4.FC5AV1C-L -t report --compiler pgiacc --walltime 00:30 --no-run -m ascent --output-root /path/to/
```

Note: Make sure /path/to/ exists are writable by the job.



### Configuration files 

The following files under the `~/E3SM/cime/config/e3sm/machines` directory controls how a test case is built and run. 
 
* `config_machines.xml`
* `config_compilers.xml`
* `config_batch.xml`
* `Depend.<machine>.<compiler>` 
* `../allactive/config_pesall.xml`

Most of the tags in the above `.xml` files are self-descriptive. 
Look for the machine names, such as cori-gpu and summit in the above config files, 
and modify as needed before a case is created. 
While the compiler options for all codes can be put in the `config_compilers.xml` file, 
the `Depend.<machine>.<compiler>` file can be used to put special rules for the files of interest. 
The `config_pesall.xml` specifies how many tasks and threads to use. 
The current setting will run with 6 MPI tasks and 6 GPUs on both Cori GPU and Summit. 

### Recompiling after source code changes

During the Hackathon, we will work on the MAM code of E3SM, which is available at 
`~/E3SM/components/cam/src/physics/cam` directory. After the source code changes, 
you can run the `case.build` command at the `$case_dir` directory to update the `e3sm.exe` binary. 
If you change the compilation flags after a case is built, 
you can do that in the `$case_dir/Depend.<machine>.<compiler>` file.

### Changing run parameters (after a case is bult)

You can change the run parameters after a case is built by directly modifying 
the `.xml` files in the `$case_dir` directory. 
For example, you can change the environment variables in the `env_mach_specific.xml` file. 
However, some run parameters, e.g., the total number of MPI tasks (NTASKS) and 
the number of OpenMP threads per MPI task (NTHRDS) in the `env_mach_pes.xml` cannot be modified 
once the case is setup without first invoking `case_setup --reset` first. 

Here are some examples:
####  Changing the number of tasks
``` shell
    vi $case_dir/env_mach_pes.xml`
    ./case.setup --reset   # you may need to run twice if fails first time
    ./case.build           # The case.setup script recommends to run case.build --clean, but it seems you don't have to
    ./case.submit
```
####  Changing the environment variables
``` shell
    vi env_mach_specific.xml    #e.g, you can add/remove PGI_ACC_TIME and PGI_ACC_NOTIFY for your runs
    ./case.submit
```
####  Changing the jsrun command line for Summit 
``` shell
    get a local copy of `/gpfs/alpine/world-shared/cli115/mpirun.summit` and modify as needed 
    update `env_mach_specific.xml` with the new location of your mpirun.summit copy 
    ./case.submit
```

### Checking correctness after code modifications

The code changes to the MAM code should not change the following output in the `$case_dir/run/atm.log.*` file. 
```shell
[zz217@login1.summit run]$ zgrep "nstep, te" atm.log.215777.200710-220630.gz
 nstep, te        0   0.32987957619845209E+10   0.32987957619845209E+10  -0.00000000000000000E+00   0.98516774202542074E+05
 nstep, te        1   0.32979837976085649E+10   0.32981390622175961E+10   0.21467254203638765E-02   0.98515947713783375E+05
 nstep, te        2   0.32982191263845558E+10   0.32982124014411197E+10  -0.92979410655164691E-04   0.98517279522859841E+05
 nstep, te        3   0.32982645460327196E+10   0.32982415383623142E+10  -0.31810640991090922E-03   0.98516912132544632E+05
 nstep, te        4   0.32981560748194294E+10   0.32981715194480414E+10   0.21353966429054141E-03   0.98516614303255526E+05
 nstep, te        5   0.32982189325749450E+10   0.32982073206766548E+10  -0.16054774822913342E-03   0.98516634832387077E+05
 nstep, te        6   0.32982553713582382E+10   0.32982417004531345E+10  -0.18901627121594125E-03   0.98516424524526432E+05
 nstep, te        7   0.32983368047160912E+10   0.32983225705329466E+10  -0.19680409139130858E-03   0.98516503767995979E+05
 nstep, te        8   0.32984366013566103E+10   0.32984145050625343E+10  -0.30550701387070613E-03   0.98516465063232987E+05
 nstep, te        9   0.32984791751837430E+10   0.32984736878417053E+10  -0.75868829867275226E-04   0.98516553442492514E+05
```

Alternatively, you can just check the md5 hash of the first 10 lines in the file
```shell
cat  atm.log.215777.200710-220630.gz| head -n 10 | md5sum
``` 



