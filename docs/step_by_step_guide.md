# Useful Aliases
Setting some aliases may be useful in running the model. You can edit your bash settings file to add aliases. Run `source` on that file to start using your new aliases. Examples:
- Chrysalis – `source ~/.bashrc`
- Compy – `source ~/.bash_profile`
- Perlmutter – `source ~/.bash_profile.ext`

Note that the specific file name may differ amongst machines. For example, it might be named `~/.bashrc`.ext.

## Batch Jobs
To check on all batch jobs: 

`alias sqa='squeue -o "%8u %.7a %.4D %.9P %7i %.2t %.10r %.10M %.10l %.8Q %j" --sort=P,-t,-p'`

To check on your batch jobs:

`alias sq='sqa -u $USER'`

The output of `sq` uses several abbreviations: ST = Status, R = running, PD = pending, CG = completing.

## Directories
You will be working in several directories. 

<run_scripts_dir>: `${HOME}/E3SM/scripts`

<code_source_dir>: `${HOME}/E3SM/code`

<simulations_dir> differs amongst machines:

- Anvil/Chrysalis (LCRC) – <simulations_dir>: `/lcrc/group/e3sm/<username>/E3SMv2`
- Compy (PNNL) – <simulations_dir>: `/compyfs/<username>/E3SMv2`
- Perlmutter (NERSC) – <simulations_dir>: `/global/cfs/cdirs/e3sm/<username>/E3SMv2`

So, it may be useful to set the following aliases:
```shell
# Model running
alias run_scripts="cd <run_scripts_dir>"
alias simulations="cd <simulations_dir>"
```

# Configuring the Model Run – Run Script
Start with an example of a run script for a low-resolution coupled simulation:  

We'll use the template script in the E3SM repository, which uses Perlmutter: https://github.com/E3SM-Project/E3SM/blob/master/run_e3sm.template.sh

Create a new run script or copy an existing one (such as the template above). The path to it should be `<run_scripts_dir>/run.<case_name>.sh`

`# Machine and project`
- `readonly MACHINE=pm-cpu`: the name of the machine you’re running on.
- `readonly PROJECT="e3sm"`: SLURM project accounting (typically `e3sm`).

`# Simulation`
- `readonly COMPSET="WCYCL1850"`: compset (configuration)
- `readonly RESOLUTION="ne30pg2_r05_IcoswISC30E3r5"`: resolution. In this example, we have:
  - `ne30pg2`: atmosphere (ne30 dynamics grid -- 30 spectral elements, pg2 physics grid)
  - `r05`: land and river on 1/2 lat/lon grid (commonly referred to as "tri-grid")
  - `IcoswISC30E3r5`: ocean and sea-ice on Icosahedral 30 km mesh with ice shelves cavities (wISC), E3SMv3 (E3) revision r5.
- `readonly CASE_NAME="your_casename"`: case name
- `# readonly CASE_GROUP=""`: This will let you mark multiple cases as part of the same group for later processing (e.g., with PACE). 

> [!IMPORTANT]
> If this is part of a simulation campaign, ask your group lead about using a `CASE_GROUP` label. Otherwise, please use a unique name to distinguish from existing `CASE_GROUP` label names, (e.g., “v2.LR“).

`# Code and compilation`
- `readonly CHECKOUT="latest"`: Date the code was checked out on, in the form `{year}{month}{day}`. The source code will be checked out in `<code_source__dir>/{year}{month}{day}`.
- `readonly BRANCH="master"`: branch the code was checked out from. Valid options include “master”, a branch name, or a git hash. For provenance purposes, it is best to specify the git hash.
- `readonly DEBUG_COMPILE=false`: option to compile with DEBUG flag (leave set to false)

> [!IMPORTANT]
> BEFORE RUNNING: Change `CHECKOUT` to a date string like 20240301. 

> [!IMPORTANT]
> A case is tied to one code base and one executable. That is, if you change `CHECKOUT` or `BRANCH`, then you should also change `CASE_NAME`.

`# Run options`
- `readonly MODEL_START_TYPE="initial"`: specify how the model should start – use initial conditions,  continue from existing restart files, branch, or hybrid (respectively: initial, continue, branch, hybrid).
- `readonly START_DATE="0001-01-01"`: model start date. Typically year 1 for simulations with perpetual (time invariant) forcing or a real year for simulations with transient forcings.

`# Set paths`
- `readonly CODE_ROOT="${HOME}/E3SMv3/code/${CHECKOUT}"`: where the E3SM code will be checked out.
- `readonly CASE_ROOT="/pscratch/sd/r/${USER}/e3sm-scratch/${CASE_NAME}"`: where the results will go. The directory will be `<simulations_dir>/${CASE_NAME}`.

`# Sub-directories`
- `readonly CASE_BUILD_DIR=${CASE_ROOT}/build`: all the compilation files, including the executable.
- `readonly CASE_ARCHIVE_DIR=${CASE_ROOT}/archive`: where short-term archived files will reside.

`# Define type of run`
- `readonly run='XS_2x5_ndays'`: type of simulation to run – i.e, a short test for verification or a long production run. (See next section for details).

`# Coupler history`
- `readonly HIST_OPTION="nyears"`
- `readonly HIST_N="5"`

`# Leave empty (unless you understand what it does)`
- `readonly OLD_EXECUTABLE=""`: this is a somewhat risky option that allows you to re-use a pre-existing executable. This is not recommended because it breaks provenance.

`# --- Toggle flags for what to do ----`

This section controls what operations the script should perform. The run_e3sm script can be invoked multiple times with the user having the option to bypass certain steps by toggling true / false.
- `do_fetch_code=true`: fetch the source code from Github.
- `do_create_newcase=true`: create new case.
- `do_case_setup=true`: case setup.
- `do_case_build=false`: compile.
- `do_case_submit=true`: submit simulation.

The first time the script is called, all the flags should be set to true. Subsequently, the user may decide to bypass code checkout (`do_fetch_code=false`) or compilation (`do_case_build=false`). A user may also prefer to manually submit the job by setting `do_case_submit=false` and then invoking `./case.submit`.

# Running the Model

## Short Tests
Before starting a long production run, it is *highly recommended* to perform a few short tests to verify:
1. The model starts without errors.
2. The model produces BFB (bit-for-bit) results after a restart.
3. The model produces BFB results when changing PE layout.

(1) can spare you from a considerable amount of frustration. Imagine submitting a large job on a Friday afternoon, only to discover Monday morning that the job started to run on Friday evening and died within seconds because of a typo in a namelist variable or input file.

Many code bugs can be caught with (2) and (3). While the E3SM nightly tests should catch such non-BFB errors, it is possible that you’ll be running a slightly different configuration (e.g., a different physics option) for which those tests have not been performed.

### Running Short Tests
The type of run to perform is controlled by the script variable `run`. You should typically perform at least two short tests (two different layouts, with and without restart). Let’s start with a short test using the 'S' (small) PE layout and running for 2x5 days: `readonly run='S_2x5_ndays'`. If you have not fetched and compiled the code, set all the toggle flags to true:
```shell
do_fetch_code=true
do_create_newcase=true
do_case_setup=true
do_case_build=true
do_case_submit=true
```
At this point, execute the run_e3sm script:
```shell
cd <run_scripts_dir>
./run.<case_name>.sh
```
Fetching the code and compiling it will take some time (30 to 45 minutes). Once the script finishes, the test job will have been submitted to the batch queue. 

You can immediately edit the script to prepare for the second short test. In this case, we will be running for 10 days (without restart) using the 'M' (medium PE layout): `readonly run='M_1x10_ndays'`. Since the code has already been fetched and compiled, change the toggle flags:
```shell
do_fetch_code=false
do_create_newcase=true
do_case_setup=true
do_case_build=false
do_case_submit=true
```
and execute the script:
```shell
cd <run_scripts_dir>
./run.<case_name>.sh
```

Since we are bypassing the code fetch and compilation (by re-using the previous executable), the script should only take a few seconds to run and submit the second test.

> [!TIP]
> The short tests use separate output directories, so it is safe to submit and run multiple tests at once. If you’d like, you could submit additional test, for example 10 days with the medium 80 nodes ('M80') layout (`M80_1x10_ndays`).

### Verrifying Results are BFB
Once the short tests are complete, we can confirm the results were bit-for-bit (BFB) the same. All the test output is located under the `tests` directory. To verify that the results are indeed BFB, we extract the global integral from the atmosphere log files (lines starting with ‘nstep, te’) and make sure that they are identical for all tests. 
```shell
cd <simulations_dir>/<case_name>/tests
for test in *
do
  zgrep -h '^ nstep, te ' ${test}/run/atm.log.*.gz | uniq > atm_${test}.txt
done
md5sum *.txt
<hash>  atm_M_1x10_ndays.txt
<matching hash>  atm_M80_1x10_ndays.txt
<matching hash>  atm_S_2x5_ndays.txt
```
If the BFB check fails, you should stop here and understand why. If they succeed, you can now start the production simulation.

## Production Simulation
To prepare for the long production simulation, edit the run_e3sm script and set `readonly run='production'`. In addition, you may need to customize some variables in the code block below to configure run options:

`# Production simulation`
- `readonly PELAYOUT="L"`: 1=single processor, S=small, M=medium, L=large, X1=very large, X2=very very large. Production simulations typically use M or L. The size determines how many nodes will be used. The exact number of nodes will differ amongst machines.
- `readonly WALLTIME="34:00:00"`: maximum wall clock time requested for the batch jobs.
- `readonly STOP_OPTION="nyears"`: see next line
- `readonly STOP_N="50"`: units and length of each segment (i.e., each batch job). E.g, the current configuration stops after 50 years.
- `readonly REST_OPTION="nyears"`: see next line
- `readonly REST_N="5"`: units and frequency for writing restart files (make sure `STOP_N` is a multiple of `REST_N`, otherwise the model will stop without writing a restart fie at the end). E.g., the current configurations saves restart files after every 5 years. 10 restart files will be saved, since `STOP_N=50`.
- `readonly RESUBMIT=”9”`: number of resubmissions beyond the original segment. This simulation would run for a total of 500 years (=inital 50 + 9x50).
- `readonly DO_SHORT_TERM_ARCHIVING=false`: leave set to false if you want to manually run the short term archive.

Since the code has already been fetched and compiled for the short tests, the toggle flags can be set to:
```shell
do_fetch_code=false
do_create_newcase=true
do_case_setup=true
do_case_build=false
do_case_submit=true
```
Finally, execute the script:
```shell
cd <run_scripts_dir>
./run.<case_name>.sh
```
The script will automatically submit the first job. New jobs will be automatically be resubmitted at the end until the total number of segments have been run.

# Looking at Results
`ls <simulations_dir>/<case_name>` explanation of directories:
- `build`: all the stuff to compile. The executable (`e3sm.exe`) is also there. 
- `case_scripts`: the files for your particular simulation. 
- `run`: where all the output will be. Most components (atmosphere, ocean, etc.) have their own log files. The coupler exchanges information between the components. The top level log file will be of the form `run/e3sm.log.*`.  Log prefixes correspond to components of the model:
  - `atm`: atmosphere
  - `cpl`: coupler
  - `ice`: sea ice
  - `lnd`: land
  - `ocn`: ocean
  - `rof`: river runoff

Run `tail -f run/<component>.log.<latest log file>` to keep up with a log in real time.

You can use the `sq` alias defined in the “Useful Aliases” section to check on the status of the job. The `NODE` in the output indicates the number of nodes used and is dependent on the `processor_config` / `PELAYOUT` size.  

> [!NOTE]
> When running on two different machines (such as Compy and Chrysalis) and/or two different compilers, the answers will not be the same, bit-for-bit. It is not possible using floating point operations to get bit-or-bit identical results across machines/compilers.

Logs being compressed to `.gz` files is one of the last steps before the job is done and will indicate successful completion of the segment. `less <log>.gz` will let you directly look at a gzipped log.

# Short Term Archiving
By default, E3SM will store all output files under the `<simulations_dir>/<case_name>/run/` directory. For long simulations, there could 10,000s to 100,000s of output files. Having so many files in a single directory can be very impractical, slowing down simple operations like `ls` to a crawl. CIME includes a short-term archiving utility that will neatly organize output files into a separate `<simulations_dir>/<case_name>/archive/` directory. Short term archiving can be accomplished with the following steps. 

> [!TIP]
> This can be done while the model is still running.

Use `--force-move` to move instead of copying, which can take a long time. Set `--last-date` to the latest date in the simulation you want to archive. You do not have to specify a beginning date.
```shell
cd <simulations_dir>/<case_name>/case_scripts
./case.st_archive --last-date 0051-01-01 --force-move --no-incomplete-logs
ls <e3sm_simulations_dir>/<case_name>/archive
```
Each component of the model has a directory under `archive/`. There are also two additional directories under `archive/`: `logs` holds the gzipped log files and `rest` holds the restart files.
| Component | Directory | File naming pattern |
| --- | --- | --- |
| Atmosphere (Earth Atmospheric Model) | `archive/atm/hist` | `*.eam.h*` |
| Coupler | `archive/cpl/hist` | `*.cpl.h*` |
| Sea Ice (MPAS-Sea-Ice) | `archive/ice/hist` | `*.mpassi.hist.*` |
| Land (Earth Land Model) | `archive/lnd/hist` | `*.elm.h*` |
| Ocean (MPAS-Ocean) | `archive/ocn/hist` | `*.mpaso.hist.*` |
| River Runoff (MOSART) | `archive/rof/hist` | `*.mosart.h*` |

# Performance Information
Model throughput is the number of simulated years per day (SYPD). You can find this with:
```shell
cd <simulations_dir>/<case_name>/case_scripts/timing
grep "simulated_years" e3sm*
```
PACE provides detailed performance information. Go to [PACE](https://pace.ornl.gov/) and enter your username to search for your jobs. You can also simply search by providing the JobID appended to log files (`NNNNN.yymmdd-hhmmss` where `NNNNN` is the SLURM job id). Click on a job ID to see its performance details. “Experiment Details” are listed at the top of the job’s page. There is also a helpful chart detailing how many processors and how much time each component (`atm`, `ocn`, etc.) used. White areas indicate time spent idle/waiting. The area of each box is essentially the "cost = simulation time * number of processors" of the corresponding component.

# Re-Submitting a Job After a Crash
If a job crashes, you can rerun with:
```shell
cd <simulations_dir>/<case_name>/case_scripts
# Make any changes necessary to avoid the crash
./case.submit
```
If you need to change a XML value, the following commands in the `case_scripts` directory are useful:
```shell
> ./xmlquery <variable>                   # Get value of a variable
> ./xmlchange -id <variable> -val <value> # Set value of a variable
```
Before re-submitting:
- Check that the rpointer files all point to the last restart. On very rare occasions, there might be some inconsistency if the model crashed at the end.. Run `head -n 1 rpointer.*` to see the restart date.
- `gzip` all the `*.log` files from the faulty segment so that they get moved during the next short-term archiving. To `gzip` log files from failed jobs, run `gzip *.log.<job ID>*` (where `<job ID>` has no periods/dots in it).
- Delete core or error files, if there are any. MPAS components will sometimes produce a large number of them. The following commands are useful for checking for these files:
  - `ls | grep -in core`
  - `ls | grep -in err`
- If you are re-submitting the _initial_ job, you will need to run `./xmlchange -id CONTINUE_RUN -val TRUE`.

# Post-Processing with zppy
To post-process a model run, do the following steps. 

> [!IMPORTANT]
> To post-process up to year _n_, then you must have short-term archived up to year _n_.

You can ask questions about `zppy` on the [zppy discussion board](https://github.com/E3SM-Project/zppy/discussions/categories/questions).

## Install zppy
Load the E3SM Unified environment. 

> [!TIP]
> The E3SM Unified environment activation commands can be found on [zppy's Getting started page](https://e3sm-project.github.io/zppy/_build/html/main/getting_started.html). Alternatively, they can be found using [Mache](https://github.com/E3SM-Project/mache/tree/main/mache/machines): click the relevant machine and find the `base_path` listed under `[e3sm_unified]` -- the activation command will be `source <base_path>/load_latest_e3sm_unified_<machine_name>.sh`.

If you need a feature in `zppy` that has not yet been included in the E3SM Unified environment, you can construct a [development environment](https://e3sm-project.github.io/zppy/_build/html/main/getting_started.html#b-development-environment).

## Configuration File
In `<run_scripts_dir>`, create a new post-processing configuration file, or copy an existing one, and call it `post.<case_name>.cfg`. 

> [!TIP]
> Good example configuration files can be found in the `zppy` [integration test directory](https://github.com/E3SM-Project/zppy/tree/main/tests/integration/generated) -- `test_complete_run_<machine_name>.cfg

Edit the file and customize as needed. The file is structured with `[section]` and `[[sub-sections]]`. There is a `[default]` section, followed by additional sections for each available zppy task (`climo`, `ts`, `e3sm_diags`, `mpas_analysis`, …). Sub-sections can be used to have multiple instances of a particular task, for example having both regridded monthly and globally averaged time series files. Refer to the `zppy` [schematics documentation](https://e3sm-project.github.io/zppy/_build/html/main/schematics.html) for more details.

The key sections of the configuration file are:

`[default]`
- `input`, `output`, `www` paths will likely need to be edited. 

> [!NOTE]
> The output of your simulation (`<simulations_dir>/<case_name>`) is the _input_ to `zppy`. You can use the same directory for `zppy` output as well, since `zppy` will generate output under `<output>/post`

`[climo]`
- `mapping_file` path may need to be edited.
- Typically you want to generate climatology files every 20,50 years: years = `begin_year:end_yr:averaging_period` – e.g., `years = "1:80:20", "1:50:50",`.

`[ts]`
- `mapping_file` path may need to be edited.
- Typically you want to generate time series files every 10 years – e.g., `years = "1:80:10"`.

`[e3sm_diags]`
- `reference_data_path` may need to be edited.
- `short_name` is a shortened version of the case_name
- `years` should match the `[climo]` section `years`

`[mpas_analysis]`
Years can be specified separately for time series, climatology, and ENSO plots. The lists must have the same lengths and each entry will be mapped to a realization of `mpas_analysis`:
```shell
climo_years ="21-50", "51-100",
enso_years = "11-50", "11-100",
ts_years = "1-50", "1-100",
```
In this particular example, MPAS Analysis will be run twice. The first realization will produce 
climatology plots averaged over years 21-50, ENSO plots for years 11 to 50, and time series plots covering years 1 to 50. The second realization will cover years 51-100 for climatologies, 11-100 for ENSO, and 1-100 for time series.

`[global_time_series]`
- `climo_years` and `ts_years` should match their equivalents in the `[mpas_analysis]` section. 

> [!TIP]
> See the `zppy` [parameters documentation](https://e3sm-project.github.io/zppy/_build/html/main/parameters.html) for more information on parameters.


## Launch zppy
Run `zppy -c post.<case_name>.cfg`. This will submit a number of jobs. Run `sq` to see what jobs are running.

`zppy` automatically handles dependencies of jobs. E.g., `e3sm_diags` jobs are dependent on `climo` and `ts` jobs, so they wait for those to finish. MPAS Analysis jobs re-use computations, so they are chained.

Most jobs run quickly, though E3SM Diags may take around an hour and MPAS Analysis may take several hours.

`zppy` creates a new directory `<simulations_dir>/<case_name>/post`. Each realization will have a shell script (typically `bash`). This is the actual file that has been submitted to the `batch` system. There will also be a log file `*.o<job ID>` as well as a `*.status` file. The status file indicates the state (WAITING, RUNNING, OK, ERROR). These files can be found in `<simulations_dir>/<case_name>/post/scripts`. Once all the jobs are complete, you can check their status.
```shell
cd <simulations_dir>/<case_name>/post/scripts
cat *.status # should be a list of "OK"
grep -v "OK" *.status # lists files without "OK"
```
If you re-run `zppy`, it will check the status of tasks and will skip a task if its status is “OK”. As your simulation progresses, you can update the post-processing years in the configuration file and re-run `zppy`. Newly added task will be submitted, while previously completed ones will be skipped.

## Tasks
If you run `ls <simulations_dir>/<case_name>/post/scripts` you’ll see files like `e3sm_diags_180x360_aave_model_vs_obs_0001-0020.status`. This is one e3sm_diags job. Parts of this file name are explained below:
| Part of File Name | Meaning |
| --- | --- |
| `e3sm_diags` | Task |
| `180x360_aave` | Grid |
| `model_vs_obs` | `model_vs_model` or `model_vs_obs` |
| `0001-0020` | First and last years |
There is also a corresponding output file. It will have the same name but end with `.o<job ID>` instead of `.status`.

## Output
The post-processing output is organized hierarchically. Examples:
- `<e3sm_simulations_dir>/<case_name>/post/atm/180x360_aave/ts/monthly/10yr` has the time series files – one variable per file, in 10 year periods as defined in `<run_scripts_dir>/post.<case_name>.cfg`.  
- `<e3sm_simulations_dir>/<case_name>/post/atm/180x360_aave/clim/20yr` similarly has climatology files for 20 year periods, as defined in `<run_scripts_dir>/post.<case_name>.cfg``.
- `<e3sm_simulations_dir>/<case_name>/post/atm/glb/ts/monthly/10yr` has globally averaged files for 10 years periods as defined in `<run_scripts_dir>/post.<case_name>.cfg`.

# Documenting the Model Run
You should create a Confluence page for your model run in the relevant Confluence space. Use the [Simulation Run Template](https://acme-climate.atlassian.net/wiki/spaces/EWCG/pages/2297299190) as a template. See below for how to fill out this template.

<!-- TODO: where should the Confluence pages be made for v3 (and for each group)? -->

## Code
`code_root_dir` and `tag_name` are defined in `<run_scripts_dir>/run.<case_name>.sh` as `CODE_ROOT` and `BRANCH` respectively.
```shell
cd <code_root_dir>/<tag_name>
git log
```
The commit hash at the top is the most recent commit. Add `“<branch name>, <commit hash>”` to this section of your page.

## Configuration
`Compset` and `Res` are specified on in the PACE “Experiment Details” section. See “Performance Information” above for how to access PACE. Choose the latest job and list these settings on your page. Custom parameters should also be listed. Find these by running:
```shell
cd <run_scripts_dir>
grep -n "EOF >> user_nl" run.<case_name>.sh # Find the line numbers to look at
```
Copy the code blocks after `cat <<EOF >> user_nl_eam`, `cat << EOF >> user_nl_elm`, and `cat << EOF >> user_nl_mosart` to your page.

## Scripts
Push your `<run_scripts_dir>/run.<case_name>.sh` to the relevant GitHub repo/directory -- likely [E3SM Data Docs](https://github.com/E3SM-Project/e3sm_data_docs/tree/main/run_scripts), under `v3/original`. Then link it on this section of your page.

## Output Files
Specify the path to your output files: `<simulations_dir>/<case_name>`.

## Jobs
Fill out a table with columns for “Job”, “Years”, “Nodes”, “SYPD”, and “Notes”.

Log file names will give you the job IDs. Logs are found in `<simulations_dir>/<case_name>/run/`. If you have done short term archiving, then they will instead be in `<simulations_dir>/<case_name>/archive/logs/`.  Use `ls` to see what logs are in the directory. The job ID will be the two-part (period-separated) number after `.log.`.

PACE’s “Experiment Details” section shows `JobID` as well. In the table, link each job ID to its corresponding PACE web page. Note that failed jobs will not have a web page on PACE, but you should still list them in the table.

Use `zgrep "DATE=" <log> | head -n 1` to find the start date. Use `zgrep "DATE=" <log> | tail -n 1` to find the end date. If you would like, you can write a bash function to make this easier:
```shell
get_dates()
{
    for f in atm.log.*.gz; do
        echo $f
        zgrep "DATE=" $f | head -n 1
        zgrep "DATE=" $f | tail -n 1
	echo ""
    done
}
```
(If `zgrep` is unavailable, use `less <log>` to look at a gzipped log file. Scroll down a decent amount to `DATE=` to find the start date. Use `SHIFT+g` to go to the end of the file. Scroll up to `DATE=` to find the end date.)

In the “Years” column specify `<start> - <end>`, with each in `year-month-day` format.

To find the number of nodes, first look at the Processor # / Simulation Time chart on PACE. The x-axis lists the highest MPI rank used, with base-0 numbering of ranks. (PE layouts often don’t fit exactly `N` nodes but instead fill `N-1` nodes and have some number of ranks left over on the final node, leaving some cores on that node unused). Then, find `MPI tasks/node` in the “Experiment Details” section. The number of nodes can then be calculated as `ceil((highest MPI rank + 1)/(MPI tasks/node))`.

The SYPD (simulated years per day) is listed in PACE’s “Experiment Details” section as `Model Throughput`.

In the “Notes” section of the table, mention if a job failed or if you changed anything before re-running a job.

## Global Time Series
> [!NOTE]
> The plots will be available online at the URL corresponding to `<www>/global_time_series/` (where `www` is specified in the `zppy` cfg). See the [E3SM Diags quick guide](https://e3sm-project.github.io/e3sm_diags/_build/html/master/quickguides/quick-guide-general.html) to find the URLs for the web portals on each E3SM machine (listed as `<web_address>`).

You can download the images and then upload them to your Confluence page.

## E3SM Diags
> [!NOTE]
> The plots will be available online at the URL corresponding to `<www>/e3sm_diags/` (where `www` is specified in the `zppy` cfg). See the [E3SM Diags quick guide](https://e3sm-project.github.io/e3sm_diags/_build/html/master/quickguides/quick-guide-general.html) to find the URLs for the web portals on each E3SM machine (listed as `<web_address>`).

Replace the baseline diagnostics in the template's table with relevant ones (e.g., diags on v3 `piControl` and another relevant `v3` run). Add your own diagnostics links in the last columns, labeling them as `<start_year>-<end_year>`.

## MPAS Analysis

> [!NOTE]
> The plots will be available online at the URL corresponding to `<www>/mpas_analysis/` (where `www` is specified in the `zppy` cfg). See the [E3SM Diags quick guide](https://e3sm-project.github.io/e3sm_diags/_build/html/master/quickguides/quick-guide-general.html) to find the URLs for the web portals on each E3SM machine (listed as `<web_address>`).

Make a bulleted list of links, e.g., for `<url_path>/ts_0001-0050_climo_0021-0050/`, create a bullet `“1-50 (time series), 21-50 (climatology)”.

## ILAMB
<!-- TODO: Add this section? -->

> [!NOTE]
> The plots will be available online at the URL corresponding to `<www>/ilamb/` (where `www` is specified in the `zppy` cfg). See the [E3SM Diags quick guide](https://e3sm-project.github.io/e3sm_diags/_build/html/master/quickguides/quick-guide-general.html) to find the URLs for the web portals on each E3SM machine (listed as `<web_address>`).

# Long Term Archiving with zstash
Simulations that are deemed sufficiently valuable should be archived using `zstash` for long-term preservation. You can ask questions about `zstash` on the [zstash discussion board](https://github.com/E3SM-Project/zstash/discussions/categories/questions).

> [!IMPORTANT]
> Compy, Anvil and Chrysalis do not have local HPSS. We rely on NERSC HPSS for long-term archiving. If you are archiving a simulation run on Compy or LCRC (Chrysalis/Anvil), do all of the following steps. If you are archiving a simulation run on NERSC (Perlmutter), skip to step 4.

## 1. Clean up directory
Log into the machine that you ran the simulation on. Remove all `eam.i` files except the latest one. Dates are of the form `<YYYY-MM-DD>`.
```shell
$ cd <simulations_dir>/<case_name>/run
$ ls | wc -l # See how many items are in this directory
$ mv <case_name>.eam.i.<Latest YYYY-MM-DD>-00000.nc tmp.nc
$ rm <case_name>.eam.i.*.nc
$ mv tmp.nc <case_name>.eam.i.<Latest YYYY-MM-DD>-00000.nc
$ ls | wc -l # See how many items are in this directory
```
There may still be more files than is necessary to archive. You can probably remove `*.err`, `*.lock`, `*debug_block*`, `*ocean_block_stats*` files.

## 2. `zstash create` & Transfer to NERSC HPSS
On the machine that you ran the simulation on:

If you don’t have one already, create a directory for utilities, e.g., `utils/`. Then, open a file in that directory called `batch_zstash_create.bash` and paste the following in it, making relevant edits:
```shell
#!/bin/bash
# Run on <machine name>
# Load E3SM Unified
<Command to load the E3SM Unified environment>
# List of experiments to archive with zstash
EXPS=(\
<case_name> \
)
# Loop over simulations
for EXP in "${EXPS[@]}"
do
    echo === Archiving ${EXP} ===
    cd <simulations_dir>/${EXP}
    mkdir -p zstash
    stamp=`date +%Y%m%d`
    time zstash create -v --hpss=globus://nersc/home/<first letter>/<username>/E3SMv2/${EXP} --maxsize 128 . 2>&1 | tee zstash/zstash_create_${stamp}.log
done
```
Load the E3SM Unified environment.
> [!TIP]
> The E3SM Unified environment activation commands can be found on [zppy's Getting started page](https://e3sm-project.github.io/zppy/_build/html/main/getting_started.html). Alternatively, they can be found using [Mache](https://github.com/E3SM-Project/mache/tree/main/mache/machines): click the relevant machine and find the `base_path` listed under `[e3sm_unified]` -- the activation command will be `source <base_path>/load_latest_e3sm_unified_<machine_name>.sh`.

Then, do the following:

```shell
$ screen # Enter screen
$ screen -ls # Output should say "Attached"
$ ./batch_zstash_create.bash 2>&1 | tee batch_zstash_create.log
# Control A D to exit screen
# DO NOT CONTROL X / CONTROL C (as for emacs). This will terminate the task running in screen!!!
$ screen -ls # Output should say "Detached"
$ hostname
# If you log in on another login node, 
# then you will need to ssh to this one to get back to the screen session.
$ tail -f batch_zstash_create.log # Check log without going into screen
# Wait for this to finish
$ screen -r # Return to screen
# Check that output ends with `real`, `user`, `sys` time information
$ exit # Terminate screen
$ screen -ls # The screen should no longer be listed
$ ls <simulations_dir>/<case_name>/zstash
# `index.db`, and a `zstash_create` log should be present
# No tar files should be listed
# If you'd like to know how much space the archive or entire simulation use, run:
$ du -sh <simulations_dir>/<case_name>/zstash
$ du -sh <simulations_dir>/<case_name>
```
Then, on NERSC/Perlmutter:
```shell
$ hsi
$ ls /home/<first letter>/<username>/E3SMv2/<case_name>
# Tar files and `index.db` should be listed.
# Note `| wc -l` doesn't work on hsi
$ exit
```

## 3. `zstash check`
On a NERSC machine (Perlmutter):
```shell
$ cd /global/homes/<first letter>/<username>
$ emacs batch_zstash_check.bash
```
Paste the following in that file, making relevant edits:
```shell
#!/bin/bash
# Run on NERSC dtn
# Load environment that includes zstash
<Command to load the E3SM Unified environment>
# List of experiments to archive with zstash
EXPS=(\
<case_name> \
)
# Loop over simulations
for EXP in "${EXPS[@]}"
do
    echo === Checking ${EXP} ===
    cd /global/cfs/cdirs/e3sm/<username>/E3SMv2
    mkdir -p ${EXP}/zstash
    cd ${EXP}
    stamp=`date +%Y%m%d`
    time zstash check -v --hpss=/home/<first letter>/<username>/E3SMv2/${EXP} --workers 2 2>&1 | tee zstash/zstash_check_${stamp}.log
done
```

> [!TIP]
> If you want to check a long simulation, you can use the `--tars` option to split the checking into more manageable pieces:
> ```shell
> # Starting at 00005a until the end
> zstash check --tars=00005a-
> # Starting from the beginning to 00005a (included)
> zstash check --tars=-00005a
> # Specific range
> zstash check --tars=00005a-00005c
> # Selected tar files
> zstash check --tars=00003e,00004e,000059
> # Mix and match
> zstash check --tars=000030-00003e,00004e,00005a-
> ```

Then, do the following:
```shell
$ ssh dtn01.nersc.gov
$ screen
$ screen -ls # Output should say "Attached"
$ cd /global/homes/<first letter>/<username>
$ ./batch_zstash_check.bash
# Control A D to exit screen
# DO NOT CONTROL X / CONTROL C (as for emacs). This will terminate the task running in screen!!!
$ screen -ls # Output should say "Detached"
$ hostname
# If you log in on another login node, 
# then you will need to ssh to this one to get back to the screen session.
# Wait for the script to finish
# (On Chrysalis, for 165 years of data, this takes ~5 hours)
$ screen -r # Return to screen
# Check that output ends with `INFO: No failures detected when checking the files.`
# as well as listing real`, `user`, `sys` time information
$ exit # Terminate screen
$ screen -ls # The screen should no longer be listed
$ exit # exit data transfer node
$ cd /global/cfs/cdirs/e3sm/<username>/E3SMv2/<case_name>/zstash
$ tail zstash_check_<stamp>.log
# Output should match the output from the screen (without the time information)
```

## 4. Document
On a NERSC machine (Perlmutter):
```shell
$ hsi
$ ls /home/<first letter>/<username>/E3SMv2
# Check that the simulation case is now listed in this directory
$ ls /home/<first letter>/<username>/E3SMv2/<case_name>
# Should match the number of files in the other machine's `<simulations_dir>/<case_name>/zstash`
$ exit
$ cd /global/cfs/cdirs/e3sm/<username>/E3SMv2/<case_name>/zstash
$ ls
# `index.db` and `zstash_check` log should be the only items listed
# https://www2.cisl.ucar.edu/resources/storage-and-file-systems/hpss/managing-files-hsi
# cput will not transfer a file if it exists.
$ hsi "cd /home/<first letter>/<username>/E3SMv2/<case_name>; cput -R <zstash_check log>"
$ hsi
$ ls /home/<first letter>/<username>/E3SMv2/<case_name>
# tar files, `index.db`, `zstash_create` log, and `zstash_check` log should be present
$ exit
```
Update the simulation Confluence page with information regarding this simulation (For Water Cycle’s v2 work, that page is [V2 Simulation Planning](https://acme-climate.atlassian.net/wiki/spaces/ED/pages/2766340117)). In the `zstash archive` column, specify: 
- `/home/<first letter>/<username>/E3SMv2/<case_name>`
- `zstash_create_<stamp>.log`
- `zstash_check_<stamp>.log`

<!-- TODO: Is there a v3 planning page to link? -->

## 5. Delete files
On a NERSC machine (Perlmutter):
```shell
$ hsi
$ ls /home/<first letter>/<username>/E3SMv2/<case_name>
# tar files, `index.db`, `zstash_create` log, and `zstash_check` log should be present
# So, we can safely delete these items on cfs
$ exit
$ ls /global/cfs/cdirs/e3sm/<username>/E3SMv2/<case_name>
# Should match output from the `ls` on `hsi` above
$ cd /global/cfs/cdirs/e3sm/<username>
$ ls E3SMv2
# Only the <case_name> you just transferred to HPSS should be listed.
# We only want to remove that one.
$ rm -rf E3SMv2
```
On the machine that you ran the simulation on:
```shell
$ cd <simulations_dir>/<case_name>
$ ls zstash
# tar files, index.db, `zstash_create` log should be present
$ rm -rf zstash # Remove the zstash directory, keeping original files
$ cd <simulations_dir>
```

## More info
Refer to [zstash's best practices for E3SM](https://e3sm-project.github.io/zstash/_build/html/master/best_practices.html) for details.

# Publishing the simulation data (Optional)
The E3SM Project has a policy to publish all official simulation campaigns once those simulations are documented in publications. Refer to step 3 in [Simulation Data Management](https://acme-climate.atlassian.net/wiki/spaces/DOC/pages/1159594096) for guidance on requesting data publication.
