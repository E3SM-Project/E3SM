Climate reproducibility testing
===============================

Requiring model changes to pass stringent tests before being accepted as part of E3SM’s main development 
branch is critical for quickly and efficiently producing a trustworthy model.  Depending on their
impacts on model output, code modifications can be classified into three types:

1. Technical changes that continue to produce bit-for-bit identical solutions
2. Changes that cause the model solution to differ,  yet produce a statistically identical climate when
   averaged over a sufficiently long time
3. Changes that lead to a different model climate

Only (3) impacts model climate, and changes of this type should only be implemented within the code 
after an in-depth demonstration of improvement. However, distinguishing between (2) and (3) requires
a comprehensive analysis of both a baseline climate and the currently produced climate.   

Through the CMDV Software project, we've provided a set of climate reproducibility tests to determine
whether or not non-bit-for-bit (nb4b) model changes are climate changing. The current tests provided are:

* **MVK** --  This tests the null hypothesis that the baseline (n) and modified (m) model Short Independent 
  Simulation Ensembles (SISE) represent the same climate state, based on the equality of distribution 
  of each variable's annual global average in the standard monthly model output between the two 
  simulations. The (per variable) null hypothesis uses the non-parametric, two-sample (n and m) 
  Kolmogorov-Smirnov test as the univariate test of of equality of distribution of global means.
 
 
Running the tests
-----------------

These tests are built into E3SM-CIME as system tests and will be launched using the `create_test` scripts. 
*However*, since these tests use high level statistics, they have additional python dependencies which
need to be installed on your system and accessible via the compute nodes (if you're on a batch machine).
Primarily, the statistical analysis of the climates is done through [EVV](https://github.com/LIVVkit/evv4esm) 
(confluence page: [EVV](https://acme-climate.atlassian.net/wiki/spaces/EIDMG/pages/774177270/EVV)) which will
generate a portable test website to describe the results (pass or fail) in detail (see the extended output 
section below). The full set of additional requirements are listed in the `requirements.txt` file in this directory. To install these
dependencies, make sure a python version `> 2.7` is loaded/installed, and then use `pip` like: 

```
pip install --user -r requirements.txt
``` 

*NOTE: Work is underway to include all of these dependencies in the* `e3sm_unified` *conda environment,
and allow CIME to optionally load and use this environment for these tests. Once this is done, you will
not have to manage these dependencies yourself. If you run into problems with getting these dependencies
working on your machine, please open an issue on E3SM's Github and tag @jhkennedy, or send 
Joseph H. Kennedy <kennedyjh@ornl.gov> an email.*

After the dependencies are installed, change to the `$E3SM/cime/scripts` directory (where `$E3SM` is the 
directory containing E3SM). Then to run one of the tests, you will use the `create_test` script like normal. 
To run the `MVK` test and generate a baseline, you would run a command like: 

```
./create_test MVK_PL.ne4_oQU240.FC5AV1C-04P2 -g --baseline-root "/PATH/TO/BASELINE" 
```

And to compare to the baseline, you would run a command like: 


```
./create_test MVK_PL.ne4_oQU240.FC5AV1C-04P2 -c --baseline-root "/PATH/TO/BASELINE" 
```

Importantly, these tests run a 20 member ensemble for at least 13 months (using the last 12 for the 
statistical tests) and, depending on the machine, may take some fiddling to execute within a particular 
queue's wallclock time limit. You may want to over-ride the requested walltime using `--walltime HH:MM:SS` 
option to `create_test`. Additionally, if messing with the wallclock time isn't enough, you can adjust the 
`STOP_N` and `RESUBMIT` setting for the tests in `$E3SM/cime/config/config_tests.xml`. For example, the
MVK test is setup like:   

```xml
  <test NAME="MVK">
    <DESC>climate reproducibility test using the multivariate K-S test</DESC>
    <INFO_DBUG>1</INFO_DBUG>
    <DOUT_S>FALSE</DOUT_S>
    <CONTINUE_RUN>FALSE</CONTINUE_RUN>
    <STOP_OPTION>nmonths</STOP_OPTION>
    <STOP_N>2</STOP_N>
    <REST_OPTION>$STOP_OPTION</REST_OPTION>
    <REST_N>$STOP_N</REST_N>
    <HIST_OPTION>$STOP_OPTION</HIST_OPTION>
    <HIST_N>$STOP_N</HIST_N>
    <RESUBMIT>6</RESUBMIT>
  </test>
```  

which will submit a job to run for two months, stop, and then resubmit itself 6 more times (continuing from 
where it left off), leading to 14 months of simulation time. These settings, combined with the 
`--walltime 02:00:00` option, allowed this test to run successfully on OLCF's Titan machine. The full set of
commands to run the MVK test used on titan are:

*Install dependencies (only need to do once):*
```
module load python
cd $E3SM/cime/scripts/climate_reproducibility
python pip --user -r requirements.txt
``` 

*Generate a baseline:*
```
module load python  # need python 2.7; default is 2.6
cd $E3SM/cime/scripts
./create_test MVK_PL.ne4_oQU240.FC5AV1C-04P2 --baseline-root "${PROJWORK}/cli115/${USER}/baselines" -g --walltime 02:00:00
```

*Compare to a baseline:*
```
module load python  # need python 2.7; default is 2.6
cd $E3SM/cime/scripts
./create_test MVK_PL.ne4_oQU240.FC5AV1C-04P2 --baseline-root "${PROJWORK}/cli115/${USER}/baselines" -c --walltime 02:00:00
```

Test pass/fail and extended output
----------------------------------

When you launch these tests, CIME will ouput the location of the case directory, which will look 
something like this:

```
# On titan:
./create_test MVK_PL.ne4_oQU240.FC5AV1C-04P2 --baseline-root "${PROJWORK}/cli115/${USER}/baselines" -c --walltime 02:00:00
    Creating test directory /ccs/home/$USER/acme_scratch/cli115/MVK_PL.ne4_oQU240.FC5AV1C-04P2.titan_pgi.C.YYYYMMDD_HHMMSS_RANDOMID
```

Let's call that directory `$CASE_DIR`. Once all the jobs are finished, navigate to that directory and 
you can `cat TestStatus` to determine if the test passed or failed by looking at the `BASELINE` status:  

```
cd $CASE_DIR
view TestStatus
    ...
    PASS MVK_PL.ne4_oQU240.FC5AV1C-04P2.titan_pgi BASELINE
    ...

```

To get some basic summary statistics about the test that was run, look in the output of the final job 
submission for EVV's analysis:

 ```
view MVK_PL.ne4_oQU240.FC5AV1C-04P2.titan_pgi.C.YYYYMMDD_HHMMSS_RANDOMID.test.oJOBID
    ...
    -------------------------------------------------------------------- 
                       ______  __      __ __      __                     
                      |  ____| \ \    / / \ \    / /                     
                      | |__     \ \  / /   \ \  / /                      
                      |  __|     \ \/ /     \ \/ /                       
                      | |____     \  /       \  /                        
                      |______|     \/         \/                         
                                                                         
        Extended Verification and Validation for Earth System Models     
    -------------------------------------------------------------------- 
     
      Current run: 2018-08-02 23:24:26 
      User: kennedy 
      OS Type: Linux 3.0.101-0.46.1_1.0502.8871-cray_gem_s 
      Machine: titan-batch7 
       
    ------------------------------------------------------------------- 
     ----------------------------------------------------------------- 
       Beginning extensions test suite  
     ----------------------------------------------------------------- 
     
        Kolmogorov-Smirnov Test: YYYYMMDD_HHMMSS_RANDOMID 
          Variables analyzed: 378 
          Rejecting: 9 
          Critical value: 13.0 
          Ensembles: identical 
     
     ----------------------------------------------------------------- 
       Extensions test suite complete  
     ----------------------------------------------------------------- 
     
    ------------------------------------------------------------------- 
     Done!  Results can be seen in a web browser at: 
       /lustre/atlas/proj-shared/cli115/$USER/MVK_PL.ne4_oQU240.FC5AV1C-04P2.titan_pgi.C.YYYYMMDD_HHMMSS_RANDOMID/run/MVK_PL.ne4_oQU240.FC5AV1C-04P2.titan_pgi.C.YYYYMMDD_HHMMSS_RANDOMID.eve/index.html 
    -------------------------------------------------------------------
    ...
```

EVV also prints the location of the output website where you can see the details of the analysis. For 
the MVK test, you will be able to view per variable Q-Q plots, P-P plots, the K-S test statistic, and 
whether it rejects or accepts the null hypothesis, as well as a description of the test itself -- you 
can see an example of the output website [here](http://livvkit.github.io/evv4esm/).

Please note: the output website uses some JavaScript to render elements of the page (especially figures), 
and opening up the `index.html` file using the `file://` protocol in a web browser will likely not work 
well (most browser have stopped allowing access to "local resources" like JavaScript through the `file://` 
protocol). You can view the website by either copying it to a hosted location (`~/WWW` which is hosted at 
`http://users.nccs.gov/~user` on Titan, for example) or copying it to your local machine and running a 
local http server (included in python!) and viewing it through an address like `http://0.0.0.0:8000/index.html`.  

**For the easiest viewing** we recommend copying the website to your local machine, and using EVV to 
view it. you can install EVV locally by running this command (will work with both python and anaconda 
environments):

```
pip install evv4esm
``` 

Then, copy the website to your local machine, and view it: 


```
# on your local machine
scp -r /lustre/atlas/proj-shared/cli115/$USER/MVK_PL.ne4_oQU240.FC5AV1C-04P2.titan_pgi.C.YYYYMMDD_HHMMSS_RANDOMID/run/MVK_PL.ne4_oQU240.FC5AV1C-04P2.titan_pgi.C.YYYYMMDD_HHMMSS_RANDOMID.eve . 
evv -o MVK_PL.ne4_oQU240.FC5AV1C-04P2.titan_pgi.C.YYYYMMDD_HHMMSS_RANDOMID.eve -s
    --------------------------------------------------------------------
                       ______  __      __ __      __                    
                      |  ____| \ \    / / \ \    / /                    
                      | |__     \ \  / /   \ \  / /                     
                      |  __|     \ \/ /     \ \/ /                      
                      | |____     \  /       \  /                       
                      |______|     \/         \/                        
                                                                        
        Extended Verification and Validation for Earth System Models    
    --------------------------------------------------------------------
    
      Current run: 2018-08-06 15:15:03
      User: fjk
      OS Type: Linux 4.15.0-29-generic
      Machine: pc0101123
      
    
    Serving HTTP on 0.0.0.0 port 8000 (http://0.0.0.0:8000/)
    
    View the generated website by navigating to:
    
        http://0.0.0.0:8000/MVK_PL.ne4_oQU240.FC5AV1C-04P2.titan_pgi.C.YYYYMMDD_HHMMSS_RANDOMID.eve/index.html
    
    Exit by pressing `ctrl+c` to send a keyboard interrupt.
    
```

