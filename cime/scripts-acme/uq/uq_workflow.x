#!/bin/bash -e

#################################################################
#                                                               #      
#                  Non-standard requirements:                   #      
#                                                               #      
# ncdump              : to read netcdf output                   #
#                       (either make sure it is in the path, or #
#                        add the path to it below)              #
#                                                               #                                  
# python's numpy      : for data manipulations                  # 
# python's matplotlib : to plot                                 # 
#         (on Titan, run 'module load python_matplotlib/1.2.1') #     
#                                                               #      
# UQTk                : for UQ analysis                         #      
#         (is included, see below to install)                   #      
#                                                               #      
#################################################################

source prepare_env.x


## Usage string
USAGE="Syntax: uq_workflow.x <RUNMODE>, where RUNMODE=-2,-1,0,1,2,3,4,5\n \
       uq_workflow.x -2  : generates input samples\n \
       uq_workflow.x -1  : reads preciously generated input samples\n \
       uq_workflow.x 0   : reads output netcdf file\n \
       uq_workflow.x 1   : diagnostic plots of the outputs\n \
       uq_workflow.x 2   : creates PC surrogate\n \
       uq_workflow.x 3   : computes sensitivities and plots results\n \
       uq_workflow.x 4   : plots sensitivities across sites for each output\n \
       uq_workflow.x 5   : plots sensitivities across outputs for each site\n"

## Argument check
if [ $# -ne 1 ]; then
    printf "$USAGE"
    exit
fi

## Get the argument (mode of run)
RUNMODE=$1

## Hardwired arguments
OUTLIST="TLAI GPP TOTVEGC" # TOTVEGC EFLX_LH_TOT TOTSOMC" ## String containing the output names to be explored
                   ## "TLAI GPP TOTVEGC EFLX_LH_TOT TOTSOMC"
MAXSITE=10 #96       ## Will analyze sites 1 to MAXSITE
                   ## (can start with a low number for testing)
NCORES=16           ## Number of cores to multithread site analysis in parallel
                   
###########################################################################
## Generate input parameter ensemble
if [ "$RUNMODE" = "-2" ]; then
    FULLPARAMFILE=Parameters_TitanUQ_ACME.csv
    REMOVEPARS="lf_flig fr_flig"
    NSAM=10000
    cd files
    ../prepr.py $FULLPARAMFILE $NSAM $REMOVEPARS
    tar -czf prepr.tar pdomain_* pnames_* input_*
    echo "Saved the generated files in files/prepr.tar"
    cd ..

## Extracts the previously generated input parameter ensemble
## for reproducibility, and to match with the default output file
elif [ "$RUNMODE" = "-1" ]; then
    cd files
    tar -xvzf prepr_saved.tar pdomain_* pnames_* input_*
    echo "Extracted input files from files/prepr_saved.tar"
    cd ..

## Run ALM as a black-box
# Input: input_clm.dat
# ALM evaluation as black-box
# Output: FLUXNETUQ_*.nc

## Read the resuling NetCDF file
elif [ "$RUNMODE" = "0" ]; then
    echo "Reading the NetCDF file"
    for OBS in $OUTLIST; do
        ./postp.py 0 $OBS
    done

## Explore/visualize the ALM runs
## (use as an example only - many hardwired parameters in postp.py)
elif [ "$RUNMODE" = "1" ]; then
    echo "Postprocessing some outputs for visual inspection"
    # Explore specific observable at a given site
    # by plotting, e.g.
    # ./postp.py 1 TLAI 64
    # or do them all together
    ./create_tasks.x 1 $MAXSITE "$OUTLIST"
    ./uqmulti.py tasks $NCORES #> uqmulti1.log

## Run surrogate construction
elif [ "$RUNMODE" = "2" ]; then
    echo "Constructing Polynomial Chaos surrogates"
    # for a specific observable at a given site
    # ./postp.py 2 GPP 38
    # or do them all together
    ./create_tasks.x 2 $MAXSITE "$OUTLIST"
    ./uqmulti.py tasks $NCORES #> uqmulti2.log

## Compute sensitivities and plot surrogate/sensitivity results
elif [ "$RUNMODE" = "3" ]; then
    echo "Computing sensitivities and plotting results"
    # for a specific observable at a given site
    # ./postp.py 3 GPP 38
    # or do them all together
    ./create_tasks.x 3 $MAXSITE "$OUTLIST"
    ./uqmulti.py tasks $NCORES #> uqmulti3.log

## Plot sensitivities across all sites for each output
elif [ "$RUNMODE" = "4" ]; then
    for OBS in $OUTLIST; do
        echo "Plotting sensitivities across all sites for $OBS"
        ./postp.py 4 $OBS $MAXSITE
    done

## Plot sensitivities across all outputs for each site
elif [ "$RUNMODE" = "5" ]; then
    for ((siteId=1;siteId<=$MAXSITE;siteId++)); do
    #for siteId in {1..$MAXSITE}; do
        echo "Plotting sensitivities across all outputs for site $siteId"
        ./postp.py 5 $siteId $OUTLIST
    done

else
    printf "$USAGE"
fi

# Cleanup
rm -rf *.pyc
