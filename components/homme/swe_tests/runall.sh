#!/bin/bash

# list of SWE tests to run
declare -a cases
cases[0]="swtc2"
cases[1]="swtc5"
cases[2]="swtc6"

# FIX THIS LIST UP
# list of output files for each SWE test to collect:
declare -a outputs
outputs[0]="swtc21.nc"
outputs[1]="swtc51.nc"
outputs[2]="swtc61.nc"

args=("$@")
if [ "$#" -lt "2" ]; then
    echo "To compile executables: "
    echo "./runall.sh NMPIRANKS /path/to/output/dir {compile-new, compile-orig}"
    echo ""
    echo "To run all tests, check status or generate plots: "
    echo "./runall.sh NMPIRANKS /path/to/output/dir {run,status,plot}"
    echo ""
    echo "As above, but for a single test (2,5,6):"
    echo "./runall.sh NMPIRANKS /path/to/output/dir {run,status,plot} testno"
    exit 1
fi
WDIR=${args[1]}
opt=${args[2]}

if [ "$#" -eq "2" ]; then
    len=${#cases[@]}
    case_start=0
    case_stop=len-1
else
    case_start=${args[3]}
    case_stop=${args[3]}
fi


#THIS MIGHT BE BROKEN??
SRCDIR="`dirname \"$0\"`"              # relative
SRCDIR="`( cd \"$SRCDIR/..\" && pwd )`"  # absolutized and normalized
echo SOURCE CODE: $SRCDIR
echo OUTPUT DIR:  $WDIR
if [ -z "$SRCDIR" ] ; then
  echo ERROR: for some reason cannot find src directory. Exiting.
  exit 1  # fail
fi
cd $SRCDIR
if [ ! -d dcmip_tests ]; then 
    echo Error: SRCDIR is missing dcmip tests case directory
    exit 1
fi
mkdir -p $WDIR
cd $WDIR
if [ -f README.cmake ] || [ -d cime ]; then 
    # prevent user from accidently installing in the git clone
    echo Error: /path/to/output/dir appears to be the source code directory. Exiting.
    exit 1
fi


if [ "$opt" = "compile" ]; then
  cd $WDIR
  rm -rf *
  ./compile.sweqx >& compile.log
fi


if [ "$opt" = "run" ]; then
#FIX THIS
    if [ ! -x src/sweqx ]; then
        echo Error: sweqx not built during compile step. Exiting.
        exit 1
    fi
#FIX THIS
   for (( i=$case_start; i<=$case_stop; i++ )); do 
       echo BUILD: ${cases[i]}
       cd $WDIR/${cases[i]}
       if [ $? -ne 0 ]; then exit $?; fi
       make install
       ./build.sh
       if [ $? -ne 0 ]; then exit $?; fi

       # remove any existing outputs:
       IFS=" " read -a ARRAY <<< "${outputs[i]}"
       for file in "${ARRAY[@]}"
       do
           \rm -f $file
       done

       echo SUBMIT: ${cases[i]}
       cd $WDIR/${cases[i]}
       sbatch $JOBSCRIPT
   done
fi

#FIX THIS- SHOULD PULL IN THE APPROPRIATE NAMELIST, ETC.
./run_case.sh ${args[0]}

#run swtc2
./run_swtc2.sh ${args[0]}

#run swtc5
./run_swtc5.sh ${args[0]}

#run swtc6
./run_swtc6.sh ${args[0]}



#FIX THIS A LITTLE
missing_count=0
if [ "$opt" = "status" ] || [ "$opt" = "plot" ]; then
   for (( i=$case_start; i<=$case_stop; i++ )); do 
       echo CASE $i: ${cases[i]}
       cd $WDIR/${cases[i]}
       if [ $? -ne 0 ]; then exit $?; fi

       IFS=" " read -a ARRAY <<< "${outputs[i]}"
       for file in "${ARRAY[@]}"
       do
           if [ -f "$file" ]; then
              echo "  found: $file"
           else
              (( missing_count++ ))
              echo "  missing(${missing_count}): $file"
           fi
       done
   done
   echo missing output count: $missing_count
fi



# GENERATE PLOTS
# FIX THIS
if [ "$opt" = "plot" ]; then
    if [ $missing_count -ne 0 ]; then 
        echo ERROR: some output files are missing.  Exiting.
        exit 1; 
    fi   
    mkdir -p $WDIR/latex
    for (( i=$case_start; i<=$case_stop; i++ )); do 
        echo CASE $i: ${cases[i]}
        cd $WDIR/${cases[i]}
        
        IFS=" " read -a ARRAY <<< "${outputs[i]}"
        for file in "${ARRAY[@]}"
        do
            cp $file $WDIR/latex
        done
    done
fi

exit 0
