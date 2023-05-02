#!/bin/bash 
shopt -s extglob

usage() {
cat << EOF | less
NAME   

 E3SM-Polar-Developer.sh - Extensible polar climate development script for E3SM


SYNOPSIS 

 E3SM-Polar-Developer.sh [-s|--sandbox]   (CODE SANBDOX NAME)

                         [-f|--fetch]     (FORK) (BRANCH)

                         [-t|--tests]     (SUITE)

                         [-n|--newcase]
                         [-d|--duration]  (MONTHS)
                         [-c|--config]    (CONFIG)
                         [-e|--emc]
                         [-k|--kombo]     (NAMEFILE)

                         [-b|--build] 
                         [-q|--qsubmit]

                         [-a|--analyze]   (COMPARISON CASE)
                         [-v|--verbose]   

                         [-h|--help]       

 This is an extensible tool to develop new polar physics and BGC in E3SM from initial
 inception through to production simulations. This script creates a structure of
 sandboxes, cases, simulations and methods to compare them. Example workflows
 are provided after first listing the available options. This script version focuses 
 on developing the sea ice model in the coupled E3SM framework.  Planned version 2 
 extensions to the full polar system are listed at the bottom of this help page. This 
 script must use one of the -s, -c, or -h options, the former in combination with 
 options listed below.  Casenames created from the -s sandbox are assigned according 
 to the combined information in the -s, -c, -d, -e, and -k options.  The script will 
 work using the following machines: ${joinedmachines%,}.


OPTIONS

 [-s|--sandbox]   - One of three possible mandatory options giving the name of the 
                    code sandbox, which will be located under the directory:
                    ${CODE_BASE}
                    When specified without any other options, this lists all cases
                    already created associated with this sandbox, or specified without
                    an argument, this lists code sandboxes already created, if any.
 
 [-f|--fetch]     - Clones to a sandbox code from an E3SM repository. You may specify 
                    arguments for the code FORK and BRANCH if different from defaults:
                    FORK:   ${FORK}, BRANCH: ${BRANCH} 
                    The branch may be specified as a git hash. If the sandbox already
                    exists, it will only be overwritten if explicitly requested when
                    the user is challenged with a Y/N prompt. 

 [-t|--tests]     - Runs preset E3SM tests for the given sandbox from testsuite(s): 

                    ${joinedtestsuites%,}
                    This option waits for the tests to complate, then provides a 
                    step-by-step PASS/FAIL summary. The test cases are located under:
                    ${CASE_BASE} 
                    with nomenclature described below.  An argument must be provided 
                    specifying which test suite is to be used, with valid arguments 
                    being 'ice', 'ocn', 'land', or 'atm'.  As an example, when you are
                    making changes to the sea ice code, it's a good idea to regularly
                    run the 'ice' test suite.

 [-n|--newcase]   - Creates or overwrites a case for a given code sandbox. If the
                    sandbox needs to be overwritten, the user will be challenged 
                    for confirmation prior to proceeding. The case name has nomenclature
                    as described below, located in:
                    ${CASE_BASE}

 [-d|--duration]  - Change the maximum duration in months with annual restarts from
                    a default run length of ${DURATION} months. A specification of more 
                    than 12 months will round up the number of months to the nearest year; 
                    A duration of 11 months will run for that amount of simulated
                    time, but 13 months will run for 2 simulated years (24 months),
                    with resubmition after the first 12-month simulation. Using this 
                    option, where the duration exceeds 12 months, one may set up a 
                    self-sustaining year-on-year production run. 

 [-c|--config]    - Specifies one of many available configurations. This option
                    does not need to be specified if the default  ${CONFIG}-Case is used. 
                    To list all availale configurations, enter 'E3SM-Polar-Developer.sh -c'
                    on the command line. Only configurations relevant to E3SM Phase 3 will 
                    be maintained.

 [-e|--emc]       - Check for energy and mass conservation using analysis members
                    in the sea ice and ocean models. Note that this option slows doen the 
                    model so may not always be desired, but is important to be used
                    intermittently when testing new code. Coming in Version 2: graphing
                    conservation properties with combined -e -a options.

 [-k|--kombo]     - Test namelist combinations for each relevant E3SM component model 
                    using the input namelist file '<kombotag>.nlk' giving options to be 
                    tested for particular models using the following format:

                    [cpl]
                    [eam]
                    NAMELIST_ENTRY1 = {false,true}
                    [eamxx]
                    [elm]
                    [mosart]
                    [mpaso]
                    NAMELIST_ENTRY2 = optionnumber1
                    [omegaxx]
                    [mpassi]
                    NAMELIST_ENTRY3 = {'string1','string2','string3'} 
                    NAMELIST_ENTRY4 = {number1,number2} 
                    [mali]
                    [ww3]
                    [datm]
                    [dlnd]
                    [docn]
                    [dice]
                    [drof]
                    [dwav]

                    Here, the square brackets indicate the model for which namelist
                    options are to be changed, and the curly brackets indicate which
                    namelist combinations are to be tested for that particular model. If 
                    no curly brackets are provided, as for NAMELIST_ENTRY2, then the 
                    namelist will be changed according to this one specification. In
                    the above example, a total of 12 different namelist combinations
                    will be tested in seperate model simulations using the same code base
                    and the same code build. In this example, the different run directories 
                    will be labeled  as run.k000 through run.k011 within the case
                    directory. Note that in the above example, all of EAM, MPAS-O, and
                    MPAS-SI must be active in the configuration specified with -c, or 
                    else the build will fail with an error message. 
           
                    This option should be entered as -k <kombotag>.nlk where
                    <kombotag> is also used in the name of the case.

 [-b|--build]     - Builds or rebuilds the case for given 

 [-q|--qsubmit]   - Submits case(s) to the queue.

 [-a|--analyze]   - Analyze the output of a given simulation, and compare it against
                    equivalent simulations from other sandboxes. At present this 
                    option provides a difference between integrations. Version 2 will 
                    provide graphical interpretation of energy, mass and state 
                    variables.

 [-v|--verbose]   - Switches on verbose mode for analysis, including information on 
                    min and max differences in fields for the same mesh, and on energy
                    conservation specifics across coupler accounting, instead of just
                    total mass and energy [CURRENTLY IN DEVELOPMENT].

 [-h|--help]      - Provides this help page.


SANDBOXES, CASES, COMBINATIONS, BRANCHES, and TEST areas.

 When a new sandbox is created, it is assigned a directory according to its name,
 situated at ${CODE_BASE}/<sandbox>. 

 When a new  case is created, it is assigned a name according to the configuration, 
 duration, energy and mass conservation, namelist combination, and sandbox branch
 <config><duration>.<kombotag>.emc.<sandbox>.<sbranch>.<sremote>.<machine> under:
 ${CASE_BASE}
 If one or other of the -k and -e options are omitted, the <kombotag>.emc modifiers
 are respectively removed from the case name. Note that the -e and -k options must
 always be specified to point to a case that includes them.

 When multiple (n>1) combinations are specified with the -k option, seperate 
 case_script and run directories are sequentially numbered under the case directory 
 as case_scripts.k000 ... case_scripts.k<n> and as run.k000 ... run.k<n> directories. 

 The <sbranch> descriptor creates a summary branch name from the full branch name.
 E.g. The sbranch name from branch "eclare108213/seaice/icepack-integration" is
 "icepack-integration", the last field seperated by back slashes. As a result, if you 
 checkout a new branch in the sandbox, a new case will be created, allowing results 
 to be compared between branches for a single sandbox.

 The <sremote> descriptor creates a summary fork name from the full remote name.
 E.g. The sremote name for git@github.com:eclare108213/E3SM is 'eclare108213'

 If the -t option is used, the tests are located under the umbrella directory of
 Test.<suite>.<sandbox>.<sbranch>.<sremote>.<machine>.


EXAMPLE-A WORKFLOW - Simple three-month baseline simulation

 1) E3SM-Polar-Developer.sh -s baseline -fnb 

 2) E3SM-Polar-Developer.sh -s baseline -q

 3) E3SM-Polar-Developer.sh -s baseline -a

 This clones ${FORK} and sets up the sandbox 'baseline' ${BRANCH} 
 in ${CODE_BASE}, then creates the default ${DURATION}-Month ${CONFIG}-Case, submits it 
 to the queue, and indicates when the run is complete. The case directory 
 ${CONFIG}${DURATION}.baseline.master.E3SM-Project.${MACHINES[$MACH]} is located in 
 ${CASE_BASE}.


EXAMPLE-B WORKFLOW - Create a comparable case for a branch with changes

 1) E3SM-Polar-Developer.sh -s change1 -f git@github.com:eclare108213/E3SM \\
                                          eclare108213/seaice/icepack-integration

 2) E3SM-Polar-Developer.sh -c

 3) E3SM-Polar-Developer.sh -s change1 -c D -d 3 -nb -k nset1.nlk

    Where nset1.nlk is a file that includes a line testing the
    column package and icepack:

    [mpassi]
    config_column_physics_type = {'column_package','icepack'}

 4) E3SM-Polar-Developer.sh -s change1 -c ${CONFIG} -d ${DURATION} -k nset1.nlk -q

 5) E3SM-Polar-Developer.sh -s change1 -c ${CONFIG} -d ${DURATION} -k nset1.nlk \\
                            -a ${CONFIG}${DURATION}.baseline.master.E3SM-Project.${MACHINES[$MACH]}

 This creates the sandbox 'change1' from the fork git@github.com:eclare108213/E3SM with
 branch eclare108213/seaice/icepack-integration. The second step lists all of the
 available configurations' compsets and resolutions, including the "D" configuration. 
 The third step generates a new case and builds for the two combinations in the file
 nset1.nlk in case D3.nset1.change1.icepack-integration.eclare108213.${MACHINES[$MACH]} in 
 ${CASE_BASE} with subdirectories *.k000 and *.k001 for 
 each of the respective namelist combinations. The fourth step queues the two cases, 
 and the final step compares these two simulations with the baseline case established
 in the EXAMPLE-A WORKFLOW above.


EXAMPLE-C WORKFLOW - set up a cold-start fully-coupled production of SORRM 

 1) E3SM-Polar-Developer.sh -c

 2) E3SM-Polar-Developer.sh -s baseline -c BSORRM -d 24 -e -nb

 3) E3SM-Polar-Developer.sh -s baseline -c BSORRM -d 24 -e -q

 4) E3SM-Polar-Developer.sh -s baseline -c BSORRM -d 24 -e -a

 After checking the available configurations in the first step, the second step
 uses the baseline sandbox created in EXAMPLE-A WORKFLOW to create and build a
 new case for the fully coupled model on the SORRM ice-ocean mesh with standard
 atmosphere and land resolution that runs for 12 months and resubmits for a second
 12 months using steps (2) and (3). Step (4) tells you if and when it is complete. 
 The -e option generates a timeseries of conservation residuals for both ice and
 ocean.


EXAMPLE-D WORKFLOW - check the sandbox passes designated E3SM test suites

 1) E3SM-Polar-Developer.sh -s baseline -t ice

 From the sandbox created in EXAMPLE-A WORKFLOW, checks code integrity and for a BFB 
 match aginst E3SM-Project/E3SM master for testsuite: 'e3sm_ice_testsuite' 
 The test case appears under Test.ice.baseline.master.E3SM-Project.${MACHINES[$MACH]} 


EXAMPLE-E WORKFLOW - Create a new branch in baseline and create a new case

 1) cd ${CODE_BASE}/code/baseline

 2) git checkout -b eclare108213/seaice/basechange

 3) ...make changes in eclare108213/seaice/basechange branch...

 4) E3SM-Polar-Developer.sh -s baseline -nb

 Assuming the changes in the new branch 'eclare108213/seaice/basechange' in the 
 baseline sandbox are not in some way erroneous, this will create and build a new
 case called ${CONFIG}${DURATION}.baseline.master.basechange.${MACHINES[$MACH]}.


EXAMPLE-F WORKFLOW - Switch remote and branch in the baseline sandbox for new case

 1) cd ${CODE_BASE}/code/baseline

 2) git remote add eclare git@github.com:eclare108213/E3SM

 3) git fetch eclare

 4) git checkout -b eclare108213/seaice/icepack-integration \\
                    eclare/eclare108213/seaice/icepack-integration

 5) E3SM-Polar-Developer.sh -s baseline -nb

 6) E3SM-Polar-Developer.sh -s baseline -q

 7) E3SM-Polar-Developer.sh -s baseline -a D3.baseline.E3SM-Project.master.${MACHINES[$MACH]}
 
 Assuming the code in the imported branch 'eclare108213/seaice/icepack-integration' 
 works, this will create and build a new case from the fork eclare108213/E3SM called 
 ${CONFIG}${DURATION}.baseline.eclare108213.icepack-integration.${MACHINES[$MACH]}.
 It will then submit the job, and compare in (7) against the master baseline in
 E3SM-Project.


EXAMPLE-G WORKFLOW - Checking energy and mass conservation 

 You can check for energy and mass conservation in the sea ice and ocean models by adding
 the -e option. For this, we will draw on the simple EXAMPLE-A WORKFLOW to demonstrate,
 what this will do.  So we follow the same steps, except for cloning, but add -e: 

 1) E3SM-Polar-Developer.sh -s baseline -e -nb 

 2) E3SM-Polar-Developer.sh -s baseline -e -q

 3) E3SM-Polar-Developer.sh -s baseline -e -a

 This will generate a new case for the baseline sandbox titled:
 ${CONFIG}${DURATION}.emc.baseline.master.E3SM-Project.${MACHINES[$MACH]} 
 Here, the '.emc.' indicates that energy, water, salt and carbon conservation 
 will be checked ('emc' stands for energy and mass conservation), and the sea ice 
 model will generate a file <casename>.mpassi.hist.am.conservationCheck.<year>.nc
 indicating whether the model is conserving internally.




AUTHORS

 Andrew Roberts, Jon Wolfe, Elizabeth Hunke, Darin Comeau, Nicole Jeffery, Erin Thomas


VERSIONS

 This Version: 1.0 March 2023

 Version 2 is being designed to include the following additional features:
 1) Extension to Perlmutter 
 2) Addition of continue, hybrid, and branch simulations for production
 3) Graphing of energy and mass conservation, extension to ocean
 4) Switch to upcoming V3 Meshes in place of V2 meshes


EOF
}

#---------------------------------
main() {

# For debugging this script, uncomment line below
#set -x

# For running E3SM with debug on, change line below to true
readonly DEBUG_COMPILE=false

# Make directories created by this script world-readable
umask 022

# Get the command line options and set defaults
get_configuration $*

# Copy script into case_script directory for provenance
copy_script $*

# Fetch code from Github
clone_code

if [ -d ${CODE_ROOT} ]; then

 if [ "${do_test_suite,,}" == "true" ]; then

  e3sm_testsuites

  if [ -d ${CASE_ROOT} ]; then

   printf "\n--- Test root: $(cd ${CASE_BASE} && dirs +0)/\n    ${CASE_NAME}\n" 

  else

   printf "\n--- No test root created \n" 
   
  fi

 else

  # Create namelist combinations
  namelist_kombo

  # Create case
  create_newcase

  if [ -d ${CASE_ROOT} ]; then

   # Build
   case_build

   if [ "${exit_script,,}" != "true" ]; then

    # Submit
    case_queue

    # Generate analysis scripts
    case_analyze
 
   fi

   printf "\n--- Case root: $(cd ${CASE_BASE} && dirs +0)/\n    ${CASE_NAME}\n" 

  else

   printf "\n--- \x1B[31mNo case yet exists for this configuration\e[0m ---\n"

  fi

 fi

 # Provide directory stemming from home with tilde shorthand
 printf "\n--- Code root: $(cd ${CODE_ROOT} && dirs +0)\n" 

else

 printf "\n--- \x1B[31mNo code yet cloned\e[0m ---\n"

fi

# Provide directory stemming from home with tilde shorthand
printf "\n--- Provenance root: $(cd `dirname ${SCRIPT_PROVENANCE_DIR}` && dirs +0)/\n"
printf "    `basename ${SCRIPT_PROVENANCE_DIR}`\n"

printf "\n--- Script saved: ${SCRIPT_PROVENANCE_NAME}\n\n"

}

#---------------------------------
get_configuration() {

    declare -Ag  MAXCOMBINATIONS
    declare -Ag  TOTALCOMBINATIONS
    declare -Ag  COMPSET
    declare -Ag  RESOLUTION
    declare -Ag  PELAYOUT
    declare -Ag  CONFIGDESCRIPTION

    # Change forks and branches here if you wish to switch defaults
    DEFAULT_FORK="git@github.com:E3SM-Project/E3SM.git" 
    DEFAULT_BRANCH="master" 

    # Change default configuration and run length here
    DEFAULT_CONFIG='D'
    DEFAULT_DURATION=3

    # --- Auto-detect the machine and set project ---
    local host_node=`uname -n | cut -f 1 -d . | rev | cut -c 2- | rev`
    MACHINES=("anvil" "chrysalis")
    printf -v joinedmachines '%s,' "${MACHINES[@]}"

    PROJECTS=("condo" "e3sm")
    SCRATCHS=("/lcrc/group/e3sm" "/lcrc/group/e3sm")
    if [ "${host_node,,}" == "blueslogin" ]; then
     MACH=0
    elif [ "${host_node,,}" == "chrlogin" ]; then
     MACH=1
    else
     printf "\n--- ERROR: Unable to detect valid machines: ${joinedmachines%,}\n\n"
     exit 1
    fi

    # --- Default Settings ---
    do_code_clone=false
    do_test_suite=false
    do_case_create=false
    do_case_build=false
    do_case_queue=false
    do_case_kombo=false
    do_energy_mass=false
    do_case_analyze=false

    exit_script=false

    readonly CASE_GROUP='E3SM-Polar'
    readonly CODE_BASE="${HOME}/${CASE_GROUP}/code"
    readonly CASE_BASE="${SCRATCHS[$MACH]}/${USER}/${CASE_GROUP}"

    SANDBOX=''
    CASE2_NAME=''

    FORK="${DEFAULT_FORK}"
    BRANCH="${DEFAULT_BRANCH}"
    HASH=''

    CONFIG="${DEFAULT_CONFIG}"
    DURATION="${DEFAULT_DURATION}"

    TESTSUITES=("e3sm_ice_developer" "e3sm_land_developer" 
                "e3sm_atm_developer" "e3sm_ocn_developer")

    printf -v joinedtestsuites '%s\n                    ' "${TESTSUITES[@]}"

    # --- Parse options ---
    while getopts :sft:d:cek:nbqah-: OPT; do
      if [ "$OPT" = "-" ]; then   # long option: reformulate OPT and OPTARG
        OPT="${OPTARG%%=*}"       # extract long option name
        OPTARG="${OPTARG#$OPT}"   # extract long option argument (may be empty)
        OPTARG="${OPTARG#=}"      # if long option argument, remove assigning `=`
      fi
      case "$OPT" in
        s | sandbox )  eval "SANDBOX=\${$((OPTIND))}" # allows one optional argument
                       if [[ $SANDBOX =~ ^-.* ]] || [ -z ${SANDBOX} ] ; then
                        SANDBOX="";
                       elif [ -n ${SANDBOX} ]; then
                        OPTIND=$((OPTIND+1));
                       fi ;;
        f | fetch )    eval "FORK=\${$((OPTIND))}" # allows two optional arguments
                       if [[ $FORK  =~ ^-.* ]] || [ -z ${FORK} ] ; then 
                        FORK="${DEFAULT_FORK}";
                       elif [ -n ${FORK} ]; then 
                        OPTIND=$((OPTIND+1)); 
                       fi
                       eval "BRANCH=\${$((OPTIND))}"
                       if [[ $BRANCH =~ ^-.* ]] || [ -z ${BRANCH} ] ; then 
                        BRANCH="${DEFAULT_BRANCH}";
                       elif [ -n ${BRANCH} ]; then 
                        OPTIND=$((OPTIND+1)) 
                       fi 
                       do_code_clone=true ;;
        t | test )     SUITE="$OPTARG";
                       do_test_suite=true ;;
        n | newcase )  do_case_create=true ;;
        d | duration ) DURATION="$OPTARG" ;;
        c | config )   eval "CONFIG=\${$((OPTIND))}" # allows one optional argument
                       if [[ $CONFIG =~ ^-.* ]] || [ -z ${CONFIG} ] ; then
                        CONFIG="";
                       elif [ -n ${CONFIG} ]; then
                        OPTIND=$((OPTIND+1));
                       fi ;;
        e | emc )      do_energy_mass=true ;;
        k | kombo )    NAMEFILE="$OPTARG";
                       do_case_kombo=true;;
        b | build )    do_case_build=true ;;
        q | qsubmit )  do_case_queue=true ;;
        a | analyze )  eval "CASE2_NAME=\${$((OPTIND))}" # allows one optional argument
                       if [[ $CASE2_NAME =~ ^-.* ]] || [ -z ${CASE2_NAME} ] ; then
                        CASE2_NAME="";
                       elif [ -n ${CASE2_NAME} ]; then
                        OPTIND=$((OPTIND+1))
                       fi
		       do_case_analyze=true ;;
        h | help )     usage; exit ;;
        : )            printf "\n--- ERROR: Argument missing from -%s. Use -h for help.\n\n" \
                              ${OPTARG};
                       exit 1 ;; # getopt error
        ? )            printf "\n--- ERROR: Option -%s is incorrect. Use -h for help.\n\n" \
                              ${OPTARG}; 
                       exit 1 ;; # getopt error
      esac >&2
    done
    shift $((OPTIND-1)) # remove parsed options and args from $@ list

    # Check and set flags and variables dependent on command line options
    if [ -z $SANDBOX ] && [ ! -z $CONFIG ]; then

     if [ -d ${CODE_BASE} ]; then
      sandboxes=($(compgen -G "${CODE_BASE}/*"))
      printf "\n--- Sandbox(es) created in $(cd ${CODE_BASE} && dirs +0):\n\n"
      for ((i=0;i<${#sandboxes[@]};i++)); do
       printf "    $(basename ${sandboxes[$i]})\n"
      done
      printf "\n--- Please specify a sandbox after -s. Use -h to get help\n\n"
      exit
     else
      printf "\n--- ERROR: No sandboxes yet created. Use -h to get help\n\n"
      exit 1
     fi

    elif [ ! -z $SANDBOX ] && \
         [ "${do_code_clone,,}" != "true" ] && \
         [ "${do_test_suite,,}" != "true" ] && \
         [ "${do_case_create,,}" != "true" ] && \
         [ "${do_case_build,,}" != "true" ] && \
         [ "${do_case_queue,,}" != "true" ] && \
         [ "${do_case_kombo,,}" != "true" ] && \
         [ "${do_case_analyze,,}" != "true" ]  ; then

     sandboxcases=($(compgen -G "${CASE_BASE}/*${SANDBOX}*"))
     sandboxexist=($(compgen -G "${CODE_BASE}/${SANDBOX}"))

     if [ ! -z ${sandboxcases[0]} ] ; then

      printf "\n--- Cases so far created for '${SANDBOX}' sandbox in\n    ${CASE_BASE}:\n"
      for ((i=0;i<${#sandboxcases[@]};i++)); do
       printf "\n    $(basename ${sandboxcases[$i]})"
      done
      printf "\n\n--- Nothing else to do for sandbox '${SANDBOX}'. See -h for help.\n\n"

     elif [ ! -z ${sandboxexist[0]} ] ; then

      printf "\n\n--- No cases yet created for sandbox '${SANDBOX}'. Use -n to create.\n\n"

     else

      printf "\n--- Sandbox '${SANDBOX}' does not yet exist. Use -f to create.\n\n"

     fi

     exit 

    elif [ "${do_case_build,,}" == "true" ] && [ "${do_case_queue,,}" == "true" ]; then
     printf "\n--- ERROR: Need to build and queue in seperate steps.\n\n"
     exit 1
    elif [ "${do_case_create,,}" == "true" ] && [ "${do_case_queue,,}" == "true" ]; then
     printf "\n--- ERROR: Need to create new case and queue in seperate steps.\n\n"
     exit 1
    elif [ "${do_test_suite,,}" == "true" ] && [ "${do_case_kombo,,}" == "true" ]; then
     printf "\n--- ERROR: Test suites only work on sandbox defaults (not -k combinations)\n\n"
     exit 1
    elif [ "${do_test_suite,,}" == "true" ] && [ "${do_case_create,,}" == "true" ]; then
     printf "\n--- ERROR: Test suites create their own case (no -n required)\n\n"
     exit 1
    elif [ "${do_test_suite,,}" == "true" ] && [ "${do_case_build,,}" == "true" ]; then
     printf "\n--- ERROR: Test suites build on their own (no -b required)\n\n"
     exit 1
    elif [ "${do_test_suite,,}" == "true" ] && [ "${do_case_analyze,,}" == "true" ]; then
     printf "\n--- ERROR: Test suites already provide summary analysis (no -a required)\n\n"
     exit 1
    fi

    re='^[0-9]+$'
    if ! [[ $DURATION =~ $re ]] ; then 
     echo "--- ERROR: duration ${DURATION} is not an integer" >&2; 
    fi

    # Set compset, resolution, pe-layout, walltime and based on specific tests
    # This is only needed if no running an E3SM test suite
    if [ "${do_test_suite,,}" != "true" ]; then

     readonly MODEL_START_TYPE="initial"  # set to 'initial' or 'continue' only
     readonly START_DATE="0001-01-01"

     ####################################################################
     # SETTINGS FOR V2 MODEL CONFIGURATIONS. THESE WILL CHANGE TO V3 WHEN 
     # AVAILABLE. MAXCOMBINATIONS=0 INDICATES NOT YET READY TO RUN

     # CASE NAMELIST COMBINATION LIMITS - DO NOT CHANGE #################

     TOTALCOMBINATIONS[D]=32
     TOTALCOMBINATIONS[G]=8
     TOTALCOMBINATIONS[F]=8
     TOTALCOMBINATIONS[B]=4


     # D-CASES ##########################################################

     # D-CASE STANDARD RESOLUTION MESH
     CONFIGDESCRIPTION[D]="D-CASE AT STANDARD RESOLUTION WITH JRA 1.5 FORCING" 
     COMPSET[D]="2000_DATM%JRA-1p5_SLND_MPASSI_DOCN%SOM_DROF%JRA-1p5_SGLC_SWAV_TEST"
     RESOLUTION[D]="TL319_EC30to60E2r2"
     PELAYOUT[D]="S"
     MAXCOMBINATIONS[D]=32

     # D-CASE STANDARD RESOLUTION MESH WITH ICE SHELVES 
     CONFIGDESCRIPTION[DWISC]="D-CASE ON SORRM MESH WITH JRA 1.5 FORCING" 
     COMPSET[DWISC]="2000_DATM%JRA-1p5_SLND_MPASSI_DOCN%SOM_DROF%JRA-1p5_SGLC_SWAV_TEST"
     RESOLUTION[DWISC]="TL319_ECwISC30to60E2r1"
     PELAYOUT[DWISC]="S"
     MAXCOMBINATIONS[DWISC]=8

     # D-CASE OLD DTESTM
     CONFIGDESCRIPTION[DOLD]="D-CASE USING OLD DTESTM CONFIGURATION" 
     COMPSET[DOLD]="2000_DATM%NYF_SLND_MPASSI_DOCN%SOM_DROF%NYF_SGLC_SWAV_TEST"
     RESOLUTION[DOLD]="T62_EC30to60E2r2"
     PELAYOUT[DOLD]="S"
     MAXCOMBINATIONS[DOLD]=32

     # D-CASE COLUMN-WISE SEA ICE MODEL 
     CONFIGDESCRIPTION[DC]="D-CASE COLUMN TEST OF SEA ICE MODEL" 
     COMPSET[DC]="DTESTM-COL"
     RESOLUTION[DC]="T62_oQU480"
     PELAYOUT[DC]="S"
     MAXCOMBINATIONS[DC]=32

     # D-CASE WC14 MESH
     CONFIGDESCRIPTION[DWC14]="D-CASE ON WC14 MESH WITH JRA 1.5 FORCING" 
     COMPSET[DWC14]="2000_DATM%JRA-1p5_SLND_MPASSI_DOCN%SOM_DROF%JRA-1p5_SGLC_SWAV_TEST"
     RESOLUTION[DWC14]="TL319_WC14to60E2r3"
     PELAYOUT[DWC14]="S"
     MAXCOMBINATIONS[DWC14]=8

     # D-CASE SORRM MESH WITH ICE SHELVES 
     CONFIGDESCRIPTION[DSORRM]="D-CASE ON SORRM MESH WITH JRA 1.5 FORCING" 
     COMPSET[DSORRM]="2000_DATM%JRA-1p5_SLND_MPASSI_DOCN%SOM_DROF%JRA-1p5_SGLC_SWAV_TEST"
     RESOLUTION[DSORRM]="TL319_SOwISC12to60E2r4"
     PELAYOUT[DSORRM]="S"
     MAXCOMBINATIONS[DSORRM]=8

     # D-CASE WAVES (CURRENTLY UNDER CONSTRUCTION)
     CONFIGDESCRIPTION[DW]="D-CASE WITH WAVES ON STANDARD MESH WITH JRA 1.5 FORCING" 
     COMPSET[DW]="DTEST-JRA1p5-WW3"
     RESOLUTION[DW]="TL319_EC30to60E2r2_wQU225EC30to60E2r2"
     PELAYOUT[DW]="S"
     MAXCOMBINATIONS[DW]=0

     # D-CASE SEA ICE WITH BGC (CURRENTLY UNDER CONSTRUCTION)
     CONFIGDESCRIPTION[DBGC]="D-CASE WITH BGC ON STANDARD MESH WITH JRA 1.5 FORCING" 
     COMPSET[DBGC]="DTESTM-BGC"
     RESOLUTION[DBGC]="T62_EC30to60E2r2"
     PELAYOUT[DBGC]="S"
     MAXCOMBINATIONS[DBGC]=0


     # G-CASES ##########################################################

     # G-CASE ON STANDARD RESOLUTION 
     CONFIGDESCRIPTION[G]="G-CASE AT STANDARD RESOLUTION WITH JRA 1.5 FORCING" 
     COMPSET[G]="GMPAS-JRA1p5"
     RESOLUTION[G]="TL319_EC30to60E2r2"
     PELAYOUT[G]="L"
     MAXCOMBINATIONS[G]=8

     # G-CASE WC14 MESH WITH JRA1p4
     CONFIGDESCRIPTION[GWC14]="G-CASE ON WC14 MESH WITH JRA 1.5 FORCING" 
     COMPSET[GWC14]="GMPAS-JRA1p5" # (or GMPAS-JRA1p5)
     RESOLUTION[GWC14]="TL319_WC14to60E2r3"
     PELAYOUT[GWC14]="" 
     MAXCOMBINATIONS[GWC14]=0

     # G-CASE SORRM MESH WITH JRA1p4
     CONFIGDESCRIPTION[GSORRM]="G-CASE ON SORRM MESH WITH JRA 1.5 FORCING" 
     COMPSET[GSORRM]="GMPAS-JRA1p5-DIB-ISMF"
     RESOLUTION[GSORRM]="TL319_SOwISC12to60E2r4"
     PELAYOUT[GSORRM]="S"
     MAXCOMBINATIONS[GSORRM]=4

     # G-CASE STANDARD RESOLUTION WITH WAVES (UNDER CONSTRUCTION)
     CONFIGDESCRIPTION[GW]="G-CASE WITH WAVES ON STANDARD MESH WITH JRA 1.5 FORCING" 
     COMPSET[GW]="GMPAS-JRA1p5-WW3"
     RESOLUTION[GW]="TL319_EC30to60E2r2_wQU225EC30to60E2r2"
     PELAYOUT[GW]="L"
     MAXCOMBINATIONS[GW]=0

     # G-CASE WITH BGC 
     CONFIGDESCRIPTION[GBGC]="G-CASE WITH BGC ON STANDARD MESH WITH JRA 1.5 FORCING" 
     COMPSET[GBGC]="GMPAS-JRA1p5"
     RESOLUTION[GBGC]="TL319_EC30to60E2r2"
     PELAYOUT[GBGC]="L"
     MAXCOMBINATIONS[GBGC]=4


     # B-CASES ##########################################################

     # B-CASE STANDARD RESOLUTION
     CONFIGDESCRIPTION[B]="B-CASE ON STANDARD MESH" 
     COMPSET[B]="WCYCL1850"
     RESOLUTION[B]="ne30pg2_EC30to60E2r2"
     PELAYOUT[B]="L"
     MAXCOMBINATIONS[B]=4

     # B-CASE WITH WC14 MESH
     CONFIGDESCRIPTION[BWC14]="B-CASE WITH STANDARD RESOLUTION ATM/LND AND WC14 OCN/ICE" 
     COMPSET[BWC14]="WCYCL1850"
     RESOLUTION[BWC14]="ne30pg2_WC14to60E2r3"
     if [ "${MACHINES[${MACH}],,}" == "chrysalis" ]; then
      PELAYOUT[BWC14]="M"
      MAXCOMBINATIONS[BWC14]=2
     else
      MAXCOMBINATIONS[BWC14]=0
     fi

     # B-CASE NARRM CONFIGURATION (UNDER CONSTRUCTION)
     CONFIGDESCRIPTION[BNARRM]="B-CASE NORTH AMERICAN REGIONALLY REFINED MODEL" 
     COMPSET[BNARRM]="WCYCL1850"
     RESOLUTION[BNARRM]="northamericax4v1pg2_WC14to60E2r3"
     if [ "${MACHINES[${MACH}],,}" == "chrysalis" ]; then
      PELAYOUT[BNARRM]="M"
      MAXCOMBINATIONS[BNARRM]=2
     else
      MAXCOMBINATIONS[BNARRM]=0
     fi

     # B-CASE SORRM WITH STANDARD RESOLUTION ATMOSPHERE
     CONFIGDESCRIPTION[BSORRM]="B-CASE WITH ICE SHELVES AND REFINED SOUTHERN OCEAN" 
     COMPSET[BSORRM]="CRYO1850"
     RESOLUTION[BSORRM]="ne30pg2_SOwISC12to60E2r4"
     if [ "${MACHINES[${MACH}],,}" == "chrysalis" ]; then
      PELAYOUT[BSORRM]="M"
      MAXCOMBINATIONS[BSORRM]=2
     else
      PELAYOUT[BSORRM]="L"
      MAXCOMBINATIONS[BSORRM]=2
     fi

     # B-CASE STANDARD RESOLUTION WITH WAVES
     CONFIGDESCRIPTION[BW]="B-CASE WITH WAVES AT STANDARD RESOLUTION" 
     COMPSET[BW]="WCYCL1850-WW3"
     RESOLUTION[BW]="ne30pg2_EC30to60E2r2_wQU225EC30to60E2r2"
     PELAYOUT[BW]="L"
     MAXCOMBINATIONS[BW]=0

     # B-CASE STANDARD RESOLUTION WITH BGC 
     CONFIGDESCRIPTION[BBGC]="B-CASE WITH BGC AT STANDARD RESOLUTION" 
     COMPSET[BBGC]="BGCEXP_CNTL_CNPECACNT_1850"
     RESOLUTION[BBGC]="ne30pg2_r05_EC30to60E2r2"
     PELAYOUT[BBGC]="L"
     MAXCOMBINATIONS[BBGC]=2

     #
     ####################################################################

     # print entire configuration options of -c is specified without an argument
     if [ -z $CONFIG ]; then

      # sort keys alphabetically
      local sortedkey=($(echo ${!MAXCOMBINATIONS[@]}| tr " " "\n" | sort -n))

      printf "\n--- Available ${MACHINES[${MACH}]} configs "
      printf "(compset, res, layout, max namelist combos):\n"
      for key in "${sortedkey[@]}"; do 
       if [[ ${MAXCOMBINATIONS[$key]} > 0 ]]; then
        if [[ ${DEFAULT_CONFIG} == ${key} ]]; then
         printf "\n    $key => \t\e[1;34mDEFAULT\e[0m ${CONFIGDESCRIPTION[$key]}\n"
        else
         printf "\n    $key => \t${CONFIGDESCRIPTION[$key]}\n"
        fi
        printf '\t%s\n' "        ${COMPSET[$key]}, " 
        printf "            \t${RESOLUTION[$key]}, "
        printf "${PELAYOUT[$key]}, "; 
        printf "${MAXCOMBINATIONS[$key]}\n"; 
       fi
      done

      printf "\n--- Draft configs not yet available on ${MACHINES[${MACH}]}:\n"
      for key in "${!MAXCOMBINATIONS[@]}"; do 
       if [[ ${MAXCOMBINATIONS[$key]} == 0 ]]; then
        printf "\n    $key => \t${CONFIGDESCRIPTION[$key]}"
       fi
      done 
      echo $'\n'
      exit

     # otherwise just provide configuration options, with extra checks
     elif [[ ! "${!MAXCOMBINATIONS[@]}" =~ "${CONFIG}" ]]; then

      printf "\n--- ERROR: \x1B[31mConfiguration '${CONFIG}' not found\e[0m\n\n" 
      exit

     else

      printf "\n--- \x1B[34m${CONFIGDESCRIPTION[$CONFIG]}\e[0m\n" 
      printf '%s\n' "    ${COMPSET[$CONFIG]}, " 
      printf "    ${RESOLUTION[$CONFIG]}, "
      printf "${PELAYOUT[$CONFIG]} layout, "; 
      printf "${MAXCOMBINATIONS[$CONFIG]} combos allowed\n"; 

      # UNDER NO CIRCUMSTANCE DISABLE THIS BLOCK OF CODE >
      readonly CASETYPE=$(echo ${CONFIG} | head -c 1)
      for key in "${!TOTALCOMBINATIONS[@]}"; do
       if [[ "${key}" == "${CASETYPE}" ]] &&
          [[ "${MAXCOMBINATIONS[$CONFIG]}" -gt "${TOTALCOMBINATIONS[$CASETYPE]}" ]] ; then
        printf "\n--- ERROR: Maximum combinations for a ${CASETYPE}-case exceeded\n\n"
        exit
       fi
      done
      # END OF NO CHANGE BLOCK

     fi

     readonly HOURS=$(((DURATION+2-1)/2)) # walltime is min(ceil(months/2),12) hours
     printf -v WALLTIME "%2.2i:00:00" $((HOURS<12 ? HOURS : 12))

     readonly STOP_OPTION="nmonths"
     readonly STOP_N="$((DURATION<12 ? DURATION : 12))"
     readonly REST_OPTION="${STOP_OPTION}"
     readonly REST_N="${STOP_N}"
     readonly HIST_OPTION="${STOP_OPTION}"
     readonly HIST_N="${STOP_N}"
     readonly RESUBMIT="$(((DURATION+12-1)/12-1))"
     readonly DO_SHORT_TERM_ARCHIVING=false


     # process namefile and prepare for the casename
     if [ "${do_case_kombo,,}" == "true" ]; then
      if [[ "${NAMEFILE}" =~ ".nlk" ]]; then
       readonly KOMBOTAG=`echo ${NAMEFILE} | cut -d . -f 1`
      else
       printf "\n--- ERROR: Specify a <kombotag>.nlk namelist file with -k\n\n"
       exit 1
      fi
     fi

    fi # if [ "${do_test_suite,,}" != "true" ]

    # code directories
    readonly CODE_ROOT="${CODE_BASE}/${SANDBOX}"

    # check that the cloned code has same fork as ${FORK} and uses ${BRANCH}
    # and post warning if not. This situation could arise if the code has
    # aleady been cloned from a different fork and branch than the default
    # or from the original specification when -f was last used.  The value
    # of FORK and BRANCH is reset accordingly so that the provenance records
    # the correct code being used in the sandbox.

    if [ "${do_code_clone,,}" != "true" ] && [ -d ${CODE_ROOT} ] ; then 

     pushd ${CODE_ROOT} 

     # check origin url
     local checkurl=($(git config --get remote.origin.url))
     if [[ ! "${checkurl}" == "${FORK}" ]]; then
      printf "\n--- Fork origin is different from the default: ${FORK}"
      printf "\n    Sandbox origin is: \x1B[34m${checkurl}\e[0m\n"
      FORK=${checkurl}        
     fi

     # check local branch
     local checkbranch=($(git rev-parse --abbrev-ref HEAD))
     if [[ ! "${checkbranch}" == "${BRANCH}" ]]; then
      printf "\n--- Branch being used is not the default: ${BRANCH}"
      printf "\n    Current branch is: \x1B[34m${checkbranch}\e[0m\n"
      BRANCH=${checkbranch}
     fi

     # get code branch hash for case naming
     HASH=($(git rev-parse HEAD))

     popd

    fi

    # obtain remote
    SREMOTE=($(echo ${FORK} | cut -d / -f 1 | cut -d : -f 2))
    [ ! -z ${SREMOTE} ] && SREMOTE=".${SREMOTE}"

    # create short branch name for case name
    SBRANCH=($(echo ${BRANCH} | rev | cut -d '/' -f 1 | rev))
    [ ! -z ${SBRANCH} ] && SBRANCH=".${SBRANCH}"

    # set case name and directories
    if [ "${do_test_suite,,}" == "true" ]; then
      local CMODIFIER="Test.${SUITE}"
    else
      local CMODIFIER="${CONFIG}$((DURATION<12 ? DURATION : 12))${TEST}"
      [ "${do_case_kombo,,}" == "true" ] && CMODIFIER="${CMODIFIER}.${KOMBOTAG}"
      [ "${do_energy_mass,,}" == "true" ] && CMODIFIER="${CMODIFIER}.emc"
    fi
    readonly CASE_NAME="${CMODIFIER}.${SANDBOX}${SBRANCH}${SREMOTE}.${MACHINES[${MACH}]}"
    readonly CASE_ROOT="${CASE_BASE}/${CASE_NAME}"
    readonly CASE_BUILD_DIR=${CASE_ROOT}/build
    readonly CASE_ARCHIVE_DIR=${CASE_ROOT}/archive
    readonly CASE_SCRIPTS_DIR=${CASE_ROOT}/case_scripts
    readonly CASE_RUN_DIR=${CASE_ROOT}/run

    # analysis checks and assignments
    if [ "${do_case_analyze,,}" == "true" ] ; then
      if [ ! -z ${CASE2_NAME} ] && [ "${CASE2_NAME}" != "${CASE_NAME}" ] ; then
       readonly CASE2_ROOT="${CASE_BASE}/${CASE2_NAME}"
       if [ ! -d ${CASE2_ROOT} ]; then
        if [ -d "${CODE_BASE}/${CASE2_NAME}" ]; then
         printf "\n--- ERROR: '${CASE2_NAME}' is a sandbox, not a case\n\n"
        else
         printf "\n--- ERROR: Case '${CASE2_NAME}' not in \n    ${CASE_BASE}\n\n" 
        fi
        exit 1
       fi
      elif [ "${CASE2_NAME}" == "${CASE_NAME}" ] && (( combinations == 1 )); then
       printf "\n--- ERROR: The comparison case and -s sandbox case are identical ---\n\n" 
       exit 1
      elif [ -z ${CASE2_NAME} ] && (( combinations == 1 )); then 
       printf "\n--- ERROR: No comparison case given with -a. Use -h for help ---\n\n" 
       exit 1
      else
       readonly CASE2_ROOT=""
      fi
    fi

}

#---------------------------------
copy_script() {

    local THIS_SCRIPT_NAME=`basename $0`
    local THIS_SCRIPT_DIR=`dirname $0`
    readonly SCRIPT_PROVENANCE_NAME=${THIS_SCRIPT_NAME}.`date +%Y%m%d-%H%M%S`
    local SCRIPT_PROVENANCE_BASE="${HOME}/${CASE_GROUP}/provenance"
    if [ -z ${CASE_NAME} ]; then
     readonly SCRIPT_PROVENANCE_DIR="${SCRIPT_PROVENANCE_BASE}/${SANDBOX}"
    else
     readonly SCRIPT_PROVENANCE_DIR="${SCRIPT_PROVENANCE_BASE}/${SANDBOX}/${CASE_NAME}"
    fi

    if [ ! -d ${SCRIPT_PROVENANCE_DIR} ]; then
     mkdir -p ${SCRIPT_PROVENANCE_DIR}
    fi

    # improve readability of provenance
    if [ "${do_energy_mass,,}" == "true" ]; then
      local energy_mass_conservation="ON"
    else
      local energy_mass_conservation="OFF"
    fi

    cat << EOF > ${SCRIPT_PROVENANCE_DIR}/${SCRIPT_PROVENANCE_NAME}
#!${SHELL}
#######################################################################
#
# Specifications for ${SCRIPT_PROVENANCE_NAME}
#
# ${THIS_SCRIPT_DIR}/${THIS_SCRIPT_NAME}
#
# Command:    ${THIS_SCRIPT_NAME} $*
# Origin:     ${FORK}
# Branch:     ${BRANCH}
# Hash:       ${HASH}
# Sandbox:    ${CODE_ROOT}
# Case:	      ${CASE_NAME}
# Casedir:    ${CASE_ROOT}
#
EOF

 if [ "${do_test_suite,,}" == "true" ]; then

    cat << EOF >> ${SCRIPT_PROVENANCE_DIR}/${SCRIPT_PROVENANCE_NAME}
# Testsuite:  e3sm_${SUITE}_developer
# Machine:    ${MACHINES[${MACH}]}
#
#######################################################################
EOF

 else

    cat << EOF >> ${SCRIPT_PROVENANCE_DIR}/${SCRIPT_PROVENANCE_NAME}
# Start Type: ${MODEL_START_TYPE}
# Start Date: ${START_DATE}
# Simulation: ${CONFIGDESCRIPTION[$CONFIG]}
# Compset:    ${COMPSET[$CONFIG]} 
# Resolution: ${RESOLUTION[$CONFIG]}
# PE Layout:  ${PELAYOUT[$CONFIG]}
# Walltime:   ${WALLTIME}
# Duration:   ${STOP_N} months
# Resubmit:   ${RESUBMIT} count
# Machine:    ${MACHINES[${MACH}]}
#
# Conservation checking is ${energy_mass_conservation}
#
#######################################################################
EOF

 fi

    # add additional information if namelist combinations are included
    if [ "${do_case_kombo,,}" == "true" ] & [ -f "${NAMEFILE}" ]; then
     echo "# " >> ${SCRIPT_PROVENANCE_DIR}/${SCRIPT_PROVENANCE_NAME}
     echo "# Sea ice namelist combinations summary for ${NAMEFILE}:" >> \
             ${SCRIPT_PROVENANCE_DIR}/${SCRIPT_PROVENANCE_NAME}
     echo "# " >> ${SCRIPT_PROVENANCE_DIR}/${SCRIPT_PROVENANCE_NAME}
     while IFS= read -r line; do
      nmline=`echo "$line" | cut -d = -f 1 | rev | cut -c 2- | rev`
      if [ ! -z $nmline ]; then # remove blank lines
       echo "# ${line}" >> ${SCRIPT_PROVENANCE_DIR}/${SCRIPT_PROVENANCE_NAME}
      fi
     done < "$NAMEFILE"
     echo "# " >> ${SCRIPT_PROVENANCE_DIR}/${SCRIPT_PROVENANCE_NAME}
     echo "#######################################################################" >> \
           ${SCRIPT_PROVENANCE_DIR}/${SCRIPT_PROVENANCE_NAME}
    fi

    # add entire script beneath this
    cat ${THIS_SCRIPT_DIR}/${THIS_SCRIPT_NAME} | sed '1d' >> \
        ${SCRIPT_PROVENANCE_DIR}/${SCRIPT_PROVENANCE_NAME} 

    # provide an ASCII call tree in the base provenance directory
    pushd ${SCRIPT_PROVENANCE_BASE}
    find . -print | sed -e 's;/*/;|;g;s;|; |;g' >| "`echo ${THIS_SCRIPT_NAME} | cut -f 1 -d .`"
    popd

}

#---------------------------------
clone_code() {

    if [ "${do_code_clone,,}" != "true" ]; then
        echo $'\n--- Skipping clone of code ---'
        return
    elif [ -e ${CODE_ROOT} ]; then
        while true; do
         read -p $'\n--- Sandbox already exists. Overwrite (Y/N)?' yn
         case $yn in
          [Yy]* ) local stop_clone=false; break;;
          [Nn]* ) local stop_clone=true; break;;
          * ) printf "    Please answer yes or no\n";;
         esac
        done

        if [ "${stop_clone,,}" == "true" ]; then
         echo $'\n--- Keeping existing clone ---\n'
         exit 
        else
         rm -rf ${CODE_ROOT}
        fi
    fi

    echo $'\n--- Cloning code ---'

    echo $'\n---' "Fork: ${FORK}"

    echo $'\n---' "Branch: ${BRANCH}"

    echo $'\n'

    mkdir -p ${CODE_ROOT}

    pushd ${CODE_ROOT}

    # This will put repository, with all code
    git clone ${FORK} . 

    if [ $? != 0 ]; then
     printf "\n--- ERROR: Problem cloning ${FORK}\n\n"
     exit 1
    fi
    
    # Check out desired branch
    git checkout ${BRANCH}

    if [ $? != 0 ]; then
     printf "\n--- ERROR: Cannot checkout ${BRANCH}\n\n"
     exit 1
    fi

    # Bring in all submodule components
    git submodule update --init --recursive

    if [ $? != 0 ]; then
     printf "\n--- ERROR: Submodule update failed\n\n"
     exit 1
    fi
    
    popd
}

#---------------------------------
e3sm_testsuites() {

    # Run standard e3sm_developer test suite on sandbox and compare with e3sm master

    if [ "${do_test_suite,,}" != "true" ]; then
        printf "\n--- Skipping E3SM test suite ---\n"
        return
    else
        printf "\n--- Starting E3SM testing on current sandbox branch ---\n"
    fi


    # set E3SM test suites to be used, check they exist, then run them

    for ((i=0;i<${#TESTSUITES[@]};i++)); do

     if [[ "e3sm_${SUITE}_developer" == ${TESTSUITES[i]} ]]; then

      if grep -q ${TESTSUITES[i]} ${CODE_ROOT}/cime_config/tests.py ; then

       # run test script
       printf "\n--- Working on ${TESTSUITES[i]} testsuite ---\n\n"
       ${CODE_ROOT}/cime/scripts/create_test \
         --output-root "${CASE_ROOT}" \
         --machine "${MACHINES[${MACH}]}" \
         --project "${PROJECTS[${MACH}]}" \
         --wait "${TESTSUITES[i]}"

       if [ $? != 0 ]; then
        local testpassfail=fail
       else
        local testpassfail=pass
       fi
  
       # run last available cs_status file
       if compgen -G ${CASE_ROOT}/cs.status.* > /dev/null; then
        printf "\n\n"
        CS_STATUS_FILES=($(ls ${CASE_ROOT}/cs.status.*))
        ${CS_STATUS_FILES[-1]}
       else
        printf "\n--- ERROR: Cannot find CS_STATUS executable for '${TESTSUITES[i]}'\n\n"
        exit 2
       fi

       if [ "${testpassfail,,}" == "fail" ]; then
        printf "\n--- Testsuite '${TESTSUITES[i]}' FAILED for sandbox '${SANDBOX}'\n"
       else
        printf "\n--- Testsuite '${TESTSUITES[i]}' PASSED for sandbox '${SANDBOX}'\n"
       fi

      else

       printf "\n--- ERROR: '${TESTSUITES[i]}' testsuite not found in sandbox '${SANDBOX}'\n\n"
       exit 2

      fi 

     fi 

    done

    if [ ! -d ${CASE_ROOT} ]; then

     printf "\n--- ERROR: 'e3sm_${SUITE}_developer' not available as a testsuite\n\n"
     exit 3

    fi

}

#---------------------------------
namelist_kombo() {

    # initialize 2D arrays of all the namelist options to be read in and
    # and the final namelists to be constructed and made globally available. 

    declare -A option
    declare -A modelt
    declare -gA namelist
    declare -gA model

    if [ "${do_case_kombo,,}" != "true" ]; then
        printf "\n--- Skipping namelist combinations ---\n"
        combinations=1
        komboname[0]="k000"
        return
    elif [ -f "${CASE_ROOT}/${NAMEFILE}" ] ; then
        printf "\n--- Reading combinations from ${CASE_NAME}/${NAMEFILE}\n"
        NAMEFILE="${CASE_ROOT}/${NAMEFILE}"
    elif [ ! -f "${NAMEFILE}" ]; then
        printf "\n--- ERROR: '${NAMEFILE}' not here in $(cd `pwd` && dirs +0),"
        printf "\n            nor in case ${CASE_NAME}\n\n"
        exit 1
    else
        printf "\n--- Reading combinations from ${NAMEFILE}\n"
    fi


    # Read in the namelist combinations from the file. These must appear in the 
    # format of one line for one namelist option, with comma-seperated values
    # to be tested, surrounded by brackets, as in this example:
    #
    # [mpassi]
    # NAMELIST_ENTRY1 = {false,true}
    # NAMELIST_ENTRY2 = {'string1','string2','string3'}
    # NAMELIST_ENTRY3 = {number1,number2}
    #
    # Note that a single namelist entry can be provided, e.g.:
    #
    # NAMELIST_ENTRY4 = 'evp'
    #
    # The namelist combinations file is read in line-by-line. Blank lines are 
    # discarded. Three main things are recorded: 1) namechangenumber, which is 
    # the total number of namelist entries to change; 2) numberofoptions, which
    # is the total number of comma-separated entries to be tested for a given
    # namelist options; 3) the list of options in an array with dimensions
    # option[namechangenumber,numberofoptions]; and, 4) the model namelist to 
    # which the changes are to be added.  In the example above, they are to be
    # added to user_nl_mpassi. 

    modeltype="";
    local i=0
    while IFS= read -r line; do
     namechange[$i]=`echo "$line" | cut -d = -f 1 | rev | cut -c 2- | rev`  
     if [ ! -z ${namechange[$i]} ]; then # remove blank lines
      case $line in 
         "["*"]") # where square brackets are provided, namelist is specified 
           modeltype=`echo "$line" | cut -d [ -f 2 | cut -d ] -f 1`
           advancei=0
           ;;
         *"{"*"}") # where brackets are provided, multiple options provided
           options=`echo "$line" | cut -d { -f 2 | cut -d } -f 1`
           commas=${options//[^,]}
           ((numberofoptions[$i]=${#commas}+1))
           if [ -z ${modeltype} ]; then
            printf "\n--- ERROR: No model specified for '${namechange[$i]}'\n\n"
            exit 1
           fi
           for ((j=0;j<${numberofoptions[$i]};j++)); do
            ((jj=j+1))
            option[$i,$j]="${namechange[$i]} = `echo ${options} | cut -d , -f $jj`"
            modelt[$i,$j]=${modeltype}
           done
           advancei=1
           ;;
         *) # where only one option is provided
           numberofoptions[$i]=1
           if [ -z ${modeltype} ]; then
            printf "\n--- ERROR: No model specified for '${namechange[$i]}'\n\n"
            exit 1
           fi
           option[$i,0]=${line}
           modelt[$i,0]=${modeltype}
           advancei=1
           ;;
      esac
      ((i=i+advancei))
     fi
    done < "$NAMEFILE"
    namechangenumber=$i

    # now calculate the total number of combinations based on the nuber of
    # options to be tested, and limit according the case type.  

    combinations=1
    for ((m=0;m<${namechangenumber};m++)); do
     index[m]=0
     combinations=$((${combinations}*${numberofoptions[m]}));
    done
    if (( combinations > ${MAXCOMBINATIONS[$CONFIG]} )) ; then
     printf "\n--- ERROR: Permutations exceed maximum for config: ${MAXCOMBINATIONS}\n\n"
     exit 1
    fi

    # now construct the individual namelists to be used for each combination
    # and write them to an output array, namelist[combinations,namechangenumber]
    # with the combination labels in komboname[combinations]

    for ((i=0;i<${combinations};i++)); do
     printf -v komboname[${i}] "k%3.3i" ${i}
     printf $'\n'"    \x1B[34m${komboname[${i}]}\e[0m ->\n" 
     for ((m=$((${namechangenumber}-1));m>=0;m--)); do
      namelist[${i},${m}]=${option[${m},${index[${m}]}]}
      model[${i},${m}]=${modelt[${m},${index[${m}]}]}
      printf "    ${model[${i},${m}]}: \t${namelist[${i},${m}]}\n"
      if [ $m == $((${namechangenumber}-1)) ]; then
       ((index[m]=${index[m]}+1)) 
      elif [ ${index[m+1]} == ${numberofoptions[m+1]} ]; then
       ((index[m]=${index[m]}+1))       
       index[m+1]=0
      fi
      if [[ $i < $((${combinations}-1)) ]] && [[ ${index[0]} == ${numberofoptions[0]} ]] 
      then 
       printf "\n---- INTERNAL ERROR: Permutations incorrect: $i $((${combinations}-1))\n\n"
       exit 125
      fi
     done
    done

}

#---------------------------------
create_newcase() {

    if [ "${do_case_create,,}" != "true" ]; then
        echo $'\n--- Skipping create newcase ---'
        return
    elif [ ! -d ${CODE_ROOT} ]; then
        printf "\n--- ERROR: No code sandbox from which to create a new case ---\n\n"
        exit 1
    elif [ -e ${CASE_ROOT} ]; then
        while true; do
         read -p $'\n--- Case root already exists. Overwrite (Y/N)?' yn
         case $yn in
          [Yy]* ) local stop_script=false; break;;
          [Nn]* ) local stop_script=true; break;;
          * ) printf "    Please answer yes or no\n";;
         esac
        done
        if [ "${stop_script,,}" == "true" ]; then
         printf "\n--- Stopping the script ---\n\n"
         exit 
        else
         rm -rf ${CASE_ROOT}
        fi
    fi

    echo $'\n--- Starting to create a new case ---'
    echo $'\n'

    # set up case
    for ((i=0;i<${combinations};i++)); do

      ${CODE_ROOT}/cime/scripts/create_newcase \
        --case ${CASE_NAME} \
        --case-group ${CASE_GROUP} \
        --output-root ${CASE_ROOT} \
        --script-root "${CASE_SCRIPTS_DIR}.${komboname[${i}]}" \
        --handle-preexisting-dirs u \
        --compset "${COMPSET[$CONFIG]}" \
        --res "${RESOLUTION[$CONFIG]}" \
        --machine "${MACHINES[${MACH}]}" \
        --project "${PROJECTS[${MACH}]}" \
        --walltime ${WALLTIME} \
        --pecount "${PELAYOUT[$CONFIG]}"

      if [ $? != 0 ]; then
       printf "\n--- ERROR: Problem creating new case for ${komboname[${i}]}\n\n"
       exit 1
      fi

    done

    # copy combinations file to case directory and check namelist files exist
    if [ "${do_case_kombo,,}" == "true" ]; then
     cp ${NAMEFILE} ${CASE_ROOT}
     chmod 644 ${CASE_ROOT}/${NAMEFILE}
    fi

    # create branch information file
    printf "Case:       ${CASE_NAME}\n" >| ${CASE_ROOT}/case_info
    printf "Origin:     ${FORK}\n" >> ${CASE_ROOT}/case_info
    printf "Branch:     ${BRANCH}\n" >> ${CASE_ROOT}/case_info
    printf "Hash:       ${HASH}\n" >> ${CASE_ROOT}/case_info
    printf "Sandbox:    ${CODE_ROOT}\n" >> ${CASE_ROOT}/case_info
    printf '%s\n' "Compset:    ${COMPSET[$CONFIG]}" >> ${CASE_ROOT}/case_info
    printf "Resolution: ${RESOLUTION[$CONFIG]}\n" >> ${CASE_ROOT}/case_info
    printf "PE Layout:  ${PELAYOUT[$CONFIG]}" >> ${CASE_ROOT}/case_info

}

#---------------------------------
case_build() {

    for ((i=0;i<${combinations};i++)); do

     if [ -d ${CASE_SCRIPTS_DIR}.${komboname[${i}]} ]; then
      pushd ${CASE_SCRIPTS_DIR}.${komboname[${i}]}
     else
      echo $'\n--- Case does not exist to build ---'
      return
     fi

     # do_case_build = false
     if [ "${do_case_build,,}" != "true" ]; then

      # Use previously built executable, make sure it exists
      if [ ! -x ${CASE_BUILD_DIR}/e3sm.exe ]; then
       if [ $i == 0 ]; then
        if [ "${do_case_queue,,}" == "true" ]; then
         do_case_queue=false
         exit_script=true
         printf "\n--- No executable yet for ${SANDBOX}, queue cancelled ---\n" 
        else
         printf "\n--- No executable yet for ${SANDBOX}, See -h for help ---\n" 
        fi
       fi
      else
        ./xmlchange BUILD_COMPLETE=TRUE > /dev/null
        echo $'\n--- Skipping build for' "${komboname[${i}]}" 
      fi

     # do_case_build = true
     elif [ "${do_case_build,,}" == "true" ]; then

      printf "\n--- Building ${komboname[${i}]} ---\n\n"

      # Setup some CIME directories
      ./xmlchange EXEROOT=${CASE_BUILD_DIR} > /dev/null
      ./xmlchange RUNDIR="${CASE_RUN_DIR}.${komboname[${i}]}" > /dev/null

      # Short term archiving
      ./xmlchange DOUT_S=${DO_SHORT_TERM_ARCHIVING^^}
      ./xmlchange DOUT_S_ROOT="${CASE_ARCHIVE_DIR}.${komboname[${i}]}"

      # Build with COSP, except for a data atmosphere (datm)
      if [ `./xmlquery --value COMP_ATM` != "datm"  ]; then 
       echo $'\nConfiguring E3SM to use the COSP simulator\n'
       ./xmlchange --id CAM_CONFIG_OPTS --append --val='-cosp'
      fi

      # Extracts input_data_dir in case it is needed for user edits to the namelist later
      local input_data_dir=`./xmlquery DIN_LOC_ROOT --value`

      # clear namelist alterations
      if compgen -G "user_nl_mpassi" > /dev/null; then
       rm user_nl_mpassi
      fi
      if compgen -G "user_nl_mpaso" > /dev/null; then
       rm user_nl_mpaso
      fi

      # Finally, run CIME case.setup
      ./case.setup --reset

      # Switch on conservation analysis members if -e option given
      if [ "${do_energy_mass,,}" == "true" ]; then
       echo $'\n--- Setting up mass and energy conservation analysis ---\n'
       if [ -f "user_nl_mpaso" ]; then
        echo "config_am_conservationcheck_enable = true" >> user_nl_mpaso
       fi
       if [ -f "user_nl_mpassi" ]; then
        echo "config_am_conservationcheck_enable = true" >> user_nl_mpassi
       fi
      fi

      # If specifying namelist settings, add them to user_nl_mpassi here
      if [ "${do_case_kombo,,}" == "true" ]; then
       for ((m=0;m<${namechangenumber};m++)); do
        if [ -f "user_nl_${model[${i},${m}]}" ]; then
         echo ${namelist[${i},${m}]} >> "user_nl_${model[${i},${m}]}"
        else
         printf "\n--- ERROR: Namelist for ${model[${i},${m}]} not used in ${CONFIG} config.\n\n"
         exit 1
        fi
       done
      fi

      # Turn on debug compilation option if requested
      if [ "${DEBUG_COMPILE^^}" == "TRUE" ]; then
        ./xmlchange DEBUG=${DEBUG_COMPILE^^}
      fi

      # Remove any existing output files so they aren't erroneously used
      pushd ${CASE_RUN_DIR}.${komboname[${i}]}
      for x in `ls ${CASE_RUN_DIR}.${komboname[${i}]}`; do
       if [[ ${x} =~ ${CASE_NAME} ]]; then
        rm ${x}
       fi
      done
      popd

      # Run CIME case.build, which only needs to be done for the first combination
      if [ ${i} == 0 ]; then
       ./case.build
      fi

      # Call preview_namelists to make sure *_in and user_nl files are consistent.
      ./preview_namelists

      # check for errors
      if [ $? != 0 ]; then
       printf "\n--- ERROR: Problems with build for ${komboname[${i}]}\n\n"
       exit 1
      fi

     else

      printf "\n--- INTERNAL ERROR: Build option nondescript\n\n"
      exit 1

     fi

     popd

    done

}

#---------------------------------
case_queue() {

    for ((i=0;i<${combinations};i++)); do

     if [ ! -d "${CASE_SCRIPTS_DIR}.${komboname[${i}]}" ]; then

      echo $'\n--- Case does not exist to queue ---'

     elif [ "${do_case_queue,,}" != "true" ]; then

      echo $'\n--- Skipping queue for' "${komboname[${i}]}"

     else

      echo $'\n--- Starting queue for' "${komboname[${i}]}"
      pushd "${CASE_SCRIPTS_DIR}.${komboname[${i}]}"

      # Set simulation start date
      ./xmlchange RUN_STARTDATE=${START_DATE}

      # Segment length
      ./xmlchange STOP_OPTION=${STOP_OPTION,,},STOP_N=${STOP_N}

      # Restart frequency
      ./xmlchange REST_OPTION=${REST_OPTION,,},REST_N=${REST_N}

      # Coupler history
      ./xmlchange HIST_OPTION=${HIST_OPTION,,},HIST_N=${HIST_N}

      # Coupler budgets (always on)
      ./xmlchange BUDGETS=TRUE

      # Set resubmissions
      if (( RESUBMIT > 0 )); then
        ./xmlchange RESUBMIT=${RESUBMIT}
      fi

      # Run type
      # Start from default of user-specified initial conditions
      if [ "${MODEL_START_TYPE,,}" == "initial" ]; then
        ./xmlchange RUN_TYPE="startup"
        ./xmlchange CONTINUE_RUN="FALSE"

      # Continue existing run
      elif [ "${MODEL_START_TYPE,,}" == "continue" ]; then
        ./xmlchange CONTINUE_RUN="TRUE"

      elif [ "${MODEL_START_TYPE,,}" == "branch" ] || \
           [ "${MODEL_START_TYPE,,}" == "hybrid" ]; then

       ./xmlchange RUN_TYPE=${MODEL_START_TYPE,,}
       ./xmlchange GET_REFCASE=${GET_REFCASE}
       ./xmlchange RUN_REFDIR=${RUN_REFDIR}
       ./xmlchange RUN_REFCASE=${RUN_REFCASE}
       ./xmlchange RUN_REFDATE=${RUN_REFDATE}
       echo 'Warning: $MODEL_START_TYPE = '${MODEL_START_TYPE} 
       echo '$RUN_REFDIR = '${RUN_REFDIR}
       echo '$RUN_REFCASE = '${RUN_REFCASE}
       echo '$RUN_REFDATE = '${START_DATE}
 
      else
       printf "\n--- ERROR: '${MODEL_START_TYPE}' is unrecognized. Exiting.\n\n"
       exit 1
      fi

      if [ "${MACHINES[${MACH}],,}" == "chrysalis" ]; then
         ./xmlchange JOB_QUEUE="compute"
      fi

      if [ $? != 0 ]; then
       printf "\n--- ERROR: Problem getting ready to submit for ${komboname[${i}]}\n\n"
       exit 1
      fi

      # Run CIME case.submit
      ./case.submit

      if [ $? != 0 ]; then
       printf "\n--- ERROR: Problem submitting for ${komboname[${i}]}\n\n"
       exit 1
      fi

      popd

      [ "${MACHINES[${MACH}],,}" == "chrysalis" ] && squeue -u ${USER}
      [ "${MACHINES[${MACH}],,}" == "anvil" ] && squeue -u ${USER}

     fi

    done
 
}

#---------------------------------
case_analyze() {

    if [ "${do_case_analyze,,}" != "true" ]; then
     echo $'\n--- Skipping analysis ---'
     return
    else
     # set up provenance
     echo $'\n--- Starting analysis ---'
     local pscript="${SCRIPT_PROVENANCE_DIR}/${SCRIPT_PROVENANCE_NAME}"
     echo ": '" >> ${pscript}
    fi

    pushd ${CASE_ROOT}

    # grab month and year
    printf -v MONTH "%2.2i" $((DURATION<12 ? DURATION : 12))
    printf -v YEAR "%4.4i" $((1+RESUBMIT))

    # set conservation files to analyze individually
    EMCVARS+=( "EnergyError" "MassError" "SaltError" "CarbonError" )
    EMCVARMODIFY="relative"
    # from ridgepack_E3SM_sea_ice_conservation_am
    EVARS+=( "netEnergyFlux" 
             "energyConsFreezingPotential" 
             "energyConsSensibleHeatFlux" 
             "energyConsOceanHeatFlux" 
             "energyConsLongwaveUp" 
             "energyConsLongwaveDown" 
             "energyConsAbsorbedShortwaveFlux" 
             "energyConsOceanShortwaveFlux" 
             "energyConsAbsorbedShortwaveFlux" 
             "energyConsOceanShortwaveFlux" 
             "energyConsLatentHea" 
             "energyConsSnowfallHeat" )

    EMODIFY="energy"
    FILE_EXTE_EMC="mpassi.hist.am.conservationCheck.${YEAR}.nc"

    # set conservation files to analyze individually

    # set fields and month to analyze
    local AVE_FREQ="Monthly"
    FIELDS+=( "iceAreaCell" "iceVolumeCell" "icePressure" "uVelocityGeo" "vVelocityGeo" )
    FILE_EXTE_FIELD="mpassi.hist.am.timeSeriesStats${AVE_FREQ}.${YEAR}-${MONTH}-01.nc"

    # set series and month to analyze
    STATS+=( "IceExtent" "IceVolume" "SnowVolume" "KineticEnergy" )
    STATSMODIFY="total"
    FILE_EXTE_STATS="mpassi.hist.am.regionalStatistics.${YEAR}.${MONTH}.nc"

    # find all available combinations in cases for analysis
    if [ -z ${CASE2_ROOT} ]; then
     local casekombos=(${CASE_ROOT}/run.k*)
     local casescript=(${CASE_ROOT}/case_scripts.k*)
    else
     local casekombos=(${CASE_ROOT}/run.k* ${CASE2_ROOT}/run.k*)
     local casescript=(${CASE_ROOT}/case_scripts.k* ${CASE2_ROOT}/case_scripts.k*)
    fi 

    # check latest simulation run log being interrogated completed
    for ((i=0;i<${#casescript[@]};i++)); do
     pushd ${casescript[i]}
     local dirname=`dirname ${casescript[${i}]}`
     local casename=`basename ${dirname}`
     local komboname=`basename ${casescript[i]} | cut -d . -f 2`
     local runoutputs=(run.${casename}.*)
     if [[ ! -f ${runoutputs[-1]} ]]; then
      printf "\n    $casename $komboname: \t\x1B[31mNot run\e[0m" 
      printf "\n    $casename $komboname: \tNot run" >> "${pscript}"
     else
      for ((j=0;j<${#runoutputs[@]};j++)); do
       local lastline=`tail -1 ${runoutputs[j]}`
       if [[ ${lastline} =~ "CASE.RUN HAS FINISHED" ]]; then
        printf "\n    ${runoutputs[j]} $komboname: \x1B[32mComplete\e[0m" 
        printf "\n    ${runoutputs[j]} $komboname: Complete" >> "${pscript}"
       else
        printf "\n    ${runoutputs[j]} $komboname: \x1B[31mIncomplete\e[0m"
        printf "\n    ${runoutputs[j]} $komboname: Incomplete" >> "${pscript}"
       fi
      done
     fi
     popd
     printf "\n" | tee -a "${pscript}"
    done

    # check to make sure the necessary files exist, and for conservation checks
    for ((i=0;i<${#casekombos[@]};i++)); do
     dirnames[i]=`dirname ${casekombos[${i}]}`
     casenames[i]=`basename ${dirnames[${i}]}`
     if [ ! -d ${casekombos[${i}]} ]; then
      printf "\n--- Missing run directory for ${casenames[i]} \n\n"
      echo "'" >> ${pscript}
      exit 1
     else
      kombonames[i]="${casenames[i]} `basename ${casekombos[${i}]} | cut -d . -f 2`"
      emccheck[i]=false
      if [[ "${casenames[$i]}" =~ ".emc." ]]; then
       if [ ! -f ${casekombos[i]}/${casenames[i]}.${FILE_EXTE_EMC} ]; then
        printf "\n--- \x1B[31mMissing ${FILE_EXTE_EMC}\e[0m\n"
        printf "    for ${kombonames[i]}\n"
       else
        emccheck[i]=true
       fi
      fi
      if [ ! -f ${casekombos[i]}/${casenames[i]}.${FILE_EXTE_FIELD} ]; then
       printf "\n--- \x1B[31mMissing ${FILE_EXTE_FIELD}\e[0m\n"
       printf "    for ${kombonames[i]}\n"
      fi
      if [ ! -f ${casekombos[i]}/${casenames[i]}.${FILE_EXTE_STATS} ]; then
        printf "\n--- \x1B[31mMissing ${FILE_EXTE_STATS}\e[0m\n"
        printf "    for ${kombonames[i]}\n"
      fi
     fi
    done

    # check energy and mass if a simulation was completed and the -e option was used
    if [[ "${emccheck[@]}" =~ "true" ]]; then 
     printf "\n--- Global Energy and Mass Conservation for Month ${MONTH}, Year ${YEAR}:\n" | \
            tee -a "${pscript}"
     for ((i=0;i<${#casekombos[@]};i++)); do
      if [ "${emccheck[$i]}" == "true" ]; then
       textoutput="\n    ${kombonames[i]}"
       printf "\x1B[34m${textoutput}\e[0m:\n\n"
       printf "${textoutput}:\n\n" >> "${pscript}"
       for ((j=0;j<${#EMCVARS[@]};j++)); do
         EMCSTATS[j]=`ncks --fmt_val=%.20f -Q -P \\
                      -d Time,$((MONTH-1)),$((MONTH-1)) -d nHemispheres,0,0 \\
                      -v "${EMCVARMODIFY}${EMCVARS[j]}" \\
                      ${casekombos[i]}/${casenames[i]}.${FILE_EXTE_EMC} | \\
                      tail -4 | grep -Eo '[+-]?[0-9]+([.][0-9]+)?'`
         printf "    ${EMCVARMODIFY}${EMCVARS[j]}:  \t%+10.5e\n" "${EMCSTATS[j]}" | \
                 tee -a "${pscript}"
       done
      fi
     done
    fi

    # check to make sure there are two simulations to compare
    if [[ ${#casekombos[@]} < 2 ]]; then

     printf "\n--- Average Global Differences for Month ${MONTH}, Year ${YEAR}:\n" | \
             tee -a "${pscript}"

     for ((i=0;i<${#casekombos[@]};i++)); do

      if [ -f ${casekombos[i]}/${casenames[i]}.${FILE_EXTE_STATS} ]; then

       # provide statistics for regional values when no comparison available with another case
       for ((k=0;k<${#STATS[@]};k++)); do

         ncra -O -d nRegions,1,1 -v ${STATSMODIFY}${STATS[k]} \
                     ${casekombos[i]}/${casenames[i]}.${FILE_EXTE_STATS} \
                     ${casekombos[i]}/${casenames[i]}.${STATS[k]}.northav.nc

         ncra -O -d nRegions,2,2 -v ${STATSMODIFY}${STATS[k]} \
                     ${casekombos[i]}/${casenames[i]}.${FILE_EXTE_STATS} \
                     ${casekombos[i]}/${casenames[i]}.${STATS[k]}.southav.nc

         NSTAT[k]=`ncks --fmt_val=%.20f -Q -P -v ${STATSMODIFY}${STATS[k]} \\
               ${casekombos[i]}/${casenames[i]}.${STATS[k]}.northav.nc | \\
               tail -4 | grep -Eo '[+-]?[0-9]+([.][0-9]+)?'`

         SSTAT[k]=`ncks --fmt_val=%.20f -Q -P -v ${STATSMODIFY}${STATS[k]} \\
               ${casekombos[i]}/${casenames[i]}.${STATS[k]}.southav.nc | \\
               tail -4 | grep -Eo '[+-]?[0-9]+([.][0-9]+)?'`

         # add header
         if [[ $k -eq 0 ]]; then
          printf "\n\033[0;36m\t\t\tNorth \t\tSouth\e[0m\n"
          printf "\n\t\t\tNorth \t\tSouth\n" >> "${pscript}"
         fi

         printf "    ${STATS[k]}:  \t%+10.5e \t%+10.5e\n" "${NSTAT[k]}" "${SSTAT[k]}" | \
                 tee -a "${pscript}"

       done

      else
  
       printf "    \x1B[31mGlobal Statistics Unavailable\e[0m\n" 
       printf "    Global Statistics Unavailable\n" >> "${pscript}"

      fi

     done

     printf "\n--- Need two cases to compare\n\n"
     echo "'" >> ${pscript}
     exit

    fi

    # cycle through cases and compare them, output to screen and provenance 
    source /lcrc/soft/climate/e3sm-unified/load_latest_e3sm_unified_${MACHINES[${MACH}]}.sh

    printf "\n--- Average Field and Global Differences for Month ${MONTH}, Year ${YEAR}:\n" | \
            tee -a "${pscript}"

    for ((i=0;i<${#casekombos[@]};i++)); do
     for ((j=i+1;j<${#casekombos[@]};j++)); do

      textoutput="\n    ${kombonames[i]} - \n    ${kombonames[j]}"
      printf "\x1B[34m${textoutput}\e[0m:\n\n" 
      printf "${textoutput}:\n\n" >> "${pscript}"

      if [ -f ${casekombos[i]}/${casenames[i]}.${FILE_EXTE_FIELD} ] && \
         [ -f ${casekombos[j]}/${casenames[j]}.${FILE_EXTE_FIELD} ]; then

       # overall flag indicating if for all tested variables the test is BFB
       local bfbflag=true

       # check mesh sizes agree for comparing cases
       meshsizea=`ncdump -h ${casekombos[i]}/${casenames[i]}.${FILE_EXTE_FIELD} | \
                    grep "nCells =" | grep -Eo '[+-]?[0-9]+([.][0-9]+)?'`
       meshsizeb=`ncdump -h ${casekombos[j]}/${casenames[j]}.${FILE_EXTE_FIELD} | \
                    grep "nCells =" | grep -Eo '[+-]?[0-9]+([.][0-9]+)?'`

       # provide field differences for equivalent meshes
       if [[ ${meshsizea} -eq ${meshsizeb} ]]; then

        for ((k=0;k<${#FIELDS[@]};k++)); do

         ncdiff -O -v time${AVE_FREQ}_avg_${FIELDS[k]} \
                     ${casekombos[i]}/${casenames[i]}.${FILE_EXTE_FIELD} \
                     ${casekombos[j]}/${casenames[j]}.${FILE_EXTE_FIELD} \
                     ${casekombos[i]}/${casenames[i]}.${casenames[j]}.diff.nc

         ncwa -O -y min ${casekombos[i]}/${casenames[i]}.${casenames[j]}.diff.nc \
                        ${casekombos[i]}/${casenames[i]}.${casenames[j]}.min.nc

         ncwa -O -y max ${casekombos[i]}/${casenames[i]}.${casenames[j]}.diff.nc \
                        ${casekombos[i]}/${casenames[i]}.${casenames[j]}.max.nc

         MINDIFF[k]=`ncks --fmt_val=%.20f -Q -P -v time${AVE_FREQ}_avg_${FIELDS[k]} \\
               ${casekombos[i]}/${casenames[i]}.${casenames[j]}.min.nc | tail -4 | \\
               grep -Eo '[+-]?[0-9]+([.][0-9]+)?'`

         MAXDIFF[k]=`ncks --fmt_val=%.20f -Q -P -v time${AVE_FREQ}_avg_${FIELDS[k]} \\
               ${casekombos[i]}/${casenames[i]}.${casenames[j]}.max.nc | tail -4 | \\
               grep -Eo '[+-]?[0-9]+([.][0-9]+)?'`

         if (( $(echo "${MINDIFF[k]} != 0" | bc -l) )) &&
            (( $(echo "${MAXDIFF[k]} != 0" | bc -l) )) ; then
          bfbflag=false

          # add header
          if [[ $k -eq 0 ]]; then
           printf "\033[0;36m\t\t\tMinimum \tMaximum\e[0m\n" 
           printf "\t\t\tMinimum \tMaximum\n" >> "${pscript}"
          fi

          printf "    ${FIELDS[k]}: \t%+10.5e \t%+10.5e\n" "${MINDIFF[k]}" "${MAXDIFF[k]}" #| \ 
                 #tee -a "${pscript}"

         else
          printf "    ${FIELDS[k]}:  \t\x1B[32mBFB\e[0m\n" 
          printf "    ${FIELDS[k]}:  \tBFB\n" >> "${pscript}"
         fi

        done

       else

        printf "    \033[0;36mDifferent Meshes\e[0m" 
        printf "    Different Meshes" >> "${pscript}"

       fi

      else
  
       printf "    \x1B[31mFields Unavailable\e[0m\n" 
       printf "    Fields Unavailable\n" >> "${pscript}"

      fi

      printf "\n"


      # provide statistics differences, including for different meshes
      if [ -f ${casekombos[i]}/${casenames[i]}.${FILE_EXTE_STATS} ] && \
         [ -f ${casekombos[j]}/${casenames[j]}.${FILE_EXTE_STATS} ]; then

       # overall flag indicating if for all tested variables the test is BFB
       bfbflag=true

       # provide difference statistics for regional values regardless of mesh
       # averaged across the last full month
       for ((k=0;k<${#STATS[@]};k++)); do

         ncdiff -O -d nRegions,1,1 -v ${STATSMODIFY}${STATS[k]} \
                     ${casekombos[i]}/${casenames[i]}.${FILE_EXTE_STATS} \
                     ${casekombos[j]}/${casenames[j]}.${FILE_EXTE_STATS} \
                     ${casekombos[i]}/${casenames[i]}.${casenames[j]}.${STATS[k]}.northdiff.nc

         ncdiff -O -d nRegions,2,2 -v ${STATSMODIFY}${STATS[k]} \
                     ${casekombos[i]}/${casenames[i]}.${FILE_EXTE_STATS} \
                     ${casekombos[j]}/${casenames[j]}.${FILE_EXTE_STATS} \
                     ${casekombos[i]}/${casenames[i]}.${casenames[j]}.${STATS[k]}.southdiff.nc

         ncra -O ${casekombos[i]}/${casenames[i]}.${casenames[j]}.${STATS[k]}.northdiff.nc \
                 ${casekombos[i]}/${casenames[i]}.${casenames[j]}.${STATS[k]}.northdiffav.nc 

         ncra -O ${casekombos[i]}/${casenames[i]}.${casenames[j]}.${STATS[k]}.southdiff.nc  \
                 ${casekombos[i]}/${casenames[i]}.${casenames[j]}.${STATS[k]}.southdiffav.nc 

         NDIFF[k]=`ncks --fmt_val=%.20f -Q -P -v ${STATSMODIFY}${STATS[k]} \\
               ${casekombos[i]}/${casenames[i]}.${casenames[j]}.${STATS[k]}.northdiffav.nc | \\
               tail -4 | grep -Eo '[+-]?[0-9]+([.][0-9]+)?'`

         SDIFF[k]=`ncks --fmt_val=%.20f -Q -P -v ${STATSMODIFY}${STATS[k]} \\
               ${casekombos[i]}/${casenames[i]}.${casenames[j]}.${STATS[k]}.southdiffav.nc | \\
               tail -4 | grep -Eo '[+-]?[0-9]+([.][0-9]+)?'`

         if (( $(echo "${NDIFF[k]} != 0" | bc -l) )) &&
            (( $(echo "${SDIFF[k]} != 0" | bc -l) )) ; then
          bfbflag=false

          # add header
          if [[ $k -eq 0 ]]; then
           printf "\033[0;36m\t\t\tNorth \t\tSouth\e[0m\n" 
           printf "\t\t\tNorth \t\tSouth\n" >> "${pscript}"
          fi

          printf "    ${STATS[k]}:  \t%+10.5e \t%+10.5e\n" "${NDIFF[k]}" "${SDIFF[k]}" | \
                 tee -a "${pscript}"

         else
          printf "    ${STATS[k]}:  \t\x1B[32mBFB\e[0m\n" 
          printf "    ${STATS[k]}:  \tBFB\n" >> "${pscript}"
         fi

       done

      else
  
       printf "    \x1B[31mGlobal Statistics Unavailable\e[0m\n" 
       printf "    Global Statistics Unavailable\n" >> "${pscript}"

      fi

     done
    done

    echo "'" >> ${pscript}

    popd

}

#---------------------------------
# Silent versions of popd and pushd
pushd() {
    command pushd "$@" > /dev/null
}
popd() {
    command popd "$@" > /dev/null
}

# Run the script
#---------------------------------
main $* 
