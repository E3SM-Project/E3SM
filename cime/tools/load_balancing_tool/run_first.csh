#!/bin/csh -f

source global_variables.csh

set casebase = ${mach}_${compset}_${res}
set bldrun = 1

mkdir -p $casedir
mkdir -p ${results_dir}

set NATM = `echo $TASK_ATM:q | sed 's/,/ /g'`
set NLND = `echo $TASK_LND:q | sed 's/,/ /g'`
set NROF = `echo $TASK_ROF:q | sed 's/,/ /g'`
set NICE = `echo $TASK_ICE:q | sed 's/,/ /g'`
set NOCN = `echo $TASK_OCN:q | sed 's/,/ /g'`
set NCPL = `echo $TASK_CPL:q | sed 's/,/ /g'`
set NWAV = `echo $TASK_WAV:q | sed 's/,/ /g'`
set NGLC = `echo $TASK_GLC:q | sed 's/,/ /g'`

set RATM = `echo $ROOT_ATM:q | sed 's/,/ /g'`
set RLND = `echo $ROOT_LND:q | sed 's/,/ /g'`
set RROF = `echo $ROOT_LND:q | sed 's/,/ /g'`
set RICE = `echo $ROOT_ICE:q | sed 's/,/ /g'`
set ROCN = `echo $ROOT_OCN:q | sed 's/,/ /g'`
set RCPL = `echo $ROOT_CPL:q | sed 's/,/ /g'`
set RWAV = `echo $ROOT_WAV:q | sed 's/,/ /g'`
set RGLC = `echo $ROOT_LND:q | sed 's/,/ /g'`


set CNT = 0

if (-e ${results_dir}/test_list.out) then
  rm -f ${results_dir}/test_list.out
endif

foreach NA (${NATM})
@ CNT = $CNT + 1
  set NT = ${NTHRDS_VAL}
  @ EXPN = $CNT

  echo "case for $NA $CNT $NT "
  @ NTASKS_A = $NATM[$CNT]
  @ NTASKS_L = $NLND[$CNT]
  @ NTASKS_R = $NROF[$CNT]
  @ NTASKS_O = $NOCN[$CNT]
  @ NTASKS_I = $NICE[$CNT]
  @ NTASKS_C = $NCPL[$CNT]
  @ NTASKS_W = $NWAV[$CNT]
  @ NTASKS_G = $NGLC[$CNT]

  @ ROOTPE_ATM = $RATM[$CNT]
  @ ROOTPE_LND = $RLND[$CNT]
  @ ROOTPE_ROF = $RROF[$CNT]
  @ ROOTPE_ICE = $RICE[$CNT]
  @ ROOTPE_OCN = $ROCN[$CNT]
  @ ROOTPE_CPL = $RCPL[$CNT]
  @ ROOTPE_WAV = $RWAV[$CNT]
  @ ROOTPE_GLC = $RGLC[$CNT]

  set go = 1
  if ($NTASKS_A < 1) set NTASKS_A = 1
  if ($NTASKS_L < 1) set NTASKS_L = 1
  if ($NTASKS_I < 1) set NTASKS_I = 1
  if ($NTASKS_C < 1) set NTASKS_C = 1

  if ($go == 1) then
    echo "setting up case for $NA $CNT $NT "

    set case = t${casestr}${EXPN}_${NT}_${casebase}
    cd ${cesmsrc}/cime/scripts
    echo ${casedir}/${case} >> $results_dir/test_list.out
    ./create_newcase -case ${casedir}/${case} -res ${res} -compset {$compset} -mach ${mach}

    cd ${casedir}/${case}
    ./case_setup -clean

    #generic stuff
    ./xmlchange -file env_run.xml -id STOP_N -val $run_len
    ./xmlchange -file env_run.xml -id STOP_OPTION -val ndays
    ./xmlchange -file env_run.xml -id REST_OPTION -val never
    ./xmlchange -file env_run.xml -id TIMER_LEVEL -val 9
    ./xmlchange -file env_run.xml -id DOUT_S -val FALSE
    ./xmlchange -file env_run.xml -id COMP_RUN_BARRIERS -val TRUE

    ./xmlchange -file env_mach_pes.xml -id NTASKS_ATM -val $NTASKS_A
    ./xmlchange -file env_mach_pes.xml -id NTHRDS_ATM -val $NT
    ./xmlchange -file env_mach_pes.xml -id ROOTPE_ATM -val ${ROOTPE_ATM}

    ./xmlchange -file env_mach_pes.xml -id NTASKS_LND -val $NTASKS_L
    ./xmlchange -file env_mach_pes.xml -id NTHRDS_LND -val $NT
    ./xmlchange -file env_mach_pes.xml -id ROOTPE_LND -val ${ROOTPE_LND}

    ./xmlchange -file env_mach_pes.xml -id NTASKS_ROF -val $NTASKS_R
    ./xmlchange -file env_mach_pes.xml -id NTHRDS_ROF -val $NT
    ./xmlchange -file env_mach_pes.xml -id ROOTPE_ROF -val ${ROOTPE_ROF}

    ./xmlchange -file env_mach_pes.xml -id NTASKS_ICE -val $NTASKS_I
    ./xmlchange -file env_mach_pes.xml -id NTHRDS_ICE -val $NT
    ./xmlchange -file env_mach_pes.xml -id ROOTPE_ICE -val ${ROOTPE_ICE}

    ./xmlchange -file env_mach_pes.xml -id NTASKS_OCN -val $NTASKS_O
    ./xmlchange -file env_mach_pes.xml -id NTHRDS_OCN -val $NT
    ./xmlchange -file env_mach_pes.xml -id ROOTPE_OCN -val ${ROOTPE_OCN}

    ./xmlchange -file env_mach_pes.xml -id NTASKS_CPL -val $NTASKS_C
    ./xmlchange -file env_mach_pes.xml -id NTHRDS_CPL -val $NT
    ./xmlchange -file env_mach_pes.xml -id ROOTPE_CPL -val ${ROOTPE_CPL}

    ./xmlchange -file env_mach_pes.xml -id NTASKS_WAV -val $NTASKS_W
    ./xmlchange -file env_mach_pes.xml -id NTHRDS_WAV -val $NT
    ./xmlchange -file env_mach_pes.xml -id ROOTPE_WAV -val ${ROOTPE_WAV}

    ./xmlchange -file env_mach_pes.xml -id NTASKS_GLC -val $NTASKS_G
    ./xmlchange -file env_mach_pes.xml -id NTHRDS_GLC -val $NT
    ./xmlchange -file env_mach_pes.xml -id ROOTPE_GLC -val ${ROOTPE_GLC}

    ./case_setup

    if ($bldrun == "1") then
       ./${case}*.build

       rm tmpsubmit >& /dev/null
cat > tmpsubmit << EOF
   ./${case}*.submit
EOF
       source tmpsubmit
    endif  #bldrun

  endif  # go

end  #NATM

