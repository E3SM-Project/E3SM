#! /bin/bash
 


print_help(){
    echo "$0 [--help|-h] readme|profile[ step]|post_process /path/to/prerun_script.tool.sh [/path/to/case]"
}

if [ "$#" -lt 1 ]; then
  print_help; exit
fi

action=$1

if [ $action = "--help" ]; then
  print_help; exit
elif [ $action = "-h" ]; then
  print_help;  exit
elif [ $action = "readme" ]; then
  if [ "$#" -lt 2 ]; then
    print_help; exit
  else
    script="$(cd "$(dirname "$2")"; pwd)/$(basename "$2")"
    bash $script README
    exit
  fi
elif [ $action = "profile" ]; then
  if [ "$#" -lt 3 ]; then
    print_help; exit
  elif [ "$#" -eq 3 ]; then
    step=1
    script="$(cd "$(dirname "$2")"; pwd)/$(basename "$2")"
    case_name=$3
  elif [ "$#" -eq 4 ]; then
    step=$2
    script="$(cd "$(dirname "$3")"; pwd)/$(basename "$3")"
    case_name=$4
  fi

  cd $case_name
  ./xmlchange --file env_mach_specific.xml run_exe=" --multi-prog ./mpmd.\${LID}.conf "
  ./xmlchange --file env_run.xml PRERUN_SCRIPT="$script PROFILE $step"
  
  #now make changes to Macros.make if needed
  PROF_CPPFLAGS=$(bash $script CPPFLAGS)
  if [[ ! -z "$PROF_CPPFLAGS" ]]
  then
    echo "CPPDEFS := \$(CPPDEFS) $PROF_CPPFLAGS" >> Macros.make
  fi
  
  PROF_LDFLAGS=$(bash $script LDFLAGS)
  if [[ ! -z "$PROF_LDFLAGS" ]]
  then
    echo "LDFLAGS := \$(LDFLAGS) $PROF_LDFLAGS" >> Macros.make
  fi

  echo "The case may now be built and run:"
  echo "cd $case_name; ./case.build --clean; ./case.build; ./case.submit"
elif [ $action = "post_process" ]; then
    script="$(cd "$(dirname "$2")"; pwd)/$(basename "$2")"
    echo $script POSTPROCESS "${@:3}" 
    bash $script POSTPROCESS "${@:3}" 
    exit
fi
