version=0.1.6

# The rest of the script should not need to be modified
if [[ $HOSTNAME = "cori"* ]] || [[ $HOSTNAME = "dtn"* ]]; then
  base_path="/global/cfs/cdirs/acme/software/anaconda_envs/base"
elif [[ $HOSTNAME = "acme1"* ]] || [[ $HOSTNAME = "aims4"* ]]; then
  base_path="/usr/local/e3sm_unified/envs/base"
elif [[ $HOSTNAME = "blueslogin"* ]]; then
  base_path="/lcrc/soft/climate/e3sm-unified/base"
elif [[ $HOSTNAME = "rhea"* ]]; then
  base_path="/ccs/proj/cli900/sw/rhea/e3sm-unified/base"
elif [[ $HOSTNAME = "cooley"* ]]; then
  base_path="/lus/theta-fs0/projects/ccsm/acme/tools/e3sm-unified/base"
elif [[ $HOSTNAME = "compy"* ]]; then
  base_path="/share/apps/E3SM/conda_envs/base"
elif [[ $HOSTNAME = "gr-fe"* ]] || [[ $HOSTNAME = "ba-fe"* ]] || \
    [[ $HOSTNAME =~ ^gr[0-9]{4}$ ]] || [[ $HOSTNAME =~ ^ba[0-9]{3}$ ]]; then
  base_path="/usr/projects/climate/SHARED_CLIMATE/anaconda_envs/base"
else
  echo "Unknown host name $HOSTNAME.  Add base_path for this machine to the script."
fi

env_name=compass_${version}
if [ -d "$base_path/envs/$env_name" ]; then
  if [ -x "$(command -v module)" ] ; then
    module unload python
  fi
  source "${base_path}/etc/profile.d/conda.sh"
  conda activate "$env_name"
else
  echo "$env_name not found."
fi

