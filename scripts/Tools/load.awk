BEGIN {
	train_to_load = ENVIRON["TRAIN_ID"];
	job_cnt = 0;
	printf "set full_jobid = ${PBS_JOBID}\n";
	printf "set jobid = `echo ${full_jobid} | awk -F. '{print $1}'`\n";
	printf "echo \"`date` - Starting Climate Train %s jobs\"\n\n", train_to_load;
}
{
	if (substr($1, 1, 1) != "#") {
		train_id=$1;
		case_name=$2;
		nprocs=$3;
		run_script=$4;
		log_dir=$5;
		email=$6;
		if (train_id == train_to_load) {
			printf "###############################################################################\n";
			printf "# %s\n", case_name;
			printf "###############################################################################\n";
			printf "mkdir -p %s\n", log_dir;
			if (email == "-") {
				printf "(time %s < /dev/null) >& %s/%s.o${jobid} ; set job_status = ${status}; echo \"`date` - %s exit status ${job_status}: `tail -1 %s/%s.o${jobid}`\" &\n", run_script, log_dir, case_name, case_name, log_dir, case_name;
			}
			else {
				printf "(time %s < /dev/null) >& %s/%s.o${jobid} ; set job_status = ${status}; echo \"`date` - %s exit status ${job_status}: `tail -1 %s/%s.o${jobid}`\"; echo \"`date` - %s exit status ${job_status}: `tail -1 %s/%s.o${jobid}`\" | /usr/bin/mailx -s \"%s in ${full_jobid}\" %s &\n", run_script, log_dir, case_name, case_name, log_dir, case_name, case_name, log_dir, case_name, case_name, email;
				printf "echo \"`date` - %s started\"\n", case_name;
				printf "echo \"`date` - %s started\" | /usr/bin/mailx -s \"%s in ${full_jobid}\" %s\n", case_name, case_name, email;
			}
			job_cnt++;
		}
	}
}
END {
	printf "\n";
	printf "###############################################################################\n";
	printf "# Wait for all jobs to finish\n";
	printf "###############################################################################\n";
	printf "echo \"`date` - Waiting for all Climate Train %s jobs to complete\"\n", train_to_load;
	printf "wait\n\n";
}
