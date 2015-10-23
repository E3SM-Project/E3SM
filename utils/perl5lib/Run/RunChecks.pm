package Run::RunChecks;
use File::Basename;
my $pkg = 'Run::RunChecks';
use Log::Log4perl qw(get_logger);
my $logger;

BEGIN{
    $logger = get_logger();
}
# 

# Check that any files in $CASEROOT/LockedFiles matches the files
# in the $CASEROOT directory 
sub checkLockedFiles
{
	my $lockedfilesdir = undef; 
	my $caseroot = shift;
	if(! defined $caseroot)
	{
		$lockedfilesdir = "./LockedFiles";
	}
	else
	{
		$lockedfilesdir = "$caseroot/LockedFiles";
	}
	$logger->info( "Locked files dir: $lockedfilesdir");
	
	my @lockedfiles = glob("$lockedfilesdir/*");
	foreach my $lockedfile(@lockedfiles)
	{
		$logger->info( "locked file $lockedfile");
		# need to get the file name minus the path and extension. 
		my ($unlockedfile, $unlockedpath, $unlockedsuffix)  = fileparse($lockedfile);
	
		open my $UNLOCKEDFILE, "<", $unlockedfile  or $logger->logdie ("cannot open $unlockedfile");
		my $unlockeddata = <$UNLOCKEDFILE>;
		close $UNLOCKEDFILE;
		open my $LOCKEDFILE, "<", $lockedfile  or $logger->logdie ("cannot open $lockedfile");
		my $lockeddata = <$LOCKEDFILE>;
		close $LOCKEDFILE;

		if($unlockeddata ne $lockeddata)
		{
			$logger->info( "locked file $lockedfile has been modified and is different than the LockedFiles version");
			if($unlockedfile =~ /env_build/)
			{
				# TODO: Run xmlchange 
			#	./xmlchange -file env_build.xml -id BUILD_COMPLETE -val FALSE
            #    ./xmlchange -file env_build.xml -id BUILD_STATUS -val 1
			}
			elsif($unlockedfile =~ /env_mach_pes/)
			{
				$logger->error( "PE count has been changed!");
				$logger->error("please invoke case_setup -clean followed by case_setup");
			}
			else
			{
				$logger->error( "Cannot change $unlockedfile, please recover the original copy from the LockedFiles directory");
			}
		}
	}
}
1;
