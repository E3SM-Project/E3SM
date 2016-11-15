package Run::RunChecks;
use File::Basename;
my $pkg = 'Run::RunChecks';

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
	print "Locked files dir: $lockedfilesdir\n";
	
	my @lockedfiles = glob("$lockedfilesdir/*");
	foreach my $lockedfile(@lockedfiles)
	{
		print "locked file $lockedfile\n";
		# need to get the file name minus the path and extension. 
		my ($unlockedfile, $unlockedpath, $unlockedsuffix)  = fileparse($lockedfile);
	
		open my $UNLOCKEDFILE, "<", $unlockedfile  or die "cannot open $unlockedfile";
		my $unlockeddata = <$UNLOCKEDFILE>;
		close $UNLOCKEDFILE;
		open my $LOCKEDFILE, "<", $lockedfile  or die "cannot open $lockedfile";
		my $lockeddata = <$LOCKEDFILE>;
		close $LOCKEDFILE;

		if($unlockeddata ne $lockeddata)
		{
			print "locked file $lockedfile has been modified and is different than the LockedFiles version\n";
			if($unlockedfile =~ /env_build/)
			{
				# TODO: Run xmlchange 
			#	./xmlchange -file env_build.xml -id BUILD_COMPLETE -val FALSE
            #    ./xmlchange -file env_build.xml -id BUILD_STATUS -val 1
			}
			elsif($unlockedfile =~ /env_mach_pes/)
			{
				print "PE count has been changed!\n";
				print "please invoke cesm_setup -clean followed by cesm_setup\n";
			}
			else
			{
				print "Cannot change $unlockedfile, please recover the original copy from the LockedFiles directory\n";
			}
		}
	}
}
