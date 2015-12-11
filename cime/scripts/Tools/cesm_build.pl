#!/usr/bin/env perl 

# specify minimum version of perl
use 5.010;

use strict;
use warnings;
use File::Path qw(mkpath);
use File::Copy;
use File::Spec;
use File::Basename;
use Data::Dumper;
use Cwd;
use POSIX qw(strftime);

use English;
no if ($PERL_VERSION ge v5.18.0), 'warnings' => 'experimental::smartmatch';

#-----------------------------------------------------------------------------------------------
# Global data. 
#-----------------------------------------------------------------------------------------------

my $CASE;
my $CASEROOT;
my $CCSMROOT;
my $BUILD_THREADED;
my $CASEBUILD;
my $CASETOOLS;
my $EXEROOT;
my $CIMEROOT;
my $INCROOT;
my $LIBROOT;
my $SHAREDLIBROOT;
my $COMP_ATM;
my $COMP_LND;
my $COMP_ICE;
my $COMP_OCN;
my $COMP_GLC;
my $COMP_WAV;
my $COMP_ROF;
my $COMPILER;
my $COMP_INTERFACE;
my $MPILIB;
my $USE_ESMF_LIB;
my $DEBUG;
my $NINST_BUILD;
my $SMP_VALUE;
my $CLM_USE_PETSC;
my $CISM_USE_TRILINOS;
my $MPASLI_USE_ALBANY;
my $CLM_CONFIG_OPTS;
my $CAM_CONFIG_OPTS;
my $MACH;
my $OS;
my $COMP_CPL;
my $OCN_SUBMODEL;
my $PROFILE_PAPI_ENABLE;
my $DIN_LOC_ROOT;
my $GET_REFCASE;
my $RUN_TYPE;
my $RUN_REFDATE;
my $RUN_REFCASE;
my $RUN_REFDIR;
my $RUNDIR;
my $CONTINUE_RUN;
my $NTHRDS_CPL;
my $NTHRDS_ATM;
my $NTHRDS_LND;
my $NTHRDS_ICE;
my $NTHRDS_OCN;
my $NTHRDS_GLC;
my $NTHRDS_WAV;
my $NTHRDS_ROF;
my $NINST_ATM;
my $NINST_LND;
my $NINST_ICE;
my $NINST_OCN;
my $NINST_GLC;
my $NINST_WAV;
my $NINST_ROF;
my $NINST_VALUE;
my $SMP_BUILD;
my $BUILD_STATUS;
my $LOGDIR;
my $SHAREDPATH;

my $sysmod;

# Stash the build log paths here..
my @bldlogs;

my $LID = strftime("%y%m%d-%H%M%S", localtime);
my $banner = "-------------------------------------------------------------------------\n";

#-----------------------------------------------------------------------------------------------
sub main {

    if ($#ARGV == -1) {
	die " ERROR: must specify a caseroot input argument";
    }
    ($CASEROOT) = @ARGV;

    chdir "$CASEROOT" or die "Could not cd to $CASEROOT: $!\n";

    $sysmod = "./Tools/check_lockedfiles";
    system($sysmod) == 0 or die "$sysmod failed: $?\n";

    ($CASE, $CCSMROOT, $BUILD_THREADED, $CASEBUILD, $CASETOOLS, $EXEROOT, $CIMEROOT, $INCROOT, $LIBROOT, $SHAREDLIBROOT, $COMP_ATM, $COMP_LND, $COMP_ICE, $COMP_OCN, $COMP_GLC, $COMP_WAV, $COMP_ROF, $COMPILER, $COMP_INTERFACE, $MPILIB, $USE_ESMF_LIB, $DEBUG, $NINST_BUILD, $SMP_VALUE, $CLM_USE_PETSC, $CISM_USE_TRILINOS, $MPASLI_USE_ALBANY, $CLM_CONFIG_OPTS, $CAM_CONFIG_OPTS, $NINST_VALUE, $MACH, $OS, $COMP_CPL, $OCN_SUBMODEL, $PROFILE_PAPI_ENABLE, $DIN_LOC_ROOT, $GET_REFCASE, $RUN_TYPE, $RUN_REFDATE, $RUN_REFCASE, $RUN_REFDIR, $RUNDIR, $CONTINUE_RUN, $NTHRDS_CPL, $NTHRDS_ATM, $NTHRDS_LND, $NTHRDS_ICE, $NTHRDS_OCN, $NTHRDS_GLC, $NTHRDS_WAV, $NTHRDS_ROF, $NINST_ATM, $NINST_LND, $NINST_ICE, $NINST_OCN, $NINST_GLC, $NINST_WAV, $NINST_ROF, $SMP_BUILD, $BUILD_STATUS, $LOGDIR) = split('JGFSEP', `./xmlquery -s JGFSEP CASE CCSMROOT BUILD_THREADED CASEBUILD CASETOOLS EXEROOT CIMEROOT INCROOT LIBROOT SHAREDLIBROOT COMP_ATM COMP_LND COMP_ICE COMP_OCN COMP_GLC COMP_WAV COMP_ROF COMPILER COMP_INTERFACE MPILIB USE_ESMF_LIB DEBUG NINST_BUILD SMP_VALUE CLM_USE_PETSC CISM_USE_TRILINOS MPASLI_USE_ALBANY CLM_CONFIG_OPTS CAM_CONFIG_OPTS NINST_VALUE MACH OS COMP_CPL OCN_SUBMODEL PROFILE_PAPI_ENABLE DIN_LOC_ROOT GET_REFCASE RUN_TYPE RUN_REFDATE RUN_REFCASE RUN_REFDIR RUNDIR CONTINUE_RUN NTHRDS_CPL NTHRDS_ATM NTHRDS_LND NTHRDS_ICE NTHRDS_OCN NTHRDS_GLC NTHRDS_WAV NTHRDS_ROF NINST_ATM NINST_LND NINST_ICE NINST_OCN NINST_GLC NINST_WAV NINST_ROF SMP_BUILD BUILD_STATUS LOGDIR`);

    if (! -f "$CASE.run") {
	die "ERROR: must invoke cesm_setup script before calling build script ";
    }

    $ENV{CIMEROOT}		= $CIMEROOT		;
    $ENV{CASETOOLS}		= $CASETOOLS		;
    $ENV{EXEROOT}		= $EXEROOT		;
    $ENV{INCROOT}		= $INCROOT		;
    $ENV{LIBROOT}		= $LIBROOT		;
    $ENV{SHAREDLIBROOT}		= $SHAREDLIBROOT	;
    $ENV{CASEROOT}		= $CASEROOT		;
    $ENV{COMPILER}		= $COMPILER		;
    $ENV{COMP_INTERFACE}	= $COMP_INTERFACE	;
    $ENV{NINST_VALUE}		= $NINST_VALUE	        ;
    $ENV{BUILD_THREADED}	= $BUILD_THREADED	;
    $ENV{MACH}			= $MACH			;
    $ENV{USE_ESMF_LIB}		= $USE_ESMF_LIB		;
    $ENV{MPILIB}		= $MPILIB		;	
    $ENV{DEBUG}			= $DEBUG		;	
    $ENV{OS}			= $OS			;
    $ENV{COMP_CPL}		= $COMP_CPL		;	
    $ENV{COMP_ATM}		= $COMP_ATM		;	
    $ENV{COMP_LND}		= $COMP_LND		;	
    $ENV{COMP_ICE}		= $COMP_ICE		;	
    $ENV{COMP_OCN}		= $COMP_OCN		;	
    $ENV{COMP_GLC}		= $COMP_GLC		;	
    $ENV{COMP_WAV}		= $COMP_WAV		;	
    $ENV{COMP_ROF}		= $COMP_ROF		;	
    $ENV{CLM_CONFIG_OPTS}       = $CLM_CONFIG_OPTS      ;
    $ENV{CAM_CONFIG_OPTS}       = $CAM_CONFIG_OPTS      ;
    $ENV{OCN_SUBMODEL}          = $OCN_SUBMODEL;
    $ENV{PROFILE_PAPI_ENABLE}   = $PROFILE_PAPI_ENABLE;
#pw    $ENV{LID}  =  "`date +%y%m%d-%H%M%S`";
#pw++
    my $lid  =  "`date +%y%m%d-%H%M%S`";
    $ENV{LID}  =  $lid;
#pw--

    # Set the overall USE_PETSC variable to TRUE if any of the
    # XXX_USE_PETSC variables are TRUE.
    # For now, there is just the one CLM_USE_PETSC variable, but in
    # the future there may be others -- so USE_PETSC will be true if
    # ANY of those are true.

    my $use_petsc = 'FALSE';
    if ($CLM_USE_PETSC eq 'TRUE') {$use_petsc = 'TRUE'};
    my $sysmod = "./xmlchange -noecho -file env_build.xml -id USE_PETSC -val ${use_petsc}";
    $ENV{USE_PETSC} = ${use_petsc};
    $ENV{CLM_USE_PETSC} = $CLM_USE_PETSC;

    # Set the overall USE_TRILINOS variable to TRUE if any of the 
    # XXX_USE_TRILINOS variables are TRUE. 
    # For now, there is just the one CISM_USE_TRILINOS variable, but in
    # the future there may be others -- so USE_TRILINOS will be true if
    # ANY of those are true.

    my $use_trilinos = 'FALSE';
    if ($CISM_USE_TRILINOS eq 'TRUE') {$use_trilinos = 'TRUE'};
    $ENV{USE_TRILINOS} = ${use_trilinos};
    $ENV{CISM_USE_TRILINOS} = $CISM_USE_TRILINOS;

    # Set the overall USE_ALBANY variable to TRUE if any of the
    # XXX_USE_ALBANY variables are TRUE.
    # For now, there is just the one MPASLI_USE_ALBANY variable, but in
    # the future there may be others -- so USE_ALBANY will be true if
    # ANY of those are true.

    my $use_albany = 'FALSE';
    if ($MPASLI_USE_ALBANY eq 'TRUE') {$use_albany = 'TRUE'};
    $ENV{USE_ALBANY} = ${use_albany};
    $ENV{MPASLI_USE_ALBANY} = $MPASLI_USE_ALBANY;

    print "    .... checking namelists (calling ./preview_namelists) \n";
    $sysmod = "./preview_namelists > /dev/null";
    system($sysmod) == 0 or die "$sysmod failed: $?\n";
    
#pw++
    my $sdate = strftime("%y%m%d-%H%M%S", localtime);
    open my $CS, ">>", "./CaseStatus";
    print $CS "check input started $sdate\n";
#pw--
    checkInputData();
    buildChecks();
#pw++
    $sdate = strftime("%y%m%d-%H%M%S", localtime);
    print $CS "build started $sdate\n";
#pw--
    buildLibraries();
    buildModel();
#pw++
    $sdate = strftime("%y%m%d-%H%M%S", localtime);
    print $CS "build complete $sdate\n";
    close $CS;
#pw--
    
#pw++
    qx($CASEROOT/Tools/cesm_postbuild -cesmid $lid -cesmroot $CCSMROOT -caseroot $CASEROOT -exedir $EXEROOT)
#pw--

}

sub checkInputData()
{
    print "    .... calling data prestaging  \n";

    chdir "$CASEROOT" or die "Could not cd to $CASEROOT: $!\n";

    my @inputdatacheck = qx(./check_input_data -inputdata $DIN_LOC_ROOT -check);

    my @unknown = grep { /unknown/ } @inputdatacheck;
    if (@unknown) {
	print "      Any files with \"status uknown\" below were not found in the expected \n";
	print "      location, and are not from the input data repository. This is for \n";
	print "      informational only; this script will not attempt to find thse files. If\n";
	print "      If CESM can find (or does not need) these files no error will result. \n";
	map {print "$_\n" } @unknown;
    }
	
    my @missing = grep { /missing/ } @inputdatacheck;
    if (@missing) {
	print "Attempting to download missing data\n";
	qx(./check_input_data -inputdata $DIN_LOC_ROOT -export);

	print "Now checking if required input data is in $DIN_LOC_ROOT \n";

	@inputdatacheck = qx(./check_input_data -inputdata $DIN_LOC_ROOT -check);
	@missing = grep { /missing/ } @inputdatacheck;
	if (@missing) {
	    print "The following files were not found, they are required\n";
	    map {print "$_\n" } @missing;
	    print "Invoke the following command to obtain them:\n";
	    print "./check_input_data -inputdata $DIN_LOC_ROOT -export";
	    print "\n";
	    die;
	}
    }
	
    if( ($GET_REFCASE eq 'TRUE') && ($RUN_TYPE ne 'startup') && ($CONTINUE_RUN eq 'FALSE') )  {
	my $refdir = "${RUN_REFDIR}/${RUN_REFCASE}/${RUN_REFDATE}";
	if (! -d "${DIN_LOC_ROOT}/${refdir}") {

	    print "***************************************************************** \n";
	    print "ccsm_prestage ERROR: $refdir is not on local disk \n";
	    print "obtain this data from the svn input data repository \n";
	    print "> mkdir -p $refdir \n";
	    print "> cd $refdir \n";
	    print "> cd ..\n";
	    print "> svn export --force https://svn-ccsm-inputdata.cgd.ucar.edu/trunk/inputdata/${refdir} \n";
	    print " or set GET_REFCASE to FALSE in env_run.xml \n";
	    print " and prestage the restart data to $RUNDIR manually \n";
            print " ***************************************************************** \n";
	    die;

     	} else {

	    print " - Prestaging REFCASE ($refdir) to $RUNDIR\n";
		
	    # prestage the reference case's files.
	    mkpath ($RUNDIR) if (! -d $RUNDIR);
		
	    my @refcasefiles = glob("${DIN_LOC_ROOT}/${refdir}/*${RUN_REFCASE}*");
	    foreach my $rcfile (@refcasefiles) {

		my $rcbaseline = basename($rcfile);
		if(! -f "${RUNDIR}/$rcbaseline") {
		    my $sysmod = "ln -s $rcfile $RUNDIR/.";
		    system($sysmod) == 0 or warn "$sysmod failed: $?\n";
		}
			
		# copy the refcases' rpointer files to the run directory
		my @rpointerfiles = glob("${DIN_LOC_ROOT}/$refdir/*rpointer*");
		foreach my $rpointerfile(@rpointerfiles) {
		    copy($rpointerfile, ${RUNDIR});
		}
		chdir "$RUNDIR" or die "Could not cd to $RUNDIR: $!\n";

		my @cam2_list = glob("*.cam2.*");
		foreach my $cam2file(@cam2_list) {
		    my $camfile = $cam2file;
		    $camfile =~ s/cam2/cam/g;
		    symlink($cam2file, $camfile);
		}
	
		my @allrundirfiles = glob("$RUNDIR/*");
		foreach my $runfile(@allrundirfiles) {
		    chmod 0755, $runfile;
		}
	    }
	}
    }
}


sub buildChecks()
{
    print "    .... calling cesm build checks \n";
	
    chdir "$CASEROOT" or die "Could not cd to $CASEROOT: $!\n";
	
    my $atmstr = 0;
    my $lndstr = 0;
    my $icestr = 0;
    my $ocnstr = 0;
    my $rofstr = 0;
    my $glcstr = 0;
    my $wavstr = 0;
    my $cplstr = 0;

    if ($NTHRDS_ATM > 1 || $BUILD_THREADED eq 'TRUE') {$atmstr = 1;}
    if ($NTHRDS_LND > 1 || $BUILD_THREADED eq 'TRUE') {$lndstr = 1;}
    if ($NTHRDS_OCN > 1 || $BUILD_THREADED eq 'TRUE') {$ocnstr = 1;}
    if ($NTHRDS_ROF > 1 || $BUILD_THREADED eq 'TRUE') {$rofstr = 1;}
    if ($NTHRDS_GLC > 1 || $BUILD_THREADED eq 'TRUE') {$glcstr = 1;}
    if ($NTHRDS_WAV > 1 || $BUILD_THREADED eq 'TRUE') {$wavstr = 1;}
    if ($NTHRDS_CPL > 1 || $BUILD_THREADED eq 'TRUE') {$cplstr = 1;}
	
    $ENV{'SMP'} = 'FALSE';
    if( $NTHRDS_ATM > 1 || $NTHRDS_CPL > 1 || $NTHRDS_ICE > 1 ||
	$NTHRDS_LND > 1 || $NTHRDS_OCN > 1 || $NTHRDS_GLC > 1 || $NTHRDS_WAV > 1) { $ENV{'SMP'} = 'TRUE';}

    my $smpstr = "a$atmstr"."l$lndstr"."r$rofstr"."i$icestr"."o$ocnstr". "g$glcstr"."w$wavstr"."c$cplstr";

    $sysmod = "./xmlchange -noecho -file env_build.xml -id SMP_VALUE -val $smpstr";
    system($sysmod) == 0 or die "$sysmod failed: $?\n";

    $ENV{'SMP_VALUE'} = $smpstr;
	
    my $inststr = "a$NINST_ATM"."l$NINST_LND"."r$NINST_ROF"."i$NINST_ICE"."o$NINST_OCN"."g$NINST_GLC"."w$NINST_WAV";

    $sysmod = "./xmlchange -noecho -file env_build.xml -id NINST_VALUE -val $inststr";
    system($sysmod) == 0 or die "$sysmod failed: $?\n";

    $ENV{'NINST_VALUE'} = $inststr;
	
    # set the overall USE_PETSC variable to TRUE if any of the XXX_USE_PETSC variables are TRUE.
    # For now, there is just the one CLM_USE_PETSC variable, but in the future, there may be others,
    # so USE_PETSC should be  true if ANY of those are true.

    $ENV{'use_petsc'} = 'FALSE';
    if ( (defined $CLM_USE_PETSC) && ($CLM_USE_PETSC eq 'TRUE')) {$ENV{'use_petsc'} = "TRUE";}

    $sysmod = "./xmlchange -noecho -file env_build.xml -id USE_PETSC -val $ENV{'use_petsc'}";
    system($sysmod) == 0 or die "$sysmod failed: $?\n";

    # set the overall USE_TRILINOS variable to TRUE if any of the XXX_USE_TRILINOS variables are TRUE. 
    # For now, there is just the one CISM_USE_TRILINOS variable, but in the future, there may be others, 
    # so USE_TRILINOS should be  true if ANY of those are true.

    $ENV{'use_trilinos'} = 'FALSE';
    if ( (defined $CISM_USE_TRILINOS) && ($CISM_USE_TRILINOS eq 'TRUE')) {$ENV{'use_trilinos'} = "TRUE";}

    $sysmod = "./xmlchange -noecho -file env_build.xml -id USE_TRILINOS -val $ENV{'use_trilinos'}";
    system($sysmod) == 0 or die "$sysmod failed: $?\n";

    # set the overall USE_ALBANY variable to TRUE if any of the XXX_USE_ALBANY variables are TRUE.
    # For now, there is just the one MPASLI_USE_ALBANY variable, but in the future, there may be others,
    # so USE_ALBANY should be  true if ANY of those are true.

    $ENV{'use_albany'} = 'FALSE';
    if ( (defined $MPASLI_USE_ALBANY) && ($MPASLI_USE_ALBANY eq 'TRUE')) {$ENV{'use_albany'} = "TRUE";}

    $sysmod = "./xmlchange -noecho -file env_build.xml -id USE_ALBANY -val $ENV{'use_albany'}";
    system($sysmod) == 0 or die "$sysmod failed: $?\n";

    if( ($NINST_BUILD ne $NINST_VALUE) && ($NINST_BUILD != 0)) {
	print " ERROR, NINST VALUES HAVE CHANGED \n";
	print " NINST_BUILD = $NINST_BUILD \n";
	print " NINST_VALUE = $NINST_VALUE \n";
	print " A manual clean of your obj directories is strongly recommended \n";
	print " You should execute the following: \n";
	print " ./$CASE.clean_build \n";
	print " Then rerun the build script interactively \n";
	print " ---- OR ---- \n";
	print " You can override this error message at your own risk by executing:  \n";
	print "./xmlchange -file env_build.xml -id NINST_BUILD -val 0 \n";
	print " Then rerun the build script interactively \n";
	die;
    }

    if ($SMP_BUILD ne $SMP_VALUE && $SMP_BUILD != 0) {
	print "  ERROR SMP STATUS HAS CHANGED \n";
	print "  SMP_BUILD = $SMP_BUILD \n";
	print "  SMP_VALUE = $SMP_VALUE \n";
	print "  A manual clean of your obj directories is strongly recommended\n";
	print "  You should execute the following: \n";
	print "    ./$CASE.clean_build\n";
	print "  Then rerun the build script interactively\n";
	print "  ---- OR ----\n";
	print "  You can override this error message at your own risk by executing\n";
	print "    ./xmlchange -file env_build.xml -id SMP_BUILD -val 0\n";
	print "  Then rerun the build script interactively\n";
	die;
    }

    if($BUILD_STATUS != 0) {
	print "  ERROR env_build HAS CHANGED \n";
	print "  A manual clean of your obj directories is strongly recommended \n";
	print "  You should execute the following:  \n";
	print "      ./$CASE.clean_build \n";
	print "  Then rerun the build script interactively \n";
	print "    ---- OR ---- \n";
	print "  You can override this error message at your own risk by executing  \n";
	print "      rm LockedFiles/env_build* \n";
	print "  Then rerun the build script interactively  \n";
	die;
    }

    if ($COMP_INTERFACE eq 'ESMF' && $USE_ESMF_LIB ne 'TRUE') {
	print " ERROR COMP_INTERFACE IS ESMF BUT USE_ESMF_LIB IS NOT TRUE \n";
	print " SET USE_ESMF_LIB to TRUE with  \n";
	print "     ./xmlchange -file env_build.xml -id USE_ESMF_LIB -value TRUE \n";
	die;
    }
	
    if($MPILIB eq 'mpi-serial' && $USE_ESMF_LIB eq 'TRUE') {
	print "  ERROR MPILIB is mpi-serial and USE_ESMF_LIB IS TRUE \n";
	print "    MPILIB can only be used with an ESMF library built with mpiuni on \n";
	print "  Set USE_ESMF_LIB to FALSE with  \n";
	print "    ./xmlchange -file env_build.xml -id USE_ESMF_LIB -val FALSE \n";
	print "  ---- OR ---- \n";
	print "  Make suer the ESMF_LIBDIR used was built with mipuni (or change it to one that was) \n";
	print "  And comment out this if block in Tools/models_buildexe \n";
	die;
    }

    $sysmod = "./xmlchange -noecho -file env_build.xml -id BUILD_COMPLETE -val FALSE";
    system($sysmod) == 0 or die "$sysmod failed: $?\n";

    my @lockedfiles = glob("LockedFiles/env_build*");
    foreach my $lf (@lockedfiles) {
	unlink($lf);
    }
}


sub buildLibraries()
{
    print "    .... calling cesm builds for utility libraries (compiler is $COMPILER) \n";

    chdir $EXEROOT;

    if ($MPILIB eq 'mpi-serial') {
	my $sysmod = "cp -p -f $CIMEROOT/externals/mct/mpi-serial/\*.h  $LIBROOT/include/.";
	system($sysmod) == 0 or die "$sysmod failed: $?\n";
    }
    
    my $debugdir = "nodebug";
    if ($DEBUG eq 'TRUE') {$debugdir = "debug";}
    
    my $threaddir = 'nothreads';
    if ($ENV{'SMP'} eq 'TRUE' or $ENV{BUILD_THREADED} eq 'TRUE')
	{
		$threaddir = 'threads';
	}
    
    $ENV{'SHAREDPATH'}  = "$SHAREDLIBROOT/$COMPILER/$MPILIB/$debugdir/$threaddir";
    $SHAREDPATH = $ENV{'SHAREDPATH'};
    if(! -e $SHAREDPATH){ mkpath $SHAREDPATH;}
    mkpath("$SHAREDPATH/lib"    ) if (! -d "$SHAREDPATH/lib"    );
    mkpath("$SHAREDPATH/include") if (! -d "$SHAREDPATH/include");

    my @libs = qw/mct gptl pio csm_share/;
    print "      build libraries: @libs\n";

    foreach my $lib(@libs) {
	
	mkpath("$SHAREDPATH/$lib") if (! -d "$SHAREDPATH/$lib");
	chdir "$SHAREDPATH/$lib" or die "Could not cd to $SHAREDPATH/$lib: $!\n";

	my $file_build = "$SHAREDPATH/$lib.bldlog.$LID";
	my $now = localtime;
	print "      $now $file_build\n";
	open my $FB, ">", $file_build or die $!;
	map { print $FB "$_: $ENV{$_}\n"} sort keys %ENV;
	close $FB;
	
	eval {system("$CASEBUILD/buildlib.$lib $SHAREDPATH $CASEROOT >> $file_build 2>&1")};
	if ($?)	{
	    print "ERROR: buildlib.$lib failed, see $file_build\n";
	    die "ERROR: cat $file_build\n";
	}
	# push the file_build path into the bldlogs array..
	push(@bldlogs, $file_build);
    }
}

sub buildModel()
{
    print "    .... calling cesm builds for component libraries  \n";

    chdir "$CASEROOT" or die "Could not cd to $CASEROOT: $!\n";

    my @modelsbuildorder = qw( atm lnd ice ocn glc wav rof );
    my %models = ( atm => $COMP_ATM, lnd => $COMP_LND, ice => $COMP_ICE,
                   ocn => $COMP_OCN, glc => $COMP_GLC, wav => $COMP_WAV,
		   rof => $COMP_ROF);
    my $model;

    foreach $model(@modelsbuildorder) {

	my $comp = $models{$model};
	my $compspec = "";
	my $objdir = "";
	my $libdir = ""; 
	my $bldroot = "";

	if ("$comp" eq "clm") {

	    my $ESMFDIR;
	    if ($USE_ESMF_LIB eq "TRUE") {
		$ESMFDIR = "esmf";
	    } else {
		$ESMFDIR = "noesmf";
            }
	    for ("$CLM_CONFIG_OPTS") {
		when (/.*clm4_0.*/) {
		    print "         - Building clm4_0 Library \n";
		    $objdir = "$EXEROOT/$model/obj" ; if (! -d "$objdir") {mkpath "$objdir"};
		    $libdir = "$EXEROOT/$model"     ; if (! -d "$libdir") {mkpath "$libdir"};
		    $compspec = "lnd";
		    print "       bldroot is $EXEROOT \n";
		    print "       objdir  is $objdir \n";
		    print "       libdir  is $libdir \n";
		} default {
		    print "         - Building clm4_5/clm5_0 shared library \n";
		    $bldroot = "$SHAREDPATH/$COMP_INTERFACE/$ESMFDIR/" ;
		    $objdir  = "$bldroot/$comp/obj" ; if (! -d "$objdir") {mkpath "$objdir"};
		    $libdir  = "$bldroot/lib"       ; if (! -d "$libdir") {mkpath "$libdir"};
		    $compspec = "clm";
		    print "       bldroot is $bldroot \n";
		    print "       objdir  is $objdir \n";
		    print "       libdir  is $libdir \n";
		}
	    }

	}  else {

	    $objdir = "$EXEROOT/$model/obj" ; if (! -d "$objdir") {mkpath -p "$objdir";}
	    $libdir = "$EXEROOT/$model"     ; if (! -d "$libdir") {mkpath -p "$libdir";}
	    $compspec = $comp;
	}

	$ENV{'MODEL'} = $model;
	my $file_build = "$EXEROOT/${model}.bldlog.$LID";
	my $now = localtime;
        print "      $now $file_build\n";

	# build the component library
	chdir "$EXEROOT/$model" or die "Could not cd to $EXEROOT/$model: $!\n";
	eval{ system("$CASEBUILD/$comp.buildlib $CASEROOT $bldroot $compspec >> $file_build 2>&1") };
	if($?) { die "ERROR: $comp.buildlib failed, see $file_build\n";	}

	#push the file_build path into the bldlogs array..
	push (@bldlogs, $file_build);
		
	#--- copy .mod files... 
	my @lcmods = glob("$objdir/*_comp_*.mod");
	my @ucmods = glob("$objdir/*_COMP_*.mod");
	foreach my $mod (@lcmods, @ucmods) {
	    copy($mod, $INCROOT);
	}
    }

    my $file_build = "$EXEROOT/cesm.bldlog.$LID";
    my $now = localtime;
    print "      $now $file_build\n";

    mkpath "$EXEROOT/cesm/obj" if (! -d "$EXEROOT/cesm/obj");
    mkpath "$EXEROOT/cesm"     if (! -d "$EXEROOT/cesm");

    # create the model executable 
    chdir "$EXEROOT/cesm" or die "Could not cd to $EXEROOT/cesm: $!\n";
    eval{ system("$CASEBUILD/model.buildexe $CASEROOT >> $file_build 2>&1") };
    if ($?) {die "ERROR: model.buildexe failed, see $file_build\n";}

    push(@bldlogs, $file_build);
	
    #--- Copy the just-built cesm.exe to cesm.exe.$LID
    copy("$EXEROOT/cesm.exe", "$EXEROOT/cesm.exe.$LID");
    chmod 0755, "$EXEROOT/cesm.exe.$LID" or warn "could not change perms on $EXEROOT/cesm.exe.$LID, $?";
	
    #copy build logs to CASEROOT/logs
    if(length($LOGDIR) > 0) {
	if(! -d "$LOGDIR/bld") {
	    mkpath "$LOGDIR/bld";
	}
	chdir "$EXEROOT" or die "Could not cd to $EXEROOT: $!\n";
	#system("gzip $EXEROOT/*.bldlog.$LID*");	
	#my @gzlogs = glob("$EXEROOT/*bldlog.$LID.*");
	
	#foreach my $gzlog(@gzlogs)
	foreach my $log(@bldlogs) {
	    system("gzip $log");
	}	
	my @gzlogs = glob("$EXEROOT/*bldlog.$LID.*");
	foreach my $gzlog(@gzlogs) {
	    copy($gzlog, "$LOGDIR/bld");
	}
    }
	
    chdir "$CASEROOT" or die "Could not cd to $CASEROOT: $!\n";
    
    $sysmod = "./xmlchange -noecho -file env_build.xml -id BUILD_COMPLETE -val TRUE";
    system($sysmod) == 0 or die "$sysmod failed: $?\n";
    
    $sysmod = "./xmlchange -noecho -file env_build.xml -id BUILD_STATUS -val 0";
    system($sysmod) == 0 or die "$sysmod failed: $?\n";
    
    $sysmod = "./xmlchange -noecho -file env_build.xml -id SMP_BUILD -val $SMP_VALUE";
    system($sysmod) == 0 or die "$sysmod failed: $?\n";
    
    $sysmod = "./xmlchange -noecho -file env_build.xml -id NINST_BUILD -val $NINST_BUILD";
    system($sysmod) == 0 or die "$sysmod failed: $?\n";
    
    my @files2unlink = glob("./LockedFiles/env_build*");
    foreach my $file2unlink(@files2unlink) {
	unlink($file2unlink);
    }
	
    foreach my $file (qw( ./env_build.xml ) ) {
	copy($file, "./LockedFiles/$file.locked");
	if ($?) { die "ERROR locking file $file, exiting.. ";}
	print " .... locking file $file\n";
    }
    
#pw    my $sdate = strftime("%y%m%d-%H%M%S", localtime);
#pw    open my $CS, ">>", "./CaseStatus";
#pw    print $CS "build complete $sdate\n";
#pw    close $CS;
}

main() unless caller; 

