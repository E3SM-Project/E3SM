#!/usr/bin/env perl
use strict;
use Cwd;
###This is a program to extract timing information based on cesm1_3_+.

if($#ARGV < 1 ){
   print "Usage: getTiming.pl -fin input_file \n";
   exit
}
my $caseroot = getcwd();   # current working directory
$ENV{CASEROOT}=$caseroot;  # put this in environment

my $dirs = "$caseroot/Tools";
unshift @INC, $dirs;
require XML::Lite;
require SetupTools;

# Retreve and expand all case xml variables
my %xmlvars=();
SetupTools::getxmlvars($caseroot, \%xmlvars);
foreach my $attr (keys %xmlvars) {
    $xmlvars{$attr} = SetupTools::expand_env_var($xmlvars{$attr}, \%xmlvars);
}



my $fin;
my $opt=shift(@ARGV);
while($#ARGV > -1){
   if($opt eq "-fin"){
      $fin= shift(@ARGV);
   }
   else {
     print "Usage: getTiming.pl -fin input_file \n";
     exit
   }
   $opt=shift(@ARGV);
}

my $m999         = -999;

my $inittype = "FALSE";
if ($xmlvars{CONTINUE_RUN} eq "FALSE" && $xmlvars{RUN_TYPE} eq "startup") {$inittype = "TRUE";}
if ($xmlvars{CONTINUE_RUN} eq "FALSE" && $xmlvars{RUN_TYPE} eq "hybrid" ) {$inittype = "TRUE";}

my $tlen = 1.0;
if ($xmlvars{NCPL_BASE_PERIOD} eq "decade") {$tlen = 3650.0;}
if ($xmlvars{NCPL_BASE_PERIOD} eq "year")   {$tlen = 365.0;}
if ($xmlvars{NCPL_BASE_PERIOD} eq "days")   {$tlen = 1.0;}
if ($xmlvars{NCPL_BASE_PERIOD} eq "hour")   {$tlen = 1.0/24.0;}
my($nprocs,$ncount);
&gettime2('CPL:CLOCK_ADVANCE ',$nprocs, $ncount);
my $nsteps = $ncount / $nprocs;
my $adays = $nsteps * $tlen / $xmlvars{ATM_NCPL};
my $odays = $adays;

if ($inittype eq "TRUE") {$odays = $adays - ($tlen / $xmlvars{OCN_NCPL}) ;}

my $atm = $xmlvars{NTASKS_ATM} * $xmlvars{NTHRDS_ATM} ;
my $lnd = $xmlvars{NTASKS_LND} * $xmlvars{NTHRDS_LND} ;
my $rof = $xmlvars{NTASKS_ROF} * $xmlvars{NTHRDS_ROF} ;
my $wav = $xmlvars{NTASKS_WAV} * $xmlvars{NTHRDS_WAV} ;
my $ice = $xmlvars{NTASKS_ICE} * $xmlvars{NTHRDS_ICE} ;
my $ocn = $xmlvars{NTASKS_OCN} * $xmlvars{NTHRDS_OCN} ;
my $glc = $xmlvars{NTASKS_GLC} * $xmlvars{NTHRDS_GLC} ;
my $cpl = $xmlvars{NTASKS_CPL} * $xmlvars{NTHRDS_CPL} ;

my $apemin = $xmlvars{ROOTPE_ATM};
my $lpemin = $xmlvars{ROOTPE_LND};
my $wpemin = $xmlvars{ROOTPE_ROF};
my $rpemin = $xmlvars{ROOTPE_WAV};
my $ipemin = $xmlvars{ROOTPE_ICE};
my $opemin = $xmlvars{ROOTPE_OCN};
my $gpemin = $xmlvars{ROOTPE_GLC};
my $cpemin = $xmlvars{ROOTPE_CPL};

my $apemax = $xmlvars{ROOTPE_ATM} + $xmlvars{NTASKS_ATM} * $xmlvars{PSTRID_ATM} - 1 ;
my $lpemax = $xmlvars{ROOTPE_LND} + $xmlvars{NTASKS_LND} * $xmlvars{PSTRID_LND} - 1 ;
my $rpemax = $xmlvars{ROOTPE_ROF} + $xmlvars{NTASKS_ROF} * $xmlvars{PSTRID_ROF} - 1 ;
my $wpemax = $xmlvars{ROOTPE_WAV} + $xmlvars{NTASKS_WAV} * $xmlvars{PSTRID_WAV} - 1 ;
my $ipemax = $xmlvars{ROOTPE_ICE} + $xmlvars{NTASKS_ICE} * $xmlvars{PSTRID_ICE} - 1 ;
my $opemax = $xmlvars{ROOTPE_OCN} + $xmlvars{NTASKS_OCN} * $xmlvars{PSTRID_OCN} - 1 ;
my $gpemax = $xmlvars{ROOTPE_GLC} + $xmlvars{NTASKS_GLC} * $xmlvars{PSTRID_GLC} - 1 ;
my $cpemax = $xmlvars{ROOTPE_CPL} + $xmlvars{NTASKS_CPL} * $xmlvars{PSTRID_CPL} - 1 ;

my $peminmax = $apemin;
if( $lpemin > $peminmax ) { $peminmax = $lpemin; }
if( $rpemin > $peminmax ) { $peminmax = $rpemin; }
if( $wpemin > $peminmax ) { $peminmax = $wpemin; }
if( $ipemin > $peminmax ) { $peminmax = $ipemin; }
if( $opemin > $peminmax ) { $peminmax = $opemin; }
if( $gpemin > $peminmax ) { $peminmax = $gpemin; }
##if( $cpemin > $peminmax ) { $peminmax = $cpemin; }
$peminmax = $peminmax + 1;   

my $maxoffset = 40;
my $extraoff  = 20;
my $aoffset = int(($maxoffset * $apemin) / $peminmax) + $extraoff;
my $loffset = int(($maxoffset * $lpemin) / $peminmax) + $extraoff;
my $roffset = int(($maxoffset * $rpemin) / $peminmax) + $extraoff;
my $woffset = int(($maxoffset * $wpemin) / $peminmax) + $extraoff;
my $ioffset = int(($maxoffset * $ipemin) / $peminmax) + $extraoff;
my $goffset = int(($maxoffset * $gpemin) / $peminmax) + $extraoff;
my $ooffset = int(($maxoffset * $opemin) / $peminmax) + $extraoff;
##$coffset = int(($maxoffset * $cpemin) / $peminmax);
my $coffset = 0;
my $lid = $ENV{lid};
my $timeroot = $ENV{timeroot};
my $date = localtime();

print "
---------------- CCSM TIMING PROFILE ---------------------

  Case        : $xmlvars{CASE}
  LID         : $lid
  Machine     : $xmlvars{MACH}
  Caseroot    : $xmlvars{CASEROOT}
  Timeroot    : $timeroot
  CCSM User   : $xmlvars{CCSMUSER}
  CCSM Tag    : $xmlvars{CCSM_REPOTAG}  (best guess)
  Curr Date   : $date

  grid        : $xmlvars{GRID}
  compset     : $xmlvars{CCSM_COMPSET}
  run_type    : $xmlvars{RUN_TYPE}, continue_run = $xmlvars{CONTINUE_RUN} (inittype = $inittype)
  stop_option : $xmlvars{STOP_OPTION}, stop_n = $xmlvars{STOP_N}
  run_length  : $adays days ($odays for ocean)


\n";

print ("  component       comp_pes    root_pe   tasks  x threads instances (stride) \n");
print ("  ---------        ------     -------   ------   ------  ---------  ------  \n");
printf("  cpl = %-8s   %-6u      %-6u   %-6u x %-6u  %-6u (%-6u) \n",$xmlvars{COMP_CPL},$cpl,$xmlvars{ROOTPE_CPL},$xmlvars{NTASKS_CPL},$xmlvars{NTHRDS_CPL},1,$xmlvars{PSTRID_CPL});
printf("  glc = %-8s   %-6u      %-6u   %-6u x %-6u  %-6u (%-6u) \n",$xmlvars{COMP_GLC},$glc,$xmlvars{ROOTPE_GLC},$xmlvars{NTASKS_GLC},$xmlvars{NTHRDS_GLC},$xmlvars{NINST_GLC},$xmlvars{PSTRID_GLC});
printf("  wav = %-8s   %-6u      %-6u   %-6u x %-6u  %-6u (%-6u) \n",$xmlvars{COMP_WAV},$wav,$xmlvars{ROOTPE_WAV},$xmlvars{NTASKS_WAV},$xmlvars{NTHRDS_WAV},$xmlvars{NINST_WAV},$xmlvars{PSTRID_WAV});
printf("  lnd = %-8s   %-6u      %-6u   %-6u x %-6u  %-6u (%-6u) \n",$xmlvars{COMP_LND},$lnd,$xmlvars{ROOTPE_LND},$xmlvars{NTASKS_LND},$xmlvars{NTHRDS_LND},$xmlvars{NINST_LND},$xmlvars{PSTRID_LND});
printf("  rof = %-8s   %-6u      %-6u   %-6u x %-6u  %-6u (%-6u) \n",$xmlvars{COMP_ROF},$rof,$xmlvars{ROOTPE_ROF},$xmlvars{NTASKS_ROF},$xmlvars{NTHRDS_ROF},$xmlvars{NINST_ROF},$xmlvars{PSTRID_ROF});
printf("  ice = %-8s   %-6u      %-6u   %-6u x %-6u  %-6u (%-6u) \n",$xmlvars{COMP_ICE},$ice,$xmlvars{ROOTPE_ICE},$xmlvars{NTASKS_ICE},$xmlvars{NTHRDS_ICE},$xmlvars{NINST_ICE},$xmlvars{PSTRID_ICE});
printf("  atm = %-8s   %-6u      %-6u   %-6u x %-6u  %-6u (%-6u) \n",$xmlvars{COMP_ATM},$atm,$xmlvars{ROOTPE_ATM},$xmlvars{NTASKS_ATM},$xmlvars{NTHRDS_ATM},$xmlvars{NINST_ATM},$xmlvars{PSTRID_ATM});
printf("  ocn = %-8s   %-6u      %-6u   %-6u x %-6u  %-6u (%-6u) \n",$xmlvars{COMP_OCN},$ocn,$xmlvars{ROOTPE_OCN},$xmlvars{NTASKS_OCN},$xmlvars{NTHRDS_OCN},$xmlvars{NINST_OCN},$xmlvars{PSTRID_OCN});
my($nmin,$nmax,$tmin,$tmax,$wmin,$wmax,$nullf);
my($fmin,$fmax,$lmin,$lmax,$rmin,$rmax,$imin,$imax,$gmin,$gmax);
my($amin,$amax,$omin,$omax,$otmin,$otmax,$wtmin,$wtmax);
&gettime(' CPL:INIT '    ,$nmin,$nmax,$nullf);
&gettime(' CPL:RUN_LOOP ',$tmin,$tmax,$nullf);
&gettime(' CPL:TPROF_WRITE ',$wtmin,$wtmax,$nullf);
&gettime(' CPL:FINAL '   ,$fmin,$fmax,$nullf);
&gettime(' CPL:LND_RUN ' ,$lmin,$lmax,$nullf);
&gettime(' CPL:ROF_RUN ' ,$rmin,$rmax,$nullf);
&gettime(' CPL:ICE_RUN ' ,$imin,$imax,$nullf);
&gettime(' CPL:GLC_RUN ' ,$gmin,$gmax,$nullf);
&gettime(' CPL:WAV_RUN ' ,$wmin,$wmax,$nullf);
&gettime(' CPL:ATM_RUN ' ,$amin,$amax,$nullf);
&gettime(' CPL:OCN_RUN ' ,$omin,$omax,$nullf);
&gettime(' CPL:OCNT_RUN ',$otmin,$otmax,$nullf);
# pick OCNT_RUN for tight coupling
if ($otmax > $omax) {
  $omin = $otmin;
  $omax = $otmax;
}
my($cmin,$cmax,$xmin,$xmax,$ocnwaittime,$null);
&gettime(' CPL:RUN ' ,$cmin,$cmax,$nullf);
&gettime(' CPL:COMM ',$xmin,$xmax,$nullf);

&gettime(' CPL:C2O_INITWAIT ',$ocnwaittime,$null,$nullf);
my $ocnrunitime;
if ( $odays != 0.0 ) {
   $ocnrunitime = ($omax) * ($adays/$odays - 1.0);
} else {
   $ocnrunitime = 0.0;
}
my $correction = $ocnrunitime - $ocnwaittime;
if ($correction < 0) {$correction = 0.0;}

$tmax = $tmax + $wtmax + $correction;
$omax = $omax + $ocnrunitime;

my $pecost = $xmlvars{TOTALPES};
if ($xmlvars{COST_PES} > 0) {
    $pecost = $xmlvars{COST_PES};
}

print("\n");
print ("  total pes active           : $xmlvars{TOTALPES} \n");
print ("  pes per node               : $xmlvars{PES_PER_NODE} \n");
print ("  pe count for cost estimate : $pecost \n");
print("\n");

print ("  Overall Metrics: \n");
printf("    Model Cost:         %10.2f   pe-hrs/simulated_year \n",($tmax*365.*$pecost)/(3600.*$adays));
printf("    Model Throughput:   %10.2f   simulated_years/day \n",(86400.*$adays)/($tmax*365.));
print("\n");

printf("    Init Time   :  %10.3f seconds \n",$nmax);
printf("    Run Time    :  %10.3f seconds   %10.3f seconds/day \n",$tmax,$tmax/$adays);
printf("    Final Time  :  %10.3f seconds \n",$fmax);

print("\n");

printf("    Actual Ocn Init Wait Time     :  %10.3f seconds \n",$ocnwaittime);
printf("    Estimated Ocn Init Run Time   :  %10.3f seconds \n",$ocnrunitime);
printf("    Estimated Run Time Correction :  %10.3f seconds \n",$correction);
printf("      (This correction has been applied to the ocean and total run times) \n");

print("\n");
print("Runs Time in total seconds, seconds/model-day, and model-years/wall-day \n");
print("CPL Run Time represents time in CPL pes alone, not including time associated with data exchange with other components \n");
print("\n");

my $tmaxr = 0.0;
my $lmaxr = 0.0;
my $rmaxr = 0.0;
my $imaxr = 0.0;
my $amaxr = 0.0;
my $omaxr = 0.0;
my $gmaxr = 0.0;
my $wmaxr = 0.0;
my $cmaxr = 0.0;
my $xmaxr = 0.0;
if ($tmax > 0.0) { $tmaxr = ($adays*86400.)/($tmax*365.); }
if ($lmax > 0.0) { $lmaxr = ($adays*86400.)/($lmax*365.); }
if ($rmax > 0.0) { $rmaxr = ($adays*86400.)/($rmax*365.); }
if ($imax > 0.0) { $imaxr = ($adays*86400.)/($imax*365.); }
if ($amax > 0.0) { $amaxr = ($adays*86400.)/($amax*365.); }
if ($omax > 0.0) { $omaxr = ($adays*86400.)/($omax*365.); }
if ($gmax > 0.0) { $gmaxr = ($adays*86400.)/($gmax*365.); }
if ($wmax > 0.0) { $wmaxr = ($adays*86400.)/($wmax*365.); }
if ($cmax > 0.0) { $cmaxr = ($adays*86400.)/($cmax*365.); }
if ($xmax > 0.0) { $xmaxr = ($adays*86400.)/($xmax*365.); }

printf("    TOT Run Time:  %10.3f seconds   %10.3f seconds/mday   %10.2f myears/wday \n",$tmax,$tmax/$adays,$tmaxr);
printf("    LND Run Time:  %10.3f seconds   %10.3f seconds/mday   %10.2f myears/wday \n",$lmax,$lmax/$adays,$lmaxr);
printf("    ROF Run Time:  %10.3f seconds   %10.3f seconds/mday   %10.2f myears/wday \n",$rmax,$rmax/$adays,$rmaxr);
printf("    ICE Run Time:  %10.3f seconds   %10.3f seconds/mday   %10.2f myears/wday \n",$imax,$imax/$adays,$imaxr);
printf("    ATM Run Time:  %10.3f seconds   %10.3f seconds/mday   %10.2f myears/wday \n",$amax,$amax/$adays,$amaxr);
printf("    OCN Run Time:  %10.3f seconds   %10.3f seconds/mday   %10.2f myears/wday \n",$omax,$omax/$adays,$omaxr);
printf("    GLC Run Time:  %10.3f seconds   %10.3f seconds/mday   %10.2f myears/wday \n",$gmax,$gmax/$adays,$gmaxr);
printf("    WAV Run Time:  %10.3f seconds   %10.3f seconds/mday   %10.2f myears/wday \n",$wmax,$wmax/$adays,$wmaxr);
printf("    CPL Run Time:  %10.3f seconds   %10.3f seconds/mday   %10.2f myears/wday \n",$cmax,$cmax/$adays,$cmaxr);
printf("    CPL COMM Time: %10.3f seconds   %10.3f seconds/mday   %10.2f myears/wday \n",$xmax,$xmax/$adays,$xmaxr);


print ("\n\n---------------- DRIVER TIMING FLOWCHART --------------------- \n\n");

my $pstrlen = 25;
my $hoffset =  1;
print ("   NOTE: min:max driver timers (seconds/day):   \n");
my $xoff = $pstrlen+$hoffset+$coffset;
printf(" %${xoff}s CPL (pes %u to %u) \n",' ',$cpemin,$cpemax);
$xoff = $pstrlen+$hoffset+$ooffset;
printf(" %${xoff}s OCN (pes %u to %u) \n",' ',$opemin,$opemax);
$xoff = $pstrlen+$hoffset+$loffset;
printf(" %${xoff}s LND (pes %u to %u) \n",' ',$lpemin,$lpemax);
$xoff = $pstrlen+$hoffset+$roffset;
printf(" %${xoff}s ROF (pes %u to %u) \n",' ',$rpemin,$rpemax);
$xoff = $pstrlen+$hoffset+$ioffset;
printf(" %${xoff}s ICE (pes %u to %u) \n",' ',$ipemin,$ipemax);
$xoff = $pstrlen+$hoffset+$aoffset;
printf(" %${xoff}s ATM (pes %u to %u) \n",' ',$apemin,$apemax);
$xoff = $pstrlen+$hoffset+$goffset;
printf(" %${xoff}s GLC (pes %u to %u) \n",' ',$gpemin,$gpemax);
$xoff = $pstrlen+$hoffset+$woffset;
printf(" %${xoff}s WAV (pes %u to %u) \n",' ',$wpemin,$wpemax);
print ("\n");

&prttime(' CPL:CLOCK_ADVANCE '   ,$coffset,$adays,$m999);
&prttime(' CPL:OCNPRE1_BARRIER ' ,$coffset,$adays,$m999);
&prttime(' CPL:OCNPRE1 '         ,$coffset,$adays,$m999);
&prttime(' CPL:ATMOCN1_BARRIER ' ,$coffset,$adays,$m999);
&prttime(' CPL:ATMOCN1 '         ,$coffset,$adays,$m999);
&prttime(' CPL:OCNPREP_BARRIER ' ,$coffset,$adays,$m999);
&prttime(' CPL:OCNPREP '         ,$coffset,$adays,$m999);
&prttime(' CPL:C2O_BARRIER '     ,$ooffset,$odays,$coffset);
&prttime(' CPL:C2O '             ,$ooffset,$odays,$coffset);
&prttime(' CPL:LNDPREP_BARRIER ' ,$coffset,$adays,$m999);
&prttime(' CPL:LNDPREP '         ,$coffset,$adays,$m999);
&prttime(' CPL:C2L_BARRIER '     ,$loffset,$adays,$coffset);
&prttime(' CPL:C2L '             ,$loffset,$adays,$coffset);
&prttime(' CPL:ICEPREP_BARRIER ' ,$coffset,$adays,$m999);
&prttime(' CPL:ICEPREP '         ,$coffset,$adays,$m999);
&prttime(' CPL:C2I_BARRIER '     ,$ioffset,$adays,$coffset);
&prttime(' CPL:C2I '             ,$ioffset,$adays,$coffset);
&prttime(' CPL:WAVPREP_BARRIER ' ,$coffset,$adays,$m999);
&prttime(' CPL:WAVPREP '         ,$coffset,$adays,$m999);
&prttime(' CPL:C2W_BARRIER '     ,$ioffset,$adays,$coffset);
&prttime(' CPL:C2W '             ,$ioffset,$adays,$coffset);
&prttime(' CPL:ROFPREP_BARRIER ' ,$coffset,$adays,$m999);
&prttime(' CPL:ROFPREP '         ,$coffset,$adays,$m999);
&prttime(' CPL:C2R_BARRIER '     ,$roffset,$adays,$coffset);
&prttime(' CPL:C2R '             ,$roffset,$adays,$coffset);
&prttime(' CPL:ICE_RUN_BARRIER ' ,$ioffset,$adays,$m999);
&prttime(' CPL:ICE_RUN '         ,$ioffset,$adays,$m999);
&prttime(' CPL:LND_RUN_BARRIER ' ,$loffset,$adays,$m999);
&prttime(' CPL:LND_RUN '         ,$loffset,$adays,$m999);
&prttime(' CPL:ROF_RUN_BARRIER ' ,$roffset,$adays,$m999);
&prttime(' CPL:ROF_RUN '         ,$roffset,$adays,$m999);
&prttime(' CPL:WAV_RUN_BARRIER ' ,$roffset,$adays,$m999);
&prttime(' CPL:WAV_RUN '         ,$roffset,$adays,$m999);
&prttime(' CPL:OCNT_RUN_BARRIER ',$ooffset,$odays,$m999);
&prttime(' CPL:OCNT_RUN '        ,$ooffset,$odays,$m999);
&prttime(' CPL:O2CT_BARRIER '    ,$ooffset,$odays,$coffset);
&prttime(' CPL:O2CT '            ,$ooffset,$odays,$coffset);
&prttime(' CPL:OCNPOSTT_BARRIER ',$coffset,$adays,$m999);
&prttime(' CPL:OCNPOSTT '        ,$coffset,$adays,$m999);
&prttime(' CPL:ATMOCNP_BARRIER ' ,$coffset,$adays,$m999);
&prttime(' CPL:ATMOCNP '         ,$coffset,$adays,$m999);
&prttime(' CPL:L2C_BARRIER '     ,$loffset,$adays,$coffset);
&prttime(' CPL:L2C '             ,$loffset,$adays,$coffset);
&prttime(' CPL:LNDPOST_BARRIER ' ,$coffset,$adays,$m999);
&prttime(' CPL:LNDPOST '         ,$coffset,$adays,$m999);
&prttime(' CPL:GLCPREP_BARRIER ' ,$coffset,$adays,$m999);
&prttime(' CPL:GLCPREP '         ,$coffset,$adays,$m999);
&prttime(' CPL:C2G_BARRIER '     ,$goffset,$adays,$coffset);
&prttime(' CPL:C2G '             ,$goffset,$adays,$coffset);
&prttime(' CPL:R2C_BARRIER '     ,$roffset,$adays,$coffset);
&prttime(' CPL:R2C '             ,$roffset,$adays,$coffset);
&prttime(' CPL:ROFPOST_BARRIER ' ,$coffset,$adays,$m999);
&prttime(' CPL:ROFPOST '         ,$coffset,$adays,$m999);
&prttime(' CPL:BUDGET1_BARRIER ' ,$coffset,$adays,$m999);
&prttime(' CPL:BUDGET1 '         ,$coffset,$adays,$m999);
&prttime(' CPL:I2C_BARRIER '     ,$ioffset,$adays,$coffset);
&prttime(' CPL:I2C '             ,$ioffset,$adays,$coffset);
&prttime(' CPL:ICEPOST_BARRIER ' ,$coffset,$adays,$m999);
&prttime(' CPL:ICEPOST '         ,$coffset,$adays,$m999);
&prttime(' CPL:FRACSET_BARRIER ' ,$coffset,$adays,$m999);
&prttime(' CPL:FRACSET '         ,$coffset,$adays,$m999);
&prttime(' CPL:ATMOCN2_BARRIER ' ,$coffset,$adays,$m999);
&prttime(' CPL:ATMOCN2 '         ,$coffset,$adays,$m999);
&prttime(' CPL:OCNPRE2_BARRIER ' ,$coffset,$adays,$m999);
&prttime(' CPL:OCNPRE2 '         ,$coffset,$adays,$m999)
&prttime(' CPL:C2O2_BARRIER '    ,$ooffset,$odays,$coffset);
&prttime(' CPL:C2O2 '            ,$ooffset,$odays,$coffset);
&prttime(' CPL:ATMOCNQ_BARRIER'  ,$coffset,$adays,$m999);
&prttime(' CPL:ATMOCNQ '         ,$coffset,$adays,$m999);
&prttime(' CPL:ATMPREP_BARRIER ' ,$coffset,$adays,$m999);
&prttime(' CPL:ATMPREP '         ,$coffset,$adays,$m999);
&prttime(' CPL:C2A_BARRIER '     ,$aoffset,$adays,$coffset);
&prttime(' CPL:C2A '             ,$aoffset,$adays,$coffset);
&prttime(' CPL:OCN_RUN_BARRIER ' ,$ooffset,$odays,$m999);
&prttime(' CPL:OCN_RUN '         ,$ooffset,$odays,$m999);
&prttime(' CPL:ATM_RUN_BARRIER ' ,$aoffset,$adays,$m999);
&prttime(' CPL:ATM_RUN '         ,$aoffset,$adays,$m999);
&prttime(' CPL:GLC_RUN_BARRIER ' ,$goffset,$adays,$m999);
&prttime(' CPL:GLC_RUN '         ,$goffset,$adays,$m999);
&prttime(' CPL:W2C_BARRIER '     ,$goffset,$adays,$coffset);
&prttime(' CPL:W2C '             ,$goffset,$adays,$coffset);
&prttime(' CPL:WAVPOST_BARRIER ' ,$coffset,$adays,$m999);
&prttime(' CPL:WAVPOST '         ,$coffset,$adays,$m999);
&prttime(' CPL:G2C_BARRIER '     ,$goffset,$adays,$coffset);
&prttime(' CPL:G2C '             ,$goffset,$adays,$coffset);
&prttime(' CPL:GLCPOST_BARRIER ' ,$coffset,$adays,$m999);
&prttime(' CPL:GLCPOST '         ,$coffset,$adays,$m999);
&prttime(' CPL:A2C_BARRIER '     ,$aoffset,$adays,$coffset);
&prttime(' CPL:A2C '             ,$aoffset,$adays,$coffset);
&prttime(' CPL:ATMPOST_BARRIER ' ,$coffset,$adays,$m999);
&prttime(' CPL:ATMPOST '         ,$coffset,$adays,$m999);
&prttime(' CPL:BUDGET2_BARRIER ' ,$coffset,$adays,$m999);
&prttime(' CPL:BUDGET2 '         ,$coffset,$adays,$m999);
&prttime(' CPL:BUDGET3_BARRIER ' ,$coffset,$adays,$m999);
&prttime(' CPL:BUDGET3 '         ,$coffset,$adays,$m999);
&prttime(' CPL:BUDGETF_BARRIER ' ,$coffset,$adays,$m999);
&prttime(' CPL:BUDGETF '         ,$coffset,$adays,$m999);
&prttime(' CPL:O2C_BARRIER '     ,$ooffset,$odays,$coffset);
&prttime(' CPL:O2C '             ,$ooffset,$odays,$coffset);
&prttime(' CPL:OCNPOST_BARRIER ' ,$coffset,$adays,$m999);
&prttime(' CPL:OCNPOST '         ,$coffset,$adays,$m999);
&prttime(' CPL:RESTART_BARRIER ' ,$coffset,$adays,$m999);
&prttime(' CPL:RESTART'          ,$coffset,$adays,$m999);
&prttime(' CPL:HISTORY_BARRIER ' ,$coffset,$adays,$m999);
&prttime(' CPL:HISTORY '         ,$coffset,$adays,$m999);
&prttime(' CPL:TSTAMP_WRITE '    ,$coffset,$adays,$m999);
&prttime(' CPL:TPROF_WRITE '     ,$coffset,$adays,$m999);
&prttime(' CPL:RUN_LOOP_BSTOP '  ,$coffset,$adays,$m999);

#print ("\n --- overall total --- \n");
#&prttime(' CPL:RUN_LOOP '    ,$coffset,$adays,$m999);

print ("\n\n");

print ("More info on coupler timing:\n");

print ("\n");
&prttime(' CPL:OCNPRE1 '         ,$coffset,$adays,$m999);
&prttime(' CPL:ocnpre1_atm2ocn ' ,$coffset,$adays,$m999);
print ("\n");
&prttime(' CPL:OCNPREP '         ,$coffset,$adays,$m999);
&prttime(' CPL:OCNPRE2 '         ,$coffset,$adays,$m999);
&prttime(' CPL:ocnprep_avg '     ,$coffset,$adays,$m999);
&prttime(' CPL:ocnprep_diagav '  ,$coffset,$adays,$m999);
print ("\n");
&prttime(' CPL:LNDPREP '         ,$coffset,$adays,$m999);
&prttime(' CPL:lndprep_atm2lnd ' ,$coffset,$adays,$m999);
&prttime(' CPL:lndprep_mrgx2l '  ,$coffset,$adays,$m999);
&prttime(' CPL:lndprep_diagav '  ,$coffset,$adays,$m999);
print ("\n");
&prttime(' CPL:ICEPREP '         ,$coffset,$adays,$m999);
&prttime(' CPL:iceprep_ocn2ice ' ,$coffset,$adays,$m999);
&prttime(' CPL:iceprep_atm2ice ' ,$coffset,$adays,$m999);
&prttime(' CPL:iceprep_mrgx2i '  ,$coffset,$adays,$m999);
&prttime(' CPL:iceprep_diagav '  ,$coffset,$adays,$m999);
print ("\n");
&prttime(' CPL:WAVPREP '         ,$coffset,$adays,$m999);
&prttime(' CPL:wavprep_atm2wav ' ,$coffset,$adays,$m999);
&prttime(' CPL:wavprep_ocn2wav ' ,$coffset,$adays,$m999);
&prttime(' CPL:wavprep_ice2wav ' ,$coffset,$adays,$m999);
&prttime(' CPL:wavprep_mrgx2w '  ,$coffset,$adays,$m999);
&prttime(' CPL:wavprep_diagav '  ,$coffset,$adays,$m999);
print ("\n");
&prttime(' CPL:ROFPREP '         ,$coffset,$adays,$m999);
&prttime(' CPL:rofprep_l2xavg '  ,$coffset,$adays,$m999);
&prttime(' CPL:rofprep_lnd2rof ' ,$coffset,$adays,$m999);
&prttime(' CPL:rofprep_mrgx2r '  ,$coffset,$adays,$m999);
&prttime(' CPL:rofprep_diagav '  ,$coffset,$adays,$m999);
print ("\n");
&prttime(' CPL:GLCPREP '         ,$coffset,$adays,$m999);
&prttime(' CPL:glcprep_avg '     ,$coffset,$adays,$m999);
&prttime(' CPL:glcprep_lnd2glc ' ,$coffset,$adays,$m999);
&prttime(' CPL:glcprep_mrgx2g '  ,$coffset,$adays,$m999);
&prttime(' CPL:glcprep_diagav '  ,$coffset,$adays,$m999);
print ("\n");
&prttime(' CPL:ATMPREP '         ,$coffset,$adays,$m999);
&prttime(' CPL:atmprep_xao2atm ' ,$coffset,$adays,$m999);
&prttime(' CPL:atmprep_ocn2atm ' ,$coffset,$adays,$m999);
&prttime(' CPL:atmprep_alb2atm ' ,$coffset,$adays,$m999);
&prttime(' CPL:atmprep_ice2atm ' ,$coffset,$adays,$m999);
&prttime(' CPL:atmprep_lnd2atm ' ,$coffset,$adays,$m999);
&prttime(' CPL:atmprep_mrgx2a '  ,$coffset,$adays,$m999);
&prttime(' CPL:atmprep_diagav '  ,$coffset,$adays,$m999);
print ("\n");
&prttime(' CPL:ATMOCNP '         ,$coffset,$adays,$m999);
&prttime(' CPL:ATMOCN1 '         ,$coffset,$adays,$m999);
&prttime(' CPL:ATMOCN2 '         ,$coffset,$adays,$m999);
&prttime(' CPL:atmocnp_ice2ocn ' ,$coffset,$adays,$m999);
&prttime(' CPL:atmocnp_wav2ocn ' ,$coffset,$adays,$m999);
&prttime(' CPL:atmocnp_fluxo '   ,$coffset,$adays,$m999);
&prttime(' CPL:atmocnp_fluxe '   ,$coffset,$adays,$m999);
&prttime(' CPL:atmocnp_mrgx2o '  ,$coffset,$adays,$m999);
&prttime(' CPL:atmocnp_accum '   ,$coffset,$adays,$m999);
&prttime(' CPL:atmocnp_ocnalb '  ,$coffset,$adays,$m999);
print ("\n");
&prttime(' CPL:ATMOCNQ '         ,$coffset,$adays,$m999);
&prttime(' CPL:atmocnq_ocn2atm ' ,$coffset,$adays,$m999);
&prttime(' CPL:atmocnq_fluxa '   ,$coffset,$adays,$m999);
&prttime(' CPL:atmocnq_atm2ocnf ',$coffset,$adays,$m999);
print ("\n");
&prttime(' CPL:OCNPOSTT '        ,$coffset,$adays,$m999);
&prttime(' CPL:OCNPOST '         ,$coffset,$adays,$m999);
&prttime(' CPL:ocnpost_diagav '  ,$coffset,$adays,$m999);
print ("\n");
&prttime(' CPL:LNDPOST '         ,$coffset,$adays,$m999);
&prttime(' CPL:lndpost_diagav '  ,$coffset,$adays,$m999);
&prttime(' CPL:lndpost_acc2lr '  ,$coffset,$adays,$m999);
&prttime(' CPL:lndpost_acc2lg '  ,$coffset,$adays,$m999);
print ("\n");
&prttime(' CPL:ROFOST '          ,$coffset,$adays,$m999);
&prttime(' CPL:rofpost_diagav '  ,$coffset,$adays,$m999);
&prttime(' CPL:rofpost_histaux ' ,$coffset,$adays,$m999);
&prttime(' CPL:rofpost_rof2lnd ' ,$coffset,$adays,$m999);
&prttime(' CPL:rofpost_rof2ice ' ,$coffset,$adays,$m999);
&prttime(' CPL:rofpost_rof2ocn ' ,$coffset,$adays,$m999);
print ("\n");
&prttime(' CPL:ICEPOST '         ,$coffset,$adays,$m999);
&prttime(' CPL:icepost_diagav '  ,$coffset,$adays,$m999);
print ("\n");
&prttime(' CPL:WAVPOST '         ,$coffset,$adays,$m999);
&prttime(' CPL:wavpost_diagav '  ,$coffset,$adays,$m999);
print ("\n");
&prttime(' CPL:GLCPOST '         ,$coffset,$adays,$m999);
&prttime(' CPL:glcpost_diagav '  ,$coffset,$adays,$m999);
&prttime(' CPL:glcpost_glc2lnd ' ,$coffset,$adays,$m999);
&prttime(' CPL:glcpost_glc2ice ' ,$coffset,$adays,$m999);
&prttime(' CPL:glcpost_glc2ocn ' ,$coffset,$adays,$m999);
print ("\n");
&prttime(' CPL:ATMPOST '         ,$coffset,$adays,$m999);
&prttime(' CPL:atmpost_diagav '  ,$coffset,$adays,$m999);
print ("\n");
&prttime(' CPL:BUDGET '          ,$coffset,$adays,$m999);
&prttime(' CPL:BUDGET1 '         ,$coffset,$adays,$m999);
&prttime(' CPL:BUDGET2 '         ,$coffset,$adays,$m999);
&prttime(' CPL:BUDGET3 '         ,$coffset,$adays,$m999);
&prttime(' CPL:BUDGETF '         ,$coffset,$adays,$m999);
print ("\n\n");

exit;

#--------- end of main ------------
#----------------------------------
sub prttime{
   my($str,$offset,$div,$coff)=@_;
   my($min,$max,$found,$zoff,$maxd,$mind);
   my($cstr,$clen,$csp);

   my $datalen = 20;

   $cstr = "<---->";
   $clen = 6;

   &gettime($str,$min,$max,$found);
   $mind = $min;
   $maxd = $max;
   if ($div >= 1.0) {
     $mind = $min/$div;
     $maxd = $max/$div;
   }

   if ($mind >= 0.0 && $maxd >= 0.0 && $found > 0.5) {
      if ($coff >= 0) {
         $zoff = $pstrlen + $coff + int(($datalen-$clen)/2);
         $csp = $offset   - $coff - int(($datalen-$clen)/2);
         printf (" %-${zoff}s%-${csp}s %8.3f:%8.3f \n",$str,$cstr,$mind,$maxd);
      }
      else {
         $zoff = $pstrlen + $offset;
         printf (" %-${zoff}s %8.3f:%8.3f \n",$str,$mind,$maxd);
      }
   }
}

#----------------------------------

sub gettime{
   my($str,$min,$max,$found)=@_;
   my(@tmp,@tmp2);

#  set max here to something nonzero but small to avoid divide by zero
   $found = 0;
   $min = 0;
   $max = 0;

   my $strw = $str;
   $strw =~ s/^\s+//;
   $strw =~ s/\s+$//;
   @tmp=`grep -w "$strw" $fin | grep -E '\\(' | grep -v hashtable`;

#   print ("tcx1 $#tmp $tmp[0]\n");
   if ($#tmp == 0) {
#       print "tcx2 $tmp[0]\n";
#pw       if ($tmp[0] =~ m/\s*${strw}\s*\d+\s*(\d*\.\d+)\s*\(.*\)\s*(\d*\.\d+)\s*\(.*\)/) {
#jpe       if ($tmp[0] =~ m/\s*${strw}\s*\d+\s*\d+\s*\S+\s*\S+\s*(\d*\.\d+)\s*\(.*\)\s*(\d*\.\d+)\s*\(.*\)/) {
       if ($tmp[0] =~ m/\s*${strw}\s*[^\(]+\s+(\d*\.\d+)\s*\(.*\)\s*(\d*\.\d+)\s*\(.*\)/) {
#	   print "tcxa $tmp[0]\n";
#           print "tcxb $1 $2 \n";
	  $max=$1;
          $min=$2;
          $found=1.0;
       }
   }

   @_[1]=$min;
   @_[2]=$max;
   @_[3]=$found;

}

#----------------------------------

sub gettime2{
   my($str,$procs,$count)=@_;
   my(@tmp,@tmp2);

#  initialize procs and count
   $procs = 0;
   $count = 0;

   my $strw = $str;
   $strw =~ s/^\s+//;
   $strw =~ s/\s+$//;
   @tmp=`grep -w "$strw" $fin | grep -E '\\(' | grep -v hashtable`;

   if ($#tmp == 0) {
       if ($tmp[0] =~ m/\s*${strw}\s*([\.e\+\d]+)\s*(\d+)\s*/) {
	  $procs=$2;
	  $count=$1;
       }
   }
#   print "Here $procs $count\n";
   @_[1]=$procs;
   @_[2]=$count;

}

#----------------------------------
