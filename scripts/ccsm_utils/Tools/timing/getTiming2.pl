#!/usr/bin/env perl

###This is a program to extract timing information based on ccsm4_0_alpha versions.

if($#ARGV < 1 ){
   print "Usage: getTiming.pl -fin input_file \n";
   exit
}

$opt=shift(@ARGV);
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

$m999         = -999;

$case         = $ENV{'CASE'};
$lid          = $ENV{'lid'};
$mach         = $ENV{'MACH'};
$caseroot     = $ENV{'caseroot'};
$timeroot     = $ENV{'timeroot'};
$date         = $ENV{'date'};
$ccsmuser     = $ENV{'ccsmuser'};
$CCSMTAG      = $ENV{'CCSM_REPOTAG'};
$CCSMCS       = $ENV{'CCSM_COMPSET'};
$GRID         = $ENV{'GRID'};
$RUN_TYPE     = $ENV{'RUN_TYPE'};
$CONTINUE_RUN = $ENV{'CONTINUE_RUN'};
$STOP_OPTION  = $ENV{'STOP_OPTION'};
$STOP_N       = $ENV{'STOP_N'};
$COMP_CPL     = $ENV{'COMP_CPL'};
$NTASKS_CPL   = $ENV{'NTASKS_CPL'};
$NTHRDS_CPL   = $ENV{'NTHRDS_CPL'};
$ROOTPE_CPL   = $ENV{'ROOTPE_CPL'};
$PSTRID_CPL   = $ENV{'PSTRID_CPL'};
$COMP_ATM     = $ENV{'COMP_ATM'};
$NTASKS_ATM   = $ENV{'NTASKS_ATM'};
$NTHRDS_ATM   = $ENV{'NTHRDS_ATM'};
$ROOTPE_ATM   = $ENV{'ROOTPE_ATM'};
$NINST_ATM    = $ENV{'NINST_ATM'};
$PSTRID_ATM   = $ENV{'PSTRID_ATM'};
$COMP_LND     = $ENV{'COMP_LND'};
$NTASKS_LND   = $ENV{'NTASKS_LND'};
$NTHRDS_LND   = $ENV{'NTHRDS_LND'};
$ROOTPE_LND   = $ENV{'ROOTPE_LND'};
$NINST_LND    = $ENV{'NINST_LND'};
$PSTRID_LND   = $ENV{'PSTRID_LND'};
$COMP_ROF     = $ENV{'COMP_ROF'};
$NTASKS_ROF   = $ENV{'NTASKS_ROF'};
$NTHRDS_ROF   = $ENV{'NTHRDS_ROF'};
$ROOTPE_ROF   = $ENV{'ROOTPE_ROF'};
$NINST_ROF    = $ENV{'NINST_ROF'};
$PSTRID_ROF   = $ENV{'PSTRID_ROF'};
$COMP_WAV     = $ENV{'COMP_WAV'};
$NTASKS_WAV   = $ENV{'NTASKS_WAV'};
$NTHRDS_WAV   = $ENV{'NTHRDS_WAV'};
$ROOTPE_WAV   = $ENV{'ROOTPE_WAV'};
$NINST_WAV    = $ENV{'NINST_WAV'};
$PSTRID_WAV   = $ENV{'PSTRID_WAV'};
$COMP_ICE     = $ENV{'COMP_ICE'};
$NTASKS_ICE   = $ENV{'NTASKS_ICE'};
$NTHRDS_ICE   = $ENV{'NTHRDS_ICE'};
$ROOTPE_ICE   = $ENV{'ROOTPE_ICE'};
$NINST_ICE    = $ENV{'NINST_ICE'};
$PSTRID_ICE   = $ENV{'PSTRID_ICE'};
$COMP_OCN     = $ENV{'COMP_OCN'};
$NTASKS_OCN   = $ENV{'NTASKS_OCN'};
$NTHRDS_OCN   = $ENV{'NTHRDS_OCN'};
$ROOTPE_OCN   = $ENV{'ROOTPE_OCN'};
$NINST_OCN    = $ENV{'NINST_OCN'};
$PSTRID_OCN   = $ENV{'PSTRID_OCN'};
$COMP_GLC     = $ENV{'COMP_GLC'};
$NTASKS_GLC   = $ENV{'NTASKS_GLC'};
$NTHRDS_GLC   = $ENV{'NTHRDS_GLC'};
$ROOTPE_GLC   = $ENV{'ROOTPE_GLC'};
$NINST_GLC    = $ENV{'NINST_GLC'};
$PSTRID_GLC   = $ENV{'PSTRID_GLC'};
$TOTALPES     = $ENV{'TOTALPES'};
$PESPERNODE   = $ENV{'PES_PER_NODE'};
$COSTPES      = $ENV{'COST_PES'};
$NCPL_BASE_PERIOD = $ENV{'NCPL_BASE_PERIOD'};
$ATM_NCPL     = $ENV{'ATM_NCPL'};
$OCN_NCPL     = $ENV{'OCN_NCPL'};

$inittype = "FALSE";
if ($CONTINUE_RUN eq "FALSE" && $RUN_TYPE eq "startup") {$inittype = "TRUE";}
if ($CONTINUE_RUN eq "FALSE" && $RUN_TYPE eq "hybrid" ) {$inittype = "TRUE";}

$tlen = 1.0;
if ($NCPL_BASE_PERIOD eq "decade") {$tlen = 3650.0;}
if ($NCPL_BASE_PERIOD eq "year")   {$tlen = 365.0;}
if ($NCPL_BASE_PERIOD eq "days")   {$tlen = 1.0;}
if ($NCPL_BASE_PERIOD eq "hour")   {$tlen = 1.0/24.0;}

&gettime2(' DRIVER_CLOCK_ADVANCE ',$nprocs, $ncount);
$nsteps = $ncount / $nprocs;
$adays = $nsteps * $tlen / $ATM_NCPL;
$odays = $adays;
if ($inittype eq "TRUE") {$odays = $adays - ($tlen / $OCN_NCPL) ;}

$atm = $NTASKS_ATM * $NTHRDS_ATM ;
$lnd = $NTASKS_LND * $NTHRDS_LND ;
$rof = $NTASKS_ROF * $NTHRDS_ROF ;
$wav = $NTASKS_WAV * $NTHRDS_WAV ;
$ice = $NTASKS_ICE * $NTHRDS_ICE ;
$ocn = $NTASKS_OCN * $NTHRDS_OCN ;
$glc = $NTASKS_GLC * $NTHRDS_GLC ;
$cpl = $NTASKS_CPL * $NTHRDS_CPL ;

$apemin = $ROOTPE_ATM;
$lpemin = $ROOTPE_LND;
$wpemin = $ROOTPE_ROF;
$rpemin = $ROOTPE_WAV;
$ipemin = $ROOTPE_ICE;
$opemin = $ROOTPE_OCN;
$gpemin = $ROOTPE_GLC;
$cpemin = $ROOTPE_CPL;

$apemax = $ROOTPE_ATM + $NTASKS_ATM * $PSTRID_ATM - 1 ;
$lpemax = $ROOTPE_LND + $NTASKS_LND * $PSTRID_LND - 1 ;
$rpemax = $ROOTPE_ROF + $NTASKS_ROF * $PSTRID_ROF - 1 ;
$wpemax = $ROOTPE_WAV + $NTASKS_WAV * $PSTRID_WAV - 1 ;
$ipemax = $ROOTPE_ICE + $NTASKS_ICE * $PSTRID_ICE - 1 ;
$opemax = $ROOTPE_OCN + $NTASKS_OCN * $PSTRID_OCN - 1 ;
$gpemax = $ROOTPE_GLC + $NTASKS_GLC * $PSTRID_GLC - 1 ;
$cpemax = $ROOTPE_CPL + $NTASKS_CPL * $PSTRID_CPL - 1 ;

$totpes = $TOTALPES;

$peminmax = $apemin;
if( $lpemin > $peminmax ) { $peminmax = $lpemin; }
if( $rpemin > $peminmax ) { $peminmax = $rpemin; }
if( $wpemin > $peminmax ) { $peminmax = $wpemin; }
if( $ipemin > $peminmax ) { $peminmax = $ipemin; }
if( $opemin > $peminmax ) { $peminmax = $opemin; }
if( $gpemin > $peminmax ) { $peminmax = $gpemin; }
##if( $cpemin > $peminmax ) { $peminmax = $cpemin; }
$peminmax = $peminmax + 1;   

$maxoffset = 40;
$extraoff  = 20;
$aoffset = int(($maxoffset * $apemin) / $peminmax) + $extraoff;
$loffset = int(($maxoffset * $lpemin) / $peminmax) + $extraoff;
$roffset = int(($maxoffset * $rpemin) / $peminmax) + $extraoff;
$woffset = int(($maxoffset * $wpemin) / $peminmax) + $extraoff;
$ioffset = int(($maxoffset * $ipemin) / $peminmax) + $extraoff;
$goffset = int(($maxoffset * $gpemin) / $peminmax) + $extraoff;
$ooffset = int(($maxoffset * $opemin) / $peminmax) + $extraoff;
##$coffset = int(($maxoffset * $cpemin) / $peminmax);
$coffset = 0;

print "
---------------- CCSM TIMING PROFILE ---------------------

  Case        : $case
  LID         : $lid
  Machine     : $mach
  Caseroot    : $caseroot
  Timeroot    : $timeroot
  CCSM User   : $ccsmuser
  CCSM Tag    : $CCSMTAG  (best guess)
  Curr Date   : $date

  grid        : $GRID
  compset     : $CCSMCS
  run_type    : $RUN_TYPE, continue_run = $CONTINUE_RUN (inittype = $inittype)
  stop_option : $STOP_OPTION, stop_n = $STOP_N
  run_length  : $adays days ($odays for ocean)


\n";

print ("  component       comp_pes    root_pe   tasks  x threads instances (stride) \n");
print ("  ---------        ------     -------   ------   ------  ---------  ------  \n");
printf("  cpl = %-8s   %-6u      %-6u   %-6u x %-6u  %-6u (%-6u) \n",$COMP_CPL,$cpl,$ROOTPE_CPL,$NTASKS_CPL,$NTHRDS_CPL,1,$PSTRID_CPL);
printf("  glc = %-8s   %-6u      %-6u   %-6u x %-6u  %-6u (%-6u) \n",$COMP_GLC,$glc,$ROOTPE_GLC,$NTASKS_GLC,$NTHRDS_GLC,$NINST_GLC,$PSTRID_GLC);
printf("  wav = %-8s   %-6u      %-6u   %-6u x %-6u  %-6u (%-6u) \n",$COMP_WAV,$wav,$ROOTPE_WAV,$NTASKS_WAV,$NTHRDS_WAV,$NINST_WAV,$PSTRID_WAV);
printf("  lnd = %-8s   %-6u      %-6u   %-6u x %-6u  %-6u (%-6u) \n",$COMP_LND,$lnd,$ROOTPE_LND,$NTASKS_LND,$NTHRDS_LND,$NINST_LND,$PSTRID_LND);
printf("  rof = %-8s   %-6u      %-6u   %-6u x %-6u  %-6u (%-6u) \n",$COMP_ROF,$rof,$ROOTPE_ROF,$NTASKS_ROF,$NTHRDS_ROF,$NINST_ROF,$PSTRID_ROF);
printf("  ice = %-8s   %-6u      %-6u   %-6u x %-6u  %-6u (%-6u) \n",$COMP_ICE,$ice,$ROOTPE_ICE,$NTASKS_ICE,$NTHRDS_ICE,$NINST_ICE,$PSTRID_ICE);
printf("  atm = %-8s   %-6u      %-6u   %-6u x %-6u  %-6u (%-6u) \n",$COMP_ATM,$atm,$ROOTPE_ATM,$NTASKS_ATM,$NTHRDS_ATM,$NINST_ATM,$PSTRID_ATM);
printf("  ocn = %-8s   %-6u      %-6u   %-6u x %-6u  %-6u (%-6u) \n",$COMP_OCN,$ocn,$ROOTPE_OCN,$NTASKS_OCN,$NTHRDS_OCN,$NINST_OCN,$PSTRID_OCN);

&gettime(' DRIVER_INIT '    ,$nmin,$nmax,$nullf);
&gettime(' DRIVER_RUN_LOOP ',$tmin,$tmax,$nullf);
&gettime(' DRIVER_TPROF_WRITE ',$wmin,$wmax,$nullf);
&gettime(' DRIVER_FINAL '   ,$fmin,$fmax,$nullf);
&gettime(' DRIVER_LND_RUN ' ,$lmin,$lmax,$nullf);
&gettime(' DRIVER_ROF_RUN ' ,$rmin,$rmax,$nullf);
&gettime(' DRIVER_ICE_RUN ' ,$imin,$imax,$nullf);
&gettime(' DRIVER_GLC_RUN ' ,$gmin,$gmax,$nullf);
&gettime(' DRIVER_WAV_RUN ' ,$wmin,$wmax,$nullf);
&gettime(' DRIVER_ATM_RUN ' ,$amin,$amax,$nullf);
&gettime(' DRIVER_OCN_RUN ' ,$omin,$omax,$nullf);
&gettime(' DRIVER_OCNT_RUN ',$otmin,$otmax,$nullf);
# pick OCNT_RUN for tight coupling
if ($otmax > $omax) {
  $omin = $otmin;
  $omax = $otmax;
}
&gettime(' DRIVER_CPL_RUN ' ,$cmin,$cmax,$nullf);
&gettime(' DRIVER_CPL_COMM ',$xmin,$xmax,$nullf);

&gettime(' DRIVER_C2O_INITWAIT ',$ocnwaittime,$null,$nullf);
if ( $odays != 0.0 ) {
   $ocnrunitime = ($omax) * ($adays/$odays - 1.0);
} else {
   $ocnrunitime = 0.0;
}
$correction = $ocnrunitime - $ocnwaittime;
if ($correction < 0) {$correction = 0.0;}

$tmax = $tmax + wmax + $correction;
$omax = $omax + $ocnrunitime;

$pecost = $totpes;
if ($COSTPES > 0) {
    $pecost = $COSTPES;
}

print("\n");
print ("  total pes active           : $totpes \n");
print ("  pes per node               : $PESPERNODE \n");
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

$tmaxr = 0.0;
$lmaxr = 0.0;
$rmaxr = 0.0;
$imaxr = 0.0;
$amaxr = 0.0;
$omaxr = 0.0;
$gmaxr = 0.0;
$wmaxr = 0.0;
$cmaxr = 0.0;
$xmaxr = 0.0;
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

$pstrlen = 25;
$hoffset =  1;
print ("   NOTE: min:max driver timers (seconds/day):   \n");
$xoff = $pstrlen+$hoffset+$coffset;
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

&prttime(' DRIVER_CLOCK_ADVANCE '   ,$coffset,$adays,$m999);
&prttime(' DRIVER_OCNPRE1_BARRIER ' ,$coffset,$adays,$m999);
&prttime(' DRIVER_OCNPRE1 '         ,$coffset,$adays,$m999);
&prttime(' DRIVER_ATMOCN1_BARRIER ' ,$coffset,$adays,$m999);
&prttime(' DRIVER_ATMOCN1 '         ,$coffset,$adays,$m999);
&prttime(' DRIVER_OCNPREP_BARRIER ' ,$coffset,$adays,$m999);
&prttime(' DRIVER_OCNPREP '         ,$coffset,$adays,$m999);
&prttime(' DRIVER_C2O_BARRIER '     ,$ooffset,$odays,$coffset);
&prttime(' DRIVER_C2O '             ,$ooffset,$odays,$coffset);
&prttime(' DRIVER_LNDPREP_BARRIER ' ,$coffset,$adays,$m999);
&prttime(' DRIVER_LNDPREP '         ,$coffset,$adays,$m999);
&prttime(' DRIVER_C2L_BARRIER '     ,$loffset,$adays,$coffset);
&prttime(' DRIVER_C2L '             ,$loffset,$adays,$coffset);
&prttime(' DRIVER_ICEPREP_BARRIER ' ,$coffset,$adays,$m999);
&prttime(' DRIVER_ICEPREP '         ,$coffset,$adays,$m999);
&prttime(' DRIVER_C2I_BARRIER '     ,$ioffset,$adays,$coffset);
&prttime(' DRIVER_C2I '             ,$ioffset,$adays,$coffset);
&prttime(' DRIVER_WAVPREP_BARRIER ' ,$coffset,$adays,$m999);
&prttime(' DRIVER_WAVPREP '         ,$coffset,$adays,$m999);
&prttime(' DRIVER_C2W_BARRIER '     ,$ioffset,$adays,$coffset);
&prttime(' DRIVER_C2W '             ,$ioffset,$adays,$coffset);
&prttime(' DRIVER_ROFPREP_BARRIER ' ,$coffset,$adays,$m999);
&prttime(' DRIVER_ROFPREP '         ,$coffset,$adays,$m999);
&prttime(' DRIVER_C2R_BARRIER '     ,$roffset,$adays,$coffset);
&prttime(' DRIVER_C2R '             ,$roffset,$adays,$coffset);
&prttime(' DRIVER_ICE_RUN_BARRIER ' ,$ioffset,$adays,$m999);
&prttime(' DRIVER_ICE_RUN '         ,$ioffset,$adays,$m999);
&prttime(' DRIVER_LND_RUN_BARRIER ' ,$loffset,$adays,$m999);
&prttime(' DRIVER_LND_RUN '         ,$loffset,$adays,$m999);
&prttime(' DRIVER_ROF_RUN_BARRIER ' ,$roffset,$adays,$m999);
&prttime(' DRIVER_ROF_RUN '         ,$roffset,$adays,$m999);
&prttime(' DRIVER_WAV_RUN_BARRIER ' ,$roffset,$adays,$m999);
&prttime(' DRIVER_WAV_RUN '         ,$roffset,$adays,$m999);
&prttime(' DRIVER_OCNT_RUN_BARRIER ',$ooffset,$odays,$m999);
&prttime(' DRIVER_OCNT_RUN '        ,$ooffset,$odays,$m999);
&prttime(' DRIVER_O2CT_BARRIER '    ,$ooffset,$odays,$coffset);
&prttime(' DRIVER_O2CT '            ,$ooffset,$odays,$coffset);
&prttime(' DRIVER_OCNPOSTT_BARRIER ',$coffset,$adays,$m999);
&prttime(' DRIVER_OCNPOSTT '        ,$coffset,$adays,$m999);
&prttime(' DRIVER_ATMOCNP_BARRIER ' ,$coffset,$adays,$m999);
&prttime(' DRIVER_ATMOCNP '         ,$coffset,$adays,$m999);
&prttime(' DRIVER_L2C_BARRIER '     ,$loffset,$adays,$coffset);
&prttime(' DRIVER_L2C '             ,$loffset,$adays,$coffset);
&prttime(' DRIVER_LNDPOST_BARRIER ' ,$coffset,$adays,$m999);
&prttime(' DRIVER_LNDPOST '         ,$coffset,$adays,$m999);
&prttime(' DRIVER_GLCPREP_BARRIER ' ,$coffset,$adays,$m999);
&prttime(' DRIVER_GLCPREP '         ,$coffset,$adays,$m999);
&prttime(' DRIVER_C2G_BARRIER '     ,$goffset,$adays,$coffset);
&prttime(' DRIVER_C2G '             ,$goffset,$adays,$coffset);
&prttime(' DRIVER_R2C_BARRIER '     ,$roffset,$adays,$coffset);
&prttime(' DRIVER_R2C '             ,$roffset,$adays,$coffset);
&prttime(' DRIVER_ROFPOST_BARRIER ' ,$coffset,$adays,$m999);
&prttime(' DRIVER_ROFPOST '         ,$coffset,$adays,$m999);
&prttime(' DRIVER_BUDGET1_BARRIER ' ,$coffset,$adays,$m999);
&prttime(' DRIVER_BUDGET1 '         ,$coffset,$adays,$m999);
&prttime(' DRIVER_I2C_BARRIER '     ,$ioffset,$adays,$coffset);
&prttime(' DRIVER_I2C '             ,$ioffset,$adays,$coffset);
&prttime(' DRIVER_ICEPOST_BARRIER ' ,$coffset,$adays,$m999);
&prttime(' DRIVER_ICEPOST '         ,$coffset,$adays,$m999);
&prttime(' DRIVER_FRACSET_BARRIER ' ,$coffset,$adays,$m999);
&prttime(' DRIVER_FRACSET '         ,$coffset,$adays,$m999);
&prttime(' DRIVER_ATMOCN2_BARRIER ' ,$coffset,$adays,$m999);
&prttime(' DRIVER_ATMOCN2 '         ,$coffset,$adays,$m999);
&prttime(' DRIVER_OCNPRE2_BARRIER ' ,$coffset,$adays,$m999);
&prttime(' DRIVER_OCNPRE2 '         ,$coffset,$adays,$m999)
&prttime(' DRIVER_C2O2_BARRIER '    ,$ooffset,$odays,$coffset);
&prttime(' DRIVER_C2O2 '            ,$ooffset,$odays,$coffset);
&prttime(' DRIVER_ATMOCNQ_BARRIER'  ,$coffset,$adays,$m999);
&prttime(' DRIVER_ATMOCNQ '         ,$coffset,$adays,$m999);
&prttime(' DRIVER_ATMPREP_BARRIER ' ,$coffset,$adays,$m999);
&prttime(' DRIVER_ATMPREP '         ,$coffset,$adays,$m999);
&prttime(' DRIVER_C2A_BARRIER '     ,$aoffset,$adays,$coffset);
&prttime(' DRIVER_C2A '             ,$aoffset,$adays,$coffset);
&prttime(' DRIVER_OCN_RUN_BARRIER ' ,$ooffset,$odays,$m999);
&prttime(' DRIVER_OCN_RUN '         ,$ooffset,$odays,$m999);
&prttime(' DRIVER_ATM_RUN_BARRIER ' ,$aoffset,$adays,$m999);
&prttime(' DRIVER_ATM_RUN '         ,$aoffset,$adays,$m999);
&prttime(' DRIVER_GLC_RUN_BARRIER ' ,$goffset,$adays,$m999);
&prttime(' DRIVER_GLC_RUN '         ,$goffset,$adays,$m999);
&prttime(' DRIVER_W2C_BARRIER '     ,$goffset,$adays,$coffset);
&prttime(' DRIVER_W2C '             ,$goffset,$adays,$coffset);
&prttime(' DRIVER_WAVPOST_BARRIER ' ,$coffset,$adays,$m999);
&prttime(' DRIVER_WAVPOST '         ,$coffset,$adays,$m999);
&prttime(' DRIVER_G2C_BARRIER '     ,$goffset,$adays,$coffset);
&prttime(' DRIVER_G2C '             ,$goffset,$adays,$coffset);
&prttime(' DRIVER_GLCPOST_BARRIER ' ,$coffset,$adays,$m999);
&prttime(' DRIVER_GLCPOST '         ,$coffset,$adays,$m999);
&prttime(' DRIVER_A2C_BARRIER '     ,$aoffset,$adays,$coffset);
&prttime(' DRIVER_A2C '             ,$aoffset,$adays,$coffset);
&prttime(' DRIVER_ATMPOST_BARRIER ' ,$coffset,$adays,$m999);
&prttime(' DRIVER_ATMPOST '         ,$coffset,$adays,$m999);
&prttime(' DRIVER_BUDGET2_BARRIER ' ,$coffset,$adays,$m999);
&prttime(' DRIVER_BUDGET2 '         ,$coffset,$adays,$m999);
&prttime(' DRIVER_BUDGET3_BARRIER ' ,$coffset,$adays,$m999);
&prttime(' DRIVER_BUDGET3 '         ,$coffset,$adays,$m999);
&prttime(' DRIVER_BUDGETF_BARRIER ' ,$coffset,$adays,$m999);
&prttime(' DRIVER_BUDGETF '         ,$coffset,$adays,$m999);
&prttime(' DRIVER_O2C_BARRIER '     ,$ooffset,$odays,$coffset);
&prttime(' DRIVER_O2C '             ,$ooffset,$odays,$coffset);
&prttime(' DRIVER_OCNPOST_BARRIER ' ,$coffset,$adays,$m999);
&prttime(' DRIVER_OCNPOST '         ,$coffset,$adays,$m999);
&prttime(' DRIVER_RESTART_BARRIER ' ,$coffset,$adays,$m999);
&prttime(' DRIVER_RESTART'          ,$coffset,$adays,$m999);
&prttime(' DRIVER_HISTORY_BARRIER ' ,$coffset,$adays,$m999);
&prttime(' DRIVER_HISTORY '         ,$coffset,$adays,$m999);
&prttime(' DRIVER_TSTAMP_WRITE '    ,$coffset,$adays,$m999);
&prttime(' DRIVER_TPROF_WRITE '     ,$coffset,$adays,$m999);
&prttime(' DRIVER_RUN_LOOP_BSTOP '  ,$coffset,$adays,$m999);

#print ("\n --- overall total --- \n");
#&prttime(' DRIVER_RUN_LOOP '    ,$coffset,$adays,$m999);

print ("\n\n");

print ("More info on coupler timing:\n");

print ("\n");
&prttime(' DRIVER_OCNPRE1 '         ,$coffset,$adays,$m999);
&prttime(' driver_ocnpre1_atm2ocn ' ,$coffset,$adays,$m999);
print ("\n");
&prttime(' DRIVER_OCNPREP '         ,$coffset,$adays,$m999);
&prttime(' DRIVER_OCNPRE2 '         ,$coffset,$adays,$m999);
&prttime(' driver_ocnprep_avg '     ,$coffset,$adays,$m999);
&prttime(' driver_ocnprep_diagav '  ,$coffset,$adays,$m999);
print ("\n");
&prttime(' DRIVER_LNDPREP '         ,$coffset,$adays,$m999);
&prttime(' driver_lndprep_atm2lnd ' ,$coffset,$adays,$m999);
&prttime(' driver_lndprep_mrgx2l '  ,$coffset,$adays,$m999);
&prttime(' driver_lndprep_diagav '  ,$coffset,$adays,$m999);
print ("\n");
&prttime(' DRIVER_ICEPREP '         ,$coffset,$adays,$m999);
&prttime(' driver_iceprep_ocn2ice ' ,$coffset,$adays,$m999);
&prttime(' driver_iceprep_atm2ice ' ,$coffset,$adays,$m999);
&prttime(' driver_iceprep_mrgx2i '  ,$coffset,$adays,$m999);
&prttime(' driver_iceprep_diagav '  ,$coffset,$adays,$m999);
print ("\n");
&prttime(' DRIVER_WAVPREP '         ,$coffset,$adays,$m999);
&prttime(' driver_wavprep_atm2wav ' ,$coffset,$adays,$m999);
&prttime(' driver_wavprep_ocn2wav ' ,$coffset,$adays,$m999);
&prttime(' driver_wavprep_ice2wav ' ,$coffset,$adays,$m999);
&prttime(' driver_wavprep_mrgx2w '  ,$coffset,$adays,$m999);
&prttime(' driver_wavprep_diagav '  ,$coffset,$adays,$m999);
print ("\n");
&prttime(' DRIVER_ROFPREP '         ,$coffset,$adays,$m999);
&prttime(' driver_rofprep_l2xavg '  ,$coffset,$adays,$m999);
&prttime(' driver_rofprep_lnd2rof ' ,$coffset,$adays,$m999);
&prttime(' driver_rofprep_mrgx2r '  ,$coffset,$adays,$m999);
&prttime(' driver_rofprep_diagav '  ,$coffset,$adays,$m999);
print ("\n");
&prttime(' DRIVER_GLCPREP '         ,$coffset,$adays,$m999);
&prttime(' driver_glcprep_avg '     ,$coffset,$adays,$m999);
&prttime(' driver_glcprep_lnd2glc ' ,$coffset,$adays,$m999);
&prttime(' driver_glcprep_mrgx2g '  ,$coffset,$adays,$m999);
&prttime(' driver_glcprep_diagav '  ,$coffset,$adays,$m999);
print ("\n");
&prttime(' DRIVER_ATMPREP '         ,$coffset,$adays,$m999);
&prttime(' driver_atmprep_xao2atm ' ,$coffset,$adays,$m999);
&prttime(' driver_atmprep_ocn2atm ' ,$coffset,$adays,$m999);
&prttime(' driver_atmprep_alb2atm ' ,$coffset,$adays,$m999);
&prttime(' driver_atmprep_ice2atm ' ,$coffset,$adays,$m999);
&prttime(' driver_atmprep_lnd2atm ' ,$coffset,$adays,$m999);
&prttime(' driver_atmprep_mrgx2a '  ,$coffset,$adays,$m999);
&prttime(' driver_atmprep_diagav '  ,$coffset,$adays,$m999);
print ("\n");
&prttime(' DRIVER_ATMOCNP '         ,$coffset,$adays,$m999);
&prttime(' DRIVER_ATMOCN1 '         ,$coffset,$adays,$m999);
&prttime(' DRIVER_ATMOCN2 '         ,$coffset,$adays,$m999);
&prttime(' driver_atmocnp_ice2ocn ' ,$coffset,$adays,$m999);
&prttime(' driver_atmocnp_wav2ocn ' ,$coffset,$adays,$m999);
&prttime(' driver_atmocnp_fluxo '   ,$coffset,$adays,$m999);
&prttime(' driver_atmocnp_fluxe '   ,$coffset,$adays,$m999);
&prttime(' driver_atmocnp_mrgx2o '  ,$coffset,$adays,$m999);
&prttime(' driver_atmocnp_accum '   ,$coffset,$adays,$m999);
&prttime(' driver_atmocnp_ocnalb '  ,$coffset,$adays,$m999);
print ("\n");
&prttime(' DRIVER_ATMOCNQ '         ,$coffset,$adays,$m999);
&prttime(' driver_atmocnq_ocn2atm ' ,$coffset,$adays,$m999);
&prttime(' driver_atmocnq_fluxa '   ,$coffset,$adays,$m999);
&prttime(' driver_atmocnq_atm2ocnf ',$coffset,$adays,$m999);
print ("\n");
&prttime(' DRIVER_OCNPOSTT '        ,$coffset,$adays,$m999);
&prttime(' DRIVER_OCNPOST '         ,$coffset,$adays,$m999);
&prttime(' driver_ocnpost_diagav '  ,$coffset,$adays,$m999);
print ("\n");
&prttime(' DRIVER_LNDPOST '         ,$coffset,$adays,$m999);
&prttime(' driver_lndpost_diagav '  ,$coffset,$adays,$m999);
&prttime(' driver_lndpost_acc2lr '  ,$coffset,$adays,$m999);
&prttime(' driver_lndpost_acc2lg '  ,$coffset,$adays,$m999);
print ("\n");
&prttime(' DRIVER_ROFOST '          ,$coffset,$adays,$m999);
&prttime(' driver_rofpost_diagav '  ,$coffset,$adays,$m999);
&prttime(' driver_rofpost_histaux ' ,$coffset,$adays,$m999);
&prttime(' driver_rofpost_rof2lnd ' ,$coffset,$adays,$m999);
&prttime(' driver_rofpost_rof2ice ' ,$coffset,$adays,$m999);
&prttime(' driver_rofpost_rof2ocn ' ,$coffset,$adays,$m999);
print ("\n");
&prttime(' DRIVER_ICEPOST '         ,$coffset,$adays,$m999);
&prttime(' driver_icepost_diagav '  ,$coffset,$adays,$m999);
print ("\n");
&prttime(' DRIVER_WAVPOST '         ,$coffset,$adays,$m999);
&prttime(' driver_wavpost_diagav '  ,$coffset,$adays,$m999);
print ("\n");
&prttime(' DRIVER_GLCPOST '         ,$coffset,$adays,$m999);
&prttime(' driver_glcpost_diagav '  ,$coffset,$adays,$m999);
&prttime(' driver_glcpost_glc2lnd ' ,$coffset,$adays,$m999);
&prttime(' driver_glcpost_glc2ice ' ,$coffset,$adays,$m999);
&prttime(' driver_glcpost_glc2ocn ' ,$coffset,$adays,$m999);
print ("\n");
&prttime(' DRIVER_ATMPOST '         ,$coffset,$adays,$m999);
&prttime(' driver_atmpost_diagav '  ,$coffset,$adays,$m999);
print ("\n");
&prttime(' DRIVER_BUDGET '          ,$coffset,$adays,$m999);
&prttime(' DRIVER_BUDGET1 '         ,$coffset,$adays,$m999);
&prttime(' DRIVER_BUDGET2 '         ,$coffset,$adays,$m999);
&prttime(' DRIVER_BUDGET3 '         ,$coffset,$adays,$m999);
&prttime(' DRIVER_BUDGETF '         ,$coffset,$adays,$m999);
print ("\n\n");

exit;

#--------- end of main ------------
#----------------------------------
sub prttime{
   local($str,$offset,$div,$coff)=@_;
   local($min,$max,$found,$zoff,$maxd,$mind);
   local($cstr,$clen,$csp);

   $datalen = 20;

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
   local($str,$min,$max,$found)=@_;
   local(@tmp,@tmp2);

#  set max here to something nonzero but small to avoid divide by zero
   $found = 0;
   $min = 0;
   $max = 0;

   $strw = $str;
   $strw =~ s/^\s+//;
   $strw =~ s/\s+$//;
   @tmp=`grep -w "$strw" $fin | grep -E '\\(' | grep -v hashtable`;

#   print ("tcx1 $#tmp $tmp[0]\n");
   if ($#tmp == 0) {
#       print "tcx2 $tmp[0]\n";
#pw       if ($tmp[0] =~ m/\s*${strw}\s*\d+\s*(\d*\.\d+)\s*\(.*\)\s*(\d*\.\d+)\s*\(.*\)/) {
       if ($tmp[0] =~ m/\s*${strw}\s*\d+\s*\d+\s*\S+\s*\S+\s*(\d*\.\d+)\s*\(.*\)\s*(\d*\.\d+)\s*\(.*\)/) {
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
   local($str,$procs,$count)=@_;
   local(@tmp,@tmp2);

#  initialize procs and count
   $procs = 0;
   $count = 0;

   $strw = $str;
   $strw =~ s/^\s+//;
   $strw =~ s/\s+$//;
   @tmp=`grep -w "$strw" $fin | grep -E '\\(' | grep -v hashtable`;

   if ($#tmp == 0) {
       if ($tmp[0] =~ m/\s*${strw}\s*(\d+)\s*\d+\s*(\S+)/) {
	  $procs=$1;
	  $count=$2;
       }
   }

   @_[1]=$procs;
   @_[2]=$count;

}

#----------------------------------
