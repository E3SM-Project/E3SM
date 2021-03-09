#!/usr/bin/env perl
#-----------------------------------------------------------------------
# Usage: perl REGRID.pl RESLN
#        where RESLN denotes the resolution desired.
#        Supported resolutions are:
#        T5, T21, T31, T42, T63, T85, T170,
#        10by15, 4by5, 2by2.5, 1by1.25, ncfilename (*new*)
#-----------------------------------------------------------------------

@ARGV;

$resln = $ARGV[0];


#-----------------------------------------------------
# Assign necessary dimensions given the resolution
#-----------------------------------------------------

if( $resln eq "T42" ) {
    $nlat = 64;
    $nlon = 128;

} elsif ( $resln eq "T63" ) {
    $nlat = 96;
    $nlon = 192;

} elsif ( $resln eq "T85" ) {
    $nlat = 128;
    $nlon = 256;

} elsif ( $resln eq "T170" ) {
    $nlat = 256;
    $nlon = 512;

} elsif ( $resln eq "T31" ) {
    $nlat = 48;
    $nlon = 96;

} elsif ( $resln eq "T21" ) {
    $nlat = 32;
    $nlon = 64;

} elsif ( $resln eq "T5" ) {
    $nlat = 8;
    $nlon = 16;

} elsif ( $resln eq "10by15" ) {
    $nlat = 19;
    $nlon = 24;

    $blat = "-90.0";
    $blon = "0.0";
    $elat = "90.0";
    $elon = "345.0";

} elsif ( $resln eq "4by5" ) {
    $nlat = 46;
    $nlon = 72;

    $blat = "-90.0";
    $blon = "0.0";
    $elat = "90.0";
    $elon = "355.0";

} elsif ( $resln eq "2by2.5" ) {
    $nlat = 91;
    $nlon = 144;
    $blat = "-90.0";
    $blon = "0.0";
    $elat = "90.0";
    $elon = "357.5";

} elsif ( $resln eq "1by1.25" ) {
    $nlat = 181;
    $nlon = 288;
    $blat = "-90.0";
    $blon = "0.0";
    $elat = "90.0";
    $elon = "358.75";

} elsif ( $resln eq "0.5by.625" ) {
    $nlat = 361;
    $nlon = 576;
    $blat = "-90.0";
    $blon = "0.0";
    $elat = "90.0";
    $elon = "359.375";

} elsif ( -e "$resln" ) { 
  print "Copying from input file $resln\n";
  $createfromfile=1;
  
} else {

  die "Usage: perl REGRID.pl <resolution> \n Supported resolutions are T5, T21, T31, T42, T63, T85, T170, 10by15, 4by5, 2by2.5, 1by1.25 or {ncfilename} where ncfilename is the name of an existing netcdf file containing variables lat and lon. \nExiting script";


}

if($createfromfile) {
  $outfile = "CREATE_DIMS_FROM_FILE.ncl";
  #------------------------------------------
  # Run NCL script
  #------------------------------------------

  print("CREATING DIMENSION FILE FOR $resln GRID\n");
  system( "ncl < CREATE_DIMS_FROM_FILE.ncl" );

}else{
  #---------------------------------------------
  # determine grid type (gaussian or regular)
  #---------------------------------------------

  $leadchar = substr($resln,0,1);

  if( $leadchar eq "T" ) {

    $indic = "GAUS";
    $infile = "CREATE_DIMS_GAU.ncl";

  }else{

    $indic = "REG";
    $infile = "CREATE_DIMS_REG.ncl";

  }

  #------------------------------------------------------------------
  # Set up NCL script to create netCDF file for requested resolution
  #------------------------------------------------------------------

  $keyword1 = "RESLN =";
  $keyword2 = "NLAT =";
  $keyword3 = "NLON =";
  $keyword4 = "BEGLAT =";
  $keyword5 = "BEGLON =";
  $keyword6 = "ENDLAT =";
  $keyword7 = "ENDLON =";

  $resexpr = qq!"$resln"!;
  
  $outfile = "CREATE_DIMS_$resln.ncl";

  open(INF, "$infile");
  open(OUTF,">$outfile");

  while(<INF>){

    if(/$keyword1/){
      print OUTF "RESLN = $resexpr\n";
    } elsif (/$keyword2/) {
      print OUTF "NLAT = $nlat\n";
    } elsif (/$keyword3/) {
      print OUTF "NLON = $nlon\n";
    } elsif (/$keyword4/) {
      print OUTF "BEGLAT = $blat\n";
    } elsif (/$keyword5/) {
      print OUTF "BEGLON = $blon\n";
    } elsif (/$keyword6/) {
      print OUTF "ENDLAT = $elat\n";
    } elsif (/$keyword7/) {
      print OUTF "ENDLON = $elon\n";
    } else {
      print OUTF "$_";
    }


  }

  close(INF);
  close(OUTF);
  #------------------------------------------
  # Run NCL script
  #------------------------------------------

  print("CREATING DIMENSION FILE FOR $resln GRID\n");
  system( "ncl < CREATE_DIMS_$resln.ncl" );
}



exit;




