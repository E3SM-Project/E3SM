package ConfigCESM;
my $pkg_nm = 'ConfigCESM';

use strict;
use English;
use IO::File;
use XML::LibXML;
use Data::Dumper;
use ConfigCase; # needed for cfg_ref reference below

# Check for the existence of XML::LibXML in whatever perl distribution happens to be in use.
# If not found, print a warning message then exit.
eval {
    require XML::LibXML;
    XML::LibXML->import();
};
if ($@)
{
    my $warning = <<END;
WARNING:
  The perl module XML::LibXML is needed for XML parsing in the CESM script system.
  Please contact your local systems administrators or IT staff and have them install it for
  you, or install the module locally.  

END
    print "$warning\n";
    exit(1);
}

#-------------------------------------------------------------------------------
sub getPrimaryComponent
{
    my ($compset) = @_;

    # TODO - this should be via parsing an xml file
    # put the list of modules in this xml file as well - as opposed to hard-coding

    # Note the primary component sets the compsets for associated with the "stand-alone" version of that
    # component via a compsets.xml file

    my $primary_component;
    if ($compset =~ m/%/) {

	# match longname
	if ($compset =~ m/_DATM/ && $compset =~ m/_DLND/ && $compset =~ m/_DICE/  && $compset =~ m/_DOCN/) {
	    $primary_component = 'drv';
	} elsif ($compset =~ m/_SATM/ && $compset =~ m/_SLND/ && $compset =~ m/_SICE/  && $compset =~ m/_SOCN/) {
	    $primary_component = 'drv';
	} elsif ($compset =~ m/_XATM/ && $compset =~ m/_XLND/ && $compset =~ m/_XICE/  && $compset =~ m/_XOCN/) {
	    $primary_component = 'drv';
	} elsif ($compset =~ m/_CAM/ && $compset =~ m/_CLM/ && $compset =~ m/_CICE/ && $compset =~ m/_POP/){
	    $primary_component = 'allactive';
	} elsif ($compset =~ m/_CAM/ && $compset =~ m/_CLM/ && $compset =~ m/_CICE/ && $compset =~ m/_DOCN%SOM/){
	    $primary_component = 'allactive';
	} elsif ($compset =~ m/_CAM/ && $compset =~ m/_CLM/ && $compset =~ m/_CICE/ && $compset =~ m/_DOCN%DOM/){
	    $primary_component = 'cam';
	} elsif ($compset =~ m/_DATM/ && $compset =~ m/_CLM/ && $compset =~ m/_SICE/ && $compset =~ m/_SOCN/){
	    $primary_component = 'clm';
	} elsif ($compset =~ m/_DATM/ && $compset =~ m/_SLND/ && $compset =~ m/_POP/ ) {
	    $primary_component = 'pop';
	} elsif ($compset =~ m/_DATM/ && $compset =~ m/_SLND/ && $compset =~ m/_CICE/  && $compset =~ m/_DOCN/) {
	    $primary_component = 'cice';
	} elsif ($compset =~ m/_SATM/ && $compset =~ m/_DLND/ &&$compset =~ m/_SICE/  && $compset =~ m/_SOCN/ && $compset =~ m/_CISM/) {
	    $primary_component = 'cism';
	}
	if (! defined $primary_component) {
	    die " create_newcase ERROR: primary_component $primary_component does not match supported primary components via longname";
	}

    } else {

	# match alias
	if ($compset =~ m/^[AXS].*/  ) {
	    $primary_component = 'drv';
	} elsif ($compset =~ m/^B.*/) {
	    $primary_component = 'allactive';
	} elsif ($compset =~ m/^E.*/) {
	    $primary_component = 'allactive';
	} elsif ($compset =~ m/^[FP].*/) {
	    $primary_component = 'cam';
	} elsif ($compset =~ m/^I.*/) {
	    $primary_component = 'clm';
	} elsif ($compset =~ m/^[CG].*/) {
	    $primary_component = 'pop';
	} elsif ($compset =~ m/^D.*/) {
	    $primary_component = 'cice';
	} elsif ($compset =~ m/^T.*/) {
	    $primary_component = 'cism';
	}
	if (! defined $primary_component) {
	    die " create_newcase ERROR: primary_component $primary_component does not match supported primary components via alias";
	}
    }
    return $primary_component;
}

#-------------------------------------------------------------------------------
1;
