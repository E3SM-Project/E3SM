use strict;
use XML::Lite;
package Depends::Checks;

sub checkLibXML()
{
	my $libXMLInstalled = 1;
	eval
	{
		require XML::LibXML;
		XML::LibXML->import();
	};
	if($@)
	{
		$libXMLInstalled = 0;
		print "XML::LibXML not installed!!\n";

	
		my $xmlmach = XML::Lite->new("../../../machines/config_machines.xml");
		
	}

}
1;
