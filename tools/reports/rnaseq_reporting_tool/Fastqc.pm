#############################################
##
## This library converts the Bustard XML data file to Latex
##
## Who:  Sergi Sayols
## When: 27-01-2014
##
#############################################
use strict;
use warnings;

our (@ISA, @EXPORT_OK);
BEGIN {
	require Exporter;
	@ISA = qw(Exporter);
	@EXPORT_OK = qw(reportFastqc);  # symbols to export on request
}

# loop over the files and print the Latex report
sub reportFastqc {

	my $dir = $_[0];
	my $templateFileName = quotemeta ".+\.zip";

	opendir(DIR, $dir) or die $!;
	for(my $f = readdir(DIR)) {

		# ignore files not matching the pattern
		print parseFastqcFile($f) if ($f =~ m/$templateFileName/);
	}
}

sub parseFastqcFile {
	return "";
}

1;
