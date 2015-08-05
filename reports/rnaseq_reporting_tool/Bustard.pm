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
	use XML::LibXSLT;
	use XML::LibXML;
	package report::xmlparser;

	require Exporter;
	push @ISA,qw(Exporter);
	push @EXPORT_OK,qw(reportBustard);  # symbols to export on request
}

# loop over the files and print the Latex report
sub reportBustard {

	my $dir = $_[0];
	my $templateFileName = "^BustardSummary.*\.xml";
	my $XSLT = "BustardSummary.toLatex.xsl";

	opendir(DIR, $dir) or die $!;
	for my $f (readdir(DIR)) {

		# ignore files not matching the pattern
		print parseBustardFile($XSLT,$f) if ($f =~ m/\Q$templateFileName/);
	}
}

sub parseBustardFile {
	# the arguments for this command are stylesheet and source files
	my ($style_file,$source_file) = @_;

	die("XSLT file $style_file doesnt exist\n") if(! -e $style_file);
	die("XML file $source_file doesnt exist\n") if(! -e $source_file);

	# initialize the parser and XSLT processor
	my $parser = XML::LibXML->new();
	my $xslt = XML::LibXSLT->new();
	my $stylesheet = $xslt->parse_stylesheet_file($style_file);

	# for each source file: parse, transform, print out result
	my $source_doc = $parser->parse_file($source_file);
	my $result = $stylesheet->transform($source_doc);
	return $stylesheet->output_string($result);
}

1;
