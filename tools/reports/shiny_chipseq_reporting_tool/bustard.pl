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
use XML::LibXSLT;
use XML::LibXML;

my $f = $ARGV[0];
my $XSLT = "BustardSummary.toMD.xsl";

print parseBustardFile($XSLT,$f);

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
