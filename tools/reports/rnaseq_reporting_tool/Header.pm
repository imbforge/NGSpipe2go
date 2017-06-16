#############################################
##
## This library prints the Latex header and the foot of the report
##
## Who:  Sergi Sayols
## When: 28-01-2014
##
#############################################
use strict;
use warnings;

our (@ISA, @EXPORT_OK);
BEGIN {
	require Exporter;
	@ISA = qw(Exporter);
	@EXPORT_OK = qw(reportHeader reportFoot);  # symbols to export on request
}

sub reportHeader {

	my $title  = $_[0];
	my $author = $_[1];

	print '\documentclass[a4paper,10pt]{article}', "\n";
	print '\usepackage[utf8]{inputenc}', "\n";
	print '\usepackage{graphicx}', "\n";
	print '\usepackage{amsmath}', "\n";
	print '\usepackage{tabularx}', "\n";
	print '\usepackage{multicol}', "\n";
	print '\usepackage{caption}', "\n";
	print '\usepackage{array}', "\n";
	print '\usepackage{tabularx}', "\n";
	print '\usepackage{rotating}', "\n";
	print '\usepackage[a4paper,hmargin=2cm,vmargin=2cm]{geometry}%margins', "\n";
	print '% Title Page', "\n";
	print '\title{', $title, '}' , "\n";
	print '\author{', $author, '}' , "\n";

	print '\begin{document}', "\n";
	print '\maketitle', "\n";

}

sub reportFoot {
	print '\end{document}', "\n";
}

1;
