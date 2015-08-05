###################################################
##
## This program generates automatix Latex reports by parsing
## the contents of the secified reports folder
##
## A parser function specific for each is written in a separate module
##
## Who:  Sergi Sayols
## When: 28-01-2014
##
###################################################
use strict;
use warnings;
use Getopt::Long;
use Switch;
use lib '/fsimb/groups/imb-bioinfocf/projects/cf_internal/imb_cf_2013_04_sayols_infrastructure_pipelines/src/';
use reports::Header  qw(reportHeader);
use reports::Header  qw(reportFoot);
use reports::Bustard qw(reportBustard);
use reports::Fastqc  qw(reportFastqc);
use reports::STAR    qw(reportSTAR);
use reports::Bowtie  qw(reportBowtie);

############################
##
## INIT TASKS: Take options from command line
##
############################
my %opt = ();
$opt{BASEDIR} = "./report";
$opt{TITLE}   = "";
$opt{AUTHOR}  = "";
$opt{HELP}    = 0;	# help

# error handle
sub usage {
	my $message = $_[0];
	my $command = $0;	# get program name from command line
	$command =~ s#^.*/##;

	print STDERR (
		$message, 
		"Usage: $command [options]\n" . 
	   	"Options:\n" . 
	   	"  --basedir=/path/to/reports   (default: ./report)\n" . 
	   	"  --title=quoted report title  (default: \'\')\n" . 
	   	"  --author=quoted authors list (default: \'\')\n" . 
	   	"  --help, this help screen"
	);

	die("\n")
}

# parse options
GetOptions(
	"basedir" => \$opt{BASEDIR},
	"title"   => \$opt{TITLE},
	"author"  => \$opt{AUTHOR},
	"help"    => \$opt{HELP}) 
	or usage("Invalid commmand line options.\n");

# validate arguments
die("BASEDIR dir $opt{BASEDIR} does not exist") unless(-d $opt{BASEDIR});
usage("") if($opt{HELP});
usage("Incompatible options.\n") if($opt{HELP} && ($opt{BASEDIR} || $opt{TITLE} || $opt{AUTHOR}));

############################
##
## MAIN PROGRAM
##
############################
opendir(my $dh, $opt{BASEDIR}) || die $!;

reportHeader($opt{TITLE},$opt{AUTHOR});

while(my $d = readdir $dh) {

	next unless(-d $opt{BASEDIR} . '/' . $d);	# skip if file is not a directory
	next if ($d =~ /^\.{1,2}$/);	# skip . and .. directories

	# dispatch the function to parse the directory
	switch ($d) {
		case "SAV"    { reportBustard($opt{BASEDIR} . '/' . $d); }
		case "fastqc" { reportFastqc ($opt{BASEDIR} . '/' . $d); }
		case "STAR"   { reportSTAR   ($opt{BASEDIR} . '/' . $d); }
		case "Bowtie" { reportBowtie ($opt{BASEDIR} . '/' . $d); }
		else { print "Don\'t know what to do with $d\n"; }
	}
}
closedir $dh;

reportFoot();
