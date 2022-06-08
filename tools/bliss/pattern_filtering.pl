## pattern_filtering.pl
##
## What: Match and trim a fastq file using a perl regex.
##       Incorporate the UMI (if any) into the read name.
## Who: Sergi Sayols
## When: 20-06-2022
##
## The script keeps + trims reads (==4 lines) in a fastq file matching a provided regex.
## If the regex contains a UMI, this is added to the read name.
## The expected call should be something like:
##   $ perl pattern_filtering.pl -1 x.R1.fastq.gz \             # R1 compressed fastq
##                               -2 x.R2.fastq.gz \             # R2 compressed fastq (optional)
##                               -o out \                       # suffix for the output file(s) name
##                               -r '^(........)(.TCACACGT' \   # perl regex to match+trim
##                               -u '$1'                        # group from the perl regex that contains the UMI
##
use strict; 
use warnings; 

# parse command line arguments
# expect: perl pattern_filtering.pl -1 x.r1.fastq.gz -2 x.r2.fastq.gz -o out_suffix -r regex -u umi_group
use Getopt::Std;
my %opts;
getopts('1:2:o:r:u:d', \%opts);

if($opts{d}) {  # debug enabled
	print "read 1    : $opts{1}\n";
	print "read 2    : $opts{2}\n";
	print "out prefix: $opts{o}\n";
	print "pattern   : $opts{r}\n";
	print "umi       : $opts{u}\n";
}

# open R1 for reading/writing
my ($R1, $R2, $R1_out, $R2_out);
if($opts{1} =~ /\.gz$/) {
	open($R1, "gunzip -c $opts{1} |") or die "gunzip $opts{1} : $!";
} else {
	open($R1, "<", $opts{1}) or die $!;
}
my $x = $opts{1};
$x =~ s/\.fastq(\.gz|)$/.$opts{o}.fastq$1/;
if($opts{1} =~ /\.gz$/) {
	open($R1_out, "| gzip -c > $x") or die "gzip $x : $!";
} else {
	open($R1_out, ">", $x) or die $!;
}

# process input files and write output if it matches the regular expression
my $matches=0;
my $reads=0;

# PAIRED END
if($opts{2}) {
	# open R2 for reading/writing
	if($opts{2} =~ /\.gz$/) {
		open($R2, "gunzip -c $opts{2} |") or die "gunzip $opts{2} : $!";
	} else {
		open($R2, "<", $opts{2}) or die $!;
	}
	my $x = $opts{2};
	$x =~ s/\.fastq(\.gz|)$/.$opts{o}.fastq$1/;
	if($opts{2} =~ /\.gz$/) {
		open($R2_out, "| gzip -c > $x") or die "gzip $x : $!";
	} else {
		open($R2_out, ">", $x) or die $!;
	}

	# process (match RE + write)
	while(!eof($R1) and !eof($R2)) {
		$reads++;
		chomp(my @r1 = (my $header=<$R1>, my $seq=<$R1>, my $plus=<$R1>, my $qual=<$R1>));
		chomp(my @r2 = ($header=<$R2>, $seq=<$R2>, $plus=<$R2>, $qual=<$R2>));

		# if pattern is in read
		if($r1[1] =~ s/$opts{r}//) {
			$matches++;
			# add umi to the end of the read name
			if($opts{u}) {
				my $umi = eval $opts{u};      # eval is expensive. Use wisely...
				$r1[0] =~ s/(\s|\n)/_$umi$1/;   # add the UMI after the first whitespace (or at the end if no whitespace)
				$r2[0] =~ s/(\s|\n)/_$umi$1/;   # add the UMI after the first whitespace (or at the end if no whitespace)
			}
			# remove pattern from qual string
			$r1[3] = substr($r1[3], 0, $-[0]) . substr($r1[3], $+[0]);
			# write output
			print $R1_out join("\n", @r1),"\n";
			print $R2_out join("\n", @r2),"\n";
		}
	}

	# close files
	close $R2;
	close $R2_out;
}
# SINGLE END
else {
	# process (match RE + write)
	while(!eof($R1)) {
		$reads++;
		chomp(my @r1 = (my $header=<$R1>, my $seq=<$R1>, my $plus=<$R1>, my $qual=<$R1>));
		# if pattern is in read
		if($r1[1] =~ s/$opts{r}//) {
			$matches++;
			# add umi to the end of the read name
			if($opts{u}) {
				my $umi = eval $opts{u};      # eval is expensive. Use wisely...
				$r1[0] =~ s/(\s|\n)/_$umi$1/;   # add the UMI after the first whitespace (or at the end if no whitespace)
			}
			# remove pattern from qual string
			$r1[3] = substr($r1[3], 0, $-[0]) . substr($r1[3], $+[0]);
			# write output
			print $R1_out join("\n", @r1),"\n";
		}
	}
}

print "Found ", $matches, " out of ", $reads, " reads matching ", $opts{r}, "\n";

# close files
close $R1;
close $R1_out;
