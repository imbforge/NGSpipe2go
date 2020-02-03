## Who:  Sergi Sayols
## When: 04-01-2014
our $CWD = $ARGV[0];
our $OUT = $ARGV[1];

## loop over the *Log.final.out files, get the name and construct the statistics
my @data;
opendir(DH,$CWD);
my @f = readdir(DH);
closedir(DH);

foreach my $f (@f) {

	# check the extension
	next unless ($f =~ /\.log$/);

	# read and build the hash to be pushed into the array
	my $reg;	# pointer to a hash
	open(FH,$CWD . $f);
	$reg->{'sample'} = $f =~ s/^Sample_richly_2014_03__//r;
	$reg->{'sample'} = $reg->{'sample'} =~ s/\.cut.*\.log$//r;
	while(<FH>) {
		chomp;
		my @l =  split /:\s/,$_;	# split fields and
		$l[1] =~ s/\t//;		# translate tabs
		$reg->{'processed'} = $l[1] if($l[0] =~ /reads processed/);
		$reg->{'aligned'}  = $l[1] if($l[0] =~ /reads with at least one reported alignment/);
		$reg->{'notaligned'} = $l[1] if($l[0] =~ /reads that failed to align/);
		$reg->{'sampled'} = $l[1] if($l[0] =~ /reads with alignments suppressed due to -m/);
	}
	close($f);

	# push the hash into the array
	push @data,$reg;
}

##
## Print output
##
if($OUT eq "tex") {
	## print header
	print('\begin{table}[h!]',"\n");
	print('\tiny',"\n");
	print('\begin{tabular}{|c|c|c|c|c|}\hline',"\n");
	print('Sample & processed & aligned & not aligned & map to more than 20 loci\\\\\\hline',"\n");

	# print the table
	for my $reg (@data) {
		print($reg->{'sample'} =~ s/_/\\_/rg, ' & ',
			  $reg->{'processed'}, ' & ',
			  $reg->{'aligned'} =~ s/%/\\%/r, ' & ',
			  $reg->{'notaligned'} =~ s/%/\\%/r, ' & ',
			  $reg->{'sampled'} =~ s/%/\\%/r, '\\\\\\hline', "\n");
	}

	## print foot
	print('\end{tabular}',"\n");
	print('\caption{Mapping statistics.}\label{qcstatistic}',"\n");
	print('\end{table}',"\n");
} else {
	## print header
	print('Sample,processed,aligned,not.aligned,more.20.loci',"\n");

	# print the table
	for my $reg (@data) {
		print($reg->{'sample'}, ',',
			  $reg->{'processed'}, ',',
			  $reg->{'aligned'} =~ s/\s.+$//r, ',',
			  $reg->{'notaligned'} =~ s/\s.+$//r, ',',
			  $reg->{'sampled'} =~ s/\s.+$//r, "\n");
	}
}
