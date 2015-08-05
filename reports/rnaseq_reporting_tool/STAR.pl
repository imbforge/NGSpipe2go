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
	next unless ($f =~ /.+\.Log\.final\.out$/);

	# read and build the hash to be pushed into the array
	my $reg;	# pointer to a hash
	open(FH,$CWD . $f);
	$reg->{'sample'} = $f =~ s/\.Log\.final\.out$//r;
	while(<FH>) {
		chomp;
		my @l =  split /\|/,$_;	# split fields and
		$l[1] =~ s/\t//;  # translate tabs
		$reg->{'input'} = $l[1] if($l[0] =~ /Number of input reads/);
		$reg->{'uniq'}  = $l[1] if($l[0] =~ /Uniquely mapped reads number/);
		$reg->{'uniqP'} = $l[1] =~ s/%//r if($l[0] =~ /Uniquely mapped reads %/);
		$reg->{'multi1'} = $l[1] if($l[0] =~ /Number of reads mapped to multiple loci/);
		$reg->{'multi2'} = $l[1] if($l[0] =~ /Number of reads mapped to too many loci/);
		$reg->{'multiP1'}= $l[1] =~ s/%//r if($l[0] =~ /% of reads mapped to multiple loci/);
		$reg->{'multiP2'}= $l[1] =~ s/%//r if($l[0] =~ /% of reads mapped to too many loci/);
		$reg->{'unmap1'}= $l[1] =~ s/%//r if($l[0] =~ /% of reads unmapped: too many mismatches/);
		$reg->{'unmap2'}= $l[1] =~ s/%//r if($l[0] =~ /% of reads unmapped: too short/);
		$reg->{'unmap3'}= $l[1] =~ s/%//r if($l[0] =~ /% of reads unmapped: other/);
	}
	close(FH);

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
	print('Sample & input reads & uniq mapped & multimapped & unmapped\\\\\\hline',"\n");

	## print the records
	for my $reg (@data) {
		print($reg->{'sample'} =~ s/_/\\_/rg, ' & ',
			  $reg->{'input'} , ' & ',
			  $reg->{'uniq'}  , ' (', $reg->{'uniqP'} , '\\%)', ' & ',
			  $reg->{'multi1'} + $reg->{'multi1'}, ' (', $reg->{'multiP1'} + $reg->{'multiP2'}, '\\%)', ' & ',
			  $reg->{'unmap1'} + $reg->{'unmap2'} + $reg->{'unmap3'}, '\\%\\\\\\hline', "\n");
	}

	## print foot
	print('\end{tabular}',"\n");
	print('\caption{Mapping statistics.}\label{qcstatistic}',"\n");
	print('\end{table}',"\n");
} else {
	## print header
	print('Sample,reads,unique,uniqueP,multimapped,multiP,unmapped',"\n");

	# print the table
	for my $reg (@data) {
		print($reg->{'sample'}, ',',
			  $reg->{'input'},  ',',
			  $reg->{'uniq'},   ',',
			  $reg->{'uniqP'},  '%,',
			  $reg->{'multi1'}+$reg->{'multi2'},  ',',
			  $reg->{'multiP1'}+$reg->{'multiP2'}, '%,',
			  $reg->{'unmap1'}+$reg->{'unmap2'}+$reg->{'unmap3'},"\%\n");
	}
}
