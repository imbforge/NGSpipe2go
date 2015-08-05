## Who:  Sergi Sayols
## When: 04-01-2014
our $CWD = $ARGV[0];

## print header
print('{\tiny',"\n");
print('\begin{longtable}{@{}cccc@{}}',"\n");
print('\centering',"\n");
print('& Duplication & read quals & sequence bias\\\\\\hline',"\n");

## loop over the *_fastqc.zip files, get the name and construct the statistics
opendir(DH,$CWD);
my @f = readdir(DH);
closedir(DH);

foreach my $f (@f) {

	# check the extension
	next unless ($f =~ /_fastqc.zip$/);

	# read and build the hash to be pushed into the array
	my $s = $f =~ s/_fastqc.zip$//r;
	system("unzip $CWD/$f ". $s ."_fastqc/Images/duplication_levels.png -d $CWD/ -qq -o >> FASTQC.log");
	die("error unzipping $f") if($? > 1);
	system("unzip $CWD/$f ". $s ."_fastqc/Images/per_base_quality.png -d $CWD/ -qq -o >> FASTQC.log");
	die("error unzipping $f") if($? > 1);
	system("unzip $CWD/$f ". $s ."_fastqc/Images/per_base_sequence_content.png -d $CWD/ -qq -o >> FASTQC.log");
	die("error unzipping $f") if($? > 1);

	# print the latex row
	print($s =~ s/_/\\_/rg," &\n",
	  	  "\\includegraphics[width=.23\\textwidth]{".$CWD."/".$s."_fastqc/Images/duplication_levels.png} &\n",
	  	  "\\includegraphics[width=.23\\textwidth]{".$CWD."/".$s."_fastqc/Images/per_base_quality.png} &\n",
	  	  "\\includegraphics[width=.23\\textwidth]{".$CWD."/".$s."_fastqc/Images/per_base_sequence_content.png} \\\\\n");
}

## print foot
print('\caption{FASTQC statistics.}\label{FASTQC}',"\n");
print('\end{longtable}}',"\n");
