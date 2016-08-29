use warnings;
use strict;

my $file=shift @ARGV;
my $path=shift @ARGV;

use FindBin qw($Bin);
chdir($Bin);

my @cmd;

push @cmd, "perl $Bin/flatAnnotatorVAR.pl ${path}/$file.vcf 15 0.1 PASS";
push @cmd, "mv ${path}/$file.var ${path}/output";
push @cmd, "python $Bin/generate_matrix.py ${path}/output";
push @cmd, "Rscript $Bin/TBpredict.R ".'"'."${path}/output/matrix.csv".'"';

my $n=0;
my @steps=('annotate','mvVar','genMatrix','Rpredict');
foreach my $cmd (@cmd) {
	my $i=system($cmd);
	if ($i>0) {
		die "died at $steps[$n] with error $i";
	}
	$n++;
}

system("python $Bin/../run_feedback.py ${path}");
