use warnings;
use strict;

my $path=shift @ARGV;

use FindBin qw($Bin);
chdir($Bin);

my $i=system("Rscript $Bin/TBpredict.R ".'"'."${path}/output/matrix.csv".'"');

if ($i>0) {
	die "1";
}

system("python $Bin/../run_feedback.py ${path}");

