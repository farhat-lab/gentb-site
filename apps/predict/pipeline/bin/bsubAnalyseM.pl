use warnings;
use strict;

my $path=shift @ARGV;

use FindBin qw($Bin);
chdir($Bin);

system("Rscript $Bin/TBpredict.R ".'"'."${path}/output/matrix.csv".'"'." >${path}/output/result.json");
my $size= -s "${path}/output/result.json";

if ($size > 0 ) {
    system("python $Bin/../run_feedback.py ${path}");
}
