use warnings;
use strict;

my $file=shift @ARGV;
my $path=shift @ARGV;

use FindBin qw($Bin);
chdir($Bin);

system("perl $Bin/flatAnnotatorVAR.pl ${path}/$file.vcf 15 0.1 PASS");
system("mv ${path}/$file.var ${path}/output");
system("python $Bin/generate_matrix.py ${path}/output");
system("Rscript $Bin/TBpredict.R ".'"'."${path}/output/matrix.csv".'"'." >${path}/output/result.json");
my $size= -s "${path}/output/result.json";

if ($size > 0 ) {
        system("python $Bin/../run_feedback.py ${path}");
}
