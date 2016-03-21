use warnings;
use strict;

my $file=shift @ARGV;
my $path=shift @ARGV;

system("perl /groups/murray/run_pipeline/bin/flatAnnotatorVAR.pl ${path}/$file.vcf 15 0.1 PASS");
system("mv ${path}/$file.var ${path}/output");
system("python /groups/murray/run_pipeline/bin/generate_matrix.py ${path}/output");
system("Rscript /groups/murray/run_pipeline/bin/TBpredict.R ".'"'."${path}/output/matrix.csv".'"'." >${path}/output/result.json");
my $size= -s ${path}/output/result.json
if ($size > 0 ) {
        system("python ${path}/gentb_status_feedback.py");
}
