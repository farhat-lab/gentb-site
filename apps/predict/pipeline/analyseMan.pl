use warnings;
use strict;

# Get path of this script.
use FindBin qw($Bin);

my $path=shift @ARGV; #path to where the input fastq files are
system("mkdir ${path}/output");

my $command = "bsub";
if(not `which bsub`) {
    $command = "$Bin/no_bsub";
}
system("$command -q short -W 2:00 -o ${path}/manual.error -J manual \"source $Bin/prepare_environment.sh; perl $Bin/bin/bsubAnalyseM.pl $path\"");

