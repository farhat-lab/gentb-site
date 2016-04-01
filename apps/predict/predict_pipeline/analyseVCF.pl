use warnings;
use strict;

# Get path of this script.
use FindBin qw($Bin);

my $path=shift @ARGV; #path to where the input fastq files are
system("mkdir ${path}/output");

opendir(DIR, $path);

#search for vcf files in the path
my @files;
@files = grep { /\.vcf$/} readdir(DIR);
closedir(DIR);

#print STDERR "@files\n";

foreach(@files)
{
    (my $name= $_) =~ s/\.vcf//g;
    my $command = "bsub";
    if(not `which bsub`) {
        $command = "$Bin/no_bsub";
    }
    system("$command -q short -W 2:00 -o ${path}/$name.error -J ".$name.' "'."source $Bin/prepare_environment.sh ; perl $Bin/bin/bsubAnalyseV.pl $name $path".'"');
}

