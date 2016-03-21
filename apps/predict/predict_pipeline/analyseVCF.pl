use warnings;
use strict;

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
	system("bsub -q short -W 2:00 -o ${path}/$name.error -J ".$name.' "'."source /groups/murray/run_pipeline/prepare_environment.sh ; perl /groups/murray/run_pipeline/bin/bsubAnalyseV.pl $name $path".'"');
}

