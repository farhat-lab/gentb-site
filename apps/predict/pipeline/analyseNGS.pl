use warnings;
use strict;

# Get path of this script.
use FindBin qw($Bin);

my $pairend= shift @ARGV; #1 if two files per strain are given or 0 if only one file per strain is given
my $pex=shift @ARGV; #_R or . this is the pair-end extension
my $path=shift @ARGV; #path to where the input fastq files are
system("mkdir ${path}/output");


my $ref="$Bin/bin/h37rv_transformed.fasta";
unless (-f $ref) {
  system("bsub -q short -J transform_rv -W 5:00 -o error "."\"source $Bin/prepare_environment.sh; perl $Bin/transform_ref.pl\"");
}

opendir(DIR, $path);

#loop that searches for fastq files in the path
my @files;
if ($pairend >0) {
        if ($pex eq ".") {
		@files =  grep { /\.1\.fastq\.?g?z?$/ } readdir(DIR);
	} else {
		@files = grep { /${pex}1\.fastq\.?g?z?$/} readdir(DIR);
	}
} else {
	@files = grep { /\.fastq\.?g?z?$/} readdir(DIR);
}
closedir(DIR);

foreach(@files)
{
	my @name=split('\.',$_);
	if ($pex ne ".") {
		@name=split($pex,$name[0]);
	}
        my $command = "bsub";
        if(not `which bsub`) {
            $command = "$Bin/no_bsub";
        }
	system("$command -q long -W 99:00 -o ${path}/$name[0].error -J ".$name[0].' "'."source $Bin/prepare_environment.sh ; perl $Bin/bin/bsubAnalyse.pl $name[0] $ref $pairend $pex $path".'"');
}

