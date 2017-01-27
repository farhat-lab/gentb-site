use warnings;
use strict;

my $file=shift @ARGV;
my $refFile=shift @ARGV; #has a fasta extension
my $pairend=shift @ARGV; #1 paired end, 0 single end
my $pex=shift @ARGV; #0 single end, _R|.
my $path=shift @ARGV;

use FindBin qw($Bin);
chdir($Bin);

my @name=split('\.',$file);
my $dataFile=$name[0];

# Remove fasta extension from call
$refFile =~ s/\.fasta$//;

my @cmd;

if ($pairend>0) {
	push @cmd, "stampy.py -g $refFile -h $refFile -o ${path}/$dataFile.sam -f sam -M ${path}/${dataFile}${pex}1.fastq ${path}/${dataFile}${pex}2.fastq;";
} else {
	push @cmd, "stampy.py -g $refFile -h $refFile -o ${path}/$dataFile.sam -f sam -M ${path}/${dataFile}.fastq;";
}
push @cmd, "samtools view -bS ${path}/$dataFile.sam > ${path}/$dataFile.bam";
push @cmd, "samtools sort ${path}/$dataFile.bam ${path}/$dataFile.sorted";
push @cmd, "samtools index ${path}/$dataFile.sorted.bam";
push @cmd, "samtools depth -b $Bin/DR_regions.BED -Q 29 ${path}/$dataFile.sorted.bam >${path}/output/$dataFile.qc";
push @cmd, "Platypus.py callVariants --bamFiles=${path}/$dataFile.sorted.bam --refFile=$refFile.fasta --output=${path}/$dataFile.h37rv.vcf";
push @cmd, "perl $Bin/flatAnnotatorVAR.pl ${path}/$dataFile.h37rv.vcf 15 0.1 PASS";
push @cmd, "mv ${path}/$dataFile.h37rv.var ${path}/output";
push @cmd, "mv ${path}/$dataFile.h37rv.vcf ${path}/output";
push @cmd, "python $Bin/generate_matrix.py ${path}/output";
push @cmd, "Rscript $Bin/TBpredict.R ".'"'."${path}/output/matrix.csv".'"';

my $n=0;
my @steps=('stampy','sam2bam','sort','index','QC','platypus','annotate','mvVar','mvVCF','genMatrix','Rpredict');
foreach my $cmd (@cmd) {
	my $i=system($cmd);
	if ($i>0) {
		die "died at $steps[$n] with error $i";
	}
	$n++;
}

system("python $Bin/../run_feedback.py ${path}");

system("rm -f ${path}/$dataFile.sam");
system("rm -f ${path}/$dataFile.sorted.bam");
system("rm -f ${path}/$dataFile.sorted.bam.bai");
system("rm -f ${path}/$dataFile.bam");
