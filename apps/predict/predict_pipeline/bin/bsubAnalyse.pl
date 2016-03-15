use warnings;
use strict;

my $file=shift @ARGV;
my $refFile=shift @ARGV; #has a fasta extension
my @name=split('\.',$refFile); 
$refFile=$name[0]; #fasta extension removed, assumes path doesnot contain a period
my $pairend=shift @ARGV; #1 paired end, 0 single end
my $pex=shift @ARGV; #0 single end, _R|.
my $path=shift @ARGV;

@name=split('\.',$file);
my $dataFile=$name[0];

if ($pairend>0) {
	system("stampy.py -g $refFile -h $refFile -o ${path}/$dataFile.sam -f sam -M ${path}/${dataFile}${pex}1.fastq ${path}/${dataFile}${pex}2.fastq;");
	#print STDERR "stampy.py -g $refFile -h $refFile -M ${path}/${dataFile}${pex}1.fastq ${path}/${dataFile}${pex}2.fastq -o ${path}/$dataFile.sam -f sam";
} else {
	system("stampy.py -g $refFile -h $refFile -o ${path}/$dataFile.sam -f sam -M ${path}/${dataFile}.fastq;");
	#print STDERR "stampy.py -g $refFile -h $refFile -M ${path}/${dataFile}.fastq -o ${path}/$dataFile.sam -f sam";
}
system("samtools view -bS ${path}/$dataFile.sam > ${path}/$dataFile.bam");
system("samtools sort ${path}/$dataFile.bam ${path}/$dataFile.sorted");
system("samtools index ${path}/$dataFile.sorted.bam");
system("Platypus.py callVariants --bamFiles=${path}/$dataFile.sorted.bam --refFile=$refFile.fasta --output=${path}/$dataFile.h37rv.vcf");
system("perl /groups/murray/run_pipeline/bin/flatAnnotatorVAR.pl ${path}/$dataFile.h37rv.vcf 15 0.1 PASS");
system("mv ${path}/$dataFile.h37rv.var ${path}/output");
system("mv ${path}/$dataFile.h37rv.vcf ${path}/output");
system("python /groups/murray/run_pipeline/bin/generate_matrix.py ${path}/output");
system("Rscript /groups/murray/run_pipeline/bin/TBpredict.R ".'"'."${path}/output/matrix.csv".'"'." >${path}/output/result.json");
my $size= -s "${path}/output/result.json";
#print STDERR "$size\n";
if ($size > 0 ) {
        system("python ${path}/gentb_status_feedback.py");
}
system("rm ${path}/$dataFile.sam");
system("rm ${path}/$dataFile.sorted.bam");
system("rm ${path}/$dataFile.sorted.bam.bai");
system("rm ${path}/$dataFile.bam");
