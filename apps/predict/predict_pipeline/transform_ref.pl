use warnings;
use strict;
use Bio::SeqIO;

my $refFile="/groups/murray/run_pipeline/bin/h37rv.fasta";

$refFile=transform($refFile);
system("stampy.py -G $refFile $refFile.fasta");
system("stampy.py -g $refFile -H $refFile");
system("samtools faidx $refFile.fasta");
system("cp $refFile* /groups/murray/run_pipeline/bin/");

sub transform
{
        my $file=shift @_;
        my @name=split('\.',$file);
        my $in  = Bio::SeqIO->new(-file => "$file",-format => 'Fasta');
        my $out = Bio::SeqIO->new(-file => ">".$name[0]."_transformed.fasta",-format => 'Fasta');
        while( my $seq = $in->next_seq() )
        {
                $out->write_seq($seq);
        }
        return $name[0]."_transformed";
}

