package SeqCoord;
use strict;
use warnings;
use Exporter;

our @ISA= qw( Exporter );

# these CAN be exported.
our @EXPORT_OK = qw( load_fasta get_coord );

# these are exported by default.
our @EXPORT = qw( load_fasta get_coord );

our $seq = '';
our $len = 0;

sub load_fasta {
    my $filename = shift @_;
    open(my $fh => $filename) || die "Cannot open $filename: $!";
    while(<$fh>) {
        chomp;
        next if /^>/;
        $seq .= $_;
        $len += length($_);
    }
    warn "Loaded $len bytes of data from fasta";
    close($fh);
}

sub get_coord {
    my $maxr = 0;
    my $sense;
    my @cr;
    for my $pair (@_) {
        my ($left, $right) = @$pair;
        $sense = ($left <= $right);
        my @newpair = $sense ? ($left, $right) : ($right,$left);
        push @cr, \@newpair;
        $maxr = $newpair[1] if ($maxr < $newpair[1]);
    }

    my $n=1;
    my $sequence;
    foreach my $exon (@cr) {
        my $eseq = substr $seq, $exon->[0] - 1, $exon->[1] - $exon->[0] + 1;
        $eseq = revcomp($eseq) unless $sense;
        $sequence .= $eseq;
        $n++;
    }
    my $ret = unpack "A60" x (int(length($sequence)/60)+1),$sequence;
    return $ret;
}

sub revcomp{
    my $seq = $_[0];
    $seq = join ('', reverse (unpack 'A' x length($_[0]), $_[0]));
    # print $seq;
    $seq =~ tr/actgACTG/tgacTGAC/;
    return $seq
}

1;

