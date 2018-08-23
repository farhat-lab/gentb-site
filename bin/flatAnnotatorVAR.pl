#!/usr/bin/perl
use strict;
use warnings;
### script that takes in .vcf file and produces a .var file. Filters and combines the mutation data in .vcf file with data from 2 genome coordinate files (with headers) #####(one for coding ###regions and one for non coding regions) to add functional data. Also requires h37rv.fasta file (there reference genome sequence file) and ###get_seq_coord.pl script to exist ###in the same folder

### example command ./flatAnnotatorVAR.pl test.vcf qual{0-255} hetero{0-1} platypusfilter{PASS|ALL} (PASS now includes badReads output will be test.var) 
use FindBin qw($Bin);


$/="\n";
$,="\t";
$\="\n";
my %tt11;
&create_translation_table11(\%tt11);

my $reference = "h37rv"; #this may need to be changed in the future
my $ref_file = shift @ARGV;
my $Creference = shift @ARGV;
my $Nreference = shift @ARGV;
my $snpfile = shift @ARGV;
#my $refPath= shift @ARGV;
my $qualThresh = (shift@ARGV)||0;
my $heteroThresh = (shift@ARGV)||0;
my $platypusFilter = (shift@ARGV)||'PASS';
if ($platypusFilter =~ m/pass/i) {
  $platypusFilter = 'PASS|badReads|alleleBias';
} else {
  #$platypusFilter =~ s/\-/\|/g;
  #$platypusFilter = "qr/(?!\A.*$platypusFilter.*\z)/"; 
  $platypusFilter = ".+?"; #match anything
}
#print STDERR $platypusFilter;
my $aasnp;
(my $strain) = ($snpfile =~ m/(strain\d+)/);
my $DEBUG=1;

#print STDERR "reading file $Creference\n";
open REF, "<$Creference" or die;
my $n=1;
my %index;
my %start;
my %end;
my %symbol;
my %strand;
my %desc;
while (<REF>) {
 chomp;
 my @w=split/\t/;
 if ($n==1) {
  foreach my $i (0..$#w) {
   if ($w[$i] =~ m/start/i) {
    $index{"start"}=$i;
    #print STDERR "start index is $i\n";
   } elsif ($w[$i] =~ m/end/i) {
    $index{"end"}=$i;
    #print STDERR "end index is $i\n";
   } elsif ($w[$i] =~ m/symbol/i) {
    $index{"symbol"}=$i;
    #print STDERR "symbol index is $i\n";
   } elsif ($w[$i] =~ m/number|name/i) {
    $index{"name"}=$i;
    #print STDERR "name index is $i\n";
   } elsif ($w[$i] =~ m/strand/i) {
    $index{"strand"}=$i;
    #print STDERR "strand index is $i\n";
   } elsif ($w[$i] =~ m/desc/i) {
    $index{"desc"}=$i;
    #print STDERR "desc index is $i\n";
   }
  }
 } else {
  $start{$w[$index{"name"}]}=$w[$index{"start"}];
  $end{$w[$index{"name"}]}=$w[$index{"end"}];
  $w[$index{"symbol"}] =~ s/_/-/g;
  $symbol{$w[$index{"name"}]}=$w[$index{"symbol"}]||'\N';
  $strand{$w[$index{"name"}]}=$w[$index{"strand"}];
  $desc{$w[$index{"name"}]}=$w[$index{"desc"}];
 }
 $n++;
}

#print STDERR "reading file $Nreference\n";
open REF, "<$Nreference" or die;
$n=1;
my %type;
my %geneBefore;
my %geneAfter;
while (<REF>) {
 chomp;
 my @w=split/\t/;
 if ($n==1) {
  foreach my $i (0..$#w) {
   if ($w[$i] =~ m/start/i) {
    $index{"start"}=$i;
   } elsif ($w[$i] =~ m/end/i) {
    $index{"end"}=$i;
   } elsif ($w[$i] =~ m/type/i) {
    $index{"type"}=$i;
   } elsif ($w[$i] =~ m/post/i) {
    $index{"geneBefore"}=$i;
   } elsif ($w[$i] =~ m/strand/i) {
    $index{"strand"}=$i;
   } elsif ($w[$i] =~ m/desc/i) {
    $index{"desc"}=$i;
   } elsif ($w[$i] =~ m/pre/i) {
    $index{"geneAfter"}=$i;
   } elsif ($w[$i] =~ m/region/i) {
    $index{"name"}=$i;
   }
  }
 } else {
  $start{$w[$index{"name"}]}=$w[$index{"start"}];
  $end{$w[$index{"name"}]}=$w[$index{"end"}];
  $type{$w[$index{"name"}]}=$w[$index{"type"}]||"";
  $strand{$w[$index{"name"}]}=$w[$index{"strand"}];
  $w[$index{"desc"}] =~ s/_/-/g;
  $desc{$w[$index{"name"}]}=$w[$index{"desc"}];
  $geneBefore{$w[$index{"name"}]}=$w[$index{"geneBefore"}]; #print STDERR $w[$index{"geneBefore"}]; 
  $geneAfter{$w[$index{"name"}]}=$w[$index{"geneAfter"}];
  $symbol{$w[$index{"name"}]}='\N';
 }
 $n++;
}
	
#print STDERR "reading file $snpfile" if $DEBUG;
open IN,"<$snpfile";
my $asnpfile  = $snpfile;
#$asnpfile =~ s/\.vcf$/.var/;
#open SNP,">$asnpfile";
print STDOUT "reference\tregionid\tgenesymbol\tvarname\tdesc\tqual\tdepth\tbidir\thqr\thqr_ref\tfq\tcodon\taltcodon\tcodpos\tplatypusfilter";
while (<IN>) {
 next if m/^#/;
 chomp;
 my $line=$_;
 my @f=split/\t/;
 my ($from,$ref_allele,$allele,$snpqual,$filter,$info)=($f[1],@f[3..7]);
 (my $depth) = ($info =~ /TC=(\d+)/);
 (my $tcf)   = ($info =~ /TCF=(\d+)/);
 (my $tcr)   = ($info =~ /TCR=(\d+)/);
 (my $nf)    = ($info =~ /NF=(\d+)/);
 (my $nr)    = ($info =~ /NR=(\d+)/);
 my $bidir = '';
 my ($dpr1,$dpr2,$dp1,$dp2) = ($tcf-$nf, $tcr-$nr, $nf, $nr);
 my ($fq) = ($info =~ /RMSMQ=([-.\d]+)/)||'\N';
 $bidir = ($dp1 && $dp2)?'Y':'N';
 next if $ref_allele eq 'N'; #ambiguous reference
 next if $allele =~ m/,/; #heterogenous allele
 #my $cnc = 'N'; #coding/non-coding = C (coding) or N (non-coding)
 #my $sns = '';  #for coding SNPs: S (synonymous) or N (non-synonymous); for everything else ''=undefined
 #my $snpstrand = '.'; # for coding SNPs: '+' or '-'; for everything else '.' = undefined; 
 $aasnp='';
 my $hqr_var = $dp1 + $dp2;
 my $hqr_ref = $dpr1 + $dpr2;
 if ($hqr_var == 0 && $hqr_var ==0) { # to protect against division by zero error
 	$hqr_ref=1;
 }
 if ($filter =~ m/;/) {
	#print STDERR $filter;
 } else {
  if ($snpqual >$qualThresh && $hqr_var/($hqr_var+$hqr_ref) >$heteroThresh && $filter =~ m/$platypusFilter/ ) { 
   my ($genename,$codon,$codpos,$altcodon,$name);#$nucpos
   my ($type,$before_gene,$after_gene);
   if (length($ref_allele)==length($allele)) {#SNP
    my $march=$from;
    while ($march<=$from+length($allele)-1) {
     while (my ($gname, $genomicPos) = each %start) {
      if ($genomicPos <$march && $end{$gname}>=$march) { #for zero based coordinates
        $genename=$gname;
      }
     }
     if ($genename =~ /^Rv/i) { #coding
      #print STDERR "march is $march before, $genename, $ref_allele, $allele, $from";
      ($name, $codon, $altcodon, $codpos, $march)=&assign_allele($genename, $ref_allele, $allele, $march, $from);
      #print STDERR "$ name, $codon, $altcodon, $codpos, march is $march after";
     } else { ## 1 noncoding variant
      my $ref_base=substr($ref_allele,$march-$from,1);
      my $base=substr($allele,$march-$from,1);
      ($name, $codon, $altcodon, $codpos)=&annotnoncoding($march, $genename,$ref_base, $base);
     }
     print STDOUT $reference,$genename,$symbol{$genename},$name,$desc{$genename},$snpqual,$depth,$bidir,$hqr_var,$hqr_ref,$fq,$codon,$altcodon,$codpos,$filter;    
     $march++;
    }
   } elsif (length($allele)!=length($ref_allele)) { #indel may or may not also have SNP
    if ($allele eq '.') {
     $allele='';
    } elsif ($ref_allele eq ".") {
     $ref_allele='';
    }
    if (length($ref_allele)>0 && substr($allele,0,length($ref_allele)) eq $ref_allele) {
     $allele=substr($allele,length($ref_allele));
     $ref_allele='';
     $from++;
    }
    if (length($allele)>0 && substr($ref_allele,0,length($allele)) eq $allele) {
     $ref_allele=substr($ref_allele,length($allele));
     $allele='';
     $from++;
    }
    if (length($allele)>0 && length($ref_allele)>0) {
     my $march=$from;
     if (length($allele)>length($ref_allele)) { #at least partly an insertion
      while ($march <=$from+length($ref_allele)-1) {
       while (my ($gname, $genomicPos) = each %start) {
        if ($genomicPos <$march && $end{$gname}>=$march) { #for zero based coordinates
         $genename=$gname;
        }
       }
       my $ref_base=substr($ref_allele,$march-$from,1);
       my $base=substr($allele,$march-$from,1);
       if ($genename =~ /^Rv/i) { # 1 gene
        ($name, $codon, $altcodon, $codpos)=&assign_allele($genename, $ref_allele, $allele, $march, $from);
       } else { ## 1 noncoding variant
        ($name, $codon, $altcodon, $codpos)=&annotnoncoding($march, $genename,$ref_base, $base);
       }
       print STDOUT $reference,$genename,$symbol{$genename},$name,$desc{$genename},$snpqual,$depth,$bidir,$hqr_var,$hqr_ref,$fq,$codon,$altcodon,$codpos,$filter;
       $march++;
      } 
     } elsif (length($ref_allele)>length($allele)) { #at least partly a deletion
      while ($march <=$from+length($allele)-1) {
       while (my ($gname, $genomicPos) = each %start) {
        if ($genomicPos <$march && $end{$gname}>=$march) { #for zero based coordinates
         $genename=$gname;
        }
       }
       my $ref_base=substr($ref_allele,$march-$from,1);
       my $base=substr($allele,$march-$from,1);
       if ($genename =~ /^Rv/i) { # 1 gene
        ($name, $codon, $altcodon, $codpos, $march)=&assign_allele($genename, $ref_allele, $allele, $march, $from);
       } else { ## 1 noncoding variant
        ($name, $codon, $altcodon, $codpos)=&annotnoncoding($march, $genename,$ref_base, $base);
       }
       print STDOUT $reference,$genename,$symbol{$genename},$name,$desc{$genename},$snpqual,$depth,$bidir,$hqr_var,$hqr_ref,$fq,$codon,$altcodon,$codpos,$filter;
       $march++;
      } 
     } #completes SNP porition of INS or DEL $march variable marks where we are
     while (my ($gname, $genomicPos) = each %start) {
      if ($genomicPos <$march && $end{$gname}>=$march) { #for zero based coordinates
       $genename=$gname;
      }
     }
     ($name,$codon,$altcodon,$codpos) = &annotindel($from,$genename,$ref_allele, $allele);
     print STDOUT $reference,$genename,$symbol{$genename},$name,$desc{$genename},$snpqual,$depth,$bidir,$hqr_var,$hqr_ref,$fq,$codon,$altcodon,$codpos,$filter; 
    } else { #completes INDEL portion of combo SNP/indel
     while (my ($gname, $genomicPos) = each %start) {
      if ($genomicPos <$from && $end{$gname}>=$from) { #for zero based coordinates
        $genename=$gname;
      }
     }
     ($name,$codon,$altcodon,$codpos) = &annotindel($from,$genename,$ref_allele, $allele);
     print STDOUT $reference,$genename,$symbol{$genename},$name,$desc{$genename},$snpqual,$depth,$bidir,$hqr_var,$hqr_ref,$fq,$codon,$altcodon,$codpos,$filter;
    }
   } #close indel
  } #close high Quality
 } #close no semicolon
} #vcf file complete

#close SNP;
close IN;

sub assign_allele{
 my $genename= shift @_;
 my $ref_allele=shift @_;
 my $allele=shift @_;
 my $march=shift @_;
 my $from =shift @_;
 my $cp; my $name; my $codon; my $altcodon; my $codpos;
 if ($strand{$genename} eq "+") {
  $cp = ($march - $start{$genename}) % 3; #within the codon 1,2
  $cp = 3 unless $cp; # 0->3;
  #print STDERR "$cp", length($allele);
  if ($cp == 3 || length($allele)== 1) {
   my $ref_base=substr($ref_allele,$march-$from,1);
   my $base=substr($allele,$march-$from,1);
   ($name, $codon, $altcodon, $codpos)=&annotcoding($march, $genename, $ref_base, $base);
  } elsif ($cp != 3 && length($allele)>1) {
   if (($cp == 1 && length($allele)<=3) || ($cp ==2 && length($allele)<=2)) {
    #print STDERR "HERE!!";
    my $ref_base=substr($ref_allele,$march-$from, length($allele));
    my $base=substr($allele,$march-$from,length($allele));
    ($name, $codon, $altcodon, $codpos)=&annotcoding($march, $genename,$ref_base, $base);
    $march=$march+length($allele)-1;
   } elsif ($cp == 1 && length($allele)>3) {
    my $ref_base=substr($ref_allele,$march-$from, 3);
    my $base=substr($allele,$march-$from, 3);
    ($name, $codon, $altcodon, $codpos)=&annotcoding($march, $genename,$ref_base, $base);
    $march=$march+2;
   } elsif ($cp == 2 && length($allele)>2) {
    my $ref_base=substr($ref_allele,$march-$from, 2);
    my $base=substr($allele,$march-$from,2);
    ($name, $codon, $altcodon, $codpos)=&annotcoding($march, $genename,$ref_base, $base);
    $march++;
   }
  }
 } elsif ($strand{$genename} eq "-") {
  $cp = ($end{$genename} - $march + 1) % 3;
  $cp = 3 unless $cp; # 0->3;
  if ($cp == 1 || length($allele)== 1) {  #codon count for negative strand is 1,3,2,1,3,2
   my $ref_base=substr($ref_allele,$march-$from,1);
   my $base=substr($allele,$march-$from,1);
   ($name, $codon, $altcodon, $codpos)=&annotcoding($march, $genename, $ref_base, $base);
  } else { #$cp!=1 and length($allele)>1
   if (($cp == 3 && length($allele)<=3) || ($cp ==2 && length($allele)<=2)) {
    #print STDERR "HERE!!";
    my $ref_base=substr($ref_allele,$march-$from, length($allele));
    my $base=substr($allele,$march-$from,length($allele));
    ($name, $codon, $altcodon, $codpos)=&annotcoding($march, $genename,$ref_base, $base);
    $march=$march+length($allele)-1;
   } elsif ($cp == 3 && length($allele)>3) {
    my $ref_base=substr($ref_allele,$march-$from, 3);
    my $base=substr($allele,$march-$from, 3);
    ($name, $codon, $altcodon, $codpos)=&annotcoding($march, $genename,$ref_base, $base);
    $march=$march+2;
   } elsif ($cp == 2 && length($allele)>2) {
    my $ref_base=substr($ref_allele,$march-$from, 2);
    my $base=substr($allele,$march-$from,2);
    ($name, $codon, $altcodon, $codpos)=&annotcoding($march, $genename,$ref_base, $base);
    $march++;
   }
  }
 }
 return ($name, $codon, $altcodon, $codpos, $march);
}

sub check_synonymous{
 my $codon = shift @_;
 my $aapos = shift @_;
 my $altcodon = shift @_;
 my $aaref; my $aavar; 
 foreach my $codnmbr (1..length($codon)/3) {
  my @bases=split //, $codon;
  my @altbases=split //,$altcodon;
  #print STDERR (@bases)[($codnmbr-1)*3..(($codnmbr-1)*3+2)];
  $aaref.=$tt11{join('',(@bases)[($codnmbr-1)*3..(($codnmbr-1)*3+2)])};
  $aavar.=$tt11{join('',(@altbases)[($codnmbr-1)*3..(($codnmbr-1)*3+2)])};
 }
 if ($aaref eq $aavar) {
  return ('S', $aaref, $aavar); #synonymous
 }else{
  return ('N', $aaref, $aavar); #non-synonymous
 }
}

sub create_translation_table11{
 my $ref = shift @_;
 my @tt11_aa = unpack 'A' x 64,'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG';
 my $i = 0;
 foreach my $base1 ('T','C','A','G') {
  foreach my $base2 ('T','C','A','G') {
   foreach my $base3 ('T','C','A','G') {
     $ref->{$base1.$base2.$base3} = $tt11_aa[$i++];
   }
  }
 }
}

sub revcomp{
 my $seq = shift @_;
 $seq = join ('', reverse (unpack 'A' x length($seq), $seq));
 $seq =~ tr/actgACTG/tgacTGAC/;
 return $seq
}

sub annotcoding{
 my $from=shift @_;
 my $genename=shift @_;
 my $txStart =$start{$genename};
 my $txEnd = $end{$genename};
 my $strand = $strand{$genename};
 my $genesymbol = $symbol{$genename};
 my $genedesc = $desc{$genename};
 my $ref_allele = shift @_;
 my $allele = shift @_;
 my $cnc = 'C';
 my $sns;
 my $aasnp;
 my $aapos;
 my $checkallele;
 my $nucpos;
 my $codpos;
 my $codon;
 my $altcodon;
 if ($genesymbol =~ /rr[sl]/) {
  $cnc = 'N';
  $sns = '';
  $codon = '\N';
  $altcodon = '\N';
  $codpos = '\N';
  if ($strand  eq '+') {
   $nucpos = $from - $txStart;
  }else{
   $nucpos = $txEnd - $from + 1;
  }
 }else{
  if ($strand eq '+') {
   $nucpos = $from - $txStart;
   $aapos = int(($from - $txStart - 1) / 3) + 1; # from start of gene
   $codpos = ($from - $txStart) % 3; #within the codon 1,2
   $codpos = 3 unless $codpos; # 0->3;
   my $codonStart = int(($from - $txStart - 1)/3)*3 + $txStart + 1; #genomic coordinate of first base in the first codon affected (int truncates anything after the point)
   my $codonEnd = int(($from + (length($allele)-1) - $txStart - 1)/3)*3 + $txStart + 3; #genomic coordinate of last base in last codon affected was=$codonStart + 2;
   $codon = `perl $Bin/get_seq_coord.pl -coord $codonStart-$codonEnd -nodefline $ref_file`; #h37rv.fasta is reference fasta file this and get_seq_coord.pl should be in the $Bin folder
   chomp $codon;
   $codon=~s/\n//g;
   #my $ref_base = substr($codon,$codpos-1,1);
   #my $input_base = substr($ref_allele,0,1);
   #unless (uc($ref_base) eq uc($input_base)) {
    #print STDERR "Reference base at $from ($ref_base) is different than then SNP reference ($input_base)!";
   #}
   $checkallele=$allele;
   $altcodon=$codon;
   substr($altcodon, $codpos-1, length($checkallele))=$checkallele;
  }else{ # '-' strand
   $nucpos = $txEnd - $from + 1;
   $aapos = int(($txEnd - $from) / 3) + 1;
   $codpos = ($txEnd - $from + 1) % 3;
   $codpos = 3 unless $codpos; # 0->3;
   my $codonStart = $txEnd - int(($txEnd - $from)/3)*3 -2;
   my $codonEnd = $txEnd - int(($txEnd - $from -length($allele)+1)/3)*3; #genomic coordinate of last base in last codon affected was=$codonStart + 2;
   $codon = `perl $Bin/get_seq_coord.pl -coord $codonStart-$codonEnd -nodefline $ref_file`; #h37rv.fasta is reference fasta file this and get_seq_coord.pl should be in the $Bin folder
   chomp $codon;
   #$codon = &revcomp($codon); #reverse complement
   my $codposrev;
   $codposrev = 3 if ($codpos==1);
   $codposrev = 1 if ($codpos==3);
   $codposrev = 2 if ($codpos==2);
   #my $ref_base = &revcomp(substr($codon,-$codposrev,1));
   #my $input_base = substr($ref_allele,0,1);
   #unless (uc($ref_base) eq uc($input_base)) {
   # print STDERR "Reference base at $from ($ref_base) is different than then SNP reference ($input_base), allele is ($allele) ref allele is ($ref_allele), codon ($codon) start at pos $codpos of $codonStart ends at $codonEnd!";
   # if (uc(substr($ref_allele,0,1)) eq &revcomp(uc($ref_base))) {
   #  print STDERR "  Assuming the SNP is reported in reverse complement.";
   #  $ref_allele = &revcomp($ref_allele);
   #  $allele = &revcomp($allele);
   # }
   #}
   $altcodon=$codon;
   $checkallele = $allele;
   substr($altcodon, $codposrev-1,length($checkallele))=$checkallele;
   $codon=&revcomp($codon);
   $altcodon=&revcomp($altcodon);
  }
  #print STDERR "altcodon is $altcodon , codon is $codon, $codpos, $allele, $checkallele\n";
  ($sns, my $aaref, my $aavar) = &check_synonymous($codon,$aapos,$altcodon);
  #print STDERR "$sns, $aaref, $aavar\n";
  $aasnp = $aaref.$aapos.$aavar; #general variable
   if ($aasnp =~ /\*/) {
	$sns = 'Z';
   }
 }
 my $name = "SNP_${cnc}${sns}_$from";
 $name .= "_${ref_allele}${nucpos}${allele}"; # rrs and rrl will have the nucleotide positions between ref and alt alleles
 $name .= "_${aasnp}" unless ($genesymbol =~ /rr[sl]/);
 my $temp;
 if ($genesymbol eq '\N') { 
    $temp=$genename;
 } else {
    $temp=$genesymbol;
 }
 $name .= "_".$temp;
 return ($name, $codon, $altcodon, $codpos);
}

sub annotnoncoding{ #snp non coding only
 my $from = shift@_;
 my $genename=shift @_;
 my $txStart =$start{$genename};
 my $txEnd = $end{$genename};
 my $strand = $strand{$genename};
 my $genedesc = $desc{$genename};
 my $before_gene = $geneBefore{$genename};
 my $after_gene = $geneAfter{$genename};
 my $type = $type{$genename};
 my $ref_allele = shift @_;
 my $allele = shift @_;
 my $cnc = 'N';
 my $aapos;
 my $nucpos;
 my $codpos;
 if($type eq 'promoter') {
  $cnc = 'P';
  $nucpos = ".";
  if ($strand eq "+") {
   $nucpos = $txEnd-$from+1;
   if (length($allele)>1) {
    my $u=$nucpos-length($allele)+1;
    $nucpos=".".$nucpos."-".$u;
   }
  } elsif ($strand eq "-") {
   $nucpos = $from-$txStart;
   if (length($allele)>1) {
    my $u=$nucpos+length($allele)-1;
    $nucpos=".".$nucpos."-".$u;
   }
  }
 } else {
  $cnc = 'I';
  $nucpos = $txEnd-$from+1; #by default look at position from next downstream gene
  if (length($allele)>1) {
   my $u=$nucpos-length($allele)+1;
   $nucpos=".".$nucpos."-".$u;
  }
 }
 my $name = "SNP_${cnc}_$from";
 $name .= "_${ref_allele}${nucpos}${allele}"; # rrs and rrl will have the nucleotide positions between ref and alt alleles
 #print STDERR ($ref_allele, $allele, $nucpos, $name); 
 $name.="_".$genedesc if ($genedesc ne '\N');
 $name.="_inter_".${before_gene}."_".$after_gene if ($genedesc eq '\N');
 #print STDERR $name;
 return ($name, '\N', '\N','\N');   
}

sub annotindel {
 my $from = shift@_;
 my $genename=shift @_;
 my $ref_allele = shift @_;
 my $allele = shift @_;
 my $txStart =$start{$genename};
 my $txEnd = $end{$genename};
 my $strand = $strand{$genename};
 my $name; 
 my $kind;
 my $nucpos='';
 if ($genename =~ m/^Rv/i) {
  my $genesymbol = $symbol{$genename};
  my $genedesc = $desc{$genename};
  my $cnc = 'C';
  my $aapos='';
  if ($genesymbol =~ /rr[sl]/) {
   $cnc = 'N';
   if ($strand  eq '+') {
    $nucpos = $from - $txStart - 1;
   }else{
    $nucpos = $txEnd - $from + 1;
   }
  }else{
   if ($strand eq '+') {
    $nucpos = $from - $txStart - 1;
    $aapos = int(($from - $txStart - 1) / 3) + 1; # from start of gene
   }else{ # '-' strand
    $aapos = int(($txEnd - $from) / 3) + 1;
    $nucpos = $txEnd - $from + 1;
   }
  }
  if (length($allele)>length($ref_allele)) {
   $kind="INS";
   if (((length($allele)-length($ref_allele)))%3 ==0) {
    $cnc.='I';
   } else {
    $cnc.='F';
   }
  } else {
   $kind="DEL";
   if (((length($ref_allele)-length($allele)))%3 ==0) {
    $cnc.='D';
   } else {
    $cnc.='F';
   }
  }
  $name = "${kind}_${cnc}_$from";
  if ($kind eq "DEL") {
   $name .= "_d${nucpos}${ref_allele}"; 
  } else {
   $name .= "_i${nucpos}${allele}"; 
  }
  $name .= "_${aapos}" unless ($genesymbol =~ /rr[sl]/);
  my $temp;
  if ($genesymbol eq '\N') {
    $temp=$genename;
  } else {
    $temp=$genesymbol;
  }
  $name .= "_".$temp;
 } else { #indel in non-coding region
  my $genedesc = $desc{$genename};
  my $before_gene = $geneBefore{$genename};
  my $after_gene = $geneAfter{$genename};
  my $type = $type{$genename};
  my $cnc = 'N';
  if($type eq 'promoter') {
   $cnc = 'P';
   if ($strand eq "+") {
    $nucpos = $txEnd-$from+1;
    if (length($allele)>1) {
     my $u=$nucpos-length($allele)+1;
     $nucpos=".".$nucpos."-".$u;
    }
   } elsif ($strand eq "-") {
    $nucpos = $from-$txStart;
    #if (length($allele)>1) {
    # my $u=$nucpos+length($allele)-1;
    # $nucpos=$nucpos."-".$u;
    #}
   }
  } else {
   $cnc = 'I';
   $nucpos = $txEnd-$from+1;
   #if (length($allele)>1) {
   # my $u=$nucpos-length($allele)+1;
   # $nucpos=$nucpos."-".$u;
   #}
  }
  if (length($allele)>length($ref_allele)) {
   $kind="INS";
  } else {
   $kind="DEL";
  }
  $name = "${kind}_${cnc}_$from";
  if ($kind eq "DEL") {
   $name .= "_d${nucpos}${ref_allele}"; 
  } else {
   $name .= "_i${nucpos}${allele}"; 
  }
  $name.="_".$genedesc if ($genedesc ne '\N' && $genedesc ne '');
  $name.="_inter_".${before_gene}."_".$after_gene if ($genedesc eq '\N' || $genedesc eq '');
 }
 return ($name, '\N','\N','\N'); 
}


