#!/usr/bin/perl
use strict;
use warnings;
### script that takes in .vcf file and produces a .var file. Filters and combines the mutation data in .vcf file with data from 2 genome coordinate files (with headers) #####(one for coding ###regions and one for non coding regions) to add functional data. Also requires h37rv.fasta file (there reference genome sequence file) and ###get_seq_coord.pl script to exist ###in the same folder

### example command ./flatAnnotatorVAR.pl test.vcf qual{0-255} hetero{0-1} platypusfilter{PASS|ALL} (output will be test.var) 
use Cwd qw(abs_path);
use FindBin qw($Bin);
use lib abs_path($Bin);
use SeqCoord;

$/="\n";
$,="\t";
$\="\n";
my %tt11;
&create_translation_table11(\%tt11);

my $reference = "h37rv"; #this may need to be changed in the future
my $ref_file = shift @ARGV; # ${h37rv}.fasta
my $Creference = shift @ARGV; # ${h37rv_coding}.txt
my $Nreference = shift @ARGV; # ${h37rv_noncoding}.txt
my $snpfile = shift @ARGV; # ${file}.vcf
#my $refPath= shift @ARGV;
my $qualThresh = (shift@ARGV)||0; # 15
my $heteroThresh = (shift@ARGV)||0; # 0.1
my $platypusFilter = (shift@ARGV)||'PASS'; # PASS
if ($platypusFilter =~ m/pass/i) {
  $platypusFilter = 'PASS'; #NOTE: can include '|badReads|alleleBias' if using stampy/platypus pipeline as calls are conservative
} else {
  #$platypusFilter =~ s/\-/\|/g;
  #$platypusFilter = "qr/(?!\A.*$platypusFilter.*\z)/"; 
  $platypusFilter = ".+?"; #match anything
}
#print STDERR $platypusFilter;
#my $aasnp;
(my $strain) = ($snpfile =~ m/(strain\d+)/);
my $DEBUG=1;

# Load the fasta file into the SeqCoord
load_fasta($ref_file);

#####################################################################################################################################################
#READING ANNOTATION TABLES----READING ANNOTATION TABLES----READING ANNOTATION TABLES----READING ANNOTATION TABLES----READING ANNOTATION TABLES----- #
#####################################################################################################################################################



#print STDERR "reading file $Creference\n";
open REF, "<$Creference" or die;
my $n=1;
my %index;
my %start;
my %end;
my %symbol;
my %strand;
my %desc;
my %region;
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
  $region{$w[$index{"start"}]}=$w[$index{"name"}];
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
  $region{$w[$index{"start"}]}=$w[$index{"name"}];
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

# Only sort once.
my @genomePositions = sort {$a <=> $b} keys %region;
	
#############################################################################################################################################
#READING VCF & Annotating----READING VCF & Annotating----READING VCF & Annotating----READING VCF & Annotating----READING VCF & Annotating---#
#############################################################################################################################################


open IN,"<$snpfile";
my $asnpfile  = $snpfile;

my $OUT = *STDOUT;
print $OUT "reference\tregionid1\tgenesymbol1\tregionid2\tgenesymbol2\tvarname\tvar_len\tdesc1\tdesc2\tqual\tdepth\tbidir\thqr\thqr_ref\tfq\tcodon\taltcodon\tcodpos\tplatypusfilter";
my $counthc;
my $source;
while (<IN>) {
if ($_ =~ m/^source/i) {
  $source = ($_ =~ /source="?(\w+)\s+/i);
} elsif ( $_ =~ m/^#/) {
  next;
}
 chomp;
 my $line=$_;
 my @f=split/\t/;
 my ($from,$ref_allele,$allele,$snpqual,$filter,$info) = ($f[1],@f[3..7]);
 my ($depth, $bidir, $hqr_var, $hqr_ref, $fq, $exclude);

 if ($source and $source =~ m/platypus/i) {
  ($depth, $bidir, $hqr_var, $hqr_ref, $fq, $exclude) = &qualityControl_platypus($line);
 } else {
  ($depth, $bidir, $hqr_var, $hqr_ref, $fq, $exclude) = &qualityControl_pilon($line);
 }
 if ($exclude ==0 ) { 
   my ($genename,$lgenename, $codon,$codpos,$altcodon,$name);#$nucpos
   my ($type,$before_gene,$after_gene);
   if (length($ref_allele)==length($allele)) {#SNP
    my $march=$from;
    #while ($march<=$from+length($allele)-1) {
    foreach my $genomePos (@genomePositions) { 
      if ($genomePos <$march && $end{$region{$genomePos}}>=$march) { #for zero based coordinates
        $genename=$region{$genomePos};
      }
      if ($genomePos <($march+length($allele)-1) && $end{$region{$genomePos}}>=($march+length($allele)-1)) { #for zero based coordinates
        $lgenename=$region{$genomePos};
      }
     }
     if ($lgenename ne $genename) {
      #print STDERR "$genename, $lgenename, $from, $march\n";
      my $ccoord=($end{$genename}-$march)|1;
      #print STDERR "$ccoord\n$end{$genename}-$march\n";
      my $ref_bases1= substr($ref_allele, 0, $ccoord);
      my $ref_bases2= substr($ref_allele, -($march + length($ref_allele)- $start{$lgenename}-1));
      #print STDERR "$ref_allele, $ref_bases1, $ref_bases2\n";
      my $bases1= substr($allele, 0, $ccoord);
      my $bases2= substr($allele, -($march + length($ref_allele)- $start{$lgenename}-1));  #assume here that gene only spans two regions
      #print STDERR "$allele, $bases1, $bases2\n";
      my $sym1=''; my $sym2='';
      if ($genename =~ /^Rv/i) { #1st coding
       (my $name1, my $codon1, my $altcodon1, my $codpos1, $march)=&assign_allele($genename, $ref_bases1, $bases1, $march, $from);
       #print STDERR $name1, $march;
        if ($symbol{$genename} ne '\N') {
                $sym1=$symbol{$genename};
        }
       if ($lgenename =~ /^Rv/i) {#2nd coding
        (my $name2, my $codon2, my $altcodon2, my $codpos2, $march)=&assign_allele($lgenename, $ref_bases2, $bases2, $start{$lgenename}+1, $start{$lgenename}+1);
	if ($symbol{$lgenename} ne '\N') {
		$sym2=$symbol{$lgenename};
	}
        print $OUT $reference,$genename, $sym1, $lgenename,$sym2,$name1."&".$name2,length($ref_allele),$desc{$genename},$desc{$lgenename},$snpqual,$depth,$bidir,$hqr_var,$hqr_ref,$fq,$codon1."&".$codon2,$altcodon1."&".$altcodon2,$codpos1."&".$codpos2,$filter;
        $march++;
       } else { #2nd noncoding
        (my $name2, my $codon2, my $altcodon2, my $codpos2)=&annotnoncoding($start{$lgenename}+1, $lgenename, $ref_bases2, $bases2);
        #print STDERR $name2, $march;
        print $OUT $reference,$genename, $sym1,$lgenename,$symbol{$lgenename},$name1."&".$name2,length($ref_allele),$desc{$genename},$desc{$lgenename},$snpqual,$depth,$bidir,$hqr_var,$hqr_ref,$fq,$codon1,$altcodon1,$codpos1,$filter;
        $march=$from+length($allele)-1;
       }
       $march++;
      } else { #1st noncoding
       (my $name1, my $codon1, my $altcodon1, my $codpos1)=&annotnoncoding($march, $genename, $ref_bases1, $bases1);
       if ($lgenename =~ /^Rv/i) {#2nd coding
        (my $name2, my $codon2, my $altcodon2, my $codpos2, $march)=&assign_allele($lgenename, $ref_bases2, $bases2, $start{$lgenename}+1, $start{$lgenename}+1);
        #my $test=$start{$lgenename}+1;  
        #print STDERR $allele, $bases1, $bases2, $lgenename, $test, $from;
        print $OUT $reference,$genename, $symbol{$genename},$lgenename,$symbol{$lgenename},$name1."&".$name2,length($ref_allele),$desc{$genename},$desc{$lgenename},$snpqual,$depth,$bidir,$hqr_var,$hqr_ref,$fq,$codon2, $altcodon2,$codpos2,$filter;
        $march++;
       } else { #2nd noncoding
        (my $name2, my $codon2, my $altcodon2, my $codpos2)=&annotnoncoding($start{$lgenename}+1, $lgenename, $ref_bases2, $bases2);
        print $OUT $reference,$genename, $symbol{$genename},$lgenename,$symbol{$lgenename},$name1."&".$name2,length($ref_allele),$desc{$genename},$desc{$lgenename},$snpqual,$depth,$bidir,$hqr_var,$hqr_ref,$fq,$codon1, $altcodon1, $codpos1, $filter;
        $march=$from+length($allele)-1;
       }
      }
     } else { #snp that only spans one region
      if ($genename =~ /^Rv/i) { #coding
       ($name, $codon, $altcodon, $codpos, $march)=&assign_allele($genename, $ref_allele, $allele, $march, $from);
      } else { ## 1 noncoding variant
       my $ref_base=$ref_allele; #substr($ref_allele,$march-$from,1);
       my $base=$allele; #substr($allele,$march-$from,1);
       ($name, $codon, $altcodon, $codpos)=&annotnoncoding($march, $genename,$ref_base, $base);
       $march=$march+length($allele)-1;
      }
      print $OUT $reference,$genename,$symbol{$genename}."\t\\N\t\\N",$name,length($ref_allele),$desc{$genename}."\t\\N",$snpqual,$depth,$bidir,$hqr_var,$hqr_ref,$fq,$codon,$altcodon,$codpos,$filter;    
      $march++;
     }
    #}
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
     my $len_var;
     #print STDERR "warning combo snp and indel at pos $from\n";
     my $march=$from;
     my ($name1, $codon1, $altcodon1, $codpos1);
     my ($name2, $codon2, $altcodon2, $codpos2);
     if (length($allele)>length($ref_allele)) { #at least partly an insertion
      $len_var=length($allele);	
      foreach my $genomePos (@genomePositions) {
       if ($genomePos <$march && $end{$region{$genomePos}}>=$march) { #for zero based coordinates
        $genename=$region{$genomePos};
       }
      }
      my $ref_base=substr($ref_allele,$march-$from,length($ref_allele));
      my $base=substr($allele,$march-$from,length($ref_allele));
      if ($genename =~ /^Rv/i) { # 1 gene
       ($name1, $codon1, $altcodon1, $codpos1, $march)=&assign_allele($genename, $ref_base, $base, $march, $from);
      } else { ## 1 noncoding variant
       ($name1, $codon1, $altcodon1, $codpos1)=&annotnoncoding($march, $genename,$ref_base, $base);
       $march=$march+length($ref_allele)-1;
      }
      $march++;
      $allele=substr($allele, length($ref_allele)); 
      $ref_allele='';
     } elsif (length($ref_allele)>length($allele)) { #at least partly a deletion
      $len_var=length($ref_allele);
      foreach my $genomePos (@genomePositions) {
       if ($genomePos <$march && $end{$region{$genomePos}}>=$march) { #for zero based coordinates
        $genename=$region{$genomePos};
       }
      }
      my $ref_base=substr($ref_allele,$march-$from,length($allele));
      my $base=substr($allele,$march-$from,length($allele));
      if ($genename =~ /^Rv/i) { # 1 gene
       ($name1, $codon1, $altcodon1, $codpos1, $march)=&assign_allele($genename, $ref_base, $base, $march, $from);
      } else { ## 1 noncoding variant
       ($name1, $codon1, $altcodon1, $codpos1)=&annotnoncoding($march, $genename,$ref_base, $base);
       $march=$march+length($allele)-1;
      }
      $march++;
      $ref_allele=substr($ref_allele, length($allele));
      $allele='';
     } #completes SNP porition of INS or DEL $march variable marks where we are
     my $lgenename;
     foreach my $genomePos (@genomePositions) {
      if ($genomePos < $march && $end{$region{$genomePos}}>= $march) { #for zero based coordinates
       $lgenename=$region{$genomePos};
      }
     }
     ($name2,$codon2,$altcodon2,$codpos2) = &annotindel($march,$lgenename,$ref_allele, $allele);
     if ($lgenename ne $genename) {
       print $OUT $reference,$genename, $symbol{$genename},$lgenename,$symbol{$lgenename},$name1.'&'.$name2,$len_var,$desc{$genename}, $desc{$lgenename},$snpqual,$depth,$bidir,$hqr_var,$hqr_ref,$fq,$codon1,$altcodon1,$codpos1,$filter;
     } else {
       print $OUT $reference,$genename,$symbol{$genename}."\t\\N\t\\N",$name1.'&'.$name2,$desc{$genename},$snpqual,$depth,$bidir,$hqr_var,$hqr_ref,$fq,$codon1,$altcodon1,$codpos1,$filter;
     }
     $march++;
    } else { #completes combo SNP/indel, below this is pure INDEL
     foreach my $genomePos (@genomePositions) {
      if ($genomePos <$from && $end{$region{$genomePos}}>=$from) { #for zero based coordinates
        $genename=$region{$genomePos};
      }
     }
     my $len_var=length($allele)||length($ref_allele);
     ($name,$codon,$altcodon,$codpos) = &annotindel($from,$genename,$ref_allele, $allele);
     print $OUT $reference,$genename,$symbol{$genename}."\t\\N\t\\N",$name,$len_var,$desc{$genename}."\t\\N",$snpqual,$depth,$bidir,$hqr_var,$hqr_ref,$fq,$codon,$altcodon,$codpos,$filter;
    }
   } #close indel
  } #close high Quality
} #vcf file complete

close SNP;
close IN;

#####################################################################################################################################################
#SUBROUTINES----SUBROUTINES----SUBROUTINES----SUBROUTINES----SUBROUTINES----SUBROUTINES----SUBROUTINES----SUBROUTINES----SUBROUTINES----SUBROUTINES #
#####################################################################################################################################################

#-----------------------------------------------------------------------------------------------------------------------
#  Loop through a coding genetic substitution of one or more bases length and split into codon level (3bp) changes
#-----------------------------------------------------------------------------------------------------------------------

sub assign_allele{
 my $genename= shift @_;
 my $ref_allele=shift @_;
 my $allele=shift @_;
 my $march=shift @_;
 my $from =shift @_;
 my $cp; my $name; my $cnc; my $sns; my $aaref=''; my $aapos=''; my $aavar=''; my $codon; my $altcodon; my $codpos; my $nucpos;
 my $fsns; my $faaref=''; my $faavar=''; my $faapos=''; my $fcodon=''; my $faltcodon=''; my $fnucpos;
 my $ref_base;
 while ($march <= $from+length($allele)-1) {
  if ($strand{$genename} eq "+") {
   $cp = ($march - $start{$genename}) % 3; #within the codon 1,2; for allele changes >1 this is the starting codon position
   $cp = 3 unless $cp; # 0->3;
   #print STDERR "$cp", length($allele);
   if ($cp == 3 || length($allele)== 1) {
    $ref_base=substr($ref_allele,$march-$from,1);
    my $base=substr($allele,$march-$from,1);
    ($cnc, $sns, $nucpos, $aaref, $aapos, $aavar, $codon, $altcodon, $codpos)=&annotcoding($march, $genename, $ref_base, $base);
   } elsif ($cp != 3 && length($allele)>1) {
    if (($cp == 1 && length($allele)<=3) || ($cp ==2 && length($allele)<=2)) {
     #print STDERR "HERE!!";
     $ref_base=substr($ref_allele,$march-$from, length($allele));
     my $base=substr($allele,$march-$from,length($allele));
     ($cnc, $sns, $nucpos, $aaref, $aapos, $aavar, $codon, $altcodon, $codpos)=&annotcoding($march, $genename,$ref_base, $base);
     $march=$march+length($allele)-1;
    } elsif ($cp == 1 && length($allele)>3) {
     $ref_base=substr($ref_allele,$march-$from, 3);
     my $base=substr($allele,$march-$from, 3);
     ($cnc, $sns, $nucpos, $aaref, $aapos, $aavar, $codon, $altcodon, $codpos)=&annotcoding($march, $genename,$ref_base, $base);
     $march=$march+2;
    } elsif ($cp == 2 && length($allele)>2) {
     $ref_base=substr($ref_allele,$march-$from, 2);
     my $base=substr($allele,$march-$from,2);
     ($cnc, $sns, $nucpos, $aaref, $aapos, $aavar, $codon, $altcodon, $codpos)=&annotcoding($march, $genename,$ref_base, $base);
     $march++;
    }
   }
  } elsif ($strand{$genename} eq "-") {
   $cp = ($end{$genename} - $march + 1) % 3;
   $cp = 3 unless $cp; # 0->3;
   if ($cp == 1 || length($allele)== 1) {  #codon count for negative strand is 1,3,2,1,3,2
    $ref_base=substr($ref_allele,$march-$from,1);
    my $base=substr($allele,$march-$from,1);
    ($cnc, $sns, $nucpos, $aaref, $aapos, $aavar, $codon, $altcodon, $codpos)=&annotcoding($march, $genename, $ref_base, $base);
   } else { #$cp!=1 and length($allele)>1
    if (($cp == 3 && length($allele)<=3) || ($cp ==2 && length($allele)<=2)) {
     #print STDERR "HERE!!";
     $ref_base=substr($ref_allele,$march-$from, length($allele));
     my $base=substr($allele,$march-$from,length($allele));
     ($cnc, $sns, $nucpos, $aaref, $aapos, $aavar, $codon, $altcodon, $codpos)=&annotcoding($march, $genename,$ref_base, $base);
     $march=$march+length($allele)-1;
    } elsif ($cp == 3 && length($allele)>3) {
     $ref_base=substr($ref_allele,$march-$from, 3);
     my $base=substr($allele,$march-$from, 3);
     ($cnc, $sns, $nucpos, $aaref, $aapos, $aavar, $codon, $altcodon, $codpos)=&annotcoding($march, $genename,$ref_base, $base);
     $march=$march+2;
    } elsif ($cp == 2 && length($allele)>2) {
     $ref_base=substr($ref_allele,$march-$from, 2);
     my $base=substr($allele,$march-$from,2);
     ($cnc, $sns, $nucpos, $aaref, $aapos, $aavar, $codon, $altcodon, $codpos)=&annotcoding($march, $genename,$ref_base, $base);
     $march++;
    }
   } 
  }
  if (length($allele) > 1) {
   $faaref .= $aaref;
   $faavar .= $aavar;
   $fcodon .= $codon;
   $faltcodon .= $altcodon;
   unless (length($faapos) > 0) {
    $faapos=$aapos;
    $fsns=$sns;
    $fnucpos=$nucpos;
   }
   if (length($fsns) >0 && $fsns ne $sns && $fsns eq 'S') {
    $fsns=$sns;
   }  
  } else {
   $faaref = $aaref;
   $faavar = $aavar;
   $fcodon = $codon;
   $faltcodon = $altcodon;
   $fsns = $sns;
  }
  $march++;
 }
 if (length($allele) <=3) {
	$name = "SNP_${cnc}${fsns}_${from}";
 } else {
	$name = "LSP_${cnc}${fsns}_${from}";
 }
 my $genesymbol=$symbol{$genename};
 if (length($allele) >1) {
  if ($nucpos >$fnucpos) {
   $nucpos=$nucpos+length($ref_base)-1;
  } else {
   $nucpos=$nucpos-length($ref_base)+1;
  }
  $name .= "_${ref_allele}${fnucpos}-${nucpos}${allele}"; # rrs and rrl will have the nucleotide positions between ref and alt alleles    	
 } else { #length==1
  $name .= "_${ref_allele}${nucpos}${allele}"; # rrs and rrl will have the nucleotide positions between ref and alt alleles
 }
 if ($genesymbol =~ /rr[sl]/) {
  $fcodon='\N';
  $faltcodon='\N';
 } else {
  if (length($allele) >1 && $faapos != $aapos) {
   $name .= "_${faaref}${faapos}-${aapos}${faavar}" 
  } else {
   $name .= "_${faaref}${aapos}${faavar}";
  }
 }
 my $temp;
 if ($genesymbol eq '\N') {
    $temp=$genename;
 } else {
    $temp=$genesymbol;
 }
 $name .= "_".$temp;
 return ($name, $fcodon, $faltcodon, $codpos, $march);
}

#---------------------------------------------------------------------------------------
#  Check if a single nucleotide subs is synonymous
#---------------------------------------------------------------------------------------

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


#---------------------------------------------------------------------------------------
#  Annotate a coding SNP
#---------------------------------------------------------------------------------------

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
 my $aapos;
 my $aaref; my $aavar;
 my $checkallele;
 my $nucpos;
 my $codpos;
 my $codon;
 my $altcodon;
 if ($genesymbol =~ /rr[sl]/) {
  $cnc = 'N';
  $sns = '';
  $codon='\N';
  $altcodon='\N';
  $codpos='\N';
  $aapos='';
  $aaref='';
  $aavar='';
  if ($strand  eq '+') {
   $nucpos = $from - $txStart;
  }else{
   $nucpos = $txEnd - $from + 1;
  }
 }else{
  if ($strand eq '+') {
   $nucpos = $from - $txStart;
   #print STDERR "$ref_allele, $allele, $strand, $genename, $from, $txStart, $nucpos\n";
   $aapos = int(($from - $txStart - 1) / 3) + 1; # from start of gene
   $codpos = ($from - $txStart) % 3; #within the codon 1,2
   $codpos = 3 unless $codpos; # 0->3;
   my $codonStart = int(($from - $txStart - 1)/3)*3 + $txStart + 1; #genomic coordinate of first base in the first codon affected (int truncates anything after the point)
   my $codonEnd = int(($from + (length($allele)-1) - $txStart - 1)/3)*3 + $txStart + 3; #genomic coordinate of last base in last codon affected was=$codonStart + 2;
   $codon = get_coord([$codonStart, $codonEnd]);
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
   $codon = get_coord([$codonStart, $codonEnd]);
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
  ($sns, $aaref, $aavar) = &check_synonymous($codon,$aapos,$altcodon);
  if (($aavar =~ /\*/ && $aaref !~ /\*/ ) || ($aaref =~/\*/ && $aavar !~ /\*/)) {
	$sns = 'Z';
  }
 }
 return ($cnc, $sns, $nucpos, $aaref, $aapos, $aavar, $codon, $altcodon, $codpos);
}


#---------------------------------------------------------------------------------------
#  Annotate a non coding SNP
#---------------------------------------------------------------------------------------

sub annotnoncoding{ 
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
    $nucpos=$nucpos."-".$u;
   }
  } elsif ($strand eq "-") {
   $nucpos = $from-$txStart;
   if (length($allele)>1) {
    my $u=$nucpos+length($allele)-1;
    $nucpos=$nucpos."-".$u;
   }
  }
 } else {
  $cnc = 'I';
  $nucpos = $txEnd-$from+1; #by default look at position from next downstream gene
  if (length($allele)>1) {
   my $u=$nucpos-length($allele)+1;
   $nucpos=$nucpos."-".$u;
  }
 }
 my $name;
 if (length($allele)>3) {
 	$name= "LSP_${cnc}_$from";
 } else {
 	$name = "SNP_${cnc}_$from";
 }
 $name .= "_${ref_allele}${nucpos}${allele}"; # rrs and rrl will have the nucleotide positions between ref and alt alleles
 #print STDERR ($ref_allele, $allele, $nucpos, $name); 
 $name.="_".$genedesc if ($genedesc ne '\N');
 $name.="_inter-".${before_gene}."-".$after_gene if ($genedesc eq '\N');
 #print STDERR $name;
 return ($name, '\N', '\N','\N');   
}


#---------------------------------------------------------------------------------------
#  Annotate an coding or noncoding INDEL variant
#---------------------------------------------------------------------------------------

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
 my $codpos=0; #where the indel starts
 my $sns;
 my $aavar; #aa sequence inserted or deleted -(0-2) bases at the begining to preserve translation frame
 my $codon; #H37Rv codon that contains the start position of the indel
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
    $codpos = ($from - $txStart) % 3; #within the codon 1,2
    $codpos = 3 unless $codpos; # 0->3;
    my $codonStart = int(($from - $txStart - 1)/3)*3 + $txStart + 1; #genomic coordinate of first base in the first codon affected (int truncates anything after the point)
    my $codonEnd = $codonStart + 2;
    $codon = get_coord([$codonStart, $codonEnd]);
   }else{ # '-' strand
    $aapos = int(($txEnd - $from + 1 ) / 3) + 1;
    $nucpos = $txEnd - $from + 1 ;
    $codpos = ($txEnd - $from + 2) % 3;
    $codpos = 3 unless $codpos; # 0->3;
    my $codonStart = $txEnd - int(($txEnd - $from +1)/3)*3 -2;
    my $codonEnd = $codonStart + 2;
    $codon = get_coord([$codonStart, $codonEnd]);
    $codon=&revcomp($codon);
    }
  }
  if (length($allele)>length($ref_allele)) {
   my $alleleTotranslate;
   $kind="INS";
   if (((length($allele)-length($ref_allele)))%3 ==0) {
    $cnc.='I';
   } else {
    $cnc.='F';
   }
   if ($strand eq "-") {
    $allele=&revcomp($allele);
   }
   if ($codpos && $codpos==1) {
	if (length($allele)<3) {
		$alleleTotranslate=$allele.substr($codon,0,3-length($allele));
	} else {
		my $remain=(length($allele))%3;
		$alleleTotranslate=substr($allele, 0, length($allele)-$remain);
	}
   } elsif ($codpos==2) {
   	$alleleTotranslate=substr($codon,0,1).$allele;
	if (length($allele)==1) {
		$alleleTotranslate=$alleleTotranslate.substr($codon,1,1);
	} else {
		my $remain=(length($alleleTotranslate))%3;
	        my $alleleTotranslate2=substr($alleleTotranslate, 0, length($alleleTotranslate)-$remain);
       		$alleleTotranslate=$alleleTotranslate2;
        }
   } else {#codpos is 3
        $alleleTotranslate=substr($codon,0,2).$allele;
        if (length($allele)>1) {
	        my $remain=(length($alleleTotranslate))%3;
       		my $alleleTotranslate2=substr($alleleTotranslate, 0, length($alleleTotranslate)-$remain);
		$alleleTotranslate=$alleleTotranslate2;
        }	
   }
   #print STDERR "strand is $strand, codpos is $codpos, codon is $codon, will translate $alleleTotranslate original is $allele\n";
   ($sns, $aavar, $aavar) = &check_synonymous($alleleTotranslate,$aapos,$alleleTotranslate);
  } else {
   $kind="DEL";
   $aavar="";
   if (((length($ref_allele)-length($allele)))%3 ==0) {
    $cnc.='D';
   } else {
    $cnc.='F';
   }
  }
  if ($aavar ne "" && $aavar =~ /\*/) {
     $cnc =~	s/I|F/Z/;
  }
  $name = "${kind}_${cnc}_$from";
  if ($kind eq "DEL") {
   $name .= "_d${nucpos}${ref_allele}"; 
  } else {
   $name .= "_i${nucpos}${allele}"; 
  }
  $name .= "_${aapos}${aavar}" unless ($genesymbol =~ /rr[sl]/);
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
  $name.="_inter-".${before_gene}."-".$after_gene if ($genedesc eq '\N' || $genedesc eq '');
 }
 return ($name, '\N','\N','\N'); 
}


#---------------------------------------------------------------------------------------
#  Check the qualtiy of a variant PILON format
#---------------------------------------------------------------------------------------

sub qualityControl_pilon
{
  my $line = shift @_;  ##VCF full line
  my $ex=0;
  my @elements = split /\t/, $line;
  my ($from,$ref_allele,$allele,$snpqual,$filter,$info)=($elements[1],@elements[3..7]);
  $ex=1 if $info =~ /IMPRECISE/i;
  (my $depth) = ($info =~ /DP=(\d+)/);
  $depth=0 unless defined $depth;
  #(my $A) = ($info =~ /BC=(\d+),/);
  #(my $C) = ($info =~ /BC=\d+,(\d+),/);
  #(my $G) = ($info =~ /BC=\d+,\d+,(\d+),/);
  #(my $T)  = ($info =~ /BC=\d+,\d+,\d+,(\d+)/);
  $ex=1 if $ref_allele =~ /N/; #ambiguous reference
  $ex=1 if $allele =~ /N/; #ambiguous/imprecise change
  $ex=1 if $allele =~ m/,|</; #heterogenous allele or <dup>
  (my $hqr)= ($info =~ /AF=((\d+(\.\d*)?)|(\.\d+))$/)||0.7666666;    ###MARTIN PLEASE CHECK THIS REGEX
  if ($hqr eq ".") { #for some indels this may happen
   $hqr=0.7666666; #to patch this error
  }
  if ($snpqual eq ".") {
    $snpqual=21;
  }
  if ($filter =~ m/;/ || $snpqual <=$qualThresh || $filter !~ m/$platypusFilter/  || $hqr <=$heteroThresh ) {
      $ex=1;
  }
  my $hqr_var = $hqr * $depth;
  my $hqr_ref = (1-$hqr) * $depth;
  return $depth, '-', $hqr_var, $hqr_ref, -1, $ex;
}

#---------------------------------------------------------------------------------------
#  Check the qualtiy of a variant Platypus format
#---------------------------------------------------------------------------------------

sub qualityControl_platypus
{
  my $line = shift @_;  ##VCF full line
  my $ex=0;
  my @elements = split /\t/, $line;
  my ($from,$ref_allele,$allele,$snpqual,$filter,$info)=($elements[1],@elements[3..7]);
  (my $depth) = ($info =~ /TC=(\d+)/);
  (my $tcf)   = ($info =~ /TCF=(\d+)/);
  (my $tcr)   = ($info =~ /TCR=(\d+)/);
  (my $nf)    = ($info =~ /NF=(\d+)/);
  (my $nr)    = ($info =~ /NR=(\d+)/);
  my ($dpr1,$dpr2,$dp1,$dp2) = ($tcf-$nf, $tcr-$nr, $nf, $nr);
  my ($fq) = ($info =~ /RMSMQ=([-.\d]+)/)||'\N';
  my $bidir = '';
  $bidir = ($dp1 && $dp2)?'Y':'N';
  $ex=1 if $ref_allele =~ /N/; #ambiguous reference
  $ex=1 if $allele =~ /N/; #ambiguous/imprecise change
  $ex=1 if $allele =~ m/,/; 
  my $hqr_var = $dp1 + $dp2;
  my $hqr_ref = $dpr1 + $dpr2;
  if ($hqr_var == 0 && $hqr_var ==0) { # to protect against division by zero error
   $hqr_ref=1;
  }
  if ($snpqual eq ".") {
    $snpqual=21;
  }
  if ($filter =~ m/;/ || $snpqual <=$qualThresh || $filter !~ m/$platypusFilter/ || $hqr_var/($hqr_var + $hqr_ref) < $heteroThresh) { 
    $ex=1;
  }
  return $depth, $bidir, $hqr_var, $hqr_ref, $fq, $ex;
}
