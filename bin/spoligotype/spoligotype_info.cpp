// $Id: spoligotype.cpp,v 1.9 2013/05/18 23:52:37 ioerger Exp $

// Copyright 2014, Thomas R. Ioerger

// compilation command: g++ -std=c++0x spoligotype_fasta.cpp -o spoligotype_fasta

#include <unordered_map>
#include <fstream>
#include <vector>
#include <string>
#include <string.h>

using namespace std;

extern "C" {

typedef vector<char*> SEQS;

#define BUFLEN 1000
#define EOL '\n'
int W;
#define L 16 // reads with <=2 mismatches have at least one contiguous stretch of 16 matches

void usage()
{
  printf("example: spoligotype <reads> [-maxmm=N] [-mincnt=N] [-maxreads=N]\n");
  exit(0);
}

SEQS* read_eland_query_sequences(char* fname,vector<char*>& headers,int max)
{
  SEQS* seqs=new SEQS();
  char buf[BUFLEN];
  ifstream fil(fname);
  while (fil.getline(buf,BUFLEN))
  {
    if (buf[0]=='#') continue;
    headers.push_back(strdup(buf+1));
    fil.getline(buf,BUFLEN);
    seqs->push_back(strdup(buf));
    if (max>0 && seqs->size()>max) break;
  }
  return seqs;
}

void reverse_complement(const char* src,char* dest)
{
  int n=strlen(src);
  const char* nucs="AGCT";
  const char* comp="TCGA";
  for (int i=0 ; i<n ; i++)
  {
    //dest[n-i-1] = 'N';                                                                                       
    dest[n-i-1] = src[i];
    for (int j=0 ; j<4 ; j++)
      if (src[i]==nucs[j]) dest[n-i-1] = comp[j];
  }
  dest[n] = '\0';
}

int count_mismatches(const char* G,int i,char* H,int j,int w,int m)
{
  int cnt=0;
  for (int k=0 ; k<w ; k++)
  {
    if (G[i+k]!=H[j+k]) cnt++;
    if (cnt>=m) return cnt;
  }
  return cnt;
}

// table of spoligotype patterns:
// http://www.biomedcentral.com/content/supplementary/1471-2148-6-107-S4.pdf
// http://www.biomedcentral.com/1471-2148/6/107
// genetic rearrangements in PGRS 17 and 18

// SpolDB4:
//  http://www.pasteur-guadeloupe.fr/tb/bd_myco.html
//  http://www.pasteur-guadeloupe.fr/tb/spoldb4/spoldb4.pdf

int main(int argc,char** argv)
{
  SEQS* seqs;
  char buf[BUFLEN];
  char seq[100];
  char* pseq;
  int i,j,k;

  int maxmm=0;
  int MIN_CNT=2;
  int max_reads=1500000;

  if (argc==1) usage();

  for (int i=2 ; i<argc ; i++)
  {
    if (strncmp(argv[i],"-maxmm=",7)==0) { sscanf(argv[i]+7,"%d",&maxmm); continue; }
    if (strncmp(argv[i],"-mincnt=",8)==0) { sscanf(argv[i]+8,"%d",&MIN_CNT); continue; }
    if (strncmp(argv[i],"-maxreads=",10)==0) { sscanf(argv[i]+10,"%d",&max_reads); continue; }
    usage();
  }

  vector<char*> headers;
  seqs = read_eland_query_sequences(argv[1],headers,max_reads);
  W = strlen(seqs->at(0));
  printf("# reads: %lu\n",seqs->size());
  printf("# read length: %d\n",W);
  if (argc>2) sscanf(argv[2],"%d",&maxmm);
  printf("# mismatches: %d\n",maxmm);
  printf("# min_count: %d\n",MIN_CNT);
  printf("# max_reads: %d\n",max_reads);
  fflush(stdout);

  printf("# note: modified oligo for spacer 2 to match H37Rv\n");
  printf("# note: modified oligo for spacer 1 (5/18/13)\n");

  // concat [1:35] onto end of genome and rev to simulate circular chromosome
  // filter out low quality seqs?

  // DR-1 Rep097 starts at 3119185 in H37Rv
  // these are in reverse order on complement strand
  // spacer 43 starts at 3119455 in H37Rv

  // changed spacer 1, 5/18/13

  const char* spacers[] = {
//"ATAGAGGGTCGCCGGTTCTGGATCA", // 1
"ATAGAGGGTCGCCGGCTCTGGATCA", // 1 changed, 5/18/13
//"CCTCATAATTGGGCGACAGCTTTTG", // 2 rc = CAAAAGCTGTCGCCCAATTATGAGG
  "CCTCATGCTTGGGCGACAGCTTTTG", // change AA->GC to match H37Rv in region 3123576 -400 0
"CCGTGCTTCCAGTGATCGCCTTCTA", // 3
"ACGTCATACGCCGACCAATCATCAG", // 4
"TTTTCTGACCACTTGTGCGGGATTA", // 5
"CGTCGTCATTTCCGGCTTCAATTTC", // 6
"GAGGAGAGCGAGTACTCGGGGCTGC", // 7
"CGTGAAACCGCCCCCAGCCTCGCCG", // 8
"ACTCGGAATCCCATGTGCTGACAGC", // 9
"TCGACACCCGCTCTAGTTGACTTCC", // 10
"GTGAGCAACGGCGGCGGCAACCTGG", // 11
"ATATCTGCTGCCCGCCCGGGGAGAT", // 12
"GACCATCATTGCCATTCCCTCTCCC", // 13
"GGTGTGATGCGGATGGTCGGCTCGG", // 14
"CTTGAATAACGCGCAGTGAATTTCG", // 15
"CGAGTTCCCGTCAGCGTCGTAAATC", // 16
"GCGCCGGCCCGCGCGGATGACTCCG", // 17
"CATGGACCCGGGCGAGCTGCAGATG", // 18
"TAACTGGCTTGGCGCTGATCCTGGT", // 19
"TTGACCTCGCCAGGAGAGAAGATCA", // 20
"TCGATGTCGATGTCCCAATCGTCGA", // 21
"ACCGCAGACGGCACGATTGAGACAA", // 22
"AGCATCGCTGATGCGGTCCAGCTCG", // 23
"CCGCCTGCTGGGTGAGACGTGCTCG", // 24
"GATCAGCGACCACCGCACCCTGTCA", // 25
"CTTCAGCACCACCATCATCCGGCGC", // 26
"GGATTCGTGATCTCTTCCCGCGGAT", // 27
"TGCCCCGGCGTTTAGCGATCACAAC", // 28
"AAATACAGGCTCCACGACACGACCA", // 29
"GGTTGCCCCGCGCCCTTTTCCAGCC", // 30
"TCAGACAGGTTCGCGTCGATCAAGT", // 31
"GACCAAATAGGTATCGGCGTGTTCA", // 32
"GACATGACGGCGGTGCCGCACTTGA", // 33
"AAGTCACCTCGCCCACACCGTCGAA", // 34
"TCCGTACGCTCGAAACGCTTCCAAC", // 35
"CGAAATCCAGCACCACATCCGCAGC", // 36
"CGCGAACTCGTCCACAGTCCCCCTT", // 37
"CGTGGATGGCGGATGCGTTGTGCGC", // 38
"GACGATGGCCAGTAAATCGGCGTGG", // 39
"CGCCATCTGTGCCTCATACAGGTCC", // 40
"GGAGCTTTCCGGCTTCTATCAGGTA", // 41
"ATGGTGGGACATGGACGAGCGCGAC", // 42
"CGCAGAATCGCACCGGGTGCGGGAG", // 43
NULL};

  int present[50];
  k = 0;  
  while (spacers[k]!=NULL)
  {
    const char* spacer=spacers[k];
    int l=strlen(spacer);
    char sp_rc[BUFLEN];
    reverse_complement(spacer,sp_rc);
    int cnt=0;

    if (maxmm==0) // use strstr() for fast string matching
    {
      for (i=0 ; i<seqs->size() ; i++)
        if (strstr(seqs->at(i),spacer)!=NULL || strstr(seqs->at(i),sp_rc)!=NULL) cnt++;
      printf("%d %s %d\n",k+1,spacer,cnt);
      if (cnt<MIN_CNT) present[k] = 0;
      else present[k] = 1;
      fflush(stdout);
      k++;
    }
    else // use count_mismatches, slower
    {
     for (i=0 ; i<seqs->size() ; i++)
     {
      //if (cnt>=MIN_CNT) break;
      int found=0;
      for (j=0 ; j<W-l ; j++)
      {
        int mm=count_mismatches(spacer,0,seqs->at(i),j,l,maxmm+1);
        if (mm<=maxmm) found = 1;
        mm = count_mismatches(sp_rc,0,seqs->at(i),j,l,maxmm+1);
        if (mm<=maxmm) found = 1;
      }
      cnt += found;
      //if (cnt>=MIN_CNT) break;
     }
     printf("%d %s %d\n",k+1,spacer,cnt);
     if (cnt==0) present[k] = 0;
     else present[k] = 1;
     fflush(stdout);
     k++;
    }
  }

  for (i=0 ; i<43 ; i++) printf("%c",present[i]+'0');
  printf("%c\n",present[42]+'0');
  
  printf("! Octal:");
  printf("[%s]", argv[1]); // pass along file name for parsing
  for (i=0 ; i<14 ; i++)
  { 
    int x=present[3*i]*4+present[3*i+1]*2+present[3*i+2];
    printf("%c",x+'0');
  }
  printf("%c\n",present[42]+'0');
}

}
