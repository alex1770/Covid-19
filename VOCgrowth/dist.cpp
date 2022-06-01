// Distance between two aligned fasta files

#include <stdio.h>
#include <getopt.h>
#include <assert.h>
#include <error.h>
#include <string>
#include <fstream>
#include <iostream>
using std::string;
typedef unsigned char UC;

int main(int ac,char **av){
  bool showall=false;
  int ignorestart=0,ignoreend=0;
  std::ios_base::sync_with_stdio(false);
  while(1)switch(getopt(ac,av,"ae:s:")){
    case 'a': showall=true;break;
    case 's': ignorestart=atoi(optarg);break;
    case 'e': ignoreend=atoi(optarg);break;
    case -1: goto ew0;
    default: goto err0;
  }
 ew0:
  if(optind+2!=ac){
  err0:
    fprintf(stderr,"Usage: dist [options] file1.fasta file2.fasta\n");
    fprintf(stderr,"       -a         Show result for each genome\n");
    fprintf(stderr,"       -s<int>    Ignore this many bases at the start (default 0)\n");
    fprintf(stderr,"       -e<int>    Ignore this many bases at the end (default 0)\n");
    exit(1);
  }
  UC upper[256];
  for(int i=0;i<256;i++)upper[i]=toupper(i);
  std::ifstream fp0(av[optind]);
  if(fp0.fail())error(1,errno,"Couldn't open %s",av[optind]);
  std::ifstream fp1(av[optind+1]);
  if(fp1.fail())error(1,errno,"Couldn't open %s",av[optind+1]);
  int n=0,n0=0,n1=0;
  string l,name0,name1,gen0,gen1;
  long long int dist=0,totsize=0;
  while(1){
    n++;
    while(std::getline(fp0,l)){
      if(l.size()==0)continue;
      if(l[0]=='>'){name0=l.substr(1);continue;}
      n0++;break;
    }
    if(n0<n&&n0!=1)break;
    if(n0==n)gen0=l;
    while(std::getline(fp1,l)){
      if(l.size()==0)continue;
      if(l[0]=='>'){name1=l.substr(1);continue;}
      n1++;break;
    }
    if(n1<n&&n1!=1)break;
    if(n1==n)gen1=l;
    if(n0<n&&n1<n)break;
    if(gen0.size()!=gen1.size())error(2,0,"Genomes %s and %s are of different lengths, %lu and %lu\n",name0.c_str(),name1.c_str(),gen0.size(),gen1.size());
    unsigned int i;
    long long int d=0;
    for(i=ignorestart;i<gen0.size()-ignoreend;i++)d+=(upper[UC(gen0[i])]!=upper[UC(gen1[i])]);
    if(showall)printf("Distance %8lld between %s and %s\n",d,name0.c_str(),name1.c_str());
    dist+=d;
    totsize+=gen0.size()-ignorestart-ignoreend;
  }
  printf("Read %d and %d genomes\n",n0,n1);
  printf("Total distance = %lld = %g per kbase\n",dist,dist/(totsize/1000.));
}
