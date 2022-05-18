// Aligns SARS-CoV-2 fasta files to reference genome Wuhan-Hu-1
// Todo: handle files with \r\n lines?

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <getopt.h>
#include <error.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <ftw.h>
#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <set>
#include <tuple>
#include <algorithm>
using std::string;
using std::pair;
using std::tie;
using std::vector;
using std::set;

typedef unsigned char UC;
typedef unsigned int UI;

// R = number of bases from which the indexes are formed
#define R 9
vector<int> refdict[1<<R*2];

// Maximum number of bases
#define MAXGS 40000

// Count threshold for offsets
#define MINOFFSETCOUNT 20

const int undefined=0x7fffffff;
const int infinity=1000000000;

double cpu(){return clock()/double(CLOCKS_PER_SEC);}
int timings=1;
#define MAXTIM 50
double ncpu[MAXTIM]={0},lcpu[MAXTIM]={0},tcpu[MAXTIM]={0};
void tick(int i){if(timings)lcpu[i]=cpu();}
void tock(int i){if(timings){ncpu[i]+=1;tcpu[i]+=cpu()-lcpu[i];}}
void prtim(){
  int i;
  double x=(ncpu[0]>0?tcpu[0]/ncpu[0]:0);
  for(i=1;i<MAXTIM;i++)if(ncpu[i]){
    double t=tcpu[i]-ncpu[i]*x;
    printf("Time %2d: CPU %12gs / %12g = %12gus\n",i,t,ncpu[i],t/ncpu[i]*1e6);
  }
}

// Split string into a sequence of substrings using any character from sep (default whitespace) as a separator.
vector<string> split(string in,string sep=" \r\t\n\f\v"){
  size_t i,j,p=0;
  vector<string> rv;
  while(1){
    i=in.find_first_not_of(sep,p);if(i==string::npos)i=in.size();
    j=in.find_first_of(sep,i);if(j==string::npos)j=in.size();
    if(i==j)return rv;
    rv.push_back(in.substr(i,j-i));
    p=j;
  }
}

// Take an input like this
// (GISAID style, unknown date)  >hCoV-19/England/PHEC-L303L83F/2021|2021|2021-06-30
// (GISAID style, known date)    >hCoV-19/Austria/CeMM11657/2021|2021-06-14|2021-07-01
// (GISAID style, extra prefix)  >hCoV-19/env/Austria/CeMM11657/2021|2021-06-14|2021-07-01
// (COG-UK style, no date)       >England/PHEC-YYF8DBE/2022
// and return (ID,sample date), e.g.,
// ("England/PHEC-L303L83F/2021","2021-06-30"), or "" means not available
pair<string,string> parseheader(string &header){
  assert(header.size()>0&&header[0]=='>');
  vector<string> hs=split(header,">|");
  string id,date;
  if(hs.size()>0){
    vector<string> ida=split(hs[0],"/");
    int n=ida.size();
    if(n>=3)id=ida[n-3]+"/"+ida[n-2]+"/"+ida[n-1];
    if(hs.size()>=3&&hs[2].size()==10&&hs[2][0]=='2'&&hs[2][1]=='0')date=hs[2];
  }
  return {id,date};
}

void readheaders(string &dir,vector<pair<string,string>>&index){
  string fname=dir+"/index";
  std::ifstream fp(fname);// Use read-lock
  if(fp.fail())return;
  string l;
  while(std::getline(fp,l))index.push_back({l.substr(0,10),l.substr(11)});
  fp.close();
  int n=index.size();
  printf("Read %d genome ID%s from %s\n",n,n==1?"":"s",fname.c_str());
}

void writeheaders(string &dir,vector<pair<string,string>>&index){
  std::ofstream fp(dir+"/index");// Use write-lock
  for(auto &ind:index)fp<<ind.first<<" "<<ind.second<<"\n";
  fp.close();
}

int main(int ac,char**av){
  string datadir;
  int compression=0;
  while(1)switch(getopt(ac,av,"c:x:")){
    case 'c': compression=atoi(optarg);break;
    case 'x': datadir=strdup(optarg);break;
    case -1: goto ew0;
    default: goto err0;
  }
 ew0:
  if(optind<ac){
  err0:
    fprintf(stderr,"Usage: align [options]\n");
    fprintf(stderr,"       -x<string> Data directory\n");
    fprintf(stderr,"       -c<int>    Compression mode (0=default=standard fasta=uncompressed ACGT)\n");
    exit(1);
  }

  // It seems you need this otherwise std::getline will be ridiculously slow because it synchronises to C stdio
  std::ios_base::sync_with_stdio(false);
  
  vector<pair<string,string>> index;
  set<string> done;
  if(datadir!=""){
    mkdir(datadir.c_str(),0777);
    readheaders(datadir,index);
    string id,date;
    for(auto &ind:index){
      tie(id,date)=parseheader(ind.second);
      done.insert(id);
    }
  }
  
  int i,t;
  
  string refgenome;
  int N;
  {
    const char*reffn="refgenome";
    std::ifstream fp(reffn);
    if(fp.fail())error(1,errno,"\nCouldn't open %s",reffn);
    std::getline(fp,refgenome);
    fp.close();
    N=refgenome.size();
    assert(N<=MAXGS);
  }

  int base2num[256];
  for(i=0;i<256;i++)base2num[i]=-1;
  base2num['A']=base2num['a']=0;
  base2num['C']=base2num['c']=1;
  base2num['G']=base2num['g']=2;
  base2num['T']=base2num['t']=3;
  for(i=0,t=0;i<N;i++){
    t=t>>2|base2num[refgenome[i]&255]<<(R-1)*2;
    if(i>=R-1&&t>=0&&t<(1<<R*2))refdict[t].push_back(i-(R-1));
  }

  int jumppen[MAXGS+1];
  for(t=0;t<=MAXGS;t++)jumppen[t]=int(floor(sqrt(t)+1e-6));

  int linenum=0,nwrite=0;
  bool last;
  string header,line;
  UC genome[MAXGS];
  last=!std::getline(std::cin,header);linenum++;
  while(!last){
    string id,date;
    tie(id,date)=parseheader(header);
    if(done.count(id)||date==""){
      // Skip genome we already have, or one for which the date isn't known
      while(1){
        last=!std::getline(std::cin,header);linenum++;
        if(last)break;
        int s=header.size();
        if(s>0&&header[0]=='>')break;
      }
      continue;
    }
    index.push_back({date,header});
    done.insert(id);
    int M=0;// Length of genome;
    while(1){
      last=!std::getline(std::cin,line);linenum++;
      if(last)break;
      int s=line.size();
      if(s>0&&line[0]=='>')break;
      if(M+s>MAXGS){fprintf(stderr,"Warning: Overlong genome at line %d\n",linenum);continue;}
      memcpy(genome+M,&line[0],s);M+=s;
    }
    
    // Offset = (index in ref genome) - (index in current genome)
    UC offsetcount[MAXGS*2]={0};
    int indexkey[MAXGS]={0};
    int badi=-1;
    for(i=0,t=0;i<M;i++){
      int b=base2num[genome[i]];
      if(b<0){badi=i;continue;}
      t=t>>2|b<<(R-1)*2;
      if(i>=badi+R){
        assert(t>=0&&t<(1<<R*2));
        indexkey[i-(R-1)]=t;
        for(int x:refdict[t]){
          int k=x-(i-(R-1));
          if(offsetcount[MAXGS+k]<MINOFFSETCOUNT)offsetcount[MAXGS+k]++;
        }
      }else if(i>=R-1)indexkey[i-(R-1)]=undefined;
    }
    //if(!ok){fprintf(stderr,"Warning: Can't find offsets for genome at lines preceding %d\n",linenum);break;}
    //for(int o=0;o<MAXGS*2;o++)if(offsetcount[o]==MINOFFSETCOUNT)printf("Offset %d\n",o-MAXGS);

    int pointoffset[MAXGS],offsets[MAXGS][2]={0};
    for(i=0;i<=M-R;i++){
      t=indexkey[i];
      pointoffset[i]=undefined;
      if(t!=undefined)for(int x:refdict[t])if(offsetcount[MAXGS+x-i]>=MINOFFSETCOUNT)pointoffset[i]=x-i;
    }
    // Approach from right
    int nearest=undefined;
    for(i=M-1;i>M-R;i--)offsets[i][1]=undefined;
    for(i=M-R;i>=0;i--){
      if(pointoffset[i]!=undefined)nearest=pointoffset[i];
      offsets[i][1]=nearest;
    }
    // Approach from left
    nearest=undefined;
    for(i=0;i<R-1;i++)offsets[i][0]=undefined;
    for(i=0;i<=M-R;i++){
      if(pointoffset[i]!=undefined)nearest=pointoffset[i];
      offsets[i+R-1][0]=nearest;
    }

    // Dyn prog on two allowable offsets: offsets[i][]
    int bp[MAXGS][2]={0};// Back pointers; bp[i][j] is defined if value is finite
    int st[2]={0,0};// State: st[j] = score (lower is better) given ended with offset offsets[i-1][j]
    for(i=0;i<M;i++){
      // Transition j -> k,  j=prev offset index, k=current offset index
      int k,newst[2]={infinity,infinity};
      for(k=0;k<2;k++){
        int j;
        int off=offsets[i][k];
        if(off!=undefined){
          int best=infinity,bestj=0;
          for(j=0;j<2;j++){
            int v=0;
            if(i>0){// Initial jump is free
              int p=offsets[i-1][j];
              if(p==undefined)continue;
              v=jumppen[abs(off-p)];
            }
            v+=st[j]+(i+off<0||i+off>=N||genome[i]!=refgenome[i+off]);
            if(v<best){best=v;bestj=j;}
          }
          newst[k]=best;
          bp[i][k]=bestj;
        }
      }
      st[0]=newst[0];
      st[1]=newst[1];
    }

    // Write aligned genome, out[]
    UC out[MAXGS+1];
    memset(out,'-',N);
    out[N]=0;
    int s=(st[1]<st[0]);
    for(i=M-1;i>=0;i--){
      int o=offsets[i][s];
      if(o!=undefined&&i+o>=0&&i+o<N)out[i+o]=genome[i];
      s=bp[i][s];
    }

    FILE*fp;
    if(datadir=="")fp=stdout; else fp=fopen((datadir+"/"+date).c_str(),"a");
    fprintf(fp,"%s|C%d\n",header.c_str(),compression);
    header=line;// Next header is the last-read line
    switch(compression){
    case 0:
      fprintf(fp,"%s\n",out);
      break;
    case 1:
      for(i=0;i<N;i++){
        int j;
        if(out[i]==refgenome[i])continue;
        if(out[i]=='N'||out[i]=='-'){
          for(j=i;j<N&&out[j]==out[i];j++);
          fprintf(fp,"%d-%d %c\n",i,j-1,out[i]);
          i=j;
        }else{
          fprintf(fp,"%d %c\n",i,out[i]);
        }
      }
      break;
    default:
      error(1,0,"Unknown compression type %d\n",compression);
    }
      
    if(datadir!="")fclose(fp);
    nwrite++;
  }
  if(datadir!=""){
    std::sort(index.begin(),index.end());
    writeheaders(datadir,index);
    printf("Wrote %d new genome%s\n",nwrite,nwrite==1?"":"s");
  }
  prtim();
}
