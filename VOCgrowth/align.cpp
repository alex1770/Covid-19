// Aligns SARS-CoV-2 fasta files to reference genome Wuhan-Hu-1
// Todo: handle files with \r\n lines
//       multithread if need for speed arises

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
#include <unordered_set>
using std::string;
using std::pair;
using std::tie;
using std::vector;
using std::unordered_set;

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

// Read headers from everything under 'path' into 'headers'
// (But actually, maybe this only gets used on my format fasta files)
// GISAID fasta headers look like this:
// >hCoV-19/England/PHEC-L303L83F/2021|2021|2021-06-30
// or occasionally like this:
// >hCoV-19/env/Austria/CeMM11657/2021|2021-06-14|2021-07-01
// COG-UK fasta headers look like this:
// >England/PHEC-YYF8DBE/2022

// Takes an input like this
// >hCoV-19/England/PHEC-L303L83F/2021|2021|2021-06-30
// >hCoV-19/env/Austria/CeMM11657/2021|2021-06-14|2021-07-01
// >England/PHEC-YYF8DBE/2022
// and returns (ID,sample date), e.g.,
// ("England/PHEC-L303L83F/2021","2021-06-30")
// Sample date is "" if this date is not available
pair<string,string> parseheader(string head){
  assert(head.size()>0&&head[0]=='>');
  vector<string> hs=split(head,">/|");
  if(hs.size()>=3){
    UI i=0;
    if(hs[i]=="hCoV-19")i++;
    if(hs[i]=="env"||hs[i]=="cat")i++;// check - alter
    if(hs.size()>=i+3){
      string id=hs[i]+"/"+hs[i+1]+"/"+hs[i+2];
      string date;
      if(hs.size()>=i+4&&hs[i+3].size()==10&&hs[i+3][0]=='2')date=hs[i+3];
      return {id,date};
    }
  }
  return {"",""};
}

unordered_set<string> headers;
void readheaders(string &path){
  auto fn=[](const char*fpath,const struct stat*sb,int typeflag)->int{
    if(typeflag==FTW_F){
      std::ifstream fp(fpath);
      if(fp.fail())error(1,errno,"\nCouldn't open %s",fpath);
      string l;
      while(std::getline(fp,l)){
        if(l.size()>0&&l[0]=='>'){
          string id,date;
          tie(id,date)=parseheader(l);
          if(id!=""&&date!="")headers.insert(id);
        }
      }
      fp.close();
    }
    return 0;
  };
  ftw(path.c_str(),fn,100);
  int n=headers.size();
  printf("Read %d genome header%s\n",n,n==1?"":"s");
}

int main(int ac,char**av){
  string datadir;
  while(1)switch(getopt(ac,av,"x:")){
    case 'x': datadir=strdup(optarg);break;
    case -1: goto ew0;
    default: goto err0;
  }
 ew0:
  if(optind<ac){
  err0:
    fprintf(stderr,"Usage: align [options]\n");
    fprintf(stderr,"       -x<string> Data directory\n");
    exit(1);
  }

  if(datadir!=""){
    mkdir(datadir.c_str(),0777);
    readheaders(datadir);
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
  string name;
  UC genome[MAXGS];
  last=!std::getline(std::cin,name);linenum++;
  while(!last){
    string id,date;
    tie(id,date)=parseheader(name);
    if(headers.count(id)||date==""){
      // Skip genome we already have, or one for which the date isn't known
      while(1){
        last=!std::getline(std::cin,name);linenum++;
        if(last)break;
        int s=name.size();
        if(s>0&&name[0]=='>')break;
      }
      continue;
    }
    FILE*fp;
    if(datadir=="")fp=stdout; else fp=fopen((datadir+"/"+date).c_str(),"a");
    fprintf(fp,"%s\n",name.c_str());
    int M=0;// Length of genome;
    while(1){
      last=!std::getline(std::cin,name);linenum++;
      if(last)break;
      int s=name.size();
      if(s>0&&name[0]=='>')break;
      if(M+s>MAXGS){fprintf(stderr,"Warning: Overlong genome at line %d\n",linenum);continue;}
      memcpy(genome+M,&name[0],s);M+=s;
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

    UC out[MAXGS+1];
    memset(out,'-',N);
    out[N]=0;
    int s=(st[1]<st[0]);
    for(i=M-1;i>=0;i--){
      int o=offsets[i][s];
      if(o!=undefined&&i+o>=0&&i+o<N)out[i+o]=genome[i];
      s=bp[i][s];
    }
    fprintf(fp,"%s\n",out);
    if(datadir!="")fclose(fp);
    nwrite++;
  }
  if(datadir!="")printf("Wrote %d new genome%s\n",nwrite,nwrite==1?"":"s");
}
