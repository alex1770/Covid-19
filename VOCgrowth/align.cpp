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
#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <tuple>
#include <algorithm>
#include <unordered_set>
#include <unordered_map>
using std::string;
using std::pair;
using std::tie;
using std::vector;
using std::unordered_set;
using std::unordered_map;

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
    fprintf(stderr,"Time %2d: CPU %12gs / %12g = %12gus\n",i,t,ncpu[i],t/ncpu[i]*1e6);
  }
}

// Split string into a sequence of substrings using any character from sep (default whitespace) as a separator.
vector<string> split(string in,string sep=" \r\t\n\f\v",bool ignoreempty=false,size_t startpos=0){
  size_t i,j=startpos;
  vector<string> rv;
  // Imagine separators at -1 and n; j points to 1 after the previous separator
  while(1){
    if(ignoreempty){
      i=in.find_first_not_of(sep,j);
      if(i==string::npos)return rv;
    }else i=j;
    j=in.find_first_of(sep,i);if(j==string::npos)j=in.size();
    rv.push_back(in.substr(i,j-i));
    if(j==in.size())return rv;
    j++;
  }
}

string getid(string gisaidname){
  vector<string> ida=split(gisaidname,"/");
  int n=ida.size();
  if(n>=3)return ida[n-3]+"/"+ida[n-2]+"/"+ida[n-1];
  return "";
}

// Take an input like this
// (GISAID style, unknown date)  >hCoV-19/England/PHEC-L303L83F/2021|2021|2021-06-30
// (GISAID style, known date)    >hCoV-19/Austria/CeMM11657/2021|2021-06-14|2021-07-01
// (GISAID style, extra prefix)  >hCoV-19/env/Austria/CeMM11657/2021|2021-06-14|2021-07-01
// (COG-UK style, no date)       >England/PHEC-YYF8DBE/2022
// and extract the ID, e.g., "England/PHEC-L303L83F/2021" or "Austria/CeMM11657/2021". "" means not available
string parseheader(string &header){
  assert(header.size()>0&&header[0]=='>');
  vector<string> hs=split(header,"|",false,1);
  if(hs.size()>0)return getid(hs[0]);
  return "";
}

// IDs_present file format:
// Space-separated. No spaces within fields.
// date ID
unordered_set<string> readIDs(string dir){
  unordered_set<string> done;
  string fname=dir+"/IDs_present";
  std::ifstream fp(fname);// Use read-lock
  if(fp.fail())return done;
  string l;
  while(std::getline(fp,l)){
    vector<string> ll=split(l);
    done.insert(ll[1]);
  }
  fp.close();
  int n=done.size();
  fprintf(stderr,"Read %d genome ID%s from %s\n",n,n==1?"":"s",fname.c_str());
  return done;
}

// metadata file format:
// Tab-separated, like GISAID.
// date ID provisionallineage location
// location is "Continent / Country[ / subsidiary [ / subsidiary ...]]"
unordered_map<string,string> readmeta(string dir){
  unordered_map<string,string> id2date;
  string fname=dir+"/metadata.tsv";
  std::ifstream fp(fname);// Use read-lock
  if(fp.fail())return id2date;
  string l;
  std::getline(fp,l);
  vector<string> headers=split(l,"\t");
  int idcol=-1,datecol=-1;
  for(UI i=0;i<headers.size();i++){
    if(headers[i]=="ID")idcol=i;
    if(headers[i]=="date")datecol=i;
  }
  assert(idcol>=0&&datecol>=0);
  while(std::getline(fp,l)){
    vector<string> ll=split(l,"\t");
    id2date[ll[idcol]]=ll[datecol];
  }
  fp.close();
  int n=id2date.size();
  fprintf(stderr,"Read %d metadata entr%s from %s\n",n,n==1?"y":"ies",fname.c_str());
  return id2date;
}

bool okdate(string date){
  return date.size()==10&&date[0]=='2'&&date[1]=='0'&&date[4]=='-'&&date[7]=='-';
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
    fprintf(stderr,"       -c<int>    Compression mode (0=default=uncompressed fasta output)\n");
    exit(1);
  }

  // It seems you need this otherwise std::getline will be ridiculously slow because it synchronises to C stdio
  std::ios_base::sync_with_stdio(false);

  unordered_set<string> done;
  unordered_map<string,string> id2date;
  if(datadir!=""){
    mkdir(datadir.c_str(),0777);
    done=readIDs(datadir);
    id2date=readmeta(datadir);
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

  // refdict[R-tuple of bases, X] = list of positions in the reference genome with that R-tuple, X
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
  int skip[3]={0};
  while(!last){// Main loop
    string id=parseheader(header);
    if(id=="")skip[0]++; else if(done.count(id))skip[1]++; else if(datadir!=""&&id2date.count(id)==0)skip[2]++; else goto ok0;
    // Skip genome we already have, or one for which the date isn't known
    while(1){
      last=!std::getline(std::cin,header);linenum++;
      if(last)break;
      int s=header.size();
      if(s>0&&header[0]=='>')break;
    }
    continue;
  ok0:
    string date;
    if(datadir!="")date=id2date[id];
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
    // offsetcount[MAXGS+offset] = number of possible uses of this offset (assuming any R-tuple in genome can match same R-tuple anywhere in reference genome)
    // indexkey[i] = R-tuple of bases at position i in the current genome
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
    //for(int o=0;o<MAXGS*2;o++)if(offsetcount[o]==MINOFFSETCOUNT)fprintf(stderr,"Offset %d\n",o-MAXGS);

    // Build two possible offsets, offsets[i][], to use at each position, i in the current genome.
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

    // Dyn prog on the two allowable offsets: offsets[i][]
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
            if(v<=best){best=v;bestj=j;}
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
    header=line;// Next header is the last-read line
    switch(compression){
    case 0:
      fprintf(fp,"%s",id.c_str());
      if(date!="")fprintf(fp,"|%s",date.c_str());
      fprintf(fp,"\n%s\n",out);
      break;
    case 1:
      fprintf(fp,"%s|C%d",id.c_str(),compression);
      for(i=0;i<N;i++){
        int j;
        if(out[i]==refgenome[i])continue;
        // Output in count-from-1 notation
        if(out[i]=='N'||out[i]=='-'){
          for(j=i;j<N&&out[j]==out[i];j++);
          fprintf(fp,"|%d-%d%c",i+1,j,out[i]);
          i=j;
        }else{
          fprintf(fp,"|%d%c",i+1,out[i]);
        }
      }
      fprintf(fp,"\n");
      break;
    case 2:
      fprintf(fp,"%s|C%d|",id.c_str(),compression);
      int p;
      p=0;
      for(i=0;i<N;i++){
        int j;
        if(out[i]==refgenome[i])continue;
        if(out[i]=='N'||out[i]=='-'){
          for(j=i;j<N&&out[j]==out[i];j++);
          fprintf(fp,"%d+%d%c",i-p,j-i,out[i]);
          p=i=j;
        }else{
          fprintf(fp,"%d%c",i-p,out[i]);
          p=i+1;
        }
      }
      fprintf(fp,"\n");
      break;
    default:
      error(1,0,"Unknown compression type %d\n",compression);
    }
      
    if(datadir!="")fclose(fp);
    nwrite++;
  }
  fprintf(stderr,"Wrote %d new genome%s\n",nwrite,nwrite==1?"":"s");
  fprintf(stderr,"Skipped %d because ID could not be read\n",skip[0]);
  fprintf(stderr,"Skipped %d because already stored\n",skip[1]);
  fprintf(stderr,"Skipped %d because date was not available\n",skip[2]);
  prtim();
}
