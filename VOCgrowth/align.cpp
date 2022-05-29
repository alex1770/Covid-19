// Aligns fasta format input to a given reference genome

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
#include <set>
#include <unordered_set>
#include <unordered_map>
using std::string;
using std::pair;
using std::tie;
using std::vector;
using std::set;
using std::unordered_set;
using std::unordered_map;
using std::min;
using std::max;

typedef unsigned char UC;
typedef unsigned int UI;
typedef long long int int64;

// R = number of bases from which the indexes are formed
#define R 10
vector<int> refdict[1<<R*2];

const int undefined=0x7f7f7f7f;
const int infinity=1000000000;

template<class T> struct array2d {
  size_t rows,cols;
  vector<T> data;
  array2d(){}
  array2d(size_t r, size_t c):rows(r),cols(c),data(r*c){}
  void size(size_t r, size_t c){rows=r;cols=c;data.resize(r*c);}
  T* operator[](size_t index){return &data[index*cols];}// First level indexing
};

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

// Not currently used
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
// and extract the ID, e.g., "hCoV-19/env/Austria/CeMM11657/2021", and possibly prepend with given prefix (to put COG-UK on the same footing as GISAID)
// "" means not available
string parseheader(const string &idprefix,const string &header){
  assert(header.size()>0&&header[0]=='>');
  vector<string> hs=split(header,"|",false,1);
  if(hs.size()>0){
    if(idprefix=="")return hs[0];
    return idprefix+"/"+hs[0];
  }
  return "";
}

// IDs_present file format - just a simple list of IDs
unordered_set<string> readIDs(string dir){
  unordered_set<string> done;
  string fname=dir+"/IDs_present";
  std::ifstream fp(fname);// Use read-lock
  if(fp.fail())return done;
  string l;
  while(std::getline(fp,l))done.insert(l);
  fp.close();
  int n=done.size();
  fprintf(stderr,"Read %d genome ID%s from %s\n",n,n==1?"":"s",fname.c_str());
  return done;
}

// IDs_present file format - just a simple list of IDs
void writeIDs(string dir,unordered_set<string> done){
  string fname=dir+"/IDs_present";
  std::ofstream fp(fname);// Use write-lock
  if(fp.fail()){fprintf(stderr,"Couldn't write %s\n",fname.c_str());return;}
  vector<string> vs;
  for(const string &s:done)vs.push_back(s);
  std::sort(vs.begin(),vs.end());
  for(string &s:vs)fp<<s<<"\n";
  fp.close();
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
  string reffn="refgenome";
  string idprefix,datadir;
  int compression=0,minoffsetcount=20;
  while(1)switch(getopt(ac,av,"c:p:m:r:x:")){
    case 'c': compression=atoi(optarg);break;
    case 'm': minoffsetcount=atoi(optarg);break;
    case 'p': idprefix=strdup(optarg);break;
    case 'r': reffn=strdup(optarg);break;
    case 'x': datadir=strdup(optarg);break;
    case -1: goto ew0;
    default: goto err0;
  }
 ew0:
  if(optind<ac){
  err0:
    fprintf(stderr,"Usage: align [options]\n");
    fprintf(stderr,"       -c<int>    Compression mode (0=default=uncompressed fasta output)\n");
    fprintf(stderr,"       -m<int>    minoffsetcount (default 20)\n");
    fprintf(stderr,"       -p<string> ID prefix (e.g., \"hCoV-19\" to put COG-UK on same footing as GISAID)\n");
    fprintf(stderr,"       -r<string> Reference genome fasta file (default \"refgenome\")\n");
    fprintf(stderr,"       -x<string> Data directory\n");
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
  
  int i,j,t;
  
  string refgenome;
  int N;
  {
    std::ifstream fp(reffn);
    if(fp.fail())error(1,errno,"\nCouldn't open %s",reffn.c_str());
    while(std::getline(fp,refgenome))if(refgenome.size()>0&&refgenome[0]!='>')break;      
    fp.close();
    N=refgenome.size();
  }

  // refdict[R-tuple of bases, X] = list of positions in the reference genome with that R-tuple, X
  int base2num[256];
  for(i=0;i<256;i++)base2num[i]=-1;
  base2num['A']=base2num['a']=0;
  base2num['C']=base2num['c']=1;
  base2num['G']=base2num['g']=2;
  base2num['T']=base2num['t']=3;
  int badj=-1;
  for(j=0,t=0;j<N;j++){
    int b=base2num[refgenome[j]&255];
    if(b<0){badj=j;continue;}
    t=t>>2|b<<(R-1)*2;
    if(j>=badj+R){assert(t>=0&&t<(1<<R*2));refdict[t].push_back(j-(R-1));}
  }

  vector<int> j2num_i(N), j2ind_i(N), list_i;

  int linenum=0,nwrite=0;
  bool last;
  string header,line;
  vector<UC> genome;
  last=!std::getline(std::cin,header);linenum++;
  int skip[3]={0};
  while(!last){// Main loop
    tick(1);
    string id=parseheader(idprefix,header);
    if(id=="")skip[0]++; else if(done.count(id))skip[1]++; else if(datadir!=""&&id2date.count(id)==0)skip[2]++; else goto ok0;
    // Skip genome we already have, or one for which the date isn't known
    while(1){
      last=!std::getline(std::cin,header);linenum++;
      if(last)break;
      int s=header.size();
      if(s>0&&header[0]=='>')break;
    }
    tock(1);
    continue;
  ok0:
    tock(1);
    tick(2);
    string date;
    if(datadir!="")date=id2date[id];
    done.insert(id);
    int M=0;// Length of genome;
    while(1){
      last=!std::getline(std::cin,line);linenum++;
      if(last)break;
      int s=line.size();
      if(s>0&&line[0]=='>')break;
      genome.resize(M+s);
      memcpy(&genome[M],&line[0],s);M+=s;
    }
    tock(2);
    
    // Offset = (index in ref genome) - (index in current genome) = j-i
    // offsetcount[M+offset] = number of possible uses of this offset (assuming any R-tuple in genome can match same R-tuple anywhere in reference genome)
    // indexkey[i] = R-tuple of bases at position i in the current genome
    vector<int> offsetcount(M+N);
    vector<int> indexkey(M);
    int badi=-1;
    tick(3);
    for(i=0,t=0;i<M;i++){
      int b=base2num[genome[i]];
      if(b<0){badi=i;continue;}
      int i1=i-(R-1);
      t=t>>2|b<<(R-1)*2;
      if(i1>badi){
        assert(t>=0&&t<(1<<R*2));
        indexkey[i1]=t;
        for(int j:refdict[t])offsetcount[M+j-i1]++;
      }else if(i>=R-1)indexkey[i1]=undefined;
    }
    //for(int o=0;o<M+N;o++)if(offsetcount[o]>=minoffsetcount)fprintf(stderr,"Offset %d %d\n",o-M,offsetcount[o]);
    tock(3);

    tick(4);
    vector<int> pointoffset_i(M,undefined),pointoffset(N,undefined);
    vector<int> best(N);
    for(i=0;i<=M-R;i++){
      t=indexkey[i];
      if(t!=undefined){
        int best_i=minoffsetcount-1;
        for(int j:refdict[t]){
          int c=offsetcount[M+j-i];
          if(c>best[j]){best[j]=c;pointoffset[j]=j-i;}
          if(c>best_i){best_i=c;pointoffset_i[i]=j-i;}
        }
      }
    }
    tock(4);
    tick(5);
    // Build up to 4 possible offsets, i2j[i][0,1], j2i[j][0,1], to use at each position (i in the current genome, j in ref)
    // i2j[][] works well for deletions
    // j2i[][] works well for insertions
    array2d<int> i2j(M,2),j2i(N,2);
    // Approach from right
    int nearest=undefined;
    for(i=M-1;i>M-R;i--)i2j[i][1]=undefined;
    for(i=M-R;i>=0;i--){
      if(pointoffset_i[i]!=undefined)nearest=pointoffset_i[i];
      if(nearest!=undefined)i2j[i][1]=i+nearest; else i2j[i][1]=undefined;
    }
    nearest=undefined;
    for(j=N-1;j>N-R;j--)j2i[j][1]=undefined;
    for(j=N-R;j>=0;j--){
      if(best[j]>=minoffsetcount)nearest=pointoffset[j];
      if(nearest!=undefined)j2i[j][1]=j-nearest; else j2i[j][1]=undefined;
    }
    // Approach from left
    nearest=undefined;
    for(i=0;i<R-1;i++)i2j[i][0]=undefined;
    for(i=0;i<=M-R;i++){
      if(pointoffset_i[i]!=undefined)nearest=pointoffset_i[i];
      if(nearest!=undefined)i2j[i+R-1][0]=i+R-1+nearest; else i2j[i+R-1][0]=undefined;
    }
    nearest=undefined;
    for(j=0;j<R-1;j++)j2i[j][0]=undefined;
    for(j=0;j<=N-R;j++){
      if(best[j]>=minoffsetcount)nearest=pointoffset[j];
      if(nearest!=undefined)j2i[j+R-1][0]=j+R-1-nearest; else j2i[j+R-1][0]=undefined;
    }
    tock(5);

    // Make antichains - make an linear order of all allowable (i,j) such that a later (i,j) is never less-in-the-partial-order than an earlier one.
    tick(9);
    memset(&j2num_i[0],0,N*sizeof(int));
    for(j=0;j<N;j++){
      int i0=j2i[j][0],i1=j2i[j][1];
      if(        i0>=0&&i0<M)j2num_i[j]++;
      if(i1!=i0&&i1>=0&&i1<M)j2num_i[j]++;
    }
    for(i=0;i<M;i++){
      int j0=i2j[i][0],j1=i2j[i][1];
      if(        j0>=0&&j0<N&&i!=j2i[j0][0]&&i!=j2i[j0][1])j2num_i[j0]++;
      if(j1!=j0&&j1>=0&&j1<N&&i!=j2i[j1][0]&&i!=j2i[j1][1])j2num_i[j1]++;
    }
    int tot=0;
    for(j=0;j<N;j++){
      j2ind_i[j]=tot;
      tot+=j2num_i[j];
    }
    fprintf(stderr,"Total %6d   Ratio=%g\n",tot,tot/double(max(M,N)));
    list_i.resize(tot);
    for(j=0;j<N;j++){
      int i0=j2i[j][0],i1=j2i[j][1];
      if(        i0>=0&&i0<M)list_i[j2ind_i[j]++]=i0;
      if(i1!=i0&&i1>=0&&i1<M)list_i[j2ind_i[j]++]=i1;
    }
    for(i=0;i<M;i++){
      int j0=i2j[i][0],j1=i2j[i][1];
      if(        j0>=0&&j0<N&&i!=j2i[j0][0]&&i!=j2i[j0][1])list_i[j2ind_i[j0]++]=i;
      if(j1!=j0&&j1>=0&&j1<N&&i!=j2i[j1][0]&&i!=j2i[j1][1])list_i[j2ind_i[j1]++]=i;
    }
    tock(9);
    
    tick(10);
    for(j=0,tot=0;j<N;j++){
      j2ind_i[j]=tot;
      tot+=j2num_i[j];
    }
    vector<int> mintree(M*2+50);
    // mintree[] is a binary tree to do min-query and range-min-update O(logn) time. (Could do O(1) time if feeling energetic.)
    // It implements an array val[0...M-1] with queries of the form val[i1] and updates of the form {val[i]=min(val[i],v0) for all i>i1}.
    // At stage j, val[i]+i+j represents the best(lowest) score achievable if you start at (i,j) and descend by legal jumps to (-1,-1),
    // but not including the mutation penalty for (i,j) itself.
    // You could also say, at stage j*, val[i]+i+j = best score achievable from (i,j) for j>=j*, using only waypoints <j*.
    // A legal jump is a move of the form (i,j) -> (i',j'), where i'<i, j'<j.
    // Such moves (with -1<=i'<i, -1<=j'<j) incur
    // (i) a skip penalty of (j-1-j')+(i-1-i')+(2 if (i',j') isn't in the list), and
    // (ii) a mutation penalty of (refgenome[j]!=genome[i])*C for some C to be decided on.
    // For these purposes, the list is deemed to include (-1,-1) and (M,N).
    vector<int> nbp0(tot),nbp1(tot);
    for(j=0;j<N;j++){
      int k;
      if(j2num_i[j]>1)std::sort(&list_i[j2ind_i[j]],&list_i[j2ind_i[j]+j2num_i[j]],std::greater<>());
      for(k=0;k<j2num_i[j];k++){
        int i1=list_i[j2ind_i[j]+k];
        assert(i1>=0&&i1<M);
        // See if waypoint (i1,j) improves val_{j+1}(i) for some i>i1, otherwise it will be left with its value based on earlier waypoints (*,<j)
        int mi=infinity;
        {
          int i=i1,m=M+1,p=0;
          do{mi=min(mi,mintree[p+i]);p+=m;m=(m+1)>>1;i=i>>1;}while(m>1&&i<m);
        }
        int v0=mi-2+(refgenome[j]!=genome[i1])*2;
        nbp0[j2ind_i[j]+k]=i1;
        nbp1[j2ind_i[j]+k]=v0;
        {
          int i=i1+1,m=M+1,p=0;
          do{mintree[p+i]=min(mintree[p+i],v0);p+=m;m=(m+1)>>1;i=(i+1)>>1;}while(m>1&&i<m);
        }
      }
    }
    tock(10);

    tick(11);
    // Write aligned genome, out[]
    vector<UC> out(N+1);
    memset(&out[0],'-',N);
    out[N]=0;
    {
      int k,vl=infinity;
      {
        int i=M,m=M+1,p=0;
        do{vl=min(vl,mintree[p+i]);p+=m;m=(m+1)>>1;i=i>>1;}while(m>1&&i<m);
      }
      int i1=M;
      for(j=N-1;j>=0;j--){
        for(k=0;k<j2num_i[j];k++){
          int p=j2ind_i[j]+k;
          int i=list_i[p];
          if(i<i1)assert(nbp1[p]>=vl);
          if(i<i1&&nbp1[p]==vl){
            vl-=-2+(refgenome[j]!=genome[i])*2;
            i1=i;
            out[j]=genome[i];
            break;
          }
        }
      }
    }
    tock(11);
    
    tick(8);
    FILE*fp;
    if(datadir=="")fp=stdout; else fp=fopen((datadir+"/"+date).c_str(),"a");
    header=line;// Next header is the last-read line
    switch(compression){
    case 0:
      fprintf(fp,">%s",id.c_str());
      if(date!="")fprintf(fp,"|%s",date.c_str());
      fprintf(fp,"\n%s\n",&out[0]);
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
    tock(8);
    nwrite++;
  }
  if(datadir!="")writeIDs(datadir,done);
  fprintf(stderr,"Wrote %d new genome%s\n",nwrite,nwrite==1?"":"s");
  fprintf(stderr,"Skipped %d because ID could not be read\n",skip[0]);
  fprintf(stderr,"Skipped %d because already stored\n",skip[1]);
  fprintf(stderr,"Skipped %d because date was not available\n",skip[2]);
  prtim();
}
