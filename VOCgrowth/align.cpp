// Aligns fasta format input to a given reference genome
// Intended to be used for genomes you expect to be somewhat similar

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
using std::fill;

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
  array2d(size_t r, size_t c, const T&val):rows(r),cols(c),data(r*c,val){}
  size_t size(void){return (rows*cols)/sizeof(T);}
  void resize(size_t r, size_t c){rows=r;cols=c;data.resize(r*c);}
  void fill(const T&val){std::fill(data.begin(),data.end(),val);}
  T* operator[](size_t index){return &data[index*cols];}// First level indexing
};

double cpu(){return clock()/double(CLOCKS_PER_SEC);}
int timings=0;
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

bool okdate(const string &date){
  return date.size()==10&&date[0]=='2'&&date[1]=='0'&&date[4]=='-'&&date[7]=='-';
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
    if(okdate(ll[datecol]))id2date[ll[idcol]]=ll[datecol];
  }
  fp.close();
  int n=id2date.size();
  fprintf(stderr,"Read %d metadata entr%s from %s\n",n,n==1?"y":"ies",fname.c_str());
  return id2date;
}

// Changes runs of <minrun to undefined
void smoothpointoffset(vector<int> &po,int minrun){
  int i,i0=0,n=po.size();
  if(n==0)return;
  int cur=po[0];
  for(i=1;i<=n;i++){
    if(i==n||po[i]!=cur){
      if(i-i0<minrun)fill(po.begin()+i0,po.begin()+i,undefined);
      i0=i;
      if(i<n)cur=po[i];
    }
  }
}

// Extend runs of hits to represent whole of R-block
void extend(vector<int> &po){
  int i,n=po.size();
  for(i=n-R;i>=0;i--){
    if(po[i]!=undefined&&po[i+1]==undefined)fill(po.begin()+i+1,po.begin()+i+R,po[i]);
  }
}

void prarr(vector<int>&vv){
  int i,n=vv.size();
  for(i=0;i<n;i++)printf("%6d  %10d\n",i,vv[i]);
}

int main(int ac,char**av){
  int deb=0;
  string reffn="refgenome";
  string idprefix,datadir;
  int compression=0,minrun=3;
  double bigthr=10,smallthr=1;
  double df=0;
  while(1)switch(getopt(ac,av,"c:d:p:M:m:r:s:tx:")){
    case 'c': compression=atoi(optarg);break;
    case 'd': df=atof(optarg);break;
    case 'M': bigthr=atof(optarg);break;
    case 'm': smallthr=atof(optarg);break;
    case 'p': idprefix=strdup(optarg);break;
    case 'r': reffn=strdup(optarg);break;
    case 's': minrun=atoi(optarg);break;
    case 't': timings=1;break;
    case 'x': datadir=strdup(optarg);break;
    case -1: goto ew0;
    default: goto err0;
  }
 ew0:
  if(optind<ac){
  err0:
    fprintf(stderr,"Usage: align [options]\n");
    fprintf(stderr,"       -c<int>    Compression mode (0=default=uncompressed fasta output)\n");
    fprintf(stderr,"       -M<float>  big threshold (default 10)\n");
    fprintf(stderr,"       -m<float>  small threshold (default 1)\n");
    fprintf(stderr,"       -p<string> ID prefix (e.g., \"hCoV-19\" to put COG-UK on same footing as GISAID)\n");
    fprintf(stderr,"       -r<string> Reference genome fasta file (default \"refgenome\")\n");
    fprintf(stderr,"       -s<int>    Min run length for smoothing offsets (default 3)\n");
    fprintf(stderr,"       -t         Enable timings\n");
    fprintf(stderr,"       -x<string> Data directory\n");
    exit(1);
  }
  if(deb)fprintf(stderr,"M=%g m=%g s=%d df=%g\n",bigthr,smallthr,minrun,df);

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
  
  // refdict[R-tuple of bases, X] = list of positions in the reference genome which start that R-tuple, X
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

  vector<int> j2num_i(N), j2ind_i(N);
  vector<float> offsetcount;
  vector<int> indexkey;
  vector<double> best_j(N);
  vector<int> pointoffset_i,pointoffset_j(N);
  array2d<int> i2j,j2i;
  vector<int> list_i;
  vector<int> mintree;
  vector<int> nbp0,nbp1;
  vector<UC> out(N+1);
  int linenum=0,nwrite=0;
  bool last;
  string header,line;
  vector<UC> genome;
  last=!std::getline(std::cin,header);linenum++;
  int skip[3]={0};
  while(!last){// Main loop
    tick(0);tock(0);
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
    tick(14);
    indexkey.resize(M);
    offsetcount.resize(M+N);
    fill(offsetcount.begin(),offsetcount.end(),0);
    fill(indexkey.begin(),indexkey.end(),undefined);
    tock(14);
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
        if(refdict[t].size()){
          float x=1./refdict[t].size();
          for(int j:refdict[t])offsetcount[M+j-i1]+=x;
        }
      }
    }
    tock(3);

    tick(4);
    // First pass: work out a "backbone" of strong offsets
    pointoffset_i.resize(M);
    fill(pointoffset_i.begin(),pointoffset_i.end(),undefined);
    fill(pointoffset_j.begin(),pointoffset_j.end(),undefined);
    fill(best_j.begin(),best_j.end(),bigthr-1e-6);
    for(i=0;i<=M-R;i++){
      t=indexkey[i];
      if(t!=undefined){
        double best_i=bigthr-1e-6;
        for(int j:refdict[t]){
          double c=offsetcount[M+j-i];
          if(c>best_j[j]){best_j[j]=c;pointoffset_j[j]=j-i;}
          if(c>best_i){best_i=c;pointoffset_i[i]=j-i;}
        }
      }
    }
    smoothpointoffset(pointoffset_i,minrun);
    smoothpointoffset(pointoffset_j,minrun);
    extend(pointoffset_i);
    extend(pointoffset_j);
    tock(4);
    if(deb){
      printf("First pass\n");
      for(i=0;i<max(M,N);i++){
        printf("%6d",i);
        if(i<M)printf("  %10d",pointoffset_i[i]); else printf("           .");
        if(i<N)printf("  %10d",pointoffset_j[i]); else printf("           .");
        printf("\n");
      }
    }
    if(deb){
      printf("\n");
      for(i=0;i<max(M,N);i++){
        printf("%6d",i);
        if(i<M)printf("  %10d",pointoffset_i[i]); else printf("           .");
        if(i<N)printf("  %10d",pointoffset_j[i]); else printf("           .");
        printf("\n");
      }
    }

    tick(5);
    // Build up to 4 possible offsets, i2j[i][0,1], j2i[j][0,1], to use at each position (i in the current genome, j in ref)
    // i2j[][] works well for deletions
    // j2i[][] works well for insertions
    i2j.resize(M,3);i2j.fill(undefined);
    j2i.resize(N,3);j2i.fill(undefined);
    // Approach from right
    int nearest=undefined;
    for(i=M-1;i>=0;i--){
      if(pointoffset_i[i]!=undefined)nearest=pointoffset_i[i];
      if(nearest!=undefined)i2j[i][1]=i+nearest;
    }
    nearest=undefined;
    for(j=N-1;j>=0;j--){
      if(pointoffset_j[j]!=undefined)nearest=pointoffset_j[j];
      if(nearest!=undefined)j2i[j][1]=j-nearest;
    }
    // Approach from left
    nearest=undefined;
    for(i=0;i<M;i++){
      if(pointoffset_i[i]!=undefined)nearest=pointoffset_i[i];
      if(nearest!=undefined)i2j[i][0]=i+nearest;
    }
    nearest=undefined;
    for(j=0;j<N;j++){
      if(pointoffset_j[j]!=undefined)nearest=pointoffset_j[j];
      if(nearest!=undefined)j2i[j][0]=j-nearest;
    }
    tock(5);
    if(deb){
      for(i=0;i<max(M,N);i++){
        int k;
        printf("%6d",i);
        if(i<M)printf("  %10d",pointoffset_i[i]); else printf("           .");
        if(i<N)printf("  %10d",pointoffset_j[i]); else printf("           .");
        for(k=0;k<2;k++){
          if(i<M)printf("  %10d",i2j[i][k]-i); else printf("           .");
        }
        for(k=0;k<2;k++){
          if(i<N)printf("  %10d",i-j2i[i][k]); else printf("           .");
        }
        printf("\n");
      }
    }
    
    tick(6);
    fill(best_j.begin(),best_j.end(),smallthr-1e-6);
    for(i=0;i<=M-R;i++){
      t=indexkey[i];
      if(t!=undefined){
        double best_i=smallthr-1e-6;
        for(int j:refdict[t]){
          double c=offsetcount[M+j-i];// -df*abs(j/double(M)-i/double(N));
          if(i!=j2i[j][0]&&i!=j2i[j][1]&&c>best_j[j]){best_j[j]=c;j2i[j][2]=i;}
          if(j!=i2j[i][0]&&j!=i2j[i][1]&&c>best_i){best_i=c;i2j[i][2]=j;}
        }
      }
    }
    for(i=M-R;i>=0;i--)if(i2j[i][2]!=undefined&&i2j[i+1][2]==undefined)for(int i1=i+1;i1<i+R;i1++)i2j[i1][2]=i2j[i][2]+i1-i;
    for(j=N-R;j>=0;j--)if(j2i[j][2]!=undefined&&j2i[j+1][2]==undefined)for(int j1=j+1;j1<j+R;j1++)j2i[j1][2]=j2i[j][2]+j1-j;
    if(0){
      for(i=0;i<max(M,N);i++){
        int k;
        printf("%6d",i);
        if(i<M)printf("  %10d",pointoffset_i[i]); else printf("           .");
        if(i<N)printf("  %10d",pointoffset_j[i]); else printf("           .");
        for(k=0;k<3;k++){
          if(i<M)printf("  %10d",i2j[i][k]-i); else printf("           .");
        }
        for(k=0;k<3;k++){
          if(i<N)printf("  %10d",i-j2i[i][k]); else printf("           .");
        }
        printf("\n");
      }
    }
    tock(6);
    
    // Make antichains - make an linear order of all allowable (i,j) such that a later (i,j) is never less-in-the-partial-order than an earlier one.
    tick(9);
    fill(j2num_i.begin(),j2num_i.end(),0);
    for(j=0;j<N;j++){
      int i0=j2i[j][0],i1=j2i[j][1],i2=j2i[j][2];
      if(                i0>=0&&i0<M)j2num_i[j]++;
      if(        i1!=i0&&i1>=0&&i1<M)j2num_i[j]++;
      if(i2!=i1&&i2!=i0&&i2>=0&&i2<M)j2num_i[j]++;
    }
    for(i=0;i<M;i++){
      int j0=i2j[i][0],j1=i2j[i][1],j2=i2j[i][2];
      if(                j0>=0&&j0<N&&i!=j2i[j0][0]&&i!=j2i[j0][1]&&i!=j2i[j0][2])j2num_i[j0]++;
      if(        j1!=j0&&j1>=0&&j1<N&&i!=j2i[j1][0]&&i!=j2i[j1][1]&&i!=j2i[j1][2])j2num_i[j1]++;
      if(j2!=j1&&j2!=j0&&j2>=0&&j2<N&&i!=j2i[j2][0]&&i!=j2i[j2][1]&&i!=j2i[j2][2])j2num_i[j2]++;
    }
    int tot=0;
    for(j=0;j<N;j++){
      j2ind_i[j]=tot;
      tot+=j2num_i[j];
    }
    list_i.resize(tot);
    if(deb)fprintf(stderr,"Total %6d   Ratio=%g\n",tot,tot/double(max(M,N)));
    for(j=0;j<N;j++){
      int i0=j2i[j][0],i1=j2i[j][1],i2=j2i[j][2];
      if(                i0>=0&&i0<M)list_i[j2ind_i[j]++]=i0;
      if(        i1!=i0&&i1>=0&&i1<M)list_i[j2ind_i[j]++]=i1;
      if(i2!=i1&&i2!=i0&&i2>=0&&i2<M)list_i[j2ind_i[j]++]=i2;
    }
    for(i=0;i<M;i++){
      int j0=i2j[i][0],j1=i2j[i][1],j2=i2j[i][2];
      if(                j0>=0&&j0<N&&i!=j2i[j0][0]&&i!=j2i[j0][1]&&i!=j2i[j0][2])list_i[j2ind_i[j0]++]=i;
      if(        j1!=j0&&j1>=0&&j1<N&&i!=j2i[j1][0]&&i!=j2i[j1][1]&&i!=j2i[j1][2])list_i[j2ind_i[j1]++]=i;
      if(j2!=j1&&j2!=j0&&j2>=0&&j2<N&&i!=j2i[j2][0]&&i!=j2i[j2][1]&&i!=j2i[j2][2])list_i[j2ind_i[j2]++]=i;
    }
    tock(9);
    
    tick(10);
    for(j=0,tot=0;j<N;j++){
      j2ind_i[j]=tot;
      tot+=j2num_i[j];
    }

    mintree.resize(M*2+50);
    fill(mintree.begin(),mintree.end(),0);
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
    nbp0.resize(tot);
    nbp1.resize(tot);
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


    // Write aligned genome, out[]
    tick(11);
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
