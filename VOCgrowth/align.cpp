// Aligns SARS-CoV-2 fasta files to reference genome Wuhan-Hu-1

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

typedef unsigned char UC;
typedef unsigned int UI;

// R = number of bases from which the indexes are formed
#define R 9
vector<int> refdict[1<<R*2];

// Maximum number of bases
#define MAXGS 40000 // 40000 // alter

// Count threshold for offsets
#define MINOFFSETCOUNT 20

const int undefined=0x7f7f7f7f;
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
  //  return gisaidname;// alter
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
  
  int i,j,t;
  
  string refgenome;
  int N;
  {
    const char*reffn="refgenome";
    std::ifstream fp(reffn);
    if(fp.fail())error(1,errno,"\nCouldn't open %s",reffn);
    while(std::getline(fp,refgenome))if(refgenome.size()>0&&refgenome[0]!='>')break;      
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
  int badi=-1;
  for(i=0,t=0;i<N;i++){
    int b=base2num[refgenome[i]&255];
    if(b<0){badi=i;continue;}
    t=t>>2|b<<(R-1)*2;
    if(i>=badi+R){assert(t>=0&&t<(1<<R*2));refdict[t].push_back(i-(R-1));}
  }

  int jumppen[MAXGS+1];
  for(t=0;t<=MAXGS;t++)jumppen[t]=int(floor(sqrt(t)+1e-6));

  vector<int> j_to_num_i(N), j_to_ind_i(N), list_i;

  int linenum=0,nwrite=0;
  bool last;
  string header,line;
  UC genome[MAXGS];
  last=!std::getline(std::cin,header);linenum++;
  int skip[3]={0};
  while(!last){// Main loop
    tick(1);
    string id=parseheader(header);
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
      if(M+s>MAXGS){fprintf(stderr,"Warning: Overlong genome at line %d\n",linenum);continue;}
      memcpy(genome+M,&line[0],s);M+=s;
    }
    tock(2);
    
    // Offset = (index in ref genome) - (index in current genome) = j-i
    // offsetcount[MAXGS+offset] = number of possible uses of this offset (assuming any R-tuple in genome can match same R-tuple anywhere in reference genome)
    // indexkey[i] = R-tuple of bases at position i in the current genome
    int offsetcount[MAXGS*2]={0};
    int indexkey[MAXGS]={0};
    int badi=-1;
    tick(3);
    for(i=0,t=0;i<M;i++){
      int b=base2num[genome[i]];
      if(b<0){badi=i;continue;}
      int j=i-(R-1);
      t=t>>2|b<<(R-1)*2;
      if(j>badi){
        assert(t>=0&&t<(1<<R*2));
        indexkey[j]=t;
        for(int x:refdict[t])offsetcount[MAXGS+x-j]++;
      }else if(i>=R-1)indexkey[j]=undefined;
    }
    //for(int o=0;o<MAXGS*2;o++)if(offsetcount[o]>=MINOFFSETCOUNT)fprintf(stderr,"Offset %d %d\n",o-MAXGS,offsetcount[o]);
    tock(3);

    tick(4);
    int pointoffset_i[MAXGS],pointoffset[MAXGS],best[MAXGS]={0};
    for(i=0;i<=M-R;i++){
      t=indexkey[i];
      pointoffset_i[i]=undefined;
      if(t!=undefined){
        int best_i=MINOFFSETCOUNT-1;
        for(int j:refdict[t]){
          int c=offsetcount[MAXGS+j-i];
          if(c>best[j]){best[j]=c;pointoffset[j]=j-i;}
          if(c>best_i){best_i=c;pointoffset_i[i]=j-i;}
        }
      }
    }
    tock(4);
    tick(5);
    // Build up to 4 possible offsets, i_to_j[i][0,1], j_to_i[j][0,1], to use at each position (i in the current genome, j in ref)
    // i_to_j[][] works well for deletions
    // j_to_i[][] works well for insertions
    int i_to_j[MAXGS][2],j_to_i[MAXGS][2];
    // Approach from right
    int nearest=undefined;
    for(i=M-1;i>M-R;i--)i_to_j[i][1]=undefined;
    for(i=M-R;i>=0;i--){
      if(pointoffset_i[i]!=undefined)nearest=pointoffset_i[i];
      if(nearest!=undefined)i_to_j[i][1]=i+nearest; else i_to_j[i][1]=undefined;
    }
    nearest=undefined;
    for(j=N-1;j>N-R;j--)j_to_i[j][1]=undefined;
    for(j=N-R;j>=0;j--){
      if(best[j]>=MINOFFSETCOUNT)nearest=pointoffset[j];
      if(nearest!=undefined)j_to_i[j][1]=j-nearest; else j_to_i[j][1]=undefined;
    }
    // Approach from left
    nearest=undefined;
    for(i=0;i<R-1;i++)i_to_j[i][0]=undefined;
    for(i=0;i<=M-R;i++){
      if(pointoffset_i[i]!=undefined)nearest=pointoffset_i[i];
      if(nearest!=undefined)i_to_j[i+R-1][0]=i+R-1+nearest; else i_to_j[i+R-1][0]=undefined;
    }
    nearest=undefined;
    for(j=0;j<R-1;j++)j_to_i[j][0]=undefined;
    for(j=0;j<=N-R;j++){
      if(best[j]>=MINOFFSETCOUNT)nearest=pointoffset[j];
      if(nearest!=undefined)j_to_i[j+R-1][0]=j+R-1-nearest; else j_to_i[j+R-1][0]=undefined;
    }
    tock(5);

    /*
    for(j=27500;j<N;j++){
      fprintf(stderr,"%6d |",j);
      if(best[j]>=MINOFFSETCOUNT)fprintf(stderr," %10d %6d |",pointoffset[j],best[j]); else fprintf(stderr," ---------- %6d |",best[j]);
      int y;
      for(y=0;y<2;y++){
        i=j_to_i[j][y];
        fprintf(stderr," %10d %s",i,i>=0&&i<M&&genome[i]==refgenome[j]?"*":".");
      }
      fprintf(stderr,"\n");
    }
    exit(0);
    */
    
    /*
    for(j=0;j<N;j++){
      fprintf(stderr,"%6d |",j);
      int y;
      for(y=0;y<2;y++){
        i=j_to_i[j][y];
        //fprintf(stderr," *%d",off);
        if(i!=undefined){
          fprintf(stderr," %6d %6d |",i,j);
          if(i>=0&&i<M){
            int x;
            for(x=0;x<2;x++){
              int j_i=i_to_j[i][x];
              fprintf(stderr," %6d %6d |",i,j_i);
            }
          }
        }
      }
      fprintf(stderr,"\n");
    }
    exit(0);
    */
    /*
    for(i=0;i<M;i++){
      int y;
      for(y=0;y<2;y++){
        j=i_to_j[i][y];
        if(j!=undefined&&(y==0||j!=i_to_j[i][y-1])){
          if(j>=0&&j<N){
            int uniq=1;
            int x;
            for(x=0;x<2;x++){
              int i_j=j_to_i[j][x];
              if(i_j!=undefined){
                if(i==i_j)uniq=0;
              }
            }
            if(uniq)fprintf(stderr,"UUU %6d %6d\n",i,j);
          }
        }
      }
    }
    exit(0);
    */

    /*
    int count[MAXGS]={0};
    for(j=0;j<N;j++){
      int y;
      for(y=0;y<2;y++){
        i=j_to_i[j][y];
        if(i!=undefined&&(y==0||i!=j_to_i[j][y-1])){
          if(i>=0&&i<M){
            int uniq=1;
            int x;
            for(x=0;x<2;x++){
              int j_i=i_to_j[i][x];
              if(j_i!=undefined){
                if(j==j_i)uniq=0;
              }
            }
            if(uniq){
              count[i]++;
              fprintf(stderr,"UUU %6d %6d\n",i,j);
            }
          }
        }
      }
    }
    for(i=0;i<M;i++)if(count[i]!=0)fprintf(stderr,"CCC %6d %6d\n",i,count[i]);
    exit(0);
    */

    tick(9);
    memset(&j_to_num_i[0],0,N*sizeof(int));
    for(i=0;i<M;i++){
      j=i_to_j[i][0];
      if(j>=0&&j<N){
        if(i!=j_to_i[j][0]&&i!=j_to_i[j][1])j_to_num_i[j]++;
        int j1=i_to_j[i][1];
        if(j1!=j&&j1>=0&&j1<N&&i!=j_to_i[j1][0]&&i!=j_to_i[j1][1])j_to_num_i[j1]++;
      }
    }
    int tot=0;
    for(j=0;j<N;j++){
      j_to_ind_i[j]=tot;
      tot+=j_to_num_i[j];
    }
    assert(tot<=2*M);
    list_i.resize(tot);
    for(i=0;i<M;i++){
      j=i_to_j[i][0];
      if(j>=0&&j<N){
        if(i!=j_to_i[j][0]&&i!=j_to_i[j][1])list_i[j_to_ind_i[j]++]=i;
        int j1=i_to_j[i][1];
        if(j1!=j&&j1>=0&&j1<N&&i!=j_to_i[j1][0]&&i!=j_to_i[j1][1])list_i[j_to_ind_i[j1]++]=i;
      }
    }
    fprintf(stderr,"Total %6d\n",tot);
    tock(9);
    
    tick(6);
    // Dyn prog on the two allowable offsets: j_to_i[i][]
    int bp[MAXGS][2]={0};// Back pointers; bp[i][j] is defined if value is finite
    int st[2]={0,0};// State: st[j] = score (lower is better) given ended with offset j_to_i[i-1][j]
    for(j=0;j<N;j++){
      // Transition x -> y,  x=prev offset index, y=current offset index
      int y,newst[2]={infinity,infinity};
      for(y=0;y<2;y++){
        int x;
        i=j_to_i[j][y];
        if(i!=undefined){
          int best=infinity,bestx=0;
          for(x=0;x<2;x++){
            int i0=-1;
            int v=0;
            if(j>0){// Initial jump is free
              i0=j_to_i[j-1][x];
              if(i0==undefined)continue;
              v=jumppen[abs(i-1-i0)];
            }
            v+=st[x]+(i0>=i||i<0||i>=M||refgenome[j]!=genome[i]);
            if(v<=best){best=v;bestx=x;}
          }
          newst[y]=best;
          bp[j][y]=bestx;
        }
      }
      st[0]=newst[0];
      st[1]=newst[1];
      //fprintf(stderr,"%6d | %10d %10d %d %10d | %10d %10d %d %10d\n",j,j_to_i[j][0],j-j_to_i[j][0],bp[j][0],st[0],j_to_i[j][1],j-j_to_i[j][1],bp[j][1],st[1]);
    }
    tock(6);

    tick(7);
    // Write aligned genome, out[]
    UC out[MAXGS+1];
    memset(out,'-',N);
    out[N]=0;
    int s=(st[1]<st[0]);
    //fprintf(stderr,"Score %d\n",st[s]);
    int pri=M,prj=N-1;
    for(j=N-1;j>=0;j--){
      i=j_to_i[j][s];
      //fprintf(stderr,"%6d (%10d %10d) %10d  %6d\n",j,j_to_i[j][0],j_to_i[j][1],o,i);
      assert(i!=undefined);
      //fprintf(stderr,"%6d %6d %6d %s %s\n",j,i,j-i,i>=0&&i<M&&refgenome[j]==genome[i]?"*":".",i!=undefined&&i>=0&&i<pri?"U":".");
      if(i!=undefined&&i>=0&&i<pri){
        out[j]=genome[i];pri=i;
      }
      s=bp[j][s];
      if(0)if(j==0||(j_to_i[j-1][s]!=undefined&&j_to_i[j-1][s]+1!=i)){
        fprintf(stderr,"%6d - %6d   %10d -->  %6d - %6d\n",j,prj,j-i,i,prj-(j-i));
        prj=j-1;
      }
    }
    tock(7);
    
    tick(8);
    FILE*fp;
    if(datadir=="")fp=stdout; else fp=fopen((datadir+"/"+date).c_str(),"a");
    header=line;// Next header is the last-read line
    switch(compression){
    case 0:
      fprintf(fp,">%s",id.c_str());
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
