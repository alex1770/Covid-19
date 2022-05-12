
// Aligns SARS-CoV-2 fasta files to reference genome Wuhan-Hu-1
// Todo: handle files with \r\n lines

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <getopt.h>
#include <error.h>
#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <unordered_map>
#include <unordered_set>
using std::string;
using std::vector;
using std::array;
using std::min;
using std::max;
using std::pair;
using std::unordered_map;
using std::unordered_set;

typedef unsigned char UC;
template<class T> struct array2d {
  size_t rows,cols;
  vector<T> data;
  array2d(){}
  array2d(size_t r, size_t c):rows(r),cols(c),data(r*c){}
  void resize(size_t r, size_t c){rows=r;cols=c;data.resize(r*c);}
  T* operator[](size_t index){return &data[index*cols];}// First level indexing
};

// R = number of bases from which the indexes are formed
#define R 9
vector<int> refdict[1<<R*2];

// Maximum number of bases
#define MAXGS 50000

// Count threshold for offsets
#define MINOFFSETCOUNT 20

const int undefined=0x7fffffff;
const int infinity=1000000000;

int main(){
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
  }

  int base2num[256];
  for(i=0;i<256;i++)base2num[i]=-1;
  base2num['A']=0;
  base2num['C']=1;
  base2num['G']=2;
  base2num['T']=3;
  for(i=0,t=0;i<N;i++){
    t=t>>2|base2num[refgenome[i]&255]<<(R-1)*2;
    if(i>=R-1&&t>=0&&t<(1<<R*2))refdict[t].push_back(i-(R-1));
  }

  int jumppen[MAXGS+1];
  for(t=0;t<=MAXGS;t++)jumppen[t]=int(floor(sqrt(t)+1e-6));

  int linenum=0;
  bool last;
  string name;
  UC genome[MAXGS];
  last=!std::getline(std::cin,name);linenum++;
  while(!last){
    printf("%s\n",name.c_str());
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
    unordered_map<int,int> offsetcount;
    vector<int> indexkey(M-(R-1));
    int badi=-1;
    for(i=0,t=0;i<M;i++){
      int b=base2num[genome[i]];
      if(b<0){badi=i;continue;}
      t=t>>2|b<<(R-1)*2;
      if(i>=badi+R){
        assert(t>=0&&t<(1<<R*2));
        indexkey[i-(R-1)]=t;
        for(int x:refdict[t])offsetcount[x-(i-(R-1))]+=1;
      }else if(i>=R-1)indexkey[i-(R-1)]=undefined;
    }
    unordered_set<int> okoffsets;
    for(auto &o:offsetcount)if(o.second>=MINOFFSETCOUNT)okoffsets.insert(o.first);
    if(okoffsets.size()==0){fprintf(stderr,"Warning: Can't find offsets for genome at lines preceding %d\n",linenum);break;}
    //for(int o:okoffsets)printf("Offset %d : %d\n",o,offsetcount[o]);

    array2d<int> offsets(M,2);

    // Approach from right
    int nearest=undefined;
    for(i=M-1;i>M-R;i--)offsets[i][1]=undefined;
    for(i=M-R;i>=0;i--){
      t=indexkey[i];
      if(t!=undefined)for(int x:refdict[t])if(okoffsets.count(x-i))nearest=x-i;
      offsets[i][1]=nearest;
    }

    // Approach from left
    nearest=undefined;
    for(i=0;i<R-1;i++)offsets[i][0]=undefined;
    for(i=R-1;i<M;i++){
      t=indexkey[i-(R-1)];
      if(t!=undefined)for(int x:refdict[t])if(okoffsets.count(x-(i-(R-1))))nearest=x-(i-(R-1));
      offsets[i][0]=nearest;
    }

    // Dyn prog on two allowable offsets: offsets[i][]
    array2d<int> bp(M,2);// Back pointers; bp[i][j] is defined if value is finite
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
            v+=st[j]+(genome[i]!=refgenome[i+off]);
            if(v<best){best=v;bestj=j;}
          }
          newst[k]=best;
          bp[i][k]=bestj;
        }
      }
      st[0]=newst[0];
      st[1]=newst[1];
    }

    vector<UC> out(N+1);
    memset(&out[0],'-',N);
    int s=(st[1]<st[0]);
    for(i=M-1;i>=0;i--){
      int o=offsets[i][s];
      if(o!=undefined){assert(i+o>=0&&i+o<N);out[i+o]=genome[i];}
      s=bp[i][s];
    }
    printf("%s\n",&out[0]);
  }
  
}
