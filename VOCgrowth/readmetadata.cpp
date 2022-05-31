// Reads GISAID and COG-UK metadata files and synthesises a reduced metadata output, in
// preparation for align.cpp reading fasta files and for subsequent variant detection.

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
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
#include <unordered_map>
using std::string;
using std::cerr;
using std::vector;
using std::unordered_map;

typedef unsigned char UC;
typedef unsigned int UI;

struct metadata {
  string date,id,lineage,location;
  bool operator<(const metadata&other){
    if(date<other.date)return true;
    if(date==other.date&&id<other.id)return true;
    return false;
  }
};

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

bool okdate(string date){
  return date.size()==10&&date[0]=='2'&&date[1]=='0'&&date[4]=='-'&&date[7]=='-';
}

// Not currently used
string getid(string gisaidname){
  vector<string> ida=split(gisaidname,"/");
  int n=ida.size();
  if(n>=3)return ida[n-3]+"/"+ida[n-2]+"/"+ida[n-1];
  return "";
}

string processlocation(string &id,string &loc){
  if(loc=="UK"){
    std::size_t f=id.find('/');
    assert(f!=std::string::npos);
    return "Europe / United Kingdom / "+id.substr(0,f);
  }
  return loc;
}

unordered_map<string,metadata> csv2map(string fn,string sep,string idprefix,string id,string date,string lineage,string location){
  unordered_map<string,metadata> ret;
  std::ifstream fp(fn);
  if(fp.fail())error(1,errno,"Couldn't open %s",fn.c_str());
  string header;
  std::getline(fp,header);
  vector<string> headers=split(header,sep);
  int idcol=-1,datecol=-1,lincol=-1,loccol=-1;
  for(UI i=0;i<headers.size();i++){
    if(headers[i]==id)idcol=i;
    if(headers[i]==date)datecol=i;
    if(headers[i]==lineage)lincol=i;
    if(headers[i]==location)loccol=i;
  }
  assert(idcol>=0&&datecol>=0&&lincol>=0&&loccol>=0);
  string line;
  while(std::getline(fp,line)){
    vector<string> data=split(line,sep);
    string date=data[datecol];
    if(okdate(date)){
      string id=idprefix+data[idcol];
      string loc=processlocation(id,data[loccol]);
      ret[id]={date,id,data[lincol],loc};
    }
  }
  fp.close();
  return ret;
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
    cerr << "Usage: readmetadata [options]\n";
    cerr << "       -x<string> Data directory\n";
    exit(1);
  }

  unordered_map<string,metadata> id2meta,cog2meta;
  
  if(datadir!="")mkdir(datadir.c_str(),0777);
  
  cerr << "Loading metadata\n";
  
  id2meta=csv2map("metadata.tsv",     "\t","",       "Virus name",   "Collection date","Pango lineage","Location");
  cerr << "Read " << id2meta.size() << " entries from GISAID metadata\n";
  
  cog2meta=csv2map("cog_metadata.csv",",","hCoV-19/","sequence_name","sample_date",    "lineage",      "country");
  cerr << "Read " << cog2meta.size() << " entries from COG-UK metadata\n";
  
  for(auto &m:cog2meta)id2meta[m.first]=m.second;

  vector<metadata> vm;
  for(const auto &m:id2meta)vm.push_back(m.second);
  std::sort(vm.begin(),vm.end());

  std::ofstream rfp;
  if(datadir!="")rfp.open(datadir+"/metadata.tsv");
  std::ostream &fp=(datadir!="")?rfp:std::cout;
  fp<<"date\tID\tlineage\tlocation\n";
  for(metadata &m:vm)fp<<m.date<<"\t"<<m.id<<"\t"<<m.lineage<<"\t"<<m.location<<"\n";
}
