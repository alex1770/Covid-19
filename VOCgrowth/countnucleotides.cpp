#include <stdio.h>
#include <string.h>
#include <assert.h>

#define REFLEN 29903

int counts[REFLEN][256];

int main(){
  int i,j,t;
  FILE *fp;
  unsigned char l[REFLEN+2];
  
  fp=fopen("fastalist","w");
  memset(counts,0,sizeof(counts));
  while(fgets((char*)l,REFLEN+2,stdin)){
    assert(l[0]=='>');
    fprintf(fp,"%s",l+1);
    assert(fgets((char*)l,REFLEN+2,stdin));
    for(i=0;i<REFLEN&&l[i];i++)counts[i][l[i]]++;
  }
  for(i=0;i<REFLEN;i++){
    printf("%5d",i);
    for(j=32;j<127;j++)if(counts[i][j])printf(" %c%d",j,counts[i][j]);
    for(j=t=0;j<256;j++)if(j<32||j>=127)t+=counts[i][j];
    if(t)printf(" ?%d",t);
    printf("\n");
  }
}
