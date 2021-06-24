wget https://cog-uk.s3.climb.ac.uk/phylogenetics/latest/cog_metadata.csv

tail -n +2 cog_metadata.csv|awk -F, '{if($7=="B.1.1.7")a[$5][0]+=1; else if($7=="B.1.617.2")a[$5][1]+=1;}END{for(x in a)if(a[x][0]>=10 && a[x][1]>=10)printf("%s  %6.3f   %6d  %6d\n",x,log(a[x][1]/a[x][0]),a[x][0],a[x][1]);}'|sort|head -n -4 > alphadelta
