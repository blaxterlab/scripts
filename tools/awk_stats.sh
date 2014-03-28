#!/bin/sh
sort -n |
awk 'BEGIN{c=0;sum=0;}\
/^[^#]/{if($1 > 0){a[c++]=$1;sum+=$1;sd+=$1*$1;}}\
END{ave=sum/c;\
stdev=sqrt((sd-ave*ave*c)/(c-1));\
if((c%2)==1){median=a[int(c/2)];}\
else{median=(a[c/2]+a[c/2-1])/2;}\
print "count\t",c,"\nmean\t",ave,"\nmedian\t",median,"\nmin\t",a[0],"\nmax\t",a[c-1],"\nsd\t",stdev}'
