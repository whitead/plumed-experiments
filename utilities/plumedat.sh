#!/bin/bash
# bash script
if(($#==0))
then
  cat << END_OF_HELP
###############
syntax:

  plumedat.sh f1 f2 ... < file
or
  plumedat.sh -l < file

file is the name of the COLVAR file
f1, f2, ... are the names of the required fields
if a required field is not available in the COLVAR file, "NA" is written in the output
with -l, the available choices are listed

example:

plumedat.sh time temp < COLVAR

prints a two-column file, with the time in the first column and the temperature
in the second column
###############
END_OF_HELP
exit 0
fi

# put the list of requested fields in this temporary file
TMPFILE=.tmp$$
test -e $TMPFILE && rm $TMPFILE
for name
do
  echo $name >> $TMPFILE
done

awk 'BEGIN{
# read the list of requested fields
  print_list=0;
  for(i=0;;i++){
    err=getline asked[i] < "'$TMPFILE'";
    if(asked[i]=="-l") print_list=1;
    if(err<=0) break;
  }
  nasked=i;
}{
  if(NF==0) next;
  if($1=="#!" && $2=="FIELDS") {
    irun++;
    if(print_list==1) print "List of available keywords (run",irun,")"
    for(i=0;i<nasked;i++) iasked[i]=0;
    for(j=3;j<=NF;j++){
      if(print_list==1) {
        printf("   %s",$j);
        if((j-1)%10==0 || j==NF) printf("\n");
      } else {
        for(i=0;i<nasked;i++) {
          if($j==asked[i]) {
            iasked[i]=j-2;
          }
        }
      }
    }
   next
  }
  if(print_list==1) next
  if(substr($0,0,1)=="#") next
  for(i=0;i<nasked;i++){
    j=iasked[i];
    if(j>0) printf("%s ",$j)
    else printf("NA")
  }
  printf("\n");
}'

rm $TMPFILE
