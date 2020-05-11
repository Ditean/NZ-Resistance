
# First Line Drugs

for i in {resistance[@]}
do
  cat $i | while read line
  do
    compare=(cat $line | cut -f 1-3)
    if grep -q $compare $vcf
    then
      echo -e "$barcode "





    #cat $vcf | grep "NC_000962.3  $compare"
