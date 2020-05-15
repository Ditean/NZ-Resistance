
# Produce resistance report from VCF file

touch resistance_profile_$i.tsv

echo -e "BARCODE\tCLASS\tDRUG\tSTATUS\tDETAILS" > resistance_profile.tsv

# Name of first line drugs
first_line=(rifampicin isoniazid)

for i in ${first_line[@]}
do
  drug_array=()
  while read line
  do
    drug_array+=($(echo $line | cut -f 1))
  done < $ROOT/reference/${i}.txt

  for snp in ${drug_array[@]}
  do
    resistance=$(cat file.vcf | grep -e "NC_000962.3\t"${snp})

done

#==========

for i in ${first_line[@]}
do
  drug_array=()
  while read line
  do
    ref=$(echo $line | awk '{print $1,$2,$3}' OFS='\t')



# Add whole line into an array -> then you can just call on each part as $1 $2 $3

while read line
do
  resistance=$(echo $line | cut -f 1-3)
  while read vcf
  do
    snp=$(echo $vcf | cut -f 2-4)
    if [[ $resistance == $snp ]]
    then
      echo -e "$BARCODE\tFIRST CLASS\t${i}\tRESISTANCE\t" # Need to find way to parse
    elif



cat vcf | grep $line | cut -f 2-4



# ======================

# Drugs by class

first_line=(isoniazid rifamplicin ethambutol pyrazinimide)
second_line=(amikacin capreomycin ethionamide fluoroquinolones kanamycin linezolid para-aminosalicylic streptomycin)

touch resistance_profile.tsv

echo -e "BARCODE\tCLASS\tDRUG\tSTATUS\tDETAILS" > resistance_profile.tsv

for i in ${first_line[@]}
do
  while read line
  do
    snp=$(echo $line | awk '{print $1,$2,$3}' OFS='\t')
    if grep -q $snp "$ROOT/vcf/file"
    then
      echo -e "$FILE\tFirst Line\tDrug Resistance Predicted\t${snp}" >> resistance_profile.tsv
    fi
  done < high/${i}.txt
done

### CONSIDERATION

while read line
do
  snp=$(echo $line | awk '{print $1,$2,$3}' OFS='\t')
  if grep -q $snp "$ROOT/vcf/file"
  then
    awk 'NR==1,/input/{sub(/CHECK/,"YES")}1'
    break
  fi
done # Have awk change a value from N/A to Resistance or No resistance predicted 


awk 'NR==1,/input/{sub(/CHECK/,"YES")}1' log_setup.txt > temp.txt && mv temp.txt log_setup.txt
