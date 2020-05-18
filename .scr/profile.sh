
# Produce resistance report from VCF file


first_line=(isoniazid rifampicin ethambutol pyrazinamide)
second_line=(amikacin capreomycin ethionamide fluoroquinolones kanamycin linezolid para-aminosalicylic streptomycin)

# Assign barcode
barcode=$(cat $ROOT/dataset/manifest.tsv | grep $SAMPLE | grep $run_stamp | cut -f 1)


for i in ${first_line[@]}
do
  while read line
  do
    snp=$(echo $line | awk '{print $1,$2,$3}' OFS='\t')
    gene=$(echo $line | awk 'BEGIN{FS=OFS=" "}{print $4,$5}')
    pmid=$(echo $line | awk 'NF && NF-1 {print ($(NF-1))}')
    # If the mutation is found
    mutation=
    if grep -q $snp "$ROOT/output/${SAMPLE}.vcf" # NEED TO LOOP TO CORRECT FILE
    then
      echo -e "${barcode}\tFirst Line\t${i}\tResistance\t${gene}\t${snp}\t${pmid}" >> $ROOT/dataset/resistance_profile.tsv
      mutation=FOUND
      break
    fi
  done < $ROOT/reference/mutations/${i}.txt

  # If no mutation is found
  if [[ ! $mutation == "FOUND" ]]
  then
    echo -e "${barcode}\tFirst Line\t${i}\tSensitive\tNo resistance predicted\tNO\tNO\tNO\tNO\tNO" >> $ROOT/dataset/resistance_profile.tsv
  fi
done

for i in ${second_line[@]}
do
  while read line
  do
    snp=$(echo $line | awk '{print $1,$2,$3}' OFS='\t')
    gene=$(echo $line | awk 'BEGIN{FS=OFS=" "}{print $4,$5}')
    pmid=$(echo $line | awk 'NF && NF-1 {print ($(NF-1))}')
    # If the mutation is found
    mutation=
    if grep -q $snp "$ROOT/output/${SAMPLE}.vcf" # NEED TO LOOP TO CORRECT FILE
    then
      echo -e "${barcode}\tSecond Line\t${i}\tResistance\t${gene}\t${snp}\t${pmid}" >> $ROOT/dataset/resistance_profile.tsv
      mutation=FOUND
      break
    fi
  done < $ROOT/reference/mutations/${i}.txt

  # If no mutation is found
  if [[ ! $mutation == "FOUND" ]]
  then
    echo -e "${barcode}\tSecond Line\t${i}\tSensitive\tNo resistance predicted\tNO\tNo\tNO\tNO\tNO" >> $ROOT/dataset/resistance_profile.tsv
  fi
done
