
# Produce resistance report from VCF file


first_line=(Isoniazid Rifampicin Ethambutol Pyrazinamide)
second_line=(Amikacin Capreomycin Ethionamide Fluoroquinolones Kanamycin Linezolid Para-aminosalicylic Streptomycin)

# Assign barcode
barcode=$(cat $ROOT/dataset/manifest.tsv | grep $SAMPLE | grep $run_stamp | cut -f 1)


for i in ${first_line[@]}
do
  while read line
  do
    snp=$(echo $line | awk '{print $1,$2,$3}' OFS='\t')
    gene=$(echo $line | awk 'BEGIN{FS=OFS=" "}{print $4,$5}')
    pmid=$(echo $line | awk '{print $6}' OFS='\t')
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
    echo -e "${barcode}\tFirst Line\t${i}\tSensitive\tNo resistance predicted\tNO\tNO\tNO\tNO" >> $ROOT/dataset/resistance_profile.tsv
  fi
done

for i in ${second_line[@]}
do
  while read line
  do
    snp=$(echo $line | awk '{print $1,$2,$3}' OFS='\t')
    gene=$(echo $line | awk 'BEGIN{FS=OFS=" "}{print $4,$5}')
    pmid=$(echo $line | awk '{print $6}' OFS='\t')
#    pmid=$(echo $line | awk 'NF && NF-1 {print ($(NF-1))}')
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
    echo -e "${barcode}\tSecond Line\t${i}\tSensitive\tNo resistance predicted\tNO\tNO\tNO\tNO" >> $ROOT/dataset/resistance_profile.tsv
  fi
done


##### Create a second file that contains all mutations (including low confidence)

for i in ${first_line[@]}
do
  while read line
  do
    snp=$(echo $line | awk '{print $1,$2,$3}' OFS='\t')
    gene=$(echo $line | awk 'BEGIN{FS=OFS=" "}{print $4,$5}')
    pmid=$(echo $line | awk '{print $6}' OFS='\t')
    confidence=$(echo $line | awk '{print $7}' OFS='\t')
    drug=$(echo $line | awk '{print $8}' OFS='\t')
    # If the mutation is found
    if [[ $i == $drug ]]
    then
      if grep -q $snp "$ROOT/output/${SAMPLE}.vcf" # NEED TO LOOP TO CORRECT FILE
      then
        echo -e "${barcode}\tFirst Line\t${i}\tResistance\t${gene}\t${snp}\t${pmid}\t${confidence}" >> $ROOT/dataset/genomic_profile.tsv
      fi
    fi
  done < $ROOT/reference/mutations/Mtb_snp_db.txt
done

for i in ${second_line[@]}
do
  while read line
  do
    snp=$(echo $line | awk '{print $1,$2,$3}' OFS='\t')
    gene=$(echo $line | awk 'BEGIN{FS=OFS=" "}{print $4,$5}')
    pmid=$(echo $line | awk '{print $6}' OFS='\t')
    confidence=$(echo $line | awk '{print $7}' OFS='\t')
    drug=$(echo $line | awk '{print $8}' OFS='\t')
    # If the mutation is found
    if [[ $i == $drug ]]
    then
      if grep -q $snp "$ROOT/output/${SAMPLE}.vcf" # NEED TO LOOP TO CORRECT FILE
      then
        echo -e "${barcode}\tSecond Line\t${i}\tResistance\t${gene}\t${snp}\t${pmid}\t${confidence}" >> $ROOT/dataset/genomic_profile.tsv
      fi
    fi
  done < $ROOT/reference/mutations/Mtb_snp_db.txt
done
