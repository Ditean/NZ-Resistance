#!/bin/bash

if [[ FASTA == TRUE ]]
then
  # Create Variable with basename + ending removed
  stub=$(basename $FORWARD)
  r1=$(echo $stub | sed 's/.fastq.gz//')

  stub=$(basename $REVERSE)
  r2=$(echo $stub | sed 's/.fastq.gz//')

  # Pre-Trimming Statistics
  fastqc $FORWARD $REVERSE -t $THREADS -o ${STAMP}/${SAMPLE} --extract

  #Move data files
  mv $STAMP/${SAMPLE}/${r1}_fastqc/fastqc_data.txt $STAMP/${SAMPLE}/pre_trimmed_${r1}_fastqc.txt
  mv $STAMP/${SAMPLE}/${r2}_fastqc/fastqc_data.txt $STAMP/${SAMPLE}/pre_trimmed_${r2}_fastqc.txt
  #Move graphs
  mv $STAMP/${SAMPLE}/${r1}_fastqc/Images/per_base_quality.png $STAMP/${SAMPLE}_pre_R1_quality.png
  mv $STAMP/${SAMPLE}/${r2}_fastqc/Images/per_base_quality.png $STAMP/${SAMPLE}_pre_R2_quality.png

  #Remove DIR from Pre-Trimming
  rm -r $STUMP/${SAMPLE}/${r1}_fastqc
  rm -r $STUMP/${SAMPLE}/${r2}_fastqc

  # Post-Trimming
  fastqc ${STAMP}/${SAMPLE}/${SAMPLE}.R1.fastq.gz ${STAMP}/${SAMPLE}/${SAMPLE}.R2.fastq.gz -t $THREADS -o $STAMP --extract

  mv $STAMP/${SAMPLE}/${r1}_fastqc/fastqc_data.txt $STAMP/${SAMPLE}/post_trimmed_${r1}_fastqc.txt
  mv $STAMP/${SAMPLE}/${r2}_fastqc/fastqc_data.txt $STAMP/${SAMPLE}/post_trimmed_${r2}_fastqc.txt
  #Move graphs
  mv $STAMP/${SAMPLE}/${r1}_fastqc/Images/per_base_quality.png $STAMP/${SAMPLE}_post_R1_quality.png
  mv $STAMP/${SAMPLE}/${r2}_fastqc/Images/per_base_quality.png $STAMP/${SAMPLE}_post_R2_quality.png

  #Remove DIR from Post-Trimming
  rm -r $STUMP/${SAMPLE}/${r1}_fastqc
  rm -r $STUMP/${SAMPLE}/${r2}_fastqc

else # If the user didn't provide PATH to files

  fastqc ${SAMPLE}.R1.*fastq.gz ${REVERSE}.R2.*fastq.gz -t $THREADS -o $STAMP --extract

  #Move data files
  mv $STAMP/${SAMPLE}/${SAMPLE}*R1_fastqc/fastqc_data.txt $STAMP/${SAMPLE}/pre_trimmed_${SAMPLE}_R1_fastqc.txt
  mv $STAMP/${SAMPLE}/${SAMPLE}*R1_fastqc/fastqc_data.txt $STAMP/${SAMPLE}/pre_trimmed_${SAMPLE}_R2_fastqc.txt
  #Move graphs
  mv $STAMP/${SAMPLE}/${SAMPLE}_R1_fastqc/Images/per_base_quality.png $STAMP/${SAMPLE}_pre_R1_quality.png
  mv $STAMP/${SAMPLE}/${SAMPLE}_R1_fastqc/Images/per_base_quality.png $STAMP/${SAMPLE}_pre_R2_quality.png

  #Remove DIR from Pre-Trimming
  rm -r $STUMP/${SAMPLE}/${SAMPLE}_R1_fastqc
  rm -r $STUMP/${SAMPLE}/${SAMPLE}_R2_fastqc

  # Post-Trimming
  fastqc ${STAMP}/${SAMPLE}/${SAMPLE}.R1.fastq.gz ${STAMP}/${SAMPLE}/${SAMPLE}.R2.fastq.gz -t $THREADS -o $STAMP --extract

  mv $STAMP/${SAMPLE}/${SAMPLE}*R1_fastqc/fastqc_data.txt $STAMP/${SAMPLE}/post_trimmed_${SAMPLE}_R1_fastqc.txt
  mv $STAMP/${SAMPLE}/${SAMPLE}*R2_fastqc/fastqc_data.txt $STAMP/${SAMPLE}/post_trimmed_${SAMPLE}_R2_fastqc.txt
  #Move graphs
  mv $STAMP/${SAMPLE}/${SAMPLE}*R1_fastqc/Images/per_base_quality.png $STAMP/${SAMPLE}_post_R1_quality.png
  mv $STAMP/${SAMPLE}/${SAMPLE}*R1_fastqc/Images/per_base_quality.png $STAMP/${SAMPLE}_post_R2_quality.png

  #Remove DIR from Post-Trimming
  rm -r $STUMP/${SAMPLE}/${SAMPLE}_R1_fastqc
  rm -r $STUMP/${SAMPLE}/${SAMPLE}_R2_fastqc

#### STATISTICS

# Info on FastQC version and encoding of sequences
i1=$(cat $STUMP/${SAMPLE}/pre_trimmed_${SAMPLE}_R1_fastqc.txt | grep "##FastQC" | awk '{print substr($0, index($0,$2)) }')
i2=$(cat $STAMP/${SAMPLE}/pre_trimmed_${SAMPLE}_R1_fastqc.txt | grep "Encoding" | awk '{print substr($0, index($0,$2)) }')

# Info on Pre-Trimmed Reads
i3=$(cat $STAMP/${SAMPLE}/pre_trimmed_${SAMPLE}_R1_fastqc.txt | grep "Total Sequences" | awk '{print substr($0, index($0,$3)) }')
i4=$(cat $STAMP/${SAMPLE}/pre_trimmed_${SAMPLE}_R1_fastqc.txt | grep "%GC" | awk '{print substr($0, index($0,$2)) }')

# Info on Post-Trimmed Reads
i5=$(cat $STAMP/${SAMPLE}/post_trimmed_${SAMPLE}_R1_fastqc.txt | grep "Total Sequences" | awk '{print substr($0, index($0,$3)) }')
i6=$(cat $STAMP/${SAMPLE}/post_trimmed_${SAMPLE}_R1_fastqc.txt | grep "%GC" | awk '{print substr($0, index($0,$2)) }')

#NEED TO ADD IN BARCODE

while read line
do
  position=$(cat $line | awk '{print $1}')
  if [[ $barcode == $line ]]
  then
    " >> temp.tsv
  else
    echo $line >> temp.tsv
  fi
done < manifest.txt

mv tempt.tsv manifest.tsv






  i1=$(cat $STAMP/${r1}_fastqc/fastqc_data.txt | grep "##FastQC" | awk '{print substr($0, index($0,$2)) }')
  i3=$(cat $STAMP/${r1}_fastqc/fastqc_data.txt | grep "Total Sequences" | awk '{print substr($0, index($0,$3)) }')
  i4=$(cat $STAMP/${r1}_fastqc/fastqc_data.txt | grep "%GC" | awk '{print substr($0, index($0,$2)) }')
  i7=$(cat $STAMP/${r2}_fastqc/fastqc_data.txt | grep "Total Sequences" | awk '{print substr($0, index($0,$3)) }')
  i8=$(cat $STAMP/${r2}_fastqc/fastqc_data.txt | grep "%GC" | awk '{print substr($0, index($0,$2)) }')

else

  fastqc ${SAMPLE}.R1.*fastq.gz ${REVERSE}.R2.*fastq.gz -t $THREADS -o $STAMP --extract
  i1=$(cat $STAMP/${SAMPLE}*R1*_fastqc/fastqc_data.txt | grep "##FastQC" | awk '{print substr($0, index($0,$2)) }')
  i2=$(cat $STAMP/${SAMPLE}*R1*_fastqc/fastqc_data.txt | grep "Encoding" | awk '{print substr($0, index($0,$2)) }')
  i3=$(cat $STAMP/${SAMPLE}*R1*_fastqc/fastqc_data.txt | grep "Total Sequences" | awk '{print substr($0, index($0,$3)) }')
  i4=$(cat $STAMP/${SAMPLE}*R1*_fastqc/fastqc_data.txt | grep "%GC" | awk '{print substr($0, index($0,$2)) }')
  i5=$(cat $STAMP/${SAMPLE}*R2*_fastqc/fastqc_data.txt | grep "##FastQC" | awk '{print substr($0, index($0,$2)) }')
  i6=$(cat $STAMP/${SAMPLE}*R2*_fastqc/fastqc_data.txt | grep "Encoding" | awk '{print substr($0, index($0,$2)) }')
  i7=$(cat $STAMP/${SAMPLE}*R2*_fastqc/fastqc_data.txt | grep "Total Sequences" | awk '{print substr($0, index($0,$3)) }')
  i8=$(cat $STAMP/${SAMPLE}*R2*_fastqc/fastqc_data.txt | grep "%GC" | awk '{print substr($0, index($0,$2)) }')

# Pre-Trimming
fastqc ${SAMPLE}.R1.fastq.gz ${SAMPLE}.R2.fastq.gz -t $THREADS --extract -o $STAMP

i1=$(cat $STAMP/${SAMPLE}._fastqc/fastqc_data.txt | grep "##FastQC" | awk '{print substr($0, index($0,$2)) }')
i2=$(cat $STAMP/${SAMPLE}_fastqc/fastqc_data.txt | grep "Encoding" | awk '{print substr($0, index($0,$2)) }')
i3=$(cat $STAMP/${SAMPLE}_fastqc/fastqc_data.txt | grep "Total Sequences" | awk '{print substr($0, index($0,$3)) }')
i4=$(cat $STAMP/${SAMPLE}_fastqc/fastqc_data.txt | grep "%GC" | awk '{print substr($0, index($0,$2)) }')

mv ${SAMPLE}_fastqc
#barcode=$(cat $ROOT/dataset/manifest.tsv | grep $SAMPLE | grep $run_stamp | cut -f 1)

#message=$(echo -e "\t$i1\t$i2\t$i3\t$i4")

#sed '/$BARCODE/$/$message/' file
#awk '{if(/^$BARCODE /) {$0=$0"$message"}}'

# Post-Trimming
fastqc ${SAMPLE}.R1.fastq.gz ${SAMPLE}.R2.fastq.gz -t $THREADS --extract -o $STAMP

i5=$(cat $STAMP/${SAMPLE}_fastqc/fastqc_data.txt | grep "##FastQC" | awk '{print substr($0, index($0,$2)) }')
i6=$(cat $STAMP/${SAMPLE}_fastqc/fastqc_data.txt | grep "Encoding" | awk '{print substr($0, index($0,$2)) }')
i7=$(cat $STAMP/${SAMPLE}_fastqc/fastqc_data.txt | grep "Total Sequences" | awk '{print substr($0, index($0,$3)) }')
i8=$(cat $STAMP/${SAMPLE}_fastqc/fastqc_data.txt | grep "%GC" | awk '{print substr($0, index($0,$2)) }')

statistics=$(echo -e "$1\t$2\t$3\t$4\t$5\t$6\t$7$\t$8")

while read line
do
  position=$(cat $line | awk '{print $1}')
  if [[ $barcode == $line ]]
  then
    echo -e "${line}\t${statistics}" >> temp.tsv
  else
    echo $line >> temp.tsv
  fi
done < manifest.txt

mv tempt.tsv manifest.tsv




#lineno=$(cat data/manifest.tsv | grep -n $barcode | cut -d':' -f 1)

#sed '${lineno}s/$/$message/' file
