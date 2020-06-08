#!/bin/bash

if [[ FASTA == TRUE ]]
then
  stub=$(basename $FORWARD)
  r1=$(echo $stub | sed 's/.fastq.gz//')

  stub=$(basename $REVERSE)
  r2=$(echo $stub | sed 's/.fastq.gz//')

  fastqc $FORWARD $REVERSE -t $THREADS -o $STAMP --extract
  mv $STAMP/${r1}_fastqc/fastqc_data.txt $STAMP/pre_trimmed_${r1}_fastqc.txt
  mv $STAMP/${r1}_fastqc/fastqc_data.txt $STAMP/pre_trimmed_${r1}_fastqc.txt
else
  
  fastqc ${SAMPLE}.R1.*fastq.gz ${REVERSE}.R2.*fastq.gz -t $THREADS -o $STAMP --extract








  i1=$(cat $STAMP/${r1}_fastqc/fastqc_data.txt | grep "##FastQC" | awk '{print substr($0, index($0,$2)) }')
  i2=$(cat $STAMP/${r1}_fastqc/fastqc_data.txt | grep "Encoding" | awk '{print substr($0, index($0,$2)) }')
  i3=$(cat $STAMP/${r1}_fastqc/fastqc_data.txt | grep "Total Sequences" | awk '{print substr($0, index($0,$3)) }')
  i4=$(cat $STAMP/${r1}_fastqc/fastqc_data.txt | grep "%GC" | awk '{print substr($0, index($0,$2)) }')
  i5=$(cat $STAMP/${r2}_fastqc/fastqc_data.txt | grep "##FastQC" | awk '{print substr($0, index($0,$2)) }')
  i6=$(cat $STAMP/${r2}_fastqc/fastqc_data.txt | grep "Encoding" | awk '{print substr($0, index($0,$2)) }')
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



#lineno=$(cat data/manifest.tsv | grep -n $barcode | cut -d':' -f 1)

#sed '${lineno}s/$/$message/' file
