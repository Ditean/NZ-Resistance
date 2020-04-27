# BWA script

# Check reference genome has been indexed

# MAKE AN ARRAY WITH THE REF GENOME SUBS

index=(amb ann bwt fai pac sa)

for i in ${index[@]}
do
  if [[ ! -f ${GENOME}.${i} ]]
  then
    bwa index $GENOME
    break
  fi
done

# Run BWA
bwa mem -t $THREADS $GENOME $STAMP/${SAMPLE}.R1.fastq.gz $STAMP/${SAMPLE}.R2.fastq.gz

exit 0
