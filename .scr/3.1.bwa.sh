# BWA script

# THINGS TO DO:
# - Check reference genome has been indexed
# - Threads max setting


# MAKE AN ARRAY WITH THE REF GENOME SUBS

# 3.1.1 - Set thread level for BWA mem
if [[ $THREADS -ge 7 ]]
then
  $THREADS=7 # (MAX SETTING)
  echo -e "WARNING: BWA thread adjusted to maximum of 7"
fi

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
bwa mem -t $THREADS $GENOME $STAMP/${SAMPLE}/${SAMPLE}.R1.fastq.gz $STAMP/${SAMPLE}/${SAMPLE}.R2.fastq.gz > $STAMP/${SAMPLE}/${SAMPLE}.sam

exit 0
