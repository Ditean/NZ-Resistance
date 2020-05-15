# Samtools script

# Set max thread of 8?

# NEED TO ADD LIMIT TO THREADS HERE?
if [[ $THREADS -ge 7 ]]
then
  THREADS=7 # (MAX SETTING)
  echo -e "WARNING: SAMtools thread adjusted to maximum of 7"
fi

samtools view -bq 10 -@ ${THREADS} $STAMP/${SAMPLE}/${SAMPLE}.sam > $STAMP/${SAMPLE}/${SAMPLE}.bam

samtools sort $STAMP/${SAMPLE}/${SAMPLE}.bam -n -O sam | samtools fixmate -m -O bam - $STAMP/${SAMPLE}/${SAMPLE}.fixmate.bam

samtools sort $STAMP/${SAMPLE}/${SAMPLE}.fixmate.bam | samtools markdup -r -S - $STAMP/${SAMPLE}/${SAMPLE}.fixmate.markdup.bam
