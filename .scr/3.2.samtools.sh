# Samtools script

samtools view -bq 10 -@ ${THREADS} $STAMP/alignment/${SAMPLE}.sam > $STAMP/alignment/${SAMPLE}.bam

samtools sort $STAMP/alignment/${SAMPLE}.bam -n -O sam | samtools fixmate -m -O bam - $STAMP/alignment/${SAMPLE}.fixmate.bam

samtools sort $STAMP/alignment/${SAMPLE}.fixmate.bam | samtools markdup -r -S - $STAMP/alignment/${SAMPLE}.fixmate.markdup.bam
