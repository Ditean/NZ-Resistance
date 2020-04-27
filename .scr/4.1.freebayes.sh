# Freebayes script

freebayes -f $GENOME --ploidy 1 $STAMP/alignment/${SAMPLE}.fixmate.markdup.bam | vcffilter -f 'QUAL > 1 & QUAL / AO > 10 & SAF > 0 & SAR > 0 & RPR > 1 & RPL > 1' > $STAMP/vcf/${SAMPLE}.vcf
