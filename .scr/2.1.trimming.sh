# Script will take input from the wrapper, and trim files using fastq-mcf
#
# - Need to activate specific conda environment based on required environment
#
# - Consider having no default for R1 and R2 so that you can define here? Or you can make
#   an array in the previous script that is curated?
#
# - Could also make a quick check for count of R1 == R2

# FUNCTIONS
function adddate(){
	date "+%d-%m-%Y %H:%M:%S %Z"
}

function date_temp(){
	date "+%d-%m-%Y"
}

# 2.1 Script has begun =========================================================

echo -e "Trimming: INITIATED"

# 2.1.1.1 - Check for reference genome
if [ ! -f $GENOME ]
then
  echo -e "WARNING: Can not detect the reference genome!\nINVALIDE PATH: $GENOME"
  exit 1
fi

# 2.1.1.2 - Check for adaptor sequence file
if [ ! -f $ADAPTOR ]
then
  echo -e "WARNING: Can not detect the adaptor file!\nINVALID PATH: $ADAPTOR"
	exit 1
fi

# 2.1.2 Running Fastq-MCF ======================================================

# 2.1.2.1 - Initial trimming
if [[ $FASTA == TRUE ]]
do
	fastq-mcf $ADAPTOR $FORWARD $REVERSE -q $QUALITY \
	-o $STAMP/${SAMPLE}${SAMPLE}.R1.fastq.gz \
	-o $STAMP/${SAMPLE}/${SAMPLE}.R2.fastq.gz
else
	fastq-mcf $ADAPTOR $ROOT/${SAMPLE}*.fastq.gz -q $QUALITY \
	-o $STAMP/${SAMPLE}/${SAMPLE}.R1.fastq.gz \
	-o $STAMP/${SAMPLE}/${SAMPLE}.R2.fastq.gz
fi

exit 0
