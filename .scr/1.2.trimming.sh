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

POSITIONAL=()
while [[ $# -gt 0 ]]
do
	key="$1"

	case $key in
		-h|--help)
		HELP=YES
		shift
		shift
		;;

		--setup) # Run 1.1.setup.sh ------------------------------------------------
		SETUP=YES
		shift
		;;

		-r1|--forward) # User provided R1 ------------------------------------------
		FORWARD="$2"
		shift
		shift
		;;

		-r2|--reverse) # User provided R2 ------------------------------------------
		REVERSE="$2"
		shift
		shift
		;;

		-g|--genome) # User provided reference genome ------------------------------
		GENOME="$2"
		shift
		shift
		;;

		ADAPTOR
		*)
		POSITIONAL+=("1") # Store unknown options ----------------------------------
		shift
		;;
	esac
done
set -- "${POSITIONAL[@]}" # restore positional parameters ----------------------

# May not need this as it is already set?

# 1.2.1 Script has begun =======================================================

echo -e "Trimming: INITIATED"

# 1.2.2.1 Checking for PE filename match ---------------------------------------
if [ ! -f $FORWARD == ${REVERSE%%R2*}R1${REVERSE#*R2} ]
then
  echo -e "WARNING: Paired reads do not match!\n$FORWARD\n$REVERSE"
  exit 0
fi

# 1.2.2.2 Checking for reference files -----------------------------------------
if [ ! -f $GENOME ]
then
  echo -e "WARNING: Can not detect the reference genome!\nINVALIDE PATH: $GENOME"
  exit 0
fi

if [ ! -f $ADAPTOR ]
then
  echo -e "WARNING: Can not detect the adaptor file!\nINVALID PATH: $ADAPTOR"
fi

# 1.2.3 Running Fastq-MCF ======================================================

mkdir -p $STAMP/fastq-mcf

# Initial trimming
if [[ $FASTA == TRUE ]]
do
	fastq-mcf $ADAPTOR $FORWARD $REVERSE -q $QUALITY \
	-o $STAMP/fastq-mcf/${sample}.R1.fastq.gz \
	-o $STAMP/fastq-mcf/${sample}.R2.fastq.gz
else
	fastq-mcf $ADAPTOR $ROOT/${i}*.fastq.gz -q $QUALITY \
	-o $STAMP/fastq-mcf/${i}.R1.fastq.gz \
	-o $STAMP/fastq-mcf/${i}.R2.fastq.gz
fi

exit 0
