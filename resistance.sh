#!/bin/bash

# Resistance wrapper

# Jordan Taylor - 23/04/2020

## NOTES:
# - re-setup script to run samples in an array. Each script runs once per pair
# [[ -z "$var" ]] && echo "Empty" || echo "Not empty" <- check if a variable is empty
# ADD A FUNCTION FOR MASKING REFERENCE GENOME
#		- User can either provide a masked genome or provide an unmasked based on options?
#		- Add a function that will do the masking of the reference genome
#				> Try to make as much a function as possible

# 0: Initial set-up=============================================================

set -x # If return = non-zero, terminate parent script

# 0.1 - FUNCTIONS---------------------------------------------------------------

# 0.1.1 - date function
function adddate(){
	date "+%d-%m-%Y %H:%M:%S %Z"
}

# 0.1.2 - temporary directory for run
function date_temp(){
	date "+%d-%m-%Y"
}

# 0.1.3 - using datasets (NCBI) to download reference genome
function datasets(){
  curl -o datasets 'https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/LATEST/linux-amd64/datasets'

  chmod +x datasets
  ./datasets download assembly GCF_000195955.2
  unzip ncbi_dataset.zip
  mv ncbi_dataset/data/GCF_*/GCF_*_genomic.fna .core/NC_000962.3.fna

  rm -r ncbi_data*
  rm -r datasets
}

# 0.1.4 - download Illumina sequencing adaptors
function container(){
  curl someaddress
  mv contam.fa $ROOT/.core/contam.fa && cp $ROOT/.core/contam.fa $ROOT/reference/adaptor.fa
}


# 0.2 - Command line user input ------------------------------------------------
POSITIONAL=()
while [[ $# -gt 0 ]] # input is greater than 0
do
	key="$1"

	case $key in
		-h|--help)
		HELP=YES
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

		-a|--adaptor)
		ADAPTOR="$2"
		shift
		shift
		;;

		-t|--threads)
		THREADS="$2"
		shift
		shift
		;;

		-q|--quality)
		QUALITY="$2"
		shift
		shift
		;;

		--hoard)
		HOARD=YES
		shift
		;;

		# Consider adding a setting for --FORCE
		#			force the samples despite wrong names of files?

		# Consider a Snakemake options?

		*)
		>&2 printf "Error: Invalid argument\n"
		exit 1
		;;

	esac
done
set -- "${POSITIONAL[@]}" # restore positional parameters ----------------------


# 1.0: Set-up and checking user input===========================================

# 1.0.1 - help message CONSIDER AN EOF OPTION
if [[ $HELP == YES ]]
then
	echo -e "TB Resistance Profiler"
	echo -e "======================\n"
	echo -e "Usage: resistance.sh [ -r1 forward ] [ -r2 reverse ] [ -g genome ] [ -a adaptor ] [ -q quality ]\n"
	echo -e "  -r1"
	exit 0
fi

# 1.0.2 - Script working directory (ROOT)
ROOT=$(dirname $(realpath $0))
export ROOT=$ROOT

### CONSIDER ADDING ALL PATH VARIABLES HERE?


# 1.1.setup.sh - Run setup sub-script-------------------------------------------

if [[ ( $SETUP == YES ) || ( ! -d $ROOT/refernce ) ]]
then
	echo "STEP 0: SETUP"
	bash ./scr/1.1.setup.sh
	echo -e "SETUP 0: COMPLETE - please check the log files for summary of setup\n======================================================\nPlease re-run the resistance pipeline"
	exit 0
fi

# 1.2.1 - User provided R1 and R2
#					> Check files exist
#					> Check names match
#					> Set FASTA=TRUE (If TRUE - scripts use user provided R1 + R2)
#					> Create basename of sample
#					> Export names for sub-scritps to use

if [[ ( -z "$FORWARD" ) && ( -n "$REVERSE" ) ]] || [[ ( -n "$FORWARD" ) && ( -z "$REVERSE" ) ]] ## IF1
then
	echo -e "ERROR: Provide both a forward and reverse sequence when providing a path"
	exit 1
else

	if [[ ( -n "$FORWARD" ) && ( -n "$REVERSE" ) ]] && [[ ( -f $"FORWARD" ) && ( -f "$REVERSE" ) ]] # Check files actually exist ## IF2
	then
		base_forward=$(basename $FORWARD) # Take only the base name for comparison
		base_reverse=$(basename $REVERSE) # If files are in separate directories, non-base will not match

		# If the names match, export the forward and reverse #### CONSIDER MAKING BOTH IN SAME ARRAY - therefore one script in future
		if [[ $base_forward == ${base_reverse//R1/R2/} ]] ## IF 3
		then
			FASTA=TRUE # If TRUE, use user provided sequences
			base_forward=$(basename $FORWARD R1.fastq.gz) # Removes either '_' / '.' from before R1
			sample_ID=${base_forward::-1} # Basename for both R1 and R2 files

			# base_forward=$(basename $FORWARD)
			# temp=${base_forward%R1*}
			# base_forward=${temp::-1}
			# temp=${base_forward%L001*}
			# base_forward=${temp::-1}

			# Export variables
			export FASTA=$FASTA
			export FORWARD=$FORWARD
			export REVERSE=$REVERSE
			export SAMPLE=$sample_ID

		else
			echo -e "ERROR: $base_forward and $base_reverse do not match"
			exit 1
		fi ## FI 3

	else
		echo -e "ERROR: Can not detect forward and reverse sequences"
		exit 1

	fi ## FI 2

fi ## FI1

# 1.2.2 - Directory for processing sequences
run_stamp=date_temp
mkdir -p tmp/$run_stamp

export STAMP=$run_stamp

# 1.2.3 - Reference genome

if [[ -z $GENOME ]]
then
	GENOME=$ROOT/reference/H37Rv.fa && export GENOME=$GENOME
else
	export GENOME=$GENOME
fi

# 1.2.4 - Adaptor reference

if [[ -z $ADAPTOR ]]
then
	ADAPTOR=$ROOT/reference/adaptor.fa && export ADAPTOR=$ADAPTOR
else
	export ADAPTOR=$ADAPTOR
fi

# 2.0: Initialise pipeline======================================================

echo -e "RESISTANCE pipeline\nINITILISING PROGRAM\n"

# 2.0.1 - Create an Array of sample names from input DIR if FASTA=FALSE
if [[ ! $FASTA == TRUE ]]
	declare -a sample_array
	sample_array=($(find $ROOT/input/ -name "*fastq.gz" -printf "%f\n" | sed -i 's/*.R..fastq.gz//g' | uniq ))
	for i in ${sample_array[@]}
	do
		[[ ( -f $ROOT/${i}.R1.fastq.gz ) && ( -f $ROOT/${i}.R2.fastq.gz ) ]] || echo -e "$i does not have matching pair"; sample_array[$sample_array[(i)$i]]=() # Remove the index for sample in array
	done
fi

# 2.1.trimming.sh - run fastq-mcf trimming

if [[ $FASTA == TRUE ]]
then
	bash ./.scr/2.1.trimming.sh
else
	for i in ${sample_array[@]}
		do
			export SAMPLE=$i
			bash ./.scr/2.1.trimming.sh # Variables should be exported - so no need to provide flags
		done
fi
# DO A QC HERE?

# 2.1.bwa.sh - Align all reads to reference genome =============================
for i in ${sample_array[@]}
do
	export sample=$i
	bash ./.scr/3.1.bwa.sh # In code set limit for
done

# 2.2.samtools.sh - Convert to BAM and remove low quality alignments
bash ./.scr/3.2.samtools.sh

# 4.1.freebayes.sh - Perform variant calling on alignments

bash ./.scr/4.1.freebayes.sh

# Shift the initial VCF results

for i in ${sample_array[@]}
do
	cp $STAMP/${i}/${i}.vcf $ROOT/output/${i}.vcf
done

#### JUNK SCRIPT PILE ##########################################################
# Old code that was wrong / obsolete
#

# Getops
#*)
#POSITIONAL+=("$1") # Store unknown options --------------------------------
#shift
#;;

# 1.2.trimming.sh - Fastq-mcf trim of R1 and R2=================================
# NEED TO FIGURE OUT A WAY TO MATCH FILES ------ MAKE THE SUB-SCRIPT PERFORM THE CHECK?

#if [[ ( -n $FORWARD)]][[ ( ! -f $FORWARD ) && ( ! -f ${FORWARD//{R,r}1/R2} ) ]
#then
#	echo "Can not detect paired files"
#	exit 0
#fi

#if [ ! -f $REVERSE ] && [ ! -f ${REVERSE//{R,r}1/R2} ]
#then
#	echo "Can not detect paired files"
#	exit 0
#fi

# establish an array
#declare -A sample_array # Add this later

#then
#	sample_array=([sample_${i}_R1]=$FORWARD [sample_${i}_R2]=$REVERSE)
#else
#	for
# fi
#fi

#if [[ ( -n "$FORWARD" ) && ( -n "$REVERSE" ) ]] # double check the names match first
#then
#	R1=$FORWARD
#	R2=${FORWARD%%R1*}R2${FORWARD#*R1}
#	[[ $R2 == $REVERSE ]] || exit 1
#else
#	j=1
#	for i in `find $ROOT/input/ -name "*.fastq.gz" | grep R1`
#	do # REWORK ARRAY: ONLY NEED TO HAVE THE BASENAME - GET RID OF BOTH R1 and R2
#		[[ -f ${i//R1/R2}]] && sample_array=([sample_${j}_R1]=$i [sample_${j}_R2]=${i//R1/R2}); j++  || echo -e "ERROR: R2 file not found: ${i//R1/R2}"
#	done
#fi
#

# Run in array?
