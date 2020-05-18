#!/bin/bash

# Resistance wrapper

# Jordan Taylor - 23/04/2020

#################################################
#	Remaining Jobs:
#		- Clean up getopts
#				> Hoard mode & Metadata?
#		- Establish resistance calling script
#		- Solution to downloading container and BED file
#		- Metadata for report
#################################################

# Settings
#set -x # Debug mode


# 0: Initial set-up=============================================================

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

function metadata(){
	while true
	do
		echo "Please provide report metadata"
		read -p "BARCODE: " barcode
		read -p "Sequencing Date: " seq_date
		read -p " " something

		echo "Provided details"
		echo "Barcode: $barcode"
		echo "Sequencing Date: $seq_date"

		read -p "Continue (Y/N)" answer
		if [[ $answer == "Y" ]]
		then
			break
		fi
	done
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

		-a|--adaptor) # Adaptors for removal in trimming - must be FASTA format ----
		ADAPTOR="$2"
		shift
		shift
		;;

		--meta) # User provided metadata for report --------------------------------
		META="$2"
		shift
		shift
		;;

		-t|--threads)
		THREADS="$2"
		shift
		shift
		;;

		-q|--quality) # Phred Score used in trimming -------------------------------
		QUALITY="$2"
		shift
		shift
		;;

		--hoard) # Do not delete any intermediary after running --------------------
		HOARD=YES
		shift
		;;

		--mask) # If yes, mask reference genome (only call resistance SNPs) --------
		MASK=YES
		shift
		;;

		--bed) # User provided BED for masking -------------------------------------
		BED="$2"  # CONSIDER HAVING THIS SET AS A DEFINED VALUE?
		shift
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

# TEMP
QUALITY=30
THREADS=8

export QUALITY=$QUALITY
export THREADS=$THREADS

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

BED=$ROOT/reference/genome.bed
export BED=$BED

### CONSIDER ADDING ALL PATH VARIABLES HERE?


# 1.1.setup.sh - Run setup sub-script-------------------------------------------

if [[ ( $SETUP == YES ) || ( ! -d $ROOT/reference ) ]] # If [[ ( -n $SETUP )]]
then
	echo "STEP 0: SETUP"
	bash ./.scr/1.1.setup.sh
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
	echo -e "ERROR: Provide both a forward and reverse sequence"
	exit 1
fi

if [[ ( -n "$FORWARD" ) && ( -n "$REVERSE" ) ]] && [[ ( -f "$FORWARD" ) && ( -f "$REVERSE" ) ]] # Check files actually exist ## IF2
then
	base_forward=$(basename $FORWARD) # Take only the base name for comparison
	base_reverse=$(basename $REVERSE) # If files are in separate directories, non-base will not match

		# If the names match, export the forward and reverse #### CONSIDER MAKING BOTH IN SAME ARRAY - therefore one script in future
# ========================================
	if [[ $base_forward == ${base_reverse//R1/R2/} ]] ## IF 3
	then
		FASTA=TRUE # If TRUE, use user provided sequences
		sample_ID=$(echo $base_forward | sed 's/.R1.*fastq.gz//')

			# Export variables
		export FASTA=$FASTA
		export FORWARD=$FORWARD
		export REVERSE=$REVERSE
		export SAMPLE=$sample_ID

	else
		echo -e "ERROR: $base_forward and $base_reverse do not match"
		exit 1
	fi ## FI 3
# ==========================================
elif [[ ( -n "$FORWARD" ) && ( -n "$REVERSE" ) ]] && [[ ( ! -f "$FORWARD" ) || ( ! -f "$REVERSE" ) ]]
then
	echo -e "ERROR: Can not detect forward and reverse sequences"
	exit 1
fi

	#statements
# 1.2.2 - Directory for processing sequences
run_stamp=$(date_temp)
mkdir -p $ROOT/tmp/$run_stamp

export STAMP=${ROOT}/tmp/${run_stamp}

# 1.2.3 - Reference genome -- need to check if masking is required

# MASK CHECK

if [[ -n $MASK ]]
then
	if [[ -z $GENOME ]]
	then
		maskfasta -fi $ROOT/reference/H37Rv.fa -bed $BED -fo $ROOT/reference/masked_H37Rv.fa && GENOME=$ROOT/reference/masked_H37Rv.fa
	else
		maskfasta -fi $GENOME -bed $BED -f $ROOT/reference/masked_user_genome.fa && GENOME=$ROOT/reference/masked_user_genome.fa
	fi
elif [[ ( -z $MASK ) && ( -z $GENOME ) ]]
	then
		GENOME=$ROOT/reference/H37Rv.fa
fi

export GENOME=$GENOME

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
then

	declare -a sample_array
	declare -a mismatch_array
	declare -a man_array

	input_array=($(find $ROOT/input/ -name "*fastq.gz" -printf "%f\n" | grep "R1"))

	for i in ${input_array[@]}
	do
		if [[ ( -f $ROOT/input/${i} ) && ( -f $ROOT/input/${i//R1/R2} ) ]]
		then
			sample_array+=($(echo $i | sed 's/.R1.*fastq.gz//'))
			man_array+=($(echo $i))
		else
			mismatch_array+=($(echo ${i//R1/R2}))
			echo "ERROR: Can not detect ${i//R1/R2}"
		fi
	done
fi

# BARCODE
if [[ $FASTA == TRUE ]]
then
	linenumber=$(sed -n '$=' $ROOT/dataset/resistance_profile.tsv)
	if [ $linenumber == 1 ]
	then
		barcode=1000
	else
		barcode=$(cat $ROOT/dataset/resistance_profile.tsv | tail -n 1 | cut -f 1)
		addup=$((barcode++))
	fi

	echo -e "${barcode}\t${base_forward}\t${base_reverse}\t${sample_ID}\t$(adddate)" >> $ROOT/dataset/manifest.tsv

elif [[ ! $FASTA == TRUE ]]
then
	linenumber=$(sed -n '$=' $ROOT/dataset/resistance_profile.tsv)
	if [ $linenumber == 1 ]
	then
		barcode=1000
	else
		barcode=$(cat $ROOT/dataset/resistance_profile.tsv | tail -n 1 | cut -f 1)
		addup=$((barcode++))
	fi
	for i in ${man_array[@]}
	do
		name=$(echo $i | sed 's/.R1.*fastq.gz//')
		echo -e "$barcode\t$i\t${i//R1/R2}\t$name\t$(adddate)" >> $ROOT/dataset/manifest.tsv
		addup=$((barcode++))
	done
fi

# 2.1.trimming.sh - run fastq-mcf trimming
# Requires the following variables
#	- FASTA
#			> R1 and R2
#			> Basename of file
#	- SAMPLE
# - Adaptor
# - Quality
# -

if [[ $FASTA == TRUE ]]
then
	mkdir ${STAMP}/$sample_ID
	bash ./.scr/2.1.trimming.sh
else
	for i in ${sample_array[@]}
		do
			export SAMPLE=$i
			mkdir ${STAMP}/$i
			bash ./.scr/2.1.trimming.sh # Variables should be exported - so no need to provide flags
		done
fi
# DO A QC HERE?

# 2.1.bwa.sh - Align all reads to reference genome =============================
for i in ${sample_array[@]}
do
	export SAMPLE=$i
	bash ./.scr/3.1.bwa.sh # In code set limit for
done

# 2.2.samtools.sh - Convert to BAM and remove low quality alignments
for i in ${sample_array[@]}
do
	export SAMPLE=$i
	bash ./.scr/3.2.samtools.sh
done

# 4.1.freebayes.sh - Perform variant calling on alignments
for i in ${sample_array[@]}
do
	export SAMPLE=$i
	bash ./.scr/4.1.freebayes.sh
done

# Shift the initial VCF results

for i in ${sample_array[@]}
do
	cp $STAMP/${i}/${i}.vcf $ROOT/output/${i}.vcf
done

for i in ${sample_array[@]}
do
	export SAMPLE=$i
	export run_stamp=$run_stamp
	bash ./.scr/profile.sh
done

echo "The Pipeline has now finished"

# Call resistance mutations
