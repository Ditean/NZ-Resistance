#!/bin/bash

# Resistance wrapper

# Jordan Taylor - 23/04/2020

## THINGS TO DO:
# - re-setup script to run samples in an array. Each script runs once per pair
# [[ -z "$var" ]] && echo "Empty" || echo "Not empty" <- check if a variable is empty

# Part 0: Initial set-up========================================================

# Part 0.1 - FUNCTIONS----------------------------------------------------------

# Part 0.1.1 - date function
function adddate(){
	date "+%d-%m-%Y %H:%M:%S %Z"
}

# Part 0.2 - Getops user input--------------------------------------------------
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

		#*)
		#POSITIONAL+=("$1") # Store unknown options --------------------------------
		#shift
		#;;
	esac
done
set -- "${POSITIONAL[@]}" # restore positional parameters ----------------------


# Part 0.2.1 - help message-----------------------------------------------------
if [[ $HELP == YES ]]
then
	echo -e "TB Resistance Profiler"
	echo -e "======================\n"
	echo -e "Usage: resistance.sh [ -r1 forward ] [ -r2 reverse ] [ -g genome ] [ -a adaptor ] [ -q quality ]\n"
	echo -e "  -r1"
	exit 0
fi

# Part 0.2.2 - check correct user input-----------------------------------------

if [[ ( -z "$FORWARD" ) && ( -n "$REVERSE" ) ]] || [[ ( -n "$FORWARD" ) && ( -z "$REVERSE" ) ]]
then
	echo -e "ERROR: Provide both a forward and reverse sequence when providing a path"
	exit 1
fi

# Part 0.3 - ENVIRONMENTAL VARIABLES--------------------------------------------

# Part 0.3.1 - Root DIR---------------------------------------------------------
ROOT=$(dirname $(realpath $0))

# Part 0.4: Exporting variables=================================================

# Part 0.4.1 - Export DIR path--------------------------------------------------
export ROOT=$ROOT

# Part 0.4.2 - Export user provided R1 and R2-----------------------------------
if [[ ( -n "$FORWARD" ) && ( -n "$REVERSE" ) ]] # -n = TRUE if string is nonzero || -z = TRUE if string is zero
then # Add a check here for matching names
	export FORWARD=$FORWARD && export REVERSE=$REVERSE
fi

echo -e "RESISTANCE pipeline\nINITILISING PROGRAM\n"


# 1.1.setup.sh - Run setup sub-script===========================================

if [[ ( $SETUP == YES ) || ( ! -d $ROOT/refernce ) ]]
then
	echo "STEP 0: SETUP"
	bash ./scr/1.1.setup.sh
	echo -e "SETUP 0: COMPLETE - please check the log files for summary of setup\n======================================================\nPlease re-run the resistance pipeline"
	exit 0
fi

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
declare -A sample_array # Add this later

# Some command to see if to use input or variable
declare -a sample_array
sample_array=($(find $ROOT/input/ -name "*fastq.gz" -printf "%f\n" | sed -i 's/*.R..fastq.gz//g' | uniq ))

# Now have the array - now check if a file does not exit

for i in ${sample_array[@]}
do
	[[ ( -f $ROOT/${i}.R1.fastq.gz ) && ( -f $ROOT/${i}.R2.fastq.gz ) ]] || echo -e "$i does not have matching pair"; sample_array[$sample_array[(i)$i]]=() # Remove the index for sample in array
done

## NEED TO EXPORT THE ARRAYS

#i=1

#if [[ ( -n "$FORWARD" ) && ( -n "$REVERSE" ) ]]
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
bash ./.scr/1.2.trimming.sh # Variables should be exported - so no need to provide flags

# DO A QC HERE?

# 2.1.bwa.sh - Align all reads to reference genome =============================

bash ./.scr/2.1.bwa.sh

# 2.2.samtools.sh - Convert to BAM and remove low quality alignments

# 3.1.freebayes.sh - Perform variant calling on alignments
