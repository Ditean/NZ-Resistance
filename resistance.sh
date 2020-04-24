#!/bin/bash

# Resistance wrapper

# Jordan Taylor - 23/04/2020

## THINGS TO DO:
# - re-setup script to run samples in an array. Each script runs once per pair
# [[ -z "$var" ]] && echo "Empty" || echo "Not empty" <- check if a variable is empty


# current date
function adddate(){
	date "+%d-%m-%Y %H:%M:%S %Z"
}

# Command line parsed arguments ================================================

# PART 1 -----------------------------------------------------------------------
# Default values for getops values - what we can do is use a 'is it empty loop'


# PART 2 -----------------------------------------------------------------------
# User input for files
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

		ADAPTOR

		*)
		>&2 printf "Error: Invalid argument\n"
		exit 1
		;;

		#*)
		#POSITIONAL+=("$1") # Store unknown options ----------------------------------
		#shift
		#;;
	esac
done
set -- "${POSITIONAL[@]}" # restore positional parameters ----------------------

# Part 1: Establish variables

# Part 1.1 - help message
if [[ $HELP == YES ]]
then
	echo -e "Help statement for the program\n- add options here\n- add other details\nverison 0.0.0 EARLY ALPHA \vCreated by: Jordan Taylor"
	exit 0
fi

# Part 1.2 - check correct user input

if [[ ( -z "$FORWARD" ) && ( -n "$REVERSE" ) ]] || [[ ( -n "$FORWARD" ) && ( -z "$REVERSE" ) ]]
then
	echo -e "ERROR: Provide both a forward and reverse sequence when providing a path"
	exit 1
fi


# Part 2: Exporting variables

# If R1 and R2
if [[ ( -n "$FORWARD" ) && ( -n "$REVERSE" ) ]] # -n = TRUE if string is nonzero || -z = TRUE if string is zero
then
	export FORWARD=$FORWARD && export REVERSE=$REVERSE
fi



echo -e "RESISTANCE pipeline\nINITILISING PROGRAM\n"

if [[ ( $SETUP == YES ) || ( ! -d $ROOT/refernce ) ]]
then
	echo "STEP 0: SETUP"
	bash ./scr/1.1.setup.sh
	echo -e "SETUP 0: COMPLETE - please check the log files for summary of setup\n======================================================\nPlease re-run the resistance pipeline"
	exit 0
fi


## Defining Environmental Paths

# ROOT - sets path to current directory
ROOT=$(dirname $(realpath $0))

# Run initial check for directories
# If variable = null and input is not present - run set-up script and exit
if && [ ! -d $ROOT/input ]
	bash ./.scr/1.1.setup.sh
fi
# 1.1.setup --------------------------------------------------------------------
# run setup of all directories and reference genomes

bash ./.scr/1.1.setup.sh

# 1.2.trimming.sh
# NEED TO FIGURE OUT A WAY TO MATCH FILES ------ MAKE THE SUB-SCRIPT PERFORM THE CHECK?

if [[ ]]

if [[ ( ! -f $FORWARD ) && ( ! -f ${FORWARD//{R,r}1/R2} ) ]
then
	echo "Can not detect paired files"
	exit 0
fi

if [ ! -f $REVERSE ] && [ ! -f ${REVERSE//{R,r}1/R2} ]
then
	echo "Can not detect paired files"
	exit 0
fi

# establish an array

declare -A sample_array

#i=1

#if [[ ( -n "$FORWARD" ) && ( -n "$REVERSE" ) ]]
#then
#	sample_array=([sample_${i}_R1]=$FORWARD [sample_${i}_R2]=$REVERSE)
#else
#	for
# fi
#fi

if [[ ( -n "$FORWARD" ) && ( -n "$REVERSE" ) ]] # double check the names match first
then
	R1=$FORWARD
	R2=${FORWARD%%R1*}R2${FORWARD#*R1}
	[[ $R1 == ]]
else
	j=1
	for i in `find $ROOT/input/ -name "*.fastq.gz" | grep R1`
	do # REWORK ARRAY: ONLY NEED TO HAVE THE BASENAME - GET RID OF BOTH R1 and R2
		[[ -f ${i//R1/R2}]] && sample_array=([sample_${j}_R1]=$i [sample_${j}_R2]=${i//R1/R2}); j++  || echo -e "ERROR: R2 file not found: ${i//R1/R2}"
	done
fi


bash ./.scr/1.2.trimming.sh -r1 $FORWARD -r2 $REVERSE -a $ADAPTOR
