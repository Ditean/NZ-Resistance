#!/bin/bash

# Resistance wrapper

# Jordan Taylor - 23/04/2020

# current date
function adddate(){
	date "+%d-%m-%Y %H:%M:%S %Z"
}

# Command line parsed arguments ================================================

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

		*)
		POSITIONAL+=("1") # Store unknown options ----------------------------------
		shift
		;;
	esac
done
set -- "${POSITIONAL[@]}" # restore positional parameters ----------------------

if [[ $HELP == YES ]]
then
	echo -e "Help statement for the program\n- add options here\n- add other details\nverison 0.0.0 EARLY ALPHA \vCreated by: Jordan Taylor"
	exit 0
fi

if [[ $SETUP == YES ]]
then
	echo "Beginning setup......."
	bash ./scr/1.1.setup.sh
	echo -e "Setup completed - please check the log files for summary of setup\n======================================================\nPlease re-run the resistance pipeline"
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
# Run itial setup of directories and scripts

bash ./.scr/1.1.setup.sh
