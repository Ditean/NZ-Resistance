#!/bin/bash

# Resistance wrapper

# Jordan Taylor - 23/04/2020

# current date
function adddate(){
	date "+%d-%m-%Y %H:%M:%S %Z"
}

## Defining Environmental Paths

# ROOT - sets path to current directory
ROOT=$(dirname $(realpath $0))

# Run itial setup of directories and scripts

bash ./.scr/1.1.setup.sh
