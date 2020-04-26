# PROJECT PROGRESS
# ================
#   - think about if you need .core files in report (maybe only reference)
#   - consider making user provided references be shifted in this script


# Script for creating all necessary directories

# Functions ====================================================================

# datasets ---------------------------------------------------------------------
# The dataset function will download the NCBI command-line program and download
# the required genome
function datasets(){
  curl -o datasets 'https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/LATEST/linux-amd64/datasets'

  chmod +x datasets
  ./datasets download assembly GCF_000195955.2
  unzip ncbi_dataset.zip
  mv ncbi_dataset/data/GCF_*/GCF_*_genomic.fna .core/NC_000962.3.fna

  rm -r ncbi_data*
  rm -r datasets
}

# Download container file and place within the .core directory + copy version into reference folder
function container(){
  curl someaddress
  mv contam.fa $ROOT/.core/contam.fa && cp $ROOT/.core/contam.fa $ROOT/reference/adaptor.fa
}

# ====================
echo -e "SETTING UP WORKFLOW \nSET-UP COMMENCED: $(adddate)\n"

# Create a log of the set-up
touch log_setup.txt

cat <<- EOF >> $ROOT/log_setup.txt
Set-up log for bovis transmission workflow

WORKFLOW SET-UP SUMMARY

USER: $USER
Date: $(adddate)

DIRECTORIES:
.core				CHECK
.scr				CHECK
input       CHECK
output      CHECK
reference   CHECK
reads				CHECK

STATUS: PENDING

FILES:
.core-
     |- NC_002945.4.fa		CHECK
     |- contam.fa		CHECK

STATUS: PENDING

.scr-
    |-script_1 			CHECK

STATUS: PENDING

reference-
	  			|- bovis.fa		CHECK
	  			|- adaptor.fa	CHECK
          |- genome.bed CHECK # still need to chek whether or not we have this

STATUS: PENDING


=========================================
END OF SET-UP FILE
=========================================
EOF

echo -e "====================\nCreating directories\n====================\n"


# CREATE CORE DIR - consider making this manditory, therefore if it isn't there ask them to repair the git clone
if [ ! -d $ROOT/.core ]
then
	echo -e "STEP 1: Core DIR created"
	mkdir $ROOT/.core
	awk 'NR==1,/.core/{sub(/CHECK/,"YES")}1' log_setup.txt > temp.txt && mv temp.txt log_setup.txt
else
	echo "STEP 1: Core DIR already present"
	awk 'NR==1,/.core/{sub(/CHECK/,"ALREADY_PRESENT")}1' log_setup.txt > temp.txt && mv temp.txt log_setup.txt
fi

# CREATE REFERENCE DIR
if [ ! -d $ROOT/reference ]
then
	echo -e "STEP 2: Reference DIR created"
        mkdir $ROOT/reference
        awk 'NR==1,/reference/{sub(/CHECK/,"YES")}1' log_setup.txt > temp.txt && mv temp.txt log_setup.txt
else
        echo "STEP 2: Reference DIR already present"
        awk 'NR==1,/reference/{sub(/CHECK/,"ALREADY_PRESENT")}1' log_setup.txt > temp.txt && mv temp.txt log_setup.txt
fi

# CREATE SCR DIR
if [ ! -d $ROOT/.scr ]
then
	echo -e "STEP 3: SCR DIR created"
        mkdir $ROOT/.scr
        awk 'NR==1,/.scr/{sub(/CHECK/,"YES")}1' log_setup.txt > temp.txt && mv temp.txt log_setup.txt
else
        echo "STEP 3: SCR DIR already present"
        awk 'NR==1,/.scr/{sub(/CHECK/,"ALREADY_PRESENT")}1' log_setup.txt > temp.txt && mv temp.txt log_setup.txt
fi

# CREATE input DIR
if [ ! -d $ROOT/input ]
then
        echo -e "STEP 4: Input DIR created\n"
        mkdir $ROOT/input
        awk 'NR==1,/input/{sub(/CHECK/,"YES")}1' log_setup.txt > temp.txt && mv temp.txt log_setup.txt
else
        echo -e "STEP 4: Input DIR already present\n"
        awk 'NR==1,/input/{sub(/CHECK/,"ALREADY_PRESENT")}1' log_setup.txt > temp.txt && mv temp.txt log_setup.txt
fi

# CHECK ALL DIRECTORIES ARE READY



# Move in the core files (references from HCS)
# The USER must have access to both HCS and have access to the Jordan directory in Htin Lab

## HCS TRANSFER
echo -e "====================\nTransfering files from HCS\n====================\n"
# H37Rv reference genome
if [ ! -f $ROOT/references/NC_000962.3.fa ]
then
  echo "STEP 5: Preparing H37Rv genome"
  if [ -f $ROOT/.core/NC_000962.3.fna ]
  then
    cp $ROOT/.core/NC_000962.3.fna $ROOT/reference/H37Rv.fa
    awk 'NR==1,/NC_000962.3.fa/{sub(/CHECK/,"YES")}1' log_setup.txt > temp.txt && mv temp.txt log_setup.txt
  else
    datasets
    awk 'NR==1,/NC_000962.3.fa/{sub(/CHECK/,"DOWNLOADED")}1' log_setup.txt > temp.txt && mv temp.txt log_setup.txt
  fi
else
	echo "Step 5: H37Rv genome detected"
	awk 'NR==1,/NC_002945.4.fa/{sub(/CHECK/,"ALREADY_PRESENT")}1' log_setup.txt > temp.txt && mv temp.txt log_setup.txt
fi

# Contam file - Need to think about a way to download this file from a repository
if [ ! -f $ROOT/.core/contam.fa ]
then
        cp /home/jordan/mount-hcs/storage.hcs-p01.otago.ac.nz/micro-shared$/Htin\ Lab/jordan/reference/container/contam.fa $ROOT/.core/
        echo -e "STEP 6: Adaptor contam downloaded\n"
        awk 'NR==1,/contam.fa/{sub(/CHECK/,"YES")}1' log_setup.txt > temp.txt && mv temp.txt log_setup.txt
else
        echo -e "Step 6: Adaptor contam already present\n"
        awk 'NR==1,/contam.fa/{sub(/CHECK/,"ALREADY_PRESENT")}1' log_setup.txt > temp.txt && mv temp.txt log_setup.txt
fi

## REFERENCE DIR SET-UP

echo -e "====================\nSetting up reference files\n====================\n"

# Bovis reference genome
if [ ! -f $ROOT/.reference/bovis.fa ]
then
	cp $ROOT/.core/NC_002945.4.fa $ROOT/.reference/bovis.fa
	echo "Step 7: Bovis reference setup"
	awk 'NR==1,/bovis.fa/{sub(/CHECK/,"YES")}1' log_setup.txt > temp.txt && mv temp.txt log_setup.txt
else
	echo "Step 7: Bovis reference already present"
	awk 'NR==1,/bovis.fa/{sub(/CHECK/,"ALREADY_PRESENT")}1' log_setup.txt > temp.txt && mv temp.txt log_setup.txt
fi

# Sequencing adaptor reference
if [ ! -f $ROOT/.reference/adaptor.fa ]
then
	cp $ROOT/.core/contam.fa $ROOT/.reference/adaptor.fa
	echo -e "Step 8: Adaptor reference setup\n"
	awk 'NR==1,/adaptor.fa/{sub(/CHECK/,"YES")}1' log_setup.txt > temp.txt && mv temp.txt log_setup.txt
else
	echo -e "Step 8: Adaptor reference already present\n"
	awk 'NR==1,/adaptor.fa/{sub(/CHECK/,"ALREADY_PRESENT")}1' log_setup.txt > temp.txt && mv temp.txt log_setup.txt
fi

echo -e "====================\nSet-up is now complete\n$(adddate)\n====================\n"
