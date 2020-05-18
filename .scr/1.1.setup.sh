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

function adddate(){
	date "+%d-%m-%Y %H:%M:%S %Z"
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
dataset     CHECK
reference   CHECK

STATUS: PENDING

FILES:
.core-
     |- NC_000962.3.fna		CHECK
     |- contam.fa		CHECK
     |- mask_genome.bed   CHECK
     |- mutations  CHECK
               |- rifampicin.txt CHECK
               |- streptomycin.txt CHECK
               |- pyrazinamide.txt CHECK
               |- para-aminosalicylic.txt  CHECK
               |- linezolid.txt  CHECK
               |- kanamycin.txt  CHECK
               |- isoniazid.txt  CHECK
               |- fluoroquinolones.txt CHECK
               |- ethionamide.txt  CHECK
               |- ethambutol.txt CHECK
               |- capreomycin.txt  CHECK
               |- amikacin.txt CHECK

STATUS: PENDING

.scr-
    |-1.1.setup.sh		  CHECK
    |-2.1.trimming.sh   CHECK
    |-3.1.bwa.sh        CHECK
    |-3.2.samtools.sh   CHECK
    |-4.1.freebayes.sh  CHECK

STATUS: PENDING

reference-
	  			|- H37Rv.fa		CHECK
	  			|- adaptor.fa	CHECK
          |- genome.bed CHECK
          |- mutations  CHECK
                    |- rifampicin.txt CHECK
                    |- streptomycin.txt CHECK
                    |- pyrazinamide.txt CHECK
                    |- para-aminosalicylic.txt  CHECK
                    |- linezolid.txt  CHECK
                    |- kanamycin.txt  CHECK
                    |- isoniazid.txt  CHECK
                    |- fluoroquinolones.txt CHECK
                    |- ethionamide.txt  CHECK
                    |- ethambutol.txt CHECK
                    |- capreomycin.txt  CHECK
                    |- amikacin.txt CHECK

STATUS: PENDING

dataset-
        |- manifest.tsv CHECK
        |- resistance_profile.tsv CHECK

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

# CREATE SCR DIR
if [ ! -d $ROOT/.scr ]
then
	echo -e "STEP 2: SCR DIR created"
        mkdir $ROOT/.scr
        awk 'NR==1,/.scr/{sub(/CHECK/,"YES")}1' log_setup.txt > temp.txt && mv temp.txt log_setup.txt
else
        echo "STEP 2: SCR DIR already present"
        awk 'NR==1,/.scr/{sub(/CHECK/,"ALREADY_PRESENT")}1' log_setup.txt > temp.txt && mv temp.txt log_setup.txt
fi

# CREATE input DIR
if [ ! -d $ROOT/input ]
then
        echo -e "STEP 3: Input DIR created\n"
        mkdir $ROOT/input
        awk 'NR==1,/input/{sub(/CHECK/,"YES")}1' log_setup.txt > temp.txt && mv temp.txt log_setup.txt
else
        echo -e "STEP 3: Input DIR already present\n"
        awk 'NR==1,/input/{sub(/CHECK/,"ALREADY_PRESENT")}1' log_setup.txt > temp.txt && mv temp.txt log_setup.txt
fi

# CREATE output DIR
if [ ! -d $ROOT/output ]
then
        echo -e "STEP 4: OUTPUT DIR created\n"
        mkdir $ROOT/output
        awk 'NR==1,/output/{sub(/CHECK/,"YES")}1' log_setup.txt > temp.txt && mv temp.txt log_setup.txt
else
        echo -e "STEP 4: OUTPUT DIR already present\n"
        awk 'NR==1,/output/{sub(/CHECK/,"ALREADY_PRESENT")}1' log_setup.txt > temp.txt && mv temp.txt log_setup.txt
fi


# CREATE REFERENCE DIR
if [ ! -d $ROOT/reference ]
then
	echo -e "STEP 5: Reference DIR created"
        mkdir $ROOT/reference
        awk 'NR==1,/reference/{sub(/CHECK/,"YES")}1' log_setup.txt > temp.txt && mv temp.txt log_setup.txt
else
        echo "STEP 5: Reference DIR already present"
        awk 'NR==1,/reference/{sub(/CHECK/,"ALREADY_PRESENT")}1' log_setup.txt > temp.txt && mv temp.txt log_setup.txt
fi

# CREATE DATASET DIR
if [ ! -d $ROOT/dataset ]
then
  echo -e "STEP X: Dataset DIR created"
  mkdir $ROOT/dataset
  awk 'NR==1,/dataset/{sub(/CHECK/,"YES")}1' log_setup.txt > temp.txt && mv temp.txt log_setup.txt
else
  echo "STEP X: Dataset DIR already present"
  awk 'NR==1,/reference/{sub(/CHECK/,"ALREADY_PRESENT")}1' log_setup.txt > temp.txt && mv temp.txt log_setup.txt

# CHECK ALL DIRECTORIES ARE READY

setup_array=(.core .scr input output reference dataset)
for i in ${setup_array[@]}
do
  if [[ ! -d ${ROOT}/${i} ]]
  then
    echo -e "ERROR: Final check of directories failed\nPlease re-run setup again"
    awk 'NR==1,/STATUS:/{sub(/PENDING/,"FAILED")}1' log_setup.txt > temp.txt && mv temp.txt log_setup.txt
    exit 1
  fi
done
awk 'NR==1,/STATUS:/{sub(/PENDING/,"PASSED")}1' log_setup.txt > temp.txt && mv temp.txt log_setup.txt
# Move in the core files (references from HCS)
# The USER must have access to both HCS and have access to the Jordan directory in Htin Lab

echo -e "====================\nEstablishing required files\n====================\n"
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

# Reference genome
if [ ! -f $ROOT/.core/NC_000962.3.fa ]
then
  echo "Downloading reference genome"
  dataset
  awk 'NR==1,/NC_000962.3.fa/{sub(/CHECK/,"YES")}1' log_setup.txt > temp.txt && mv temp.txt log_setup.txt
else
  awk 'NR==1,/NC_000962.3.fa/{sub(/CHECK/,"ALREADY_PRESENT")}1' log_setup.txt > temp.txt && mv temp.txt log_setup.txt
fi

# Contams
if [ ! -f $ROOT/.core/contam.fa ]
then
  echo "Can not detect contam file - please download"
  awk 'NR==1,/contam.fa/{sub(/CHECK/,"ERROR")}1' log_setup.txt > temp.txt && mv temp.txt log_setup.txt
else
echo "Contam file present"
  awk 'NR==1,/contam.fa/{sub(/CHECK/,"YES")}1' log_setup.txt > temp.txt && mv temp.txt log_setup.txt
fi

# Mask
if [ ! -f $ROOT/.core/mask_genome.bed ]
then
  echo "ERROR: Can not detect mask_genome.bed"
  awk 'NR==1,/mask_genome.bed/{sub(/CHECK/,"ERROR")}1' log_setup.txt > temp.txt && mv temp.txt log_setup.txt
else
  echo "Masking file detected"
  awk 'NR==1,/mask_genome.bed/{sub(/CHECK/,"YES")}1' log_setup.txt > temp.txt && mv temp.txt log_setup.txt
fi

#### ADD IN A STATUS PENDING CHECK

echo -e "====================\nSetting up reference files\n====================\n"

# H37Rv reference genome
if [ ! -f $ROOT/reference/H37Rv.fa ]
then
	cp $ROOT/.core/NC_000962.3.fna $ROOT/reference/H37Rv.fa
	echo "Step 7: H37Rv reference setup"
	awk 'NR==1,/H37Rv.fa/{sub(/CHECK/,"YES")}1' log_setup.txt > temp.txt && mv temp.txt log_setup.txt
else
	echo "Step 7: H37Rv reference already present"
	awk 'NR==1,/H37Rv.fa/{sub(/CHECK/,"ALREADY_PRESENT")}1' log_setup.txt > temp.txt && mv temp.txt log_setup.txt
fi

# Sequencing adaptor reference
if [ ! -f $ROOT/reference/adaptor.fa ]
then
	cp $ROOT/.core/contam.fa $ROOT/reference/adaptor.fa
	echo -e "Step 8: Adaptor reference setup\n"
	awk 'NR==1,/adaptor.fa/{sub(/CHECK/,"YES")}1' log_setup.txt > temp.txt && mv temp.txt log_setup.txt
else
	echo -e "Step 8: Adaptor reference already present\n"
	awk 'NR==1,/adaptor.fa/{sub(/CHECK/,"ALREADY_PRESENT")}1' log_setup.txt > temp.txt && mv temp.txt log_setup.txt
fi

# Genome bed
if [ ! -f $ROOT/reference/genome.bed ]
then
  if [ -f $ROOT/.core/mask_genome.bed ]
  then
    cp $ROOT/.core/mask_genome.bed $ROOT/reference/genome.bed
    echo "Masking BED file set up"
    awk 'NR==1,/genome.bed/{sub(/CHECK/,"YES")}1' log_setup.txt > temp.txt && mv temp.txt log_setup.txt
  fi
else
  echo "Masking BED file already present"
  awk 'NR==1,/genome.bed/{sub(/CHECK/,"ALREADY_PRESENT")}1' log_setup.txt > temp.txt && mv temp.txt log_setup.txt
fi

# mutation files
if [ ! -d $ROOT/reference/mutations ]
then
  cp -R $ROOT/.core/mutations $ROOT/reference/
fi

# Resistance profile
if [ ! -f $ROOT/dataset/resistance_profile.tsv ]
then
  touch $ROOT/dataset/resistance_profile.tsv
  echo -e "BARCODE\tCLASS\tDRUG\tSTATUS\tGENE\tPOSITION\tREFERENCE\tALTERNATE\tPMID" > $ROOT/dataset/resistance_profile.tsv # Move to set-up
  awk 'NR==1,/resistance_profile.tsv/{sub(/CHECK/,"YES")}1' log_setup.txt > temp.txt && mv temp.txt log_setup.txt
else
  echo "Resistance profile detected"
  awk 'NR==1,/resistance_profile.tsv/{sub(/CHECK/,"ALREADY_PRESENT")}1' log_setup.txt > temp.txt && mv temp.txt log_setup.txt

# Manifest File
if [ ! -f $ROOT/dataset/manifest.tsv ]
then
  touch $ROOT/dataset/manifest.tsv
  echo -e "BARCODE\tFORWARD SEQ\tREVERSE SEQ\tSAMPLE ID\tPROCESSING DATE" > $ROOT/dataset/manifest.tsv
  awk 'NR==1,/manifest.tsv/{sub(/CHECK/,"YES")}1' log_setup.txt > temp.txt && mv temp.txt log_setup.txt
else
  echo "Manifest detected"
  awk 'NR==1,/manifest.tsv/{sub(/CHECK/,"ALREADY_PRESENT")}1' log_setup.txt > temp.txt && mv temp.txt log_setup.txt
fi

# Contam file - Need to think about a way to download this file from a repository
#if [ ! -f $ROOT/.core/contam.fa ]
#then
#        cp /home/jordan/mount-hcs/storage.hcs-p01.otago.ac.nz/micro-shared$/Htin\ Lab/jordan/reference/container/contam.fa $ROOT/.core/
#        echo -e "STEP 6: Adaptor contam downloaded\n"
#        awk 'NR==1,/contam.fa/{sub(/CHECK/,"YES")}1' log_setup.txt > temp.txt && mv temp.txt log_setup.txt
#else
#        echo -e "Step 6: Adaptor contam already present\n"
#        awk 'NR==1,/contam.fa/{sub(/CHECK/,"ALREADY_PRESENT")}1' log_setup.txt > temp.txt && mv temp.txt log_setup.txt
#fi

## REFERENCE DIR SET-UP


echo -e "====================\nSet-up is now complete\n$(adddate)\n====================\n"
