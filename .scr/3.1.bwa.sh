# BWA script

# Check reference genome has been indexed

# MAKE AN ARRAY WITH THE REF GENOME SUBS

for i in ${index[@]}
do
  if [[ ! -f $i ]]
  then
    bwa index $GENOME
    break
  fi
done

# Run BWA
bwa mem -t $THREADS
