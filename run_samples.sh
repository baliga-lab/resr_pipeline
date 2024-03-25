#!/bin/bash

#infile="test-accessions.tsv"
#infile="test-accessions-20240117"
if [ "$#" -ne 1 ]; then
    echo "usage: run_samples.sh <inputfile>"
    exit 1
fi
infile=$1

while read l; do
  path="/proj/omics4tb2/wwu/SRA_DOWNLOADS/SRALite_REAL/$l*"
  ls $path
  if [ $? -eq 0 ]
  then
    full=`ls $path | head -n 1`
    if [ -f $full ];
    then
	batchfile="$l""_snp_calling.sh"
	echo "$full exists"
	mkdir -p $l && cd $l
	fasterq-dump --split-3 $full
	# handle fasterq-dump error !!!!!
	if [ $? -eq 0 ]; then
	    # same input dir and result dir for easier management
	    cd .. && ./make_snp_calling.py $l "$l" && chmod u+x $batchfile && sbatch $batchfile
	else
	    # cleanup and log the error
	    cd .. && rm -rf $l
	    echo "$l" >> "fasterq_dump_failed_conversion.txt"
	    echo "FAILED SRA -> FASTQ conversion - skipping this sample and add to error list"
	fi
    else
	echo "$full not exists"
    fi
  else
    echo "skipping sample $l"
  fi

done < $infile
