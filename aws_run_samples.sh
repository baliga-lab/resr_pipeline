#!/bin/bash

if [ "$#" -ne 1 ]; then
    echo "usage: aws_run_samples.sh <inputfile>"
    exit 1
fi
infile=$1

while read l; do
    read -r sample url <<< "$l"
    resultsfile="$sample""-results.tar.gz"
    if test -f $resultsfile; then
        echo "Result exists, next"
        continue
    fi
    echo "Downloading Sample: $sample"
    mkdir -p $sample && cd $sample
    wget $url
    if [ $? -eq 0 ]
    then
	batchfile="$sample""_snp_calling.sh"
	full=`ls | head -n 1`
	echo "fasterq-dump --split-3 $full"
	fasterq-dump --split-3 $full
  	# handle fasterq-dump error !!!!!
  	if [ $? -eq 0 ]; then
  	    # same input dir and result dir for easier management
  	    cd .. && ./make_snp_calling.py $sample "$sample" && chmod u+x $batchfile && "./$batchfile"
	    aws s3 cp "$sample""-results.tar.gz" s3://baliga-bucket1
        rm -rf $resultsfile $sample "$sample""_snp_calling.sh"
  	else
  	    # cleanup and log the error
  	    cd .. && rm -rf $sample
  	    echo "$sample" >> "aws_failed_sample.txt"
  	    echo "FAILED SRA -> FASTQ conversion - skipping this sample and add to error list"
  	fi
    else
  	echo "$sample could not be downloaded"
  	echo "$sample" >> "aws_failed_sample.txt"
    fi
done < $infile
