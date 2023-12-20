# resR SNP Pipeline - A SNP based analysis pipeline

## Description

This is a data analysis pipeline based on the unfixed and fixed SNP 
calling pipelines of https://github.com/MtbEvolution/resR_Project
It was modified and augmented for easier configuration, running in various environments
and uniform results.

## Dependencies

  * pandas
  * globalsearch
  * all the dependencies of resR project
    * VarScan 2.4.0
    * samtools
    * BWA
    * sickle

## Usage

```
./make_snp_calling.py <path to FASTQ files>
``` 
