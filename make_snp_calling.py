#!/usr/bin/env python3

import argparse
import jinja2
import os
from globalsearch.rnaseq.find_files import find_fastq_files

DESCRIPTION = "make_fixed_snp_calling"

FIXED_TEMPLATE = """#!/bin/bash
sickle pe -l 35 -f {{fastq1}} -r {{fastq2}} -t sanger -o {{fixed_result_path}}/{{trimmed1}} -p {{fixed_result_path}}/{{trimmed2}} -s {{fixed_result_path}}/{{trimmedS}}
bwa aln -R 1 {{fasta_path}} {{fixed_result_path}}/{{trimmed1}} > {{fixed_result_path}}/{{sai1}}
bwa aln -R 1 {{fasta_path}} {{fixed_result_path}}/{{trimmed2}} > {{fixed_result_path}}/{{sai2}}
bwa aln -R 1 {{fasta_path}} {{fixed_result_path}}/{{trimmedS}} > {{fixed_result_path}}/{{saiS}}
bwa sampe -a 1000 -n 1 -N 1 {{fasta_path}} {{fixed_result_path}}/{{sai1}} {{fixed_result_path}}/{{sai2}} {{fixed_result_path}}/{{trimmed1}} {{fixed_result_path}}/{{trimmed2}} > {{fixed_result_path}}/{{paired_sam}}
bwa samse -n 1 {{fasta_path}} {{fixed_result_path}}/{{saiS}} {{fixed_result_path}}/{{trimmedS}} > {{fixed_result_path}}/{{single_sam}}
samtools view -bhSt {{fasta_path}}.fai {{fixed_result_path}}/{{paired_sam}} -o {{fixed_result_path}}/{{paired_bam}}
samtools view -bhSt {{fasta_path}}.fai {{fixed_result_path}}/{{single_sam}} -o {{fixed_result_path}}/{{single_bam}}
samtools merge {{fixed_result_path}}/{{merged_bam}} {{fixed_result_path}}/{{paired_bam}} {{fixed_result_path}}/{{single_bam}}
samtools sort {{fixed_result_path}}/{{merged_bam}} -o {{fixed_result_path}}/{{sorted_bam}}
samtools index {{fixed_result_path}}/{{sorted_bam}}
samtools mpileup -q 30 -Q 20 -ABOf {{fasta_path}} {{fixed_result_path}}/{{sorted_bam}} > {{fixed_result_path}}/{{pileup_file}}
java -jar {{varscan_path}} mpileup2snp {{fixed_result_path}}/{{pileup_file}} --min-coverage 3 --min-reads2 2 --min-avg-qual 20 --min-var-freq 0.01 --min-freq-for-hom 0.9 --p-value 99e-02 --strand-filter 1 > {{fixed_result_path}}/{{varscan_file}}
java -jar {{varscan_path}} mpileup2cns {{fixed_result_path}}/{{pileup_file}} --min-coverage 3 --min-avg-qual 20 --min-var-freq 0.75 --min-reads2 2 --strand-filter 1 > {{fixed_result_path}}/{{cns_file}}
perl {{resr_script_path}}/PPE_PE_INS_filt.pl {{resr_db_path}}/PPE_INS_loci.list {{fixed_result_path}}/{{varscan_file}} > {{fixed_result_path}}/{{ppe_file}}
perl {{resr_script_path}}/fixed_format_trans.pl {{fixed_result_path}}/{{ppe_file}} > {{fixed_result_path}}/{{for_file}}
perl {{resr_script_path}}/fix_extract.pl {{fixed_result_path}}/{{for_file}} > {{fixed_result_path}}/{{fix_file}}
perl {{resr_script_path}}/unfix_pileup_match.pl {{fixed_result_path}}/{{for_file}} {{fixed_result_path}}/{{pileup_file}} > {{fixed_result_path}}/{{forup_file}}
cut -f2-4 {{fixed_result_path}}/{{fix_file}} > {{fixed_result_path}}/{{snp_file}}
python {{resr_script_path}}/annotate_mtb_results.py {{fixed_result_path}}/{{snp_file}} {{fixed_result_path}}/{{varscan_file}} {{resr_db_path}} > {{fixed_result_path}}/{{result_file}}
"""

UNFIXED_TEMPLATE="""#!/bin/bash

sickle pe -l 35 -f {{fastq1}} -r {{fastq2}} -t sanger -o {{unfixed_result_path}}/{{trimmed1}} -p {{unfixed_result_path}}/{{trimmed2}} -s {{unfixed_result_path}}/{{trimmedS}}
bwa aln -R 1 {{fasta_path}} {{unfixed_result_path}}/{{trimmed1}} > {{unfixed_result_path}}/{{sai1}}
bwa aln -R 1 {{fasta_path}} {{unfixed_result_path}}/{{trimmed2}} > {{unfixed_result_path}}/{{sai2}}
bwa aln -R 1 {{fasta_path}} {{unfixed_result_path}}/{{trimmedS}} > {{unfixed_result_path}}/{{saiS}}
bwa sampe -a 1000 -n 1 -N 1 {{fasta_path}} {{unfixed_result_path}}/{{sai1}} {{unfixed_result_path}}/{{sai2}} {{unfixed_result_path}}/{{trimmed1}} {{unfixed_result_path}}/{{trimmed2}} > {{unfixed_result_path}}/{{paired_sam}}
bwa samse -n 1 {{fasta_path}} {{unfixed_result_path}}/{{saiS}} {{unfixed_result_path}}/{{trimmedS}} > {{unfixed_result_path}}/{{single_sam}}
samtools view -bhSt {{fasta_path}}.fai {{unfixed_result_path}}/{{paired_sam}} -o {{unfixed_result_path}}/{{paired_bam}}
samtools view -bhSt {{fasta_path}}.fai {{unfixed_result_path}}/{{single_sam}} -o {{unfixed_result_path}}/{{single_bam}}
samtools merge {{unfixed_result_path}}/{{merged_bam}} {{unfixed_result_path}}/{{paired_bam}} {{unfixed_result_path}}/{{single_bam}}
samtools sort {{unfixed_result_path}}/{{merged_bam}} -o {{unfixed_result_path}}/{{sorted_bam}}
samtools index {{unfixed_result_path}}/{{sorted_bam}}

#generate mpileup file for each sample
samtools mpileup -q 30 -Q 20 -ABOf {{fasta_path}} {{unfixed_result_path}}/{{sorted_bam}} > {{unfixed_result_path}}/{{pileup_file}}

#using varscan for SNP calling
java -jar {{varscan_path}} mpileup2snp {{unfixed_result_path}}/{{pileup_file}} --min-coverage 3 --min-reads2 2 --min-avg-qual 20 --min-var-freq 0.01 --min-freq-for-hom 0.9 --p-value 99e-02 --strand-filter 1 > {{unfixed_result_path}}/{{varscan_file}}
java -jar {{varscan_path}} mpileup2cns {{unfixed_result_path}}/{{pileup_file}} --min-coverage 3 --min-avg-qual 20 --min-var-freq 0.75 --min-reads2 2 --strand-filter 1 > {{unfixed_result_path}}/{{cns_file}}

#exclude regions belonging to PPE/PE and insertion sequences, and also exclude the regions that were recently marked as error-prone (Marin, Bioinformatics, 2022)
perl {{resr_script_path}}/PE_IS_filt.pl {{resr_db_path}}/Excluded_loci_mask.list {{unfixed_result_path}}/{{varscan_file}} > {{unfixed_result_path}}/{{ppe_file}}
perl {{resr_script_path}}/unfixed_format_trans.pl {{unfixed_result_path}}/{{ppe_file}} > {{unfixed_result_path}}/{{for_file}}

#extract read location from mpileup file (where does a mutation allele locate on a seqeuncing read), for further filtering based on tail distribution
perl {{resr_script_path}}/mix_pileup_merge.pl {{unfixed_result_path}}/{{for_file}} {{unfixed_result_path}}/{{pileup_file}} > {{unfixed_result_path}}/{{forup_file}}

#average sequencing depth, only include samples with genome coverage rate >0.9 and sequencing depth >20X
# the python script has the same effect as the sed/awk line, but is easier to understand
#sed 's/:/\\t/g' {{cns_file}}|awk '{if ($6 >= 3){n++;sum+=$6}} END {print \"\\t\",n/4411532,\"\\t\",sum/n}' > {{dep_file}}
python3 {{resr_script_path}}/avg_sequencing_depth.py {{unfixed_result_path}}/{{cns_file}} > {{unfixed_result_path}}/{{dep_file}}

#extract unfixed SNPs from forup files, this will create two files: "markdisc" and "markkept"; the suspected false positives(such as mutations with tail region enrichment) will be moved to markdisc file
perl {{resr_script_path}}/mix_extract_0.95.pl {{unfixed_result_path}}/{{forup_file}} > {{unfixed_result_path}}/{{mix_file}}
perl {{resr_script_path}}/forup_format.pl {{unfixed_result_path}}/{{mix_file}} > {{unfixed_result_path}}/{{mixfor_file}}
perl {{resr_script_path}}/info_mark.pl {{unfixed_result_path}}/{{mixfor_file}} > {{unfixed_result_path}}/{{mixmark_file}}
perl {{resr_script_path}}/redepin_filt.pl {{resr_db_path}}/Excluded_loci_mask.list {{unfixed_result_path}}/{{dep_file}} {{unfixed_result_path}}/{{mixmark_file}}

#filter list of highly repeated mutations with similar mutational frequency
#for those unfixed mutations that arise >=5 times in the 50K isolates, further check their reliability based on 1) the ratio in "markkept"; 2) the distribution of the mutational frequency.
cat {{unfixed_result_path}}/*mixmarkkept > {{unfixed_result_path}}/all_KEPT.txt
perl {{resr_script_path}}/loci_freq_count.pl {{unfixed_result_path}}/all_KEPT.txt > {{unfixed_result_path}}/kept_repeat.txt
cat {{unfixed_result_path}}/*mixmark > {{unfixed_result_path}}/all_MIX.txt
perl {{resr_script_path}}/loci_freq_count.pl {{unfixed_result_path}}/all_MIX.txt > {{unfixed_result_path}}/mix_repeat.txt
perl {{resr_script_path}}/repeat_number_merge.pl {{unfixed_result_path}}/mix_repeat.txt {{unfixed_result_path}}/kept_repeat.txt > {{unfixed_result_path}}/merge_kept_mix.txt
perl {{resr_script_path}}/ratio.pl {{unfixed_result_path}}/merge_kept_mix.txt > {{unfixed_result_path}}/merge_kept_mix_ratio.txt
# WW: Comparing to 5 can lead to 0 *per5up.txt files !!!
awk '$4>=5' {{unfixed_result_path}}/merge_kept_mix_ratio.txt |awk '$6>0.6'|cut -f1|while read i;do echo $i > {{unfixed_result_path}}/$i.per5up.txt;grep -w $i {{unfixed_result_path}}/all_KEPT.txt|cut -f12 >> {{unfixed_result_path}}/$i.per5up.txt;done
paste {{unfixed_result_path}}/*per5up.txt > {{unfixed_result_path}}/5up_0.6_paste.txt
perl {{resr_script_path}}/stdv.pl {{unfixed_result_path}}/5up_0.6_paste.txt |awk '$2<0.25'|cut -f1 > {{unfixed_result_path}}/5up_0.6_0.25.list
perl {{resr_script_path}}/freq_extract.pl {{unfixed_result_path}}/5up_0.6_0.25.list {{unfixed_result_path}}/5up_0.6_paste.txt > {{unfixed_result_path}}/5up_0.6_0.25.txt
awk '$4>=5' {{unfixed_result_path}}/merge_kept_mix_ratio.txt|cut -f1 > {{unfixed_result_path}}/5up.list
perl {{resr_script_path}}/repeat_loci.pl {{unfixed_result_path}}/5up_0.6_0.25.list {{unfixed_result_path}}/5up.list > {{unfixed_result_path}}/5up_remove_loc.list
perl {{resr_script_path}}/repeatloci_filter.pl {{unfixed_result_path}}/5up_remove_loc.list {{unfixed_result_path}}/{{markkept_file}} > {{unfixed_result_path}}/{{keptfilt_file}}

#annotation of unfixed SNPs
cut -f9-11 {{unfixed_result_path}}/{{keptfilt_file}} > {{unfixed_result_path}}/{{keptsnp_file}}
#perl {{resr_script_path}}/1_MTBC_Annotation_mtbc_4411532.pl {{keptsnp_file}} > {{keptanofilt_file}}
python {{resr_script_path}}/annotate_mtb_results.py {{unfixed_result_path}}/{{keptsnp_file}} {{unfixed_result_path}}/{{varscan_file}} {{resr_db_path}} > {{unfixed_result_path}}/{{result_file}}
tb-profiler profile --bam {{unfixed_result_path}}/ERR1023302.sorted.bam --dir {{tbprofiler_result_path}} --prefix {{sample_id}} --csv
"""

if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description=DESCRIPTION)
    parser.add_argument('--fastq_pattern', default="*_{{readnum}}.fastq", help="FASTQ pattern")
    parser.add_argument('input_path', help="input path to sample")
    parser.add_argument('result_path', help="result path")

    args = parser.parse_args()
    fastq_files = find_fastq_files(args.input_path, [args.fastq_pattern])
    if len(fastq_files) == 0:
        print("no FASTQ files found")
        exit(1)
    fq1, fq2 = fastq_files[0]
    basename1 = os.path.basename(fq1)
    basename2 = os.path.basename(fq2)
    stem1 = basename1.replace(".fastq", "").replace("gz", "")
    stem2 = basename2.replace(".fastq", "").replace("gz", "")

    # depending on the readnum format
    stem0 = stem1.replace("_1", "")
    fasta_path = "/bwa_pipeline/reference/MTB_ancestor_reference.fasta"
    resr_script_path = "scripts"
    resr_db_path = "databases"
    varscan_path = "/jarfiles/VarScan.v2.4.0.jar"

    fixed_templ = jinja2.Template(FIXED_TEMPLATE)
    fixed_result_path = os.path.join(args.result_path, "fixed")
    unfixed_result_path = os.path.join(args.result_path, "unfixed")
    tbprofiler_result_path = os.path.join(args.result_path, "tbprofiler")

    if not os.path.exists(args.result_path):
        os.makedirs(args.result_path)
        os.makedirs(unfixed_result_path)
        os.makedirs(fixed_result_path)

    config = {
        "sample_id": stem0,
        "fasta_path": fasta_path,
        "resr_script_path": resr_script_path,
        "resr_db_path": resr_db_path,
        "varscan_path": varscan_path,
        "fixed_result_path": fixed_result_path,
        "unfixed_result_path": unfixed_result_path,
        "tbprofiler_result_path": tbprofiler_result_path,
        "fastq1": fq1, "fastq2": fq2,
        "trimmed1": "%s_trimmed.fq" % stem1, "trimmed2": "%s_trimmed.fq" % stem2,
        "trimmedS": "%s_trimmedS.fq" % stem0,
        "sai1": "%s.sai" % stem1, "sai2": "%s.sai" % stem2, "saiS": "%s_S.sai" % stem0,
        "paired_sam": "%s.paired.sam" % stem0, "single_sam": "%s.single.sam" % stem0,
        "paired_bam": "%s.paired.bam" % stem0, "single_bam": "%s.single.bam" % stem0,
        "merged_bam": "%s.merged.bam" % stem0, "sorted_bam": "%s.sorted.bam" % stem0,
        "pileup_file": "%s.pileup" % stem0, "varscan_file": "%s.varscan" % stem0,
        "cns_file": "%s.cns" % stem0, "ppe_file": "%s.ppe" % stem0,
        "for_file": "%s.for" % stem0, "fix_file": "%s.fix" % stem0,
        "forup_file": "%s.forup" % stem0, "snp_file": "%s.snp" % stem0,
        "dep_file": "%s.dep" % stem0, "mix_file": "%s.mix" % stem0,
        "mixfor_file": "%s.mixfor" % stem0, "mixmark_file": "%s.mixmark" % stem0,
        "markkept_file": "%s.mixmarkkept" % stem0,
        "keptfilt_file": "%s.keptfilt" % stem0, "keptsnp_file": "%s.keptsnp" % stem0,
        "keptanofilt_file": "%s.keptanofilt" % stem0,
        "result_file": "%s.finalresult" % stem0
    }
    fixed_out = fixed_templ.render(config)
    with open("fixed_snp_calling.sh", "w") as outfile:
        outfile.write(fixed_out)

    unfixed_templ = jinja2.Template(UNFIXED_TEMPLATE)
    unfixed_out = unfixed_templ.render(config)
    with open("unfixed_snp_calling.sh", "w") as outfile:
        outfile.write(unfixed_out)
