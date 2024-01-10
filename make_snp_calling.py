#!/usr/bin/env python3

import argparse
import jinja2
import os
import json
from globalsearch.rnaseq.find_files import find_fastq_files


DESCRIPTION = "make_fixed_snp_calling"

SCRIPT_TEMPLATE_PRE = """#!/bin/bash

if [ ! -d "{{common_result_path}}" ]; then
  mkdir -p {{common_result_path}}
fi
if [ ! -d "{{fixed_result_path}}" ]; then
  mkdir -p {{fixed_result_path}}
fi
if [ ! -d "{{unfixed_result_path}}" ]; then
  mkdir -p {{unfixed_result_path}}
fi
if [ ! -d "{{tbprofiler_result_path}}" ]; then
  mkdir -p {{tbprofiler_result_path}}
fi
"""

SCRIPT_TEMPLATE_PAIRED = """
sickle pe -l 35 -f {{fastq1}} -r {{fastq2}} -t sanger -o {{common_result_path}}/{{trimmed1}} -p {{common_result_path}}/{{trimmed2}} -s {{common_result_path}}/{{trimmedS}}
bwa aln -R 1 {{fasta_path}} {{common_result_path}}/{{trimmed1}} > {{common_result_path}}/{{sai1}}
bwa aln -R 1 {{fasta_path}} {{common_result_path}}/{{trimmed2}} > {{common_result_path}}/{{sai2}}
bwa aln -R 1 {{fasta_path}} {{common_result_path}}/{{trimmedS}} > {{common_result_path}}/{{saiS}}
bwa sampe -a 1000 -n 1 -N 1 {{fasta_path}} {{common_result_path}}/{{sai1}} {{common_result_path}}/{{sai2}} {{common_result_path}}/{{trimmed1}} {{common_result_path}}/{{trimmed2}} > {{common_result_path}}/{{paired_sam}}
bwa samse -n 1 {{fasta_path}} {{common_result_path}}/{{saiS}} {{common_result_path}}/{{trimmedS}} > {{common_result_path}}/{{single_sam}}
samtools view -bhSt {{fasta_path}}.fai {{common_result_path}}/{{paired_sam}} -o {{common_result_path}}/{{paired_bam}}
samtools view -bhSt {{fasta_path}}.fai {{common_result_path}}/{{single_sam}} -o {{common_result_path}}/{{single_bam}}
samtools merge {{common_result_path}}/{{merged_bam}} {{common_result_path}}/{{paired_bam}} {{common_result_path}}/{{single_bam}}
samtools sort {{common_result_path}}/{{merged_bam}} -o {{common_result_path}}/{{sorted_bam}}
"""

SCRIPT_TEMPLATE_SINGLE = """#!/bin/bash
sickle se -l 35 -f {{fastq}} -t sanger -o {{common_result_path}}/{{trimmed}}
bwa aln -R 1 {{fasta_path}} {{common_result_path}}/{{trimmed}} > {{common_result_path}}/{{sai}}
bwa samse -n 1 {{fasta_path}} {{common_result_path}}/{{sai}} {{common_result_path}}/{{trimmed}} > {{common_result_path}}/{{single_sam}}
samtools view -bhSt {{fasta_path}}.fai {{common_result_path}}/{{single_sam}} -o {{common_result_path}}/{{single_bam}}
samtools sort {{common_result_path}}/{{single_bam}} -o {{common_result_path}}/{{sorted_bam}}
"""

SCRIPT_TEMPLATE_COMMON = """

# Start here with sorted BAM
samtools index {{common_result_path}}/{{sorted_bam}}
samtools mpileup -q 30 -Q 20 -ABOf {{fasta_path}} {{common_result_path}}/{{sorted_bam}} > {{common_result_path}}/{{pileup_file}}
java -jar {{varscan_path}} mpileup2snp {{common_result_path}}/{{pileup_file}} --min-coverage 3 --min-reads2 2 --min-avg-qual 20 --min-var-freq 0.01 --min-freq-for-hom 0.9 --p-value 99e-02 --strand-filter 1 > {{common_result_path}}/{{varscan_file}}
java -jar {{varscan_path}} mpileup2cns {{common_result_path}}/{{pileup_file}} --min-coverage 3 --min-avg-qual 20 --min-var-freq 0.75 --min-reads2 2 --strand-filter 1 > {{common_result_path}}/{{cns_file}}

# FIXED specific part
perl {{resr_script_path}}/PPE_PE_INS_filt.pl {{resr_db_path}}/PPE_INS_loci.list {{common_result_path}}/{{varscan_file}} > {{fixed_result_path}}/{{ppe_file}}
perl {{resr_script_path}}/fixed_format_trans.pl {{fixed_result_path}}/{{ppe_file}} > {{fixed_result_path}}/{{for_file}}
perl {{resr_script_path}}/fix_extract.pl {{fixed_result_path}}/{{for_file}} > {{fixed_result_path}}/{{fix_file}}
perl {{resr_script_path}}/unfix_pileup_match.pl {{fixed_result_path}}/{{for_file}} {{common_result_path}}/{{pileup_file}} > {{fixed_result_path}}/{{forup_file}}
cut -f2-4 {{fixed_result_path}}/{{fix_file}} > {{fixed_result_path}}/{{snp_file}}
python {{resr_script_path}}/annotate_mtb_results.py {{fixed_result_path}}/{{snp_file}} {{common_result_path}}/{{varscan_file}} {{resr_db_path}} > {{fixed_result_path}}/{{result_file}}

# Unfixed part
# use cns and varscan file from fixed

#exclude regions belonging to PPE/PE and insertion sequences, and also exclude the regions that were recently marked as error-prone (Marin, Bioinformatics, 2022)
#perl {{resr_script_path}}/PE_IS_filt.pl {{resr_db_path}}/Excluded_loci_mask.list {{common_result_path}}/{{varscan_file}} > {{unfixed_result_path}}/{{ppe_file}}
python {{resr_script_path}}/varscan_drop_excluded.py {{resr_db_path}}/Excluded_loci_mask.list {{common_result_path}}/{{varscan_file}} {{unfixed_result_path}}/{{ppe_file}}
perl {{resr_script_path}}/unfixed_format_trans.pl {{unfixed_result_path}}/{{ppe_file}} > {{unfixed_result_path}}/{{for_file}}

#extract read location from mpileup file (where does a mutation allele locate on a seqeuncing read), for further filtering based on tail distribution
perl {{resr_script_path}}/mix_pileup_merge.pl {{unfixed_result_path}}/{{for_file}} {{common_result_path}}/{{pileup_file}} > {{unfixed_result_path}}/{{forup_file}}

#average sequencing depth, only include samples with genome coverage rate >0.9 and sequencing depth >20X
# the python script has the same effect as the sed/awk line, but is easier to understand
#sed 's/:/\\t/g' {{cns_file}}|awk '{if ($6 >= 3){n++;sum+=$6}} END {print \"\\t\",n/4411532,\"\\t\",sum/n}' > {{dep_file}}
python3 {{resr_script_path}}/avg_sequencing_depth.py {{common_result_path}}/{{cns_file}} > {{unfixed_result_path}}/{{dep_file}}

#extract unfixed SNPs from forup files, this will create two files: "markdisc" and "markkept"; the suspected false positives(such as mutations with tail region enrichment) will be moved to markdisc file
perl {{resr_script_path}}/mix_extract_0.95.pl {{unfixed_result_path}}/{{forup_file}} > {{unfixed_result_path}}/{{mix_file}}
perl {{resr_script_path}}/forup_format.pl {{unfixed_result_path}}/{{mix_file}} > {{unfixed_result_path}}/{{mixfor_file}}
perl {{resr_script_path}}/info_mark.pl {{unfixed_result_path}}/{{mixfor_file}} > {{unfixed_result_path}}/{{mixmark_file}}

# this script implicitly generates ".mixmarkkept" and ".mixmarkfilt"
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
#awk '$4>=5' {{unfixed_result_path}}/merge_kept_mix_ratio.txt |awk '$6>0.6'|cut -f1|while read i;do echo $i > {{unfixed_result_path}}/$i.per5up.txt;grep -w $i {{unfixed_result_path}}/all_KEPT.txt|cut -f12 >> {{unfixed_result_path}}/$i.per5up.txt;done
paste {{unfixed_result_path}}/*per5up.txt > {{unfixed_result_path}}/5up_0.6_paste.txt
perl {{resr_script_path}}/stdv.pl {{unfixed_result_path}}/5up_0.6_paste.txt |awk '$2<0.25'|cut -f1 > {{unfixed_result_path}}/5up_0.6_0.25.list
perl {{resr_script_path}}/freq_extract.pl {{unfixed_result_path}}/5up_0.6_0.25.list {{unfixed_result_path}}/5up_0.6_paste.txt > {{unfixed_result_path}}/5up_0.6_0.25.txt
awk '$4>=5' {{unfixed_result_path}}/merge_kept_mix_ratio.txt|cut -f1 > {{unfixed_result_path}}/5up.list
perl {{resr_script_path}}/repeat_loci.pl {{unfixed_result_path}}/5up_0.6_0.25.list {{unfixed_result_path}}/5up.list > {{unfixed_result_path}}/5up_remove_loc.list
perl {{resr_script_path}}/repeatloci_filter.pl {{unfixed_result_path}}/5up_remove_loc.list {{unfixed_result_path}}/{{markkept_file}} > {{unfixed_result_path}}/{{keptfilt_file}}

#annotation of unfixed SNPs
cut -f9-11 {{unfixed_result_path}}/{{keptfilt_file}} > {{unfixed_result_path}}/{{keptsnp_file}}
#perl {{resr_script_path}}/1_MTBC_Annotation_mtbc_4411532.pl {{keptsnp_file}} > {{keptanofilt_file}}
python {{resr_script_path}}/annotate_mtb_results.py {{unfixed_result_path}}/{{keptsnp_file}} {{common_result_path}}/{{varscan_file}} {{resr_db_path}} > {{unfixed_result_path}}/{{result_file}}
"""

SCRIPT_TEMPLATE_TBPROFILER1 = """
tb-profiler profile --bam {{common_result_path}}/{{sample_id}}.sorted.bam --dir {{tbprofiler_result_path}} --prefix {{sample_id}} --csv
"""

SCRIPT_TEMPLATE_TBPROFILER2 = """
conda run -n {{conda_env}} tb-profiler profile --bam {{common_result_path}}/{{sample_id}}.sorted.bam --dir {{tbprofiler_result_path}} --prefix {{sample_id}} --csv
"""

if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description=DESCRIPTION)
    parser.add_argument('--paired_fastq_pattern', default="*_{{readnum}}.fastq", help="paired FASTQ pattern")
    parser.add_argument('--single_fastq_pattern', default="*.fastq", help="single end FASTQ pattern")
    parser.add_argument('--config', default="config.json", help="config file")
    parser.add_argument('input_path', help="input path to sample")
    parser.add_argument('result_path', help="result path")

    args = parser.parse_args()
    paired_end = True
    abs_input_path = os.path.abspath(args.input_path)
    fastq_files = find_fastq_files(abs_input_path, [args.paired_fastq_pattern])
    if len(fastq_files) == 0:
        print("no paired FASTQ files found - generating single read pipeline")
        paired_end = False
        fastq_files = find_fastq_files(abs_input_path, [args.single_fastq_pattern])
        if len(fastq_files) == 0:
            print("NO single FASTQ files found - FAILURE !!!")
            exit(1)
        else:
            fq1 = fastq_files[0][0]
            basename = os.path.basename(fq1)
            stem0 = basename.replace(".fastq", "").replace("gz", "")
    else:
        print("paired FASTQ files found - generating paired end pipeline")
        fq1, fq2 = fastq_files[0]
        basename1 = os.path.basename(fq1)
        basename2 = os.path.basename(fq2)
        stem1 = basename1.replace(".fastq", "").replace("gz", "")
        stem2 = basename2.replace(".fastq", "").replace("gz", "")

        # depending on the readnum format
        stem0 = stem1.replace("_1", "")

    with open(args.config) as infile:
        config_file = json.load(infile)
        print(config_file)

    #fasta_path = "/bwa_pipeline/reference/MTB_ancestor_reference.fasta"
    resr_script_path = "scripts"
    resr_db_path = "databases"
    #varscan_path = "/jarfiles/VarScan.v2.4.0.jar"

    common_result_path = os.path.join(args.result_path, "common")
    fixed_result_path = os.path.join(args.result_path, "fixed")
    unfixed_result_path = os.path.join(args.result_path, "unfixed")
    tbprofiler_result_path = os.path.join(args.result_path, "tbprofiler")

    # TODO: check if the BWA indexes exist, otherwise generate them with
    # bwa index -a bwtsw <fasta>
    config = {
        "sample_id": stem0,
        "fasta_path": config_file['fasta_path'],
        "resr_script_path": resr_script_path,
        "resr_db_path": resr_db_path,
        "varscan_path": config_file['varscan_path'],
        "fixed_result_path": fixed_result_path,
        "unfixed_result_path": unfixed_result_path,
        "common_result_path": common_result_path,
        "tbprofiler_result_path": tbprofiler_result_path,

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
    if "tb_profiler_env" in config_file:
        tb_profiler_templ = SCRIPT_TEMPLATE_TBPROFILER2
        config['conda_env'] = config_file['tb_profiler_env']
    else:
        tb_profiler_templ = SCRIPT_TEMPLATE_TBPROFILER1

    if paired_end:
        config.update({
            "fastq1": fq1, "fastq2": fq2,
            "trimmed1": "%s_trimmed.fq" % stem1, "trimmed2": "%s_trimmed.fq" % stem2,
            "trimmedS": "%s_trimmedS.fq" % stem0,
            "sai1": "%s.sai" % stem1, "sai2": "%s.sai" % stem2, "saiS": "%s_S.sai" % stem0
        })
        script_templ = jinja2.Template(SCRIPT_TEMPLATE_PRE + SCRIPT_TEMPLATE_PAIRED +
                                       SCRIPT_TEMPLATE_COMMON + tb_profiler_templ)
    else:
        config.update({
            "fastq": fq1, "trimmed": "%s_trimmedS.fq" % stem0,
            "sai": "%s.sai" % stem0
        })
        script_templ = jinja2.Template(SCRIPT_TEMPLATE_PRE + SCRIPT_TEMPLATE_SINGLE +
                                       SCRIPT_TEMPLATE_COMMON + tb_profiler_templ)

    script_out = script_templ.render(config)
    with open("%s_snp_calling.sh" % stem0, "w") as outfile:
        outfile.write(script_out)
