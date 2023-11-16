#!/usr/bin/env python3

import argparse
import jinja2
import os
from globalsearch.rnaseq.find_files import find_fastq_files

DESCRIPTION = "make_fixed_snp_calling"

TEMPLATE = """#!/bin/bash
sickle pe -l 35 -f {{fastq1}} -r {{fastq2}} -t sanger -o {{trimmed1}} -p {{trimmed2}} -s {{trimmedS}}
bwa aln -R 1 {{fasta_path}} {{trimmed1}} > {{sai1}}
bwa aln -R 1 {{fasta_path}} {{trimmed2}} > {{sai2}}
bwa aln -R 1 {{fasta_path}} {{trimmedS}} > {{saiS}}
bwa sampe -a 1000 -n 1 -N 1 {{fasta_path}} {{sai1}} {{sai2}} {{trimmed1}} {{trimmed2}} > {{paired_sam}}
bwa samse -n 1 {{fasta_path}} {{saiS}} {{trimmedS}} > {{single_sam}}
samtools view -bhSt {{fasta_path}}.fai {{paired_sam}} -o {{paired_bam}}
samtools view -bhSt {{fasta_path}}.fai {{single_sam}} -o {{single_bam}}
samtools merge {{merged_bam}} {{paired_bam}} {{single_bam}}
samtools sort {{merged_bam}} -o {{sorted_bam}}
samtools index {{sorted_bam}}
samtools mpileup -q 30 -Q 20 -ABOf {{fasta_path}} {{sorted_bam}} > {{pileup_file}}
java -jar {{varscan_path}} mpileup2snp {{pileup_file}} --min-coverage 3 --min-reads2 2 --min-avg-qual 20 --min-var-freq 0.01 --min-freq-for-hom 0.9 --p-value 99e-02 --strand-filter 1 > {{varscan_file}}
java -jar {{varscan_path}} mpileup2cns {{pileup_file}} --min-coverage 3 --min-avg-qual 20 --min-var-freq 0.75 --min-reads2 2 --strand-filter 1 > {{cns_file}}
perl {{resr_path}}/PPE_PE_INS_filt.pl {{resr_path}}/PPE_INS_loci.list {{varscan_file}} > {{ppe_file}}
perl {{resr_path}}/format_trans.pl {{ppe_file}} > {{for_file}}
perl {{resr_path}}/fix_extract.pl {{for_file}} > {{fix_file}}
perl {{resr_path}}/unfix_pileup_match.pl {{for_file}} {{pileup_file}} > {{forup_file}}
cut -f2-4 {{fix_file}} > {{snp_file}}
"""

UNFIXED_TEMPLATE="""#!/bin/bash
TB_FASTA=/bwa_pipeline/reference/MTB_ancestor_reference.fasta
RESR=/resR_scripts

# parameters filename
i=Sample_R1.fastq
j=Sample_R2.fastq
fq1=Sample_tr_1.fq
fq2=Sample_tr_2.fq
fq3=Sample_tr_S.fq
sai1=Sample_R1.sai
sai2=Sample_R2.sai
sai3=Sample_R.sai
samp=Sample.piared.sam
sams=Sample.single.sam
bamp=Sample.paired.bam
bams=Sample.single.bam
bamm=Sample.merge.bam
sortbam=Sample.sort.bam
pileup=Sample.pileup
var=Sample.varscan
cns=Sample.cns
ppe=Sample.ppe
format=Sample.for
forup=Sample.forup
dep=Sample.dep
mix=Sample.mix
mixfor=Sample.mixfor
mixmark=Sample.mixmark
markkept=Sample.markkept
keptfilt=Sample.keptfilt
keptsnp=Sample.keptsnp
keptanofilt=Sample.keptanofilt

sickle pe -l 35 -f $i -r $j -t sanger -o $fq1 -p $fq2 -s $fq3
bwa aln -R 1 $TB_FASTA $fq1 > $sai1
bwa aln -R 1 $TB_FASTA $fq2 > $sai2
bwa aln -R 1 $TB_FASTA $fq3 > $sai3
bwa sampe -a 1000 -n 1 -N 1 $TB_FASTA $sai1 $sai2 $fq1 $fq2 > $samp
bwa samse -n 1 $TB_FASTA $sai3 $fq3 > $sams
samtools view -bhSt $TB_FASTA.fai $samp -o $bamp
samtools view -bhSt $TB_FASTA.fai $sams -o $bams
samtools merge $bamm $bamp $bams
samtools sort $bamm -o $sortbam
samtools index $sortbam
#generate mpileup file for each sample
samtools mpileup -q 30 -Q 20 -ABOf $TB_FASTA $sortbam > $pileup
#using varscan for SNP calling
java -jar /jarfiles/VarScan.v2.4.0.jar mpileup2snp $pileup --min-coverage 3 --min-reads2 2 --min-avg-qual 20 --min-var-freq 0.01 --min-freq-for-hom 0.9 --p-value 99e-02 --strand-filter 1 > $var
java -jar /jarfiles/VarScan.v2.4.0.jar mpileup2cns $pileup --min-coverage 3 --min-avg-qual 20 --min-var-freq 0.75 --min-reads2 2 --strand-filter 1 > $cns
#exclude regions belonging to PPE/PE and insertion sequences, and also exclude the regions that were recently marked as error-prone (Marin, Bioinformatics, 2022)
perl $RESR/PE_IS_filt.pl $RESR/Excluded_loci_mask.list $var > $ppe
perl $RESR/format_trans.pl $ppe > $format
#extract read location from mpileup file (where does a mutation allele locate on a seqeuncing read), for further filtering based on tail distribution
perl $RESR/mix_pileup_merge.pl $format $pileup > $forup
#average sequencing depth, only include samples with genome coverage rate >0.9 and sequencing depth >20X
sed 's/:/\t/g' $cns|awk '{if (\$6 >= 3){n++;sum+=\$6}} END {print \"\t\",n/4411532,\"\t\",sum/n}' > $dep
#extract unfixed SNPs from forup files, this will create two files: "markdisc" and "markkept"; the suspected false positives(such as mutations with tail region enrichment) will be moved to markdisc file
perl $RESR/mix_extract_0.95.pl $forup > $mix
perl $RESR/forup_format.pl $mix > $mixfor
perl $RESR/info_mark.pl $mixfor > $mixmark
perl $RESR/redepin_filt.pl $RESR/Excluded_loci_mask.list $dep $mixmark
#filter list of highly repeated mutations with similar mutational frequency
#for those unfixed mutations that arise >=5 times in the 50K isolates, further check their reliability based on 1) the ratio in "markkept"; 2) the distribution of the mutational frequency.
cat *mixmarkkept > all_KEPT.txt; perl $RESR/loci_freq_count.pl all_KEPT.txt > kept_repeat.txt
cat *mixmark > all_MIX.txt;perl $RESR/loci_freq_count.pl all_MIX.txt > mix_repeat.txt
perl $RESR/repeat_number_merge.pl mix_repeat.txt kept_repeat.txt > merge_kept_mix.txt
perl $RESR/ratio.pl merge_kept_mix.txt > merge_kept_mix_ratio.txt
awk '$4>=5' merge_kept_mix_ratio.txt |awk '$6>0.6'|cut -f1|while read i;do echo $i > $i.per5up.txt;grep -w $i all_KEPT.txt|cut -f12 >> $i.per5up.txt;done
paste *per5up.txt > 5up_0.6_paste.txt
perl $RESR/stdv.pl 5up_0.6_paste.txt |awk '$2<0.25'|cut -f1 > 5up_0.6_0.25.list
perl $RESR/freq_extract.pl 5up_0.6_0.25.list 5up_0.6_paste.txt > 5up_0.6_0.25.txt
awk '$4>=5' merge_kept_mix_ratio.txt|cut -f1 > 5up.list
perl $RESR/repeat_loci.pl 5up_0.6_0.25.list 5up.list > 5up_remove_loc.list
perl $RESR/unfix_scripts/repeatloci_filter.pl 5up_remove_loc.list $markkept > $keptfilt
#annotation of unfixed SNPs
cut -f9-11 $keptfilt > $keptsnp
perl $RESR/1_MTBC_Annotation_mtbc_4411532.pl $keptsnp > $keptanofilt
"""

if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description=DESCRIPTION)
    parser.add_argument('--fastq_pattern', default="*_{{readnum}}.fastq", help="FASTQ pattern")
    parser.add_argument('input_path', help="input path to sample")
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
    resr_path = "/resR_scripts"
    varscan_path = "/jarfiles/VarScan.v2.4.0.jar"

    fixed_templ = jinja2.Template(TEMPLATE)
    config = {
        "fasta_path": fasta_path, "resr_path": resr_path, "varscan_path": varscan_path,
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
        "forup_file": "%s.forup" % stem0, "snp_file": "%s.snp" % stem0
    }
    fixed_out = fixed_templ.render(config)
    with open("fixed_snp_calling.sh", "w") as outfile:
        outfile.write(fixed_out)
