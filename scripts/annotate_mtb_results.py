#!/usr/bin/env python3

import pandas
import argparse
import os

MTB_GENOME_SIZE = 4411532

def read_databases(dbpath):
    # 1. read Tuberculist information
    header = ['GeneID', 'GeneName', 'Strand', 'Start', 'Stop', 'Description', 'Category']
    tuberculist_df = pandas.read_csv(os.path.join(dbpath, '2_Tuberculist_new_20150307'),
                                     sep='\t',
                                     header=None, names=header)
    gene_infos = {}
    gene_ids = []

    igr = {}
    igr_keys = []

    for index, row in tuberculist_df.iterrows():
        gene_ids.append(row['GeneID'])
        gene_infos[row['GeneID']] = {
            "name": row["GeneName"],
            "strand": row["Strand"],
            "start": row["Start"],
            "stop": row["Stop"],
            "description": row["Description"],
            "category": row["Category"]
        }
    # we have the gene_infos map, now build igr
    for index, row in tuberculist_df.iterrows():
        if index > 0:
            if row["Start"] > tuberculist_df.iloc[index - 1]["Stop"]:
                prev_gene = tuberculist_df.iloc[index - 1]["GeneID"]
                igr_key = "%s-%s" % (prev_gene, row["GeneID"])
                igr[igr_key] = {
                    "start": gene_infos[prev_gene]["stop"] + 1,
                    "stop": row['Start'] - 1
                }
                igr_keys.append(igr_key)

    # 2. read Codon -> AA mapping
    code = {}
    genetic_codes_df = pandas.read_csv(os.path.join(dbpath, '3_genetic_codes'),
                                       sep='\t',
                                       header=None, names=['Codon', 'Result'])
    for index, row in genetic_codes_df.iterrows():
        code[row['Codon']] = row['Result']

    # 3. Read genome
    genome = ""
    with open(os.path.join(dbpath, '4_mtbc_sequence.fasta')) as infile:
        for line in infile:
            if not line.startswith(">"):
                genome += line.strip()

    return gene_ids, gene_infos, code, genome, igr_keys, igr

DESCRIPTION = """annotate_mtb_results - Annotate SNP results with Tuberculist information"""

def is_in_gene(pos, gene_info):
    return pos >= gene_info['start'] and pos <= gene_info['stop']


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description=DESCRIPTION)
    parser.add_argument('snpfile', help="SNP result file")
    parser.add_argument('varscan', help="VarScan result file")
    parser.add_argument('dbpath', help="path for database files")
    args = parser.parse_args()

    special_map = {}
    varscan_df = pandas.read_csv(args.varscan, sep='\t')
    for index, row in varscan_df.iterrows():
        chrom = row["Chrom"]
        pos = row["Position"]
        cons, cov, reads1, reads2, freq, pval = row["Cons:Cov:Reads1:Reads2:Freq:P-value.1"].split(":")
        special_map[pos] = {
            "cons": cons,
            "cov": cov,
            "reads1": reads1,
            "reads2": reads2,
            "freq": freq,
            "pval": pval
        }

    gene_ids, gene_infos, code, genome, igr_keys, igr = read_databases(args.dbpath)
    snp_df = pandas.read_csv(args.snpfile, sep='\t', header=None,
                             names=["VarscanPos", "Ref", "Alt"])
    header = ["VarscanPosition", "Ref", "Alt", "CodonPos", "Type_WTAA_MutAA",
              "WTCodon_MutCodon", "Gene ID", "Name",
              "Description", "Type", "Cons", "Cov", "Reads1", "Reads2", "Freq", "Pval"]
    print('\t'.join(header))
    for index, row in snp_df.iterrows():
        special = special_map[row['VarscanPos']]
        for gene_id in gene_ids:
            gene_info = gene_infos[gene_id]
            if is_in_gene(row['VarscanPos'], gene_info):
                if gene_id.startswith('MTB'):
                    print("%s\t%s\t%s\t-\t---\t---\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" %
                          (row["VarscanPos"], row["Ref"], row["Alt"],
                           gene_id,
                           gene_info['name'], gene_info['description'], gene_info['category'],
                           # varscan special
                           special["cons"],
                           special["cov"],
                           special["reads1"],
                           special["reads2"],
                           special["freq"],
                           special["pval"]))
                else:
                    if gene_info['strand'] == '+':
                        #print("NON-MTB GENE (+) %s\n" % gene_id);
                        length = gene_info['stop'] - gene_info['start'] + 1
                        seq = genome[gene_info['start'] - 1:(gene_info['start'] -  1) + length]
                        loci = row["VarscanPos"] - gene_info['start'] + 1
                        ct = int(loci / 3)
                        remain = loci % 3
                        if remain == 0:
                            codon = ct
                            wildtype = seq[loci - 3: (loci - 3) + 3]
                            mutation = wildtype[:2] + row['Alt']
                        elif remain == 1:
                            codon = ct + 1
                            wildtype = seq[loci - 1: (loci - 1) + 3]
                            mutation = row['Alt'] + wildtype[1:3]
                        elif remain == 2:
                            codon = ct + 1
                            wildtype = seq[loci - 2: (loci - 2) + 3]
                            mutation = wildtype[0] + row['Alt'] + wildtype[2]
                    elif gene_info['strand'] == '-':
                        #print("NON-MTB GENE (-): %s\n" % gene_id);
                        length = gene_info['stop'] - gene_info['start'] + 1
                        seq = genome[gene_info['start'] - 1:(gene_info['start'] -  1) + length]
                        seq = seq[::-1]  # reverse the sequence
                        loci = gene_info['stop'] - row['VarscanPos'] + 1
                        ct = int(loci / 3)
                        remain = loci % 3
                        if remain == 0:
                            codon = ct
                            wildtype = seq[loci - 3: (loci - 3) + 3]
                            mutation = wildtype[0:2] + row['Alt']
                        elif remain == 1:
                            codon = ct + 1
                            wildtype = seq[loci - 1: (loci - 1) + 3]
                            mutation = row['Alt'] + wildtype[1:3]
                        elif remain == 2:
                            codon = ct + 1
                            wildtype = seq[loci - 2: (loci - 2) + 3]
                            mutation = wildtype[0] + row['Alt'] + wildtype[2]

                        # translate bases
                        wildtype = wildtype.translate(str.maketrans("ATGC", "TACG"))
                        mutation = mutation.translate(str.maketrans("ATGC", "TACG"))

                    if code[wildtype] == code[mutation]:
                        restype = "Synonymous"
                    else:
                        restype = "Nonsynonymous"

                    code_info = '%s-%s-%s' % (restype, code[wildtype], code[mutation])
                    triplets = '%s-%s' % (wildtype, mutation)
                    print("%s\t%s\t%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" %
                          (row["VarscanPos"], row["Ref"], row["Alt"],
                           codon, code_info, triplets, gene_id, gene_info['name'],
                           gene_info['description'], gene_info['category'],
                           # varscan special
                           special["cons"],
                           special["cov"],
                           special["reads1"],
                           special["reads2"],
                           special["freq"],
                           special["pval"]))

            else:
                #print("not a gene location")
                pass

        # II. Iterate over joined genes
        for j in igr_keys:
            j_info = igr[j]
            if row['VarscanPos'] >= j_info['start'] and row['VarscanPos'] <= j_info['stop']:
                gene1, gene2 = j.split('-')
                if gene1 == 'Rv3924c' and gene2 == 'Rv0001':
                    left = row['VarscanPos'] - j_info['start']
                    right = MTB_GENOME_SIZE - row['VarscanPos'];
                else:
                    left = row['VarscanPos'] - gene_infos[gene1]['stop']
                    right= gene_infos[gene2]['start'] - row['VarscanPos']
                gene_info1 = gene_infos[gene1]
                gene_info2 = gene_infos[gene2]
                loc_str = "%s%d-%d%s" % (gene_info1['strand'], left, right, gene_info2['strand'])
                name_str = "%s-%s" % (gene_info1['name'], gene_info2['name'])
                desc_str = "%s##%s" % (gene_info1['description'], gene_info2['description'])
                cat_str = "%s##%s" % (gene_info1['category'], gene_info2['category'])
                print("%d\t%s\t%s\t-\t---\t%s\t%s\t%s\t%s\t%s\t%s" %
                      (row['VarscanPos'], row['Ref'], row['Alt'],
                       loc_str, j, name_str, desc_str, cat_str, freq))
