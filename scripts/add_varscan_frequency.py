#!/usr/bin/env python3

import argparse
import pandas as pd

DESCRIPTION = "add_varscan_frequency.py - add varscan information"

if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description=DESCRIPTION)
    parser.add_argument('result', help="preliminary result file")
    parser.add_argument('varscan', help="varscan file")
    parser.add_argument('outfile', help="output file")
    parser.add_argument('errfile', help="error file")
    args = parser.parse_args()

    error_lines = []
    result_rows = []
    with open(args.result) as infile:
        header = infile.readline().strip().split('\t')
        line_num = 2
        for line in infile:
            comps = line.strip().split('\t')
            if line_num == 2:
                refline = comps
            num_comps = len(comps)
            if num_comps != 10:
                error_lines.append(line.strip())
            else:
                result_rows.append(comps)
            line_num += 1

    header = ["VarscanPosition", "Ref", "Alt", "CodonPos", "Type_WTAA_MutAA",
              "WTCodon_MutCodon", "Gene ID", "Name",
              "Description", "Type"]

    freq_map = {}
    varscan_df = pd.read_csv(args.varscan, sep='\t')
    #print(varscan_df)
    for index, row in varscan_df.iterrows():
        chrom = row["Chrom"]
        pos = row["Position"]
        special = row["Cons:Cov:Reads1:Reads2:Freq:P-value.1"]
        freq_map[pos] = special.split(":")[4]

    frequencies = []
    for row in result_rows:
        position = int(row[0])  # at index 0
        freq = freq_map[position]
        row.append(freq)

    with open(args.outfile, 'w') as outfile:
        outfile.write("%s\n" % '\t'.join(header))
        for row in result_rows:
            outfile.write("%s\n" % '\t'.join(row))

    with open(args.errfile, 'w') as errfile:
        for line in error_lines:
            errfile.write("%s\n" % line)
