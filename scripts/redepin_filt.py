#!/usr/bin/env python3

"""
Reimplementation of redepin_filt.pl for better understanding what is actually
happening
"""
import argparse
import pandas

DESCRIPTION = "redepin_filt.py - repeat loci filter"

if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description=DESCRIPTION)
    parser.add_argument('excluded_loci', help="excluded loci file")
    parser.add_argument('seq_depth', help="sequencing depth file")
    parser.add_argument('mixmark', help="mixmark file")
    args = parser.parse_args()

    mixmarkkept_path = args.mixmark + "kept"
    mixmarkdisc_path = args.mixmark + "disc"

    excluded_loci = set()
    with open(args.excluded_loci) as infile:
        for line in infile:
            excluded_loci.add(line.strip())
    with open(args.seq_depth) as infile:
        ratio, avg_seq_depth = infile.readline().strip().split('\t')
        ratio = float(ratio)
        depth = float(avg_seq_depth)
        depth_upper = depth * 1.5
        depth_lower = depth * 0.5
    mixmark_df = pandas.read_csv(args.mixmark, sep='\t', header=None)
    mixmark_kept_df = pandas.DataFrame()  # kept
    mixmark_disc_df = pandas.DataFrame()  # discarded
    for index, row in mixmark_df.iterrows():
        loci = row[8]
        seq_depth = row[12]
        avg_reads_loc = row[0]
        if (not loci in excluded_loci and
            seq_depth > depth_lower and seq_depth < depth_upper and
            avg_reads_loc > 0.25 and avg_reads_loc < 0.75 and
            "0.1" in row[4]):
            b0, b1 = map(int, row[6].split(":"))
            row_df = pandas.DataFrame([row])
            if b1 < 50 and row[5] < 0.5:
                mixmark_kept_df = pandas.concat([mixmark_kept_df, row_df], axis=0, ignore_index=True)
            elif b1 > 50:
                mixmark_kept_df = pandas.concat([mixmark_kept_df, row_df], axis=0, ignore_index=True)
            else:
                row_df.insert(0, '', 'Filter')
                mixmark_disc_df = pandas.concat([mixmark_disc_df, row_df], axis=0, ignore_index=True)
    mixmark_kept_df.to_csv(mixmarkkept_path, sep='\t', header=None, index=False)
    mixmark_disc_df.to_csv(mixmarkdisc_path, sep='\t', header=None, index=False)
