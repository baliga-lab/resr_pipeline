#!/usr/bin/env python3

"""
Replacement script for PE_IS_filt.py. This rewrite is to make the
intention more clear.
"""
import argparse
import pandas

DESCRIPTION = """varscan_drop_excluded.py - drop excluded regions from Varscan file"""

if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description=DESCRIPTION)
    parser.add_argument('exclude_list', help="exclude list")
    parser.add_argument('varscan', help="input varscan file")
    parser.add_argument('outfile', help="output varscan file")
    args = parser.parse_args()
    excluded = set()
    with open(args.exclude_list) as infile:
        for line in infile:
            excluded.add(int(line.strip()))
    df = pandas.read_csv(args.varscan, sep='\t')
    result_indices = []
    for index, row in df.iterrows():
        if not row['Position'] in excluded:
            result_indices.append(index)
    result_df = df[df.index.isin(result_indices)]
    result_df.to_csv(args.outfile, sep='\t', index=False)
