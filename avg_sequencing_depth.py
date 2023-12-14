#!/usr/bin/python3

"""
Calculate average sequencing depth, with Python instead of SED/AWK
to make it more maintainable/readable
"""
import sys

MAGIC_NUMBER = 4411532.0
if __name__ == '__main__':
    with open(sys.argv[1]) as infile:
        infile.readline()  # skip header
        n = 0
        the_sum = 0
        for line in infile:
            comps = line.strip().split('\t')
            coverage = int(comps[4].split(':')[1])
            if coverage >= 3:
                n += 1
                the_sum += coverage
        col0 = n / MAGIC_NUMBER
        col1 = the_sum / n
        print("\t%f\t%f" % (col0, col1))
