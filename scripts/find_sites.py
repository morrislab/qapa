#!/usr/bin/env python
#
# Find overlapping regions of BED file and report position of the highest peak.
# Input BED file must be sorted first by strand, chr, start, end
# 
#import pdb
import os.path
import sys
import fileinput
import argparse
import pybedtools
from pybedtools import featurefuncs

class Interval:
    def __init__(self, line):
        l = line.rstrip().split("\t")
        self.chr = l[0]
        self.start = int(l[1])
        self.end = int(l[2])
        self.name = l[3]
        self.score = int(l[4])
        self.strand = l[5]
        self.peak =  self.end if self.strand == "+" else self.start
        self.read_count = self.score

    def merge(self, b):
        '''Merge this interval with another interval
        '''
        self.start = min(self.start, b.start)
        self.end = max(self.end, b.end)
        self.read_count += b.score

        if self.strand == "+" and b.score >= self.score:
            self.peak = b.end
            self.score = b.score
        elif self.strand == "-" and b.score > self.score:
            self.peak = b.start
            self.score = b.score

    def finalize(self):
        '''Print the interval using the highest peak'''
        (finalstart, finalend) = (None, None)
        if self.strand == "+":
            finalstart = self.peak - 1
            finalend = self.peak
        else:
            finalstart = self.peak
            finalend = self.peak + 1
        self.name = self.chr + ":" + str(self.start) + "-" + str(self.end)
        print("\t".join([self.chr, str(finalstart), str(finalend), 
            self.name, str(self.score), self.strand]))

    def finalize_interval(self):
        print("\t".join([self.chr, str(self.start), str(self.end), 
            self.peak, str(self.read_count), self.strand]))

def overlaps(a, b):
    '''Determine if two intervals overlap in genomic coordinates'''
    return a.chr == b.chr \
        and a.strand == b.strand \
        and a.start <= b.end \
        and a.end >= b.start

def getoptions():
    desc = "Process BED file of read alignments to find poly(A) site clusters" \
           " and report position and count with the highest coverage"
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('bedfile', nargs=1, 
                        help='reads as BED file')
    args = parser.parse_args()
    return args

def sort_bed(bedobj):
    # sort by strand first
    tmpbed = bedobj._tmp()
    os.system('sort {0} -k 6,6 -k1,1 -k2,2n -k3,3n > {1}'.format(bedobj.fn, tmpbed))
    return pybedtools.BedTool(tmpbed)

def sort_by_strand(bedobj, strand):
    tmpbed = bedobj._tmp()
    if strand == "+":
        os.system('sort {0} -k1,1 -k3,3n > {1}'.format(bedobj.fn, tmpbed))
    else:
        os.system('sort {0} -k1,1 -k2,2n > {1}'.format(bedobj.fn, tmpbed))
    return pybedtools.BedTool(tmpbed)


def group_reads(bedfile):
    reads = pybedtools.BedTool(bedfile)
   
    forward = reads.filter(lambda x: x.strand == "+").saveas()
    if len(forward) > 0:
        forward = sort_by_strand(forward, strand="+")\
                .groupby(g=[1,3,6], c=[2,4], o=['min','count'], full=True)\
                .cut([0,6,2,3,7,5])\
                .saveas()

    reverse = reads.filter(lambda x: x.strand == "-").saveas()
    if len(reverse) > 0:
        reverse = sort_by_strand(reverse, strand="-")\
                .groupby(g=[1,2,6], c=[3,4], o=['max','count'], full=True)\
                .cut([0,1,6,3,7,5])\
                .saveas()
 
    grouped_reads = forward.cat(reverse, postmerge=False)
    grouped_reads = sort_bed(grouped_reads).saveas()

    return grouped_reads

 
def merge_bed(bedobj):
    curInterval = None
    
    for line in bedobj:
        # use our own custom Interval class
        myInterval = Interval(str(line))

        if curInterval is None:
            curInterval = myInterval
        elif overlaps(curInterval, myInterval):
            curInterval.merge(myInterval)
        else:
            # finalize current interval
            curInterval.finalize()
            # start a new one
            curInterval = myInterval

    curInterval.finalize()

if __name__ == '__main__':
    args = getoptions()
    pybedtools.set_tempdir(os.path.dirname(os.path.abspath(args.bedfile[0])))
    pybedtools.cleanup(remove_all=False)
    # group redundant alignments
    bedobj = group_reads(args.bedfile[0])
    # merge
    merge_bed(bedobj)
