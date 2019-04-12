# Process pre-sorted BED-like file and collapse 3' UTRs that have similar
# 3' ends.
# BED file must be sorted by 3' coordinate, strand, followed by start
# coordinate.

from __future__ import print_function
import sys
import re
import pandas as pd
import numpy as np
import warnings
import logging

logging.basicConfig(level=logging.INFO, stream=sys.stderr,
                    format='%(asctime)s - %(levelname)-8s - %(message)s')
logger = logging.getLogger('collapse')
logging.captureWarnings(True)

_TAG = 'collapse' 

class Interval:
    def __init__(self, l, sp=None):
        #l = line.rstrip().split("\t")
        self.chrom = l[0]
        self.start = int(l[1])
        self.end = int(l[2])
        self.name = l[3]
        self.score = int(l[4])
        self.strand = l[5]
        # Last exon coordinates
        self.start2 = int(l[6])
        self.end2 = int(l[7])
        self.gene_id = l[8]
        #self.utr_id = l[9]
        self.species = self._guess_species(sp)

    def is_forward(self):
        return self.strand == "+"

    def merge(self, b):
        '''Merge this interval with another interval
        '''
        self.start = min(self.start, b.start)
        self.end = max(self.end, b.end)
        if not re.search(b.name, self.name):
            self.name = self.name + "," + b.name
        self.start2 = min(self.start2, b.start2)
        self.end2 = max(self.end2, b.end2)

    def set_score(self):
        '''Set the score column'''
        if self.is_forward():
            self.score = self.end - self.end2
        else:
            self.score = self.start2 - self.start

    def finalize(self, species=None):
        '''Print the interval using the highest peak'''
        if self.is_forward():
            utr_co = [self.end2, self.end]
        else:
            utr_co = [self.start, self.start2]
        new_name = [self.name, self.species, self.chrom, 
                    self.start, self.end, self.strand, 'utr'] + utr_co
        new_name = "_".join([str(x) for x in new_name])
        self.set_score()
        return [self.chrom, self.start, self.end, new_name, self.score,
                self.strand, self.gene_id]

    def _guess_species(self, species=None):
        if re.match(r'ENST\d+', self.name):
            return 'hg19'
        elif re.match(r'ENSMUST\d+', self.name):
            return 'mm10'
        elif species is not None:
            return species
        warnings.warn('Could not guess species from gene name!' +
            ' To disable this warning, use --species option', Warning)
        return 'unk'


def overlaps(a, b, dist5, dist3):
    '''Determine if two intervals have close 3' ends'''

    if a.chrom == b.chrom and a.strand == b.strand:
        if a.is_forward() and abs(b.end - a.end) <= dist3:
            return True
        elif not a.is_forward() and abs(a.start - b.start) <= dist3:
            return True
    return False


def same_gene(a, b):
    '''Get Ensembl Gene ID and check if they are the same'''
    return a.gene_id == b.gene_id


# def getoptions():
#     usage = "usage: python %prog [options] sorted_grouped_reads.bed annotation.db"
#     desc = "Process sorted BED file and collapse consecutive intervals that" + \
#     " have similar 3' ends."
#     parser = OptionParser(usage=usage, description=desc)
#     parser.add_option("-d", type = "int", default = 24,
#             dest = "dist3", metavar = "DISTANCE",
#             help = "Maximum distance between 3' ends to merge [%default]")
#     parser.add_option("-f", type = "int", default = 3,
#             dest = "dist5", metavar = "DISTANCE",
#             help = "Maximum distance between 5' ends to merge [%default]")
#     (opts, args) = parser.parse_args()
#     if len(args) < 1:
#         parser.error("Missing arguments")
#     return (opts, args)


def merge_bed(args, inputfile):
    '''Go through a sorted BED file and merge intervals together'''

    if inputfile == '-':
        df = pd.read_table(sys.stdin)
    else:
        df = pd.read_table(inputfile)

    # Remove 3' UTRs with length == 0
    df = df[df.utr_length > 0]

    # Sort by three prime coordinate
    logger.info("Sorting data frame by 3' end", tag=_TAG)
    df['three_prime'] = np.where(df['strand'] == '+', df['end'], df['start'])
    df = df.sort_values(['strand', 'seqnames', 'three_prime'])

    prev_interval = None
    collapsed_three_prime = []
    overlap_diff_genes = set()

    logger.info("Iterating and merging intervals by 3' end", tag=_TAG)
    for index, line in df.iterrows():
        my_interval = Interval(line, args.species)

        if prev_interval is None:
            prev_interval = my_interval
        elif overlaps(prev_interval, my_interval, args.dist5, args.dist3):
            if same_gene(prev_interval, my_interval):
                prev_interval.merge(my_interval)
            else:
                # print("Skipping overlapping but different "
                #       "genes %s and %s" %
                #       (prev_interval.gene_id, my_interval.gene_id),
                #       file=sys.stderr)
                overlap_diff_genes.add(prev_interval.gene_id)
                overlap_diff_genes.add(my_interval.gene_id)
                prev_interval = None
        else:
            # finalize current interval
            collapsed_three_prime.append(prev_interval.finalize())
            # start a new one
            prev_interval = my_interval

    collapsed_three_prime.append(prev_interval.finalize())

    # After collapsing 3' ends, update 5' ends so that they match
    three_prime_df = pd.DataFrame(collapsed_three_prime,
                                  columns=['chr', 'start', 'end', 'name',
                                           'score', 'strand', 'gene_id'])

    three_prime_df = \
        three_prime_df[~three_prime_df['gene_id'].isin(overlap_diff_genes)]

    logger.info("Updating 5' end for each gene", tag=_TAG)

    # Filter by forward and reverse strand
    forward = three_prime_df[three_prime_df.strand == '+']
    five_prime_pos = forward.groupby('gene_id')['start'].min()
    forward = forward.join(five_prime_pos, on='gene_id', rsuffix='_r')
    forward['start'] = forward['start_r']
    forward = forward.drop('start_r', 1)

    reverse = three_prime_df[three_prime_df.strand == '-']
    five_prime_pos = reverse.groupby('gene_id')['end'].max()
    reverse = reverse.join(five_prime_pos, on='gene_id', rsuffix='_r')
    reverse['end'] = reverse['end_r']
    reverse = reverse.drop('end_r', 1)

    assert three_prime_df.shape[0] == forward.shape[0] + reverse.shape[0]

    # Join back together
    return pd.concat([forward, reverse])


if __name__ == '__main__':
    pass
