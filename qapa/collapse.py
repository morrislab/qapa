# Process pre-sorted BED-like file and collapse 3' UTRs that have similar
# 3' ends.
# BED file must be sorted by 3' coordinate, strand, followed by start
# coordinate.

import sys
import pandas as pd
import numpy as np
import logging

logger = logging.getLogger(__name__)

class Interval:
    def __init__(self, l, sp=None):
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
        self.species = self._guess_species(sp)

    def is_forward(self):
        return self.strand == "+"

    def merge(self, b):
        '''Merge this interval with another interval
        '''
        self.start = min(self.start, b.start)
        self.end = max(self.end, b.end)
        if b.name not in self.name:
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
        if self.name.startswith('ENST0'):
            return 'hsa'
        elif self.name.startswith('ENSMUST0'):
            return 'mmu'
        elif species is not None:
            return species
        return 'unk'


def overlaps(a, b, dist3):
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


def merge_bed(args, inputfile):
    '''Go through a sorted BED file and merge intervals together'''
    if args.debug:
        logger.setLevel(logging.DEBUG)

    if inputfile == '-':
        df = pd.read_table(sys.stdin)
    else:
        df = pd.read_table(inputfile)

    # Remove 3' UTRs with length == 0
    df = df[df.utr_length > 0]

    # Sort by three prime coordinate
    logger.info("Sorting data frame by 3' end")
    df['three_prime'] = np.where(df['strand'] == '+', df['end'], df['start'])
    df = df.sort_values(['strand', 'seqnames', 'three_prime'])

    prev_interval = None
    collapsed_three_prime = []
    overlap_diff_genes = set()

    logger.info("Iterating and merging intervals by 3' end")
    for index, line in df.iterrows():
        my_interval = Interval(line, args.species)

        if prev_interval is None:
            prev_interval = my_interval
        elif overlaps(prev_interval, my_interval, args.dist3):
            if same_gene(prev_interval, my_interval):
                prev_interval.merge(my_interval)
            else:
                logger.debug("Skipping overlapping but different "
                      "genes %s and %s" %
                      (prev_interval.gene_id, my_interval.gene_id))
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

    logger.info("Updating 5' end for each gene")

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

