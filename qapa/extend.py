# For each gene, find the common 5' end shared by all last exons. In most cases,
# the last exons already share a 5' end and nothing needs to be done. For
# other cases, such as when there is a last exon overlapping an internal exon:
# e.g.
#
# 1)  3'----------------[    ]     [    ]
# 2)          3'--------[    ]     [    ]
# 3)                           3'--[    ]
#
# We need to extend the 5' ends of the #1 and #2 to match that of #3, but also
# need to take into account the intron. We can use BED12 format to keep track of
# the "blocks", which will need later to extract the proper intronless sequence.

import sys
import pandas as pd
import logging

logger = logging.getLogger(__name__)


def is_plus(strands):
    return all([x == "+" for x in strands])


def _split_int(string):
    return [int(x) for x in string.split(",")]


def _extend(feature, most_five_prime, forward, num_extends=0):
    '''
    For each feature, check if the 5' end of the 3' UTR interval is in front
    of the most proximal 3' site.

    If yes, then we want to extend it by including additional upstream exons.
    If no, then just create the block* columns.
    '''
    stack = []
    block_sizes = []
    block_starts = []
    exon_starts = _split_int(feature.exonStarts)
    exon_ends = _split_int(feature.exonEnds)

    """
    Traverse exon list and record exons appear before most_five_prime position

    Use a stack to keep track of exons being considered.
    Start with the 3' most exon. If the next exon still succeeds the 5'-most
    coordinate, add this exon and move on to the next.

    The number of times you can move on is controlled by num_extends or the
    number of exons, whichever comes first.

    Then, iterate through the stack to come up with the block* fields.
    """
    if forward:
        min_index = max(0, len(exon_starts) - 1 - num_extends)

        # base case
        i = len(exon_starts) - 1
        stack.append(i)
        i -= 1

        if most_five_prime >= exon_starts[min_index]:
            while i >= min_index \
                    and exon_starts[i] >= most_five_prime:
                stack.append(i)
                i -= 1

        stack.reverse()
        for j in stack:
            block_sizes.append(exon_ends[j] - exon_starts[j])
            block_starts.append(exon_starts[j] - exon_starts[stack[0]])

        interval_start = exon_starts[stack[0]]
        interval_end = feature.end
    else:
        min_index = min(len(exon_ends) - 1, num_extends)
        # base case
        stack.append(0)
        i = 1

        if most_five_prime <= exon_ends[min_index]:
            while i <= min_index \
                    and exon_ends[i] <= most_five_prime:
                stack.append(i)
                i += 1

        for j in stack:
            block_sizes.append(exon_ends[j] - exon_starts[j])
            block_starts.append(exon_starts[j] - exon_starts[0])

        interval_start = feature.start
        interval_end = exon_ends[stack[-1]]

    # Create new block columns
    block_cols = pd.Series([len(block_sizes),
                           ",".join([str(x) for x in block_sizes]),
                           ",".join([str(x) for x in block_starts])],
                           index=['blockCount', 'blockSizes', 'blockStarts'])
    newfeature = pd.concat([feature, block_cols])
    newfeature['start'] = interval_start
    newfeature['end'] = interval_end
    return newfeature


def extend_5prime(feature_group, numextends=0):
    '''
    For a set of features grouped by gene_id, find the most 5' end and extend
    the 5' ends of all features to match it
    '''

    newgroup = []
    forward = is_plus(feature_group.strand)
    if forward:
        most_five_prime = min(feature_group.start)
    else:
        most_five_prime = max(feature_group.end)

    for index, row in feature_group.iterrows():
        newgroup.append(_extend(row, most_five_prime, forward,
                                numextends))

    newgroup = pd.DataFrame(newgroup)\
                 .drop(['exonStarts', 'exonEnds'], axis=1)

    if forward:
        newgroup = newgroup.groupby('name2')\
                        .apply(lambda g: g[g['end'] >= g['start'].max()])
    else:
        newgroup = newgroup.groupby('name2')\
                        .apply(lambda g: g[g['start'] <= g['end'].min()])
    return newgroup


def main(args, input_filename):
    if args.debug:
        logger.setLevel(logging.DEBUG)

    if input_filename == '-':
        df = pd.read_table(sys.stdin)
    else:
        df = pd.read_table(input_filename)

    # Process each group
    newdf = []

    from concurrent.futures import ProcessPoolExecutor, as_completed

    with ProcessPoolExecutor(args.cores) as executor:
        futures = {executor.submit(extend_5prime, group, args.numextends): name2 \
                    for name2, group in df.groupby('name2')}
        for future in as_completed(futures):
            name2 = futures[future]
            try:
                newdf.append(future.result())
            except Exception as exc:
                logger.exception("Error extending %s: %s" % (name2, exc))

    # for name2, group in tqdm(df.groupby('name2'), desc="extend"):
    #     newdf.append(extend_5prime(group, args.numextends))

    return pd.concat(newdf).sort_values(['seqnames', 'start'])

