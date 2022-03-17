# Update the 3' ends of GENCODE annotations with more
# complete poly(A) site coordinates from PolyAsite database and GENCODE poly(A)
# site track

import sys
import os
import pybedtools
from pybedtools import featurefuncs
import re
import logging

logger = logging.getLogger(__name__)


def extend_feature(feature, length=24):
    """Extend the 3' end by length
    """
    if feature.strand == "+":
        feature.end = feature.end + length
    else:
        if feature.start >= length:
            feature.start = feature.start - length
        else:
            feature.start = 0
    return feature


def restore_feature(feature, length=24):
    """Restore the original 3' end coordinate
    """
    if feature.strand == "+":
        feature.end = feature.end - length
    else:
        if feature.start == 0:
            exonStarts = feature[9].split(",")
            feature.start = int(exonStarts[0])
        else:
            feature.start = feature.start + length
    return feature


def update_3prime(feature, min_distance=24, min_intermediate_pas=4, custom=False):
    """
    Use the overlapping poly(A) site features and update the 3' end coordinate
    of the interval.
    """
    # 0-based indices!
    site_start = 12
    site_end = site_start + 1
    site_name = site_end + 1
    site_numsamples = site_name + 1
    start2 = 6
    end2 = start2 + 1

    if feature[site_start] == '-1':
        return feature

    # Update cluster poly(A) site coordinates if match is from PolyAsite
    match = re.match(r'(chr)?.*:(\d+):.*', feature[site_name])
    if match and not custom:
        if feature.strand == "+":
            feature[site_end] = int(match.group(2))
            feature[site_start] = int(feature[site_end]) - 1
        else:
            feature[site_start] = int(match.group(2))
            feature[site_end] = int(feature[site_start]) + 1

    # Update 3' end coordinate
    if feature.strand == "+":
        dist_from_three_prime = abs(int(feature[site_end]) - feature.end)
        feature.end = int(feature[site_end])
        # Remove feature if match is outside of 3' UTR
        if feature.end < int(feature[end2]):
            return None
        feature.score = feature.end - int(feature[end2])
    else:
        dist_from_three_prime = abs(int(feature[site_start]) - feature.start)
        feature.start = int(feature[site_start])
        if feature.start > int(feature[start2]):
            return None
        feature.score = int(feature[start2]) - feature.start

    if dist_from_three_prime > min_distance and \
            not custom and \
            int(feature[site_numsamples]) < min_intermediate_pas:
        #logger.debug("Skipping {}. Distance from 3' end: {}"\
                     #.format(feature[0:6], dist_from_three_prime))
        return None

    return feature


def resolve_overlaps(feature):
    """
    For overlapping features, choose the shorter isoform
    """
    names = feature.name.split("|")
    if len(names) > 1:
        lengths = [int(x) for x in feature.score.split("|")]
        idx = lengths.index(min(lengths))
        feature.score = lengths[idx]
        feature.name = names[idx]
        feature[9] = feature[9].split("|")[idx]
        feature[10] = feature[10].split("|")[idx]

    return feature


def _add_chr(feature):
    if feature.chrom.startswith("chr"):
        return feature
    feature.chrom = 'chr' + feature.chrom
    return feature

def _move_num_exp(feature):
    # Move number of experiments field to score field
    feature.score = feature[7]
    return feature

def sort_bed(bedobj):
    """
    Use GNU sort
    """
    tmpbed = bedobj._tmp()
    os.system('sort {0} -k1,1 -k2,2n -k3,3n > {1}'.format(bedobj.fn, tmpbed))
    return pybedtools.BedTool(tmpbed)


def preprocess_gencode_polya(gencode_polya_file):
    """
    Preprocess GENCODE polyA track
    """
    logger.info("Preprocessing %s" % gencode_polya_file)
    gencode = pybedtools.BedTool(gencode_polya_file)\
        .filter(lambda x: x.name == 'polyA_site')\
        .saveas()
    gencode = sort_bed(gencode)
    validate(gencode, gencode_polya_file)
    return gencode


def preprocess_polyasite(polyasite_file, min_polyasite):
    """
    Pre-proecess polyasite file
    """
    logger.info("Preprocessing %s" % polyasite_file)
    # add an empty filter step as hack to read gzipped files
    polyasite = pybedtools.BedTool(polyasite_file)\
        .filter(lambda x: x)\
        .saveas()
    polyasite = sort_bed(polyasite)
    validate(polyasite, polyasite_file)

    pas_filter = re.compile("(DS|TE)$")
    is_v2 = polyasite.field_count() == 11

    if not is_v2:
        logger.info("Detected PolyASite version 1")
        field_index = 3
    else:
        logger.info("Detected PolyASite version 2")
        field_index = 9
        polyasite = polyasite.each(_add_chr)\
                             .each(_move_num_exp)\
                             .saveas()

    polyasite_te = polyasite\
        .filter(lambda x: int(x[4]) >= min_polyasite)\
        .filter(lambda x: pas_filter.search(x[field_index]))\
        .cut(range(0, 6))\
        .saveas()
    validate(polyasite_te, polyasite_file)
    return polyasite_te


def validate(bedobj, filename):
    """
    Validate the input BED file using solution from
    https://github.com/daler/pybedtools/issues/252
    """
    try:
        ft = bedobj.file_type
        if ft == 'empty':
            raise IOError(f"BED file {filename} is empty!")
    except IndexError as err:
        raise err(f"Error reading the BED file {filename}. Is the file properly "
                  "formatted?")


def main(args, input_filename, fout=sys.stdout):
    if args.debug:
        logger.setLevel(logging.DEBUG)

    # Load intervals
    fin = sys.stdin if input_filename == '-' else input_filename
    utrs = pybedtools.BedTool(fin).each(extend_feature).saveas()
    custom_mode = False

    # Load databases with pybedtools
    if args.other:
        logger.info("Annotating with %s" % args.other)
        custom_mode = True
        custom = pybedtools.BedTool(args.other)
        validate(custom, args.other)
        sites = sort_bed(custom)
    elif not args.no_annotation:
        gencode = preprocess_gencode_polya(args.gencode_polya)

        polyasite_te = preprocess_polyasite(args.polyasite, args.min_polyasite)

        sites = gencode.cat(polyasite_te, postmerge=False)

        # Downstream 1kb PAS
        # pas_filter = re.compile("DS$")
        # polyasite_ds = polyasite\
        #                 .filter(lambda x: pas_filter.search(x.name))\
        #                 .saveas()

    # Procedure:
    #   - Intersect with databases
    #   - Update 3' coordinate with overlapping poly(A) feature
    #   - Sort
    #   - Group by features and collapse feature name and scores
    #   - Restore BED format with cut()
    #   - Use custom function to resolve overlapping features from groupby

    if args.no_annotation:
        overlap_utrs = utrs.each(restore_feature)\
                           .saveas()
        logger.info("Skipping annotation step")
    else:
        logger.info("Finding intersect")
        overlap_utrs = utrs.intersect(sites, s=True, wa=True, wb=True)\
                           .each(restore_feature)\
                           .each(update_3prime,
                                 min_intermediate_pas=args.intermediate_polyasite,
                                 custom=custom_mode)\
                           .saveas()
        logger.debug("Annotated %d UTRs", len(overlap_utrs))

    if len(overlap_utrs) == 0:
        raise RuntimeError("Failed to find any overlap between UTRs and "
                "annotation. Please check your input files.")


    overlap_utrs = sort_bed(overlap_utrs)\
        .groupby(g=[1, 2, 3, 6, 7, 8, 9], c=[4, 5, 10, 11],
                 o=['collapse'] * 4, delim="|")\
        .cut([0, 1, 2, 7, 8, 3, 4, 5, 6, 9, 10])\
        .each(resolve_overlaps)

    # overlap_utrs = utrs.intersect(sites, s=True, loj=True, stream=True)
    # annotated_utrs = overlap_utrs.filter(lambda x: x[12] != "-1").saveas()
    # unannotated_utrs = overlap_utrs.filter(lambda x: x[12] == "-1").saveas()

    # annotated_utrs = annotated_utrs.each(update_3prime)

    # annotated_utrs = pybedtools.BedTool(tmpbed)\
    #                    .groupby(g=[1, 2, 3, 6, 7, 8, 9],
    #                             c=[4, 5, 10, 11],
    #                             o=['collapse'] * 4, delim="|")\
    #                    .cut([0, 1, 2, 7, 8, 3, 4, 5, 6, 9, 10])\
    #                    .each(resolve_overlaps)\
    #                    .saveas()

    # Now look for downstream sites (DS) for annotated UTRs (add new feature if
    # found)
    #   find closest
    #   update 3prime
    #   cut events
    # max_ds_distance = 1100
    # ds_annotated_utrs = annotated_utrs.closest(polyasite_ds, s=True, D='a',
    #                                            fd=True)\
    #                         .filter(lambda x: 0 < x[19] < max_ds_distance)\
    #                         .each(update_3prime)\
    #                         .sort()\
    #                         .groupby(g=[1,2,3,6,7,8,9], c=[4,6,10,11])
    #                         .saveas()

    # Repeat for unannotated UTRs (update if found)
    # ds_unannotated_utrs = overlap_utrs.closest(polyasite_ds, s=True, D='a', fd=True)\
    #                 .filter(lambda x: 0 < x[19] < max_ds_distance)\
    #                 .each(update_3prime)\
    #                 .saveas()

    # Concat
    # final_utrs = annotated_utrs.cat(*[ds_annotated_utrs, ds_unannotated_utrs],
    #                                 postmerge=False)

    header = ["seqnames", "start", "end", "name", "utr_length", "strand",
              "lastexon_cds_start", "lastexon_cds_end", "name2",
              "exonStarts", "exonEnds"]
    fout.write("\t".join(header) + "\n")
    fout.write(str(overlap_utrs))

