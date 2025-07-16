import sys
from Bio import SeqIO
import pybedtools
import logging

logger = logging.getLogger(__name__)


def get_sequences(bed_file, genome):
    bed = pybedtools.BedTool(bed_file)
    logger.info("Extract sequences from %s" % genome)
    return bed.sequence(genome, s=True, name=True, fullHeader=False,
                        split=True)


def filter_sequences(fasta_file, min_length=100, fout=sys.stdout):
    seqs = set()
    ids = set()
    skipped = 0
    handle = open(fasta_file, 'r')

    for record in SeqIO.parse(handle, "fasta"):
        sequence = str(record.seq)

        if len(sequence) > min_length and \
                not (sequence in seqs or record.id in ids):
            seqs.add(sequence)
            ids.add(record.id)
            SeqIO.write(record, fout, "fasta")
        else:
            logger.debug("Skipping %s" % record.id)
            skipped += 1

    handle.close()

    if skipped > 0:
        logger.info("Skipped %d sequences (too short or is a duplicate)" % skipped)


def add_decoys(genome, decoys_output_file, fout=sys.stdout):
    '''
    Append genome sequence to output 3'UTRome to generate a 'decoy-aware' 3'UTRome for use with Salmon >= v1.0.0's selective alignment procedure
    Follows recipe outlined here: https://combine-lab.github.io/alevin-tutorial/2019/selective-alignment/
    '''

    ids = []
    handle = open(genome, 'r')

    for record in SeqIO.parse(genome, "fasta"):
        ids.append(record.id)
        SeqIO.write(record, fout, "fasta")

    handle.close()

    # Write seqnames to a text file (one per-line)
    with open(decoys_output_file, "w") as decoys_out:
        for id in ids:
            decoys_out.write(f"{id}\n")

    logger.info(f"Outputted decoy seqnames to {decoys_output_file}")


def main(args):
    if args.debug:
        logger.setLevel(logging.DEBUG)

    seqs = get_sequences(args.bed_file[0], args.genome)

    fout = open(args.output_file[0], "w")
    filter_sequences(seqs.seqfn, fout=fout)

    if args.decoys:
        logger.info(f"Append genome sequence - {args.genome} - as decoys to the output 3'UTR FASTA - {args.output_file[0]}")
        add_decoys(args.genome, args.decoys_output_file, fout=fout)

    fout.close()
