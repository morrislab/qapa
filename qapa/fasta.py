from __future__ import print_function
import fileinput
import sys
from Bio import SeqIO
import pybedtools


def get_sequences(bed_file, genome):
    bed = pybedtools.BedTool(bed_file)
    return bed.sequence(genome, s=True, name=True, fullHeader=False,
                        split=True)


def filter_sequences(fasta_file, min_length=100, fout=sys.stdout):
    seqs = set()
    ids = set()
    skipped = 0
    handle = fileinput.input(fasta_file)

    for record in SeqIO.parse(handle, "fasta"):
        sequence = str(record.seq)

        if len(sequence) > min_length and \
                not (sequence in seqs or record.id in ids):
            seqs.add(sequence)
            seqs.add(record.id)
            SeqIO.write(record, fout, "fasta")
        else:
            skipped += 1

    handle.close()

    if skipped > 0:
        print("[%s] Skipped %d sequences" % ('filter_fasta', skipped),
              file=sys.stderr)


def main(args):
    seqs = get_sequences(args.bed_file[0], args.genome)

    fout = open(args.output_file[0], "w")
    filter_sequences(seqs.seqfn, fout=fout)
    fout.close()
