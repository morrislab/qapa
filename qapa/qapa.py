from __future__ import print_function
import sys
import os
import os.path
import argparse
import tempfile

from . import extract
from . import annotate
from . import extend
from . import collapse
from . import fasta

__version__ = '1.0.0'


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


def getoptions(args=None):
    desc = """
RNA-seq Quantification of Alternative Polyadenylation (QAPA)

Choose one of the sub-commands below. Include '-h' for more information.
Note: unless otherwise specified, all input files can be in compressed
(.gz) format."""
    # url = "URL: https://www.github.com/morrislab/QAPA"
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description=desc,
                                     fromfile_prefix_chars='@')
    parser.add_argument('-v', '--version', action='version',
                        version='%(prog)s ' + __version__)
    parser.add_argument('-t', '--temp', type=str, help='Set temp directory')
    optional = parser._action_groups.pop()
    subparsers = parser.add_subparsers(title="sub-commands")
    parser._action_groups.append(optional)

    # build utrs
    desc = """Extract 3' UTRs from GENCODE annotation table in genePred format,
           followed by annotation with GENCODE poly(A) track and PolyASite
           (http://polyasite.unibas.ch). Output is a BED file sent to STDOUT.
           """
    build_parser = subparsers.add_parser('build', description=desc,
                                         help="Build 3' UTR reference library")
    optional = build_parser._action_groups.pop()
    build_parser.add_argument('annotation_file', nargs=1,
                              help="Input annotation table")
    # build_parser.add_argument('output_file', nargs=1, help='Output filename')

    required = build_parser.add_argument_group('required named arguments')
    required.add_argument("--db", type=str, required=True,
                          help="Ensembl gene identifier table")
    required.add_argument('-g', '--gencode_polya', dest="gencode_polya",
                          required=True,
                          help="GENCODE poly(A) site track")
    required.add_argument('-p', '--polyasite', dest="polyasite", required=True,
                          help="Gruber et a. POLYASITE database")
    optional.add_argument('-m', '--min_polyasite', dest="min_polyasite",
                          type=int, default=3,
                          help="Minimum score in POLYASITE [%(default)s]")
    optional.add_argument('-i', '--intermediate_polyasite',
                          dest="intermediate_polyasite", type=int, default=4,
                          help="Minimum score in POLYASITE for creating "
                          "intermedate PAS entries [%(default)s]")

    optional.add_argument("-d", type=int, default=24,
                          dest="dist3", metavar="DISTANCE",
                          help="Maximum distance between 3' ends to merge [%(default)s]")
    optional.add_argument("-f", type=int, default=3,
                          dest="dist5", metavar="DISTANCE",
                          help="Maximum distance between 5' ends to merge [%(default)s]")
    optional.add_argument("-s", "--save", action='store_true',
                          help="Don't automatically delete intermediate files")
    build_parser.set_defaults(func=build)
    build_parser._action_groups.append(optional)

    # fetch FASTA
    desc = """Calls bedtools getfasta to extract sequences from BED file,
           followed by removal of duplicates and short sequences (which
           can cause issues with Sailfish indices)"""
    fasta_parser = subparsers.add_parser('fasta', description=desc,
                                         help="Extract FASTA "
                                         "sequences from from BED file")
    fasta_parser.add_argument('bed_file', nargs=1, help='Input BED filename')
    fasta_parser.add_argument('output_file', nargs=1, help='Output filename')
    optional = fasta_parser._action_groups.pop()
    required = fasta_parser.add_argument_group('required named arguments')
    required.add_argument('-f', '--fi', type=str, dest='genome', required=True,
                          help='Genome FASTA file (*uncompressed)')
    fasta_parser.set_defaults(func=fetch_sequences)
    fasta_parser._action_groups.append(optional)

    # quantify APA
    desc = """Collects 3' UTR expression quantifications from one or more
          samples and merges them into a single table, followed by Poly(A) Site
          Usage (PAU) calculation. This step requires 3' UTR expression
          quantification to be already carried out (e.g. via Sailfish)."""
    quant_parser = subparsers.add_parser('quant', description=desc,
                                         help='Compute PAU for one or more '
                                              'samples')
    quant_parser.add_argument('quant_files', nargs='+',
                              help="Filepaths of one or more 3' UTR "
                                   "quantification files. Expects each file "
                                   "to be inside its own directory (e.g. "
                                   "./path/to/samples_*/quant.sf). The "
                                   "directory base name will be used as the "
                                   "sample name.")
    optional = quant_parser._action_groups.pop()
    required = quant_parser.add_argument_group('required named arguments')
    required.add_argument("--db", type=str, required=True,
                          help="Ensembl gene identifier table")
    optional.add_argumenet('-f', '--field', type=str, metavar='FIELD',
                           default='TPM',
                           help='Field to merge [%(default)s]')
    optional.add_argument('-s', '--save', type=str, metavar='FILE',
                          help='Save intermediate file of merged samples as '
                               'FILE')
    quant_parser.set_defaults(func=quant)
    quant_parser._action_groups.append(optional)

    args = parser.parse_args(args=args)

    if args.temp:
        eprint("[Setting temporary directory to {}]".format(args.temp))
        tempfile.tempdir = args.temp

    return args


###############################################################################


def build(args):
    tf1 = tempfile.NamedTemporaryFile(mode='w', prefix='qapa_extract_',
                                      delete=False)
    tf2 = tempfile.NamedTemporaryFile(mode='w', prefix='qapa_anno_',
                                      delete=False)
    tf3 = tempfile.NamedTemporaryFile(mode='w', prefix='qapa_extend_',
                                      delete=False)
    try:
        # 1) get last exon from table
        eprint("[Extracting 3' UTRs from table]")
        extract.main(args, tf1)
        eprint(tf1.name)
        tf1.close()

        # 2) annotate 3' ends
        eprint("[Annotating 3' UTRs]")
        annotate.main(args, tf1.name, tf2)
        eprint(tf2.name)
        tf2.close()

        # 3) extend 5'
        eprint("[Checking 5' ends]")
        result = extend.main(args, tf2.name)
        result.to_csv(tf3, sep="\t", index=False, header=True)
        eprint(tf3.name)
        tf3.close()

        # # 4) collapse 3' ends
        eprint("[Collapsing 3' ends]")
        # fout = open(args.output_file[0], 'w')
        result = collapse.merge_bed(args, tf3.name)
        result.to_csv(sys.stdout, sep="\t", index=False, header=False)
        # fout.close()
    except Exception as e:
        eprint("Error: {}".format(e))
    finally:
        if not args.save:
            os.unlink(tf1.name)
            os.unlink(tf2.name)
            os.unlink(tf3.name)


def fetch_sequences(args):
    fasta.main(args)
    eprint("[Sequences written to {}]".format(args.output_file[0]))


def quant(args):
    intermediate_name = ''

    if not args.save:
        merged_tmp = tempfile.NamedTemporaryFile(prefix='qapa_',
                                                 suffix='_merge',
                                                 delete=False)
        intermediate_name = merged_tmp.name
    else:
        intermediate_name = args.save

    try:
        cmd = "create_merged_data.R --ensembl {} -f {} {} > {}".format(
            args.db, args.field, " ".join(args.quant_files), intermediate_name)
        # eprint(cmd)
        os.system(cmd)

        cmd = "compute_pau.R -e {}".format(intermediate_name)
        # eprint(cmd)
        os.system(cmd)

    finally:
        if not args.save:
            os.unlink(intermediate_name)


def main():
    args = getoptions()
    args.func(args)


if __name__ == '__main__':
    main()
