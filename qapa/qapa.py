import sys
import os
import os.path
import re
import argparse
import tempfile
import logging

from .version import __version__

logging.basicConfig(level=logging.INFO, stream=sys.stderr,
                    format='%(name)-13s - %(asctime)s - %(levelname)-8s - %(message)s')
logger = logging.getLogger(__name__)

def _check_input_files(inputs, parser):
    for input_file in inputs:
        if input_file and not os.path.exists(input_file):
            parser.error("No such file: {}".format(input_file))


def getoptions(args=None):
    desc = """
RNA-seq Quantification of Alternative Polyadenylation (QAPA)

Choose one of the sub-commands below. Include '-h' for more information.
Note: unless otherwise specified, all input files can be in compressed
(.gz) format."""
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description=desc,
                                     fromfile_prefix_chars='@')
    parser.add_argument('-v', '--version', action='version',
                        version='%(prog)s ' + __version__)
    optional = parser._action_groups.pop()
    subparsers = parser.add_subparsers(title="sub-commands", dest='subcommand')
    parser._action_groups.append(optional)

    # common args
    common = argparse.ArgumentParser(add_help=False)
    common.add_argument('-t', '--temp', type=str,
                        help="Set temp directory [{}]".
                        format(tempfile.gettempdir()))

    # common args for build and fasta
    common_build = argparse.ArgumentParser(add_help=False)
    common_build.add_argument('--debug', action="store_true",
                        help="Increase verbosity of log messages for"
                        " troubleshooting.")

    # build utrs
    desc = """
Extract 3' UTRs from GENCODE annotation table in genePred format,
followed by annotation with GENCODE poly(A) track and PolyAsite.
Alternatively, the second step can be carried out using a custom BED
file using option -o.

Output is in BED format plus additional gene symbol column
(to STDOUT).
           """
    build_parser = subparsers.add_parser('build', description=desc,
                                         formatter_class=argparse.RawDescriptionHelpFormatter,
                                         help="build 3' UTR reference library",
                                         parents=[common, common_build])
    optional = build_parser._action_groups.pop()
    build_parser.add_argument('annotation_file', nargs=1,
                              help="Input annotation table in genePred format."
                              " Recognizes UCSC tables or genePred files "
                              " generated from GTF format using 'gtfToGenePred"
                              " -genePredExt' command.")

    required = build_parser.add_argument_group('required named arguments')
    required.add_argument("--db", type=str, required=True,
                          help="Ensembl gene identifier table")
    required_polya = build_parser.add_argument_group(
        'required named arguments when using public poly(A) site annotations'
    )
    required_polya.add_argument('-g', '--gencode_polya', dest="gencode_polya",
                                help="GENCODE poly(A) site track")
    required_polya.add_argument('-p', '--polyasite', dest="polyasite",
                                help="PolyAsite database")
    optional.add_argument('-m', '--min_polyasite', dest="min_polyasite",
                          type=int, default=3,
                          help="Minimum score in PolyAsite [%(default)s]")
    optional.add_argument('-i', '--intermediate_polyasite',
                          dest="intermediate_polyasite", type=int, default=4,
                          help="Minimum score in PolyAsite for creating "
                          "intermedate PAS entries [%(default)s]")
    optional.add_argument("-e", type=int, default=0,
                          dest="numextends", metavar="DISTANCE",
                          help="Number of exons to extend in 5' direction. "
                          "Used for resolving different 5' ends of "
                          "overlapping 3' UTRs. Setting 0 will exclude "
                          "internal 3' UTRs. [%(default)s]")
    optional.add_argument("-d", type=int, default=24,
                          dest="dist3", metavar="DISTANCE",
                          help="Maximum distance between 3' ends to merge [%(default)s]")
    optional.add_argument("-o", "--other", default=None,
                          help="Use this option to annotate 3' UTRs with a "
                          "custom BED file of poly(A) sites. This option "
                          "cannot be used in conjunction with -g and -p.")
    optional.add_argument("-C", "--no_skip_random_chromosomes", default=True,
                          action="store_true",
                          help="Disable skipping of chromosomes that don't "
                          "match the regular expression pattern "
                          " '^(chr)*[0-9XYM]+$'")
    optional.add_argument("-s", "--save", action='store_true',
                          help="Don't automatically delete intermediate files")
    optional.add_argument("-N", "--no_annotation", action='store_true',
                          help="Skip annotation step. Use this option if you "
                          "only have a gene model annotation file to build "
                          "a 3' UTR library.")
    optional.add_argument("--species", type=str,
                          help="Set species descriptor. Useful if not using 'hsa' or "
                          "'mmu', otherwise 'unk' is used.")
    optional.add_argument("-c", "--cores", type=int,
                          help="The number of cores for multiprocessing. "
                          "Defaults to all available cores.")
    build_parser.set_defaults(func=build)
    build_parser._action_groups.append(optional)

    # fetch FASTA
    desc = """Calls bedtools getfasta to extract sequences from BED file,
           followed by removal of duplicates and short sequences (which
           can cause issues with Sailfish indices)"""
    fasta_parser = subparsers.add_parser('fasta', description=desc,
                                         help="extract FASTA "
                                         "sequences from BED file",
                                         parents=[common, common_build])
    fasta_parser.add_argument('bed_file', nargs=1, help='Input BED filename')
    fasta_parser.add_argument('output_file', nargs=1, help='Output filename')
    optional = fasta_parser._action_groups.pop()
    required = fasta_parser.add_argument_group('Required named arguments')
    required.add_argument('-f', '--fi', type=str, dest='genome', required=True,
                          help='Genome FASTA file (*uncompressed)')
    fasta_parser.set_defaults(func=fetch_sequences)
    fasta_parser._action_groups.append(optional)

    # quantify APA
    desc = """Merge 3' UTR expression quantifications from one or more
          samples into a single table, followed by Poly(A) Site
          Usage (PAU) calculation. This step requires 3' UTR expression
          quantification to be already carried out."""
    quant_parser = subparsers.add_parser('quant', description=desc,
                                         help='compute PAU for one or more '
                                              'samples',
                                         parents=[common])
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
    optional.add_argument('-f', '--field', type=str, metavar='FIELD',
                          default='TPM',
                          help='Field to merge [%(default)s]')
    optional.add_argument('-F', '--format', type=str, metavar='FORMAT',
                          choices=['sailfish', 'salmon'], default='salmon',
                          help="Specify transcript quantification method. For "
                               "Sailfish v0.8 or earlier, use 'sailfish'. "
                               "Otherwise, use 'salmon'. [%(default)s]")
    optional.add_argument('-s', '--save', type=str, metavar='FILE',
                          help='Save intermediate file of merged samples as '
                               'FILE')
    quant_parser.set_defaults(func=quant)
    quant_parser._action_groups.append(optional)

    args = parser.parse_args(args=args)

    if args.subcommand == 'build':
        if args.other is None and \
            (args.gencode_polya is None or args.polyasite is None) and \
            not args.no_annotation:
                parser.error("Missing arguments: -g and -p flags are mutually "
                             "inclusive")

        if args.other and (args.gencode_polya or args.polyasite):
            logger.info("Option -o (custom BED) will be used for build phase and "
                   "-g and -p will be ignored")

        if not (args.species is None or \
                re.match(r'^[a-zA-Z0-9]+$', args.species)):
            parser.error("Species must be alphanumeric.")

        _check_input_files([args.polyasite, args.gencode_polya, args.db,
                            args.other, args.annotation_file[0]], build_parser)

    elif args.subcommand == 'fasta':
        _check_input_files([args.bed_file[0], args.genome], fasta_parser)
    elif args.subcommand == 'quant':
        _check_input_files([args.db], quant_parser)
    else:
        parser.error("Specify sub-command: 'build', 'fasta', 'quant'")

    if args.temp:
        if not os.path.exists(args.temp):
            parser.error("No such directory: {}".format(args.temp))
        logger.info("Setting temporary directory to {}".format(args.temp))
        tempfile.tempdir = args.temp

    return args


###############################################################################


def build(args):
    from . import extract, extend, annotate, collapse
    if args.debug:
        logger.setLevel(logging.DEBUG)

    tf1 = tempfile.NamedTemporaryFile(mode='w', prefix='qapa_extract_',
                                      delete=False)
    tf2 = tempfile.NamedTemporaryFile(mode='w', prefix='qapa_anno_',
                                      delete=False)
    tf3 = tempfile.NamedTemporaryFile(mode='w', prefix='qapa_extend_',
                                      delete=False)

    try:
        # 1) get last exon from table
        logger.info("Extracting 3' UTRs from %s" % args.annotation_file[0])
        extract.main(args, tf1)
        logger.debug(tf1.name)
        tf1.close()

        # 2) annotate 3' ends
        logger.info("Annotating 3' UTRs")
        annotate.main(args, tf1.name, tf2)
        logger.debug(tf2.name)
        tf2.close()

        # 3) extend 5'
        logger.info("Checking 5' ends")
        result = extend.main(args, tf2.name)
        result.to_csv(tf3, sep="\t", index=False, header=True)
        logger.debug(tf3.name)
        tf3.close()

        # # 4) collapse 3' ends
        logger.info("Collapsing 3' ends")
        result = collapse.merge_bed(args, tf3.name)
        result.to_csv(sys.stdout, sep="\t", index=False, header=False)
    except Exception as e:
        logger.exception("Error occurred in build:")
        exit(1)
    finally:
        if not args.save:
            os.unlink(tf1.name)
            os.unlink(tf2.name)
            os.unlink(tf3.name)


def fetch_sequences(args):
    from . import fasta
    if args.debug:
        logger.setLevel(logging.DEBUG)
    fasta.main(args)
    logger.info("Sequences written to {}".format(args.output_file[0]))


def quant(args):
    intermediate_name = ''

    if not args.save:
        merged_tmp = tempfile.NamedTemporaryFile(prefix='qapa_merge_',
                                                 delete=False)
        intermediate_name = merged_tmp.name
    else:
        intermediate_name = args.save

    try:
        cmd = "create_merged_data.R --db {} -f {} -F {} {} > {}".format(
            args.db, args.field, args.format, " ".join(args.quant_files),
            intermediate_name)
        logger.info(cmd)
        os.system(cmd)

        cmd = "compute_pau.R -e {}".format(intermediate_name)
        logger.info(cmd)
        os.system(cmd)

    finally:
        if not args.save:
            os.unlink(intermediate_name)


def main():
    args = getoptions()
    logger.info("Version %s" % __version__)

    args.func(args)
    logger.info("Finished!")


if __name__ == '__main__':
    main()
