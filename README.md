# RNA-seq Quantification of Alternative Polyadenylation (QAPA)

Analysis of alternative polyadenylation (APA) from RNA-seq
data (human and mouse). QAPA consists of two main components:

  1. Extraction and annotation of 3′ UTR sequences from gene models
  1. Calculation of relative usage of alternative 3′ UTR isoforms based on
     transcript-level abundance.

Note that QAPA itself does not perform transcript quantification. It relies on
other tools such as [Sailfish](https://github.com/kingsfordgroup/sailfish) and
[Salmon](https://github.com/COMBINE-lab/salmon).

---

## Installation

QAPA consists of both Python and R scripts.  A conda virtual environment
can be created using the provided `environment.yml` file.

1. Clone the repo:

        git clone https://github.com/morrislab/qapa.git
        cd qapa

2. (Optional) Install `mamba` for faster Conda management

        conda install -c conda-forge mamba

3.  Create the environment

        mamba env create -f environment.yml
        conda activate qapa

4. Test that `qapa` command is available

        qapa --help
---

# Usage

QAPA has three sub-commands:

  1. [`build`](#build-3-utrs-from-annotation-build): Generate a 3′ UTR
     library from annotations
  2. [`fasta`](#extract-3-utr-sequences-fasta): Extract sequences for
     indexing by transcript quantification tools
  3. [`quant`](#quantify-3-utr-isoform-usage-quant): Calculate relative
     usage of alternative 3′ UTR isoforms

## 1) Build 3′ UTRs from annotation (`build`)

### Download pre-compiled libraries

Pre-compiled libraries for human and mouse are available for download below.
Otherwise, continue reading to build from scratch.

_Updated in [v1.3.0](https://github.com/morrislab/qapa/releases/tag/v1.3.0),
the following libraries are pre-compiled with Polyasite V2_:

  - [Human (hg38)](https://github.com/morrislab/qapa/releases/download/v1.3.0/qapa_3utrs.gencode_V31.hg38.bed.gz)
  - [Mouse (mm10)](https://github.com/morrislab/qapa/releases/download/v1.3.0/qapa_3utrs.gencode_VM22.mm10.bed.gz)

(Old) Versions v1.2.3 and earlier:

  - [Human (hg19)](https://zenodo.org/record/1222196/files/qapa_3utrs.gencode.hg19.tar.gz)
  - [Mouse (mm10)](https://zenodo.org/record/1222196/files/qapa_3utrs.gencode.mm10.tar.gz)

### Prepare annotation files

To run `build`, gene and poly(A) annotation sources need to be prepared:

**A. Gene annotation**

1. Ensembl gene metadata table from
   [Biomart](http://www.ensembl.org/biomart/martview).

   Human and mouse tables are provided in the `examples` folder.  To obtain a fresh
   copy, download a table of Ensembl Genes from Biomart with the following
   attributes:

   1. Ensembl Gene ID
   1. Ensembl Transcript ID
   1. Gene Type
   1. Transcript Type
   1. Gene Name

   Alternatively, download via MySQL (see
   [here](http://www.ensembl.org/info/data/mysql.html) for more details):

        mysql --user=anonymous --host=martdb.ensembl.org --port=5316 -A ensembl_mart_89 \
            -e "select stable_id_1023 as 'Gene stable ID', stable_id_1066 as 'Transcript stable ID', \
            biotype_1020 as 'Gene type', biotype_1064 as 'Transcript type', \
            display_label_1074 as 'Gene name' from mmusculus_gene_ensembl__transcript__main" \
            > ensembl_identifiers.txt

   To change the species, replace the table name (e.g. for human, use
   `hsapiens_gene_ensembl__transcript__main`).

2. GENCODE gene prediction annotation table

   Download from [UCSC Table Browser](https://genome.ucsc.edu/cgi-bin/hgTables)
   or alternatively via MySQL (see
   [here](https://genome.ucsc.edu/goldenpath/help/mysql.html) for more details).
   For example, to download mm10 gene predictions:

        mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A \
            -e "select * from wgEncodeGencodeBasicVM9" mm10 > gencode.basic.txt

   Alternatively, if you are starting from a GTF/GFF file, you can convert
   it to genePred format using the UCSC tool
   [`gtfToGenePred`](http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/gtfToGenePred):

        gtfToGenePred -genePredExt custom_genes.gtf custom_genes.genePred

   Note that it is important to include the `-genePredExt` option!

**B. Poly(A) site annotation (optional)**

This step can be skipped, otherwise continue reading. Two options are available:

**Option 1**: standard approach using PolyASite and GENCODE poly(A) track (as described in the [paper](#citation))

1. PolyASite database

   Download BED files (human or mouse) from http://polyasite.unibas.ch/.

2. GENCODE poly(A) sites track

   Download from [UCSC Table Browser](https://genome.ucsc.edu/cgi-bin/hgTables)
   or alternatively via MySQL. For example, to download mm10 annotations:

        mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A \
            -e "select chrom, txStart, txEnd, name2, score, strand \
            from wgEncodeGencodePolyaVM9 where name2 = 'polyA_site'" -N mm10 \
            > gencode.polyA_sites.bed

**Option 2**: use custom BED track of poly(A) sites

A custom BED file of poly(A) sites can be used to annotate 3′ UTRs.
Each entry must contain the start (0-based) and end coordinate of a poly(A)
site.

### Commands

Once the data files have been prepared, we can then use `build` to create the 3'
UTR library. See `qapa build -h` for usage details. The following
describes several example use cases:

1. To extract 3′ UTRs from annotation, run:

       qapa build --db ensembl_identifiers.txt -g gencode.polyA_sites.bed -p clusters.mm10.bed gencode.basic.txt > output_utrs.bed

2. If using a custom BED file, replace the `-g` and `-p` options with `-o`:

       qapa build --db ensembl_identifiers.txt -o custom_sites.bed gencode.basic.txt > output_utrs.bed

3. If using a custom genePred file converted from GTF, supply the file as in 1.
   (e.g. the first positional argument):

       qapa build --db ensembl_identifiers.txt -o custom_sites.bed custom_genes.genePred > output_utrs.bed

4. If bypassing the poly(A) annotation step, include the `-N` option:

       qapa build -N --db ensembl_identifiers.txt gencode.basic.txt > output.utrs.bed

Results will be saved in the file `output_utrs.bed` (default is STDOUT).
It is important that the sequence IDs are not modified as it will be parsed by
`quant` below.

Additional notes:

  - 3' UTRs that contain introns will be skipped.
  - Chromosome names that contain underscores are currently not supported and will be
    skipped.
  - Only genes with `Gene Type = 'protein_coding'` will be considered.

### Troubleshooting tips

  - Use `--debug` option to produce more verbose logging messages
  - Use `--save` and `--temp <dir_path>` to save intermediate files generated
  by `qapa build`. `<dir_path>` should be a user accessible directory.


## 2) Extract 3′ UTR sequences (`fasta`)

To extract sequences from the BED file prepared by `build`, a reference genome in
FASTA format is required. e.g. http://hgdownload.soe.ucsc.edu/downloads.html.

Then, run the command:

    qapa fasta -f genome.fa output_utrs.bed output_sequences.fa

Essentially `fasta` is a wrapper that calls `bedtools getfasta`. Note that
`genome.fa` must be uncompressed. Sequences will be saved in
`output_sequences.fa`.

Since the initial publication, Salmon added a 'selective alignment' procedure ([Srivastava, Malik et al., Genome Biology, 2020](https://doi.org/10.1186/s13059-020-02151-8)). The 'full' selective alignment jointly aligns the reads to the input transcriptome and the genome, discarding alignments that map better to the genome than the transcriptome. To generate a fasta file compatible with Salmon's selective alignment procedure (a 'decoy-aware transcriptome'), run `qapa fasta` with the `--decoys` flag:

    qapa fasta -f genome.fa --decoys output_utrs.bed output_sequences.fa

This will save the transcript sequences appended with the genome sequence from `genome.fa` to `output_sequences.fa`. This also generates a file named `decoys.txt`, which contains the sequence identifiers from the genome sequence to differentiate target transcripts from decoy sequences. The name of the text file can be controlled using the optional `-d` or `--decoys-output-file` argument.

## 3) Quantify 3′ UTR isoform usage (`quant`)

Expression quantification of 3′ UTR isoforms must be carried out first using the
FASTA file prepared by `fasta` as the index. For example, to index the sequences
using Salmon:

    salmon index -t output_sequences.fa -i utr_library

If you generated a genome-augmented fasta file with `qapa fasta --decoys`, additionally pass the `decoys.txt` file to `salmon index`:

    salmon index -t output_sequences.fa -i utr_library --decoys decoys.txt

Following expression quantification, QAPA expects the results to be located inside its
own sub-directory. For example, typical Sailfish/Salmon results may appear with
the following directory structure:

    project/
      |-- sample1/quant.sf
      |-- sample2/quant.sf
      |-- (etc.)

The `quant` sub-command calls two R scripts:
`create_merged_table.R` and `compute_pau.R`. The first script combines the
quantifications from each sample into a single table. The second script computes
the relative proportion of each isoform in a gene, measured as Poly(A) Usage
(PAU) (`compute_pau.R`).

    qapa quant --db ensembl_identifiers.mm10.txt project/sample*/quant.sf > pau_results.txt

Results will be saved in the file `pau_results.txt` (default is STDOUT).

For advanced usage, the R scripts can be run on its own. Run
`create_merged_table.R -h` and `compute_pau.R -h` for usage details.

The final output format is as follows:

Column | Description
------ | -----------
APA_ID | unique identifier consisting of in the format `<Ensembl Gene ID>_<number>_<(P\|D\|S)>`, where P = proximal, D = distal, and S = single
Transcript | one or more Ensembl Transcript IDs
Gene | Ensembl Gene ID
Gene_Name | gene symbol
Chr | chromosome
LastExon.Start | start coordinate of last exon
LastExon.End | end coordinate of last exon
Strand | + or -
UTR3.Start | start coordinate of 3′ UTR
UTR3.End | end coordinate of 3′ UTR
Length | length of the 3′ UTR
Num_Events | number of PAS per gene
*sample1*.PAU | PAU estimate for *sample1*
*sample2*.PAU | PAU estimate for *sample2*
*sample1*.TPM | TPM estimate for *sample1*
*sample2*.TPM | TPM estimate for *sample2*

# Citation

Ha, K.C.H., Blencowe, B.J., Morris, Q. (2018). QAPA: a new method for the
systematic analysis of alternative polyadenylation from RNA-seq data. Genome
Biol. 19, 45.
[[link]](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-018-1414-4)
