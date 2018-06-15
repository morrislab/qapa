# RNA-seq Quantification of Alternative Polyadenylation (QAPA)

Analysis of alternative polyadenylation (APA) from RNA-seq
data (human and mouse). QAPA consists of two main components:

  1. Extraction and annotation of 3' UTR sequences from gene models
  1. Calculation of relative usage of alternative 3' UTR isoforms based on
     transcript-level abundance.

Note that QAPA itself does not perform transcript quantification. It relies on
other tools such as [Sailfish](https://github.com/kingsfordgroup/sailfish) and
[Salmon](https://github.com/COMBINE-lab/salmon).

# Installation

QAPA consists of both Python (2.7+ or 3.5+) and R scripts.

1. Install the latest development version of
   [bedtools](https://github.com/arq5x/bedtools2). *Note: do not use 2.26.0
   stable release as there is a
   [bug](https://github.com/arq5x/bedtools2/issues/435) with the groupBy tool*.

2. Install the latest development version (or the latest
   [release](https://github.com/morrislab/qapa/releases/latest)) of QAPA from
   GitHub . QAPA requires the Python packages pandas, numpy, pybedtools, and
   biopython. These will be automatically installed, if necessary.

        git clone https://github.com/morrislab/qapa.git
        cd qapa
        python setup.py install

3. Install [R](https://www.r-project.org/) and the R packages optparse, dplyr,
   data.table, and stringr. In the R console, the packages can be installed in
   one line:

        install.packages(c("optparse", "dplyr", "data.table", "stringr"))

   Alternatively, run the provided `install_packages.R` helper script from
   command line:

        Rscript scripts/install_packages.R

4. To test if installation is working:

        cd          # change to root directory
        which qapa  # should return path of qapa executable
        qapa -h     # should display help message

# Usage

QAPA has three sub-commands: `build`, `fasta`, and `quant`. Unless specified,
all input files can be in compressed (gzip) format.

## 3' UTR library extraction (`build` and `fasta`)

Pre-defined libraries for human and mouse are available for download below. Otherwise, continue to the next section to build from scratch.

  - [Human (hg19)](https://zenodo.org/record/1222196/files/qapa_3utrs.gencode.hg19.tar.gz)
  - [Mouse (mm10)](https://zenodo.org/record/1222196/files/qapa_3utrs.gencode.mm10.tar.gz)

The following data sources are required:

**A. Gene annotation**

1. Ensembl gene metadata table from [Biomart](http://www.ensembl.org/biomart).
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

        mysql --user=anonymous --host=martdb.ensembl.org --port=5316 -A ensembl_mart_89
            -e "select stable_id_1023 as 'Gene stable ID', stable_id_1066 as 'Transcript
            stable ID', biotype_1020 as 'Gene type', biotype_1064 as 'Transcript type',
            display_label_1074 as 'Gene name' from mmusculus_gene_ensembl__transcript__main"
            > ensembl_identifiers.txt

   To change the species, replace the table name (e.g. for human, use
   `hsapiens_gene_ensembl__transcript__main`).

2. GENCODE gene prediction annotation table

   Download from [UCSC Table Browser](https://genome.ucsc.edu/cgi-bin/hgTables)
   or alternatively via MySQL. For example, to download mm10 gene predictions:

        mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -e "select * from
            wgEncodeGencodeBasicVM9" mm10 > gencode.basic.txt

   Alternatively, if you are starting from a GTF/GFF file, you can convert
   it to genePred format using the UCSC tool
   [`gtfToGenePred`](http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/gtfToGenePred):

        gtfToGenePred -genePredExt custom_genes.gtf custom_genes.genedPred

**B. Poly(A) site annotation**

As of v1.2.0, this step is optional. Otherwise, two options are available:

Option 1: standard approach (as described in the [paper](#citation))

1. PolyAsite database

   Download BED files (human or mouse) from http://polyasite.unibas.ch/.

2. GENCODE poly(A) sites track

   Download from [UCSC Table Browser](https://genome.ucsc.edu/cgi-bin/hgTables)
   or alternatively via MySQL (see
   [here](https://genome.ucsc.edu/goldenpath/help/mysql.html) for more details).
   For example, to download mm10 annotations:

        mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -e "select chrom, txStart,
        txEnd, name2, score, strand from wgEncodeGencodePolyaVM9 where name2 =
        'polyA_site'" -N mm10 > gencode.polyA_sites.bed

Option 2: use custom BED track

1. Custom BED track of poly(A) sites

    A custom BED file of poly(A) can be used to annotate 3' UTRs.
    Each entry must contain the start (0-based) and end coordinate of a poly(A)
    site.

**C. Reference genome**

A reference genome in FASTA format is required for extracting sequences from
   BED files. e.g. 
   http://hgdownload.soe.ucsc.edu/downloads.html.

### Run

To extract 3' UTRs from annotation, run:

    qapa build --db ensembl_identifiers.txt -g gencode.polyA_sites.bed -p clusters.mm10.bed 
        gencode.basic.txt > output_utrs.bed

If using a custom BED file, replace the `-g` and `-p` options with `-o`:

    qapa build --db ensembl_identifiers.txt -o custom_sites.bed
        gencode.basic.txt > output_utrs.bed

If using a custom genePred file converted from GTF, include the `-H`
option:

    qapa build -H --db ensembl_identifiers.txt -o custom_sites.bed
        custom_genes.genePred > output_utrs.bed
 
If bypassing the poly(A) annotation step, include the `-N` option:

    qapa build -N --db ensembl_identifiers.txt gencode.basic.txt > output.utrs.bed

Results will be saved in the file `output_utrs.bed` (default is STDOUT).

To extract sequences from the resulting BED file, use the `fasta` sub-command
(calls `bedtools getfasta`) along with the downloaded reference genome file:

    qapa fasta -f genome.fa output_utrs.bed output_sequences.fa

Note that `genome.fa` must be uncompressed. Sequences will be saved in
`output_sequences.fa`. This file can then be indexed by transcript
quantification tools like Sailfish/Salmon:
    
    salmon index -t output_sequences.fa -i utr_library
     
It is important that the sequence IDs
are not modified as it will be parsed by the `quant` sub-command.

## Quantification of 3' UTR isoform usage (`quant`)

Expression quantification of 3' UTR isoforms must be carried out first. Once
this is done, QAPA expects the quantification results to be located inside its
own sub-directory. For example, typical Sailfish/Salmon results may appear with
the following directory structure:

    project/
      |-- sample1/quant.sf
      |-- sample2/quant.sf
      |-- (etc.)

The `quant` sub-command is a wrapper that calls two R scripts:
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
UTR3.Start | start coordinate of 3' UTR
UTR3.End | end coordinate of 3' UTR
Length | length of the 3' UTR
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
