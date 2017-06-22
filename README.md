## RNA-seq Quantification of Alternative Polyadenylation (QAPA)

Tools for analyzing alternative polyadenylation (APA) from RNA-seq
data (human and mouse). QAPA consists of two main components:

  1. Extraction and annotation of 3' UTR sequences from gene models
  1. Calculation of relative usage of alternative 3' UTR isoforms based on
     transcript-level abundance.

Note that QAPA itself does not perform transcript quantification. It relies on
other tools such as [Sailfish](https://github.com/kingsfordgroup/sailfish) and
[Salmon](https://github.com/COMBINE-lab/salmon).

## Installation

QAPA consists of both Python (2.7+ or 3.5+) and R scripts.

1. Install [R](https://www.r-project.org/) and the R packages optparse, dplyr,
   data.table, and stringr. In the R console, the packages can be installed in
   one line:

        install.packages(c("optparse", "dplyr", "data.table", "stringr"))

2. Install the latest development version of
   [bedtools](https://github.com/arq5x/bedtools2). *Note: do not use 2.26.0
   stable release as there is a
   [bug](https://github.com/arq5x/bedtools2/issues/435) with the groupBy tool*.

3. Install the latest QAPA source code from GitHub. QAPA requires the Python
   packages pandas, numpy, pybedtools, and biopython. These will be
   automatically installed, if necessary.

        git clone git@github.com:morrislab/qapa.git
        cd qapa
        python setup.py install

        # To test if installation is working:
        cd
        which qapa
        qapa -h

## Usage

QAPA has three sub-commands: `build`, `fasta`, and `quant`.

### 3' UTR library extraction (`build` and `fasta`)

#### Pre-requisites

The following data sources are required:

1. Ensembl gene metadata table from [Biomart](http://www.ensembl.org/biomart/).
   Human and mouse tables are provided in the `examples` folder.  To obtain a fresh
   copy, download a table of Ensembl Genes from Biomart with the following
   attributes: Ensembl Gene ID, Ensembl Transcript ID, Gene Type, Transcript Type,
   and Gene Name. Alternatively, download via MySQL (see
   [here](http://www.ensembl.org/info/data/mysql.html) for more details):

        mysql --user=anonymous --host=martdb.ensembl.org --port=5316 -A ensembl_mart_89
            -e "select stable_id_1023 as 'Gene stable ID', stable_id_1066 as 'Transcript
            stable ID', biotype_1020 as 'Gene type', biotype_1064 as 'Transcript type',
            display_label_1074 as 'Gene name' from mmusculus_gene_ensembl__transcript__main"
            > ensembl_identifiers.txt

   To change the species, replace the table name (e.g. for human, use
   `hsapiens_gene_ensembl__transcript__main`).

2. PolyASite database

   Download BED files (human or mouse) from http://polyasite.unibas.ch/.

3. GENCODE poly(A) sites track

   Download from UCSC Table Browser or alternatively via MySQL (see
   [here](https://genome.ucsc.edu/goldenpath/help/mysql.html) for more details):

        mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -e "select chrom, txStart,
        txEnd, name2, score, strand from wgEncodeGencodePolyaVM9 where name2 =
        'polyA_site'" -N mm10 > gencode.polyA_sites.bed

4. GENCODE gene prediction annotation table

   Download from UCSC Table Browser or alternatively via MySQL:

        mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -e "select * from
            wgEncodeGencodeBasicVM9" mm10 > gencode.basic.txt

   Note that the `-N` option (suppress column headings) is not used here.

5. A reference genome in FASTA format is required for extracting sequences from BED files. Can be downloaded from http://hgdownload.soe.ucsc.edu/downloads.html.

#### Run

To extract 3' UTRs from annotation, run:

    qapa build --db ensembl_identifiers.txt -g gencode.polyA_sites.bed
        -p polyASite.bed gencode.basic.txt > output_utrs.bed

Results will be saved in the file `output_utrs.bed` (default is STDOUT).

To extract sequences from the resulting BED file, use the `fasta` sub-command:

    qapa fasta -f genome.fa output_utrs.bed output_sequences.fa

Sequences will be saved in `output_sequences.fa`. This file can then be indexed
by transcript quantification tools like Sailfish/Salmon.

### Quantification of 3' UTR isoform usage (`quant`)

Expression quantification of 3' UTR isoforms must be carried out first. Once
this is done, QAPA expects the quantification results to be located inside its
own directory. For example, Sailfish/Salmon results may appear like the following:

    project/
      |-- sample1/quant.sf
      |-- sample2/quant.sf

The `quant` sub-command calls two R scripts that merge the results from each
sample into a single table, followed by computing the relative proportion of
each isoform in a gene, measured as Poly(A) Usage (PAU).

    qapa quant --db ensembl_identifiers.mm10.txt project/sample*/quant.sf > apa_results.txt

Results will be saved in the file `apa_results.txt` (default is STDOUT).

The output format is following:

1. APA_ID: unique identifier consisting of in the format `<Ensembl Gene ID>_<number>_<(P|D|S)>`, where P = proximal, D = distal, and S = single
1. Transcript: one or more Ensembl Transcript IDs
1. Gene: Ensembl Gene ID
1. Gene_Name: gene symbol
1. Chr: chromosome
1. LastExon.Start: start coordinate of last exon
1. LastExon.End: end coordinate of last exon
1. Strand: `+` or `-`
1. UTR3.Start: start coordinate of 3' UTR
1. UTR3.End: end coordinate of 3' UTR
1. Length: length of the 3' UTR
1. Num_Events: number of PAS per gene
1. *sample1*.PAU: PAU estimate for *sample1*
1. *sample2*.PAU: PAU estimate for *sample2*
1. *sample1*.TPM: TPM estimate for *sample1*
1. *sample2*.TPM: TPM estimate for *sample2*

