#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(optparse))

args <- commandArgs(TRUE)

option.list <- list(
  make_option(c("--db"), type="character", default=NULL,
              help="Ensembl identifiers database file. Required unless -m is
              specified. [%default]"),
  make_option(c("-f", "--field"), type="character",
              default="TPM", help="Field to merge [%default]"),
  make_option(c("-m", "--merge_only"), action="store_true", default=FALSE,
              help="Perform merge of fields only. Don't get additional metadata.
              [%default]"),
  make_option(c("-F", "--format"), type="character", default="salmon",
              help="Specify transcript quantification method. For Sailfish
              v0.8 or earlier, use 'sailfish'. Otherwise, use 'salmon'.
              [%default]"),
  make_option(c("-a", "--all_genes"), action="store_true", default=FALSE,
              help="If set to TRUE, do NOT filter for only protein-coding
              genes [%default]"),
  make_option(c("-n", "--non_standard"), action="store_true", default=FALSE,
              help="If set to TRUE, disables QAPA parsing of name column. A
              valid Ensembl Transcript ID is still required to extract
              additional metadata, unless -m is specified. [%default]")
  )
desc <- paste("\nMerge multlple quantification runs into a single summary table.",
              "Output is sent to STDOUT.")
parser <- OptionParser(option_list=option.list,
                       description = desc,
                       usage="usage: %prog [options] <quant files to merge>")
opt <- parse_args(parser, args=args, positional_arguments=TRUE)

if (length(opt$args) < 1) {
  write("Missing arguments", stderr())
  stop(print_help(parser))
}

if (opt$options$merge_only) {
  write("Merge-only mode enabled.", stderr())
} else if (is.null(opt$options$db)) {
  stop("Ensembl identifiers database is required.")
}

if (!opt$options$format %in% c("sailfish", "salmon")) {
  stop("QAPA currently supports sailfish or salmon formats")
}

suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dplyr))
#library(assertthat)

#### Functions #################################################################
join_iterative <- function(files, by = NULL, ...) {
  # Modified from plyr::join_all
  #
  # Join file by file to avoid loading everythig into memory
  pb <- txtProgressBar(style = 3, file = stderr())
  N = length(files)
  if (N == 1) {
    return(load_data(files[1], ...))
  }
  setTxtProgressBar(pb, 1/N)
  x <- load_data(files[1], ...)
  setkeyv(x, by)
  for (i in 2:length(files)) {
    setTxtProgressBar(pb, i/N)
    y <- load_data(files[i], ...)
    setkeyv(y, by)
    #assert_that(all(x$Transcript == y$Transcript) == TRUE)
    x <- y[x]
    #assert_that(all(x$Transcript == y$Transcript) == TRUE)
  }
  setTxtProgressBar(pb, 1)
  close(pb)
  #write("\n", stderr())
  return(x)
}

load_data <- function(path, format, field) {
  # Load data and keep only the selected field

  if (file.exists(path)) {

    if (format == "sailfish") {
      column.names <- c("Transcript", "Length", "TPM", "NumReads")

      if (! field %in% column.names) {
        stop(sprintf("The specified field (%s) cannot be found!", field))
      }

      m <- read.table(path, sep="\t", check.names = FALSE, quote = NULL,
                      stringsAsFactors=FALSE, col.names = column.names)
      m <- data.table(m[,c(1,2, which(colnames(m) == field))])
    } else {
      m <- read.table(path, sep="\t", check.names = FALSE, quote = NULL,
                      header = TRUE, stringsAsFactors=FALSE)
      if (! field %in% colnames(m)) {
        stop("The specified field cannot be found!")
      }
      m <- data.table(m[,c(1,2, which(colnames(m) == field))])
      setnames(m, "Name", "Transcript")
    }

    setnames(m, field, basename(dirname(path)))
    return(m)
  } else {
    stop(paste(path, "does not exist."))
  }
  return(NULL)
}

#### Split Ensembl field ####
format_multi_ensembl_ids <- function(ids) {
  # Format Ensembl Transcript and Ensembl Gene IDs if there are multiple
  # e.g.
  # ENSMUST00000111043_ENSMUSG00000048482,ENSMUST00000111044_ENSMUSG00000048482_mm9_chr1
  # becomes
  # ENSMUST00000111043,ENSMUST00000111044_ENSMUSG00000048482_mm9_chr1
  # Test regex: https://regex101.com/r/zuDsy1/1
  split_ids <- str_match(ids, "^(([^_]+_[^_,]+)(,[^_]+_[^_,]+)*)_(([^_]+).+)")
  
  if (is.na(split_ids[1])) {
      stop("Unable to format Ensembl ID by regex")
  }
  
  # Separate multiple Transcript_Gene name
  ens <- strsplit(split_ids[,2], ",")
  # Split transcript and gene names, then re-arrange to combine transcripts and genes
  ens <- lapply(ens, function(e) {
    strsplit(e, "_") %>% do.call("rbind", .) %>%
      apply(., 2, function(y) paste(unique(y), collapse=",")) %>%
      paste(., collapse="_")
  })
  stopifnot(length(ens) == length(ids))
  apply(cbind(ens, split_ids[,5]), 1, paste, collapse="_")
}

separate_ensembl_field <- function(df) {
  tx_pattern <- "^([^_]+_[^_,]+)(,[^_]+_[^_,]+)*_[^_]+_(chr)*\\w+_\\d+_\\d+_[-+]_utr_\\d+_\\d+"
  if (grepl(tx_pattern, df$Transcript[1], perl = TRUE)) {
    # Format Ensembl Transcript and Ensembl Gene IDs if there are multiple
    # Remove utr tag
    df[, Transcript := str_extract(Transcript, tx_pattern) %>%
           format_multi_ensembl_ids() %>%
           #str_replace("_(hg\\d+|mm\\d+|unk)", "") %>%
           str_replace("_utr", "")]

    # Split by underscore
    df[, c("Transcript", "Gene", "Species", "Chr", "LastExon.Start",
        "LastExon.End", "Strand", "UTR3.Start", "UTR3.End") :=
         tstrsplit(Transcript, "_", fixed=TRUE)]
    df[, Species := NULL]

    df[, ':=' (
      LastExon.Start = as.numeric(as.character(LastExon.Start)),
      LastExon.End = as.numeric(as.character(LastExon.End)),
      UTR3.Start = as.numeric(as.character(UTR3.Start)),
      UTR3.End = as.numeric(as.character(UTR3.End))
    )]
    df[, Length := abs(UTR3.End - UTR3.Start)]
  } else {
    warning("Unable to find Ensembl IDs by regex")
    df
  }
}

#### Add Ensembl Gene ID column ####
extract_one_transcript <- function(ids) {
  str_extract(ids, "^[^,]+")
}

add_ensembl_metadata <- function(df, dbfile, all_genes = FALSE,
                                 non_standard = FALSE) {
  # Add additional Ensembl metadata
  # df - merged_data frame
  # dbfile - path of Ensembl metadata file

  db <- read.table(dbfile, header = TRUE, sep = "\t", stringsAsFactors = FALSE, 
                   quote = NULL) %>%
    data.table()

  df[, tid := extract_one_transcript(Transcript)]
  
  if (!all_genes) {
      db <- db[Gene.type == "protein_coding"]
  }
  gid <- db[, .(tid=Transcript.stable.ID,
                Gene=Gene.stable.ID,
                Gene_Name=Gene.name)] %>%
      unique()
  if ("Gene" %in% colnames(df)) {
      df[, Gene := NULL]   # Remove Gene column so it doesn't conflict with gid results
  }
  setkey(df, tid)
  setkey(gid, tid)
  df <- gid[df]
  df[, tid := NULL]

  c <- length(which(!is.na(df$Gene_Name)))
  pct <- c/nrow(df)*100
  write(sprintf("Found %d/%d (%0.2f%%) matches", c, nrow(df), pct), stderr())
  if (c == 0) {
    warning("No annotation matches were found. Are you using the correct database?")
  }
  if (pct < 90.0) {
    warning("A large proportion of events could not be matched to the annotation.",
            " This could suggest a possible issue with the database used.",
            " In a typical run, all of the events have a match since",
            " the 3' UTR library should have been built using the same database file.")
  }
  if (!non_standard) {
    meta_cols <- c("Transcript", "Gene", "Gene_Name", "Chr",
                   "LastExon.Start", "LastExon.End", "Strand", "UTR3.Start",
                   "UTR3.End", "Length")
    setcolorder(df, c(meta_cols, 
                      sort(colnames(df)[which(!colnames(df) %in% meta_cols)])))
  }
  return(unique(df))
}
################################################################################

#### Join data ####
write(paste("Merging samples by", opt$options$field), stderr())
merged_data <- join_iterative(opt$args, by = c("Transcript", "Length"),
                              field = opt$options$field,
                              format = opt$options$format)

# Remove random chromosomes
# if (!opt$options$merge_only && all(grepl("chr[0-9XY]+", merged_data$Transcript))) {
#   write("Removing random chromosomes", stderr())
#   merged_data <- merged_data[grep("chr[0-9XY]+_\\d+_\\d+", Transcript),]
# }

if (nrow(merged_data) < 100) {
  warning("Less than 100 rows in final table.")
}

if (!opt$options$merge_only) {
  write("Separating Ensembl IDs", stderr())
  separate_ensembl_field(merged_data)
  write("Adding Ensembl metadata", stderr())
  merged_data <- add_ensembl_metadata(merged_data, opt$options$db,
                                      opt$options$all_genes,
                                      opt$options$non_standard)
}

#### Data adjustment ####
# first_sample <- which(colnames(merged_data) == "Length") + 1
# samples_ix <- seq(first_sample, ncol(merged_data))

#### Write output ####
write.table(merged_data, file="", quote=F, row.names=F, sep="\t")
write("\nFinished merging data", stderr())
