#!/usr/bin/env Rscript
#
# Copyright (c) 2014-2016 Kevin Ha
# University of Toronto

suppressPackageStartupMessages(library(optparse))

args <- commandArgs(TRUE)

option.list <- list(
  make_option(c("--ensembl"), type="character", default=NULL,
              help="Ensembl identifiers database file [%default]"),
  make_option(c("-f", "--field"), type="character",
              default="TPM", help="Field to merge [%default]"),
  make_option(c("-m", "--merge_only"), action="store_true", default=FALSE,
              help="Perform merge of fields only. Don't get additional metadata
              [%default]")
  )
desc <- paste("\nMerge multlple quantification runs into a single summary table.",
              "Output is sent to STDOUT.")
parser <- OptionParser(option_list=option.list,
                       description = desc,
                       usage="usage: %prog [options] <quant directories to merge>")
opt <- parse_args(parser, args=args, positional_arguments=TRUE)

if (length(opt$args) < 1) {
  write("Missing arguments", stderr())
  stop(print_help(parser))
}

if (opt$options$merge_only) {
  write("Merge-only mode enabled.", stderr())
}

if (grepl(".db$", opt$options$ensembl)) {
  warning(paste("Did you supply a sqlite3 database to --ensembl?",
                 "We now use a tab-delimited file now!"))
}

suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dplyr))
library(assertthat)

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
    assert_that(all(x$Transcript == y$Transcript) == TRUE)
    x <- y[x]
    assert_that(all(x$Transcript == y$Transcript) == TRUE)
  }
  setTxtProgressBar(pb, 1)
  close(pb)
  write("\n", stderr())
  return(x)
}

load_data <- function(path, field = "TPM") {
  # Load data and keep only the selected field

  if (file.exists(path)) {
    column.names <- c("Transcript", "Length", "TPM", "NumReads")

    if (! field %in% column.names) {
      stop("The specified field cannot be found!")
    }

    m <- read.table(path, sep="\t", check.names = FALSE, stringsAsFactors=FALSE)
    stopifnot(ncol(m) == length(column.names))
    colnames(m) <- column.names
    m <- data.table(m[,c(1,2, which(column.names == field))])
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
  split_ids <- str_match(ids, "^(ENS.+)_((hg19|mm9|mm10).+)")
  # Separate multiple Transcript_Gene name
  ens <- strsplit(split_ids[,2], ",")
  # Split transcript and gene names, then re-arrange to combine transcripts and genes
  ens <- lapply(ens, function(e) {
    strsplit(e, "_") %>% do.call("rbind", .) %>%
      apply(., 2, function(y) paste(unique(y), collapse=",")) %>%
      paste(., collapse="_")
  })
  stopifnot(length(ens) == length(ids))
  apply(cbind(ens, split_ids[,3]), 1, paste, collapse="_")
}

separate_ensembl_field <- function(df) {
  tx_pattern <- "^ENS(MUS)*T.*_(hg\\d+|mm\\d+)_chr[0-9XY]+_\\d+_\\d+_[-+]_utr_\\d+_\\d+"
  if (grepl(tx_pattern, df$Transcript[1], perl = TRUE)) {
    # Format Ensembl Transcript and Ensembl Gene IDs if there are multiple
    # Remove hg19 info
    # Remove utr tag
    df[, Transcript := str_extract(Transcript, tx_pattern) %>%
           format_multi_ensembl_ids() %>%
           str_replace("_(hg\\d+|mm\\d+)", "") %>%
           str_replace("_utr", "")]

    # Split by underscore
    df[, c("Transcript", "Gene", "Chr", "LastExon.Start", "LastExon.End",
           "Strand", "UTR3.Start", "UTR3.End") :=
         tstrsplit(Transcript, "_", fixed=TRUE)]

    df[, ':=' (
      LastExon.Start = as.numeric(as.character(LastExon.Start)),
      LastExon.End = as.numeric(as.character(LastExon.End)),
      UTR3.Start = as.numeric(as.character(UTR3.Start)),
      UTR3.End = as.numeric(as.character(UTR3.End))
    )]
    df[, Length := abs(UTR3.End - UTR3.Start)]
    #return(df)
  } else {
    stop("Unable to separate Ensembl IDs by regex")
  }
  #return(df)
}

#### Add Ensembl Gene ID column ####
extract_one_transcript <- function(ids) {
  str_extract(ids, "ENS(MUS)*T\\d+")
}

add_ensembl_metadata <- function(df, dbfile = opt$options$ensembl) {
  # Add additional Ensembl metadata
  # df - merged_data frame
  # dbfile - path of Ensembl metadata file

  db <- read.table(dbfile, header = TRUE, sep = "\t", stringsAsFactors=FALSE) %>%
    data.table()

  if (all(grepl("ENS(MUS)*T\\d+.*", df$Transcript[1:1000], perl=TRUE))) {
    df[, tid := extract_one_transcript(Transcript)]

    gid <- db[Gene.type == "protein_coding",
              .(tid=Transcript.stable.ID,
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
  } else {
    # should never enter this condition if input is properly prepared
    stop("Can't identify gene ID column type")
  }

  c <- length(which(!is.na(df$Gene_Name)))
  pc <- paste0("(", round(c / nrow(df)*100, digits=7), "%)")
  write(paste("Found", c, "/", nrow(df), pc, "matches"), stderr())
  if (c == 0) {
    warning("No annotation matches were found. Are you using the correct database?")
  }
  meta_cols <- c("Transcript", "Gene", "Gene_Name", "Chr",
                 "LastExon.Start", "LastExon.End", "Strand", "UTR3.Start",
                 "UTR3.End", "Length")
  setcolorder(df, c(meta_cols, sort(colnames(df)[which(!colnames(df) %in% meta_cols)])))
  return(unique(df))
}
################################################################################

#### Join data ####
write(paste("Merging samples by", opt$options$field), stderr())
merged_data <- join_iterative(opt$args, by = c("Transcript", "Length"),
                              field = opt$options$field)

# Remove random chromosomes
if (!opt$options$merge_only && all(grepl("chr[0-9XY]+", merged_data$Transcript))) {
  write("Removing random chromosomes", stderr())
  merged_data <- merged_data[grep("chr[0-9XY]+_\\d+_\\d+", Transcript),]
}

if (nrow(merged_data) < 100) {
  warning("Less than 100 rows in final table.")
}

if (!opt$options$merge_only) {
  write("Separating Ensembl IDs", stderr())
  separate_ensembl_field(merged_data)
  write("Adding Ensembl metadata", stderr())
  merged_data <- add_ensembl_metadata(merged_data, opt$options$ensembl)
}

#### Data adjustment ####
first_sample <- which(colnames(merged_data) == "Length") + 1
samples_ix <- seq(first_sample, ncol(merged_data))

#### Write output ####
write.table(merged_data, file="", quote=F, row.names=F, sep="\t")
write("\nFinished merging data", stderr())
