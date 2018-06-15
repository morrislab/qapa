#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(optparse))

# Options ####
args <- commandArgs(TRUE)

option.list <- list(
    make_option(c("-e", "--expr"), action="store_true", default=FALSE,
        help="Include expression values in output [%default]"),
    make_option(c("-q", "--qualityscores"), action="store_true",
        default=FALSE, help="*EXPERIMENTAL* Include quality scores (e.g.
         standard deviation from beta distribution [%default])")
)
desc <- "Calculate PAU given TPM values. Output written to STDOUT."
parser <- OptionParser(option_list=option.list,
                       description = desc,
                       usage="usage: %prog [options] COUNTS.tab")
opt <- parse_args(parser, args=args, positional_arguments=TRUE)

if (length(opt$args) == 0)
    stop(print_help(parser))

suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(data.table))


# Functions ---------------------------------------------------------------


calculate <- function(dt, qscore) {
  # Calculate UTR Index given a data.table
  #
  # Return as data table or NULL if dt is empty
  write("Calculating Poly(A) Usage", stderr())
  write(paste("\t", nrow(dt), "rows,", length(unique(dt$Gene)), "genes"), stderr())

  if (nrow(dt) == 0) {
      return(NULL)
  }
  
  dt[, c("APA_ID", "pau") :=
            list(update_apa_id(APA_ID, UTR3.Start, UTR3.End),
                 calc_pau(value)
            ),
     by=key(dt)]

  # Cast to wide format and add suffixes
  ui <- dcast(dt, APA_ID + Transcript + Gene + Gene_Name + Chr + LastExon.Start +
                LastExon.End + Strand + UTR3.Start + UTR3.End + Length ~ sample,
              value.var = "pau")
  measured.cols <- (which(colnames(ui) == "Length") + 1):ncol(ui)
  colnames(ui)[measured.cols] <- str_replace_all(colnames(ui)[measured.cols],
                                              "$", ".PAU")

  if (qscore) {
    dt[, c("q") := list(calc_beta_sd(value)), by=key(dt)]
    q <- dcast(dt, APA_ID + Transcript + Gene + Gene_Name + Chr + LastExon.Start +
                 LastExon.End + Strand + UTR3.Start + UTR3.End + Length ~ sample,
               value.var = "q")
    colnames(q)[measured.cols] <- str_replace_all(colnames(q)[measured.cols],
                                                  "$", ".Q")

    return(
      merge(data.table(ui), data.table(q),
            by = c("APA_ID", "Transcript", "Gene", "Gene_Name", "Chr", "LastExon.Start",
                   "LastExon.End", "Strand", "UTR3.Start", "UTR3.End", "Length")
      )
    )
  }
  return(data.table(ui))
}


beta_sd <- function(incl, excl) {
  incl <- incl + 1
  excl <- excl + 1
  s <- (incl*excl) / ((incl+excl)^2*(incl+excl+1))
  return(sqrt(s))
}

calc_pau <- function(x) {
  round(x/sum(x)*100, digits = 3)
}

calc_beta_sd <- function(x) {
  beta <- abs(x - sum(x))
  round(beta_sd(x, beta), digits = 3)
}

apa_id <- function(x) {
  # Generate an APA_ID that combines Ensembl Gene ID with a number
  # Assume that input is sorted by x
  eg.rle <- rle(as.character(x))
  paste0(rep(eg.rle$values, times=eg.rle$lengths), "_",
         unlist(sapply(eg.rle$lengths, seq_len)))
}

create_apa_id <- function(d) {
  # Create custom Ensembl Gene APA IDs
  d.new <- d[order(d$Gene, d$Length),]
  setkey(d.new, d)
  d.new$APA_ID <- apa_id(d.new$Gene)
  return(d.new)
}

update_apa_id <- function(x, start, end) {
  # Update APA_ID with proximal, distal or single status - added as another suffix
  if (length(x) == 1) {
    x <- paste(x, "S", sep = "_")
  } else {
    px <- which.min(abs(end - start))
    x[px] <- paste(x[px], "P", sep = "_")
    x[-px] <- paste(x[-px], "D", sep = "_")
  }
  return(x)
}


get_first_sample_ix <- function(df) {
  return(which(colnames(df) == "Length") + 1)
}

get_num_samples <- function(df) {
  return(ncol(df) - (get_first_sample_ix(df) - 1))
}
###############################################################################

# Input ####
using_stdin <- FALSE
file <- opt$args[1]
if (file == "-") {
  file <- file('stdin')
  using_stdin <- TRUE
} else if (!file.exists(file)) {
  stop(paste("Input file", file, "doesn't exist!"))
}
m <- data.table(read.csv(file, sep="\t", header=T, check.names=FALSE, stringsAsFactors=FALSE))
setkey(m, Gene)
m[, APA_ID := apa_id(Gene)]

# Calculate PAU values ####
write("Melting data frame", stderr())
dt <-  melt(m,
       id.vars = c("APA_ID",
                   colnames(m)[1:(get_first_sample_ix(m) - 1)]),
       variable.name = "sample")


write("Operating on forward strand", stderr())
dt.plus <- dt[dt$Strand == "+"]
setkey(dt.plus, sample, Gene, Chr, LastExon.Start)
#dt.plus[, "label" := find_distances_to_merge(LastExon.End, opt$options$width), by=key(dt.plus)]
#dt.plus <- merge_intervals(dt.plus)
#setorder(dt.plus, Chr, LastExon.Start, LastExon.End)
#setkey(dt.plus, sample, Gene, Chr, LastExon.Start)
wide.plus <- calculate(dt.plus, opt$options$qualityscores)

write("\nOperating on reverse strand", stderr())
dt.minus <- dt[Strand == "-"]
setkey(dt.minus, sample, Gene, Chr, LastExon.End)
#dt.minus[, "label" := find_distances_to_merge(LastExon.Start, opt$options$width), by=key(dt.minus)]
#dt.minus <- merge_intervals(dt.minus)
#setorder(dt.minus, Chr, LastExon.End, -LastExon.Start)
#setkey(dt.minus, sample, Gene, Chr, LastExon.End)
wide.minus <- calculate(dt.minus, opt$options$qualityscores)

pau <- rbind(wide.plus, wide.minus)

# Count the number of events per gene
setkey(pau, Gene_Name)
n_sites <- pau[, .(Num_Events = .N), by=Gene_Name]
pau <- n_sites[pau]

if (opt$options$expr) {
  # Append a suffix (e.g. .TPM) to the end to distinguish from the pau columns
  write("\nAdding input expression values", stderr())
  if (nrow(dt.plus) == 0) {
      dt.plus <- NULL
  }
  if (nrow(dt.minus) == 0) {
      dt.minus <- NULL
  }
  expr <- rbind(dt.plus, dt.minus)
  expr[,sample := paste0(sample, ".TPM")]
  expr <- dcast(expr, Transcript + Gene + Gene_Name + Chr + LastExon.Start +
                  LastExon.End + Strand + UTR3.Start + UTR3.End + Length ~ sample,
                value.var = "value")
  meta <- grep(".TPM$", colnames(expr), value = TRUE, invert = TRUE)
  pau_samples <- grep(".PAU$", colnames(pau), value = TRUE)
  tpm_samples <- str_replace(pau_samples, ".PAU$", ".TPM")
  setcolorder(expr, c(meta, tpm_samples))
  pau <- merge(pau, data.table(expr),
               by = c("Transcript", "Gene", "Gene_Name", "Chr", "LastExon.Start",
                      "LastExon.End", "Strand", "UTR3.Start", "UTR3.End", "Length"))
}
setcolorder(pau, c("APA_ID", grep("APA_ID", colnames(pau), invert = TRUE, value = TRUE)))

# Output ####
write.table(pau, file="", sep="\t", quote=F, row.names=F)
write("\nFinished computing PAU!", stderr())
