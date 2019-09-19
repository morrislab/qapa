#!/usr/bin/env Rscript

# https://stackoverflow.com/a/29132294
source_funcs <- function(x) {
    cmds <- parse(x)
    assign.funs <- sapply(cmds, function(x) {
        if(x[[1]]=="<-") {
            if(x[[3]][[1]]=="function") {
                return(TRUE)
            }
        }
        return(FALSE)
    })
    return(cmds[assign.funs])
}

files <- c("../scripts/create_merged_data.R",
           "../scripts/compute_pau.R")
for (f in files) {
    cmds <- source_funcs(f)
    eval(cmds)
}

suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dplyr))
library(testthat)
test_dir("R/")
