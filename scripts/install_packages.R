#!/usr/bin/env Rscript

writeLines(paste("Using", R.Version()$version.string))

local.lib <- .libPaths()[1]
pkgs <- c("stringr", "dplyr", "data.table", "optparse")

for (p in pkgs) {
  if (!p %in% installed.packages(lib.loc = local.lib)) {
    install.packages(p, dependencies=TRUE,
                     repos='http://cran.us.r-project.org', lib=local.lib)
  }
}

writeLines("\n-------------------------------------")
writeLines("INSTALL SUMMARY")
writeLines("-------------------------------------")
check <- vector(length = length(pkgs))
for (p in pkgs) {
  check[p] <- p %in% installed.packages(lib.loc = local.lib)
  writeLines(sprintf("%s: %s", p, ifelse(check[p] == TRUE, "OK", "FAIL")))  
}
writeLines("-------------------------------------\n")
status <- !(all(check))
q(status=status)
