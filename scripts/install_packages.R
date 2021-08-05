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
writeLines("INSTALL SUMMARY OF REQUIRED R PACKAGES")
writeLines("-------------------------------------")
check <- vector()
for (p in pkgs) {
  status <- p %in% installed.packages(lib.loc = local.lib)
  check <- append(check, status)
  writeLines(sprintf("%s: %s", p, ifelse(status == TRUE, "OK", "FAIL")))
}
writeLines("-------------------------------------\n")
status <- !(all(check))
q(status=status)
