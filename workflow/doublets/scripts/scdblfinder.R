# scDblFinder on each batch of a group of batches, called by scdblfinder.py as
#   Rscript --vanilla scdblfinder.R <params.json>
# One scDblFinder() call per batch, with the random seed set before each call. Writes the score
# and the call (1 = doublet) of every cell, in the order of the input cells.
suppressPackageStartupMessages({
  library(scDblFinder)
  library(SingleCellExperiment)
  library(Matrix)
  library(hdf5r)
  library(jsonlite)
})

params <- fromJSON(commandArgs(trailingOnly = TRUE)[1])

h5 <- H5File$new(params$input, mode = "r")
read <- function(name) h5[[name]][]
shape <- read("counts/shape")
genes <- make.unique(read("var_names"))
# cells x genes CSR in the file = genes x cells CSC in R
counts <- new("dgCMatrix", i = read("counts/indices"), p = read("counts/indptr"),
              x = read("counts/data"), Dim = c(shape[2], shape[1]),
              Dimnames = list(genes, paste0("cell", seq_len(shape[1]))))
batch <- read("batch")
batch_names <- read("batch_names")
h5$close_all()
stopifnot(length(genes) == shape[2], length(batch) == shape[1])
cat(ncol(counts), "cells,", nrow(counts), "genes,", length(batch_names), "batches\n")

score <- rep(NA_real_, ncol(counts))
doublet <- rep(NA_integer_, ncol(counts))
for (code in sort(unique(batch))) {
  name <- batch_names[code + 1]
  cells <- which(batch == code)
  sce <- SingleCellExperiment(list(counts = counts[, cells, drop = FALSE]))
  # a wrapper keeps the object out of the call that do.call() builds
  run <- function(...) scDblFinder(sce, returnType = "sce", ...)
  set.seed(params$seed)
  sce <- tryCatch(
    do.call(run, as.list(params$args)),
    error = function(e) stop("scDblFinder failed on batch '", name, "': ", conditionMessage(e))
  )
  if (is.null(sce$scDblFinder.class)) stop("scDblFinder returned no doublet calls (threshold = FALSE?)")
  score[cells] <- sce$scDblFinder.score
  doublet[cells] <- as.integer(sce$scDblFinder.class == "doublet")
  threshold <- metadata(sce)$scDblFinder.threshold
  cat(sprintf("batch '%s': %d cells, %d doublets (%.1f%%), score threshold %s\n",
              name, length(cells), sum(doublet[cells]), 100 * mean(doublet[cells]),
              if (is.null(threshold)) "NA" else format(threshold, digits = 3)))
  rm(sce)
  gc(verbose = FALSE)
}
stopifnot(!anyNA(score), !anyNA(doublet))

h5 <- H5File$new(params$output, mode = "w")
h5[["score"]] <- score
h5[["doublet"]] <- doublet
h5$close_all()
