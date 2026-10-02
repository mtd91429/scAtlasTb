# STACAS integration, called by stacas.py as
#   Rscript --vanilla stacas.R <params.json>
# Splits the object by batch and runs Run.STACAS(): anchors between pairs of batches (on the
# integration features, or on genes STACAS selects itself), rescored by STACAS and, in
# semi-supervised mode, filtered by cell type labels, then integrated along STACAS's guide tree.
# Returns the PCA that Run.STACAS computes on the integrated data (`dims` dimensions).
params <- jsonlite::fromJSON(commandArgs(trailingOnly = TRUE)[1])
source(params$utils)
suppressPackageStartupMessages(library(STACAS))
set.seed(params$seed)

obj <- read_seurat_input(params$input, scale_factor = params$scale_factor)
cells <- colnames(obj)
cat(ncol(obj), "cells,", nrow(obj), "genes,", length(unique(obj$batch)), "batches\n")

objects <- SplitObject(obj, split.by = "batch")
rm(obj)
invisible(gc())

# Run.STACAS leaves out batches below min.sample.size, which would drop their cells
min_size <- params$stacas_args$min.sample.size
if (is.null(min_size)) min_size <- formals(Run.STACAS)$min.sample.size
small <- sapply(objects, ncol) < min_size
if (any(small)) {
  stop("STACAS does not integrate batches with fewer than min_sample_size = ", min_size,
       " cells: ", paste(names(objects)[small], collapse = ", "),
       ". Lower min_sample_size or remove these batches.")
}

# default: the exchanged genes, i.e. the integration features; a number: STACAS selects its own
anchor_features <- params$anchor_features
if (is.null(anchor_features)) {
  anchor_features <- rownames(objects[[1]])
  cat("anchor features:", length(anchor_features), "integration features\n")
} else {
  cat("anchor features: STACAS selects", anchor_features, "of", nrow(objects[[1]]), "genes\n")
}

reference <- NULL
if (length(params$reference) > 0) {
  reference <- match(params$reference, names(objects))
  if (anyNA(reference)) {
    stop("reference batches not found: ", paste(params$reference[is.na(reference)], collapse = ", "))
  }
}

cell_labels <- NULL
if (isTRUE(params$semi_supervised)) {
  cell_labels <- "cell_labels"
  n <- sum(sapply(objects, function(o) sum(!is.na(o$cell_labels))))
  cat("semi-supervised:", n, "of", length(cells), "cells labelled\n")
}

# a wrapper keeps the objects out of the call that do.call() builds
run <- function(...) {
  Run.STACAS(objects, dims = params$dims, anchor.features = anchor_features,
             reference = reference, cell.labels = cell_labels, seed = params$seed,
             verbose = TRUE, ...)
}
integrated <- do.call(run, as.list(params$stacas_args))
cat("integration features:", length(rownames(integrated[["integrated"]])), "\n")

embedding <- Embeddings(integrated, "pca")
stopifnot(setequal(rownames(embedding), cells))
embedding <- embedding[cells, , drop = FALSE]
write_embedding(params$output, embedding)
cat("embedding:", nrow(embedding), "cells x", ncol(embedding), "dimensions\n")
