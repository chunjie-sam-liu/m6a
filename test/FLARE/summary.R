# Metainfo ----------------------------------------------------------------

# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: Mon May 15 16:52:42 2023
# @DESCRIPTION: filename

# Library -----------------------------------------------------------------

library(magrittr)
library(ggplot2)
library(patchwork)
library(rlang)

# args --------------------------------------------------------------------


anno_file <- "/mnt/isilon/xing_lab/liuc9/projdata/m6a/flare/sailor/APOBEC1YTH.sorted.md.bam.combined.readfiltered.formatted.varfiltered.snpfiltered.ranked.bed.anno"

# src ---------------------------------------------------------------------


# header ------------------------------------------------------------------

# future::plan(future::multisession, workers = 10)

# function ----------------------------------------------------------------


# load data ---------------------------------------------------------------
col_names <- c(
  "chrom",
  "start",
  "end",
  "name",
  "score",
  "strand",
  "gene_id",
  "gene_name",
  "gene_region",
  "overlapping"
)
anno <- data.table::fread(
  file=anno_file,
  sep="\t",
  col.names = col_names
)

# body --------------------------------------------------------------------

anno |> 
  dplyr::group_by(gene_region) |> 
  dplyr::count() |> 
  dplyr::ungroup() |> 
  dplyr::arrange(-n)

# footer ------------------------------------------------------------------

future::plan(future::sequential)

# save image --------------------------------------------------------------