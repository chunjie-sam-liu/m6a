# Metainfo ----------------------------------------------------------------

# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: Fri May 12 20:04:34 2023
# @DESCRIPTION: filename

# Library -----------------------------------------------------------------

library(magrittr)
library(ggplot2)
library(patchwork)
library(rlang)

# args --------------------------------------------------------------------


# src ---------------------------------------------------------------------


# header ------------------------------------------------------------------

future::plan(future::multisession, workers = 10)

# function ----------------------------------------------------------------


# load data ---------------------------------------------------------------

anno_filename <- "/scr1/users/liuc9/m6a/merge_output.vcf.avinput.hg38_multianno.txt"


# body --------------------------------------------------------------------
anno <- data.table::fread(
  file = anno_filename,
  sep = "\t"
)


anno |> 
  dplyr::mutate(
    mut = glue::glue("{Ref}>{Alt}")
  ) |> 
  # filter indel
  dplyr::filter(
    Ref %in% c("A", "G", "C", "T")
  ) |> 
  dplyr::filter(
    Alt %in% c("A", "G", "C", "T")
  ) |> 
  # filter snp in dbSNP
  dplyr::filter(
    !grepl(
      pattern = "rs",
      x = avsnp150
    )
  ) |> 
  # filter c-to-t
  dplyr::filter(
    mut %in% c("C>T", "G>A")
  ) ->
  anno_nors

anno_nors |> 
  dplyr::group_by(
    Func.refGene
  ) |> 
  dplyr::count() |> 
  dplyr::ungroup() |> 
  dplyr::arrange(-n)

# anno_nors |> 
#   dplyr::group_by(
#     Func.knownGene
#   ) |> 
#   dplyr::count() |> 
#   dplyr::ungroup() |> 
#   dplyr::arrange(-n)
# 
# anno_nors |> 
#   dplyr::group_by(
#     Func.ensGene
#   ) |> 
#   dplyr::count() |> 
#   dplyr::ungroup() |> 
#   dplyr::arrange(-n)

anno_nors |> 
  dplyr::select(
    Chr, Start, End, 
    Func.refGene,
    Gene.refGene,
    GeneDetail.refGene,
    ExonicFunc.refGene,
    AAChange.refGene,
    mut
  ) ->
  anno_nors_sel


anno_nors_sel |> 
  head() |> 
  View()

anno_nors_sel |> 
  dplyr::group_by(
    Func.refGene
  ) |> 
  dplyr::count() |> 
  dplyr::ungroup() |> 
  dplyr::arrange(-n)

# footer ------------------------------------------------------------------

future::plan(future::sequential)

# save image --------------------------------------------------------------