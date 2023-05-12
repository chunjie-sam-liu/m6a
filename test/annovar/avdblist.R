# Metainfo ----------------------------------------------------------------

# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: Fri May 12 15:31:11 2023
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


# body --------------------------------------------------------------------
hg38dblist <- "/mnt/isilon/xing_lab/liuc9/refdata/annovar/humandb/hg38_avdblist.txt"



db <- readr::read_tsv(
  file = hg38dblist,
  col_names = c("filename", "date", "records")
)

db |> 
  dplyr::filter(!grepl(
    pattern = "idx",
    x = filename
  )) |> 
  dplyr::mutate(date = as.character(date)) |> 
  dplyr::mutate(date = lubridate::as_date(date)) |> 
  dplyr::arrange(dplyr::desc(date)) |> 
  dplyr::mutate(
    name = gsub(
      pattern = "hg38_|\\..*$",
      replacement = "",
      x = filename
    )
  ) ->
  adb

adb |> 
  dplyr::arrange(
    name,
    dplyr::desc(date),
  ) |> 
  print(n = Inf)


# footer ------------------------------------------------------------------

future::plan(future::sequential)

# save image --------------------------------------------------------------
