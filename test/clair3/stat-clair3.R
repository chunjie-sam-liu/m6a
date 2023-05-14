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

pcc <- readr::read_tsv(file = "https://raw.githubusercontent.com/chunjie-sam-liu/chunjie-sam-liu.life/master/public/data/pcc.tsv")

# header ------------------------------------------------------------------

future::plan(future::multisession, workers = 10)

# function ----------------------------------------------------------------


# load data ---------------------------------------------------------------

anno_filename <- "/mnt/isilon/xing_lab/liuc9/projdata/m6a/clair3/merge_output.vcf.gz.avinput.hg38_multianno.txt"


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

# 
# anno_nors_sel |> 
#   head() |> 
#   View()

anno_nors_sel |> 
  dplyr::group_by(
    Func.refGene
  ) |> 
  dplyr::count() |> 
  dplyr::ungroup() |> 
  dplyr::mutate(
    Func.refGene = factor(
      Func.refGene,
      levels = Func.refGene
    )
  ) |> 
  dplyr::mutate(csum = rev(cumsum(rev(n)))) |> 
  dplyr::mutate(pos = n/2 + dplyr::lead(csum, 1)) |> 
  dplyr::mutate(pos = dplyr::if_else(is.na(pos), n/2, pos)) %>% 
  dplyr::mutate(percentage = n/sum(n)) |> 
  ggplot(aes(x = "", y = n, fill = Func.refGene)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  scale_fill_manual(
    name = NULL,
    values = pcc$color,
  ) +
  # scale_color_manual(
  #   name = NULL,
  #   values = pcc$color
  # ) +
  ggrepel::geom_label_repel(
    aes(
      y = pos,
      # label = glue::glue("{Func.refGene}\n{n} ({scales::percent(percentage)})"), 
      label = Func.refGene,
      fill = Func.refGene,
    ),
    size = 6,
    # fill = "white",
    nudge_x = 1,
    show.legend = FALSE,
  ) +
  coord_polar(theta = "y", start = 0) +
  theme_void() +
  theme(
    plot.title = element_text(
      # vjust = -2,
      hjust = 0.5,
      size = 22,
    ),
    legend.position = "none"
  )


anno_nors_sel |> nrow()

anno_nors_sel |> 
  dplyr::filter(
    !grepl(
      pattern = ";",
      x = GeneDetail.refGene
    )
  ) |> 
  dplyr::mutate(
    dist = gsub(
      pattern = "dist=",
      replacement = "",
      x = GeneDetail.refGene
    )
  ) |> 
  # dplyr::filter(dist != ".") |> 
  dplyr::mutate(dis = as.numeric(dist)) |> 
  dplyr::filter(!is.na(dis)) ->
  a

a |> 
  dplyr::filter(mut  == "C>T") |>
  ggplot(aes(x = dis)) +
  geom_density(fill = "blue", alpha = 0.5)


flare_anno <- "/scr1/users/liuc9/m6a/test4m.sorted.md.bam.combined.readfiltered.formatted.varfiltered.snpfiltered.ranked.bed.anno"

d <- readr::read_tsv(
  file = flare_anno,
  col_names = F
)



d |> 
  dplyr::group_by(
    X9
  ) |> 
  dplyr::count() |> 
  dplyr::ungroup() |> 
  dplyr::arrange(-n) |> 
  dplyr::mutate(
    X9 = factor(
      X9,
      levels = X9
    )
  ) |> 
  dplyr::mutate(csum = rev(cumsum(rev(n)))) |> 
  dplyr::mutate(pos = n/2 + dplyr::lead(csum, 1)) |> 
  dplyr::mutate(pos = dplyr::if_else(is.na(pos), n/2, pos)) %>% 
  dplyr::mutate(percentage = n/sum(n)) |> 
  ggplot(aes(x = "", y = n, fill = X9)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  scale_fill_manual(
    name = NULL,
    values = pcc$color,
  ) +
  # scale_color_manual(
  #   name = NULL,
  #   values = pcc$color
  # ) +
  ggrepel::geom_label_repel(
    aes(
      y = pos,
      # label = glue::glue("{X9}\n{n} ({scales::percent(percentage)})"), 
      label = X9,
      fill = X9,
    ),
    size = 6,
    # fill = "white",
    nudge_x = 1,
    show.legend = FALSE,
  ) +
  coord_polar(theta = "y", start = 0) +
  theme_void() +
  theme(
    plot.title = element_text(
      # vjust = -2,
      hjust = 0.5,
      size = 22,
    ),
    legend.position = "none"
  )


d |> 
  dplyr::filter(
    X9 %in% c("five_prime_utr", "CDS", "three_prime_utr")
  ) ->
  dd

dd


# footer ------------------------------------------------------------------

future::plan(future::sequential)

# save image --------------------------------------------------------------