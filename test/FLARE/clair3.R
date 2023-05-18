# Metainfo ----------------------------------------------------------------

# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: Wed May 17 16:50:16 2023
# @DESCRIPTION: 


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

anno_file <- "/mnt/isilon/xing_lab/liuc9/projdata/m6a/clair3/merge_output.vcf.gz.avinput.hg38_multianno.txt"

anno <- vroom::vroom(
  file = anno_file,
) |> 
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
  dplyr::filter(
    !grepl(
      pattern = "rs",
      x = avsnp150
    )
  ) |>
  dplyr::filter(
    mut %in% c("C>T", "G>A")
  ) |> 
  dplyr::select(
    Chr, Start, End,
    Func.refGene,
    Gene.refGene,
    GeneDetail.refGene,
    ExonicFunc.refGene,
    AAChange.refGene,
    mut
  )


# body --------------------------------------------------------------------

anno |> 
  dplyr::group_by(Func.refGene) |> 
  dplyr::count() |> 
  dplyr::ungroup() |> 
  dplyr::arrange(-n) |> 
  dplyr::mutate(
    gene_region = factor(
      x = Func.refGene,
      levels = Func.refGene
    )
  ) |> 
  dplyr::mutate(csum = rev(cumsum(rev(n)))) |> 
  dplyr::mutate(pos = n/2 + dplyr::lead(csum, 1)) |> 
  dplyr::mutate(pos = dplyr::if_else(is.na(pos), n/2, pos)) %>% 
  dplyr::mutate(percentage = n/sum(n)) |> 
  ggplot(aes(x = "", y = n)) +
  geom_bar(
    aes(fill = Func.refGene),
    stat = "identity", 
    width = 1, 
    color = "white"
  ) + 
  # ggsci::scale_fill_npg() +
  ggrepel::geom_label_repel(
    aes(
      y = pos,
      label = glue::glue("{Func.refGene}\n{n} ({scales::percent(percentage)})"),
      fill = Func.refGene,
    ),
    color = "white",
    size = 5,
    # angle = 45,
    # alpha = 0.5,
    fontface = "bold",
    lineheight = 1.2,
    # hjust = 0.3,
    # vjust = 1,
    segment.color = "black",
    nudge_x = 0.8,
    nudge_y = 1,
    box.padding = 0.5,
    # segment.curvature = 0.1,
    # segment.angle = 20,
    max.overlaps = Inf,
    show.legend = FALSE,
    # direction = "y"
    force = 0.2
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
  ) ->
  p_pie;p_pie

ggsave(
  filename = "pie-plot-c-to-u-distribution.pdf",
  plot = p_pie,
  device = "pdf",
  path = "/home/liuc9/github/m6a/test/FLARE/plots",
  width = 9,
  height = 9
)


# footer ------------------------------------------------------------------

future::plan(future::sequential)

# save image --------------------------------------------------------------

