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

fn_dis <- function(.start, .strand, .gene_id, .gene_region, .overlapping) {
  # .start <- anno_region$start[[6]]
  # .strand <- anno_region$strand[[6]]
  # .gene_id <- anno_region$gene_id[[6]]
  # .gene_region <- anno_region$gene_region[[6]]
  # .overlapping <-anno_region$overlapping[[6]]
  
  strsplit(
    x = .overlapping,
    split = "\\|"
  )[[1]] |> 
    strsplit(
      split = ":"
    ) |> 
    purrr::map(.f = function(.x) {
      data.frame(t(.x))
    }) |> 
    dplyr::bind_rows() ->
    .d
  
  colnames(.d) <- c("transcript_id","region_start", "region_stop", "strand", "region", "gene_id", "gene_name", "transcript_type", "gene_type", "feature")
  
  
  .d |> 
    dplyr::filter(
      feature == "feature_contains_query",
      # transcript_type == "protein_coding"
    ) |> 
    dplyr::filter(
      strand == .strand,
      region == .gene_region
      # region == "transcript"
    ) |> 
    dplyr::mutate(
      region_start = as.numeric(region_start),
      region_stop = as.numeric(region_stop)
    ) ->
    .dd
  
  
  .dd |> 
    dplyr::mutate(
      dis = abs(.start - region_start)
    ) |> 
    dplyr::mutate(
      ratio = dis / abs(region_stop - region_start)
    )

  
  # .dd
  # if(.strand == "+") {
  #   .r1 <- min(.dd$region_start)
  #   .r2 <- max(.dd$region_stop)
  #   # .dis <- .start - .r1
  # } else {
  #   .r1 <- max(.dd$region_start)
  #   .r2 <- min(.dd$region_stop)
  #   # .dis <- .start - .r1
  # }
  # 
  # .dis <- abs(.start - .r1) / abs(.r1 - .r2)
  # .dis
}

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
  dplyr::arrange(-n) |> 
  dplyr::mutate(
    gene_region = factor(
      x = gene_region,
      levels = gene_region
    )
  ) |> 
  dplyr::mutate(csum = rev(cumsum(rev(n)))) |> 
  dplyr::mutate(pos = n/2 + dplyr::lead(csum, 1)) |> 
  dplyr::mutate(pos = dplyr::if_else(is.na(pos), n/2, pos)) %>% 
  dplyr::mutate(percentage = n/sum(n)) |> 
  ggplot(aes(x = "", y = n)) +
  geom_bar(
    aes(fill = gene_region),
    stat = "identity", 
    width = 1, 
    color = "white"
  ) + 
  ggsci::scale_fill_npg() +
  ggrepel::geom_label_repel(
    aes(
      y = pos,
      label = glue::glue("{gene_region}\n{n} ({scales::percent(percentage)})"),
      fill = gene_region,
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
  p_pie

ggsave(
  filename = "pie-plot-c-to-u-distribution.pdf",
  plot = p_pie,
  device = "pdf",
  path = "/home/liuc9/github/m6a/test/FLARE/plots",
  width = 9,
  height = 9
)


# dis ---------------------------------------------------------------------

anno |> 
  dplyr::filter(
    gene_region %in% c("CDS", "three_prime_utr", "five_prime_utr")
  )  ->
  anno_region

future::plan(future::multisession, workers = 10)
anno_region |> 
  dplyr::mutate(
    # dis = purrr::pmap(
    dis = furrr::future_pmap(
      .l = list(
        .start = start,
        .strand = strand,
        .gene_id = gene_id,
        .gene_region = gene_region,
        .overlapping = overlapping
      ),
      .f = fn_dis
    )
  ) ->
  anno_region_dis
future::plan(future::sequential)

# anno_region_dis |> 
#   dplyr::filter(gene_region == "CDS") ->
#   .m
# .m$dis |> summary()

anno_region_dis |> 
  # dplyr::filter(gene_region == "CDS") |> 
  as.data.frame() |> 
  ggplot(aes(
    x = dis
  )) +
  geom_histogram()

# anno_region_dis |> 
#   dplyr::select(gene_region, dis) |> 
#   ggplot(aes(x = dis, y = gene_region)) +
#   ggridges::geom_density_ridges()
#   geom_boxplot()

  
anno_region_dis |> 
  dplyr::select(gene_region, dis) |> 
  tibble::as_tibble() |> 
  tidyr::unnest(cols = dis) |> 
  dplyr::mutate(
    ratio_new = dplyr::case_match(
      gene_region,
      "five_prime_utr" ~ ratio + 0,
      "CDS" ~ ratio + 1,
      "three_prime_utr" ~ ratio + 2
    )
  ) ->
  anno_region_dis_new

anno_region_dis_new |> 
  as.data.frame() |> 
  # dplyr::filter(gene_region != "three_prime_utr") |> 
  ggplot(aes(
    x = ratio_new
  )) +
  geom_density()




# footer ------------------------------------------------------------------

future::plan(future::sequential)

# save image --------------------------------------------------------------