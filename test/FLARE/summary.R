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
# anno_file <- "/scr1/users/liuc9/m6a/test4m.sorted.md.bam.combined.readfiltered.formatted.varfiltered.snpfiltered.ranked.bed.anno"

# src ---------------------------------------------------------------------


# header ------------------------------------------------------------------

# future::plan(future::multisession, workers = 10)

# function ----------------------------------------------------------------

fn_dis <- function(.start, .strand, .gene_id, .gene_region, .overlapping) {
  # .start <- anno_region$start[[1]]
  # .strand <- anno_region$strand[[1]]
  # .gene_id <- anno_region$gene_id[[1]]
  # .gene_region <- anno_region$gene_region[[1]]
  # .overlapping <-anno_region$overlapping[[1]]
  
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
    dplyr::bind_rows() |> 
    dplyr::mutate(
      X2 = as.numeric(X2),
      X3 = as.numeric(X3)
    ) ->
    .d
  
  colnames(.d) <- c("transcript_id","region_start", "region_stop", "strand", "region", "gene_id", "gene_name", "transcript_type", "gene_type", "feature")
  
  .d |> 
    dplyr::filter(
      region == "transcript"
    ) |> 
    dplyr::mutate(
      transcript_length = abs(region_start - region_stop)
    ) |> 
    dplyr::select(
      transcript_id,
      transcript_length
    ) ->
    .dt
  
  
  .d |> 
    dplyr::filter(
      feature == "feature_contains_query",
      strand == .strand
    ) |> 
    dplyr::filter(
      region == .gene_region
    ) |> 
    dplyr::select(
      transcript_id,
      region_start,
      region_stop,
      strand,
      region
    ) |> 
    dplyr::distinct() |> 
    dplyr::left_join(
      .dt,
      by = "transcript_id"
    ) ->
    .dd
  
  
  .dd |> 
    dplyr::mutate(
      rel_locl = abs(.start - region_start),
      feature_size = abs(region_stop - region_start)
    ) |> 
    dplyr::mutate(
      ratio = rel_locl / feature_size
    ) |> 
    dplyr::mutate(
      ratio = dplyr::case_match(
        region,
        "five_prime_utr" ~ ratio + 0,
        "CDS" ~ ratio + 1,
        "three_prime_utr" ~ ratio + 2
      )
    ) |> 
    dplyr::select(
      transcript_id,
      rel_locl,
      feature_size,
      ratio,
      transcript_length
    )
}

# load data ---------------------------------------------------------------
col_names <- c(
  "chrom",
  "start",
  "end",
  "confidence",
  "edit",
  "strand",
  "gene_id",
  "gene_name",
  "gene_region",
  "overlapping"
)
anno <- readr::read_tsv(
  file=anno_file,
  col_names = col_names,
  col_types = readr::cols(
    edit = "c"
  )
) |> 
  tidyr::separate(
    col = edit,
    into = c("t", "c"),
    sep = ",",
    convert = T
  ) |> 
  dplyr::mutate(
    edit_ct_ratio = t/c
  )



# body --------------------------------------------------------------------



# dis ---------------------------------------------------------------------

anno |> 
  dplyr::filter(
    gene_region %in% c("CDS", "three_prime_utr", "five_prime_utr")
  )  ->
  anno_region

future::plan(future::multisession, workers = 20)
anno_region |> 
  dplyr::mutate(
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
  ) |> 
  dplyr::select(-overlapping) |> 
  tibble::as_tibble() ->
  anno_region_dis
future::plan(future::sequential)


# plot --------------------------------------------------------------------


anno |> 
  # dplyr::filter(
  #   t > 2,
  #   edit_ct_ratio > 0.05
  # ) |>
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
  p_pie;p_pie

ggsave(
  filename = "pie-plot-c-to-u-distribution.pdf",
  plot = p_pie,
  device = "pdf",
  path = "/home/liuc9/github/m6a/test/FLARE/plots",
  width = 9,
  height = 9
)

future::plan(future::multisession, workers = 10)
anno_region_dis |> 
  dplyr::mutate(
    dis = furrr::future_map(
      .x = dis,
      .f = function(.x) {
        .x |> 
          dplyr::slice_max(
            transcript_length
          )
      }
    )
  ) |> 
  tidyr::unnest(cols = dis) ->
  anno_region_dis_unnest
future::plan(future::sequential)


anno_region_dis_unnest |> 
  dplyr::group_by(
    gene_region
  ) |> 
  dplyr::summarise(s = mean(feature_size))

# annotation error
anno_region_dis_unnest |> 
  # dplyr::filter(
  #   # confidence > 0,
  #   t > 2,
  #   edit_ct_ratio > 0.05
  # ) |> 
  ggplot(aes(
    x = ratio
  )) +
  geom_density()


anno_region_dis_unnest |> 
  ggplot(aes(
    x = ratio
  )) +
  geom_density()



# footer ------------------------------------------------------------------

future::plan(future::sequential)

# save image --------------------------------------------------------------