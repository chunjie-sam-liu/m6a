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
library(plyranges)
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg38)

# args --------------------------------------------------------------------


anno_file <- "/mnt/isilon/xing_lab/liuc9/projdata/m6a/flare/sailor/APOBEC1YTH.sorted.md.bam.combined.readfiltered.formatted.varfiltered.snpfiltered.ranked.bed.anno"
# anno_file <- "/scr1/users/liuc9/m6a/test4m.sorted.md.bam.combined.readfiltered.formatted.varfiltered.snpfiltered.ranked.bed.anno"

gtf <- readr::read_rds(
  "/mnt/isilon/xing_lab/liuc9/refdata/ensembl/Homo_sapiens.GRCh38.107.gtf.plyranges.rds"
)

# src ---------------------------------------------------------------------


# header ------------------------------------------------------------------

# future::plan(future::multisession, workers = 10)

# function ----------------------------------------------------------------

fn_extract_gtf <- function(.tid, .rst, .rsp, .start) {
  # .tid <- .dd$transcript_id[[1]]
  # .rst <- .dd$region_start[[1]]
  # .rsp <- .dd$region_stop[[1]]
  
  gtf |> 
    plyranges::filter(
      type == "CDS",
      transcript_id == .tid
    ) |> 
    plyranges::select(exon_number) |> 
    tibble::as_tibble() |> 
    dplyr::arrange(exon_number) ->
    .cdss
  
  .std <- .cdss$strand[[1]] |> as.character()
  
  # .cdss |>
  #   tibble::rowid_to_column(var = "y") |>
  #   ggplot(aes(
  #     x = start,
  #     xend = end,
  #     y = y,
  #     yend = y,
  #     color = exon_number,
  #     label = exon_number
  #   )) +
  #   geom_segment() +
  #   geom_text(
  #     hjust = 'outside', nudge_x = -0.2
  #   ) +
  #   geom_vline(xintercept = .start) +
  #   geom_vline(xintercept = .rst, color = "red") +
  #   geom_vline(xintercept = .rsp, color = "blue")
  
  .cdss_n <- if(.std == "+") {
    .cdss |> 
      dplyr::mutate(
        in_cds = end > .start & start < .start,
        pre_cds = .start > end
      ) |> 
      dplyr::mutate(
        pre_in_cds = in_cds | pre_cds
      ) |> 
      dplyr::mutate(
        dis1 = .start - start
      ) |> 
      dplyr::mutate(
        dis1 = ifelse(
          pre_in_cds == FALSE,
          0,
          dis1
        )
      ) |> 
      dplyr::mutate(
        dis2 = ifelse(
          pre_cds == TRUE,
          width,
          dis1
        )
      ) 
  } else {
    .cdss |> 
      dplyr::mutate(
        in_cds = end > .start & start < .start,
        pre_cds = start > .start
      ) |> 
      dplyr::mutate(
        pre_in_cds = in_cds | pre_cds
      ) |> 
      dplyr::mutate(
        dis1 = end - .start
      ) |> 
      dplyr::mutate(
        dis1 = ifelse(
          pre_in_cds == FALSE,
          0,
          dis1
        )
      ) |> 
      dplyr::mutate(
        dis2 = ifelse(
          pre_cds == TRUE,
          width,
          dis1
        )
      ) 
  }
  
  
  .dis = sum(.cdss_n$dis2)
  .feature_size = sum(.cdss_n$width)
  .ratio = .dis / .feature_size
  
  tibble::tibble(
    rel_loc = .dis,
    feature_size = .feature_size,
    ratio = .ratio
  )
  
}

fn_dis <- function(.start, .strand, .gene_id, .gene_region, .overlapping) {
  
  # mm <- anno_region |>
  #   dplyr::filter(
  #     strand == "+",
  #     gene_region == "CDS"
  #   )
  # .start <- mm$start[[1]]
  # .strand <- mm$strand[[1]]
  # .gene_id <- mm$gene_id[[1]]
  # .gene_region <- mm$gene_region[[1]]
  # .overlapping <-mm$overlapping[[1]]
  
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
  
  # .d |>
  #   dplyr::slice(1:3) |> 
  #   tibble::rowid_to_column(var = "y") |>
  #   ggplot(aes(
  #     x = region_start,
  #     xend = region_stop,
  #     y = y,
  #     yend = y,
  #     color = region,
  #     label = region
  #   )) +
  #   geom_segment() +
  #   geom_text(
  #     hjust = 'outside', nudge_x = -0.2
  #   )

  
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
    
  if(.gene_region == "CDS") {
    
    .dd |> 
      dplyr::mutate(start = .start) |> 
      dplyr::mutate(
        a = purrr::pmap(
          .l = list(
            .tid = transcript_id,
            .rst = region_start,
            .rsp = region_stop,
            .start = start
          ),
          .f = fn_extract_gtf
        )
      ) |> 
      tidyr::unnest(cols = a) |> 
      dplyr::mutate(
        ratio_new = ratio + 1
      ) |> 
      dplyr::select(
        transcript_id,
        rel_loc,
        feature_size,
        ratio,
        ratio_new,
        transcript_length
      )
    
  } else {
    .ddd <- if(.strand == "+") {
      .dd |> 
        dplyr::mutate(
          rel_loc = abs(.start - region_start)
        )
    } else {
      .dd |> 
        dplyr::mutate(
          rel_loc = abs(.start - region_stop)
        )
    }
    
    .ddd |> 
      dplyr::mutate(
        feature_size = abs(region_stop - region_start)
      ) |> 
      dplyr::mutate(
        ratio = rel_loc / feature_size
      ) |> 
        dplyr::mutate(
          ratio_new = dplyr::case_match(
            region,
            "five_prime_utr" ~ ratio + 0,
            "CDS" ~ ratio + 1,
            "three_prime_utr" ~ ratio + 2
          )
        ) |> 
        dplyr::select(
          transcript_id,
          rel_loc,
          feature_size,
          ratio,
          ratio_new,
          transcript_length
        )
  }
  
  
  
  
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

future::plan(future::multisession, workers = 50)
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

anno_region_dis_unnest |> 
  # dplyr::filter(
  #   t > 1
  # ) |> 
  ggplot(aes(
    x = ratio_new
  )) +
  geom_histogram(
    aes(y = ..density..),
    fill = "blue",
    color = "white",
    alpha = 0.3,
    linewidth = 0.1,
    bins = 30,
    binwidth = 0.1,
    boundary=0
  ) +
  geom_density(
    alpha = 0.1,
    linetype = 2,
    linewidth = 1.2,
    colour = "red",
    fill = "red"
  ) +
  geom_vline(
    xintercept = c(1, 2),
    color = "red",
    linetype = 3,
    linewidth = 1.2,
  ) +
  scale_x_continuous(
    expand = expansion(mult = c(0.005), add = 0),
    limits = c(0, 3),
    breaks = seq(0, 3, 0.1)
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.05), add = 0)
  ) +
  theme(
    panel.background = element_rect(
      color = "black",
      fill = NA
    ),
    panel.grid = element_blank(),
    axis.text.y = element_text(
      color = "black",
      size = 12,
      face = "bold"
    ),
    axis.title.y = element_text(
      color = "black",
      size = 14,
      face = "bold"
    ),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  labs(
    y = "Density"
  ) ->
  p_hist_den;p_hist_den

ax <- tibble::tibble(
  xmin = c(0, 1, 2),
  xmax = c(1, 2, 3),
  ymin = c(2.5, 2, 2.5),
  ymax = c(3, 3.5, 3),
  x = c(0.5, 1.5, 2.5),
  y = 1,
  label = c("5'UTR", "CDS", "3'UTR")
)

ax |> 
  ggplot(aes(
    xmin = xmin,
    xmax = xmax,
    ymin = ymin,
    ymax = ymax,
    x = x,
    y = y,
    label = label
  )) +
  geom_rect() +
  geom_text(
    color = "black",
    size = 10,
    fontface = "bold"
  ) +
  geom_hline(
    yintercept = -0.5,
    color = "white"
  ) +
  scale_x_continuous(
    expand = expansion(mult = c(0.005), add = 0),
    limits = c(0, 3),
    breaks = seq(0, 3, 0.1)
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.1, 0), add = 0)
  ) +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank()
  ) ->
  ax_plot;ax_plot


cowplot::plot_grid(
  plotlist = list(p_hist_den, ax_plot),
  align = 'v',
  ncol = 1,
  rel_heights = c(1, 0.2)
)  ->
  p_all

ggsave(
  filename = "histo-plot-c-to-u-distribution.pdf",
  plot = p_all,
  device = "pdf",
  path = "/home/liuc9/github/m6a/test/FLARE/plots",
  width = 10,
  height = 5
)


anno |> 
  dplyr::filter(
    t > 1,
  ) |>
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
  p_pie_t1;p_pie_t1

ggsave(
  filename = "pie-plot-t1-c-to-u-distribution.pdf",
  plot = p_pie_t1,
  device = "pdf",
  path = "/home/liuc9/github/m6a/test/FLARE/plots",
  width = 9,
  height = 9
)

anno_region_dis_unnest |> 
  dplyr::filter(
    t > 1
  ) |>
  ggplot(aes(
    x = ratio_new
  )) +
  geom_histogram(
    aes(y = ..density..),
    fill = "blue",
    color = "white",
    alpha = 0.3,
    linewidth = 0.1,
    bins = 30,
    binwidth = 0.1,
    boundary=0
  ) +
  geom_density(
    alpha = 0.1,
    linetype = 2,
    linewidth = 1.2,
    colour = "red",
    fill = "red"
  ) +
  geom_vline(
    xintercept = c(1, 2),
    color = "red",
    linetype = 3,
    linewidth = 1.2,
  ) +
  scale_x_continuous(
    expand = expansion(mult = c(0.005), add = 0),
    limits = c(0, 3),
    breaks = seq(0, 3, 0.1)
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.05), add = 0)
  ) +
  theme(
    panel.background = element_rect(
      color = "black",
      fill = NA
    ),
    panel.grid = element_blank(),
    axis.text.y = element_text(
      color = "black",
      size = 12,
      face = "bold"
    ),
    axis.title.y = element_text(
      color = "black",
      size = 14,
      face = "bold"
    ),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  labs(
    y = "Density"
  ) ->
  p_hist_den_t1;p_hist_den_t1


cowplot::plot_grid(
  plotlist = list(p_hist_den_t1, ax_plot),
  align = 'v',
  ncol = 1,
  rel_heights = c(1, 0.2)
)  ->
  p_all_t1;p_all_t1

ggsave(
  filename = "histo-plot-t1-c-to-u-distribution.pdf",
  plot = p_all_t1,
  device = "pdf",
  path = "/home/liuc9/github/m6a/test/FLARE/plots",
  width = 10,
  height = 5
)



# logo --------------------------------------------------------------------

seqlevels(Hsapiens) |> head(30)


future::plan(future::multisession, workers = 10)
anno_region_dis_unnest |> 
  dplyr::mutate(
    seq = furrr::future_pmap(
      .l = list(
        .chrom = chrom,
        .end = end
      ),
      .f = function(.chrom, .end) {
        .chr <- glue::glue("chr{gsub('T', '', .chrom)}")
        .point <- .end
        .point_up <- .point - 4
        .point_end <- .point + 4
        .seq <- getSeq(
          Hsapiens,
          .chr,
          .point_up,
          .point_end
        )
      }
    )
  ) ->
  anno_region_dis_unnest_seq
future::plan(future::sequential)

future::plan(future::multisession, workers = 20)
anno_region_dis_unnest_seq |> 
  dplyr::mutate(
    s_seq = purrr::map2_chr(
      .x = strand,
      .y = seq,
      .f = function(.x, .y) {
        if(.x == "+") {
          toString(
            RNAString(.y)
          )
        } else {
          toString(Biostrings::reverseComplement(
            RNAString(.y)
          ))
        }
      }
    )
  ) ->
  anno_region_dis_unnest_seq_s
future::plan(future::sequential)


ggplot() +
  ggseqlogo::geom_logo(
    data = anno_region_dis_unnest_seq_s |> 
      dplyr::pull(s_seq_p), 
    method = "prob"
  ) +
  theme(
    panel.background = element_blank(),
    axis.line = element_line(
      color = "black"
    ),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text = element_text(
      color = "black",
      size = 12,
      face = "bold"
    ),
    axis.title = element_text(
      color = "black",
      size = 14,
      face = "bold"
    )
  ) ->
  p_seqlogo;p_seqlogo


ggsave(
  filename = "seqlogo-plot-motif.pdf",
  plot = p_seqlogo,
  device = "pdf",
  path = "/home/liuc9/github/m6a/test/FLARE/plots",
  width = 10,
  height = 4
)


ggplot() +
  ggseqlogo::geom_logo(
    data = anno_region_dis_unnest_seq_s |>
      dplyr::filter(t > 1) |> 
      dplyr::pull(s_seq), 
    method = "prob"
  ) +
  theme(
    panel.background = element_blank(),
    axis.line = element_line(
      color = "black"
    ),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text = element_text(
      color = "black",
      size = 12,
      face = "bold"
    ),
    axis.title = element_text(
      color = "black",
      size = 14,
      face = "bold"
    )
  ) ->
  p_seqlogo_t1


ggsave(
  filename = "seqlogo-plot-motif-t1.pdf",
  plot = p_seqlogo_t1,
  device = "pdf",
  path = "/home/liuc9/github/m6a/test/FLARE/plots",
  width = 10,
  height = 4
)



ggplot() +
  ggseqlogo::geom_logo(
    data = anno_region_dis_unnest_seq_s |>
      dplyr::mutate(
        ac = substr(x = s_seq, start = 4, stop = 5)
      ) |> 
      dplyr::filter(ac == "AC") |> 
      dplyr::filter(t > 1) |>
      dplyr::pull(s_seq), 
    method = "prob"
  ) +
  theme(
    panel.background = element_blank(),
    axis.line = element_line(
      color = "black"
    ),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text = element_text(
      color = "black",
      size = 12,
      face = "bold"
    ),
    axis.title = element_text(
      color = "black",
      size = 14,
      face = "bold"
    )
  ) ->
  p_seqlogo_ac;p_seqlogo_ac


ggsave(
  filename = "seqlogo-plot-motif-ac-t1.pdf",
  plot = p_seqlogo_ac,
  device = "pdf",
  path = "/home/liuc9/github/m6a/test/FLARE/plots",
  width = 10,
  height = 4
)

# AC plot -----------------------------------------------------------------


anno_region_dis_unnest_seq_s |>
  dplyr::mutate(
    ac = substr(x = s_seq, start = 4, stop = 5)
  ) |> 
  dplyr::filter(ac == "AC") ->
  anno_region_dis_unnest_seq_s_ac



# test --------------------------------------------------------------------



anno_region_dis_unnest |> 
  dplyr::arrange(-edit_ct_ratio, -t)

anno_region_dis_unnest |> 
  dplyr::group_by(
    gene_region
  ) |> 
  dplyr::summarise(m = median(feature_size)) ->
  anno_region_dis_unnest_m

anno_region_dis_unnest_m

utr5_sf <- anno_region_dis_unnest_m$m[[2]] / anno_region_dis_unnest_m$m[[1]]
utr3_sf <- anno_region_dis_unnest_m$m[[3]] / anno_region_dis_unnest_m$m[[1]]

anno_region_dis_unnest |>
  dplyr::filter(
    gene_region == "CDS",
    strand == "-"
  ) |>
  dplyr::select(
    start,
    transcript_id,
    rel_loc,
    feature_size,
    ratio_new
  )
  

anno_region_dis_unnest |> 
  dplyr::filter(gene_name == "ACTB")

anno_region_dis_unnest |> 
  dplyr::filter(
    t > 1,
  ) ->
  anno_region_dis_unnest_t1

anno_region_dis_unnest_t1 |> 
  ggplot(aes(
    x = ratio_new 
  )) +
  geom_histogram(
    data = transform(anno_region_dis_unnest_t1, chrom = NULL),
    fill = "grey85"
  ) +
  geom_histogram() +
  facet_wrap(
    facets = dplyr::vars(chrom),
    scales = "free"
  )

# 
# anno_region_dis_unnest_t1 |> 
#   dplyr::group_by(gene_id) |> 
#   dplyr::count() |> 
#   dplyr::arrange(-n)

anno_region_dis_unnest_t1 |> 
  dplyr::filter(
    gene_id == "ENSG00000198492"
  ) |> 
  dplyr::select(strand, rel_loc, feature_size, ratio, ratio_new)


anno_region_dis_unnest |> 
  dplyr::filter(gene_region == "three_prime_utr") |> 
  dplyr::select(strand, rel_loc, feature_size, ratio, ratio_new) |> 
  dplyr::mutate(
    dis = ifelse(
      strand == "-",
      -rel_loc,
      rel_loc
    )
  ) |> 
  ggplot(aes(dis)) +
  geom_histogram(bins = 100)


anno_region_dis_unnest |> 
  dplyr::filter(gene_region == "three_prime_utr") |> 
  dplyr::arrange(-feature_size) |> 
  dplyr::select(feature_size)



anno_region_dis_unnest |> 
  dplyr::filter(
    t > 1,
    confidence > 0.
  ) |> 
  ggplot(aes(
    x = ratio_new 
  )) +
  geom_density()


anno_region_dis_unnest |> 
  dplyr::mutate(
    ratio_new_scale = dplyr::case_match(
      gene_region,
      "CDS" ~ ratio_new,
      "five_prime_utr" ~ scales::rescale(ratio_new, to = c(1-utr5_sf, 1), from = c(0,1)),
      "three_prime_utr" ~ scales::rescale(ratio_new, to = c(2, 2+utr3_sf), from = c(2,3))
    )
  ) ->
  anno_region_dis_unnest_scale

anno_region_dis_unnest_scale |> 
  dplyr::filter(
    t > 1,
    confidence > 0.
  ) |> 
  ggplot(aes(
    x = ratio_new_scale 
  )) +
  geom_density()

# footer ------------------------------------------------------------------

future::plan(future::sequential)

# save image --------------------------------------------------------------

save.image(file = "/home/liuc9/data/projdata/m6a/flare/sailor/summary.rda")
