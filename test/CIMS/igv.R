library(ggbio)

p.ideo <- Ideogram(genome = "hg38")
p.ideo

library(GenomicRanges)
p.ideo + xlim(GRanges("chr2", IRanges(1e8, 1e8 + 10000000)))


library(ggbio)
library(Homo.sapiens)


data(genesymbol, package = "biovizBase")

wh <- genesymbol[c("BRCA1", "NBR1")]
wh <- range(wh, ignore.strand = TRUE)
wh
p.txdb <- autoplot(Homo.sapiens, which = wh)
p.txdb

autoplot(Homo.sapiens, which = wh, gap.geom = "chevron")

autoplot(Homo.sapiens,
  which = wh, label.color = "black", color = "brown",
  fill = "brown"
)


library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
autoplot(txdb, which = wh)

library(BSgenome.Hsapiens.UCSC.hg19)

bg <- BSgenome.Hsapiens.UCSC.hg19
p.bg <- autoplot(bg, which = wh)
p.bg + zoom(1 / 100)
p.bg + zoom(1 / 1000)
p.bg + zoom(1 / 2500)
library(BSgenome.Hsapiens.UCSC.hg19)
bg <- BSgenome.Hsapiens.UCSC.hg19
## force to use geom 'segment' at this level
autoplot(bg, which = resize(wh, width = width(wh) / 2000), geom = "segment")

fl.bam <- system.file("extdata", "wg-brca1.sorted.bam", package = "biovizBase")
wh <- keepSeqlevels(wh, "chr17")
autoplot(fl.bam, which = wh)

library(BSgenome.Hsapiens.UCSC.hg19)
bg <- BSgenome.Hsapiens.UCSC.hg19
p.mis <- autoplot(fl.bam, bsgenome = bg, which = wh, stat = "mismatch")
p.mis

library(VariantAnnotation)
fl.vcf <- system.file("extdata", "17-1409-CEU-brca1.vcf.bgz", package = "biovizBase")
vcf <- readVcf(fl.vcf, "hg19")
vr <- as(vcf[, 1:3], "VRanges")
vr <- renameSeqlevels(vr, value = c("17" = "chr17"))
## small region contains data
gr17 <- GRanges("chr17", IRanges(41234400, 41234530))
p.vr <- autoplot(vr, which = wh)
## none geom
p.vr
## rect geom
p.vr + xlim(gr17)
## text geom
p.vr + xlim(gr17) + zoom()
gr17 <- GRanges("chr17", IRanges(41234415, 41234569))
tks <- tracks(p.ideo,
  mismatch = p.mis, dbSNP = p.vr, ref = p.bg, gene = p.txdb,
  heights = c(2, 3, 3, 1, 4)
) + xlim(gr17) + theme_tracks_sunset()
tks


tks + zoom()
## zoom in with scale
p.txdb + zoom(1 / 8)
## zoom out
p.txdb + zoom(2)
## next view page
p.txdb + nextView()
## previous view page
p.txdb + prevView()


library(gmoviz)



bam <- Rsamtools::BamFile(file = "/mnt/isilon/xing_lab/liuc9/projdata/m6a/test4m_result/alignment/M6A/0/aligned.sort.bam")


.coverage <- gmoviz::getCoverage(
  regions_of_interest = "1",
  bam_file = "/mnt/isilon/xing_lab/liuc9/projdata/m6a/test4m_result/alignment/M6A/0/aligned.sort.bam",
  window_size = 1000
)



library(Gviz)

afrom <- 2960000
ato <- 3160000
alTrack <- AlignmentsTrack(
  system.file(package = "Gviz", "extdata", "gapped.bam"),
  isPaired = TRUE
)

bm <- biomaRt::useEnsembl(
  host = "https://grch37.ensembl.org",
  biomart = "ENSEMBL_MART_ENSEMBL",
  dataset = "hsapiens_gene_ensembl"
)
bmt <- BiomartGeneRegionTrack(
  genome = "hg19",
  chromosome = "chr12",
  start = afrom,
  end = ato,
  filter = list(with_refseq_mrna = TRUE),
  stacking = "dense",
  biomart = bm
)

plotTracks(c(bmt, alTrack),
  from = afrom, to = ato,
  chromosome = "chr12"
)

# Chromosome 1: 58,776,845-58,784,048
afrom <- 58776845 - 10000
ato <- 58784048 + 10000
alTrack <- AlignmentsTrack(
  "/mnt/isilon/xing_lab/liuc9/projdata/m6a/test4m_result/alignment/M6A/0/aligned.sort.bam"
)
bm <- biomaRt::useEnsembl(
  host = "https://www.ensembl.org",
  # version = 107,
  GRCh = "GRCh38",
  biomart = "ENSEMBL_MART_ENSEMBL",
  dataset = "hsapiens_gene_ensembl"
)
bmt <- BiomartGeneRegionTrack(
  genome = "hg38",
  chromosome = 1,
  start = afrom,
  end = ato,
  filter = list(with_refseq_mrna = TRUE),
  stacking = "dense",
  biomart = bm
)
plotTracks(c(bmt),
  from = afrom, to = ato,
  chromosome = 1
)
plotTracks(c(alTrack, bmt),
  from = afrom, to = ato,
  chromosome = 1,
  type = "coverage"
)


afrom <- 44945200
ato <- 44947200
alTrack <- AlignmentsTrack(
  system.file(package = "Gviz", "extdata", "snps.bam"), isPaired = TRUE)
plotTracks(alTrack, chromosome = "chr21", from = afrom, to = ato)


afrom <- 58776845 - 10000
ato <- 58784048 + 10000
alTrack <- AlignmentsTrack(
  "/mnt/isilon/xing_lab/liuc9/projdata/m6a/test4m_result/alignment/M6A/0/aligned.sort.bam"
)
plotTracks(alTrack, chromosome = "1", from = afrom, to = ato)
