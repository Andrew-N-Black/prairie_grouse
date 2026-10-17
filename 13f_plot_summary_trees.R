#!/usr/bin/env Rscript
# =============================================================================
# 13f_plot_summary_trees.R — Figure 18: autosomal and Z summary trees beside the mitogenome tree
#
# Inputs : autosomes.gcf.cf.tree, chrZ.gcf.cf.tree   (13e_summary_trees.sh; gCF on nodes)
#          mito_tree.contree                         (12d_mito_tree.sh; UFBoot on nodes)
# Output : Figure18_summary_trees.pdf / .png
#
# Install once:
#   install.packages(c("BiocManager", "ape", "ggplot2", "patchwork"))
#   BiocManager::install(c("ggtree", "treeio"))
#
# Run:  Rscript 13f_plot_summary_trees.R autosomes.gcf.cf.tree chrZ.gcf.cf.tree mito_tree.contree
#
# Nuclear panels are drawn as cladograms (ASTRAL gives no terminal branch lengths); node labels
# are gene concordance factors, the % of 100-kb window trees containing that branch. The
# mitogenome panel is drawn as in Figure 16 (phylogram, UFBoot). Bars mark every species, and
# every species pair, that is monophyletic in that tree.
# =============================================================================
suppressPackageStartupMessages({
  library(ape)
  library(ggtree)
  library(ggplot2)
  library(patchwork)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) stop("usage: Rscript 13f_plot_summary_trees.R autosomes.gcf.cf.tree chrZ.gcf.cf.tree mito_tree.contree")
NUC_OUTGROUP  <- "Lmuta"
MITO_OUTGROUP <- "PV239887.1"
GCF_MIN  <- 20        # show gCF on nodes with at least this value (and >= 3 descendant tips)
UFB_MIN  <- 70        # mitogenome: show UFBoot >= this
COLS     <- c(LEPC = "#8B4513", GRPC = "#DAA520", STGR = "black")
FONT     <- "Times"   # "serif" if Times is not registered

species_of <- function(tips) {
  ifelse(grepl("^LEPC", tips), "LEPC",
  ifelse(grepl("^STGR", tips), "STGR", "GRPC"))      # MW574394.1 is GRPC
}
display_of <- function(tips) {
  d <- sub("^(LEPC|GRPC|STGR)_(F[0-9]+)_hap([12])$", "\\2 h\\3", tips)
  d <- sub("^(LEPC|GRPC|STGR)_", "", sub("_mitogenome$", "", d))
  d[tips == "MW574394.1"] <- "MW574394.1 (ref.)"
  d
}

prep_tree <- function(file, outgroup) {
  tr <- read.tree(file)
  # IQ-TREE writes "gCF" or "support/gCF" labels; keep the last number
  tr$node.label <- suppressWarnings(as.numeric(sub(".*/", "", tr$node.label)))
  tr <- root(tr, outgroup = outgroup, resolve.root = TRUE, edgelabel = TRUE)
  tr <- drop.tip(tr, outgroup)
  ladderize(tr, right = FALSE)
}

clade_bars <- function(tr) {
  sp <- species_of(tr$tip.label)
  groups <- list(LEPC = "LEPC", GRPC = "GRPC", STGR = "STGR",
                 "LEPC + GRPC" = c("LEPC", "GRPC"), "GRPC + STGR" = c("GRPC", "STGR"),
                 "LEPC + STGR" = c("LEPC", "STGR"))
  out <- list()
  for (g in names(groups)) {
    tt <- tr$tip.label[sp %in% groups[[g]]]
    if (length(tt) > 1 && length(tt) < Ntip(tr) && is.monophyletic(tr, tt))
      out[[g]] <- list(node = getMRCA(tr, tt), pair = length(groups[[g]]) == 2,
                       col = if (length(groups[[g]]) == 1) COLS[[g]] else "black")
  }
  out
}

panel <- function(tr, title, cladogram, min_lab, scale_w = NULL) {
  meta <- data.frame(label = tr$tip.label, species = factor(species_of(tr$tip.label), levels = names(COLS)),
                     display = display_of(tr$tip.label), stringsAsFactors = FALSE)
  p <- ggtree(tr, linewidth = 0.45, branch.length = if (cladogram) "none" else "branch.length") %<+% meta
  xmax <- max(p$data$x)
  p$data$support <- suppressWarnings(as.numeric(p$data$label))
  p$data$support[p$data$isTip] <- NA
  ndesc <- sapply(p$data$node, function(n) if (n <= Ntip(tr)) 1 else length(extract.clade(tr, n)$tip.label))
  p$data$show <- !p$data$isTip & !is.na(p$data$support) & p$data$support >= min_lab & ndesc >= 3
  p <- p +
    geom_tippoint(aes(fill = species), shape = 21, size = 1.9, colour = "black", stroke = 0.25) +
    geom_tiplab(aes(label = display), size = 2.3, family = FONT, offset = xmax * 0.02) +
    geom_text2(aes(subset = show, label = round(support)), hjust = 1.3, vjust = -0.4,
               size = 2.0, colour = "grey25", family = FONT)
  for (g in names(b <- clade_bars(tr))) {
    p <- p + geom_cladelab(node = b[[g]]$node, label = g, align = TRUE,
                           offset = xmax * (if (b[[g]]$pair) 0.62 else 0.42),
                           offset.text = xmax * 0.02, barsize = 0.7, fontsize = 2.8, family = FONT,
                           textcolour = b[[g]]$col, barcolour = b[[g]]$col)
  }
  if (!is.null(scale_w)) p <- p + geom_treescale(x = 0, y = -0.8, width = scale_w, fontsize = 2.4,
                                                 linesize = 0.5, family = FONT)
  p + scale_fill_manual(values = COLS, drop = FALSE, name = NULL) + hexpand(1.0) +
    ggtitle(title) + theme(plot.title = element_text(family = FONT, size = 10, hjust = 0.5))
}

options(scipen = 10)
pa <- panel(prep_tree(args[1], NUC_OUTGROUP), "Autosomes", TRUE, GCF_MIN)
pz <- panel(prep_tree(args[2], NUC_OUTGROUP), "Z chromosome", TRUE, GCF_MIN)
pm <- panel(prep_tree(args[3], MITO_OUTGROUP), "Mitogenome", FALSE, UFB_MIN, scale_w = 0.0005)

fig <- (pa | pz | pm) + plot_layout(guides = "collect") + plot_annotation(tag_levels = "A") &
  theme(legend.position = "bottom", legend.text = element_text(family = FONT, size = 10),
        plot.tag = element_text(family = FONT, face = "bold"))
fig <- fig & guides(fill = guide_legend(override.aes = list(size = 3.5)))

ggsave("Figure18_summary_trees.pdf", fig, width = 11, height = 8.5)
ggsave("Figure18_summary_trees.png", fig, width = 11, height = 8.5, dpi = 300)
message("wrote Figure18_summary_trees.pdf / .png")
