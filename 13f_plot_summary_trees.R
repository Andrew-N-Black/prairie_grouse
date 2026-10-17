#!/usr/bin/env Rscript
# =============================================================================
# 13f_plot_summary_trees.R — Figure 17: autosomal and Z summary trees beside the mitogenome tree
#
# Inputs : autosomes.gcf.cf.tree, chrZ.gcf.cf.tree   (13e_summary_trees.sh; gCF on nodes)
#          mito_tree.contree                         (12d_mito_tree.sh; UFBoot on nodes)
# Output : Figure17_summary_trees.pdf / .png
#
# Install once:
#   install.packages(c("BiocManager", "ape", "ggplot2", "patchwork"))
#   BiocManager::install(c("ggtree", "treeio"))
#
# Run:  Rscript 13f_plot_summary_trees.R autosomes.gcf.cf.tree chrZ.gcf.cf.tree mito_tree.contree
#
# Nuclear panels are drawn as cladograms (ASTRAL gives no terminal branch lengths); node labels
# are gene concordance factors, the % of 100-kb window trees containing that branch; they are
# always shown on species and species-pair nodes, because ASTRAL resolves a clade that most
# windows contradict whenever it is the most frequent of the alternatives. The
# mitogenome panel is drawn as in 12e (phylogram, UFBoot). Bars mark every species, and
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
GCF_MIN  <- 50        # nuclear: gCF is always shown on species / species-pair nodes; elsewhere only if >= this
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

# Support values belong to splits (branches), not nodes. Rather than trusting edge-label
# handling during re-rooting, record each branch's value under its tip set (the side without
# the outgroup) in the tree as read, then look every node of the rooted tree up by its tips.
split_key <- function(tips) paste(sort(tips), collapse = ",")
prep_tree <- function(file, outgroup) {
  tr0 <- read.tree(file)
  val <- suppressWarnings(as.numeric(sub(".*/", "", tr0$node.label)))   # "gCF" or "x/gCF"
  lut <- c()
  for (i in seq_len(tr0$Nnode)) {
    tips <- extract.clade(tr0, Ntip(tr0) + i)$tip.label
    if (outgroup %in% tips) tips <- setdiff(tr0$tip.label, tips)
    if (length(tips) > 1 && !is.na(val[i])) lut[split_key(tips)] <- val[i]
  }
  tr <- root(tr0, outgroup = outgroup, resolve.root = TRUE)
  tr <- ladderize(drop.tip(tr, outgroup), right = FALSE)
  tr$node.label <- sapply(seq_len(tr$Nnode), function(i) {
    v <- lut[split_key(extract.clade(tr, Ntip(tr) + i)$tip.label)]
    if (length(v) && !is.na(v)) v else NA })
  tr
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
  p <- ggtree(tr, branch.length = if (cladogram) "none" else "branch.length") %<+% meta
  xmax <- max(p$data$x)
  p$data$support <- suppressWarnings(as.numeric(p$data$label))
  p$data$support[p$data$isTip] <- NA
  b <- clade_bars(tr)
  bar_nodes <- sapply(b, `[[`, "node")
  ndesc <- sapply(p$data$node, function(n) if (n <= Ntip(tr)) 1 else length(extract.clade(tr, n)$tip.label))
  p$data$show <- !p$data$isTip & !is.na(p$data$support) &
    (p$data$node %in% bar_nodes | (p$data$support >= min_lab & ndesc >= 3))
  # labels of neighbouring nodes (e.g. two short mitogenome branches) would overlap:
  # push the upper one up and the lower one down
  p$data$ly <- p$data$y
  sh <- which(p$data$show)
  sh <- sh[order(p$data$y[sh])]
  for (k in seq_along(sh)[-1]) {
    i <- sh[k - 1]; j <- sh[k]
    if (abs(p$data$y[j] - p$data$y[i]) < 0.9 && abs(p$data$x[j] - p$data$x[i]) < xmax * 0.15) {
      p$data$ly[i] <- p$data$y[i] - 0.35
      p$data$ly[j] <- p$data$y[j] + 0.35
    }
  }
  p$data$lab <- if (cladogram) paste0(round(p$data$support), "%") else as.character(round(p$data$support))
  # Layout. The nuclear panels have a FIXED x-range (coord_cartesian below), so one tree unit is
  # a known physical width and the gaps no longer shrink when anything else changes:
  # range = -0.25 .. 3.33 x tree width over a ~4-inch panel, i.e. tree width ~ 1.1 in.
  # Measured from the tips:  name 0.055 | names end ~0.42 | species bar 0.50, its label to ~0.85 |
  # pair bar 0.98, its label to ~1.7. The mitogenome panel keeps the 12e layout.
  tip_off  <- xmax * (if (cladogram) 0.055 else 0.03)
  off_sp   <- xmax * (if (cladogram) 0.50 else 0.30)
  off_pair <- xmax * (if (cladogram) 0.98 else 0.62)
  p <- p +
    geom_tippoint(aes(fill = species), shape = 21, size = 1.9, colour = "black", stroke = 0.25) +
    geom_tiplab(aes(label = display), size = 2.3, family = FONT,
                offset = tip_off) +
    geom_label2(aes(subset = show, label = lab, y = ly), size = 2.1, family = FONT, colour = "grey15",
                label.size = 0, label.padding = unit(0.08, "lines"), fill = "white", hjust = 1.1, vjust = -0.15)
  for (g in names(b)) {
    p <- p + geom_cladelab(node = b[[g]]$node, label = g, align = TRUE,
                           offset = if (b[[g]]$pair) off_pair else off_sp,
                           offset.text = xmax * 0.03, barsize = 0.7, fontsize = 2.8, family = FONT,
                           textcolour = b[[g]]$col, barcolour = b[[g]]$col)
  }
  if (!is.null(scale_w)) p <- p + geom_treescale(x = 0, y = -0.8, width = scale_w, fontsize = 2.4,
                                                 linesize = 0.5, family = FONT)
  p + scale_fill_manual(values = COLS, drop = FALSE, name = NULL) + (if (cladogram) coord_cartesian(xlim = c(-0.25, 3.33) * xmax, clip = "off") else hexpand(1.0)) +
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

ggsave("Figure17_summary_trees.pdf", fig, width = 13, height = 8.5)
ggsave("Figure17_summary_trees.png", fig, width = 13, height = 8.5, dpi = 300)
message("wrote Figure17_summary_trees.pdf / .png")
