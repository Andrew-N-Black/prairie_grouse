#!/usr/bin/env Rscript
# =============================================================================
# 12e_plot_mito_tree_ggtree.R — consensus mitogenome tree (IQ-TREE 3, UFBoot) with ggtree
#
# Input : mito_tree_rooted.contree  (or mito_tree.contree; the script roots it)
# Output: mito_tree_ggtree.pdf / .png
#
# Install once (Bioconductor):
#   install.packages(c("BiocManager", "ape", "ggplot2"))
#   BiocManager::install(c("ggtree", "treeio"))
#
# Run:  Rscript 12e_plot_mito_tree_ggtree.R mito_tree_rooted.contree
# =============================================================================
suppressPackageStartupMessages({
  library(ape)
  library(ggtree)
  library(ggplot2)
})

args      <- commandArgs(trailingOnly = TRUE)
TREE_FILE <- if (length(args)) args[1] else "mito_tree_rooted.contree"
OUTGROUP  <- "PV239887.1"     # Lagopus muta
KEEP_OUT  <- FALSE            # TRUE draws the ptarmigan too (its long branch squashes the ingroup)
HIC_BIRDS <- c("F5485", "F5502")   # Hi-C read consensus mitogenomes, marked with *
MIN_BS    <- 70                    # show UFBoot support >= this
COLS      <- c(LEPC = "#8B4513", GRPC = "#DAA520", STGR = "black", Outgroup = "grey40")
FONT      <- "Times"               # use "serif" if Times is not registered on your system

# ---- tree: root on the outgroup, optionally drop it, ladderize ----------------
tr <- read.tree(TREE_FILE)
tr <- root(tr, outgroup = OUTGROUP, resolve.root = TRUE)
if (!KEEP_OUT) tr <- drop.tip(tr, OUTGROUP)
tr <- ladderize(tr, right = FALSE)

# ---- tip metadata -------------------------------------------------------------
options(scipen = 10)               # scale bar reads 0.0005, not 5e-04
tips <- tr$tip.label
meta <- data.frame(label = tips, stringsAsFactors = FALSE)
meta$species <- ifelse(grepl("^LEPC", tips), "LEPC",
                ifelse(grepl("^STGR", tips), "STGR",
                ifelse(tips == OUTGROUP, "Outgroup", "GRPC")))   # MW574394.1 is GRPC
meta$species <- factor(meta$species, levels = c("LEPC", "GRPC", "STGR", "Outgroup"))
# tip label = sample ID only; species is shown by the circle
meta$display <- sub("^(LEPC|GRPC|STGR)_", "", sub("_mitogenome$", "", tips))
meta$display[tips == "MW574394.1"] <- "MW574394.1 (reference)"
meta$display[tips == OUTGROUP]     <- "Lagopus muta (PV239887.1)"
is_hic <- grepl(paste(HIC_BIRDS, collapse = "|"), tips)
meta$display[is_hic] <- paste0(meta$display[is_hic], " *")

# ---- plot ---------------------------------------------------------------------
p <- ggtree(tr, linewidth = 0.6) %<+% meta
xmax <- max(p$data$x)

# bootstrap values live in the internal-node labels
p$data$support <- suppressWarnings(as.numeric(p$data$label))
p$data$support[p$data$isTip] <- NA

# clade bars: LEPC, and GRPC + STGR
n_lepc <- getMRCA(tr, tips[meta$species == "LEPC"])
n_gs   <- getMRCA(tr, tips[meta$species %in% c("GRPC", "STGR")])

# CLADE_OFFSET is measured from the longest tip (align = TRUE), so both bars
# line up in one column to the right of all labels. Increase it if a bar still
# touches a label; hexpand() must grow with it so the bar label is not cut off.
CLADE_OFFSET <- xmax * 0.55

p <- p +
  geom_tippoint(aes(fill = species), shape = 21, size = 2.6, colour = "black", stroke = 0.3) +
  geom_tiplab(aes(label = display), size = 3.2, family = FONT, colour = "black",
              offset = xmax * 0.025) +
  geom_text2(aes(subset = !isTip & !is.na(support) & support >= MIN_BS, label = support),
             hjust = 1.35, vjust = -0.45, size = 2.4, colour = "grey25", family = FONT) +
  geom_cladelab(node = n_lepc, label = "LEPC", align = TRUE, offset = CLADE_OFFSET,
                offset.text = xmax * 0.02, barsize = 0.8, fontsize = 3.4, family = FONT,
                textcolour = COLS[["LEPC"]], barcolour = COLS[["LEPC"]]) +
  geom_cladelab(node = n_gs, label = "GRPC + STGR", align = TRUE, offset = CLADE_OFFSET,
                offset.text = xmax * 0.02, barsize = 0.8, fontsize = 3.4, family = FONT) +
  geom_treescale(x = 0, y = -0.5, width = 0.0005, fontsize = 3, linesize = 0.6,
                 family = FONT) +
  scale_fill_manual(values = COLS, drop = TRUE, name = NULL) +
  hexpand(0.9) +
  theme(legend.position = "bottom",
        legend.text = element_text(family = FONT, size = 10),
        plot.margin = margin(5, 5, 5, 5),
        plot.caption = element_text(family = FONT, size = 8, hjust = 0)) +
  guides(fill = guide_legend(override.aes = list(size = 4))) +
  labs(caption = "* Hi-C read consensus; node labels = UFBoot support (>= 70); scale = substitutions/site")

ggsave("mito_tree_ggtree.pdf", p, width = 6.5, height = 6.5)
ggsave("mito_tree_ggtree.png", p, width = 6.5, height = 6.5, dpi = 300)
message("wrote mito_tree_ggtree.pdf / .png")
