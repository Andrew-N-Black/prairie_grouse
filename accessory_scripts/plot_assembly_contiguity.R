# =============================================================================
# plot_assembly_contiguity.R
# USFWS report Objective 3, Figure S5: contiguity of the 46 haplotype-resolved,
# pseudo-chromosome-scale prairie grouse assemblies (QUAST, final RagTag-ordered
# assemblies after removal of unplaced scaffolds < 50 kb).
#
# Three panels: number of scaffolds, total assembly length (Mb) and scaffold
# N50 (Mb). One row per individual, grouped by species; hap1 and hap2 are
# joined by a line. W-bearing haplotypes of the six females are drawn as open
# symbols (they lack the ~75 Mb Z chromosome, hence shorter length and N50).
#
# Run in RStudio (needs ggplot2; patchwork optional but recommended):
#   install.packages(c("ggplot2", "patchwork"))
# Outputs (working directory): FigS_assembly_contiguity.pdf / .png
# =============================================================================

library(ggplot2)

# ---- data (Table S3) -------------------------------------------------------
asm <- read.table(header = TRUE, stringsAsFactors = FALSE, text = "
sample species scaf_h1 scaf_h2 len_h1 len_h2 n50_h1 n50_h2
F5457 STGR 232 332 1040.3 1051.8 70.30 71.43
F5462 STGR 223 251 1049.2 1058.7 71.04 71.54
F5463 STGR 192 290 1038.4 1052.2 70.69 71.28
F5468 STGR 211 249 1044.5 1054.5 71.67 71.01
F5472 STGR 169 231 1045.5 1058.8 71.30 71.15
F5478 STGR 249 254 1040.2 1048.1 69.05 71.75
F5480 STGR 212 270 1042.9 1051.1 70.76 70.84
F5485 STGR 238 206 1047.7 1050.8 70.98 70.64
F5502 STGR 274 191 1058.2 1039.9 71.08 71.06
F5503 STGR 207 332 975.3 1075.3 65.51 68.21
F5540 LEPC 205 299 1045.4 1064.8 71.35 71.13
F5541 LEPC 316 248 1063.4 1048.1 70.53 70.97
F5542 LEPC 279 265 1051.3 1055.7 71.02 71.28
F5543 LEPC 322 275 1059.8 1046.6 71.46 71.02
F5544 LEPC 234 272 1053.4 1060.3 71.28 70.75
F5545 LEPC 206 331 1042.4 1056.7 70.69 71.17
F5546 LEPC 234 343 1050.2 1069.6 71.34 71.88
F5595 GRPC 322 191 1094.8 980.2 71.19 65.49
F5596 GRPC 320 190 1090.7 981.6 71.35 65.62
F5597 GRPC 288 213 1063.7 1057.1 71.87 71.13
F5598 GRPC 226 314 973.3 1076.2 67.04 70.40
F5599 GRPC 258 368 972.9 1080.3 65.51 70.54
F5600 GRPC 247 389 978.9 1081.7 66.15 70.95
")

# W-bearing haplotypes of the six females (Table S3; CHD1 PCR, Figure S7)
w_hap <- c("F5503_hap1", "F5595_hap2", "F5596_hap2",
           "F5598_hap1", "F5599_hap1", "F5600_hap1")

# ---- reshape to long format (base R, no tidyr needed) -----------------------
long <- do.call(rbind, lapply(c("h1", "h2"), function(h) {
  data.frame(sample  = asm$sample,
             species = asm$species,
             hap     = ifelse(h == "h1", "hap1", "hap2"),
             Scaffolds = asm[[paste0("scaf_", h)]],
             Length    = asm[[paste0("len_",  h)]],
             N50       = asm[[paste0("n50_",  h)]])
}))
long$id      <- paste(long$sample, long$hap, sep = "_")
long$W       <- long$id %in% w_hap
long$species <- factor(long$species, levels = c("LEPC", "GRPC", "STGR"))
# order individuals: by species, then sample ID (top to bottom)
ord <- with(asm, sample[order(factor(species, levels = c("LEPC", "GRPC", "STGR")), sample)])
long$sample <- factor(long$sample, levels = rev(ord))

sp_cols <- c(LEPC = "#A52A2A", GRPC = "#DAA520", STGR = "#000000")

# ---- one panel: W-bearing haplotypes drawn with open symbols ----------------
long$hap_W <- ifelse(long$W, paste0(long$hap, "_W"), long$hap)
panel <- function(var, xlab, show_y = TRUE, show_strip = FALSE) {
  ggplot(long, aes(x = .data[[var]], y = sample)) +
    geom_line(aes(group = sample, colour = species), linewidth = 0.4, alpha = 0.6) +
    geom_point(aes(colour = species, shape = hap_W), size = 2.2, stroke = 0.9) +
    scale_shape_manual(values = c(hap1 = 16, hap2 = 17, hap1_W = 1, hap2_W = 2),
                       breaks = c("hap1", "hap2", "hap1_W"),
                       labels = c("hap1", "hap2", "W-bearing (open)"),
                       name = NULL) +
    scale_colour_manual(values = sp_cols, name = NULL) +
    scale_x_continuous(expand = expansion(mult = c(0.05, 0.08))) +
    facet_grid(species ~ ., scales = "free_y", space = "free_y") +
    labs(x = xlab, y = NULL) +
    theme_bw(base_size = 10) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.y = element_line(colour = "grey92"),
          axis.text.y = if (show_y) element_text(size = 8) else element_blank(),
          axis.ticks.y = if (show_y) element_line() else element_blank(),
          strip.text.y.right = if (show_strip) element_text(angle = 0, face = "bold") else element_blank(),
          strip.background = if (show_strip) element_rect(fill = "grey95") else element_blank(),
          legend.position = "bottom")
}

pA <- panel("Scaffolds", "Scaffolds (n)") + labs(tag = "A")
pB <- panel("Length",    "Total length (Mb)", show_y = FALSE) + labs(tag = "B")
pC <- panel("N50",       "Scaffold N50 (Mb)", show_y = FALSE, show_strip = TRUE) + labs(tag = "C")

# ---- combine and save ---------------------------------------------------------
if (requireNamespace("patchwork", quietly = TRUE)) {
  library(patchwork)
  fig <- (pA | pB | pC) +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom")
  print(fig)
  ggsave("FigS_assembly_contiguity.pdf", fig, width = 9, height = 6.5)
  ggsave("FigS_assembly_contiguity.png", fig, width = 9, height = 6.5, dpi = 600)
} else {
  message("patchwork not installed: saving the three panels separately")
  for (p in list(A = pA, B = pB, C = pC)) print(p)
  ggsave("FigS_assembly_contiguity_A_scaffolds.png", pA, width = 4, height = 6.5, dpi = 600)
  ggsave("FigS_assembly_contiguity_B_length.png",    pB, width = 3.5, height = 6.5, dpi = 600)
  ggsave("FigS_assembly_contiguity_C_N50.png",       pC, width = 3.5, height = 6.5, dpi = 600)
}

# quick summary for the text
cat(sprintf("N50: %.2f-%.2f Mb (median %.2f); length: %.1f-%.1f Mb (median %.1f); scaffolds: %d-%d (median %d)\n",
            min(long$N50), max(long$N50), median(long$N50),
            min(long$Length), max(long$Length), median(long$Length),
            min(long$Scaffolds), max(long$Scaffolds), as.integer(median(long$Scaffolds))))
