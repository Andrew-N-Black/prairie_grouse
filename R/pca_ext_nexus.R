#Load libraries
library(readxl)
library(ggplot2)

#Read in metadata
metadata <- read_xlsx("/Users/andrewblack/Documents/Research/GROUSE/sarek_nexus_new_plus_shotguns/heterozygosity_extended_nexus.xlsx")
#Read in covariation matrix
cov<-as.matrix(read.table("~/final.cov"))

#Extract and calculate eplained variation
axes<-eigen(cov)
head(axes$values/sum(axes$values)*100)
#[1] 22.6180474  2.9970624  2.1391689  0.6050842  0.5466266
#[6]  0.4565761


#Bind vectors with metadata and plot
PC1_3<-as.data.frame(axes$vectors[,1:3])
x<-cbind(PC1_3,metadata)
 #By species and group
ggplot(data=x, aes(y=V2, x=V1)) +
    geom_point(size=6, color="black", aes(shape=GROUP, fill=SPECIES)) +
    theme_classic() +
    xlab("PC1 (22.6%)") + ylab("PC2 (3.0%)") +
    geom_hline(yintercept=0, linetype="dashed") +
    geom_vline(xintercept=0, linetype="dashed") +
    scale_fill_manual("Species", values=c("goldenrod","brown","black","grey")) +
    scale_shape_manual("Group", values=c(25,21)) +
    theme(legend.position = "right") +
    guides(
        fill = guide_legend(override.aes = list(shape=21, size=5, stroke=0.5)),
        shape = guide_legend(override.aes = list(fill="grey50", size=5))
    )+guides(
        fill = guide_legend(override.aes = list(shape=21, size=5, stroke=0.5)),
        shape = guide_legend(override.aes = list(fill=NA, size=5))
    )


-----------------------------------------------------------
#Or 3-D plot
#### Load libraries ####
library(readxl)
library(ggplot2)
library(plotly)
library(reticulate)

#### Read in metadata ####
metadata <- read_xlsx("/Users/andrewblack/Documents/Research/GROUSE/sarek_nexus_new_plus_shotguns/heterozygosity_extended_nexus.xlsx")
#### Read in covariance matrix ####
cov <- as.matrix(read.table("~/final.cov"))

#### Extract and calculate explained variation ####
axes <- eigen(cov)
head(axes$values / sum(axes$values) * 100)
# [1] 22.6180474  2.9970624  2.1391689  0.6050842  0.5466266
# [6]  0.4565761

#### Bind PC vectors with metadata ####
PC1_3 <- as.data.frame(axes$vectors[, 1:3])
x <- cbind(PC1_3, metadata)

#### 2D ggplot: PC1 vs PC2, colored by SPECIES, shaped by GROUP ####
ggplot(data = x, aes(y = V2, x = V1)) +
  geom_point(
    size = 6, color = "black",
    aes(shape = GROUP, fill = SPECIES),
    position = position_jitter(width = 0, height = 0.02)
  ) +
  theme_classic() +
  xlab("PC1 (22.6%)") + ylab("PC2 (3.0%)") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_fill_manual("SPECIES", values = c("goldenrod", "brown", "black", "grey")) +
  scale_shape_manual("GROUP", values = c(25, 21)) +
  theme(legend.position = "top") +
  guides(
    fill = guide_legend(override.aes = list(shape = 21, size = 5, stroke = 0.5)),
    shape = guide_legend(override.aes = list(fill = NA, size = 5))
  )

#### 3D plotly: PC1 vs PC2 vs PC3, colored by SPECIES, symbol by GROUP ####
species_colors <- c(
  "Tympanuchus cupido" = "goldenrod",
  "Tympanuchus pallidicinctus" = "brown",
  "Tympanuchus phasianellus" = "black",
  "Tympanuchus phasianellus/Tympanuchus cupido" = "grey"
)

group_symbols <- c("Allopatric" = "circle", "Sympatric" = "diamond")

p <- plot_ly()

# SPECIES traces
for (sp in names(species_colors)) {
  df_sub <- x[x$SPECIES == sp, ]
  p <- p %>% add_trace(
    data = df_sub,
    x = ~V1, y = ~V2, z = ~V3,
    type = 'scatter3d', mode = 'markers',
    marker = list(
      size = 8,
      color = species_colors[[sp]],
      symbol = group_symbols[df_sub$GROUP],
      line = list(color = 'black', width = 1)
    ),
    name = sp,
    legendgroup = "species",
    legendgrouptitle = list(text = "Species", font = list(color = 'black')),
    showlegend = TRUE
  )
}

# Dummy GROUP/symbol traces (legend-only, white fill, black font)
for (grp in names(group_symbols)) {
  p <- p %>% add_trace(
    x = 0, y = 0, z = 0,
    type = 'scatter3d', mode = 'markers',
    marker = list(
      size = 8, color = 'white', symbol = group_symbols[[grp]],
      line = list(color = 'black', width = 1)
    ),
    name = grp,
    legendgroup = "group",
    legendgrouptitle = list(text = "Group", font = list(color = 'black')),
    visible = "legendonly",
    showlegend = TRUE
  )
}

p <- p %>% layout(
  scene = list(
    xaxis = list(title = 'PC1 (22.6%)'),
    yaxis = list(title = 'PC2 (2.99%)'),
    zaxis = list(title = 'PC3 (2.13%)'),
    camera = list(eye = list(x = 1.5, y = 1.5, z = 1.5)),
    aspectmode = "cube"
  ),
  legend = list(
    groupclick = "toggleitem",
    font = list(color = 'black')
  ),
  margin = list(l = 50, r = 50, b = 50, t = 50)
)

p

#### Export 3D plot as SVG ####
save_image(p, path.expand("~/PCA_3D_plot.svg"), width = 1200, height = 900)
