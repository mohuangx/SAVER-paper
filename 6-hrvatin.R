# Mo Huang, mohuang@wharton.upenn.edu
# Hrvatin analysis

# x <- readRDS("SAVER-data/hrvatin.rds")
# x.saver <- readRDS("SAVER-data/hrvatin_saver.rds")
# 
# ct.dat <- read.csv("SAVER-data/hrvatin_celltypes.csv", header = TRUE,
#                    row.names = 1, stringsAsFactors = FALSE)
# 
# ct.dat2 <- ct.dat[!is.na(ct.dat$celltype), ]
# ct.dat2$celltype[!is.na(ct.dat2$subtype)] <- 
#   ct.dat2$subtype[!is.na(ct.dat2$subtype)]
# table(ct.dat2$celltype)
# ct.dat3 <- ct.dat2[ct.dat2$celltype != "ExcL23", ]
# 
# ident <- ct.dat3[colnames(x), 4]
# 
# x.sub <- x[, !is.na(ident)]
# x.saver.sub <- x.saver[, !is.na(ident)]
# 
# ident2 <- ct.dat3[colnames(x.sub), 4]
# 
# non.exc <- grep("ExcL|Hip|RSP|Sub", ident2, invert = TRUE)
# 
# ident2[non.exc] <- "Other"
# 
# ident2 <- as.factor(ident2)
# 
# library(Seurat)
# library(ggplot2)
# library(cowplot)
# 
# 
# obs <- CreateSeuratObject(x.sub, project = "obs")
# saver <- CreateSeuratObject(x.saver.sub, project = "saver")
# 
# 
# obs <- NormalizeData(obs)
# saver <- NormalizeData(saver)
# 
# obs <- FindVariableGenes(obs, do.plot = FALSE)
# saver <- FindVariableGenes(saver, do.plot = FALSE)
# 
# obs <- ScaleData(obs)
# saver <- ScaleData(saver)
# 
# obs <- RunPCA(obs, pc.genes = obs@var.genes, do.print = FALSE, 
#               pcs.compute = 50)
# PCElbowPlot(obs, num.pc = 50)
# 
# saver <- RunPCA(saver, pc.genes = saver@var.genes, do.print = FALSE,
#                 pcs.compute = 50)
# PCElbowPlot(saver, num.pc = 50)
# 
# 
# obs <- JackStraw(object = obs, num.pc = 40, num.replicate = 100,
#                  do.print = TRUE)
# 
# p1 <- JackStrawPlot(obs, PCs = 1:40)
# levels(p1$data$PC.Score)
# 
# saver <- JackStraw(object = saver, num.pc = 40, num.replicate = 100,
#                    do.print = TRUE)
# 
# p2 <- JackStrawPlot(saver, PCs = 1:40)
# levels(p2$data$PC.Score)
# 
# dim.obs <- 35
# dim.saver <- 30
# 
# obs <- RunTSNE(obs, dims.use = 1:dim.obs, check_duplicates = FALSE,
#                do.fast = TRUE)
# 
# saver <- RunTSNE(saver, dims.use = 1:dim.saver, check_duplicates = FALSE,
#                  do.fast = TRUE)
# 
# ident3 <- factor(ident2, levels = levels(ident2)[c(3, 4, 5, 6, 7, 8, 9, 1, 2, 10, 12, 13, 11)])
# 
# 
# obs.tsne <- data.frame(obs@dr$tsne@cell.embeddings, ident2 = ident2,
#                        ident3 = ident3)
# 
# saver.tsne <- data.frame(saver@dr$tsne@cell.embeddings, ident2 = ident2,
#                          ident3 = ident3)
# 
# saveRDS(list(obs.tsne, saver.tsne), "SAVER-data/hrvatin_tsne.rds")

###############################################################################
## Plotting

tsne <- readRDS("SAVER-data/hrvatin_tsne.rds")

obs.tsne <- tsne[[1]]
saver.tsne <- tsne[[2]]

library(ggplot2)
library(cowplot)

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

my.cols <- c(gg_color_hue(12), "gray50")

p3 <- ggplot(obs.tsne) + geom_point(aes(x = tSNE_1, y = tSNE_2,
                                        colour = ident3), size = 0.5) +
  scale_colour_manual(name = "Cell type", values = my.cols) +
  guides(colour = FALSE) + 
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(size = 20, face = "plain")) +
  labs(title = "Observed")

p4 <- ggplot(saver.tsne) + geom_point(aes(x = tSNE_1, y = tSNE_2,
                                          colour = ident3), size = 0.5) +
  scale_colour_manual(name = "Cell type", values = my.cols) +
  guides(colour = FALSE) + 
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(size = 20, face = "plain")) +
  labs(title = "SAVER")


p5 <- ggplot(saver.tsne) + geom_point(aes(x = tSNE_1, y = tSNE_2,
                                          colour = ident3), size = 0.5) +
  scale_colour_manual(name = "Cell type", values = my.cols) +
  guides(colour = guide_legend(override.aes = list(size=5))) +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(size = 14)) +
  labs(title = "SAVER")

gglegend <- function(x){ 
  tmp <- ggplot_gtable(ggplot_build(x)) 
  leg <- which(sapply(tmp$grobs, function(y) y$name) == "guide-box") 
  tmp$grobs[[leg]]
}

png("plots/fig2e_hrvatin.png", 6, 3, units = "in", res = 300)
plot_grid(p3, p4)
dev.off()

library(grid)
library(gridExtra)

pdf("plots/fig2e_hrvatin_legend.pdf", 8, 4)
grid.draw(gglegend(p5))
dev.off()


