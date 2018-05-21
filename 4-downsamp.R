# Mo Huang, mohuang@wharton.upenn.edu
# Perform down-sampling analysis


ref.files <- list.files("SAVER-data", pattern = "ref", full.names = TRUE)

samp.files <- list.files("SAVER-data", pattern = "samp.rds", full.names = TRUE)

saver.files <- list.files("SAVER-data", pattern = "samp_saver.rds", 
                          full.names = TRUE)

magic.files <- list.files("SAVER-data", pattern = "samp_magic.rds", 
                          full.names = TRUE)

scimpute.files <- list.files("SAVER-data", pattern = "samp_scimpute.rds",
                             full.names = TRUE)

dat <- vector("list", 4)
names(dat) <- c("baron", "chen", "manno", "zeisel")

normalizeData <- function(x, y = x) {
  sf <- colSums(y)/mean(colSums(y))
  return(sweep(x, 2, sf, "/"))
}

for (i in 1:4) {
  dat[[i]] <- vector("list", 7)
  names(dat[[i]]) <- c("X", "Y", "saver", "magic", "scimpute", "X.norm", "Y.norm")
  dat[[i]][[1]] <- readRDS(ref.files[i])
  dat[[i]][[2]] <- readRDS(samp.files[i])
  dat[[i]][[3]] <- readRDS(saver.files[i])
  dat[[i]][[4]] <- readRDS(magic.files[i])
  dat[[i]][[5]] <- readRDS(scimpute.files[i])
  dat[[i]][[5]][is.na(dat[[i]][[5]])] <- 0
}

n.cells <- sapply(dat, function(x) ncol(x[[1]]))
n.genes <- sapply(dat, function(x) nrow(x[[1]]))

cell.names <- sapply(dat, function(x) colnames(x[[1]]))
gene.names <- sapply(dat, function(x) rownames(x[[1]]))


for (i in 1:4) {
  dat[[i]][[6]] <- normalizeData(dat[[i]][[1]])
  dat[[i]][[7]] <- normalizeData(dat[[i]][[2]])
  dat[[i]][[5]] <- normalizeData(dat[[i]][[5]])
}

###############################################################################
## Correlation with reference plots

get.cor.gene <- function(X, Y) {
  sapply(1:nrow(X), function(i) cor(X[i, ], Y[i, ]))
}

get.cor.cell <- function(X, Y) {
  sapply(1:ncol(X), function(i) cor(X[, i], Y[, i]))
}

cor.dat <- vector("list", 2)
names(cor.dat) <- c("gene", "cell")

for (i in 1:2) {
  cor.dat[[i]] <- vector("list", 4)
  names(cor.dat[[i]]) <- names(dat)
  for (j in 1:4) {
    cor.dat[[i]][[j]] <- vector("list", 4)
    names(cor.dat[[i]][[j]]) <- c("Obs", "SAVER", "MAGIC", "scImpute")
  }
}

for (i in 1:4) {
  for (j in 1:4) {
    ind <- c(7, 3, 4, 5)
    if (j != 2) {
      cor.dat[[1]][[i]][[j]] <- get.cor.gene(dat[[i]][[6]], dat[[i]][[ind[j]]])
      cor.dat[[2]][[i]][[j]] <- get.cor.cell(dat[[i]][[6]], dat[[i]][[ind[j]]])
    } else {
      cor.dat[[1]][[i]][[j]] <- get.cor.gene(dat[[i]][[6]], dat[[i]][[ind[j]]]$estimate)
      cor.dat[[2]][[i]][[j]] <- get.cor.cell(dat[[i]][[6]], dat[[i]][[ind[j]]]$estimate)
    }
    
  }
}

cor.chg.dat <- vector("list", 2)
names(cor.chg.dat) <- c("gene", "cell")

for (i in 1:2) {
  cor.chg.dat[[i]] <- vector("list", 4)
  names(cor.chg.dat[[i]]) <- names(dat)
  for (j in 1:4) {
    cor.chg.dat[[i]][[j]] <- vector("list", 3)
    names(cor.chg.dat[[i]][[j]]) <- c("SAVER", "MAGIC", "scImpute")
  }
}

calc.cor.chg <- function(x, y) {
  (x-y)/abs(y)*100
}

for (i in 1:4) {
  for (j in 1:3) {
    ind <- c(2, 3, 4)
    cor.chg.dat[[1]][[i]][[j]] <- calc.cor.chg(cor.dat[[1]][[i]][[ind[j]]],
                                               cor.dat[[1]][[i]][[1]])
    cor.chg.dat[[2]][[i]][[j]] <- calc.cor.chg(cor.dat[[2]][[i]][[ind[j]]],
                                               cor.dat[[2]][[i]][[1]])
  }
}


source("bin/violin_plot.R")

pdf("plots/fig2a_corplot.pdf", 8, 8)
par(mfrow = c(2, 2), cex.main = 1.5, mar = c(5, 3, 0, 0) + 0.1, oma = c(2, 3, 3, 1), 
    mgp = c(3.5, 1, 0),
    cex.lab = 1.3, font.lab = 2, cex.axis = 1.3, bty = "n", las = 1)

x <- c(1, 2, 3, 4)
plot(x, c(-10, -10, -10, -10), type = "p", ylab = " ", xlab = " ", cex = 1.5, 
     ylim = c(-0.2, 1), xlim = c(1, 9), lwd = 2, pch = 5, axes = FALSE, main = " ")
# axis(1, at = c(1.85, 3.35, 4.85, 6.35, 7.85, 9.35), labels = c("Obs", "SAVER", "Reg",
#                                                                "KNN", "SVD", "RF"))
axis(1, at = c(2.1, 4.1, 6.1, 8.1), labels = FALSE)
text(c(2.1, 4.1, 6.1, 8.1), par()$usr[3]-0.07,
     labels = c("Baron", "Chen", "La Manno", "Zeisel"), srt = 45, adj = 1, 
     xpd = TRUE, cex = 1.4)
axis(2, pos = 1.1)
par(las = 0)
mtext("Correlation with Reference", side = 2, line = 3, cex = 1.4)

fill <- c("white", "#cb181d", "#6baed6", "#bae4b3")

xloc <- seq(1.5, 2.7, 0.4)

for (i in 1:4) {
  for (j in 1:4) {
    boxplot.ej(cor.dat[[1]][[i]][[j]], xloc = xloc[j]+2*(i-1), 
               cex.boxpoint = ps, fill = fill[j])
  }
}


x <- c(1, 2, 3, 4, 5)
plot(x, c(-10, -10, -10, -10, -10), type = "p", ylab = " ", xlab = " ", cex = 1.5, 
     ylim = c(-0.2, 1), xlim = c(1, 9), lwd = 2, pch = 5, axes = FALSE, main = " ")
# axis(1, at = c(1.85, 3.35, 4.85, 6.35, 7.85, 9.35), labels = c("Obs", "SAVER", "Reg",
#                                                                "KNN", "SVD", "RF"))
axis(1, at = c(2.1, 4.1, 6.1, 8.1), labels = FALSE)
text(c(2.1, 4.1, 6.1, 8.1), par()$usr[3]-0.07,
     labels = c("Baron", "Chen", "La Manno", "Zeisel"), srt = 45, adj = 1, 
     xpd = TRUE, cex = 1.4)
par(las = 1)
axis(2, pos = 1.1)
par(las = 0)


for (i in 1:4) {
  for (j in 1:4) {
    boxplot.ej(cor.dat[[2]][[i]][[j]], xloc = xloc[j]+1.85*(i-1), 
               cex.boxpoint = ps, fill = fill[j])
  }
}

x <- c(1, 2, 3, 4, 5)
plot(0, type = "n", ylab = " ", xlab = " ", cex = 1.5, 
     ylim = c(-150, 150), xlim = c(1, 7), lwd = 2, pch = 5, axes = FALSE, main = " ")
# axis(1, at = c(1.85, 3.35, 4.85, 6.35, 7.85, 9.35), labels = c("Obs", "SAVER", "Reg",
#                                                                "KNN", "SVD", "RF"))
axis(1, at = c(1.85, 3.35, 4.85, 6.35), labels = FALSE)
text(c(1.85, 3.35, 4.85, 6.35), par()$usr[3]-14,
     labels = c("Baron", "Chen", "La Manno", "Zeisel"), srt = 45, adj = 1, 
     xpd = TRUE, cex = 1.4)
par(las = 1)
axis(2, pos = 1.1)
par(las = 0)
mtext("% Improvement", side = 2, line = 3, cex = 1.4)

fill <- c("#cb181d", "#6baed6", "#bae4b3")

segments(1, 0, 8.5, 0, col = "gray", lty = 2, lwd = 2)

xloc <- seq(1.5, 2.2, 0.35)

for (i in 1:4) {
  for (j in 1:3) {
    boxplot.ej(cor.chg.dat[[1]][[i]][[j]], xloc = xloc[j]+1.5*(i-1), 
               cex.boxpoint = ps, fill = fill[j])
  }
}




x <- c(1, 2, 3, 4, 5)
plot(0, type = "n", ylab = " ", xlab = " ", cex = 1.5, 
     ylim = c(-100, 100), xlim = c(1, 7), lwd = 2, pch = 5, axes = FALSE, main = " ")
# axis(1, at = c(1.85, 3.35, 4.85, 6.35, 7.85, 9.35), labels = c("Obs", "SAVER", "Reg",
#                                                                "KNN", "SVD", "RF"))
axis(1, at = c(1.85, 3.35, 4.85, 6.35), labels = FALSE)
text(c(1.85, 3.35, 4.85, 6.35), par()$usr[3]-10,
     labels = c("Baron", "Chen", "La Manno", "Zeisel"), srt = 45, adj = 1, 
     xpd = TRUE, cex = 1.4)
par(las = 1)
axis(2, pos = 1.1)
par(las = 0)

xloc <- seq(1.5, 2.2, 0.35)

segments(1, 0, 8.5, 0, col = "gray", lty = 2, lwd = 2)


for (i in 1:4) {
  for (j in 1:3) {
    boxplot.ej(cor.chg.dat[[2]][[i]][[j]], xloc = xloc[j]+1.5*(i-1), 
               cex.boxpoint = ps, fill = fill[j])
  }
}


mtext("Gene", outer = TRUE, cex = 1.6, at = 0.28)
mtext("Cell", outer = TRUE, cex = 1.6, at = 0.78)

plot(0, type = "n", ylab = " ", xlab = " ", cex = 1.5, 
     ylim = c(-0.2, 1), xlim = c(1, 7), lwd = 2, pch = 5, axes = FALSE, main = " ")
rect(5.25, 0.90, 5.5, 0.95, col = "white", lwd = 2)
text(6.5, 0.93, "Observed", cex = 1.1, font = 1)
rect(5.25, 0.82, 5.5, 0.87, col = "#cb181d", lwd = 2)
text(6.5, 0.85, "SAVER", cex = 1.1, font = 1)
rect(5.25, 0.74, 5.5, 0.79, col = "#6baed6", lwd = 2)
text(6.5, 0.77, "MAGIC", cex = 1.1, font = 1)
rect(5.25, 0.66, 5.5, 0.71, col = "#bae4b3", lwd = 2)
text(6.5, 0.69, "scImpute", cex = 1.1, font = 1)

dev.off()

###############################################################################
## Correlation matrix distance (takes a while to run)

# library(SAVER)
# 
# cor.gene.dat <- vector("list", 4)
# 
# calc_cmd <- function(R1, R2) {
#   traceR1R2 <- sum(diag(crossprod(R1, R2)))
#   R1.norm <- norm(R1, type = "F")
#   R2.norm <- norm(R2, type = "F")
#   return(1-traceR1R2/(R1.norm*R2.norm))
# }
# 
# cmd.gene <- matrix(0, 4, 4, 
#                    dimnames = list(c("Baron", "Chen", "Manno", "Zeisel"),
#                                    c("Observed", "SAVER", "MAGIC", "scImpute")))
# 
# for (i in 1:4) {
#   ref.cor <- cor(t(dat[[i]]$X.norm))
#   samp.cor <- cor(t(dat[[i]]$Y.norm))
#   cmd.gene[i, 1] <- calc_cmd(samp.cor, ref.cor)
#   print(1)
#   
#   samp.cor <- cor.genes(dat[[i]]$saver)
#   cmd.gene[i, 2] <- calc_cmd(samp.cor, ref.cor)
#   print(2)
#   
#   samp.cor <- cor(t(dat[[i]]$magic))
#   cmd.gene[i, 3] <- calc_cmd(samp.cor, ref.cor)
#   print(3)
#   
#   samp.cor <- cor(t(dat[[i]]$scimpute))
#   cmd.gene[i, 4] <- calc_cmd(samp.cor, ref.cor)
#   print(4)
# }
# 
# 
# cmd.cell <- matrix(0, 4, 4, 
#                    dimnames = list(c("Baron", "Chen", "Manno", "Zeisel"),
#                                    c("Observed", "SAVER", "MAGIC", "scImpute")))
# 
# for (i in 1:4) {
#   ref.cor <- cor(dat[[i]]$X.norm)
#   samp.cor <- cor(dat[[i]]$Y.norm)
#   cmd.cell[i, 1] <- calc_cmd(samp.cor, ref.cor)
#   print(1)
#   
#   samp.cor <- cor.cells(dat[[i]]$saver)
#   cmd.cell[i, 2] <- calc_cmd(samp.cor, ref.cor)
#   print(2)
#   
#   samp.cor <- cor(dat[[i]]$magic)
#   cmd.cell[i, 3] <- calc_cmd(samp.cor, ref.cor)
#   print(3)
#   
#   samp.cor <- cor(dat[[i]]$scimpute)
#   cmd.cell[i, 4] <- calc_cmd(samp.cor, ref.cor)
#   print(4)
# }
# 
# cmd <- list(gene = cmd.gene, cell = cmd.cell)
# 
# saveRDS(cmd, "SAVER-data/fig2b_cmd.rds")

cmd <- readRDS("SAVER-data/fig2b_cmd.rds")

pdf("plots/fig2b_cmd.pdf", 7.75, 3.25)

par(mfrow = c(1, 2), cex.main = 1.5, mar = c(4, 2, 0, 1.5) + 0.1, oma = c(1.5, 3, 2, 0), 
    mgp = c(3.5, 1, 0),
    cex.lab = 1.5, font.lab = 2, cex.axis = 1.3, bty = "n", las = 1)
x <- c(1, 2, 3, 4)
plot(x, c(-10, -10, -10, -10), type = "p", ylab = " ", xlab = " ", cex = 1.5, 
     ylim = c(0, 1), xlim = c(1, 7), lwd = 2, pch = 5, axes = FALSE, main = " ")
# axis(1, at = c(1.85, 3.35, 4.85, 6.35, 7.85, 9.35), labels = c("Obs", "SAVER", "Reg",
#                                                                "KNN", "SVD", "RF"))
axis(1, at = c(2, 3.5, 5, 6.5), labels = FALSE, pos = 0, lwd.ticks = 0)
axis(1, at = c(1.1, 7), labels = FALSE, lwd.ticks = 0, pos = 0)
text(c(2, 3.5, 5, 6.5), par()$usr[3]-0.02,
     labels = c("Baron", "Chen", "La Manno", "Zeisel"), srt = 45, adj = 1, 
     xpd = TRUE, cex = 1.4)
par(las = 1)
axis(2, pos = 1.1)
par(las = 0)
mtext("CMD", side = 2, line = 3, cex = 1.4)
mtext("Gene-to-Gene", cex = 1.6, line = 0.2)

xloc <- seq(1.5, 2.5, 0.25)

fill <- c("white", "#cb181d", "#6baed6", "#bae4b3")

for (i in 1:4) {
  for (j in 1:4) {
    rect(xloc[j]+1.5*(i-1), 0, xloc[j+1]+1.5*(i-1), cmd[[1]][i, j],
         col = fill[j], lwd = 2)
  }
}

x <- c(1, 2, 3, 4)
plot(x, c(-10, -10, -10, -10), type = "p", ylab = " ", xlab = " ", cex = 1.5, 
     ylim = c(0, 0.1), xlim = c(1, 7), lwd = 2, pch = 5, axes = FALSE, main = " ")
# axis(1, at = c(1.85, 3.35, 4.85, 6.35, 7.85, 9.35), labels = c("Obs", "SAVER", "Reg",
#                                                                "KNN", "SVD", "RF"))
axis(1, at = c(2, 3.5, 5, 6.5), labels = FALSE, pos = 0, lwd.ticks = 0)
axis(1, at = c(1.1, 7), labels = FALSE, lwd.ticks = 0, pos = 0)
text(c(2, 3.5, 5, 6.5), par()$usr[3]-0.002,
     labels = c("Baron", "Chen", "La Manno", "Zeisel"), srt = 45, adj = 1, 
     xpd = TRUE, cex = 1.4)
par(las = 1)
axis(2, pos = 1.1)
par(las = 0)
mtext("Cell-to-Cell", cex = 1.6, line = 0.2)

xloc <- seq(1.5, 2.5, 0.25)

fill <- c("white", "#cb181d", "#6baed6", "#bae4b3")

for (i in 1:4) {
  for (j in 1:4) {
    rect(xloc[j]+1.5*(i-1), 0, xloc[j+1]+1.5*(i-1), cmd[[2]][i, j],
         col = fill[j], lwd = 2)
  }
}

dev.off()

###############################################################################
## Cell clustering (takes a while to run)

# ### calculate PCs using jackstraw
# library(Seurat)
# 
# calc.jackstraw <- function(x, jack = TRUE) {
#   x.seurat <- CreateSeuratObject(raw.data = x)
#   x.seurat <- NormalizeData(x.seurat)
#   x.seurat <- ScaleData(x.seurat)
#   x.seurat <- FindVariableGenes(x.seurat, do.plot = FALSE)
#   x.seurat <- RunPCA(x.seurat, pc.genes = x.seurat@var.genes,
#                      pcs.compute = 25, do.print = FALSE)
#   if (jack) {
#     x.seurat <- JackStraw(x.seurat, num.pc = 25, do.print = TRUE)
#   }
#   return(x.seurat)
# }
# 
# dat.seurat <- vector("list", 4)
# names(dat.seurat) <- names(dat)
# 
# for (i in 1:4) {
#   dat.seurat[[i]] <- vector("list", 5)
#   names(dat.seurat[[i]]) <- c("Ref", "Obs", "SAVER", "MAGIC", "scImpute")
#   for (j in 1:5) {
#     dat.seurat[[i]][[j]] <- calc.jackstraw(dat[[i]][[j]], jack = FALSE)
#     print(j)
#   }
# }
# 
# 
# for (i in 1:4) {
#   for (j in 1:5) {
#     dat.seurat[[i]][[j]]@raw.data <- NULL
#   }
# }
# 
# jackstraw.results <- vector("list", 4)
# 
# for (i in 1:4) {
#   jackstraw.results[[i]] <- sapply(dat.seurat[[i]], function(x) 
#     levels(JackStrawPlot(x, PCs = 1:25)$data$PC.Score))
# }
# 
# ### run t-SNE at different resolutions
# pcs.chosen <- matrix(c(15, 5, 15, 10, 10,
#                        15, 10, 15, 10, 10,
#                        15, 10, 10, 5, 15,
#                        20, 10, 15, 10, 25), 4, 5, byrow = TRUE)
# rownames(pcs.chosen) <- names(dat)
# colnames(pcs.chosen) <- names(dat.seurat[[1]])
# 
# calc.tsne <- function(x, pcs, res, plot.tsne = FALSE,
#                       do.return = FALSE) {
#   x <- FindClusters(x, dims.use = 1:pcs, resolution = res,
#                     print.output = FALSE, save.SNN = TRUE)
#   if (is.null(x@dr$tsne)) {
#     x <- RunTSNE(x, dims.use = 1:pcs, check_duplicates = FALSE,
#                  do.fast = TRUE)
#   }
#   if (do.return) {
#     return(x)
#   }
#   if (plot.tsne) {
#     TSNEPlot(x, do.label = TRUE)
#   }
# }
# 
# sapply(dat.seurat, function(x) sapply(x, function(y) length(y@var.genes)))
# 
# dat.seurat[[1]][[1]] <- calc.tsne(dat.seurat[[1]][[1]], pcs = pcs.chosen[1, 1],
#                                   res = 0.7, do.return = TRUE)
# 
# dat.seurat[[2]][[1]] <- calc.tsne(dat.seurat[[2]][[1]], pcs = pcs.chosen[2, 1],
#                                   res = 0.6, do.return = TRUE)
# 
# dat.seurat[[3]][[1]] <- calc.tsne(dat.seurat[[3]][[1]], pcs = pcs.chosen[3, 1],
#                                   res = 1.1, do.return = TRUE)
# 
# dat.seurat[[4]][[1]] <- calc.tsne(dat.seurat[[4]][[1]], pcs = pcs.chosen[4, 1],
#                                   res = 0.8, do.return = TRUE)
# 
# res <- seq(0.4, 1.4, 0.1)
# ident <- vector("list", 4)
# tsne <- vector("list", 4)
# names(ident) <- names(dat)
# names(tsne) <- names(dat)
# for (i in 1:4) {
#   ident[[i]] <- vector("list", 5)
#   tsne[[i]] <- vector("list", 5)
#   names(ident[[i]]) <- c("Ref", "Obs", "SAVER", "MAGIC", "scImpute")
#   names(tsne[[i]]) <- c("Ref", "Obs", "SAVER", "MAGIC", "scImpute")
#   for (j in 2:5) {
#     ident[[i]][[j]] <- vector("list", 11)
#     for (k in 1:11) {
#       dat.seurat[[i]][[j]] <- calc.tsne(dat.seurat[[i]][[j]], pcs = pcs.chosen[i, j],
#                                         res = res[k], do.return = TRUE)
#       ident[[i]][[j]][[k]] <- dat.seurat[[i]][[j]]@ident
#     }
#     tsne[[i]][[j]] <- data.frame(dat.seurat[[i]][[j]]@dr$tsne@cell.embeddings)
#     print(paste(i, j))
#   }
# }
# 
# for (i in 1:4) {
#   ident[[i]][[1]] <- dat.seurat[[i]][[1]]@ident
#   tsne[[i]][[1]] <- data.frame(dat.seurat[[i]][[1]]@dr$tsne@cell.embeddings)
# }
# 
# ### remove singletons in Baron data
# ident <- info[[2]]
# 
# bad.ind <- which(as.numeric(ident[[1]][[1]]) > 7)
# ident[[1]][[1]] <- factor(ident[[1]][[1]][-bad.ind])
# for (i in 2:5) {
#   for (j in 1:11) {
#     ident[[1]][[i]][[j]] <- factor(ident[[1]][[i]][[j]][-bad.ind])
#   }
# }
# 
# for (i in 1:5) {
#   tsne[[1]][[i]] <- tsne[[1]][[i]][-bad.ind, ]
# }
# 
# saveRDS(list(tsne, ident), "SAVER-data/fig2d_tsne.rds")

library(clusteval)
library(scales)

info <- readRDS("SAVER-data/fig2d_tsne.rds")
tsne <- info[[1]]
ident <- info[[2]]
max.jaccard <- matrix(0, 4, 4)

for (i in 1:4) {
  for (j in 1:4) {
    max.jaccard[i, j] <- max(sapply(ident[[i]][[j+1]], function(x) 
      cluster_similarity(ident[[i]][[1]], x))[1:11])
  }
}

rownames(max.jaccard) <- names(dat)
colnames(max.jaccard) <- c("Obs", "SAVER", "MAGIC", "scImpute")

sapply(ident, function(x) length(levels(x[[1]])))

cols1 <-  c("#8dd3c7", "#ffd92f", "#e78ac3", "#bebada", "#80b1d3", "#fc8d62", 
            "#b3de69")

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

cols2 <- gg_color_hue(15)

alp <- c(0.4, 0.1, 0.4, 0.4)

m <- matrix(c(6, 6, 2, 2, 3, 3,
              1, 1, 2, 2, 3, 3, 
              1, 1, 4, 4, 5, 5,
              7, 7, 4, 4, 5, 5), 4, 6, byrow = TRUE)

i <- 4
png("plots/fig2d_tsne.png", 5.5, 4, units = "in", res = 300)
layout(m)
par(mar = c(0, 0, 0, 0), oma = c(4, 6, 4, 1))
cols <- cols1
plot(tsne[[i]][[1]], col = alpha(cols[as.numeric(ident[[i]][[1]])], alp[i]), pch = 19, axes = FALSE, 
     frame.plot = TRUE, ann = FALSE, cex = 0.6)
plot(-tsne[[i]][[2]], col = alpha(cols[as.numeric(ident[[i]][[1]])], alp[i]), pch = 19, axes = FALSE,
     frame.plot = TRUE, ann = FALSE, cex = 0.6)
legend("bottomright", legend = round(max.jaccard[i, 1], 2), bty = "n",
       cex = 1.4, y.intersp = 0.5)
plot(tsne[[i]][[3]], col = alpha(cols[as.numeric(ident[[i]][[1]])], alp[i]), pch = 19, axes = FALSE,
     frame.plot = TRUE, ann = FALSE, cex = 0.6)
legend("bottomright", legend = round(max.jaccard[i, 2], 2), bty = "n",
       cex = 1.4, y.intersp = 0.5)
plot(tsne[[i]][[4]], col = alpha(cols[as.numeric(ident[[i]][[1]])], alp[i]), pch = 19, axes = FALSE,
     frame.plot = TRUE, ann = FALSE, cex = 0.6)
legend("bottomright", legend = round(max.jaccard[i, 3], 2), bty = "n",
       cex = 1.4, y.intersp = 0.5)
plot(tsne[[i]][[5]], col = alpha(cols[as.numeric(ident[[i]][[1]])], alp[i]), pch = 19, axes = FALSE,
     frame.plot = TRUE, ann = FALSE, cex = 0.6)
legend("bottomright", legend = round(max.jaccard[i, 4], 2), bty = "n",
       cex = 1.4, y.intersp = 0.5)
mtext("Reference", outer = TRUE, at = 0.16, line = -4.5, cex = 1.2)
mtext("Observed", outer = TRUE, at = 0.5, line = 1, cex = 1.2)
mtext("SAVER", outer = TRUE, at = 0.83, line = 1, cex = 1.2)
mtext("MAGIC", outer = TRUE, side = 1, at = 0.5, line = 1.5, cex = 1.2)
mtext("scImpute", outer = TRUE, at = 0.84, side = 1, line = 1.5, cex = 1.2)
par(las = 1)
mtext("Zeisel", outer = TRUE, side = 2, at = 0.52, line = 1, cex = 1.2)
dev.off()

### supp t-SNE
png("plots/suppfig8_tsne.png", 11, 6, units = "in", res = 150)
par(mfrow = c(3, 5), mar = c(0, 0, 0, 0), oma = c(1, 6, 4, 1))
for (i in 1:3) {
  if (i == 2) {
    cols <- cols2
  } else {
    cols <- cols1
  }
  plot(tsne[[i]][[1]], col = alpha(cols[as.numeric(ident[[i]][[1]])], alp[i]), pch = 19, axes = FALSE, 
       frame.plot = TRUE, ann = FALSE, cex = 0.6)
  plot(tsne[[i]][[2]], col = alpha(cols[as.numeric(ident[[i]][[1]])], alp[i]), pch = 19, axes = FALSE,
       frame.plot = TRUE, ann = FALSE, cex = 0.6)
  legend("bottomright", legend = round(max.jaccard[i, 1], 2), bty = "n",
         cex = 1.4, y.intersp = 0.5)
  plot(tsne[[i]][[3]], col = alpha(cols[as.numeric(ident[[i]][[1]])], alp[i]), pch = 19, axes = FALSE,
       frame.plot = TRUE, ann = FALSE, cex = 0.6)
  legend("bottomright", legend = round(max.jaccard[i, 2], 2), bty = "n",
         cex = 1.4, y.intersp = 0.5)
  plot(tsne[[i]][[4]], col = alpha(cols[as.numeric(ident[[i]][[1]])], alp[i]), pch = 19, axes = FALSE,
       frame.plot = TRUE, ann = FALSE, cex = 0.6)
  legend("bottomright", legend = round(max.jaccard[i, 3], 2), bty = "n",
         cex = 1.4, y.intersp = 0.5)
  plot(tsne[[i]][[5]], col = alpha(cols[as.numeric(ident[[i]][[1]])], alp[i]), pch = 19, axes = FALSE,
       frame.plot = TRUE, ann = FALSE, cex = 0.6)
  legend("bottomright", legend = round(max.jaccard[i, 4], 2), bty = "n",
         cex = 1.4, y.intersp = 0.5)
}
mtext("Reference", outer = TRUE, at = 0.1, line = 1, cex = 1.2)
mtext("Observed", outer = TRUE, at = 0.3, line = 1, cex = 1.2)
mtext("SAVER", outer = TRUE, at = 0.5, line = 1, cex = 1.2)
mtext("MAGIC", outer = TRUE, at = 0.7, line = 1, cex = 1.2)
mtext("scImpute", outer = TRUE, at = 0.9, line = 1, cex = 1.2)
par(las = 1)
mtext("Baron", outer = TRUE, side = 2, at = 0.84, line = 1, cex = 1.2)
mtext("Chen", outer = TRUE, side = 2, at = 0.5, line = 1, cex = 1.2)
mtext("La", outer = TRUE, side = 2, at = 0.2, line = 2.4, cex = 1.2)
mtext("Manno", outer = TRUE, side = 2, at = 0.16, line = 1, cex = 1.2)
dev.off()




