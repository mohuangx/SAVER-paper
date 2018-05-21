# Mo Huang, mohuang@wharton.upenn.edu
# DE analysis

X <- readRDS("SAVER-data/zeisel_ref.rds")

Y <- readRDS("SAVER-data/zeisel_samp.rds")

n.cells <- ncol(X)
n.genes <- nrow(X)

cell.names <- colnames(X)
gene.names <- rownames(X)

normalizeData <- function(x, y = x) {
  sf <- colSums(y)/mean(colSums(y))
  return(sweep(x, 2, sf, "/"))
}

X.norm <- normalizeData(X)
Y.norm <- normalizeData(Y)

header <- read.table("SAVER-data/expression_mRNA_17-Aug-2014.txt", sep="\t", 
                     nrows=10, header=FALSE, fill = TRUE, comment.char = "", 
                     stringsAsFactors = FALSE)

header.df <- data.frame(t(header[9:10, -c(1, 2)]), stringsAsFactors = FALSE)
colnames(header.df) <- c("class1", "class2")
rownames(header.df) <- header[8, -c(1, 2)]

class1.cells <- header.df[colnames(X), 1]
class2.cells <- header.df[colnames(X), 2]

cell.subtype <- header.df[colnames(X), 2]
ca1pyr1 <- which(cell.subtype == "CA1Pyr1")
ca1pyr2 <- which(cell.subtype == "CA1Pyr2")

###############################################################################
### SAVER DE (takes a while to run)

# X.saver <- readRDS("SAVER-data/zeisel_ref_saver.rds")
# 
# Y.saver <- readRDS("SAVER-data/zeisel_samp_saver.rds")
# 
# get.sig.genes <- function(x, alpha, beta) {
#   cell.subtype <- header.df[colnames(x), 2]
#   ind <- which(cell.subtype %in% c("CA1Pyr1", "CA1Pyr2"))
#   group1 <- x[, which(cell.subtype == "CA1Pyr1")]
#   group2 <- x[, which(cell.subtype == "CA1Pyr2")]
#   w <- vector("list", 10)
#   for (i in 1:10) {
#     x.samp <- t(sapply(1:n.genes, function(i) rgamma(length(ind), alpha[i, ind],
#                                                      beta[i, ind])))
#     w[[i]] <- apply(x.samp, 1,
#                     FUN = function(z)
#                       unlist(
#                         wilcox.test(z ~ cell.subtype[
#                           ind])[1]))
#     print(i)
#   }
#   n.x <- ncol(group1)
#   n.y <- ncol(group2)
#   mean.w <- Reduce("+", w)/length(w) - n.x*n.y/2
#   var.w <- n.x*n.y/12*(n.x+n.y+1) + (1+1/10)/9*apply(simplify2array(w), 1, var)
#   mean.z <- sapply(1:n.genes, function(i) (mean.w[i]-sign(mean.w[i])*0.5)/sqrt(var.w[i]))
#   mean.p <- sapply(mean.z, function(x) 2*min(pnorm(x), pnorm(x, lower.tail = FALSE)))
#   pvaluesrk2 <- cbind(mean.w, mean.p)
#   pvaluesrk.adj <- p.adjust(pvaluesrk2[, 2], method = "BH")
#   group1.mean <- rowMeans(group1)
#   group2.mean <- rowMeans(group2)
#   fold.change <- log2(group1.mean/group2.mean)
#   df <- data.frame(genes = rownames(x),
#                    statistic = pvaluesrk2[, 1],
#                    variance = var.w,
#                    pvalue.rk = pvaluesrk2[, 2],
#                    pvalue.adj = pvaluesrk.adj,
#                    log2.fold.chg = fold.change, stringsAsFactors = FALSE)
#   return(df)
# }
# 
# X.saver.de <- get.sig.genes(X.saver$estimate, X.saver$alpha, X.saver$beta)
# Y.saver.de <- get.sig.genes(Y.saver$estimate, Y.saver$alpha, Y.saver$beta)
# 
# saveRDS(list(X.saver.de, Y.saver.de), "SAVER-data/de_saver.rds")

###############################################################################
### MAST

# de.cells <- which(class2.cells %in% c("CA1Pyr1", "CA1Pyr2"))
# 
# X.sub <- X.norm[, de.cells]
# Y.sub <- Y.norm[, de.cells]
# 
# cdat <- data.frame(class2 = class2.cells[de.cells])
# 
# library(MAST)
# library(data.table)
# 
# X.sca <- FromMatrix(log2(X.sub+1), cdat)
# cdr <- colSums(assay(X.sca) > 0)
# colData(X.sca)$cngeneson <- scale(cdr)
# zlm.X <- zlm.SingleCellAssay(~ class2 + cngeneson, X.sca)
# summary.X <- summary(zlm.X, doLRT='class2CA1Pyr2')
# summary.XDt <- summary.X$datatable
# fcHurdle.X <- merge(summary.XDt[contrast=='class2CA1Pyr2' & component=='H', 
#                                 .(primerid, `Pr(>Chisq)`)], #hurdle P values
#                     summary.XDt[contrast=='class2CA1Pyr2' & component=='logFC', 
#                                 .(primerid, coef, ci.hi, ci.lo)], by='primerid')
# fcHurdle.X[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
# X.sig <- data.frame(fcHurdle.X)
# 
# Y.sca <- FromMatrix(log2(Y.sub+1), cdat)
# cdr <- colSums(assay(Y.sca) > 0)
# colData(Y.sca)$cngeneson <- scale(cdr)
# zlm.Y <- zlm.SingleCellAssay(~ class2 + cngeneson, Y.sca)
# summary.Y <- summary(zlm.Y, doLRT='class2CA1Pyr2')
# summary.YDt <- summary.Y$datatable
# fcHurdle.Y <- merge(summary.YDt[contrast=='class2CA1Pyr2' & component=='H', 
#                                 .(primerid, `Pr(>Chisq)`)], #hurdle P values
#                     summary.YDt[contrast=='class2CA1Pyr2' & component=='logFC', 
#                                 .(primerid, coef, ci.hi, ci.lo)], by='primerid')
# fcHurdle.Y[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
# Y.sig <- data.frame(fcHurdle.Y)
# 
# saveRDS(list(X.sig, Y.sig), "SAVER-data/de_mast.rds")

###############################################################################
### scDD

# library(scDD)
# library(DESeq2)
# library(BiocParallel)
# library(Biobase)
# 
# register(MulticoreParam(16))
# 
# X <- readRDS("SAVER-data/zeisel_ref.rds")
# 
# Y <- readRDS("SAVER-data/zeisel_samp.rds")
# 
# header <- read.table("data/expression_mRNA_17-Aug-2014.txt", sep="\t", 
#                      nrows=10, header=FALSE, fill = TRUE, comment.char = "", 
#                      stringsAsFactors = FALSE)
# 
# header.names <- header[, 2]
# 
# header <- header[, -c(1,2)]
# 
# cell.subtype <- as.character(header[10, which(header[8, ] %in% colnames(X))])
# 
# # only selecting CA1Pyr1 and CA1Pyr2
# ind <- which(cell.subtype %in% c("CA1Pyr1", "CA1Pyr2"))
# 
# normalizeData <- function(x, y = x) {
#   sf <- colSums(y)/mean(colSums(y))
#   return(sweep(x, 2, sf, "/"))
# }
# 
# X.norm <- normalizeData(X)
# Y.norm <- normalizeData(Y)
# 
# X1 <- X.norm[, ind]
# Y1 <- Y.norm[, ind]
# 
# cell.subtype1 <- as.character(header[10, which(header[8, ] %in% colnames(X1))])
# 
# columnData <- data.frame(group = cell.subtype1, row.names = colnames(X1))
# 
# phenoData <- data.frame(condition = ifelse(columnData$group == "CA1Pyr1", 1, 2),
#                         row.names = colnames(X1))
# phenoData <- as(phenoData, "AnnotatedDataFrame")
# 
# X.se <- ExpressionSet(assayData = as.matrix(X1),
#                       phenoData = phenoData)
# 
# Y.se <- ExpressionSet(assayData = as.matrix(Y1),
#                       phenoData = phenoData) 
# 
# prior_param <- list(alpha=0.01, mu0=0, s0=0.01, a0=0.01, b0=0.01)
# 
# X.sig <- scDD(X.se, prior_param=prior_param, n.cores = cores)
# Y.sig <- scDD(Y.se, prior_param=prior_param, n.cores = cores)
# 
# saveRDS(list(X.sig, Y.sig), "SAVER-data/de_scdd.rds")

###############################################################################
### SCDE

# library(scde)
# library(methods)
# 
# X.sub <- X[, c(ca1pyr1, ca1pyr2)]
# Y.sub <- Y[, c(ca1pyr1, ca1pyr2)]
# 
# sg <- factor(header.df[colnames(X.sub), "class2"], 
#              levels = c("CA1Pyr1", "CA1Pyr2"))
# names(sg) <- colnames(X.sub)
# 
# X.err <- scde.error.models(counts = X.sub, groups = sg, n.cores = 1, 
#                            threshold.segmentation = TRUE, 
#                            save.crossfit.plots = FALSE, 
#                            save.model.plots = FALSE,
#                            verbose = 1)
# 
# Y.err <- scde.error.models(counts = Y.sub, groups = sg, n.cores = 1, 
#                            threshold.segmentation = TRUE, 
#                            save.crossfit.plots = FALSE, 
#                            save.model.plots = FALSE,
#                            verbose = 1)
# 
# valid.cells.X <- X.err$corr.a > 0
# valid.cells.Y <- Y.err$corr.a > 0
# 
# X.err <- X.err[valid.cells.X, ]
# Y.err <- Y.err[valid.cells.Y, ]
# 
# X.prior <- scde.expression.prior(models = X.err, counts = X.sub, 
#                                  length.out = 400, show.plot = FALSE)
# 
# Y.prior <- scde.expression.prior(models = Y.err, counts = Y.sub, 
#                                  length.out = 400, show.plot = FALSE)
# 
# 
# X.diff <- scde.expression.difference(X.err, X.sub, X.prior, groups  =  sg, 
#                                      n.randomizations  =  100, n.cores  =  1,
#                                      verbose = 1)
# 
# Y.diff <- scde.expression.difference(Y.err, Y.sub, Y.prior, groups  =  sg, 
#                                      n.randomizations  =  100, n.cores  =  1,
#                                      verbose = 1)
# 
# 
# X.diff$p.value <- 2*pnorm(abs(X.diff$cZ), lower.tail = FALSE)
# Y.diff$p.value <- 2*pnorm(abs(Y.diff$cZ), lower.tail = FALSE)
# 
# saveRDS(list(X.diff, Y.diff), "SAVER-data/de_scde.rds")

###############################################################################
### FDR estimation

# library(doParallel)
# 
# cl <- makeCluster(20)
# registerDoParallel(cl)
# 
# X <- readRDS("SAVER-data/zeisel_ref.rds")
# 
# Y <- readRDS("SAVER-data/zeisel_samp.rds")
# 
# n.cells <- ncol(X)
# n.genes <- nrow(X)
# 
# cell.names <- colnames(X)
# gene.names <- rownames(X)
# 
# normalizeData <- function(x, y = x) {
#   sf <- colSums(y)/mean(colSums(y))
#   return(sweep(x, 2, sf, "/"))
# }
# 
# X.norm <- normalizeData(X)
# Y.norm <- normalizeData(Y)
# 
# 
# header <- read.table("data/expression_mRNA_17-Aug-2014.txt", sep="\t", 
#                      nrows=10, header=FALSE, fill = TRUE, comment.char = "", 
#                      stringsAsFactors = FALSE)
# 
# header.df <- data.frame(t(header[9:10, -c(1, 2)]), stringsAsFactors = FALSE)
# colnames(header.df) <- c("class1", "class2")
# rownames(header.df) <- header[8, -c(1, 2)]
# 
# class1.cells <- header.df[colnames(X), 1]
# class2.cells <- header.df[colnames(X), 2]
# 
# cell.subtype <- header.df[colnames(X), 2]
# ca1pyr1 <- which(cell.subtype == "CA1Pyr1")
# ca1pyr2 <- which(cell.subtype == "CA1Pyr2")
# 
# de.cells <- which(class2.cells %in% c("CA1Pyr1", "CA1Pyr2"))
# 
# X.sub <- X[, de.cells]
# X.sub.norm <- X.norm[, de.cells]
# 
# Y.sub <- Y[, de.cells]
# Y.sub.norm <- Y.norm[, de.cells]
# 
# package.list <- c("MAST", "data.table", "scde", "methods", "scDD", "DESeq2",
#                   "Biobase")
# X.fp <- foreach(i = 1:20, .packages = package.list,
#               .errorhandling = 'pass') %dopar% {
#   set.seed(i)
#   de.cells.perm <- sample(de.cells, length(de.cells))
#   x.samp <- t(sapply(1:n.genes, function(k) rgamma(length(de.cells), 
#                                                    X.saver$alpha[k, de.cells],
#                                                    X.saver$beta[k, de.cells])))
#   w <- apply(x.samp, 1,
#              FUN = function(z)
#                unlist(
#                  wilcox.test(z ~ cell.subtype[
#                    de.cells.perm])[1]))
#   n.x <- ncol(group1)
#   n.y <- ncol(group2)
#   mean.w <- w-n.x*n.y/2
#   var.w <- rep(n.x*n.y/12*(n.x+n.y+1), n.genes)
#   mean.z <- sapply(1:n.genes, function(i) (mean.w[i]-sign(mean.w[i])*0.5)/sqrt(var.w[i]))
#   mean.p <- sapply(mean.z, function(x) 2*min(pnorm(x), pnorm(x, lower.tail = FALSE)))
#   pvaluesrk2 <- cbind(mean.w, mean.p)
#   pvaluesrk.adj <- p.adjust(pvaluesrk2[, 2], method = "BH")
#   group1.mean <- rowMeans(group1)
#   group2.mean <- rowMeans(group2)
#   fold.change <- log2(group1.mean/group2.mean)
#   X.sig.saver <- data.frame(genes = rownames(x),
#                             statistic = pvaluesrk2[, 1],
#                             variance = var.w,
#                             pvalue.rk = pvaluesrk2[, 2],
#                             pvalue.adj = pvaluesrk.adj,
#                             log2.fold.chg = fold.change, stringsAsFactors = FALSE)
#   
#   cdat <- data.frame(class2 = class2.cells[de.cells.perm])
#   Y.sca <- FromMatrix(log2(X.sub.norm+1), cdat)
#   cdr <- colSums(assay(Y.sca) > 0)
#   colData(Y.sca)$cngeneson <- scale(cdr)
#   zlm.Y <- zlm.SingleCellAssay(~ class2 + cngeneson, Y.sca)
#   summary.Y <- summary(zlm.Y, doLRT='class2CA1Pyr2')
#   summary.YDt <- summary.Y$datatable
#   fcHurdle.Y <- merge(summary.YDt[contrast=='class2CA1Pyr2' & component=='H', 
#                                   .(primerid, `Pr(>Chisq)`)], #hurdle P values
#                       summary.YDt[contrast=='class2CA1Pyr2' & component=='logFC', 
#                                   .(primerid, coef, ci.hi, ci.lo)], by='primerid')
#   fcHurdle.Y[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
#   X.sig.mast <- data.frame(fcHurdle.Y)
#   
#   cell.subtype1 <- class2.cells[de.cells.perm]
#   columnData <- data.frame(group = cell.subtype1, 
#                            row.names = colnames(X.sub.norm))
#   phenoData <- data.frame(condition = ifelse(columnData$group == "CA1Pyr1", 1, 2),
#                           row.names = colnames(X.sub.norm))
#   phenoData <- as(phenoData, "AnnotatedDataFrame")
#   X.se <- ExpressionSet(assayData = as.matrix(X.sub.norm),
#                         phenoData = phenoData)
#   prior_param <- list(alpha=0.01, mu0=0, s0=0.01, a0=0.01, b0=0.01)
#   X.sig.scdd <- scDD(X.se, prior_param=prior_param, n.cores = 1)
#   
#   sg <- factor(class2.cells[de.cells.perm], 
#                levels = c("CA1Pyr1", "CA1Pyr2"))
#   names(sg) <- colnames(X.sub)
#   X.err <- scde.error.models(counts = X.sub, groups = sg, n.cores = 1, 
#                              threshold.segmentation = TRUE, 
#                              save.crossfit.plots = FALSE, 
#                              save.model.plots = FALSE)
#   valid.cells.X <- X.err$corr.a > 0
#   X.err <- X.err[valid.cells.X, ]
#   X.prior <- scde.expression.prior(models = X.err, counts = X.sub, 
#                                    length.out = 400, show.plot = FALSE)
#   X.sig.scde <- scde.expression.difference(X.err, X.sub, X.prior, groups  =  sg, 
#                                            n.randomizations  =  100, n.cores  =  1,
#                                            verbose = 1)
#   X.sig.scde$p.value <- 2*pnorm(abs(X.sig.scde$cZ), lower.tail = FALSE)
#   return(list(saver = X.sig.saver, mast = X.sig.mast, 
#               scdd = X.sig.scdd,
#               scde = X.sig.scde))
# }
# 
# Y.fp <- foreach(i = 1:20, .packages = package.list,
#                 .errorhandling = 'pass') %dopar% {
#   set.seed(i)
#   de.cells.perm <- sample(de.cells, length(de.cells))
#   x.samp <- t(sapply(1:n.genes, function(k) rgamma(length(de.cells), 
#                                                    Y.saver$alpha[k, de.cells],
#                                                    Y.saver$beta[k, de.cells])))
#   w <- apply(x.samp, 1,
#              FUN = function(z)
#                unlist(
#                  wilcox.test(z ~ cell.subtype[
#                    de.cells.perm])[1]))
#   n.x <- ncol(group1)
#   n.y <- ncol(group2)
#   mean.w <- w-n.x*n.y/2
#   var.w <- rep(n.x*n.y/12*(n.x+n.y+1), n.genes)
#   mean.z <- sapply(1:n.genes, function(i) (mean.w[i]-sign(mean.w[i])*0.5)/sqrt(var.w[i]))
#   mean.p <- sapply(mean.z, function(x) 2*min(pnorm(x), pnorm(x, lower.tail = FALSE)))
#   pvaluesrk2 <- cbind(mean.w, mean.p)
#   pvaluesrk.adj <- p.adjust(pvaluesrk2[, 2], method = "BH")
#   group1.mean <- rowMeans(group1)
#   group2.mean <- rowMeans(group2)
#   fold.change <- log2(group1.mean/group2.mean)
#   Y.sig.saver <- data.frame(genes = rownames(x),
#                             statistic = pvaluesrk2[, 1],
#                             variance = var.w,
#                             pvalue.rk = pvaluesrk2[, 2],
#                             pvalue.adj = pvaluesrk.adj,
#                             log2.fold.chg = fold.change, stringsAsFactors = FALSE)
#   
#   cdat <- data.frame(class2 = class2.cells[de.cells.perm])
#   Y.sca <- FromMatrix(log2(Y.sub.norm+1), cdat)
#   cdr <- colSums(assay(Y.sca) > 0)
#   colData(Y.sca)$cngeneson <- scale(cdr)
#   zlm.Y <- zlm.SingleCellAssay(~ class2 + cngeneson, Y.sca)
#   summary.Y <- summary(zlm.Y, doLRT='class2CA1Pyr2')
#   summary.YDt <- summary.Y$datatable
#   fcHurdle.Y <- merge(summary.YDt[contrast=='class2CA1Pyr2' & component=='H', 
#                                   .(primerid, `Pr(>Chisq)`)], #hurdle P values
#                       summary.YDt[contrast=='class2CA1Pyr2' & component=='logFC', 
#                                   .(primerid, coef, ci.hi, ci.lo)], by='primerid')
#   fcHurdle.Y[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
#   Y.sig.mast <- data.frame(fcHurdle.Y)
#   
#   cell.subtype1 <- class2.cells[de.cells.perm]
#   columnData <- data.frame(group = cell.subtype1, 
#                            row.names = colnames(Y.sub.norm))
#   phenoData <- data.frame(condition = ifelse(columnData$group == "CA1Pyr1", 1, 2),
#                           row.names = colnames(Y.sub.norm))
#   phenoData <- as(phenoData, "AnnotatedDataFrame")
#   Y.se <- ExpressionSet(assayData = as.matrix(Y.sub.norm),
#                         phenoData = phenoData)
#   prior_param <- list(alpha=0.01, mu0=0, s0=0.01, a0=0.01, b0=0.01)
#   Y.sig.scdd <- scDD(Y.se, prior_param=prior_param, n.cores = 1)
#   
#   sg <- factor(class2.cells[de.cells.perm], 
#                levels = c("CA1Pyr1", "CA1Pyr2"))
#   names(sg) <- colnames(Y.sub)
#   Y.err <- scde.error.models(counts = Y.sub, groups = sg, n.cores = 1, 
#                              threshold.segmentation = TRUE, 
#                              save.crossfit.plots = FALSE, 
#                              save.model.plots = FALSE)
#   valid.cells.Y <- Y.err$corr.a > 0
#   Y.err <- Y.err[valid.cells.Y, ]
#   Y.prior <- scde.expression.prior(models = Y.err, counts = Y.sub, 
#                                    length.out = 400, show.plot = FALSE)
#   Y.sig.scde <- scde.expression.difference(Y.err, Y.sub, Y.prior, groups  =  sg, 
#                                            n.randomizations  =  100, n.cores  =  1,
#                                            verbose = 1)
#   Y.sig.scde$p.value <- 2*pnorm(abs(Y.sig.scde$cZ), lower.tail = FALSE)
#   return(list(saver = Y.sig.saver, mast = Y.sig.mast, 
#               scdd = Y.sig.scdd,
#               scde = Y.sig.scde))
# }
# 
# saveRDS(list(X.fp, Y.fp), "SAVER-data/de_fdr.rds")

###############################################################################
### DE analysis

de.saver <- readRDS("SAVER-data/de_saver.rds")
de.mast <- readRDS("SAVER-data/de_mast.rds")
de.scdd <- readRDS("SAVER-data/de_scdd.rds")
de.scde <- readRDS("SAVER-data/de_scde.rds")
perm <- readRDS("SAVER-data/de_fdr.rds")

X.sig <- c(sum(de.saver[[1]]$pvalue.adj < 0.01), 
           sum(de.mast[[1]]$fdr < 0.01), 
           sum(de.scde[[1]]$p.values < 0.01), 
           sum(de.scdd[[1]]$Genes$nonzero.pvalue.adj < 0.01, na.rm = TRUE),
           sum(de.scdd[[1]]$Genes$zero.pvalue.adj < 0.01, na.rm = TRUE))
names(X.sig) <- c("SAVER", "MAST", "SCDE", "SCDD1", "SCDD2")

Y.sig <- c(sum(de.saver[[2]]$pvalue.adj < 0.01), 
           sum(de.mast[[2]]$fdr < 0.01), 
           sum(de.scde[[2]]$p.value < 0.01), 
           sum(de.scdd[[2]]$Genes$nonzero.pvalue.adj < 0.01, na.rm = TRUE),
           sum(de.scdd[[2]]$Genes$zero.pvalue.adj < 0.01, na.rm = TRUE))
names(Y.sig) <- c("SAVER", "MAST", "SCDE", "SCDD1", "SCDD2")

total.sig <- rbind(X.sig, Y.sig)

scdd.n2 <- c(sum(!is.na(de.scdd[[1]]$Genes$zero.pvalue)), 
             sum(!is.na(de.scdd[[2]]$Genes$zero.pvalue)))
pvalue.cut <- 0.01/3529*total.sig
pvalue.cut[, 5] <- pvalue.cut[, 5]*3529/scdd.n2

# remove bad runs in down-sampled
bad.run <- which(sapply(perm[[2]], function(x) is.null(x[[2]])))
perm[[2]] <- perm[[2]][-bad.run]

fp <- matrix(0, 2, 5)
for (i in 1:2) {
  fp[i, 1] <- mean(sapply(perm[[i]], function(x) sum(x[[1]]$pvalue.rk < pvalue.cut[i, 1])))
  fp[i, 2] <- mean(sapply(perm[[i]], function(x) sum(x[[2]]$Pr..Chisq. < pvalue.cut[i, 2])))
  fp[i, 3] <- mean(sapply(perm[[i]], function(x) sum(abs(x[[4]]$Z) > qnorm(pvalue.cut[i, 3]/2, 
                                                                           lower.tail = FALSE))))
  fp[i, 4] <- mean(sapply(perm[[i]], function(x) sum(x[[3]]$Genes$nonzero.pvalue < pvalue.cut[i, 4], 
                                                     na.rm = TRUE)))
  fp[i, 5] <- mean(sapply(perm[[i]], function(x) sum(x[[3]]$Genes$zero.pvalue < pvalue.cut[i, 5],
                                                     na.rm = TRUE)))
}

fdr <- fp/total.sig
fdr[, 4] <- (fdr[, 4]*total.sig[, 4] + fdr[, 5]*total.sig[, 5])/(total.sig[, 4] + total.sig[, 5])

total.sig[, 4] <- total.sig[, 4] + total.sig[, 5]

pdf("plots/fig2c_de.pdf", 6.5, 3.25)

par(mfrow = c(1, 2), cex.main = 1.5, mar = c(2, 3, 0, 1.5) + 0.1, oma = c(2, 2, 2, 0), 
    mgp = c(3.5, 1, 0),
    cex.lab = 1.5, font.lab = 2, cex.axis = 1.3, bty = "n", las = 1)

plot(0, type = "n", ylab = " ", xlab = " ", cex = 1.5, 
     ylim = c(0, 3500), xlim = c(1, 7.5), lwd = 2, pch = 5, axes = FALSE, main = " ")
# axis(1, at = c(1.85, 3.35, 4.85, 6.35, 7.85, 9.35), labels = c("Obs", "SAVER", "Reg",
#                                                                "KNN", "SVD", "RF"))
axis(1, at = c(3, 6.5), labels = FALSE, pos = 0, lwd.ticks = 0)
axis(1, at = c(1.1, 7.5), labels = FALSE, lwd.ticks = 0, pos = 0)
text(c(2.75, 5.75), par()$usr[3]-100,
     labels = c("Reference", "Observed"), xpd = TRUE, cex = 1.2)
par(las = 1)
axis(2, pos = 1.1)
par(las = 0)
mtext("Significant genes", side = 2, line = 3.5, cex = 1.3)


xloc <- seq(1.75, 3.75, 0.5)

fill <- c("#a50f15", "#3182bd", "#74c476", "#cbc9e2")

for (i in 1:2) {
  for (j in 1:4) {
    j. <- c(1, 2, 4, 3)[j]
    rect(xloc[j]+3*(i-1), 0, xloc[j+1]+3*(i-1), total.sig[i, j.],
         col = fill[j], lwd = 2)
  }
}

plot(0, type = "n", ylab = " ", xlab = " ", cex = 1.5, 
     ylim = c(0, 0.03), xlim = c(1, 7.5), lwd = 2, pch = 5, axes = FALSE, main = " ")
# axis(1, at = c(1.85, 3.35, 4.85, 6.35, 7.85, 9.35), labels = c("Obs", "SAVER", "Reg",
#                                                                "KNN", "SVD", "RF"))
axis(1, at = c(2, 6), labels = FALSE, pos = 0, lwd.ticks = 0)
axis(1, at = c(1.1, 7.5), labels = FALSE, lwd.ticks = 0, pos = 0)
text(c(2.75, 5.75), par()$usr[3]-0.001,
     labels = c("Reference", "Observed"), xpd = TRUE, cex = 1.2)
par(las = 1)
axis(2, pos = 1.1)
par(las = 0)
mtext("Estimated FDR", side = 2, line = 4, cex = 1.3)


xloc <- seq(1.75, 3.75, 0.5)

fill <- c("#a50f15", "#3182bd", "#74c476", "#cbc9e2")

segments(1.1, 0.01, 7.75, 0.01, lwd = 2, col = "darkgray", lty = 2)

for (i in 1:2) {
  for (j in 1:4) {
    j. <- c(1, 2, 4, 3)[j]
    rect(xloc[j]+3*(i-1), 0, xloc[j+1]+3*(i-1), fdr[i, j.],
         col = fill[j], lwd = 2)
  }
}


par(mfrow = c(1, 1))
plot(0, type = "n", xlab = "", ylab = "", axes = FALSE)

legend("center", c("SAVER", "MAST", "scDD", "SCDE"), pch = 22,
       pt.bg = fill, cex = 1.2, ncol = 3, text.font = 1,
       pt.cex = 2, pt.lwd = 2, xjust = 0.5, yjust = 0,
       xpd = TRUE, bty = "n")

dev.off()







