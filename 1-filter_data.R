# Mo Huang, mohuang@wharton.upenn.edu
# Clean and filter datasets

###############################################################################
## Melanoma Drop-seq data

melanoma.raw <- as.matrix(read.table("SAVER-data/GSE99330_dropseqUPM.txt.gz"))

# convert upm to counts
melanoma <- sweep(melanoma.raw, 2, apply(melanoma.raw, 2, 
                                         function(x) min(x[x!= 0])), "/")
is.wholenumber <- function(x, tol = .Machine$double.eps^0.5) {
  abs(x - round(x)) < tol
}
sum(!is.wholenumber(df))

rm(melanoma.raw)
gc()

melanoma <- round(melanoma)

# filter data
melanoma.filt <- melanoma[which(rowMeans(melanoma) > 0.01), 
                          which(colSums(melanoma) >= 500 & 
                                  colSums(melanoma) <= 20000)]

saveRDS(melanoma.filt, "SAVER-data/melanoma_dropseq.rds")
write.csv(melanoma.filt, "SAVER-data/melanoma_dropseq.csv", quote = FALSE)

###############################################################################
## Baron data

x <- read.csv("SAVER-data/GSM2230757_human1_umifm_counts.csv.gz",
              header = TRUE, check.names = FALSE)
x.dat <- t(as.matrix(x[, 4:ncol(x)]))
colnames(x.dat) <- x$barcode
x <- x.dat

ercc <- which(grepl("ERCC", rownames(x), ignore.case = TRUE))
rgenes <- which("r_" == substring(rownames(x), 1, 2))
mt <- which("mt-" == substring(rownames(x), 1, 3))
x1 <- x[-c(ercc, rgenes, mt), ]

plot(density(log10(colSums(x1))))
summary(colSums(x1))
dim(x1)

# no need to filter by cells

# look at mean gene expression
plot(density(log10(rowMeans(x1))))
abline(v = log10(0.001))
# Remove genes with mean expression less than 0.001

x2 <- x1[rowMeans(x1) >= 0.001, ]
dim(x2)

# Look at nonzero cells
plot(density(log10(rowSums(x2 != 0))))
abline(v = log10(3))
# Remove genes with less than 3 non-zero cells

x3 <- x2[rowSums(x2 != 0) >= 3, ]
dim(x3)


# build reference
lib.size <- colSums(x3)
non.zero.prop <- apply(x3, 1, function(x) sum(x != 0)/length(x))

cells.filt <- which(lib.size > 5000)
genes.filt <- which(non.zero.prop > 0.25)

data.filt <- x3[genes.filt, cells.filt]

# Down-samp

set.seed(5)
alpha <- rgamma(ncol(data.filt), 10, 100)

data.samp <- t(apply(sweep(data.filt, 2, alpha, "*"), 1, function(x)
  rpois(length(x), x)))

colnames(data.samp) <- colnames(data.filt)

saveRDS(data.filt, "SAVER-data/baron_human_ref.rds")
saveRDS(data.samp, "SAVER-data/baron_human_samp.rds")

write.csv(data.samp, "SAVER-data/baron_human_samp.csv", quote = FALSE)

###############################################################################
## Chen data

x <- read.table("SAVER-data/GSE87544_Merged_17samples_14437cells_count.txt.gz",
                header = TRUE, check.names = FALSE, row.names = 1)
x <- as.matrix(x)

ercc <- which(grepl("ERCC", rownames(x), ignore.case = TRUE))
rgenes <- which("r_" == substring(rownames(x), 1, 2))
mt <- which("mt-" == substring(rownames(x), 1, 3))
x1 <- x[-c(ercc, rgenes, mt), ]

plot(density(log10(colSums(x1))))
summary(colSums(x1))
abline(v = log10(16000))
dim(x1)

# Filter out library size greater than 15000
x2 <- x1[, which(colSums(x1) <= 15000)]
summary(colSums(x2))

# look at mean gene expression
plot(density(log10(rowMeans(x2))))
abline(v = log10(0.0002))
# Remove genes with mean expression less than 0.0002

x3 <- x2[rowMeans(x2) >= 0.0002, ]
dim(x3)

# Look at nonzero cells
plot(density(log10(rowSums(x3 != 0))))
abline(v = log10(5))
# Remove genes with less than 5 non-zero cells

x4 <- x3[rowSums(x3 != 0) >= 5, ]
dim(x4)

# build reference
lib.size <- colSums(x4)
non.zero.prop <- apply(x4, 1, function(x) sum(x != 0)/length(x))

cells.filt <- which(lib.size > 2000)
genes.filt <- which(non.zero.prop > 0.2)

data.filt <- x4[genes.filt, cells.filt]

# Down-samp

set.seed(10)
alpha <- rgamma(ncol(data.filt), 10, 100)

data.samp <- t(apply(sweep(data.filt, 2, alpha, "*"), 1, function(x)
  rpois(length(x), x)))

colnames(data.samp) <- colnames(data.filt)

saveRDS(data.filt, "SAVER-data/chen_ref.rds")
saveRDS(data.samp, "SAVER-data/chen_samp.rds")

write.csv(data.samp, "SAVER-data/chen_samp.csv", quote = FALSE)


###############################################################################
## La Manno data

x <- read.table("SAVER-data/GSE76381_EmbryoMoleculeCounts.cef.txt.gz", skip = 5,
                header = FALSE, row.names = 1, check.names = FALSE)

cellnames <- read.table("SAVER-data/GSE76381_EmbryoMoleculeCounts.cef.txt.gz",
                        skip = 1, nrows = 1, row.names = 1, 
                        stringsAsFactors = FALSE)

x <- as.matrix(x)
colnames(x) <- cellnames

ercc <- which(grepl("ERCC", rownames(x), ignore.case = TRUE))
rgenes <- which("r_" == substring(rownames(x), 1, 2))
mt <- which("mt-" == substring(rownames(x), 1, 3))
x1 <- x[-c(ercc, rgenes, mt), ]


# look at library size

plot(density(log10(colSums(x1))))
summary(colSums(x1))
dim(x1)

# no need to filter by cells

# look at mean gene expression
plot(density(log10(rowMeans(x1))))
abline(v = log10(0.01))
# Remove genes with mean expression less than 0.001

x2 <- x1[rowMeans(x1) >= 0.001, ]
dim(x2)

# Look at nonzero cells
plot(density(log10(rowSums(x2 != 0))))
abline(v = log10(3))
# Remove genes with less than 3 non-zero cells

x3 <- x2[rowSums(x2 != 0) >= 3, ]
dim(x3)

lib.size <- colSums(x3)
non.zero.prop <- apply(x3, 1, function(x) sum(x != 0)/length(x))

cells.filt <- which(lib.size > 5000)
genes.filt <- which(non.zero.prop > 0.3)

data.filt <- x3[genes.filt, cells.filt]

set.seed(15)
alpha <- rgamma(ncol(data.filt), 10, 100)

data.samp <- t(apply(sweep(data.filt, 2, alpha, "*"), 1, function(x)
  rpois(length(x), x)))

colnames(data.samp) <- colnames(data.filt)

saveRDS(data.filt, "SAVER-data/manno_human_ref.rds")
saveRDS(data.samp, "SAVER-data/manno_human_samp.rds")

write.csv(data.samp, "SAVER-data/manno_human_samp.csv", quote = FALSE)


###############################################################################
## Zeisel data

x <- read.table("SAVER-data/expression_mRNA_17-Aug-2014.txt", skip = 11,
                header = FALSE, row.names = 1, check.names = FALSE)

x <- as.matrix(x[, -1])

cellnames <- read.table("SAVER-data/expression_mRNA_17-Aug-2014.txt",
                        skip = 7, nrows = 1, row.names = 1, 
                        stringsAsFactors = FALSE)
colnames(x) <- cellnames[-1]

# build reference

lib.size <- colSums(x)
non.zero.prop <- apply(x, 1, function(y) sum(y != 0)/length(y))

cells.filt <- which(lib.size > 10000)
genes.filt <- which(non.zero.prop > 0.4)

data.filt <- x[genes.filt, cells.filt]

# compare library size before and after
lib.size.filt <- colSums(data.filt)
plot(lib.size.filt, lib.size[cells.filt], xlab = "Filtered Library Size",
     ylab = "Original Library Size", 
     main = "Comparing Library Size Before and After Filtering")

# remove one outlier cell
data.filt2 <- data.filt[, -which.min(lib.size.filt)]

# downsample

set.seed(50)

n.cells <- ncol(data.filt2)

alpha1 <- rgamma(n.cells, 10, 40)
alpha2 <- rgamma(n.cells, 10, 100)
alpha3 <- rgamma(n.cells, 10, 200)


Y1 <- t(apply(sweep(data.filt2, 2, alpha1, "*"), 1, function(x) rpois(length(x), x)))
colnames(Y1) <- colnames(data.filt2)

Y2 <- t(apply(sweep(data.filt2, 2, alpha2, "*"), 1, function(x) rpois(length(x), x)))
colnames(Y2) <- colnames(data.filt2)

Y3 <- t(apply(sweep(data.filt2, 2, alpha3, "*"), 1, function(x) rpois(length(x), x)))
colnames(Y3) <- colnames(data.filt2)

saveRDS(data.filt2, "SAVER-data/zeisel_ref.rds")
saveRDS(Y3, "SAVER-data/zeisel_samp.rds")

write.csv(Y3, "SAVER-data/zeisel_samp.csv", quote = FALSE)



###############################################################################
## Hrvatin data

x <- read.csv("SAVER-data/GSE102827_merged_all_raw.csv.gz",
              header = TRUE, row.names = 1, check.names = FALSE)

x <- as.matrix(x)

# look at mean gene expression
plot(density(log10(rowMeans(x))))
# Filter out genes with expression less than 0.00003
x1 <- x[rowMeans(x) >= 0.00003, ]

# look at non-zero expression
# Filter out genes with non-zero expression in less than 4 cells

x2 <- x1[rowSums(x1 != 0) >= 4, ]
dim(x2)

set.seed(011118)
samp.cells <- sample(1:ncol(x2), 10000)

x.sub <- x2[, samp.cells]

saveRDS(x.sub, "SAVER-data/hrvatin.rds")

write.csv(x.sub, "SAVER-data/hrvatin.csv", quote = FALSE)


