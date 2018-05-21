# Mo Huang, mohuang@wharton.upenn.edu
# Run SAVER on datasets

library(SAVER)
library(doParallel)

cl <- makeCluster(16, outfile = "")
registerDoParallel(cl)

# Melanoma Drop-seq

dat <- readRDS("SAVER-data/melanoma_dropseq.rds")
fish <- read.table("SAVER-data/fishSubset.txt", header = TRUE, row.names = 1)

gene.names <- rownames(dat)
genes <- which(gene.names %in% colnames(fish))

out <- saver(dat, pred.genes = genes, pred.genes.only = TRUE,
             do.fast = FALSE)

saveRDS(out, "SAVER-data/melanoma_dropseq_saver.rds")

# Baron

dat <- readRDS("SAVER-data/baron_human_samp.rds")

out <- saver(dat, do.fast = FALSE)

saveRDS(out, "SAVER-data/baron_human_samp_saver.rds")

# Chen
dat <- readRDS("SAVER-data/chen_samp.rds")

out <- saver(dat, do.fast = FALSE)

saveRDS(out, "SAVER-data/chen_samp_saver.rds")

# La Manno
dat <- readRDS("SAVER-data/manno_human_samp.rds")

out <- saver(dat, do.fast = FALSE)

saveRDS(out, "SAVER-data/manno_human_samp_saver.rds")

# Zeisel
dat <- readRDS("SAVER-data/zeisel_samp.rds")

out <- saver(dat, do.fast = FALSE)

saveRDS(out, "SAVER-data/zeisel_samp_saver.rds")

# Hrvatin
dat <- readRDS("SAVER-data/hrvatin.rds")

out <- saver(dat, do.fast = FALSE)

saveRDS(out, "SAVER-data/hrvatin_saver.rds")


