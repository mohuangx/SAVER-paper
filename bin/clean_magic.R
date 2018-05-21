# Mo Huang, mohuang@wharton.upenn.edu
# Clean MAGIC output

# Melanoma Drop-seq
dat <- as.matrix(read.csv("SAVER-data/melanoma_dropseq_magic.csv", 
                          row.names = 1, header = TRUE, 
                          check.names = FALSE))
fish <- read.table("data/fishSubset.txt", header = TRUE, row.names = 1)

gene.names <- colnames(dat)
genes <- which(gene.names %in% colnames(fish))

saveRDS(t(dat[, genes]), "SAVER-data/melanoma_dropseq_magic.rds")

# Baron
dat <- t(as.matrix(read.csv("SAVER-data/baron_human_samp_magic.csv", 
                            row.names = 1, header = TRUE, 
                            check.names = FALSE)))

saveRDS(dat, "SAVER-data/baron_human_samp_magic.rds")

# Chen
dat <- t(as.matrix(read.csv("SAVER-data/chen_samp_magic.csv", 
                            row.names = 1, header = TRUE, 
                            check.names = FALSE)))

saveRDS(dat, "SAVER-data/chen_samp_magic.rds")

# La Manno
dat <- t(as.matrix(read.csv("SAVER-data/manno_human_samp_magic.csv", 
                            row.names = 1, header = TRUE, 
                            check.names = FALSE)))

saveRDS(dat, "SAVER-data/manno_human_samp_magic.rds")

# Zeisel
dat <- t(as.matrix(read.csv("SAVER-data/zeisel_samp_magic.csv", 
                            row.names = 1, header = TRUE, 
                            check.names = FALSE)))

saveRDS(dat, "SAVER-data/zeisel_samp_magic.rds")



