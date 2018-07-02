# Mo Huang, mohuang@wharton.upenn.edu
# Run scImpute on datasets

# Install version 0.0.2
devtools::install_github("Vivianstats/scImpute", 
                         ref = "287535eca99044bcaa9930caf3e92a7ffa348726")
library(scImpute)

infile <- "csv"
outfile <- "csv"
out_dir <- "SAVER-data"
drop_thre <- 0.5
ncores <- 16

# Melanoma Drop-seq

count_path <- ("SAVER-data/melanoma_dropseq.csv")
fish <- read.table("SAVER-data/fishSubset.txt", header = TRUE, row.names = 1)

gene.names <- rownames(dat)
genes <- which(gene.names %in% colnames(fish))

scimpute(count_path = count_path, infile = infile, outfile = outfile,
         out_dir = out_dir, drop_thre = drop_thre, ncores = ncores)

out <- as.matrix(read.csv("SAVER-data/scimpute_count.csv", row.names = 1, 
                          check.names = FALSE))

saveRDS(out[genes, ], "SAVER-data/melanoma_dropseq_scimpute.rds")

# Baron

count_path <- ("SAVER-data/baron_human_samp.csv")

scimpute(count_path = count_path, infile = infile, outfile = outfile,
         out_dir = out_dir, drop_thre = drop_thre, ncores = ncores)

out <- as.matrix(read.csv("SAVER-data/scimpute_count.csv", row.names = 1,
                          check.names = FALSE))
saveRDS(out, "SAVER-data/baron_human_samp_scimpute.rds")

# Chen
count_path <- ("SAVER-data/chen_samp.csv")

scimpute(count_path = count_path, infile = infile, outfile = outfile,
         out_dir = out_dir, drop_thre = drop_thre, ncores = ncores)

out <- as.matrix(read.csv("SAVER-data/scimpute_count.csv", row.names = 1,
                          check.names = FALSE))
saveRDS(out, "SAVER-data/chen_samp_scimpute.rds")

# La Manno
count_path <- ("SAVER-data/manno_human_samp.csv")

scimpute(count_path = count_path, infile = infile, outfile = outfile,
         out_dir = out_dir, drop_thre = drop_thre, ncores = ncores)

out <- as.matrix(read.csv("SAVER-data/scimpute_count.csv", row.names = 1,
                          check.names = FALSE))
saveRDS(out, "SAVER-data/manno_human_samp_scimpute.rds")

# Zeisel
count_path <- ("SAVER-data/zeisel_samp.csv")

scimpute(count_path = count_path, infile = infile, outfile = outfile,
         out_dir = out_dir, drop_thre = drop_thre, ncores = ncores)

out <- as.matrix(read.csv("SAVER-data/scimpute_count.csv", row.names = 1,
                          check.names = FALSE))
saveRDS(out, "SAVER-data/zeisel_samp_scimpute.rds")
