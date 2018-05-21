# Mo Huang, mohuang@wharton.upenn.edu
# Run scImpute on datasets

library(scImpute)

infile <- "rds"
outfile <- "rds"
out_dir <- "SAVER-data"
drop_thre <- 0.5
k <- 1
ncores <- 16

# Melanoma Drop-seq

count_path <- ("SAVER-data/melanoma_dropseq.rds")
fish <- read.table("SAVER-data/fishSubset.txt", header = TRUE, row.names = 1)

gene.names <- rownames(dat)
genes <- which(gene.names %in% colnames(fish))

scimpute(count_path = count_path, infile = infile, outfile = outfile,
         out_dir = out_dir, labeled = FALSE, drop_thre = drop_thre,
         Kcluster = k, ncores = ncores)

out <- readRDS("SAVER-data/scimpute_count.rds")

saveRDS(out[genes, ], "SAVER-data/melanoma_dropseq_scimpute.rds")

# Baron

count_path <- ("SAVER-data/baron_human_samp.rds")

scimpute(count_path = count_path, infile = infile, outfile = outfile,
         out_dir = out_dir, labeled = FALSE, drop_thre = drop_thre,
         Kcluster = k, ncores = ncores)

out <- readRDS("SAVER-data/scimpute_count.rds")
saveRDS(out, "SAVER-data/baron_human_samp_scimpute.rds")

# Chen
count_path <- ("SAVER-data/chen_samp.rds")

scimpute(count_path = count_path, infile = infile, outfile = outfile,
         out_dir = out_dir, labeled = FALSE, drop_thre = drop_thre,
         Kcluster = k, ncores = ncores)

out <- readRDS("SAVER-data/scimpute_count.rds")
saveRDS(out, "SAVER-data/chen_samp_scimpute.rds")

# La Manno
count_path <- ("SAVER-data/manno_human_samp.rds")

scimpute(count_path = count_path, infile = infile, outfile = outfile,
         out_dir = out_dir, labeled = FALSE, drop_thre = drop_thre,
         Kcluster = k, ncores = ncores)

out <- readRDS("SAVER-data/scimpute_count.rds")
saveRDS(out, "SAVER-data/manno_human_samp_scimpute.rds")

# Zeisel
count_path <- ("SAVER-data/zeisel_samp.rds")

scimpute(count_path = count_path, infile = infile, outfile = outfile,
         out_dir = out_dir, labeled = FALSE, drop_thre = drop_thre,
         Kcluster = k, ncores = ncores)

out <- readRDS("SAVER-data/scimpute_count.rds")
saveRDS(out, "SAVER-data/zeisel_samp_scimpute.rds")
