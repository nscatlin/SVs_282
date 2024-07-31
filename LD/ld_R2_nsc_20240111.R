#### arguments for command line ----

args <- commandArgs(trailingOnly = TRUE)

# # help
# if (all(length(args) == 1 & args == "-h" | args == "--help")) {
#   cat("
# Description: this script merges SNPs and TEs hapmap files from usda parents to be used in
#              Tassel 5 when projecting TEs into RILs
#
# Usage: ")
#   quit()
# }

# make sure the correct number of arguments are used
# you should provide 3 arguments
if (length(args) != 3) {
  stop("incorrect number of arguments provided.

       Usage: Rscript
       ")
}

# assign arguments to variables
plink.file <- args[1]
out.dir.ld <- args[2]
chr <- args[3]
# 
# plink.file <- "~/Downloads/plink_results_SNPs-highest-LD-TE_R2_chr3.ld"
# out.dir.ld <- "~/Downloads/ld_chr3_test"
# chr <- "3"
# 



# setwd("~/projects/ld_snp-te/analysis/")
# plink.file <- "WiDiv508_ref.B73_SNPs_not.in.TEs_allchr.TEs-combined.SNP-TE-only.chr1.ld"
# out.dir.ld <- "ld_distribution"
# chr <- 10

#### libraries ----

if(!require("data.table")) install.packages("data.table")
if(!require("doParallel")) install.packages("doParalell")
if(!require("tidyverse")) install.packages("tidyverse")


# if (detectCores() > 10) {
#   num.cores <- 10
# } else {
#   num.cores <- detectCores()
# }

if (!dir.exists(out.dir.ld)) dir.create(out.dir.ld)

#### subsample ----

# load one chr at a time
plink.file.chr <- gsub("chr[0-9]+", paste0("chr", chr), plink.file, perl = TRUE)
# open table with LD among markers
LD_results_dups <- fread(plink.file.chr, header = TRUE, data.table = FALSE)

LD_results <- LD_results_dups %>%
  filter(!(str_detect(SNP_B, "_ins_")))

LD_results$dist_to_te <- LD_results[, "BP_B"] - LD_results[, "BP_A"]

LD_results_R2 <- LD_results %>%
  select(SNP_A, R2, dist_to_te)

# LD_results_max_R2_nodups_old <- LD_results %>%
#   group_by(SNP_A) %>%
#   slice_max(R2) %>%
#   slice_min(abs(dist_to_te)) %>%
#   ungroup()

LD_results_max_R2_nodups <- LD_results %>%
  arrange(SNP_A, dist_to_te) %>%
  group_by(SNP_A) %>%
  slice_max(R2) %>%
  slice_min(dist_to_te) %>%
  group_by(SNP_A) %>%
  sample_n(1) %>%
  ungroup()

LD_results_min_R2_nodups <- LD_results %>% 
  arrange(SNP_A, dist_to_te) %>%
  group_by(SNP_A) %>%
  slice_min(R2) %>%
  slice_min(dist_to_te) %>%
  group_by(SNP_A) %>%
  sample_n(1) %>%
  ungroup()

# LD_results_min_R2_nodups <- LD_results %>% 
#   group_by(SNP_A) %>% 
#   slice_min(R2) %>% 
#   slice_min(abs(dist_to_te)) %>% 
#   ungroup()



# LD_results_max_R2_nodups_old[LD_results_max_R2_nodups_old$SNP_A == "B73_ins_chr3:206305393-206306871_Oh43:205103668_36532",]
# LD_results_max_R2_nodups_new[LD_results_max_R2_nodups_new$SNP_A == "B73_ins_chr3:206305393-206306871_Oh43:205103668_36532",]



outfile.highest <- paste0(out.dir.ld, "/plink_results_SNPs-highest-LD-TE_R2_chr", chr, ".ld")
fwrite(LD_results_max_R2_nodups, outfile.highest, sep = "\t", quote = FALSE, row.names = FALSE, na = NA)

outfile.lowest <- paste0(out.dir.ld, "/plink_results_SNPs-lowest-LD-TE_R2_chr", chr, ".ld")
fwrite(LD_results_min_R2_nodups, outfile.lowest, sep = "\t", quote = FALSE, row.names = FALSE, na = NA)

outfile.R2 <- paste0(out.dir.ld, "/plink_results_SNPs_all_R2_chr", chr, ".tsv")
fwrite(LD_results_R2, outfile.R2, sep = "\t", quote = FALSE, row.names = FALSE, na = NA)
