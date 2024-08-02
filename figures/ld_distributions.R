library(ggplot2)
library(ggpubr)
library(tidyverse)
library(rstudioapi)
library(ggforce)


setwd(dirname(getActiveDocumentContext()$path))
df <- read.table("../LD/plink_results_SNPs-highest-LD-TE_R2_allchrs.ld", header = TRUE)

sig_hapmap_1 <- read.table("../LD/hapmap321/1_hapmap321_peiffer_intersection_plink_all_lmm_options_maf0.1.assoc.txt_sig_results.txt", header = TRUE)
sig_hapmap_1 <- sig_hapmap_1 %>% 
  rename(snp = rs)
sig_hapmap_2 <- read.table("../LD/hapmap321/2_hapmap321_peiffer_intersection_plink_all_lmm_options_maf0.1.assoc.txt_sig_results.txt", header = TRUE)
sig_hapmap_2 <- sig_hapmap_2 %>% 
  rename(snp = rs)

sv_gwas_bimbam <- read.table("../associations/282_Peiffer2014Genetics_w_header_bimbam.csv", header=TRUE)

df_r2_0.1 <- df[df$R2>0.1,]
df_r2_0.5 <- df[df$R2>0.5,]
df_r2_1 <- df[df$R2==1,]
df %>% 
  filter(R2 == 1) %>% 
  nrow()

# perc_0.1 <- (nrow(df_r2_0.1)/nrow(df) * 100)
# perc_0.5 <- (nrow(df_r2_0.5)/nrow(df) * 100)
# 
# perc_ld <- nrow(df)/nrow(sv_gwas_bimbam)*100
# 
# nrow(df) - nrow(df_r2_0.5)
# write.table(df_r2_0.1, "./ld_SVs_hapmapSNPS_1Mb_r2_greater_than_0.1.csv", sep = ',', row.names=FALSE, col.names=TRUE, quote=FALSE)

# write.table(df_r2_0.5, "./ld_SVs_hapmapSNPS_1Mb_r2_greater_than_0.5.csv", sep = ',', row.names=FALSE, col.names=TRUE, quote=FALSE)

sig_snps_in_ld_1 <- df_r2_0.1[df_r2_0.1$SNP_B %in% sig_hapmap_1$snp, ]
sig_snps_in_ld_2 <- df_r2_0.1[df_r2_0.1$SNP_B %in% sig_hapmap_2$snp, ]



# Significant SVs
sig_svs <- c("B73_ins_chr7:181347609-181348458_Oh43:175922047_29423",
             "B73_ins_chr10:31349992-31394134_Oh43:31339343_5713",
             "Oh43_ins_chr3:90730550-90783144_B73:89951120_16536",
             "Oh43_ins_chr3:90730550-90783144_B73:89951120_16536", 
             "Oh43_ins_chr4:234896964-234898413_B73:230228965_32318")


sig_sv_in_ld_w_sig_snp1 <- sig_snps_in_ld_1[sig_snps_in_ld_1$SNP_A %in% sig_svs,] 
sig_sv_in_ld_w_sig_snp2 <- sig_snps_in_ld_2[sig_snps_in_ld_2$SNP_A %in% sig_svs,] 

df %>% 
  filter(SNP_A %in% sig_svs) %>% 
  select(SNP_A, R2)


plot_ld <- function(df, variable, bins = 30, title = ''){
  ggplot(df, aes(x={{variable}}))+
    geom_histogram(col=I("black"), fill = "grey", linewidth = 0.4)+
    ylab("Number of SVs with SNP within 1Mb")+
    xlab(expression(italic(r)^2))+
    ggtitle(title)+
    scale_x_continuous(breaks = seq(0, 1, by = 0.1), expand = c(0,0))+
    scale_y_continuous(expand = c(0,0))+#seq(0, 3500, by = 500), expand = c(0,0))+
    # scale_fill_discrete(name = "Reference Genotype")+
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          plot.background = element_rect(fill = "white"),  # Set plot background color to white
          panel.background = element_rect(fill = "white"),
          panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          axis.ticks = element_line(size = 0.5),  # Customize tick marks
          axis.ticks.length = unit(0.2, "cm"),
          axis.text.x = element_text(color="black"),
          axis.text.y = element_text(color="black"),
          axis.text = element_text(size=12), 
          axis.title = element_text(size =14),
          legend.title = element_text(size = 14),
          legend.text = element_text(size = 12),
          plot.margin = margin(1,1,1.5,1.2, "cm"))
}

plot_ld_0.1 <- plot_ld(df=df_r2_0.1, variable=R2, title=expression("SVs with "~r^2~" greater than 0.1 for SNPs within 1Mb"))
plot_ld_0.5 <- plot_ld(df=df_r2_0.5, variable=R2, title=expression("SVs with "~r^2~" greater than 0.5 fo r SNPs within 1Mb"))
plot_highest_ld <- plot_ld(df=df, variable=R2)


png("./highest_ld_between_SVs_SNPs_within_1Mb.png", width = 8, height = 6 , res = 400, units = "in")
plot_highest_ld
dev.off()

# If you are interested in looking at a zoomed in view
library(ggbreak)
plot_highest_ld_zoomed <- plot_highest_ld+
  facet_zoom(ylim = c(0, 35), 
             zoom.size = 0.5)
png("./highest_ld_between_SVs_SNPs_within_1Mb_zoomed_ymax35.png", width = 8, height = 6 , res = 400, units = "in")
plot_highest_ld_zoomed
dev.off()





# TE = SV (cateogry 5)


sv_categories <- read.csv("../data/sv_ids_sv_cateogries.csv")
sv_categories <- sv_categories %>% 
  mutate(category = case_when(
    category == "Category_1"~ "No TE SV",
    category == "Category_2"~ "Incomplete TE SV",
    category == "Category_3"~ "TE = SV",
    category == "Category_4"~ "Multi TE SV",
    category == "Category_5"~ "TE Within SV",
    TRUE~category
  ))


ld_w_sv_categories_df <- merge(df, sv_categories[, c("sv_id", "category")], by.x = "SNP_A", by.y = "sv_id", all.x = TRUE)


te_equal_sv <- ld_w_sv_categories_df %>% 
  filter(te_category == "Category_5")

te_equal_sv_ld <- plot_ld(df=te_equal_sv, variable=R2)

png("./TE_equal_SV_highest_ld_between_SVs_SNPs_within_1Mb.png", width = 8, height = 6 , res = 400, units = "in")
te_equal_sv_ld
dev.off()


plot_te_equal_sv_ld_zoomed <- te_equal_sv_ld+
  facet_zoom(ylim = c(0, 15), 
             zoom.size = 0.5)

png("./TE_equal_SV_highest_ld_between_SVs_SNPs_within_1Mb_zoomed_ymax15.png", width = 8, height = 6 , res = 400, units = "in")
plot_te_equal_sv_ld_zoomed
dev.off()


