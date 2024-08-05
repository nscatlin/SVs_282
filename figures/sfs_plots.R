library(ggplot2)
library(MetBrewer)
library(tidyr)
library(tidyverse)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Input: vector of genotype calls for a given snp (0,1,2,...,ploidy), NAs are encoded as NA
# Output: total number of genotype calls made (number of individual genotypes*ploidy), correcting for missing data
num_geno_calls = function(x, ploidy){
  length(x[!is.na(x)])*ploidy
}

# Input: genotype table with genotype calls varying from 0,1,2,...,ploidy
# Output: vector of allele frequencies for alternate alleles
sfs = function(someTable, ploidy){
  # calculate alternate allele frequency
  alt_counts = rowSums(someTable[,c(-1:-3)], na.rm = T)
  
  # calculate number of genotype calls for each variant
  num_geno_calls = apply(someTable[,c(-1:-3)], MARGIN = 1, FUN = num_geno_calls, ploidy = ploidy)
  
  # calculate alternate allele frequency
  print("Calculating frequency of alternate alleles...")
  alt_freq = alt_counts/num_geno_calls
  result_df = data.frame(key = someTable[, 1], alt_freq = alt_freq)
  
  # Return the resulting dataframe
  return(result_df)
}


df_w_header <- read.csv("../associations/282_Peiffer2014Genetics_w_header_bimbam.csv", header = TRUE)

# df_header_sfs <- sfs(df_w_header, 2)

# Let's look at the SFS for the SVs that are analyzed in GEMMA gwas
svs_in_peiffer_gwas <- read.table("../associations/gwas_SVs.txt", header = FALSE )
# svs_in_peiffer_gwas <- read.table("./redownload_B73_Oh43/svs_in_peiffer_gwas_unimputed.txt", header = FALSE )

df_bimbam_svs_in_peiffer_gwas <- df_w_header[df_w_header$id %in% svs_in_peiffer_gwas$V1, ]

df_bimbam_svs_in_peiffer_gwas_sfs <- sfs(df_bimbam_svs_in_peiffer_gwas, 2)

sv_categories <- read.csv("../data/sv_ids_sv_cateogries.csv", header = TRUE)

df_bimbam_svs_in_peiffer_gwas_sfs_cats <- merge(df_bimbam_svs_in_peiffer_gwas_sfs, sv_categories, by.x = "key", by.y = "sv_id", all.x = TRUE)

df_bimbam_svs_in_peiffer_gwas_sfs_cats <- df_bimbam_svs_in_peiffer_gwas_sfs_cats %>% 
  mutate(category = case_when(
    category == "Category_1"~ "No TE SV",
    category == "Category_2"~ "Incomplete TE SV",
    category == "Category_3"~ "TE = SV",
    category == "Category_4"~ "Multi TE SV",
    category == "Category_5"~ "TE Within SV",
    TRUE~category
  ))

df_bimbam_svs_in_peiffer_gwas_sfs_cats <- df_bimbam_svs_in_peiffer_gwas_sfs_cats %>% 
  mutate(bin = cut(alt_freq, breaks = c(0.1, 0.19, 0.29, 0.39, 0.49, 0.59, 0.69, 0.79, 0.9), labels = FALSE))

for (i in unique(df_bimbam_svs_in_peiffer_gwas_sfs_cats$bin)){
  print(i)
  print(min(df_bimbam_svs_in_peiffer_gwas_sfs_cats$alt_freq[df_bimbam_svs_in_peiffer_gwas_sfs_cats$bin == i]))
  print(max(df_bimbam_svs_in_peiffer_gwas_sfs_cats$alt_freq[df_bimbam_svs_in_peiffer_gwas_sfs_cats$bin == i]))
}



df_bimbam_svs_in_peiffer_gwas_sfs_cats_summary <- df_bimbam_svs_in_peiffer_gwas_sfs_cats %>%
  group_by(bin, category) %>%
  summarise(total_value = sum(alt_freq)) %>%
  ungroup()

df_bimbam_svs_in_peiffer_gwas_sfs_cats_summary <- df_bimbam_svs_in_peiffer_gwas_sfs_cats_summary %>%
  mutate(category = factor(category)) %>%
  bind_rows(df_bimbam_svs_in_peiffer_gwas_sfs_cats_summary %>% 
              group_by(bin) %>% 
              summarise(category = "All SVs Combined", 
                        total_value = sum(total_value)))



df_bimbam_svs_in_peiffer_gwas_sfs_cats_summary_combined <- df_bimbam_svs_in_peiffer_gwas_sfs_cats_summary %>%
  group_by(category) %>%
  mutate(sum_total_value = sum(total_value)) %>%  # Calculate sum of total_value within each category
  ungroup()

df_bimbam_svs_in_peiffer_gwas_sfs_cats_summary_combined <- df_bimbam_svs_in_peiffer_gwas_sfs_cats_summary_combined %>%
  mutate(proportion = total_value / sum_total_value)


df_bimbam_svs_in_peiffer_gwas_sfs_cats_summary_combined_ggplot2 <- drop_na(df_bimbam_svs_in_peiffer_gwas_sfs_cats_summary_combined) %>% 
  ggplot(aes(x = as.factor(bin), y = proportion, fill = category)) +
  geom_bar(stat = "identity", position = "dodge",col=I("black"), linewidth = 0.4) +
  scale_x_discrete(expand = c(0,0), 
                   labels=c("0.1 - 0.19","0.2 - 0.29","0.3 - 0.39","0.4 - 0.49",
                            "0.5 - 0.59","0.6 - 0.69","0.7 - 0.79","0.8 - 0.9"))+
  scale_y_continuous(expand = c(0,0))+
  expand_limits(x = 0, y = 0)+
  scale_fill_manual(name = "SV Category", values = as.vector(met.brewer("Java", n=length(unique(df_bimbam_svs_in_peiffer_gwas_sfs_cats_summary_combined$category)))))+
  labs(x ="Allele Frequency", y = "Proportion of SVs ")+
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.background = element_rect(fill = "white"),  # Set plot background color to white
        panel.background = element_rect(fill = "white"),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.ticks = element_line(size = 0.8),  # Customize tick marks
        axis.ticks.length = unit(0.2, "cm"),
        axis.text.x = element_text(color="black", angle=45, vjust = 0.6),
        axis.text.y = element_text(color="black"),
        axis.text = element_text(size=15), 
        axis.title = element_text(size =17),
        legend.title = element_text(size = 17),
        legend.text = element_text(size = 15),
        plot.margin = margin(t = 1,r = 1, b = 1, l = 1, "cm"))
df_bimbam_svs_in_peiffer_gwas_sfs_cats_summary_combined_ggplot2

# png("./sfs_gwas_svs_proportions_ggplot2_unimputed.png", width = 8, height = 6 , res = 400, units = "in")
# df_bimbam_svs_in_peiffer_gwas_sfs_cats_summary_combined_ggplot2
# dev.off()



