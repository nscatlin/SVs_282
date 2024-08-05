library(ggplot2)
library(MetBrewer)
library(tidyr)
library(tidyverse)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


df_w_header <- read.csv("../associations/282_Peiffer2014Genetics_w_header_bimbam.csv", header = TRUE)

svs_in_peiffer_gwas <- read.table("../associations/gwas_SVs.txt", header = FALSE )
# svs_in_peiffer_gwas <- read.table("./redownload_B73_Oh43/svs_in_peiffer_gwas_unimputed.txt", header = FALSE )

# df_bimbam_svs_in_peiffer_gwas <- df_w_header[df_w_header$id %in% svs_in_peiffer_gwas$V1, ]
# 
# df_bimbam_svs_in_peiffer_gwas_sfs <- sfs(df_bimbam_svs_in_peiffer_gwas, 2)

sv_categories <- read.csv("../data/sv_ids_sv_cateogries.csv", header = TRUE)

sv_categories <- sv_categories %>% 
  mutate(category = case_when(
  category == "Category_1"~ "No TE SV",
  category == "Category_2"~ "Incomplete TE SV",
  category == "Category_3"~ "TE = SV",
  category == "Category_4"~ "Multi TE SV",
  category == "Category_5"~ "TE Within SV",
  TRUE~category
))

df_w_header_w_sv_cats <- merge(df_w_header, sv_categories, by.x = "id", by.y = "sv_id", all.x = TRUE)


te_equal_sv <- df_w_header_w_sv_cats$id[sv_categories$category == 'TE = SV']
b73_oh43_intersect_tes <- read.table("../data/B73_Oh43_bedtools_intersect_TEs.bed", header = FALSE)

df_te_equal_sv <- df_w_header %>% 
  filter(id %in% te_equal_sv)

df_te_equal_sv_w_intersect_tes <- merge(df_te_equal_sv, b73_oh43_intersect_tes[, c("V4", "V8")], by.x = "id", by.y = "V8", all.x = TRUE)

df_te_equal_sv_w_intersect_tes <- df_te_equal_sv_w_intersect_tes[!is.na(df_te_equal_sv_w_intersect_tes$V4), ]
df_te_equal_sv_w_intersect_tes_sfs <- sfs(df_te_equal_sv_w_intersect_tes[,-ncol(df_te_equal_sv_w_intersect_tes)], 2)
df_te_equal_sv_w_intersect_tes_sfs <- merge(df_te_equal_sv_w_intersect_tes_sfs,b73_oh43_intersect_tes[, c("V4", "V8")], by.x = "key", by.y = "V8", all.x = TRUE)


# SFSs
df_te_equal_sv_w_intersect_tes_sfs <- df_te_equal_sv_w_intersect_tes_sfs %>% 
  rename(TE_superfamily = V4)


df_w_header_sfs <- sfs(df_w_header, 2)


# TE = SV SFS Plot (all TE=SV based on Manisha's analysis, regardless if they intersect TE in bedtools) --> 1247 TEs couldn't be found through bedtools intersect
df_te_equal_sv_sfs <- sfs(df_te_equal_sv, 2)




######################
## Figure for paper ##
######################
# png("./sfs_te_equal_sv.png", width = 8, height = 6 , res = 400, units = "in")
ggplot(df_te_equal_sv_sfs, aes(x=alt_freq))+
  geom_histogram(col=I("black"), fill = "grey", linewidth = 0.4)+
  ylim(c(0,3000))+
  scale_x_continuous(breaks = seq(0, 1, by = 0.1), expand = c(0,0))+
  scale_y_continuous(breaks = seq(0, 3000, by = 500), expand = c(0,0))+#seq(0, 3500, by = 500), expand = c(0,0))+
  xlab("Allele Frequency")+
  ylab("Number of SVs ")+
  # ggtitle("Site-frequency Spectrum of SV-present alleles")+
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
        plot.margin = margin(t = 1,r = 1, b = 1, l = 1, "cm"))
# dev.off()






## Change the TE Superfamily strings to something nicer for images
# TE_superfamily_strings<- list(
#   "Gypsy_LTR_retrotransposon" = "Ty3/Gypsy",
#   "Tc1_Mariner_TIR_transposon" = "Tc1/Mariner",
#   "hAT_TIR_transposon" = "hAT",
#   "RTE_LINE_retrotransposon" = "RTE",
#   "PIF_Harbinger_TIR_transposon" = "Pif/Harbinger",
#   "CACTA_TIR_transposon" = "CACTA",
#   "Copia_LTR_retrotransposon" = "Ty1/Copia",
#   "LINE_element" = "Unknown LINE",
#   "L1_LINE_retrotransposon" = "L1 LINE",
#   "Mutator_TIR_transposon" = "Mutator",
#   "helitron " = "Helitron",
#   "LTR_retrotransposon" = "Unknown LTR"
# )

df_te_equal_sv_w_intersect_tes_sfs <- df_te_equal_sv_w_intersect_tes_sfs %>%
  mutate(TE_superfamily = case_when(
    TE_superfamily ==   "Gypsy_LTR_retrotransposon" ~ "Ty3/Gypsy",
    TE_superfamily == "Tc1_Mariner_TIR_transposon" ~ "Tc1/Mariner",
    TE_superfamily == "hAT_TIR_transposon" ~ "hAT",
    TE_superfamily == "RTE_LINE_retrotransposon" ~ "RTE",
    TE_superfamily == "PIF_Harbinger_TIR_transposon" ~ "Pif/Harbinger",
    TE_superfamily == "CACTA_TIR_transposon" ~ "CACTA",
    TE_superfamily == "Copia_LTR_retrotransposon" ~ "Ty1/Copia",
    TE_superfamily == "LINE_element" ~ "Unknown LINE",
    TE_superfamily == "L1_LINE_retrotransposon" ~ "L1 LINE",
    TE_superfamily == "Mutator_TIR_transposon" ~ "Mutator",
    TE_superfamily == "helitron " ~ "Helitron",
    TE_superfamily == "LTR_retrotransposon" ~ "Unknown LTR",
    TRUE~TE_superfamily
  ))


df_te_equal_sv_w_intersect_tes_sfs_btw_0.25_0.75 <- df_te_equal_sv_w_intersect_tes_sfs %>% 
  filter(alt_freq >= 0.25) %>% 
  filter(alt_freq <= 0.75)
df_te_equal_sv_w_intersect_tes_sfs_less_0.25 <- df_te_equal_sv_w_intersect_tes_sfs %>% 
  filter(alt_freq < 0.25)
df_te_equal_sv_w_intersect_tes_sfs_more_0.75 <- df_te_equal_sv_w_intersect_tes_sfs %>% 
  filter(alt_freq > 0.75)



# min_value <- min(df_te_equal_sv_w_intersect_tes_sfs$alt_freq)
# max_value <- max(df_te_equal_sv_w_intersect_tes_sfs$alt_freq)

custom_breaks <- seq(0, 1, by = 0.1)

sorted_categories <- df_te_equal_sv_w_intersect_tes_sfs %>%
  arrange(desc(alt_freq)) %>%
  pull(TE_superfamily)

java_colors <-  c("Ty3/Gypsy" = "#663171",
                  "Tc1/Mariner" = "#8C345B",
                  "Unknown LTR" = "#B23746",
                  "hAT" = "#D13F34",
                  "RTE" = "#DB542F",
                  "Pif/Harbinger" = "#E5692A",
                  "CACTA" = "#E87A39",
                  "Ty1/Copia" = "#E5885D",
                  "Unknown LINE" = "#E29581",
                  "L1 LINE" = "#A78E7B",
                  "Mutator" = "#597F68",
                  "helitron" = "#0C7156")

# I also like the following color palettes from MetBrewer:
#  - Johnson
#  - OKeeffe1



# Ingres= list(c("#041d2c", "#06314e", "#18527e", "#2e77ab", "#d1b252", "#a97f2f", "#7e5522", "#472c0b"), c(4, 5, 3, 6, 2, 7, 1, 8), colorblind=TRUE)
# Java = list(c("#663171", "#cf3a36", "#ea7428", "#e2998a", "#0c7156"), c(1, 4, 2, 5, 3), colorblind=TRUE)
# Johnson = list(c("#a00e00", "#d04e00", "#f6c200", "#0086a8", "#132b69"), c(3, 1, 4, 2, 5), colorblind=TRUE)
# OKeeffe1 = list(c("#6b200c", "#973d21", "#da6c42", "#ee956a", "#fbc2a9", "#f6f2ee", "#bad6f9", "#7db0ea", "#447fdd", "#225bb2", "#133e7e"), c(8, 6, 1, 4, 10, 3, 11, 5, 2, 7, 9), colorblind=TRUE)
# Veronese = list(c("#67322e", "#99610a", "#c38f16", "#6e948c", "#2c6b67", "#175449", "#122c43"), c(5, 1, 7, 2, 3, 6, 4), colorblind=TRUE)





# Stacked plot of all TE=SV TE superfamlies
#######################
## Figures for paper ##
#######################
te_equal_sv_stacked_sfs <- ggplot(df_te_equal_sv_w_intersect_tes_sfs, aes(alt_freq, fill=TE_superfamily))+
  geom_histogram(col=I("black"), linewidth = 0.4 ,breaks = seq(0, 1, by = 0.05))+
  # geom_histogram(col=I("black"), fill = "grey", linewidth = 0.4)+
  scale_x_continuous(breaks = seq(0, 1, by = 0.1), expand = c(0, 0))+
  scale_y_continuous(breaks = seq(0, 2000, by = 250), expand = c(0, 0))+
  # xlab("Allele Frequency")+
  # ylab("Number of TEs ")+
  labs(fill = "TE Superfamily", x="Allele Frequency", y="Number of TEs")+
  scale_fill_manual(values = java_colors)+
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
        plot.margin = margin(t = 1,r = 1, b = 1, l = 1, "cm"))

# png('./redownload_B73_Oh43/te_equal_sv_stacked_sfs.png', width = 8, height = 6 , res = 400, units = "in")
te_equal_sv_stacked_sfs
# dev.off()


sfs_te_superfamily <- function(df,title){ # returns ggplot geom_bar plot
  category_counts <- table(df$TE_superfamily)
  
  # Order categories by frequency counts (lowest to highest)
  ordered_categories <- names(sort(category_counts))
  df$TE_superfamily <- factor(df$TE_superfamily, levels = ordered_categories)
  plot <-ggplot(df, aes(TE_superfamily, fill=TE_superfamily)) +
    geom_bar(col=I("black"), linewidth = 0.4, show.legend = FALSE) +
    labs(fill = "TE Superfamily", x="TE Superfamily", y="Number of TEs",  title=title)+
    scale_fill_manual(values = java_colors)+
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          plot.background = element_rect(fill = "white"),  # Set plot background color to white
          panel.background = element_rect(fill = "white"),
          panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          axis.ticks = element_line(size = 0.5),  # Customize tick marks
          axis.ticks.length = unit(0.2, "cm"),
          axis.text.x = element_text(color="black",angle = 45, hjust = 1),
          axis.text.y = element_text(color="black"),
          axis.text = element_text(size=12), 
          axis.title = element_text(size =14),
          legend.title = element_text(size = 14),
          legend.text = element_text(size = 12),
          plot.margin = margin(t = 1,r = 1, b = 1, l = 1, "cm"))
  return(plot)
  
}

pcat2_func <-sfs_te_superfamily(df_te_equal_sv_w_intersect_tes_sfs_btw_0.25_0.75,"0.25 \u2264 MAF \u2264 0.75 for TE = SV")
pcat3_func <- sfs_te_superfamily(df_te_equal_sv_w_intersect_tes_sfs_less_0.25, "MAF < 0.25 TE=SV")
pcat4_func <- sfs_te_superfamily(df_te_equal_sv_w_intersect_tes_sfs_more_0.75, "MAF \u2265 0.75 TE = SV")

######################
## Figure for paper #
######################
grouped_te_superfamily_plots <- ggarrange(te_equal_sv_stacked_sfs+labs(tag = 'A.', title="SFS for TE = SV by TE Superfamily")+theme(legend.key.size = unit(0.5, "cm"), legend.title = element_text(size = 11,),legend.text = element_text(size = 10)),pcat2_func+labs(tag = 'B.'), pcat3_func+labs(tag = 'C.'), pcat4_func+labs(tag = 'D.'))
# png('./redownload_B73_Oh43/te_equal_sv_stacked_sfs_w_sfs_0.25_0.75.png', width = 15, height = 10 , res = 400, units = "in")
grouped_te_superfamily_plots
# dev.off()

