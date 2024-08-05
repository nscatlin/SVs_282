library(ggplot2)

# Freq. of missing data per SV-present allele

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

bimbam_ordered_peiffer_intersection <- read.table("../associations/282_Peiffer2014Genetics_w_header_bimbam.csv", sep = ',', header = TRUE)

na_percentage <- rowSums(is.na(bimbam_ordered_peiffer_intersection[,c(-1:-3)])) / ncol(bimbam_ordered_peiffer_intersection[,c(-1:-3)]) * 100

na_df <- data.frame(na_percentage)

non_na_percentage <- rowSums(!is.na(bimbam_ordered_peiffer_intersection[,c(-1:-3)])) / ncol(bimbam_ordered_peiffer_intersection[,c(-1:-3)]) * 100

non_na_df <- data.frame(non_na_percentage)

# png("./perc_genotype_called_per_sv.png", width = 10, height = 8 , res = 400, units = "in")
ggplot(non_na_df, aes(x = non_na_percentage)) +
  geom_histogram(binwidth = 2,aes(fill = non_na_df > 90), boundary = 0,col=I("black"), linewidth = 0.4,show.legend = FALSE) +
  scale_x_continuous(breaks = seq(0, 100, by = 10), limits = c(0, 100),expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))+#seq(0, 3500, by = 500), expand = c(0,0))+
  # scale_y_discrete(expand = c(0, 0)) +
  scale_fill_manual(values = c("FALSE" = "grey", "TRUE" = "green")) +
  labs(x = "Percentage of Genotypes Called",
       y = "Frequency") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.background = element_rect(fill = "white"),  # Set plot background color to white
        panel.background = element_rect(fill = "white"),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.ticks = element_line(size = 0.8),  # Customize tick marks
        axis.ticks.length = unit(0.2, "cm"),
        axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black"),
        axis.text = element_text(size=15), 
        axis.title = element_text(size =17),
        legend.title = element_text(size = 17),
        legend.text = element_text(size = 15),
        plot.margin = margin(t = 1,r = 1, b = 1, l = 1, "cm"))
# dev.off()

length(missing_counts[missing_counts <=0.1])
length(missing_counts[missing_counts <=0.2])
length(missing_counts[missing_counts <=0.3])
length(missing_counts[missing_counts <=0.4])