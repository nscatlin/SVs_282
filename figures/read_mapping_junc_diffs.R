library(ggplot2)
library(ggpubr)
library(grid)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

Oh43_query_df <- read.csv("../data/Oh43_reads_to_pseudoref_alleles.csv", header = TRUE)
B73_query_df <- read.csv("../data/B73_reads_to_pseudoref_alleles.csv", header = TRUE)
master_B73_reads_Oh43_reads_to_pseudoref <- rbind(B73_query_df,Oh43_query_df)

# Left
master_B73_reads_Oh43_reads_to_pseudoref_left_juncs_stack <- master_B73_reads_Oh43_reads_to_pseudoref %>%
  ggplot(aes(x = as.numeric(ref_minus_query_diff_left))) +
  geom_bar(aes(fill = sv_present_genotype), position = "stack", width=1.0)+
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))+
  # scale_y_break(c(15000, 100000))+
  # color_palette(pal)+
  # scale_fill_manual(values = josephs_palette("josephs_true"))+
  xlim(-25,40)+
  labs(x = "Difference in reads", y = "Number of junctions", title = "Left junctions")+
  # xlab("Difference in reads")+
  # ylab("Number of junctions")+
  geom_vline(xintercept = 0, linetype="dotted", 
             color = "red", linewidth=0.75)+
  # ggtitle("Left junctions")+
  # scale_fill_discrete(name = "Reference Genotype")++
  scale_fill_manual(values=josephs_palette(n=2, name="josephs_colorblind"), name = "Reference Genotype")+
  # theme_bw()+
  # theme(text = element_text(family = "Arial"),
  #       axis.text.x.top = element_blank(),
  #       axis.ticks.x.top = element_blank(),
  #       plot.title = element_text(hjust = 0.5, face = "bold"),
  #       axis.text = element_text(size=12), 
  #       axis.title = element_text(size =16),
  #       legend.title = element_text(size = 14),
  #       legend.text = element_text(size = 12))
  #
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
        plot.margin = margin(1,1,1.5,1.2, "cm"))+
  scale_fill_grey(name = "Reference Genotype")

# Right
master_B73_reads_Oh43_reads_to_pseudoref_right_juncs_stack <- master_B73_reads_Oh43_reads_to_pseudoref %>%
  ggplot(aes(x = as.numeric(ref_minus_query_diff_right))) +
  geom_bar(aes(fill = sv_present_genotype), position = "stack", width=1.0)+
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))+
  # scale_y_break(c(15000, 100000))+
  # color_palette(pal)+
  # scale_fill_manual(values = josephs_palette("josephs_true"))+
  xlim(-25,40)+
  labs(x = "Difference in reads", y = "Number of junctions", title = "Right junctions")+
  # xlab("Difference in reads")+
  # ylab("Number of junctions")+
  geom_vline(xintercept = 0, linetype="dotted", 
             color = "red", linewidth=0.75)+
  # ggtitle("Right junctions")+
  # scale_fill_discrete(name = "Reference Genotype")++
  scale_fill_manual(values=josephs_palette(n=2, name="josephs_colorblind"), name = "Reference Genotype")+
  # theme_bw()+
  # theme(text = element_text(family = "Arial"),
  #       axis.text.x.top = element_blank(),
  #       axis.ticks.x.top = element_blank(),
  #       plot.title = element_text(hjust = 0.5, face = "bold"),
  #       axis.text = element_text(size=12), 
  #       axis.title = element_text(size =16),
  #       legend.title = element_text(size = 14),
  #       legend.text = element_text(size = 12))
  #
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
        plot.margin = margin(1,1,1.5,1.2, "cm"))+
  scale_fill_grey(name = "Reference Genotype")

fig_master_left_right <- ggarrange(master_B73_reads_Oh43_reads_to_pseudoref_left_juncs_stack + rremove("ylab") + rremove("xlab"),
                                   master_B73_reads_Oh43_reads_to_pseudoref_right_juncs_stack + rremove("ylab") + rremove("xlab"),
                                   nrow=1, ncol=2, common.legend = TRUE, legend = "right")
annotate_figure(fig_master_left_right, left = textGrob("Number of junctions", rot = 90, vjust = 1, gp = gpar(cex = 1.3)),
                bottom = textGrob("Difference in read mapping", gp = gpar(cex = 1.3), hjust = 0.8))

