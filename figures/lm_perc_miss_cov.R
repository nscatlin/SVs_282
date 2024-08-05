library(ggplot2)
library(tidyr)
library(tidyverse)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#####################################################
### Linear Model of perc. missing data ~ coverage ###
#####################################################

coverage_stats <- read.csv("../data/coverages_alignments_to_B73_V5.csv", header = FALSE)
df_w_header <- read.csv("../associations/282_Peiffer2014Genetics_w_header_bimbam.csv", header = TRUE)
names(coverage_stats) <- c("genotype","avg_coverage")
coverage_stats$genotype <- gsub('(^[0-9])', 'X\\1', coverage_stats$genotype)
coverage_stats$genotype <-gsub('-', '.',coverage_stats$genotype)
coverage_stats$genotype <- tolower(coverage_stats$genotype)

missing_data <- vector("list", length = ncol(df_w_header[-1:-3]))
df_subset <- df_w_header[, -(1:3)]

# Loop through each column
for (i in seq_along(df_subset)) {
  # Calculate percentage of missing data for the current column
  missing_percentage <- mean(is.na(df_subset[[i]])) * 100
  
  # Store column name and percentage of missing data
  missing_data[[i]] <- c(colnames(df_subset)[i], missing_percentage)
}

# Convert the list to a dataframe
missing_data_df <- do.call(rbind.data.frame, missing_data)

# Rename columns
colnames(missing_data_df) <- c("Column Name", "Percentage of Missing Data")

merged_sv_miss_coverage <- merge(missing_data_df, coverage_stats, by.x="Column Name", by.y = "genotype", all = TRUE)

# Fill missing values with NA
merged_sv_miss_coverage[is.na(merged_sv_miss_coverage)] <- NA



colnames(merged_sv_miss_coverage) <- c("genotypes", "perc_missing_data", "avg_coverage")

merged_sv_miss_coverage$perc_missing_data <- as.numeric(as.character(merged_sv_miss_coverage$perc_missing_data))
# plot(merged_sv_miss_coverage$Percentage_of_Missing_Data~merged_sv_miss_coverage$avg_coverage)
# Since B73 and Oh43 do not have nay missing data, let's remove those from the dataset

merged_sv_miss_coverage_wo_B73_Oh43 <- merged_sv_miss_coverage %>% 
  filter(genotypes != 'Oh43' & genotypes != 'B73' )
# plot(model)
# write.table(merged_sv_miss_coverage_wo_B73_Oh43, "./282_genotypes_perc_miss_avg_cov.csv", sep = ',', row.names=FALSE, col.names=TRUE, quote=FALSE)



lm_model <- lm(perc_missing_data~avg_coverage, data = merged_sv_miss_coverage_wo_B73_Oh43)





# png("./redownload_B73_Oh43/lm_perc_missing_data_coverage_all_genos_wo_B73_Oh43.png", width = 8, height = 6 , res = 400, units = "in")
ggplot(data = merged_sv_miss_coverage_wo_B73_Oh43, aes(x = avg_coverage, y = as.numeric(perc_missing_data)))+
  geom_point()+
  geom_smooth(method = "lm")+
  scale_x_continuous(limits = c(0,10), expand = c(0,0))+
  scale_y_continuous(limits = c(0,100),expand = c(0,0))+#seq(0, 3500, by = 500), expand = c(0,0))+
  labs(x ="Realized Coverage", y = "Percent Missing Data ")+
  ylim(0,100)+
  xlim(0,max(merged_sv_miss_coverage_wo_B73_Oh43$avg_coverage) +5)+
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
# + annotate("text", x = 22, y = 75, label =expression(Adj.~R^2 == 0.067), parse = TRUE, size = 5)
# dev.off()
summary(lm_model)



# Let's look for outliers
plot(lm_model)

# Looks like row 265 is an outlier (vaw6)
wo_vaw6 <- merged_sv_miss_coverage_wo_B73_Oh43 %>% 
  slice(-265)

model_wo_vaw6 <- lm(perc_missing_data~avg_coverage, data = wo_vaw6)
plot(model_wo_vaw6)

# Looks like row 31 is an outlier (b57)
wo_vaw6_b57 <- wo_vaw6 %>% 
  slice(-31)

model_wo_vaw6_b57 <- lm(perc_missing_data~avg_coverage, data = wo_vaw6_b57)

ggplot(data = wo_vaw6_b57, aes(x = avg_coverage, y = as.numeric(perc_missing_data)))+
  geom_point()+
  geom_smooth(method = "lm")+
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))+#seq(0, 3500, by = 500), expand = c(0,0))+
  labs(x ="Realized Coverage", y = "Percent Missing Data ")+
  ylim(0,100)+
  xlim(0,10)+
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


