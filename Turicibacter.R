setwd("~/Documents/Model-AD/5xFAD Manuscript/20220801_mBio_R1/github_R1")

#clear objects from the workspace
rm(list=ls())

library(tidyverse)
library(dplyr)
library(tibble) 
library(ggplot2)
library(ggpubr)

#Load and format data
turicibacter_df_meta <-read.table("./Turicibacter_data.txt", 
                      header = TRUE, sep = '\t', quote = "", row.names = 1)

#Order
turicibacter_df_meta$Cohort = factor(turicibacter_df_meta$Cohort, levels = c('4mo','8mo','12mo','18mo'))

###Plot
ggplot(data = turicibacter_df_meta) +
  aes(x = Var2, y = (value*100), fill = Genotype) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values=c('#d94801', '#2171b5')) +
  geom_point(size = 1.5, colour = "black", alpha = 0.5, 
             position = position_jitterdodge(jitter.width = 0.1),
             aes(fill=Genotype, shape=Sex)) +   
  facet_wrap(Type+Cohort~., nrow = 1) +
  scale_y_log10() +
  theme_bw() +
  theme(title =  element_text(size = 10,face = "bold"))+
  theme(legend.text =  element_text(hjust = 0, size = 10,face = "bold"))+
  theme(legend.title = element_text(hjust = 0, size = 10,face = "bold"))+ 
  theme(axis.text.y = element_text(hjust = 0, size = 8, face = "bold"))+
  theme(axis.text.x = element_text(hjust = 1, size = 8, face = "bold", angle = 45))+
  theme(axis.title.y = element_text(size = 10, face = "bold"))+
  theme(axis.title.x = element_text(size = 10, face = "bold"))+
  labs(x = NULL, y = "Relative Abundance (Log10 x 100)", fill = "Genotype", 
       title = "")+
  stat_compare_means(method = "wilcox.test", label = "p.signif",
                     label.x = 1.5, vjust = 0.5, size = 4, hide.ns = T) 

#ggsave("./TRB.pdf", height = 4, width = 6.25)

