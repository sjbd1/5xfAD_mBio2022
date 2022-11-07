#Unload all loaded packages

(.packages())
invisible(lapply(paste0("package:", names(sessionInfo()$otherPkgs)),   # Unload add-on packages
                 detach,
                 character.only = TRUE, unload = TRUE))
(.packages())

#Assign the working directory
setwd("~/Documents/Model-AD/5xFAD Manuscript/20220801_mBio_R1/github_R1")


#Note: each section of the script is designed to work independently and has 
#variable names in common with other sections. Safest to "clear objects from the workspace" 
#with rm(list=ls()) before runnning each section.

##################### Heatmap of 100 Most Abundant Species #######################
#Adapted from Dave Tang's pheatmap example: https://davetang.org/muse/2018/05/15/making-a-heatmap-in-r-with-the-pheatmap-package/

#clear objects from the workspace
rm(list=ls())

# load package
library(dplyr)
library(tibble)
library(pheatmap)
library(RColorBrewer)

#Load and format data
data_meta <- read.table("./Taxonomy_data.txt", 
                      header = TRUE, sep = '\t', quote = "", row.names = 1)

#filter for: Cecal/Fecal samples only
data_meta <- filter(data_meta, Type == "CECAL")

#Change "HET" to "5xfAD"
data_meta$Genotype <- replace(data_meta$Genotype, data_meta$Genotype == "HET", "5xfAD")

#Order data by "Sex"
data_meta <- data_meta[order(factor(data_meta$Sex, 
                                    levels = c("Female", "Male"))),]
#Order data by "Genotype"
data_meta <- data_meta[order(factor(data_meta$Genotype, 
                                    levels = c("WT","5xfAD"))),]
#Order data by "Cohort"
data_meta <- data_meta[order(factor(data_meta$Cohort, 
                                    levels = c("18mo", "12mo", "8mo","4mo")), decreasing = TRUE),]
#Order data by "Type"
data_meta <- data_meta[order(factor(data_meta$Type, 
                                    levels = c("CECAL", "FECAL"))),]
#select for metadata
meta <- data_meta[,1:8]

###copy rownames to a column and add to new "group" column with other metadata
meta <- rownames_to_column(meta, var = "Sample_ID")

#Define group (e.g., Cohort_Line_SampleType)
meta$Group <- paste(meta$Sample_ID, meta$Genotype, sep="_")

#Select for data
data <- data_meta[,-(1:8)]

##selecting most abundant (e.g., Top10)
n_otu = 100
top_otus <- as.data.frame(colSums(data)) %>%
  rownames_to_column() %>%
  top_n(n_otu)

top <- data[colnames(data) %in% top_otus$rowname]

#Join with metadata and select for data & metadata subset
top_meta <- cbind(meta, top)
top_meta <- rownames_to_column(top_meta, var = "rownames")
top_meta <- column_to_rownames(top_meta, var = "Group")

top_group <- (top_meta[,-(1:10)])

#transpose for better plotting
top_trans <- t(top_group)

#test plot
pheatmap(top_trans)

#Standardized normal distribution
normfunction <- function(x){
  (x - mean(x)) / sd(x)
}

top_norm <- as.data.frame(t(apply(top_trans, 1, normfunction)))

#test plot
pheatmap(top_norm)

#Design a key for column annotation
key <- as.data.frame((cbind( top_meta$Sex, top_meta$Genotype, top_meta$Cohort, top_meta$Type)), 
                     row.names = row.names(top_meta))
colnames(key) <- c("Sex", "Genotype","Cohort", "Sample Type") #assign colnames
key$Cohort <- factor(key$Cohort, levels = c("4mo","8mo", "12mo", "18mo"))

#https://colorbrewer2.org/#type=diverging&scheme=PuOr&n=10
my_colour = list(
  'Sample Type' = c(CECAL = "#fc8d59", FECAL = "#91bfdb"),
  Cohort = c('4mo' = "#d7191c", '8mo' = "#fdae61", '12mo'="#abd9e9", '18mo'='#2c7bb6'),
  #Cohort = c('18mo'='#2c7bb6'),
  Genotype = c('5xfAD' = "#d94801", 'WT' = "#2171b5"),
  Sex = c(Female = "#80cdc1", Male="#000000" )
)

#define color pallet, bias is adjusted to use maximum amount of range and center white @ ~ z = 0
cp <- colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu")), #name = "RdYlBu" name = "Spectral"
                       bias = 2)(100)

#plot
pheatmap(top_norm, clustering_method = 'average', cluster_rows = TRUE, cluster_cols = FALSE,
         gaps_col = FALSE, show_colnames = FALSE, show_rownames = TRUE, 
         fontsize = 8, fontsize_row = 6, fontsize_col = 6, 
         cutree_cols = 2, cutree_rows = 10, angle_col = 270, 
         annotation_legend = TRUE, annotation_col = key, annotation_colors = my_colour, 
         color = cp,legend = TRUE, legend_breaks = c(-1,2,5,8,11,14), legend_labels = NA, 
         annotation_names_row = TRUE, annotation_names_col = TRUE,
         main = "")

##################### Species Richness #################################
#clear objects from the workspace
rm(list=ls())

#Load packages
library(dplyr)
library(ggplot2)
library(ggpubr)

#Load and format data
data_meta <- read.table("./Taxonomy_data.txt", 
                      header = TRUE, sep = '\t', quote = "", row.names = 1)

## Prepare data
#Extract metadata
meta <- data_meta[,(1:8)]
    
#Change "HET" to "5xfAD"
meta$Genotype <- replace(meta$Genotype, meta$Genotype=="HET", "5xfAD")

#Define group
meta$Group <- paste(meta$Cohort, meta$Genotype, meta$Type, sep="_")
    
#Extract data
OTUs <- data_meta[,-(1:8)]

#calculate richness
spp_richness <- data.frame(apply(OTUs[,-1]>0,1,sum))
names(spp_richness)[1]<-paste("richness")

#Join with metadata
spp_richness_meta <-cbind(meta, spp_richness)
      
#Make metadata factors to enable stats
sapply(spp_richness_meta, is.factor)
spp_richness_meta$Group <- as.factor(spp_richness_meta$Group)
spp_richness_meta$Genotype <- as.factor(spp_richness_meta$Genotype)
spp_richness_meta$Cohort <- as.factor(spp_richness_meta$Cohort)
spp_richness_meta$Sex <- as.factor(spp_richness_meta$Sex)
spp_richness_meta$Type <- as.factor(spp_richness_meta$Type)
      
#Test significance of richness by genotype (seperates ages and sample types)
data_summary <- spp_richness_meta[,-(1:8)] %>%
                group_by(Group) %>%
                summarise(mean = mean(richness), 
                          median = median(richness),
                          se = (sd(richness)) / sqrt(length(richness)))

#Shapiro-Wilk normality test; 
    #if p < 0.05 the data is not normally distributed
    #if p > 0.05 the data is normally distributed
shapiro.test(data_summary$mean)
# W = 0.93799, p-value = 0.4724
# data is normally distributed
      
#anova by genotype, age, and sample type (*assumes a synergistic effect)
anova_test <- aov(richness  ~ Genotype*Cohort*Type, data = spp_richness_meta)
summary(anova_test)
#                       Df Sum Sq Mean Sq F value   Pr(>F)    
#Genotype               1    350     350   0.142   0.7068    
#Cohort                 3  70759   23586   9.567 5.71e-06 ***
#Type                   1 152899  152899  62.016 1.46e-13 ***
#Genotype:Cohort        3    129      43   0.018   0.9968    
#Genotype:Type          1    174     174   0.071   0.7906    
#Cohort:Type            1      8       8   0.003   0.9552    
#Genotype:Cohort:Type   1   7672    7672   3.112   0.0791 .  
#Residuals            224 552269    2465                     
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1             
      
#Faceted Plot
#Arrange metadata
spp_richness_meta <- spp_richness_meta %>%
                    mutate(Genotype =  factor(Genotype, levels = c("5xfAD","WT"))) %>%
                    arrange(Genotype)
      
spp_richness_meta <- spp_richness_meta %>%
                    mutate(Cohort =  factor(Cohort, levels = c("4mo", "8mo", "12mo", "18mo"))) %>%
                    arrange(Cohort)
      
ggplot(data = spp_richness_meta) +
  aes(x = Cohort, y = richness, fill = Genotype) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values=c('#d94801', '#2171b5')) +
  #ylim(0, 3e+08)+
  geom_point(size = 2, colour = "black", alpha = 0.5, 
             position = position_jitterdodge(jitter.width = 0.1),
             aes(fill=Genotype, shape=Sex)) +
  facet_wrap(Type~., nrow = 1, scales = "free_x") +
  theme_bw() +
  theme(axis.text = element_text(size = 10,face = "bold"),
        axis.title = element_text(size = 10,face = "bold"))+
  theme(legend.title = element_text(hjust = 0, size = 10,face = "bold"))+ 
  theme(axis.text.y = element_text(hjust = 0, size = 10, face = "bold"))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "bold"))+
  labs(x = "", y = "Number of Species Detected", fill = "Genotype") +
  ggtitle("")+
  theme(plot.title = element_text(hjust = 0, size = 12,face = "bold")) +
  stat_compare_means(method = "wilcox.test",label = "p.signif",label.x = 1.5,vjust = 0, size = 4,hide.ns = F) ## package ggpubr lets us put significance tests directly on the plot!
      
#ggsave("./Richness.pdf", height = 3.5, width = 4.5)

##################### Shannon Diversity #################################
#clear objects from the workspace
rm(list=ls())

#Load packages
library(dplyr)
library(ggplot2)
library(ggpubr)
library(vegan)

#import data
#Load and format data
data_meta <- read.table("./Taxonomy_data.txt", 
                        header = TRUE, sep = '\t', quote = "", row.names = 1)
      
## Prepare data
#Extract metadata
meta <- data_meta[,(1:8)]
      
#Change "HET" to "5xfAD"
meta$Genotype <- replace(meta$Genotype, meta$Genotype=="HET", "5xfAD")
      
#Define group
meta$Group <- paste(meta$Cohort, meta$Genotype, meta$Type, sep="_")

#Extract data
OTUs <- data_meta[,-(1:8)]

shannon <- diversity(OTUs, "shannon")
is.numeric(shannon)

#Shapiro-Wilk normality test
shapiro.test(shannon)
#if p < 0.05 the data is not normally distributed
#if p > 0.05 the data is normally distributed
# W = 0.95011, p-value = 2.977e-07 --> data is NOT normally distributed
      
shan.alpha <- data.frame(cbind(meta$Group, meta$Genotype,
                               meta$Cohort, meta$Sex, meta$Type, shannon))
colnames(shan.alpha) <- c('Group', 'Genotype', 'Cohort', 'Sex', 'Type', 'shannon')
      
sapply(shan.alpha, is.factor)
shan.alpha$Group <- as.factor(shan.alpha$Group)
shan.alpha$Genotype <- as.factor(shan.alpha$Genotype)
shan.alpha$Cohort <- as.factor(shan.alpha$Cohort)
shan.alpha$Sex <- as.factor(shan.alpha$Sex)
shan.alpha$Type <- as.factor(shan.alpha$Type)
      
is.numeric(shan.alpha$shannon)
shan.alpha$shannon <- as.numeric(shan.alpha$shannon)
      
##Kruskal Wallis rank sum test:
    # p < 0.05 some group medians are significantly different
    # p > 0.05 no group medians are significantly different
kruskal.test(shannon ~ Group, data = shan.alpha)
#Kruskal-Wallis chi-squared = 52.345, df = 11, p-value = 2.361e-07
kruskal.test(shannon ~ Genotype, data = shan.alpha)
#Kruskal-Wallis chi-squared = 0.00023441, df = 1, p-value = 0.9878
kruskal.test(shannon ~ Cohort, data = shan.alpha)
#Kruskal-Wallis chi-squared = 21.042, df = 3, p-value = 0.0001032
kruskal.test(shannon ~ Sex, data = shan.alpha)
#Kruskal-Wallis chi-squared = 2.7859, df = 1, p-value = 0.0951
kruskal.test(shannon ~ Type, data = shan.alpha)
#Kruskal-Wallis chi-squared = 10.33, df = 1, p-value = 0.001309

##Faceted by Sample Type 
shan.alpha <- shan.alpha %>%
  mutate(Cohort =  factor(Cohort, levels = c("4mo", "8mo",  "12mo", "18mo"))) %>%
  arrange(Cohort)   
  
ggplot(data = shan.alpha) +
  aes(x = Cohort, y = shannon, fill = Genotype) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values=c('#d94801', '#2171b5')) +
  #ylim(0, 3e+08)+
  geom_point(size = 2, colour = "black", alpha = 0.5, 
             position = position_jitterdodge(jitter.width = 0.1),
             aes(fill=Genotype, shape=Sex)) +
  facet_wrap(Type~., nrow = 1, scales = "free_x") +
  theme_bw() +
  theme(axis.text = element_text(size = 10,face = "bold"),
        axis.title = element_text(size = 10,face = "bold"))+
  theme(legend.title = element_text(hjust = 0, size = 10,face = "bold"))+ 
  theme(axis.text.y = element_text(hjust = 0, size = 10, face = "bold"))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "bold"))+
  labs(x = "", y = "Shannon Diversity", fill = "Genotype") +
  ggtitle("")+
  theme(plot.title = element_text(hjust = 0, size = 12,face = "bold")) +
  stat_compare_means(method = "wilcox.test",label = "p.signif",label.x = 1.5,vjust = 0, size = 4,hide.ns = T) ## package ggpubr lets us put significance tests directly on the plot!

#ggsave("./Shannon_Type.pdf", height = 3.5, width = 4.5)
      
##Faceted by Sample Type & Sex
ggplot(data = shan.alpha) +
  #aes(x = Genotype, y = value, fill = Var2) +
  aes(x = Cohort, y = shannon, fill = Genotype) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values=c('#d94801', '#2171b5')) +
  geom_point(size = 2, colour = "black", alpha = 0.5, 
             position = position_jitterdodge(jitter.width = 0.1),
             aes(fill=Genotype, shape=Sex)) +
  #geom_point(size = 2, colour = "black", alpha = 0.5, position = position_jitterdodge(jitter.width = 0)) +
  facet_wrap(Type+Sex~., nrow = 1, scales = "free_x") +
  theme(axis.text = element_text(size = 10,face = "bold"),
        axis.title = element_text(size = 10,face = "bold"))+
  theme(legend.title = element_text(hjust = 0, size = 10,face = "bold"))+ 
  theme(axis.text.y = element_text(hjust = 0, size = 10, face = "bold"))+
  theme(axis.text.x = element_text(hjust = 0, size = 10, face = "bold"))+
  theme_bw() +
  labs(x = "Age", y = "Shannon Diversity", fill = "Genotype") +
  stat_compare_means(method = "wilcox.test",label = "p.signif",label.x = 1.5,vjust = 0, size = 4, hide.ns = T) ## package ggpubr lets us put significance tests directly on the plot!
      
#ggsave("./Shannon_Sex&Type.pdf", height = 3, width = 7)
      
##################### PCA ###########################################
#clear objects from the workspace
rm(list=ls())

#Load packages
library(vegan)
library(ggfortify)
library(dplyr)
  
#Load and format data
data_meta <- read.table("./Taxonomy_data.txt", 
                        header = TRUE, sep = '\t', quote = "", row.names = 1)

# To select for subsets (e.g., FECAL or CECAL)
#data_meta <- data_meta[grep("FECAL", data_meta$Type),]
      
## Prepare data
#Extract metadata
meta <- data_meta[,(1:8)]
      
#Change "HET" to "5xfAD"
meta$Genotype <- replace(meta$Genotype, meta$Genotype=="HET", "5xfAD")
      
#Define group
meta$Group <- paste(meta$Genotype, meta$Cohort, meta$Type, sep="_")
      
#Extract data
OTUs <- data_meta[,-(1:8)]
is.numeric(OTUs)
sapply(OTUs, mode)

#Create relative abundance matrix
OTUs_rel <- OTUs/rowSums(OTUs)
#rownames(OTUs_rel) <- meta[,1]

#PCA of relative abundance data
pca_rel <- prcomp(OTUs_rel, scale. = FALSE)
      
autoplot(pca_rel, data = meta, shape = 'Cohort', colour = 'Genotype',
         loadings = F, loadings.colour = 'blue',
         loadings.label = F, loadings.label.size = 3)
      
#PCA of raw abundance data
pca_raw <- prcomp(OTUs, scale. = FALSE)
      
autoplot(pca_raw, data = meta, shape = 'Cohort', colour = 'Genotype',
         loadings = F, loadings.colour = 'blue',
         loadings.label = F, loadings.label.size = 3)
    ###Observed separation by Cohort but not Genotype
      
#PCA of betadiversity
beta <-vegdist(OTUs_rel, method="bray")
      
pca_rel_beta <- prcomp(beta, scale. = FALSE)
      
autoplot(pca_rel_beta, data = meta, shape = 'Cohort', colour = 'Genotype',
         loadings = F, loadings.colour = 'blue',
         loadings.label = F, loadings.label.size = 3)
    ###Observed separation by Cohort but not Genotype
      
##################### PERMANOVA #################################
#Load packages
library(vegan)
library(ggfortify)
library(dplyr)
      
#Load and format data
data_meta <- read.table("./Taxonomy_data.txt", 
                        header = TRUE, sep = '\t', quote = "", row.names = 1)

# To select for subsets (e.g., FECAL or CECAL)
#data_meta <- data_meta[grep("CECAL", data_meta$Type),]
#data_meta <- data_meta[grep("18mo", data_meta$Cohort),]

## Prepare data
#Extract metadata
meta <- data_meta[,(1:8)]
      
#Change "HET" to "5xfAD"
meta$Genotype <- replace(meta$Genotype, meta$Genotype=="HET", "5xfAD")
      
#Define group
meta$Group <- paste(meta$Genotype, meta$Cohort, meta$Type, sep="_")
      
#Extract data
OTUs <- data_meta[,-(1:8)]
is.numeric(OTUs)
sapply(OTUs, mode)
      
#Create relative abundance matrix
OTUs_rel <- OTUs/rowSums(OTUs)
#rownames(OTUs_rel) <- meta[,1]
  
###For All sample types
meta_test_1 <- data.frame(subset(meta, select = c(Genotype, Sex, Type, Cohort, Housing.ID)))
      
meta_test_1$Genotype <- as.factor(meta_test_1$Genotype)
meta_test_1$Sex <- as.factor(meta_test_1$Sex)
meta_test_1$Type <- as.factor(meta_test_1$Type)
meta_test_1$Cohort <- as.factor(meta_test_1$Cohort)
meta_test_1$Housing.ID <- as.factor(meta_test_1$Housing.ID)

res_1 <- adonis(formula = OTUs_rel ~ ., 
              data = meta_test_1, perm = 999, 
              method = "bray") 
      
hist(res_1$f.perms)
res_1

#All ages      
#  Call:
#  Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
#  Genotype     1    0.1614 0.16139   5.667 0.00596  0.001 ** 
#  Sex          1    0.7444 0.74441  26.140 0.02748  0.001 ***
#  Type         1    1.2861 1.28607  45.160 0.04748  0.001 ***
#  Cohort       3    4.3819 1.46063  51.290 0.16177  0.001 ***
#  Housing.ID  51   15.4449 0.30284  10.634 0.57018  0.001 ***
#  Residuals  178    5.0690 0.02848         0.18713           
#  Total      235   27.0877                 1.00000  

#18 mo only      
#  Call:
#Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
#Genotype    1    0.5544 0.55436  17.999 0.07287  0.001 ***
#Sex         1    0.3418 0.34183  11.098 0.04493  0.001 ***
#Type        1    0.7952 0.79525  25.820 0.10453  0.001 ***
#Housing.ID 12    4.2222 0.35185  11.424 0.55500  0.001 ***
#Residuals  55    1.6940 0.03080         0.22267           
#Total      70    7.6077                 1.00000 

#4 mo only      
#  Call:
#Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
#Genotype     1    0.0182 0.01818   0.840 0.00145  0.497    
#Sex          1    0.7257 0.72565  33.527 0.05796  0.001 ***
#Type         1    1.1997 1.19967  55.428 0.09583  0.001 ***
#Housing.ID  21    8.9090 0.42424  19.601 0.71164  0.001 ***
#Residuals   77    1.6666 0.02164         0.13312           
#Total      101   12.5191                 1.00000

###For CECAL/FECAL only
meta_test_2 <- data.frame(subset(meta, select = c(Genotype, Sex, Cohort, Housing.ID)))     

meta_test_2$Genotype <- as.factor(meta_test_2$Genotype)
meta_test_2$Sex <- as.factor(meta_test_2$Sex)
meta_test_2$Cohort <- as.factor(meta_test_2$Cohort)
meta_test_2$Housing.ID <- as.factor(meta_test_2$Housing.ID)
    
res_2 <- adonis(formula = OTUs_rel ~ ., 
                data = meta_test_2, perm = 999, 
                method = "bray") 
hist(res_2$f.perms)
res_2

#All Ages      
#  Call: CECAL
#              Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
#  Genotype    1    0.1314 0.13142   6.027 0.02311  0.001 ***
#  Sex         1    0.2619 0.26192  12.012 0.04605  0.001 ***
#  Cohort      1    0.9925 0.99254  45.517 0.17450  0.001 ***
#  Housing.ID 21    3.4078 0.16228   7.442 0.59915  0.001 ***
#  Residuals  41    0.8940 0.02181         0.15719           
#  Total      65    5.6878                 1.00000           

#18 mo only      
#  Call:
#Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
#Genotype    1   0.24614 0.24614  9.8614 0.09036  0.001 ***
#Sex         1   0.13137 0.13137  5.2631 0.04822  0.004 ** 
#Housing.ID 12   1.84744 0.15395  6.1680 0.67817  0.001 ***
#Residuals  20   0.49920 0.02496         0.18325           
#Total      34   2.72415                 1.00000 

#4 mo only      
#  Call:
#Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
#Genotype    1   0.01749 0.01749  1.1621 0.00886  0.343    
#Sex         1   0.42283 0.42283 28.0974 0.21429  0.001 ***
#Housing.ID  8   1.23190 0.15399 10.2326 0.62432  0.001 ***
#Residuals  20   0.30098 0.01505         0.15253           
#Total      30   1.97319                 1.00000  

#All Ages
#  Call: FECAL
#              Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
#  Genotype     1    0.0782 0.07823   2.489 0.00389  0.027 *  
#  Sex          1    0.5058 0.50583  16.090 0.02515  0.001 ***
#  Cohort       3    3.3976 1.13253  36.025 0.16892  0.001 ***
#  Housing.ID  51   12.5791 0.24665   7.846 0.62541  0.001 ***
#  Residuals  113    3.5524 0.03144         0.17662           
#  Total      169   20.1132                 1.00000 

#18 mo only      
#  Call:
#Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
#Genotype    1    0.3393 0.33929  8.5487 0.08296  0.001 ***
#Sex         1    0.2236 0.22362  5.6341 0.05467  0.001 ***
#Housing.ID 12    2.6937 0.22447  5.6557 0.65859  0.001 ***
#Residuals  21    0.8335 0.03969         0.20378           
#Total      35    4.0901                 1.00000 

#12 mo only      
#  Call:
#Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)    
#Genotype    1   0.07202 0.072017  2.8352 0.03826  0.007 ** 
#Sex         1   0.12227 0.122271  4.8136 0.06496  0.001 ***
#Housing.ID  6   1.00207 0.167011  6.5749 0.53239  0.001 ***
#Residuals  27   0.68583 0.025401         0.36438           
#Total      35   1.88219                  1.00000 

#8 mo only      
#  Call:
#Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)    
#Genotype    1    0.2061 0.206108  5.6352 0.02909  0.003 ** 
#Sex         1    0.2705 0.270503  7.3958 0.03817  0.001 ***
#Housing.ID 23    5.2563 0.228533  6.2483 0.74176  0.001 ***
#Residuals  37    1.3533 0.036575         0.19098           
#Total      62    7.0862                  1.00000 

#4 mo only      
#  Call:
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
#Genotype    1    0.0276 0.02761   0.986 0.00295  0.393    
#Sex         1    0.3926 0.39256  14.020 0.04197  0.001 ***
#Housing.ID 21    7.6169 0.36271  12.954 0.81437  0.001 ***
#Residuals  47    1.3160 0.02800         0.14070           
#Total      70    9.3530                 1.00000    

##################### NMDS #######################################
#clear objects from the workspace
rm(list=ls())

#Load packages
library(vegan)
library(ggfortify)
library(dplyr)
      
#Load and format data
data_meta <- read.table("./Taxonomy_data.txt", 
                        header = TRUE, sep = '\t', quote = "", row.names = 1)

# Use grep to select for subsets (e.g., FECAL or CECAL)
#data_meta <- data_meta[grep("CECAL", data_meta$Type),]
      
## Prepare data
#Extract metadata
meta <- data_meta[,(1:8)]
      
#Change "HET" to "5xfAD"
meta$Genotype <- replace(meta$Genotype, meta$Genotype=="HET", "5xfAD")

#Define group
meta$Group <- paste(meta$Genotype, meta$Cohort, meta$Type, sep="_")
      
#Extract data
OTUs <- data_meta[,-(1:8)]
is.numeric(OTUs)
sapply(OTUs, mode)
      
#Create relative abundance matrix
OTUs_rel <- OTUs/rowSums(OTUs)
#rownames(OTUs_rel) <- meta[,1]
      
#Get grouping information
grouping_info<- meta
      
#Get MDS stats
sol<-metaMDS(OTUs_rel, distance = "bray", k = 2, trymax = 100)
sol
    #Dimensions: 2 
    #Stress:     0.1969434 
    #Stress type 1, weak ties
    #Two convergent solutions found after 26 tries
#Define stress for plot
stress <- stress <- signif(sol$stress,digits=2)

#basic plots to test 
stressplot(sol)
plot(sol)
      
#Make a new data frame for coloring, and shape of points
NMDS=data.frame(x=sol$point[,1],y=sol$point[,2],
                Genotype=as.factor(grouping_info$Genotype),
                Cohort=as.factor(grouping_info$Cohort),
                Sex=as.factor(grouping_info$Sex),
                Type=as.factor(grouping_info$Type),
                Housing.ID=as.factor(grouping_info$Housing.ID),
                Group=as.factor(grouping_info$Group))
      
NMDS$Cohort <- factor(NMDS$Cohort, levels = c("18mo", "12mo", "8mo","4mo"))
shape_values<-c(21,24,22,3)
gwsC <- c('#d94801', '#2171b5') 

#Plot
ggplot(data=NMDS,aes(x,y, colour = Genotype))+
  theme_bw() +
  stat_ellipse(aes(group = Genotype), level = 0.95, show.legend = F, linetype = 2, size=1,  color = "black") +
  theme(plot.title = element_text(hjust = 0.5, size = 14,face = "bold"))+ 
  scale_x_continuous(name = "NMDS1") + scale_y_continuous(name = "NMDS2")+
  annotate("text", x=-0.75, y=-0.85, size = 4, label=stress)+ #add stress indicators
  geom_point(size=3, aes(shape=Cohort, fill=Genotype))+
  scale_shape_manual(values=shape_values) + 
  scale_color_manual(values=gwsC) + scale_fill_manual(values=gwsC)+
  theme(axis.text = element_text(size = 10,face = "bold"),
        axis.title = element_text(size = 10,face = "bold"))+
  theme(axis.text.x=element_text( size = 10, hjust = 0.5))

#ggsave("./NMDS.pdf", height = 3.5, width = 4.5)

##################### NMDS_Subsets #######################################
#clear objects from the workspace
rm(list=ls())

#Load packages
library(vegan)
library(ggfortify)
library(dplyr)

#Load and format data
data_meta <- read.table("./Taxonomy_data.txt", 
                        header = TRUE, sep = '\t', quote = "", row.names = 1)

# Use grep to select for subsets (e.g., FECAL or CECAL)
data_meta <- data_meta[grep("CECAL", data_meta$Type),]
data_meta <- data_meta[grep("18mo", data_meta$Cohort),]

## Prepare data
#Extract metadata
meta <- data_meta[,(1:8)]

#Change "HET" to "5xfAD"
meta$Genotype <- replace(meta$Genotype, meta$Genotype=="HET", "5xfAD")

#Define group
meta$Group <- paste(meta$Genotype, meta$Cohort, meta$Type, sep="_")
#meta$Group <- paste(meta$Genotype, meta$Type, sep="_")

#Extract data
OTUs <- data_meta[,-(1:8)]
is.numeric(OTUs)
sapply(OTUs, mode)

#Create relative abundance matrix
OTUs_rel <- OTUs/rowSums(OTUs)
#rownames(OTUs_rel) <- meta[,1]

#Get grouping information
grouping_info<- meta

#Get MDS stats
sol<-metaMDS(OTUs_rel, distance = "bray", k = 2, trymax = 100)
sol
sol$stress
stress <- signif(sol$stress,digits=2)
#Cecal Stress:     0.1788469 
#Fecal Stress:     0.1925399

#basic plots to test 
stressplot(sol)
plot(sol)

#Make a new data frame for coloring, and shape of points
NMDS=data.frame(x=sol$point[,1],y=sol$point[,2],
                Genotype=as.factor(grouping_info$Genotype),
                Cohort=as.factor(grouping_info$Cohort),
                Sex=as.factor(grouping_info$Sex),
                Type=as.factor(grouping_info$Type),
                Housing.ID=as.factor(grouping_info$Housing.ID),
                Group=as.factor(grouping_info$Group))

NMDS$Cohort <- factor(NMDS$Cohort, levels = c("18mo", "12mo", "8mo","4mo"))

#Plot (for all Cecal)
shape_values<-c(21, 3)
gwsC <- c('#d94801', '#2171b5')

ggplot(data=NMDS,aes(x,y, colour = Genotype))+
  theme_bw() +
  stat_ellipse(aes(group = Group), level = 0.95, show.legend = F, linetype = 2, size=1,  color = "black") +
  theme(plot.title = element_text(hjust = 0.5, size = 14,face = "bold"))+ 
  scale_x_continuous(name = "NMDS1") + scale_y_continuous(name = "NMDS2")+
  annotate("text", x=-0.7, y=-0.75, size = 4, label=stress)+ #add stress indicators
  geom_point(size=3, aes(shape=Cohort, fill=Genotype))+
  scale_shape_manual(values=shape_values) + 
  scale_color_manual(values=gwsC) + scale_fill_manual(values=gwsC)+
  theme(axis.text = element_text(size = 10,face = "bold"),
        axis.title = element_text(size = 10,face = "bold"))+
  theme(axis.text.x=element_text( size = 10, hjust = 0.5))+
  annotate("text", x=-0.2, y=0.75, size = 3, label="4 mo WT")+ #annotate groups
  annotate("text", x=0.7, y=0.6, size = 3, label="4 mo 5xfAD")+ #annotate groups
  annotate("text", x=-0.5, y=-0.5, size = 3, label="18 mo WT")+ #annotate groups
  annotate("text", x=0.15, y=-0.75, size = 3, label="18 mo 5xfAD") #annotate groups
  
#ggsave("./NMDS_Cecal_Group.pdf", height = 3.5, width = 4.5)

#Plot (for all Fecal)

shape_values<-c(21,24,22,3)
gwsC <- c('#d94801', '#2171b5') 
ggplot(data=NMDS,aes(x,y, colour = Genotype))+
  theme_bw() +
  stat_ellipse(aes(group = Group), level = 0.95, show.legend = F, linetype = 2, size=1,  color = "black") +
  theme(plot.title = element_text(hjust = 0.5, size = 14,face = "bold"))+ 
  scale_x_continuous(name = "NMDS1") + scale_y_continuous(name = "NMDS2")+
  annotate("text", x=-0.9, y=-1.5, size = 4, label=stress)+ #add stress indicators
  geom_point(size=3, aes(shape=Cohort, fill=Genotype))+
  scale_shape_manual(values=shape_values) + 
  scale_color_manual(values=gwsC) + scale_fill_manual(values=gwsC)+
  theme(axis.text = element_text(size = 10,face = "bold"),
        axis.title = element_text(size = 10,face = "bold"))+
  theme(axis.text.x=element_text( size = 10, hjust = 0.5))+
  annotate("text", x=-1, y=0.9, size = 3, label="4 mo WT")+ #annotate groups
  annotate("text", x=-1.1, y=-0.6, size = 3, label="4 mo 5xfAD")+ #annotate groups
  annotate("text", x=0.5, y=1, size = 3, label="8 mo WT")+ #annotate groups
  annotate("text", x=0.6, y=0.6, size = 3, label="8 mo 5xfAD")+ #annotate groups
  annotate("text", x=0.8, y=0, size = 3, label="12 mo WT")+ #annotate groups
  annotate("text", x=0.8, y=0.2, size = 3, label="12 mo 5xfAD")+ #annotate groups
  annotate("text", x=0, y=-1.4, size = 3, label="18 mo 5xfAD")+ #annotate groups
  annotate("text", x=-0.5, y=-1, size = 3, label="18 mo WT") #annotate groups

  #ggsave("./NMDS_Fecal_Group.pdf", height = 3.5, width = 4.5)

##Plot (for 18mo Fecal & Cecal)
ggplot(data=NMDS,aes(x,y, colour = Genotype))+
  theme_bw() +
  stat_ellipse(aes(group = Group), level = 0.95, show.legend = F, linetype = 2, size=1,  color = "black") +
  theme(plot.title = element_text(hjust = 0.5, size = 14,face = "bold"))+ 
  scale_x_continuous(name = "NMDS1") + scale_y_continuous(name = "NMDS2")+
  annotate("text", x=-0.9, y=-0.7, size = 4, label=stress)+ #add stress indicators
  geom_point(size=3, aes(shape=Type, fill=Genotype))+
  scale_shape_manual(values=shape_values) + 
  scale_color_manual(values=gwsC) + scale_fill_manual(values=gwsC)+
  theme(axis.text = element_text(size = 10,face = "bold"),
        axis.title = element_text(size = 10,face = "bold"))+
  theme(axis.text.x=element_text( size = 10, hjust = 0.5))+
#annotate("text", x=-1, y=1, size = 3, label="4 mo WT")+ #annotate groups
#annotate("text", x=1, y=1, size = 3, label="4 mo 5xfAD")+ #annotate groups
#annotate("text", x=-0.5, y=-0.5, size = 3, label="8 mo WT")+ #annotate groups
#annotate("text", x=0.15, y=-0.75, size = 3, label="8 mo 5xfAD")+ #annotate groups
#annotate("text", x=-0.2, y=0.75, size = 3, label="12 mo WT")+ #annotate groups
#annotate("text", x=0.7, y=0.6, size = 3, label="12 mo 5xfAD")+ #annotate groups
annotate("text", x=0.7, y=-0.6, size = 3, label="18 mo WT")+ #annotate groups
annotate("text", x=-1.1, y=0.7, size = 3, label="18 mo 5xfAD") #annotate groups

#ggsave("./NMDS_18mo_Fecal.pdf", height = 3.5, width = 4.5)

#ggsave("./NMDS_18mo_Cecal&Fecal_Group.pdf", height = 3.5, width = 4.5)

##Plot (for 4mo Fecal & Cecal)
ggplot(data=NMDS,aes(x,y, colour = Genotype))+
  theme_bw() +
  stat_ellipse(aes(group = Group), level = 0.95, show.legend = F, linetype = 2, size=1,  color = "black") +
  theme(plot.title = element_text(hjust = 0.5, size = 14,face = "bold"))+ 
  scale_x_continuous(name = "NMDS1") + scale_y_continuous(name = "NMDS2")+
  annotate("text", x=-0.8, y=-1, size = 4, label=stress)+ #add stress indicators
  geom_point(size=3, aes(shape=Type, fill=Genotype))+
  scale_shape_manual(values=shape_values) + 
  scale_color_manual(values=gwsC) + scale_fill_manual(values=gwsC)+
  theme(axis.text = element_text(size = 10,face = "bold"),
        axis.title = element_text(size = 10,face = "bold"))+
  theme(axis.text.x=element_text( size = 10, hjust = 0.5))
#annotate("text", x=-1, y=1, size = 3, label="4 mo WT")+ #annotate groups
#annotate("text", x=1, y=1, size = 3, label="4 mo 5xfAD")+ #annotate groups
#annotate("text", x=-0.5, y=-0.5, size = 3, label="8 mo WT")+ #annotate groups
#annotate("text", x=0.15, y=-0.75, size = 3, label="8 mo 5xfAD")+ #annotate groups
#annotate("text", x=-0.2, y=0.75, size = 3, label="12 mo WT")+ #annotate groups
#annotate("text", x=0.7, y=0.6, size = 3, label="12 mo 5xfAD")+ #annotate groups
#annotate("text", x=-0.5, y=-1.5, size = 3, label="18 mo WT")+ #annotate groups
#annotate("text", x=0.15, y=-1.5, size = 3, label="18 mo 5xfAD") #annotate groups

#ggsave("./NMDS_4mo_Cecal&Fecal_Group.pdf", height = 3.5, width = 4.5)

#Plot (for 18mo Cecal)
ggplot(data=NMDS,aes(x,y, colour = Genotype))+
  theme_bw() +
  stat_ellipse(aes(group = Genotype), level = 0.95, show.legend = F, linetype = 2, size=1,  color = "black") +
  theme(plot.title = element_text(hjust = 0.5, size = 14,face = "bold"))+ 
  scale_x_continuous(name = "NMDS1") + scale_y_continuous(name = "NMDS2")+
  annotate("text", x=-0.7, y=-0.7, size = 4, label=stress)+ #add stress indicators
  geom_point(size=3, aes(fill=Genotype))+
  scale_shape_manual(values=shape_values) + 
  scale_color_manual(values=gwsC) + scale_fill_manual(values=gwsC)+
  theme(axis.text = element_text(size = 10,face = "bold"),
        axis.title = element_text(size = 10,face = "bold"))+
  theme(axis.text.x=element_text( size = 10, hjust = 0.5))+
  annotate("text", x=-0.7, y=0.55, size = 3, label="18 mo WT")+ #annotate groups
  annotate("text", x=0.5, y=-0.65, size = 3, label="18 mo 5xfAD")+ #annotate groups
  ggtitle("NMDS of 18 month cecal samples")+
  theme(plot.title = element_text(hjust = 0, size = 12,face = "bold")) 
#ggsave("./NMDS_Cecal_18mo_Genotype.pdf", height = 3.5, width = 4.5)

#Plot (for 4mo Cecal)
ggplot(data=NMDS,aes(x,y, colour = Genotype))+
  theme_bw() +
  stat_ellipse(aes(group = Genotype), level = 0.95, show.legend = F, linetype = 2, size=1,  color = "black") +
  theme(plot.title = element_text(hjust = 0.5, size = 14,face = "bold"))+ 
  scale_x_continuous(name = "NMDS1") + scale_y_continuous(name = "NMDS2")+
  annotate("text", x=-0.5, y=-0.7, size = 4, label=stress)+ #add stress indicators
  geom_point(size=3, aes(fill=Genotype))+
  scale_shape_manual(values=shape_values) + 
  scale_color_manual(values=gwsC) + scale_fill_manual(values=gwsC)+
  theme(axis.text = element_text(size = 10,face = "bold"),
        axis.title = element_text(size = 10,face = "bold"))+
  theme(axis.text.x=element_text( size = 10, hjust = 0.5))
  #annotate("text", x=-0.7, y=0.55, size = 3, label="18 mo WT")+ #annotate groups
  #annotate("text", x=0.5, y=-0.65, size = 3, label="18 mo 5xfAD") #annotate groups

#ggsave("./NMDS_Cecal_4mo_Genotype.pdf", height = 3.5, width = 4.5)

#Plot (for 18mo Fecal)
ggplot(data=NMDS,aes(x,y, colour = Genotype))+
  theme_bw() +
  stat_ellipse(aes(group = Genotype), level = 0.95, show.legend = F, linetype = 2, size=1,  color = "black") +
  theme(plot.title = element_text(hjust = 0.5, size = 14,face = "bold"))+ 
  scale_x_continuous(name = "NMDS1") + scale_y_continuous(name = "NMDS2")+
  annotate("text", x=-1, y=-0.65, size = 4, label=stress)+ #add stress indicators
  geom_point(size=3, aes(fill=Genotype))+
  scale_shape_manual(values=shape_values) + 
  scale_color_manual(values=gwsC) + scale_fill_manual(values=gwsC)+
  theme(axis.text = element_text(size = 10,face = "bold"),
        axis.title = element_text(size = 10,face = "bold"))+
  theme(axis.text.x=element_text( size = 10, hjust = 0.5))+
  annotate("text", x=-1, y=0.7, size = 3, label="18 mo 5xfAD")+ #annotate groups
  annotate("text", x=0.6, y=-0.65, size = 3, label="18 mo WT")+ #annotate groups
  ggtitle("NMDS of 18 month fecal samples")+
  theme(plot.title = element_text(hjust = 0, size = 12,face = "bold")) 

#ggsave("./NMDS_Fecal_18mo_Genotype.pdf", height = 3.5, width = 4.5)

#Plot (for 12mo Fecal)
ggplot(data=NMDS,aes(x,y, colour = Genotype))+
  theme_bw() +
  stat_ellipse(aes(group = Genotype), level = 0.95, show.legend = F, linetype = 2, size=1,  color = "black") +
  theme(plot.title = element_text(hjust = 0.5, size = 14,face = "bold"))+ 
  scale_x_continuous(name = "NMDS1") + scale_y_continuous(name = "NMDS2")+
  annotate("text", x=1, y=-0.3, size = 4, label=stress)+ #add stress indicators
  geom_point(size=3, aes(fill=Genotype))+
  scale_shape_manual(values=shape_values) + 
  scale_color_manual(values=gwsC) + scale_fill_manual(values=gwsC)+
  theme(axis.text = element_text(size = 10,face = "bold"),
        axis.title = element_text(size = 10,face = "bold"))+
  theme(axis.text.x=element_text( size = 10, hjust = 0.5))
  #annotate("text", x=-1, y=0.7, size = 3, label="18 mo 5xfAD")+ #annotate groups
  #annotate("text", x=0.6, y=-0.65, size = 3, label="18 mo WT") #annotate groups

#ggsave("./NMDS_Fecal_12mo_Genotype.pdf", height = 3.5, width = 4.5)

#Plot (for 8mo Fecal)
ggplot(data=NMDS,aes(x,y, colour = Genotype))+
  theme_bw() +
  stat_ellipse(aes(group = Genotype), level = 0.95, show.legend = F, linetype = 2, size=1,  color = "black") +
  theme(plot.title = element_text(hjust = 0.5, size = 14,face = "bold"))+ 
  scale_x_continuous(name = "NMDS1") + scale_y_continuous(name = "NMDS2")+
  annotate("text", x=-1, y=-0.65, size = 4, label=stress)+ #add stress indicators
  geom_point(size=3, aes(fill=Genotype))+
  scale_shape_manual(values=shape_values) + 
  scale_color_manual(values=gwsC) + scale_fill_manual(values=gwsC)+
  theme(axis.text = element_text(size = 10,face = "bold"),
        axis.title = element_text(size = 10,face = "bold"))+
  theme(axis.text.x=element_text( size = 10, hjust = 0.5))
  #annotate("text", x=-1, y=0.7, size = 3, label="18 mo 5xfAD")+ #annotate groups
  #annotate("text", x=0.6, y=-0.65, size = 3, label="18 mo WT") #annotate groups

#ggsave("./NMDS_Fecal_8mo_Genotype.pdf", height = 3.5, width = 4.5)

#Plot (for 4mo Fecal)
ggplot(data=NMDS,aes(x,y, colour = Genotype))+
  theme_bw() +
  stat_ellipse(aes(group = Genotype), level = 0.95, show.legend = F, linetype = 2, size=1,  color = "black") +
  theme(plot.title = element_text(hjust = 0.5, size = 14,face = "bold"))+ 
  scale_x_continuous(name = "NMDS1") + scale_y_continuous(name = "NMDS2")+
  annotate("text", x=-1, y=-0.65, size = 4, label=stress)+ #add stress indicators
  geom_point(size=3, aes(fill=Genotype))+
  scale_shape_manual(values=shape_values) + 
  scale_color_manual(values=gwsC) + scale_fill_manual(values=gwsC)+
  theme(axis.text = element_text(size = 10,face = "bold"),
        axis.title = element_text(size = 10,face = "bold"))+
  theme(axis.text.x=element_text( size = 10, hjust = 0.5))
  #annotate("text", x=-1, y=0.7, size = 3, label="18 mo 5xfAD")+ #annotate groups
  #annotate("text", x=0.6, y=-0.65, size = 3, label="18 mo WT") #annotate groups

#ggsave("./NMDS_Fecal_4mo_Genotype.pdf", height = 3.5, width = 4.5)

##################### PCOA #######################################
#clear objects from the workspace
rm(list=ls())

#Load packages      
library(vegan)
library(ggfortify)
library(dplyr)
library(tibble)

#import data
data_meta <- as.data.frame(read.table("./Taxonomy_data.txt",
                                         row.names =1, 
                                         quote=NULL, 
                                         sep="\t", 
                                         stringsAsFactors=FALSE, 
                                         header= TRUE))
      
#Change "HET" to "5xfAD"
data_meta$Genotype <- replace(data_meta$Genotype, data_meta$Genotype=="HET", "5xfAD")

# Use grep to select for subsets (e.g., FECAL or CECAL)
data_meta <- data_meta[grep("CECAL", data_meta$Type),]
data_meta <- data_meta[grep("Male", data_meta$Sex),]
      
#Prepare relative abundance matrices
OTUs <- data_meta[,-(1:8)]
OTUs_rel <- OTUs/rowSums(OTUs)
meta <- data_meta[,(1:8)]
OTUs_rel_meta <- cbind(meta, OTUs_rel)
      
#Select for metadata
meta <- OTUs_rel_meta[,(1:8)]
meta$Genotype <- replace(meta$Genotype, meta$Genotype=="HET", "5xfAD")
meta$Group <- paste(meta$Cohort, meta$Genotype, sep="_")
      
#Create a distance matrix
bray = vegdist(OTUs_rel_meta[,9:ncol(OTUs_rel_meta)], method = "bray") 
      
#Here, the number of dimensions (k) is equal to n-1, where n is the number of OTUS
pcoa = cmdscale(bray, eig = T, k = nrow(OTUs_rel_meta[,9:ncol(OTUs_rel_meta)])-1, add = T) 
    
#Calculate eigen values
eig = eigenvals(pcoa) 
      
#To get the proportion of variance for PCO1 / PCO2
pcoa_var = eig/sum(eig) 
pcoa_var
      
#Grab only the points from the PCO object, merge it with metadata and remove samples that don't make sense.
points = as.data.frame(pcoa$points[,1:3]) %>% 
  merge(., meta, by = "row.names") #%>% 
      
#Add Loadings
#ord is some ordination object like an nmds or pco. 
#env is the abundance matrix. env is the OTUs_rel_meta excluding metadata
#The two have to match in samples.
      
# create an abundance matrix the excludes the metadata
env1 <- OTUs_rel_meta[,-(1:ncol(meta))]
      
#calculate: This step takes a long time. May want to save test2 as a seperate file for callup    
test2 = envfit(ord = pcoa, env = env1) 
      
#make the loadings
test3 = as.data.frame(scores(test2, display = "vectors")) 
test3 = test3 %>% rownames_to_column()
      
#To order by cohort
points$Cohort <- factor(points$Cohort, levels = c("18mo", "12mo", "8mo", "4mo")) # may not be needed
      
shape_values<-c(21,24,22,3)
gwsC <- c('#d94801', '#2171b5') 
      
#Plot For Fecal
ggplot(data = points)+
  aes(x = V1, y = V2, colour = Genotype) +
  theme_bw() +
  geom_segment(data = test3, aes(x = 0, xend = Dim1, y = 0, yend = Dim2), 
               arrow = arrow(length = unit(0.25, "cm")), colour = "grey90")+  #This line draws the loadings
  theme(plot.title = element_text(hjust = 0.5, size = 14,face = "bold"))+ 
  stat_ellipse(aes(group = Group), linetype = 2, size = 1,show.legend = F, color = "black") +
  #scale_x_continuous(name = "NMDS1") + scale_y_continuous(name = "NMDS2")+
  geom_point(size=3, aes(shape=Cohort, fill=Genotype))+
  scale_shape_manual(values=shape_values) + 
  scale_color_manual(values=gwsC) + 
  scale_fill_manual(values=gwsC)+
  labs(x = bquote(bold("PCo1 ("~.(round(pcoa_var[1]*100, digits = 1))~"%) ")), 
       y = bquote(bold("PCo2 ("~.(round(pcoa_var[2]*100, digits = 1))~"%) ")))+ #for labeling the axes with the proportion of variance
  theme(axis.text = element_text(size = 10,face = "bold"),
        axis.title = element_text(size = 10,face = "bold"))+
  annotate("text", x=0.85, y=0.6, size = 3, label="4 mo 5xfAD")+
  annotate("text", x=1.1, y=0.2, size = 3, label="4 mo 
  WT")+ 
  annotate("text", x=0.3, y=-0.8, size = 3, label="18 mo 5xfAD")+ 
  annotate("text", x=-0.5, y=-0.8, size = 3, label="18 mo 
  WT")+
  annotate("text", x=-0.5, y=-0.2, size = 3, label="12 mo
      5xfAD")+ 
  annotate("text", x=0.2, y=0.1, size = 3, label="12 mo 
  WT")+ 
  annotate("text", x=-0.7, y=0.9, size = 3, label="8 mo
  WT")+
  annotate("text", x=-0, y=0.8, size = 3, label="8 mo 5xfAD")
   
#ggsave("./PCoA_Fecal_Male_loadings.pdf", height = 3.5, width = 4.5)

shape_values<-c(21, 3)
gwsC <- c('#d94801', '#2171b5')

#Plot For Cecal
ggplot(data = points)+
  aes(x = V1, y = V2, colour = Genotype) +
  theme_bw() +
  geom_segment(data = test3, aes(x = 0, xend = Dim1, y = 0, yend = Dim2), 
               arrow = arrow(length = unit(0.25, "cm")), colour = "grey90")+  #This line draws the loadings
  theme(plot.title = element_text(hjust = 0.5, size = 14,face = "bold"))+ 
  stat_ellipse(aes(group = Group), linetype = 2, size = 1,show.legend = F, color = "black") +
  #scale_x_continuous(name = "NMDS1") + scale_y_continuous(name = "NMDS2")+
  geom_point(size=3, aes(shape=Cohort, fill=Genotype))+
  scale_shape_manual(values=shape_values) + 
  scale_color_manual(values=gwsC) + 
  scale_fill_manual(values=gwsC)+
  labs(x = bquote(bold("PCo1 ("~.(round(pcoa_var[1]*100, digits = 1))~"%) ")), 
       y = bquote(bold("PCo2 ("~.(round(pcoa_var[2]*100, digits = 1))~"%) ")))+ #for labeling the axes with the proportion of variance
  theme(axis.text = element_text(size = 10,face = "bold"),
        axis.title = element_text(size = 10,face = "bold"))+
  annotate("text", x=0.3, y=-0.6, size = 3, label="18 mo WT ")+
  annotate("text", x=-0.6, y=-0.6, size = 3, label="4 mo WT")+ 
  annotate("text", x=0.25, y=0.6, size = 3, label="18 mo 5xfAD")+
  annotate("text", x=-0.7, y=-0.2, size = 3, label="4 mo 
  5xfAD")

#ggsave("./PCoA_Cecal_Male_Loadings.pdf", height = 3.5, width = 4.5)

####3D plot (Just for fun)
library(plotly)
shape_values3D<-c(1,2,3,4)
gwsC3D <- c('#d94801', '#2171b5') 
plot_ly(x=points$V1, y=points$V2, z=points$V3, type="scatter3d", mode="markers", symbol = points$Cohort) 
#color=points$Genotype, colors = gwsC3D, symbol = points$Cohort, symbols = shape_values3D)
      
##################### Random Forest ##################################
#clear objects from the workspace
rm(list=ls())

#Load packages
library(dplyr)
library(tibble)    
library(rfPermute)
library(ggplot2)
library(ggpubr)
      
#import data
data_meta <- as.data.frame(read.table("./Taxonomy_data.txt",
                                           row.names =1, 
                                           quote=NULL, 
                                           sep="\t", 
                                           stringsAsFactors=FALSE, 
                                           header= TRUE))
      
## Prepare data
#Extract metadata
meta <- data_meta[,(1:8)]
      
#Change "HET" to "5xfAD"
meta$Genotype <- replace(meta$Genotype, meta$Genotype=="HET", "5xfAD")

# Select for subsets (e.g., FECAL or CECAL)
meta = filter(meta, Cohort == "18mo")
#meta = filter(meta, Type == "FECAL")
#meta = filter(meta, Genotype == "WT")

#Define group (e.g., Cohort_Line_SampleType)
#By changing the grouping categories RF can be performed on different 
#separation categories. Among other things, we evaluated differences 
#with respect to: 
#meta$Group <- paste(meta$Type, sep="_")#sample type: 
#meta$Group <- paste(meta$Sex, sep="_")#sex: 
#meta$Group <- paste(meta$Cohort, sep="_")#age: 
meta$Group <- paste(meta$Genotype, sep="_")#genotype: 
#meta$Group <- paste(meta$Genotype, meta$Type, sep="_")#genotype & sample type
#meta$Group <- paste(meta$Cohort, meta$Genotype, sep="_")#genotype & cohort 

#Extract data
OTUs <- data_meta[,-(1:8)]
      
#Coverted to relative abundance
OTUs <- OTUs/rowSums(OTUs) 
      
#RF doesn't like having more features than samples, so remove low abundance taxa.
OTUs <- OTUs[,colMeans(OTUs) >= 0.001]
      
# Random forest only uses on meta factor. This line merges it but you can do that which ever way you want.
OTUs <- merge(meta, OTUs, by = "row.names")
      
# Re-establish ronames
OTUs <- column_to_rownames(OTUs, var = "Row.names")
      
# Remove unneeded meta
OTUs <- OTUs[,-(1:8)]
      
# Save and re-create OTU table to make numeric
write.table(OTUs, "./OTU_num.txt", sep = "\t", row.names = TRUE, quote=FALSE)
OTU_num <- as.data.frame(read.table("./OTU_num.txt",
                                    row.names =1, 
                                    quote=NULL,
                                    check.names=FALSE,
                                    sep="\t", 
                                    stringsAsFactors=TRUE, 
                                    header= TRUE))
      
#Order Group
#OTU_num$Group <- factor(OTU_num$Group, levels=c('4mo','8mo','12mo','18mo'))
#OTU_num$Group <- factor(OTU_num$Group, levels=c('4mo','18mo'))

#set seed
set.seed(77)

# Perform Random Forest
rf_out <- rfPermute(formula = Group ~ ., 
                    data = OTU_num, proximity = T, importance = F, 
                    ntree = 601, num.cores = 32)
      
### View Results
# Plot Confusion Matrix heatmap
confplot <- plotConfMat(rf_out, title= NULL, plot = TRUE)
#ggsave("./RF_18mo_Cecal&Fecal_Geno_ConfPlot.pdf", height = 2, width = 3)
      
# Create a plot of Random Forest proximity scores using multi-dimensional scaling. 
proxplot <- proximityPlot(rf_out, legend.loc = "bottom",   dim.x = 1, dim.y = 2)
      
safe_colorblind_palette <- c("#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
                             "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")

proximityPlot( rf_out, dim.x = 1, dim.y = 2,
               c('#d94801', '#2171b5'),#class.cols = c('#d94801','#8c2d04', '#2171b5','#084594'), # =NULL , #NULL or c('#d94801', '#2171b5') or safe_colorblind_palette
               legend.type = "none", #"legend", "label", "none"
               legend.loc = "right", #"top", "bottom", "left", "right"
               point.size = 2,
               circle.size = 0,
               circle.border = 0,
               group.type = "ellipse", #"ellipse", "hull", "contour", "none"
               group.alpha = 0.1,
               ellipse.level = 0.95,
               n.contour.grid = 100,
               label.size = 4,
               label.alpha = 0.7,
               plot = TRUE
)
#ggsave("./RF_18mo_Cecal&Fecal_Geno_ProxPlot.pdf", height = 3, width = 4)

# Dotchart of variable importance as measured by a Random Forest
varplot <- varImpPlot(rf_out, type = 1, n.var = 20)

#For barplot of variable importance
imp <- as.data.frame(varImpPlot(rf_out)) # save the varImp object as data frame
imp$variable <- rownames(imp)
rownames(imp) <- NULL  
imp <- arrange(imp, desc(MeanDecreaseAccuracy))
imp_20 <- imp[(1:20),]
imp_20$variable <- factor(imp_20$variable, levels=imp_20$variable)

ggplot(imp_20, aes(x=variable, weight=MeanDecreaseAccuracy, fill=variable))+ 
  geom_bar() + theme_bw() + coord_flip() +  
  aes(x = reorder(variable, desc(variable)))+
  xlab(NULL) + ylab("Mean Decrease Accuracy")+ 
  scale_fill_discrete(guide="none")+
  theme(axis.text.x=element_text(size=7, face = "bold"),
        axis.text.y=element_text(size=6, face = "bold"),
        axis.title=element_text(size=7, face = "bold"))

#ggsave("./RF_18mo_Cecal&Fecal_Geno_varplot.pdf", height = 2.25, width = 3)

#end RF
##################### Linear Mixed Effects Model ############################
#clear objects from the workspace
rm(list=ls())

#Load packages
library(dplyr)
library(tibble)   
library(lmerTest)
library(reshape2)
library(tidyverse)
      
#Read in your OTU table. If you've already got relative abundances, skip this part.
data_meta<-as.data.frame(read.table("./Taxonomy_data.txt",
                                    row.names =1, 
                                    quote=NULL, 
                                    sep="\t", 
                                    stringsAsFactors=FALSE, 
                                    header= TRUE))

# Select for subsets (e.g., FECAL or CECAL)
data_meta <- filter(data_meta, Cohort == "18mo")
data_meta <- filter(data_meta, Type == "FECAL")
#data_meta <- filter(data_meta, Sex == "Female")
      
#Change "HET" to "5xfAD"
data_meta$Genotype <- replace(data_meta$Genotype, data_meta$Genotype == "HET", "5xfAD")
      
#Select for data & transpose
OTUs <- t(data_meta[,-(1:8)])

#Make relative abundances.
#OTUs<- OTUs[complete.cases(OTUs), ]
column_sums <- colSums(OTUs)
OTUs_rel <- apply(OTUs, 1, '/', column_sums)
      
#Read in metadata.
var<-data_meta[,(1:8)]
      
attach(var)
      
#Transpose the OTU table & filter out low abundance taxa. 
OTUs_rel <- as.data.frame(t(OTUs_rel))
OTUs_rel$mean<-rowMeans(OTUs_rel)
OTUs_rel <- OTUs_rel %>%
  rownames_to_column("Sample") %>%
  filter(OTUs_rel$mean > 0.001) %>%
  subset(., select = -c(mean))
OTUs_rel <- setNames(data.frame(t(OTUs_rel[,-1])),OTUs_rel[,1])
      
#Make a longform data table with melt.
melted<-melt(OTUs_rel)

#Run a Linear Mixed Effects Model (LMER). 
#It's ~lmer(OTU abundance ~ Metadata category for comparison + (Random effect)
lmertest<-melted %>%
  split(.$variable) %>%
  map(~lmer(value ~ Genotype + (1|as.factor(Housing.ID)), data = .)) %>%
  map(summary) %>% 
  map(coef)
      
lmertest[["Turicibacter.sp..H121"]]
      
#Pull out p values from your test. The "10" is just the 10th position in the output. Yours might be different.
pvalues<-map_dbl(lmertest,10)
p_ordered<-(sort(pvalues))
padj<-p.adjust(p_ordered, method = "fdr")
padj
padj_df<-as.data.frame(padj)
      
#Filter out insignificant taxa.
padj_df<-padj_df %>%
  filter(padj < 0.05)
padj_df<-rownames_to_column(padj_df, "names")
padj_df
      
#Set up a dataframe for plotting.
plotting<-melted %>%
  filter(variable %in% padj_df$names)
      
detach(var)
  
#significant by lmer for 18mo all (no sample type separation)
#names         padj
#1          Turicibacter.sp..H121 2.784729e-07
#2             Romboutsia.ilealis 1.070107e-04
#3          Lactobacillus.sp..P38 1.921330e-04
#4      Ligilactobacillus.murinus 1.921330e-04
#5     Ligilactobacillus.animalis 1.921330e-04
#6 Burkholderiales.bacterium.YL45 3.493773e-03
      
#Significant species @ 18 mo (FECAL) by LMER test
#names        padj
#1       Turicibacter.sp..H121 0.001274444
#2   Ligilactobacillus.murinus 0.029149724
#3       Lactobacillus.sp..P38 0.029149724
#4  Ligilactobacillus.animalis 0.029149724
#5          Romboutsia.ilealis 0.029149724
#6 Adlercreutzia.equolifaciens 0.032350436      
      
#Significant species @ 18 mo (CECAL) by LMER test
#names        padj
#1       Turicibacter.sp..H121 0.001022314
#2          Romboutsia.ilealis 0.007205026
#3          Alistipes.megaguti 0.023318701
#4            Alistipes.dispar 0.023318701
#5       Alistipes.onderdonkii 0.023318701
#6            Alistipes.shahii 0.023318701
#7   Ligilactobacillus.murinus 0.023318701
#8       Lactobacillus.sp..P38 0.023318701
#9        Alistipes.finegoldii 0.023318701
#10         Duncaniella.sp..C9 0.023318701
#11 Ligilactobacillus.animalis 0.023893540
#12         Alistipes.communis 0.024475364
#13         Duncaniella.sp..B8 0.025939592
#14    Muribaculum.intestinale 0.040784296
#15    Lactobacillus.johnsonii 0.040784296
      
##################### Relative Abundance Plots #############################
#clear objects from the workspace
rm(list=ls())

#Load packages
library(dplyr)
library(tibble) 
library(ggplot2)
library(ggpubr)
      
#Read in OTU table
data_meta<-as.data.frame(read.table("./Taxonomy_data.txt",
                                    row.names =1, 
                                    quote=NULL, 
                                    sep="\t", 
                                    stringsAsFactors=FALSE, 
                                    header= TRUE))
      
#Metadata
metadata <- data_meta[,(1:8)]
      
#Select for subsets
#metadata = filter(metadata, Sex == "Female") 
      
#Change "HET" to "5xfAD"
metadata$Genotype <- replace(metadata$Genotype, metadata$Genotype=="HET", "5xfAD")
      
#Order by age
metadata$Cohort = factor(metadata$Cohort, levels=c('4mo','8mo','12mo','18mo'))
      
#Data
rawdata <- data_meta[,-(1:8)]
    
#Generate relative abundance data
reldata <- rawdata/rowSums(rawdata)
      
#merge reldata with metadata
reldata_meta <- merge(metadata, reldata, by=0)    
      
#=====Species Plots=======#
#Test plot
MOI = "Turicibacter.sp..H121"
reldata_meta[[MOI]]
reldata_meta[[MOI]] <- as.numeric(reldata_meta[[MOI]])
ggplot(data = reldata_meta) + 
  aes(x = Cohort, y = .data[[MOI]]*100, fill = Genotype, label = Mouse_ID) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values=c('#d94801', '#2171b5')) +
  geom_point(size = 2, colour = "black", alpha = 1, 
             position = position_jitterdodge(jitter.width = 0.1),
             aes(fill=Genotype, shape=Sex)) + 
  facet_wrap(Type + Sex~., nrow = 1, scale = "free_x") +
  theme_bw() +
  theme(title =  element_text(size = 8,face = "bold"))+
  theme(legend.text =  element_text(hjust = 0, size = 8,face = "bold"))+
  theme(legend.title = element_text(hjust = 0, size = 8,face = "bold"))+ 
  theme(axis.text.y = element_text(hjust = 0, size = 7, face = "bold"))+
  theme(axis.text.x = element_text(hjust = 0.5, size = 7, face = "bold"))+
  theme(axis.title.y = element_text(size = 6, face = "bold"))+
  theme(axis.title.x = element_text(size = 6, face = "bold"))+
  #theme(legend.position = "none")+
  #geom_text(hjust=0, vjust=0)+
  labs(x = "Age", y = "Relative Abundance (%)", fill = "Genotype", 
       title = MOI)+
  stat_compare_means(method = "wilcox.test", label = "p.signif",
                     label.x = 1.5, vjust = 0.5, size = 4, hide.ns = T) ## package ggpubr lets us put significance tests directly on the plot!
    
#ggsave(path = "./", filename=paste(MOI,".pdf",sep=""), 
#       device = "pdf", height = 3, width = 4.5,)
      
#plot all "selects" using for loop and save to file. 
#"selects" is a data frame containing the list of analytes in column 1 
selects = as.data.frame(c("Turicibacter.sp..H121", "Romboutsia.ilealis", "Alistipes.megaguti", 
                          "Alistipes.dispar", "Alistipes.onderdonkii", "Alistipes.shahii", 
                          "Ligilactobacillus.murinus", "Lactobacillus.sp..P38", "Alistipes.finegoldii", 
                          "Duncaniella.sp..C9", "Ligilactobacillus.animalis", "Alistipes.communis", 
                          "Duncaniella.sp..B8","Muribaculum.intestinale", "Lactobacillus.johnsonii", 
                          "Adlercreutzia.equolifaciens", "Burkholderiales.bacterium.YL45"))
colnames(selects) <- "moi"
      
for (i in 1:nrow(selects)) {
  MOI = selects[i,1]
  reldata_meta[[MOI]] <- as.numeric(reldata_meta[[MOI]])
  plot <- ggplot(data = reldata_meta) + 
    aes(x = Cohort, y = .data[[MOI]]*100, fill = Genotype, label = Mouse_ID) +
    geom_boxplot(outlier.shape = NA) +
    scale_fill_manual(values=c('#d94801', '#2171b5')) +
    geom_point(size = 2, colour = "black", alpha = 1, 
               position = position_jitterdodge(jitter.width = 0.1),
               aes(fill=Genotype, shape=Sex)) + 
    facet_wrap(Type + Sex~., nrow = 1, scale = "free_x") +
    theme_bw() +
    theme(title =  element_text(size = 8,face = "bold"))+
    theme(legend.text =  element_text(hjust = 0, size = 8,face = "bold"))+
    theme(legend.title = element_text(hjust = 0, size = 8,face = "bold"))+ 
    theme(axis.text.y = element_text(hjust = 0, size = 7, face = "bold"))+
    theme(axis.text.x = element_text(hjust = 0.5, size = 7, face = "bold"))+
    theme(axis.title.y = element_text(size = 6, face = "bold"))+
    theme(axis.title.x = element_text(size = 6, face = "bold"))+
    theme(legend.position = "none")+
    #geom_text(hjust=0, vjust=0)+
    labs(x = "Age", y = "Relative Abundance (%)", fill = "Genotype", 
         title = MOI)+
    stat_compare_means(method = "wilcox.test", label = "p.signif",
                       label.x = 1.5, vjust = 0.5, size = 4, hide.ns = T) ## package ggpubr lets us put significance tests directly on the plot!
  plot
  ggsave(plot, path = "./", filename=paste(MOI,".pdf",sep=""), 
         device = "pdf", height = 3, width = 4.5)
  }
