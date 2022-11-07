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

######## Heatmap of Most Abundant Fecal Metabolites ################
#Adapted from Dave Tang's pheatmap example
# https://davetang.org/muse/2018/05/15/making-a-heatmap-in-r-with-the-pheatmap-package/

#clear objects from the workspace
rm(list=ls())

# load package
library(pheatmap)
library(tibble)
library(viridis)
library(dplyr)

#Load and format data
data_meta <- as.data.frame(read.table("./FecalMX_data.txt",
                                      fill = TRUE,
                                      row.names = 1,
                                      quote=NULL, 
                                      sep="\t", 
                                      stringsAsFactors=FALSE, 
                                      header= TRUE))

#Extract metadata
meta <- data_meta[,(1:8)]

#Order data.
data <- data_meta %>%
  mutate(Sex =  factor(Sex, levels = c("Female", "Male"))) %>%
  arrange(Sex) %>%
  mutate(Genotype =  factor(Genotype, levels = c("WT", "5xfAD"))) %>%
  arrange(Genotype) %>%
  mutate(Cohort = factor(Cohort, levels = c("4mo", "8mo", "12mo", "18mo"))) %>%
  arrange(Cohort)

data <- data[,-(1:8)]

##selecting most abundant (e.g., Top100)
top_metabs <- as.data.frame(colSums(data)) %>%
  rownames_to_column() %>%
  top_n(100)

top <- data[colnames(data) %in% top_metabs$rowname]

#transpose for better plotting
top_trans <- t(top)

#normalize using std normal distribution
normfunction <- function(x){
  (x - mean(x)) / sd(x)
}

top_norm <- as.data.frame(t(apply(top_trans, 1, normfunction)))

#Design a key for column annotation
key <- as.data.frame((cbind(meta$Sex, meta$Genotype, meta$Cohort)), 
                     row.names = row.names(meta))

colnames(key) <- c("Sex", "Genotype", "Cohort")

#Arrange key
key <- key %>%
  mutate(Cohort =  factor(Cohort, levels = c("4mo", "8mo", "12mo", "18mo"))) %>%
  arrange(Cohort)

key <- key %>%
  mutate(Genotype =  factor(Genotype, levels = c("WT", "5xfAD"))) %>%
  arrange(Genotype)

#plot
pheatmap(top_norm)
pheatmap(top_norm, annotation_col =  key)
p2 <- pheatmap(top_norm, clustering_method = 'average', cluster_rows = TRUE, cluster_cols = FALSE,
               gaps_col = TRUE, main = '', show_colnames = FALSE, show_rownames = TRUE, 
               fontsize = 10, cutree_cols = 0, cutree_rows = 5, angle_col = 270, annotation_col =  key,
               legend = TRUE)
               #, filename = "./pheatmap_5xFAD_ClusterRow_Tep684.pdf", cellwidth = 8, cellheight = 8)

######## Correlation between Most Abundant Fecal Metabolites & Selected Microbes ################

#clear objects from the workspace
rm(list=ls())

library(tidyverse)
library(reshape)
library(ggdendro)
library(viridis)
library(grid)
library(Hmisc)
library(pheatmap)
library(vegan)

otus <- as.data.frame(read.table("./Taxonomy_data.txt",
                                 row.names =1, 
                                 quote=NULL, 
                                 sep="\t", 
                                 stringsAsFactors=FALSE, 
                                 header= TRUE))

#Select for data only
otu_data <- otus[,-(1:8)]

#Load selected species
selects<-read.table("./selects.txt", 
                    header = TRUE, sep = '\t', quote = "", row.names = 1)

select_otus <- otu_data[colnames(otu_data) %in% selects$names]

###Loading & processing the metabolite data
metab_meta <- as.data.frame(read.table("./FecalMX_data.txt",
                                       fill = TRUE,
                                       row.names = 1,
                                       quote=NULL, 
                                       sep="\t", 
                                       stringsAsFactors=FALSE, 
                                       header= TRUE))

# create a metadata table that contains only the samples present in both metabolite & microbe data
metab_meta <- rownames_to_column(metab_meta, var = "SampleID")
select_otus <- rownames_to_column(select_otus, var = "SampleID")
merged_data <- merge(metab_meta, select_otus,
              by.x="SampleID", by.y="SampleID", 
              all.x = FALSE, all.y = FALSE, sort = TRUE)
merged_data <- column_to_rownames(merged_data, var = "SampleID")

#Select for metabolite data only
mx <- merged_data[,(9:(ncol(merged_data)-17))]

##selecting most abundant metabolites (e.g., top 20)
top_mx <- as.data.frame(colSums(mx)) %>%
  rownames_to_column() %>%
  top_n(100)

top_mx <- mx[colnames(mx) %in% top_mx$rowname]

#Select for OTUData
select_otus_comp <- merged_data[,((ncol(merged_data)-16):ncol(merged_data))]

#Select for metadata
map <- merged_data[,(2:8)]

#correlation between select otus and 100 most abundant metabolites
correlation <- cor(top_mx,select_otus_comp, method = c('spearman'))

pheatmap(correlation, cluster_cols = TRUE, clustering_method = 'average', cluster_rows = TRUE,
         gaps_col = FALSE, color = viridis(10), show_colnames = TRUE, show_rownames = TRUE, 
         fontsize = 8, fontsize_col = 8, cutree_cols = 5, angle_col = 45, 
         cutree_rows = 6, fontsize_row = 5, silent = FALSE)
           #, width = 10, height = 12, filename = "./Outputs/spearman_selects.pdf")

#Mantel significance test
mx.dist <- vegdist(top_mx) # Bray-Curtis
otu.dist <- vegdist(select_otus_comp)
mantel(mx.dist, otu.dist)
    #Mantel statistic r: 0.4716 
    #Significance: 0.001 
mantel(mx.dist, otu.dist, method = "spearman", permutations = 9999, na.rm = TRUE)
    #Mantel statistic r: 0.2446 
    #Significance: 0.0057 

######## Total abundance plot ####################################################

#clear objects from the workspace
rm(list=ls())

library(dplyr)
library(ggplot2)
library(ggpubr)

#Load and format data
data_meta <- as.data.frame(read.table("./FecalMX_data.txt",
                                      fill = TRUE,
                                      row.names = 1,
                                      quote=NULL, 
                                      sep="\t", 
                                      stringsAsFactors=FALSE, 
                                      header= TRUE))

#Create data matrix
metabolites <- data_meta[,-(1:8)]

#Create metadata matrix 
meta <- data_meta[,(1:8)]

#Define Group
meta$Group <- paste(meta$Cohort, meta$Genotype, sep="_")

#generate a count of metabolites present in each sample and name column "metabcount"
metab.count <- data.frame(rowSums(metabolites))
#metab.count <- rename(metab.count, metabcount = rowSums.metabolites.)
colnames(metab.count) <- "metabcount"

#Join metab.count & metadata
metab.count.meta <-cbind(meta, metab.count)
metab.count.group <- data.frame(metab.count.meta[,-(1:8)])

# calculate summary statistics
data_summary <- metab.count.group %>%
  group_by(Group) %>%
  summarise(mean = mean(metabcount), 
            median = median(metabcount),
            se = (sd(metabcount)) / sqrt(length(metabcount)))

#Shapiro-Wilk normality test
  #if p < 0.05 the data is not normally distributed
  #if p > 0.05 the data is normally distributed
shapiro.test(data_summary$mean)
#W = 0.85584, p-value = 0.249
# data is normally distributed

#2-way anova by genotype and age ("*" assumes a synergistic effect)
seq.test <- aov(metabcount  ~ Cohort*Genotype, data = metab.count.meta)
summary(seq.test)
#Results
#             Df    Sum Sq   Mean Sq F value Pr(>F)
#Cohort           2 2.105e+14 1.053e+14   1.435  0.252
#Genotype         1 2.349e+13 2.349e+13   0.320  0.575
#Cohort:Genotype  2 3.043e+14 1.522e+14   2.074  0.141
#Residuals       34 2.494e+15 7.335e+13                          

#order
metab.count.meta <- metab.count.meta %>%
  mutate(Genotype =  factor(Genotype, levels = c("WT", "5xfAD"))) %>%
  arrange(Genotype) %>%
  mutate(Cohort =  factor(Cohort, levels = c("4mo", "8mo", "12mo"))) %>%
  arrange(Cohort)

#Plot
ggplot(data = metab.count.meta) +
  aes(x = Cohort, y = metabcount, fill = Genotype) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values=c('#d94801', '#2171b5')) +
  geom_point(size = 2, colour = "black", alpha = 1, 
             position = position_jitterdodge(jitter.width = 0.1),
             aes(fill=Genotype, shape=Sex)) + 
  #facet_wrap(Line~., nrow = 1, scales = "free_x") +
  #facet_wrap(Sex+Type+Cohort~., nrow = 1) +
  scale_y_log10() +
  theme_bw() +
  theme(title =  element_text(size = 10,face = "bold"))+
  theme(legend.text =  element_text(hjust = 0, size = 10,face = "bold"))+
  theme(legend.title = element_text(hjust = 0, size = 10,face = "bold"))+ 
  theme(axis.text.y = element_text(hjust = 0, size = 10, face = "bold"))+
  theme(axis.text.x = element_text(hjust = 0.5, size = 10, face = "bold"))+
  theme(axis.title.y = element_text(size = 10, face = "bold"))+
  theme(axis.title.x = element_text(size = 10, face = "bold"))+
  labs(x = NULL, y = "Total Counts (iTIC; Log10)", fill = "Genotype", 
       title = "Total Metabolite Abundance (iTIC)")+ 
  stat_compare_means(method = "wilcox.test", label = "p.signif",
                     label.x = 1.5, vjust = 0.5, size = 4, hide.ns = F)

#ggsave("./Outputs/Counts_5xfAD_iTIC.pdf", height = 3.5, width = 4.5)

######## Shannon Diversity #######################################################

#clear objects from the workspace
rm(list=ls())

library(dplyr)
library(ggplot2)
library(ggpubr)
library(vegan)

#import data
data_meta <- as.data.frame(read.table("./FecalMX_data.txt",
                                      fill = TRUE,
                                      row.names = 1,
                                      quote=NULL, 
                                      sep="\t", 
                                      stringsAsFactors=FALSE, 
                                      header= TRUE))

#Create data matrix
metabolites <- data_meta[,-(1:8)]

#Create metadata matrix 
meta <- data_meta[,(1:8)]

#Define Group
meta$Group <- paste(meta$Cohort, meta$Genotype, sep="_")

#generate a count of metabolites present in each sample and name column "metabcount"
metab.count <- data.frame(rowSums(metabolites))
#metab.count <- rename(metab.count, metabcount = rowSums.metabolites.)
colnames(metab.count) <- "metabcount"

#Join metab.count & metadata
metab.count.meta <-cbind(meta, metab.count)

#Calculate shannon diversity
shannon <- diversity(metabolites, "shannon")
is.numeric(shannon)

#Shapiro-Wilk normality test
#if p < 0.05 the data is not normally distributed
#if p > 0.05 the data is normally distributed
shapiro.test(shannon)
#data:  shannon
#W = 0.87133, p-value = 0.0003098
    #the data is NOT normally distributed

shan.alpha <- data.frame(cbind(meta$Group, meta$Genotype, meta$Cohort,
                               meta$Housing.ID, meta$Sex, meta$Sample_Type, shannon))

colnames(shan.alpha) <- c('Group', 'Genotype', 'Cohort', 'HousingID', 'Sex', 'SampleType', 'shannon')
sapply(shan.alpha, is.factor)
shan.alpha$Group <- as.factor(shan.alpha$Group)
shan.alpha$Genotype <- as.factor(shan.alpha$Genotype)
shan.alpha$Cohort <- as.factor(shan.alpha$Cohort)
shan.alpha$HousingID <- as.factor(shan.alpha$HousingID)
shan.alpha$Sex <- as.factor(shan.alpha$Sex)
shan.alpha$SampleType <- as.factor(shan.alpha$SampleType)

is.numeric(shan.alpha$shannon)
shan.alpha$shannon <- as.numeric(shan.alpha$shannon)

kruskal.test(shannon ~ Group, data = shan.alpha)
#Kruskal-Wallis chi-squared = 7.8359, df = 5, p-value = 0.1655
kruskal.test(shannon ~ Genotype, data = shan.alpha)
#Kruskal-Wallis chi-squared = 1.6837, df = 1, p-value = 0.1944
kruskal.test(shannon ~ Cohort, data = shan.alpha)
#Kruskal-Wallis chi-squared = 0.25171, df = 2, p-value = 0.8817
kruskal.test(shannon ~ HousingID, data = shan.alpha)
#Kruskal-Wallis chi-squared = 20.955, df = 22, p-value = 0.5235
kruskal.test(shannon ~ Sex, data = shan.alpha)
#Kruskal-Wallis chi-squared = 0.97726, df = 1, p-value = 0.3229

# p < 0.05 some group medians are significantly different
# p > 0.05 no group medians are significantly differnent

#order
shan.alpha <- shan.alpha %>%
  mutate(Genotype =  factor(Genotype, levels = c("5xfAD", "WT"))) %>%
  arrange(Genotype) %>%
  mutate(Cohort =  factor(Cohort, levels = c("4mo", "8mo", "12mo"))) %>%
  arrange(Cohort)

#plot
ggplot(data = shan.alpha) +
  aes(x = Cohort, y = shannon, fill = Genotype) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values=c('#d94801', '#2171b5')) +
  geom_point(size = 2, colour = "black", alpha = 0.5, 
             position = position_jitterdodge(jitter.width = 0.1),
             aes(fill=Genotype, shape=Sex)) +
  #geom_point(size = 2, colour = "black", alpha = 0.5, position = position_jitterdodge(jitter.width = 0)) +
  #facet_wrap(Sex~., nrow = 1, scales = "free_x") +
  theme(axis.text = element_text(size = 10,face = "bold"),
        axis.title = element_text(size = 10,face = "bold"))+
  theme(legend.title = element_text(hjust = 0, size = 10,face = "bold"))+ 
  theme(axis.text.y = element_text(hjust = 0, size = 10, face = "bold"))+
  theme(axis.text.x = element_text(hjust = 0, size = 10, face = "bold"))+
  theme_bw() +
  labs(x = "Age", y = "Shannon Diversity", fill = "Genotype") + 
  stat_compare_means(method = "wilcox.test", label = "p.signif", label.x = 1.5,
                     vjust = 0, size = 4, hide.ns = F) ## package ggpubr lets us put significance tests directly on the plot!

#ggsave("./Outputs/Shannon_5xfAD_iTIC.pdf", height = 3.5, width = 4.5)

######## PCA #####################################################################

#clear objects from the workspace
rm(list=ls())

library(vegan)
library(ggfortify)
library(dplyr)

#import data
data_meta <- as.data.frame(read.table("./FecalMX_data.txt",
                                      fill = TRUE,
                                      row.names = 1,
                                      quote=NULL, 
                                      sep="\t", 
                                      stringsAsFactors=FALSE, 
                                      header= TRUE))

#Create data matrix
metabolites <- data_meta[,-(1:8)]

#Create metadata matrix 
meta <- data_meta[,(1:8)]

#PCA of raw abundance data
pca_res <- prcomp(metabolites, scale. = FALSE)

autoplot(pca_res, data = meta, shape = 'Cohort', colour = 'Genotype',
         loadings = FALSE, loadings.colour = 'blue',
         loadings.label = FALSE, loadings.label.size = 3)

#replot with betadiversity
beta <-vegdist(metabolites, method="bray")

pca_res_beta <- prcomp(beta, scale. = FALSE)

autoplot(pca_res_beta, data = meta, shape = 'Cohort', colour = 'Genotype',
         loadings = FALSE, loadings.colour = 'blue',
         loadings.label = FALSE, loadings.label.size = 3)


######## PERMANOVA ###############################################################

#clear objects from the workspace
rm(list=ls())

#Load packages
library(vegan)
library(ggfortify)
library(dplyr)

#import data
data_meta <- as.data.frame(read.table("./FecalMX_data.txt",
                                      fill = TRUE,
                                      row.names = 1,
                                      quote=NULL, 
                                      sep="\t", 
                                      stringsAsFactors=FALSE, 
                                      header= TRUE))

#Create data matrix
metabolites <- data_meta[,-(1:8)]

#Create metadata matrix 
meta <- data_meta[,(1:8)]

#ADONIS
meta_test <- data.frame(subset(meta, 
                               select = c(Genotype, Cohort, Sex, Housing.ID)))

meta_test$Genotype <- as.factor(meta_test$Genotype)
meta_test$Sex <- as.factor(meta_test$Sex)
meta_test$Cohort <- as.factor(meta_test$Cohort)
meta_test$Housing.ID <- as.factor(meta_test$Housing.ID)

res <- adonis(formula = metabolites ~ ., 
              data = meta_test, perm = 999, 
              method = "bray") 
res

#Call
#           Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)   
#Genotype    1    0.0575 0.057478 0.59921 0.01386  0.905   
#Cohort      2    0.4028 0.201422 2.09984 0.09711  0.007 **
#Sex         1    0.1008 0.100777 1.05061 0.02429  0.368   
#Housing.ID 19    2.0523 0.108016 1.12607 0.49475  0.290   
#Residuals  16    1.5348 0.095923         0.36999          
#Total      39    4.1482                  1.00000          

######## NMDS ####################################################################

#clear objects from the workspace
rm(list=ls())

#Load packages
library(vegan)
library(ggfortify)
library(dplyr)

#import data
data_meta <- as.data.frame(read.table("./FecalMX_data.txt",
                                      fill = TRUE,
                                      row.names = 1,
                                      quote=NULL, 
                                      sep="\t", 
                                      stringsAsFactors=FALSE, 
                                      header= TRUE))

#Create data matrix
metabolites <- data_meta[,-(1:8)]

#Create metadata matrix 
meta <- data_meta[,(1:8)]

#Get MDS stats
sol<-metaMDS(metabolites, distance = "bray", k = 2, trymax = 100)
sol
stress <- "stress = 0.089"      
#Dimensions: 2 
#Stress:     0.08886375 
#Stress type 1, weak ties
#Two convergent solutions found after 20 tries

#basic plots to test 
stressplot(sol)
plot(sol)

#Make a new data frame, and put Group and Factor information there, to be 
#useful for coloring, and shape of points
NMDS=data.frame(x=sol$point[,1],y=sol$point[,2],Genotype=as.factor(meta$Genotype),
                Cohort=as.factor(meta$Cohort))

#order
NMDS$Cohort <- factor(NMDS$Cohort, levels = c("12mo", "8mo","4mo"))
NMDS$Genotype <- factor(NMDS$Genotype, levels = c("WT","5xfAD")) 

shape_values<-c(21,24,22,3)
gwsC <- c('#d94801', '#2171b5') 

#scale_shape_manual(values = c(21,24,22,23), name = "Cohort")

ggplot(data=NMDS,aes(x,y, colour=Genotype))+
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 14,face = "bold"))+ 
  scale_x_continuous(name = "NMDS1") + scale_y_continuous(name = "NMDS2")+
  annotate("text", x=0.25, y=-0.1, size = 4, label=stress)+ #add stress indicators
  geom_point(size=3, aes(shape=Cohort, fill=Genotype))+
  stat_ellipse(linetype = 2, aes(group = Cohort), show.legend = F, color = "black") +
  scale_shape_manual(values=shape_values) + 
  scale_color_manual(values=gwsC) + scale_fill_manual(values=gwsC)+
  theme(axis.text = element_text(size = 10,face = "bold"),
        axis.title = element_text(size = 10,face = "bold"))+
  theme(axis.text.x=element_text( size = 10, hjust = 0.5))

#p <- p + labs(title = "NMDS of Batch 1 Unfiltered Relative Abundance")
#p <- p + theme(panel.background = element_rect(fill = 'grey'))

shape_values<-c(21,24,22,3)
gwsC <- c('#d94801', '#2171b5') 

#Probable Convergence Failure

#ggsave("./Outputs/NMDS.pdf", height = 3.5, width = 4.5)

######## PCOA ####################################################################

#clear objects from the workspace
rm(list=ls())

#Load Libraries
library(vegan)
library(ggfortify)
library(dplyr)
library(randomForest)
library(tibble)

#import data
data_meta <- as.data.frame(read.table("./FecalMX_data.txt",
                                      fill = TRUE,
                                      row.names = 1,
                                      quote=NULL, 
                                      sep="\t", 
                                      stringsAsFactors=FALSE, 
                                      header= TRUE))

#Create data matrix
metabolites <- data_meta[,-(1:8)]

#Create metadata matrix 
meta <- data_meta[,(1:8)]
meta$Group <- paste(meta$Cohort, meta$Genotype, sep="_")

#Create a distance matrix
bray = vegdist(metabolites, method = "bray") 

#Here, the number of dimensions (k) is equal to n-1, where n is the number of metabolites
pcoa = cmdscale(bray, eig = T, k = nrow(metabolites[,1:ncol(metabolites)])-1, add = T) 

#Calculate eigen values
eig = eigenvals(pcoa) 

#Compute the proportion of variance for PCO1 / PCO2
pcoa_var = eig/sum(eig) 
pcoa_var

#Grab only the points from the PCO object and merge with metadata.
points = as.data.frame(pcoa$points[,1:3]) %>% 
  merge(., meta, by = "row.names") #%>% 

#Add Loadings
#ord is some ordination object like an nmds or pco. 
#env is the abundance matrix. env is the metabolites excluding metadata
#The two have to match in samples.

#create an abundance matrix the excludes the metadata
env1 <- metabolites

#calculate: This step takes a long time.   
test2 = envfit(ord = pcoa, env = env1) 

#make the loadings
test3 = as.data.frame(scores(test2, display = "vectors")) 
test3 = test3 %>% rownames_to_column()

#To order data
points$Cohort <- factor(points$Cohort, levels = c("12mo", "8mo", "4mo")) 
points$Genotype <- factor(points$Genotype, levels = c("5xfAD", "WT")) 

gwsC <- c('#d94801', '#2171b5')
shape_values<-c(24,22,3)

#Plot with loadings & ellipses
ggplot(data = points) +
  aes(x = V1, y = V2) +
  theme_bw() +
  geom_segment(data = test3, aes(x = 0, xend = Dim1, y = 0, yend = Dim2), arrow = arrow(length = unit(0.25, "cm")), colour = "grey90") +  #This line draws the loadings
  #geom_text(data = test3, aes(x = Dim1, y = Dim2, label = rowname), size = 4, alpha = .3, inherit.aes = FALSE) + 
  geom_point(aes(pch = Cohort, fill = Genotype), size = 3, alpha = 1) + #sets shape
  scale_shape_manual(values = c(21,22,23,24,3,4,8,9), name = "Cohort") + # from https://ggplot2.tidyverse.org/reference/scale_shape.html
  guides(color = F, fill = guide_legend(override.aes = list(shape = 21))) +
  stat_ellipse(linetype = 2, aes(group = Group), show.legend = F, color = "black") +
  scale_fill_manual(values = gwsC) +
  labs(color = "Sample type", x = bquote("PCo1 ("~.(round(pcoa_var[1]*100, digits = 1))~"%) "), 
       y = bquote("PCo2 ("~.(round(pcoa_var[2]*100, digits = 1))~"%) "))
#title = "PCoA analyis of 5xfAD and WT fecal microbiome")  #This line is important for labeling the axes with the proportion of variance.

#ggsave("./PCOA_5xfAD_iTIC_Loadings&Ellipses.pdf", height = 3.5, width = 4.5)

#3D plot
library(plotly)
plot_ly(x=points$V1, y=points$V2, z=points$V3, type="scatter3d", mode="markers", color=points$Cohort)

######## Random Forest ###########################################################

#clear objects from the workspace
rm(list=ls())

library(rfPermute) #https://cran.r-project.org/web/packages/rfPermute/rfPermute.pdf
library(vegan)
library(ggfortify)
library(dplyr)
library(tidyr)
library(tibble)

#import data
data_meta <- as.data.frame(read.table("./FecalMX_data.txt",
                                      fill = TRUE,
                                      row.names = 1,
                                      quote=NULL, 
                                      sep="\t", 
                                      stringsAsFactors=FALSE, 
                                      header= TRUE))

#Extract metadata
meta <- data_meta[,(1:8)]

#Define group (e.g., Cohort_Line_SampleType)
#meta$Group <- paste(meta$Cohort, meta$Genotype, sep="_")
meta$Group <- paste(meta$Cohort, meta$Genotype, sep="_")

#Extract data
metabolites <- data_meta[,-(1:8)]

#RF doesnt like having more features than samples, so remove low abundance data
metabolites <- metabolites[,colMeans(metabolites) >= 1000]

# Random forest only uses on metadata factor.
metabolites <- merge(meta, metabolites, by = "row.names")

# Re-establish ronames
metabolites <- column_to_rownames(metabolites, var = "Row.names")

# Remove unneeded metadata
metabolites <- metabolites[,-(1:8)]

# Save and re-create OTU table to make numeric
write.table(metabolites, "./temp.txt", sep = "\t", row.names = TRUE, quote=FALSE)
OTU_num <- as.data.frame(read.table("./temp.txt",
                                    row.names =1, 
                                    quote=NULL,
                                    check.names=FALSE,
                                    sep="\t", 
                                    stringsAsFactors=TRUE, 
                                    header= TRUE))

# Perform Random Forest
rf_out <- rfPermute(formula = Group ~ ., 
                    data = OTU_num, proximity = T, importance = F, 
                    ntree = 601, num.cores = 32)

# View Results
# Plot Confusion Matrix heatmap
confplot <- plotConfMat(rf_out, title= NULL, plot = TRUE)
#ggsave("./RF_5xfAD_12mo_ConfPlot.pdf", height = 3.5, width = 4.5)

# Create a plot of Random Forest proximity scores using multi-dimensional scaling. 

safe_colorblind_palette <- c("#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
                             "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")

proximityPlot( rf_out, dim.x = 1, dim.y = 2,
               #class.cols =c('#d94801', '#2171b5') , #NULL or c('#d94801', '#2171b5') or safe_colorblind_palette
               legend.type = "legend", #"legend", "label", "none"
               legend.loc = "top", #"top", "bottom", "left", "right"
               point.size = 2,
               circle.size = 5,
               circle.border = 1,
               group.type = "ellipse", #"ellipse", "hull", "contour", "none"
               group.alpha = 0.25,
               ellipse.level = 0.95,
               n.contour.grid = 100,
               label.size = 4,
               label.alpha = 0.7,
               plot = TRUE
)
#ggsave("./RF_5xfAD_12mo_ProxPlot.pdf", height = 3.5, width = 4.5)

# Dotchart of variable importance as measured by a Random Forest
varplot <- varImpPlot(rf_out, type = 1, n.var = 20)

#pdf("./RF_5xfAD_12mo_VarPlot.pdf", height = 4.5, width = 4.5)  
#varImpPlot(rf_out, type = 1, n.var = 20, main = NULL, cex = 0.8)
#dev.off()

######## LMER ####################################################################

#clear objects from the workspace
rm(list=ls())

library(vegan)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(lmerTest)
library(reshape2)

#Read in your OTU table. 
data_meta <- as.data.frame(read.table("./FecalMX_data.txt",
                                      fill = TRUE,
                                      row.names = 1,
                                      quote=NULL, 
                                      sep="\t", 
                                      stringsAsFactors=FALSE, 
                                      header= TRUE))

#Select for 12mo Samples Only
data_meta <- filter(data_meta, Cohort == "12mo")

#Select for Female Samples Only
data_meta <- filter(data_meta, Sex == "Female")

#Select for data
data <- data_meta[,-(1:8)]
data<-t(data)

#Make relative abundances
data<- data[complete.cases(data), ]
column_sums <- colSums(data)
data_norm <- apply(data, 1, '/', column_sums)

#Read in your metadata.
var<-data_meta[,(1:8)]

attach(var)

#Transpose the OTU table and filter out low abundance mx 
test <- t(data_norm)
test <- as.data.frame(test)
test$mean <- rowMeans(test)
test <- test %>%
  rownames_to_column("Sample") %>%
  filter(test$mean > 0.0001) %>%
  subset(., select = -c(mean))
test = setNames(data.frame(t(test[,-1])),test[,1])

#Make a longform data table with melt.
melted<-melt(test)

#Run a Linear Mixed Effects Model (LMER). 
#It's ~lmer(OTU abundance ~ Metadata category for comparison + (Random effect)
lmertest<-melted %>%
  split(.$variable) %>%
  map(~lmer(value ~ Genotype + (1|as.factor(Housing.ID)), data = .)) %>%
  map(summary) %>% 
  map(coef)

lmertest[["Cytosine"]]

#Pull out pvalues 
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

#Save padj_df for downstream analysis
write.table(padj_df, "./MXselects_12moFemale.txt", sep = "\t", row.names = TRUE, quote=FALSE)

detach(var)

######## Abundance Plots #########################################################

#clear objects from the workspace
rm(list=ls())

library("tidyverse")
library("readxl")
library("tibble")
library("ggpubr")
library("ggplot2")

#Read in OTU table
data_meta <- as.data.frame(read.table("./FecalMX_data.txt",
                                      fill = T,
                                      row.names = 1,
                                      check.names = F,
                                      quote=NULL, 
                                      sep="\t", 
                                      stringsAsFactors=F, 
                                      header= T))

#Order by age
data_meta$Cohort = factor(data_meta$Cohort, levels=c('4mo','8mo','12mo'))

#Order by Genotype
data_meta$Genotype = factor(data_meta$Genotype, levels=c("5xfAD","WT"))

#Order by Sex
data_meta$Sex = factor(data_meta$Sex, levels=c("Female","Male"))

#Species Plots

#Significant metabolites for 12 mo females by LMER test
#    names	padj
#    1	Nialamide	9.4223188715275e-15
#    2	Fluorene	0.000955530044829756
#    3	Ala.Pro.Lys	0.00501756101899522
#    4	Ferulic.acid	0.0165791639835986
#    5	X4.Acetylbutyric.acid	0.0165791639835986
#    6	Disaccharide.16	0.0311825180328338

#Nialamide
#Simple by time
ggplot(data = data_meta) +
  aes(x = Cohort, y = Nialamide, fill = Genotype) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values=c('#d94801', '#2171b5')) +
  geom_point(size = 2, colour = "black", alpha = 1, 
             position = position_jitterdodge(jitter.width = 0.1),
             aes(fill=Genotype, shape=Sex)) +
  #facet_wrap(Sex~., nrow = 1) +
  theme_bw() +
  theme(title =  element_text(size = 10,face = "bold"))+
  theme(legend.text =  element_text(hjust = 0, size = 10,face = "bold"))+
  theme(legend.title = element_text(hjust = 0, size = 10,face = "bold"))+ 
  theme(axis.text.y = element_text(hjust = 0, size = 10, face = "bold"))+
  theme(axis.text.x = element_text(hjust = 0, size = 10, face = "bold"))+
  theme(axis.title.y = element_text(size = 10, face = "bold"))+
  theme(axis.title.x = element_text(size = 10, face = "bold"))+
  #theme(legend.position = "none")+
  labs(x = "Age", y = "Abundance (iTIC)", fill = "Genotype", 
       title = "Nialamide")+
  stat_compare_means(method = "wilcox.test", label = "p.signif",
                     label.x = 1.5, vjust = 0.5, size = 4, hide.ns = F)

#ggsave("./Nialamide.pdf", height = 3, width = 3.5)          

#Fluorene
ggplot(data = data_meta) +
  aes(x = Cohort, y = Fluorene, fill = Genotype) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values=c('#d94801', '#2171b5')) +
  geom_point(size = 2, colour = "black", alpha = 1, 
             position = position_jitterdodge(jitter.width = 0.1),
             aes(fill=Genotype, shape=Sex)) +
  #facet_wrap(Sex~., nrow = 1) +
  theme_bw() +
  theme(title =  element_text(size = 10,face = "bold"))+
  theme(legend.text =  element_text(hjust = 0, size = 10,face = "bold"))+
  theme(legend.title = element_text(hjust = 0, size = 10,face = "bold"))+ 
  theme(axis.text.y = element_text(hjust = 0, size = 10, face = "bold"))+
  theme(axis.text.x = element_text(hjust = 0, size = 10, face = "bold"))+
  theme(axis.title.y = element_text(size = 10, face = "bold"))+
  theme(axis.title.x = element_text(size = 10, face = "bold"))+
  theme(legend.position = "none")+
  labs(x = "Age", y = "Abundance (iTIC)", fill = "Genotype", 
       title = "Fluorene")+
  stat_compare_means(method = "wilcox.test", label = "p.signif",
                     label.x = 1.5, vjust = 0.5, size = 4, hide.ns = F)

#ggsave("./Outputs/mxPlots/0104/Fluorene.pdf", height = 3, width = 3.5)


#Ala-Pro-Lys (col 144)
ggplot(data = data_meta) +
  aes(x = Cohort, y = data_meta[,144], fill = Genotype) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values=c('#d94801', '#2171b5')) +
  geom_point(size = 2, colour = "black", alpha = 1, 
             position = position_jitterdodge(jitter.width = 0.1),
             aes(fill=Genotype, shape=Sex)) +
  #facet_wrap(Sex~., nrow = 1) +
  theme_bw() +
  theme(title =  element_text(size = 10,face = "bold"))+
  theme(legend.text =  element_text(hjust = 0, size = 10,face = "bold"))+
  theme(legend.title = element_text(hjust = 0, size = 10,face = "bold"))+ 
  theme(axis.text.y = element_text(hjust = 0, size = 10, face = "bold"))+
  theme(axis.text.x = element_text(hjust = 0, size = 10, face = "bold"))+
  theme(axis.title.y = element_text(size = 10, face = "bold"))+
  theme(axis.title.x = element_text(size = 10, face = "bold"))+
  theme(legend.position = "none")+
  labs(x = "Age", y = "Abundance (iTIC)", fill = "Genotype", 
       title = "Ala-Pro-Lys")+
  stat_compare_means(method = "wilcox.test", label = "p.signif",
                     label.x = 1.5, vjust = 0.5, size = 4, hide.ns = F)

#ggsave("./Ala-Pro-Lys.pdf", height = 3, width = 3.5)

#Ferulic.acid (col 244)
ggplot(data = data_meta) +
  aes(x = Cohort, y = data_meta[,244], fill = Genotype) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values=c('#d94801', '#2171b5')) +
  geom_point(size = 2, colour = "black", alpha = 1, 
             position = position_jitterdodge(jitter.width = 0.1),
             aes(fill=Genotype, shape=Sex)) +
  #facet_wrap(Sex~., nrow = 1) +
  theme_bw() +
  theme(title =  element_text(size = 10,face = "bold"))+
  theme(legend.text =  element_text(hjust = 0, size = 10,face = "bold"))+
  theme(legend.title = element_text(hjust = 0, size = 10,face = "bold"))+ 
  theme(axis.text.y = element_text(hjust = 0, size = 10, face = "bold"))+
  theme(axis.text.x = element_text(hjust = 0, size = 10, face = "bold"))+
  theme(axis.title.y = element_text(size = 10, face = "bold"))+
  theme(axis.title.x = element_text(size = 10, face = "bold"))+
  theme(legend.position = "none")+
  labs(x = "Age", y = "Abundance (iTIC)", fill = "Genotype", 
       title = "Ferulic acid")+
  stat_compare_means(method = "wilcox.test", label = "p.signif",
                     label.x = 1.5, vjust = 0.5, size = 4, hide.ns = F)

#ggsave("./Ferulic.acid.pdf", height = 3, width = 3.5)   

#4.Acetylbutyric.acid (col 81)
ggplot(data = data_meta) +
  aes(x = Cohort, y = data_meta[,81], fill = Genotype) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values=c('#d94801', '#2171b5')) +
  geom_point(size = 2, colour = "black", alpha = 1, 
             position = position_jitterdodge(jitter.width = 0.1),
             aes(fill=Genotype, shape=Sex)) +
  #facet_wrap(Sex~., nrow = 1) +
  theme_bw() +
  theme(title =  element_text(size = 10,face = "bold"))+
  theme(legend.text =  element_text(hjust = 0, size = 10,face = "bold"))+
  theme(legend.title = element_text(hjust = 0, size = 10,face = "bold"))+ 
  theme(axis.text.y = element_text(hjust = 0, size = 10, face = "bold"))+
  theme(axis.text.x = element_text(hjust = 0, size = 10, face = "bold"))+
  theme(axis.title.y = element_text(size = 10, face = "bold"))+
  theme(axis.title.x = element_text(size = 10, face = "bold"))+
  theme(legend.position = "none")+
  labs(x = "Age", y = "Abundance (iTIC)", fill = "Genotype", 
       title = "4-Acetylbutyric acid")+
  stat_compare_means(method = "wilcox.test", label = "p.signif",
                     label.x = 1.5, vjust = 0.5, size = 4, hide.ns = F)

#ggsave("./4.Acetylbutyric.acid.pdf", height = 3, width = 3.5)    

#Disaccharide.16
ggplot(data = data_meta) +
  aes(x = Cohort, y = Disaccharide.16, fill = Genotype) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values=c('#d94801', '#2171b5')) +
  geom_point(size = 2, colour = "black", alpha = 1, 
             position = position_jitterdodge(jitter.width = 0.1),
             aes(fill=Genotype, shape=Sex)) +
  #facet_wrap(Sex~., nrow = 1) +
  theme_bw() +
  theme(title =  element_text(size = 10,face = "bold"))+
  theme(legend.text =  element_text(hjust = 0, size = 10,face = "bold"))+
  theme(legend.title = element_text(hjust = 0, size = 10,face = "bold"))+ 
  theme(axis.text.y = element_text(hjust = 0, size = 10, face = "bold"))+
  theme(axis.text.x = element_text(hjust = 0, size = 10, face = "bold"))+
  theme(axis.title.y = element_text(size = 10, face = "bold"))+
  theme(axis.title.x = element_text(size = 10, face = "bold"))+
  theme(legend.position = "none")+
  labs(x = "Age", y = "Abundance (iTIC)", fill = "Genotype", 
       title = "Disaccharide.16")+
  stat_compare_means(method = "wilcox.test", label = "p.signif",
                     label.x = 1.5, vjust = 0.5, size = 4, hide.ns = F)

#ggsave("./Disaccharide.16.pdf", height = 3, width = 3.5) 

#Tryptophan
ggplot(data = data_meta) +
  aes(x = Cohort, y = Tryptophan, fill = Genotype) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values=c('#d94801', '#2171b5')) +
  geom_point(size = 2, colour = "black", alpha = 1, 
             position = position_jitterdodge(jitter.width = 0.1),
             aes(fill=Genotype, shape=Sex)) +
  #facet_wrap(Sex~., nrow = 1) +
  theme_bw() +
  theme(title =  element_text(size = 10,face = "bold"))+
  theme(legend.text =  element_text(hjust = 0, size = 10,face = "bold"))+
  theme(legend.title = element_text(hjust = 0, size = 10,face = "bold"))+ 
  theme(axis.text.y = element_text(hjust = 0, size = 10, face = "bold"))+
  theme(axis.text.x = element_text(hjust = 0, size = 10, face = "bold"))+
  theme(axis.title.y = element_text(size = 10, face = "bold"))+
  theme(axis.title.x = element_text(size = 10, face = "bold"))+
  theme(legend.position = "none")+
  labs(x = "Age", y = "Abundance (iTIC)", fill = "Genotype", 
       title = "Tryptophan (NS)")+
  stat_compare_means(method = "wilcox.test", label = "p.signif",
                     label.x = 1.5, vjust = 0.5, size = 4, hide.ns = F)

#ggsave("./Tryptophan.pdf", height = 3, width = 3.5) 

#5-Hydroxy-3-indoleacetic acid (col 104)
ggplot(data = data_meta) +
  aes(x = Cohort, y = data_meta[,104], fill = Genotype) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values=c('#d94801', '#2171b5')) +
  geom_point(size = 2, colour = "black", alpha = 1, 
             position = position_jitterdodge(jitter.width = 0.1),
             aes(fill=Genotype, shape=Sex)) +
  #facet_wrap(Sex~., nrow = 1) +
  theme_bw() +
  theme(title =  element_text(size = 10,face = "bold"))+
  theme(legend.text =  element_text(hjust = 0, size = 10,face = "bold"))+
  theme(legend.title = element_text(hjust = 0, size = 10,face = "bold"))+ 
  theme(axis.text.y = element_text(hjust = 0, size = 10, face = "bold"))+
  theme(axis.text.x = element_text(hjust = 0, size = 10, face = "bold"))+
  theme(axis.title.y = element_text(size = 10, face = "bold"))+
  theme(axis.title.x = element_text(size = 10, face = "bold"))+
  theme(legend.position = "none")+
  labs(x = "Age", y = "Abundance (iTIC)", fill = "Genotype", 
       title = "5-Hydroxy-3-indoleacetic acid")+
  stat_compare_means(method = "wilcox.test", label = "p.signif",
                     label.x = 1.5, vjust = 0.5, size = 4, hide.ns = F)

#ggsave("./5-Hydroxy-3-indoleacetic acid.pdf", height = 3, width = 3.5) 
