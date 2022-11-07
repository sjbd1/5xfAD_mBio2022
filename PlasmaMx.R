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

##################### Heatmap of Most Abundant Plasma Metabolites ##############
#Adapted from Dave Tang's pheatmap example
# https://davetang.org/muse/2018/05/15/making-a-heatmap-in-r-with-the-pheatmap-package/

#clear objects from the workspace
rm(list=ls())

# load package
library(pheatmap)
library(tibble)
library(viridis)
library(dplyr)
library(RColorBrewer)

#Load and format data
data_meta <- as.data.frame(read.table("./PlasmaMX_data.txt",
                                      fill = TRUE,
                                      row.names = 1,
                                      quote=NULL, 
                                      sep="\t", 
                                      stringsAsFactors=FALSE, 
                                      header= TRUE))

#Extract metadata
meta <- data_meta[,(1:10)]

#Order and select for data.
data <- data_meta %>%
  mutate(Sex =  factor(Sex, levels = c("Female", "Male"))) %>%
  arrange(Sex)

data <- data %>%
  mutate(Genotype =  factor(Genotype, levels = c("WT", "5xfAD"))) %>%
  arrange(Genotype)

data <- data %>% 
  mutate(Cohort = factor(Cohort, levels = c("12mo", "18mo"))) %>%
  arrange(Cohort)

data <- data[,-(1:10)]

#set NA = 0
data[is.na(data)] <- 0

##selecting most abundant (e.g. Top100)
top_metabs <- as.data.frame(colSums(data)) %>%
  rownames_to_column() %>%
  top_n(137)

top <- data[colnames(data) %in% top_metabs$rowname]

#transpose for better plotting
top_trans <- t(top)

#normalize
#the heatmap.2() function has a parameter for scaling the rows; this can be easily implemented.
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
  mutate(Cohort =  factor(Cohort, levels = c("12mo", "18mo"))) %>%
  arrange(Cohort)

key <- key %>%
  mutate(Genotype =  factor(Genotype, levels = c("WT", "5xfAD"))) %>%
  arrange(Genotype)

#plot
p1 <- pheatmap(top_norm, 
               clustering_method = 'average', cluster_rows = TRUE, cluster_cols = FALSE,
               main ='Plasma Metabolites', show_colnames = FALSE, show_rownames = TRUE, 
               annotation_col =  key, gaps_col = F, border_color = NA,
               fontsize = 6, fontsize_col = 4, fontsize_row = 4, 
               cutree_cols = 5, cutree_rows = 10, angle_col = 270, 
               legend = TRUE)#, filename = "./Outputs/pheatmap_P180Metabs_ClusterRow_All_Norm.pdf", cellheight = 7, cellwidth = 3.25)

cp <- colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu")), #name = "RdYlBu" name = "Spectral"
                       bias = 1)(100)

library(RColorBrewer)
#https://colorbrewer2.org/#type=diverging&scheme=PuOr&n=10
my_colour = list(
  #'Sample_Type' = c(FECAL = "#91bfdb"),
  #Cohort = c('4mo' = "#d7191c", '8mo' = "#fdae61", '12mo'="#abd9e9", '12mo'="#fdae61"),
  Cohort = c('12mo'="#abd9e9", '18mo'="#fdae61"),
  Genotype = c('5xfAD' = "#d94801", 'WT' = "#2171b5"),
  Sex = c(Female = "#80cdc1", Male="#000000" )
)
# define color pallet
cp <- colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu")), #name = "RdYlBu" name = "Spectral"
                       bias = 1)(100)

#plot
pheatmap(top_norm, clustering_method = 'average', cluster_rows = TRUE, cluster_cols = FALSE,
         gaps_col = FALSE, show_colnames = FALSE, show_rownames = TRUE, border_color = NA,
         fontsize = 6, fontsize_row = 4, fontsize_col = 4, 
         cutree_cols = 5, cutree_rows = 10, angle_col = 270, 
         annotation_legend = TRUE, annotation_col = key, annotation_colors = my_colour, 
         color = cp,legend = TRUE, legend_breaks = c(-3,-1,1,3,5,7,9), legend_labels = NA, 
         annotation_names_row = TRUE, annotation_names_col = TRUE,
         main = "Plasma Metabolites")

##################### Total Concentration of Validated Metabolites ##############

#clear objects from the workspace
rm(list=ls())

library(vegan)
library(ggfortify)
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggpubr)

#import data
data_meta <- as.data.frame(read.table("./PlasmaMX_data.txt",
                                      fill = TRUE,
                                      row.names = 1,
                                      quote=NULL, 
                                      sep="\t", 
                                      stringsAsFactors=FALSE, 
                                      header= TRUE))

#Create data matrix
metabolites <- data_meta[,-(1:10)]

#Create metadata matrix 
meta <- data_meta[,(1:10)]

#Define Group
meta$Group <- paste(meta$Cohort, meta$Genotype, sep="_")

#generate a count of metabolites present in each sample and name column "metabcount"
metab.count <- data.frame(rowSums(metabolites, na.rm = T))
colnames(metab.count) <- c("metabcount")

#Join metab.count & metadata
metab.count.meta <-cbind(meta, metab.count)

# calculate summary statistics
data_summary <- metab.count.meta %>%
  group_by(Group) %>%
  summarise(mean = mean(metabcount), 
            median = median(metabcount),
            se = (sd(metabcount)) / sqrt(length(metabcount)))

#Shapiro-Wilk normality test
  #if p < 0.05 the data is not normally distributed
  #if p > 0.05 the data is normally distributed
shapiro.test(data_summary$mean)
#W = 0.97617, p-value = 0.8792
# data is normally distributed

#2-way anova by genotype and age ("*" assumes a synergistic effect)
seq.test <- aov(metabcount  ~ Cohort*Genotype, data = metab.count.meta)
summary(seq.test)
#Results
#                 Df   Sum Sq Mean Sq F value Pr(>F)
#Cohort           1    26910   26910   0.045 0.8336  
#Genotype         1  1807432 1807432   3.012 0.0923 .
#Cohort:Genotype  1   919331  919331   1.532 0.2248  
#Residuals       32 19203055  600095                                
#Metab count group means are not significantly different

#order data
metab.count.meta <- metab.count.meta %>%
  mutate(Genotype =  factor(Genotype, levels = c("5xfAD","WT"))) %>%
  arrange(Genotype) %>%
  mutate(Cohort =  factor(Cohort, levels = c("12mo", "18mo"))) %>%
  arrange(Cohort)

#Plot
ggplot(data = metab.count.meta) +
  aes(x = Cohort, y = metabcount, fill = Genotype) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values=c('#d94801', '#2171b5')) +
  geom_point(size = 2, colour = "black", alpha = 1, 
             position = position_jitterdodge(jitter.width = 0.1),
             aes(fill=Genotype, shape=Sex)) + 
  #geom_point(position = position_jitterdodge(jitter.width = .25), alpha = 0.25) +
  #facet_wrap(Line~., nrow = 1, scales = "free_x") +
  #facet_wrap(Sex+Type+Cohort~., nrow = 1) +
  #scale_y_log10() +
  theme_bw() +
  theme_bw() +
  theme(title =  element_text(size = 10,face = "bold"))+
  theme(legend.text =  element_text(hjust = 0, size = 10,face = "bold"))+
  theme(legend.title = element_text(hjust = 0, size = 10,face = "bold"))+ 
  theme(axis.text.y = element_text(hjust = 0, size = 10, face = "bold"))+
  theme(axis.text.x = element_text(hjust = 0.5, size = 10, face = "bold"))+
  theme(axis.title.y = element_text(size = 10, face = "bold"))+
  theme(axis.title.x = element_text(size = 10, face = "bold"))+
  labs(x = NULL, y = "Summed Concentration [uM]", fill = "Genotype", 
       title = "Total Validated Metabolite Concentration")+ 
  stat_compare_means(method = "wilcox.test", label = "p.signif",
                     label.x = 1.5, vjust = 0.5, size = 4, hide.ns = F)

#ggsave("./Counts_plasma.pdf", height = 3.5, width = 4.5)

##################### Shannon Diversity of Validated Plasma Metabolites ########

#clear objects from the workspace
rm(list=ls())

library(vegan)
library(ggfortify)
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggpubr)

#import data
data_meta <- as.data.frame(read.table("./PlasmaMX_data.txt",
                                      fill = TRUE,
                                      row.names = 1,
                                      quote=NULL, 
                                      sep="\t", 
                                      stringsAsFactors=FALSE, 
                                      header= TRUE))

#Create data matrix
metabolites <- data_meta[,-(1:10)]

#Create metadata matrix 
meta <- data_meta[,(1:10)]

#Define Group
meta$Group <- paste(meta$Cohort, meta$Genotype, sep="_")

#Set NA=0
metabolites[is.na(metabolites)] <- 0

#Calculate shannon diversity
shannon <- diversity(metabolites, "shannon")

#Shapiro-Wilk normality test
  #if p < 0.05 the data is not normally distributed
  #if p > 0.05 the data is normally distributed
shapiro.test(shannon)
#data:  shannon
#W = 0.87808, p-value = 0.0009052
#the data is NOT normally distributed

#Create dataframe containing Shannon diversity and relevant metadata
shan.alpha <- data.frame(cbind(meta$Group, meta$Genotype, meta$Cohort, 
                               meta$Housing_ID, meta$Sex, shannon))

colnames(shan.alpha) <- c('Group', 'Genotype', 'Cohort', 'HousingID', 'Sex', 'shannon')
sapply(shan.alpha, is.factor)
shan.alpha$Group <- as.factor(shan.alpha$Group)
shan.alpha$Genotype <- as.factor(shan.alpha$Genotype)
shan.alpha$Cohort <- as.factor(shan.alpha$Cohort)
shan.alpha$HousingID <- as.factor(shan.alpha$HousingID)
shan.alpha$Sex <- as.factor(shan.alpha$Sex)

is.numeric(shan.alpha$shannon)
shan.alpha$shannon <- as.numeric(shan.alpha$shannon)

kruskal.test(shannon ~ Group, data = shan.alpha)
#Kruskal-Wallis chi-squared = 5.7189, df = 3, p-value = 0.1261
kruskal.test(shannon ~ Genotype, data = shan.alpha)
#Kruskal-Wallis chi-squared = 1.6827, df = 1, p-value = 0.1946
kruskal.test(shannon ~ Cohort, data = shan.alpha)
#Kruskal-Wallis chi-squared = 3.6486, df = 1, p-value = 0.05611
kruskal.test(shannon ~ HousingID, data = shan.alpha)
#Kruskal-Wallis chi-squared = 21.498, df = 16, p-value = 0.1602
kruskal.test(shannon ~ Sex, data = shan.alpha)
#Kruskal-Wallis chi-squared = 10.142, df = 1, p-value = 0.00145

# p < 0.05 some group medians are significantly different
# p > 0.05 no group medians are significantly different

#Plot
shan.alpha <- shan.alpha %>%
  mutate(Genotype =  factor(Genotype, levels = c("5xfAD", "WT"))) %>%
  arrange(Genotype) %>%
  mutate(Cohort =  factor(Cohort, levels = c("12mo", "18mo"))) %>%
  arrange(Cohort)

ggplot(data = shan.alpha) +
  aes(x = Cohort, y = shannon, fill = Genotype) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values=c('#d94801', '#2171b5')) +
  geom_point(size = 2, colour = "black", alpha = 1, 
             position = position_jitterdodge(jitter.width = 0.1),
             aes(fill=Genotype, shape=Sex)) + 
  #facet_wrap(Sex~., nrow = 1, scales = "free_x") +
  theme_bw() +
  theme(title =  element_text(size = 10,face = "bold"))+
  theme(legend.text =  element_text(hjust = 0, size = 10,face = "bold"))+
  theme(legend.title = element_text(hjust = 0, size = 10,face = "bold"))+ 
  theme(axis.text.y = element_text(hjust = 0, size = 10, face = "bold"))+
  theme(axis.text.x = element_text(hjust = .5, size = 10, face = "bold"))+
  theme(axis.title.y = element_text(size = 10, face = "bold"))+
  theme(axis.title.x = element_text(size = 10, face = "bold"))+
  labs(x = NULL, y = "Shannon Diversity", fill = "Genotype",
       title = "Shannon Diversity")+
  stat_compare_means(method = "wilcox.test", label = "p.signif",
                   label.x = 1.5, vjust = 0.5, size = 4, hide.ns = F)

#ggsave("./Shannon_plasma.pdf", height = 3.5, width = 4.5)

##################### PCA of Validated Metabolite ##############################
#protips available here: https://cran.r-project.org/web/packages/ggfortify/vignettes/plot_pca.html

#clear objects from the workspace
rm(list=ls())

library(vegan)
library(ggfortify)

#import data
data_meta <- as.data.frame(read.table("./PlasmaMX_data.txt",
                                      fill = TRUE,
                                      row.names = 1,
                                      quote=NULL, 
                                      sep="\t", 
                                      stringsAsFactors=FALSE, 
                                      header= TRUE))

#Create data matrix
metabolites <- data_meta[,-(1:10)]

#Set NA values to 0
metabolites[is.na(metabolites)] <- 0

#Create metadata matrix 
meta <- data_meta[,(1:10)]

#Define Group
meta$Group <- paste(meta$Cohort, meta$Genotype, sep="_")
#meta$Group <- paste(meta$Cohort, meta$Sex, sep="_")

#Create relative abundance matrix
metabolites_rel <- metabolites/rowSums(metabolites)
#rownames(metabolites_rel) <- meta[,1]

#PCA of raw abundance data
pca_res <- prcomp(metabolites, scale. = FALSE)

autoplot(pca_res, data = meta, colour = 'Genotype', shape = 'Cohort', size = 2,
         x = 1, y = 2,
         loadings = F, loadings.colour = 'blue',
         loadings.label = F, loadings.label.size = 3)

#ggsave("./PCA_Raw_PC1vsPC2.pdf", height = 3, width = 4.5)

#PCA of relative abundance data
pca_res_rel <- prcomp(metabolites_rel, scale. = FALSE)

autoplot(pca_res_rel, data = meta, colour = 'Genotype', shape = 'Cohort', size = 2,
         x = 1, y = 2,
         loadings = F, loadings.colour = 'blue',
         loadings.label = F, loadings.label.size = 3)

#ggsave("./PCA_Rel_PC1vsPC2_NoLoad.pdf", height = 3, width = 4.5)

##################### PERMANOVA ################################################
library(vegan)
library(ggfortify)

#clear objects from the workspace
rm(list=ls())

#import data
data_meta <- as.data.frame(read.table("./PlasmaMX_data.txt",
                                      fill = T,
                                      row.names = 1,
                                      check.names = F,
                                      quote=NULL, 
                                      sep="\t", 
                                      stringsAsFactors=F, 
                                      header= T))

#Create data matrix
metabolites <- data_meta[,-(1:10)]
#metabolites[is.na(metabolites)] <- 0


#Create metadata matrix 
meta <- data_meta[,(1:10)]

# ADONIS
meta_test <- data.frame(subset(meta, 
                               select = c(Genotype, Sex, Cohort, Housing_ID)))

meta_test$Genotype <- as.factor(meta_test$Genotype)
meta_test$Sex <- as.factor(meta_test$Sex)
meta_test$Cohort <- as.factor(meta_test$Cohort)
meta_test$Housing_ID <- as.factor(meta_test$Housing_ID)

res <- adonis(formula = metabolites ~ ., 
              data = meta_test, perm = 999, 
              method = "bray",
              na.rm = T) 
res

#            Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)  
#Genotype    1   0.03013 0.030134  2.9289 0.06945  0.025 *
#Sex         1   0.02873 0.028732  2.7927 0.06622  0.029 *
#Cohort      1   0.01072 0.010724  1.0423 0.02471  0.399  
#Housing_ID 14   0.17913 0.012795  1.2437 0.41283  0.192  
#Residuals  18   0.18519 0.010288         0.42679         
#Total      35   0.43391                  1.00000    

##################### NMDS #####################################################
library(vegan)
library(ggfortify)

#clear objects from the workspace
rm(list=ls())

#import data
data_meta <- as.data.frame(read.table("./PlasmaMX_data.txt",
                                      fill = T,
                                      row.names = 1,
                                      check.names = F,
                                      quote=NULL, 
                                      sep="\t", 
                                      stringsAsFactors=F, 
                                      header= T))

#Create data matrix
metabolites <- data_meta[,-(1:10)]

#Set NA = 0
metabolites[is.na(metabolites)] <- 0

#Create metadata matrix 
meta <- data_meta[,(1:10)]

#Define Group
meta$Group <- paste(meta$Cohort, meta$Genotype, sep="_")
 
#Get MDS stats
sol<-metaMDS(metabolites, distance = "bray", k = 2, trymax = 100)
sol
stress <- "stress = 8.206837e-05"      
#Dimensions: 2 
#Stress:     8.206837e-05 
#Stress type 1, weak ties
#No convergent solutions - best solution after 100 tries
#Scaling: centring, PC rotation, halfchange scaling 
#Species: expanded scores based on ‘wisconsin(sqrt(metabolites))’

#basic plots to test 
stressplot(sol)
plot(sol)

#Make a new data frame, and put Group and Factor information there, to be 
#useful for coloring, and shape of points
NMDS=data.frame(x=sol$point[,1],y=sol$point[,2],Genotype=as.factor(meta$Genotype),
                Cohort=as.factor(meta$Cohort))

#order
NMDS$Cohort <- factor(NMDS$Cohort, levels = c("18mo", "12mo","4mo"))
NMDS$Genotype <- factor(NMDS$Genotype, levels = c("WT","5xfAD")) 

shape_values<-c(21,24,22,3)
gwsC <- c('#d94801', '#2171b5', 'grey') 

ggplot(data=NMDS,aes(x,y, colour=Genotype))+
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 14,face = "bold"))+ 
  scale_x_continuous(name = "NMDS1") + scale_y_continuous(name = "NMDS2")+
  #annotate("text", x=0.25, y=-0.1, size = 4, label=stress)+ #add stress indicators
  geom_point(size=3, aes(shape=Cohort, fill=Genotype))+
  stat_ellipse(linetype = 2, aes(group = Cohort), show.legend = F, color = "black") +
  scale_shape_manual(values=shape_values) + 
  scale_color_manual(values=gwsC) + scale_fill_manual(values=gwsC)+
  theme(axis.text = element_text(size = 10,face = "bold"),
        axis.title = element_text(size = 10,face = "bold"))+
  theme(axis.text.x=element_text( size = 10, hjust = 0.5))

#ggsave("./Outputs/NMDS_All_Ellipses.pdf", height = 3.5, width = 4.5)

#Warning messages:
#  1: In MASS::cov.trob(data[, vars]) : Probable convergence failure
#2: In MASS::cov.trob(data[, vars]) : Probable convergence failure

##################### PCOA #####################################################

#clear objects from the workspace
rm(list=ls())

#Load Libraries
library(vegan)
library(ggfortify)
library(dplyr)
library(randomForest)
library(tibble)

#import data
data_meta <- as.data.frame(read.table("./PlasmaMX_data.txt",
                                      fill = TRUE,
                                      row.names = 1,
                                      quote=NULL, 
                                      sep="\t", 
                                      stringsAsFactors=FALSE, 
                                      header= TRUE))

#Create data matrix
metabolites <- data_meta[,-(1:10)]

#Set NA values to 0
metabolites[is.na(metabolites)] <- 0

#Create metadata matrix 
meta <- data_meta[,(1:10)]
meta$Group <- paste(meta$Cohort, meta$Genotype, sep="_")

#Create a distance matrix
bray = vegdist(metabolites, method = "bray", na.rm = T) 

#Here, the number of dimensions (k) is equal to n-1, where n is the number of OTUS
pcoa = cmdscale(bray, eig = T, k = nrow(metabolites[,1:ncol(metabolites)])-1, add = T) 

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
#env is the abundance matrix. env is the metabolites excluding metadata

#calculate: This step takes a long time. 
test2 = envfit(ord = pcoa, env = metabolites, na.rm = TRUE) 

#make the loadings
test3 = as.data.frame(scores(test2, display = "vectors")) 
test3 = test3 %>% rownames_to_column()

#order data
points$Cohort <- factor(points$Cohort, levels = c("18mo", "12mo")) 
points$Genotype <- factor(points$Genotype, levels = c("5xfAD","WT")) 

#define color scheme
gwsC <- c('#d94801', '#2171b5')

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

  #ggsave("./P180_PCOA_5xfAD_Valid_Loadings&Ellipses.pdf", height = 3.5, width = 4.5)

#Plot with ellipses only
ggplot(data = points) +
  aes(x = V1, y = V2) +
  theme_bw() +
  #geom_segment(data = test3, aes(x = 0, xend = Dim1, y = 0, yend = Dim2), arrow = arrow(length = unit(0.25, "cm")), colour = "grey90") +  #This line draws the loadings
  #geom_text(data = test3, aes(x = Dim1, y = Dim2, label = rowname), size = 4, alpha = .3, inherit.aes = FALSE) + 
  geom_point(aes(pch = Cohort, fill = Genotype), size = 3, alpha = 1) + #sets shape
  scale_shape_manual(values = c(21,22,23,24,3,4,8,9), name = "Cohort") + # from https://ggplot2.tidyverse.org/reference/scale_shape.html
  guides(color = F, fill = guide_legend(override.aes = list(shape = 21))) +
  stat_ellipse(linetype = 2, aes(group = Group), show.legend = F, color = "black", level = 0.95) +
  scale_fill_manual(values = gwsC) +
  labs(color = "Sample type", x = bquote("PCo1 ("~.(round(pcoa_var[1]*100, digits = 1))~"%) "), 
       y = bquote("PCo2 ("~.(round(pcoa_var[2]*100, digits = 1))~"%) "))


#ggsave("./P180_PCOA_5xfAD_Valid_Ellipses.pdf", height = 3.5, width = 4.5)

#3D plot
library(plotly)
plot_ly(x=points$V1, y=points$V2, z=points$V3, type="scatter3d", mode="markers", 
        color=points$Group)

##################### Random Forest ############################################

#clear objects from the workspace
rm(list=ls())

#Load Libraries
library(rfPermute) #https://cran.r-project.org/web/packages/rfPermute/rfPermute.pdf
library(vegan)
library(ggfortify)
library(dplyr)
library(tidyr)
library(tibble)

#import data
data_meta <- as.data.frame(read.table("./PlasmaMX_data.txt",
                                      fill = TRUE,
                                      row.names = 1,
                                      quote=NULL, 
                                      sep="\t", 
                                      stringsAsFactors=FALSE, 
                                      header= TRUE))

#Select for subsets
#data_meta <- filter(data_meta, Cohort == "12mo")

#Extract metadata
meta <- data_meta[,(1:10)]

#Define group (e.g., Cohort_Line_SampleType)
#meta$Group <- paste(meta$Cohort, meta$Genotype, sep="_")
#meta$Group <- paste(meta$Cohort)
meta$Group <- paste(meta$Genotype)
#meta$Group <- paste(meta$Sex)

#Extract data
metabolites <- data_meta[,-(1:10)]

#RF doesnt like having more features than samples, so remove low abundance mx
#metabolites <- metabolites[,colMeans(metabolites) >= 1000]

# Random forest only uses one metadata factor
metabolites <- merge(meta, metabolites, by = "row.names")

# Re-establish rownames
metabolites <- column_to_rownames(metabolites, var = "Row.names")

# Remove unneeded metadata
metabolites <- metabolites[,-(1:10)]

# Save and re-create OTU table to make numeric
write.table(metabolites, "./temp.txt", sep = "\t", row.names = TRUE, quote=FALSE)
OTU_num <- as.data.frame(read.table("./temp.txt",
                                    row.names =1, 
                                    quote=NULL,
                                    check.names=FALSE,
                                    sep="\t", 
                                    stringsAsFactors=TRUE, 
                                    header= TRUE))
OTU_num[is.na(OTU_num)] <- 0 #set NA=0

# Perform Random Forest
rf_out <- rfPermute(formula = Group ~ ., 
                    data = OTU_num, proximity = T, importance = F, 
                    ntree = 601, num.cores = 32)

rf_out$Group <- factor(rf_out$Group, levels = c("18mo", "12mo")) 

# View Results
# Plot Confusion Matrix heatmap
confplot <- plotConfMat(rf_out, title= NULL, plot = TRUE)
#ggsave("./RF_5xfAD_18&12mo_Cohort_ConfPlot.pdf", height = 3.5, width = 4.5)

# Create a plot of Random Forest proximity scores using multi-dimensional scaling. 
#proxplot <- proximityPlot(rf_out, legend.loc = "bottom",   dim.x = 1, dim.y = 2, hull.alpha = NULL)
proxplot <- proximityPlot(rf_out, legend.loc = "bottom",   dim.x = 1, dim.y = 2)

safe_colorblind_palette <- c("#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
                             "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")

proximityPlot( rf_out, dim.x = 1, dim.y = 2,
               class.cols = NULL , #NULL or c('#d94801', '#2171b5') or safe_colorblind_palette
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

#ggsave("./RF_5xfAD_12_Genotype_ProxPlot.pdf", height = 3.5, width = 4.5)

# Dotchart of variable importance as measured by a Random Forest
varplot <- varImpPlot(rf_out, type = 1, n.var = 20)

#pdf("./RF_5xfAD_12_Genotype_VarPlot.pdf", height = 4.5, width = 4.5)  
#varImpPlot(rf_out, type = 1, n.var = 20, main = NULL, cex = 0.8)
#dev.off()

##################### Linear Mixed Effects Model ############################################

#clear objects from the workspace
rm(list=ls())

#Load Libraries

library(vegan)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(lmerTest)
library(reshape2)

#Read in your OTU table. 
data_meta <- as.data.frame(read.table("./PlasmaMX_data.txt",
                                      fill = TRUE,
                                      row.names = 1,
                                      quote=NULL, 
                                      sep="\t", 
                                      stringsAsFactors=FALSE, 
                                      header= TRUE))

#Select for subsets
#data_meta <- filter(data_meta, Cohort == "12mo")
#data_meta <- data_meta[grep("^18mo|^12mo",data_meta$Cohort),]
#data_meta <- filter(data_meta, Sex == "Male")

#Select for data
data <- data_meta[,-(1:10)]

#Set NA = 0
data[is.na(data)] <- 0

data <- t(data)

#Make relative abundances. You can ignore the complete.cases part if you want.
#data<- data[complete.cases(data), ]
column_sums <- colSums(data, na.rm = TRUE)
data_norm <- apply(data, 1, '/', column_sums)

#Read in your metadata.
var<-data_meta[,(1:10)]

attach(var)

#Transpose the OTU table and filter out low abundance mx 
data_norm <- t(data_norm)
data_norm <- as.data.frame(data_norm)
data_norm$mean <- rowMeans(data_norm)
data_norm <- data_norm %>%
  rownames_to_column("Sample") %>%
  #filter(data_norm$mean > 0.00001) %>%
  subset(., select = -c(mean))
data_norm = setNames(data.frame(t(data_norm[,-1])),data_norm[,1])

#Make a longform data table with melt.
melted<-melt(data_norm)

#Run a Linear Mixed Effects Model (LMER). 
#It's ~lmer(OTU abundance ~ Metadata category for comparison + (Random effect)
lmertest<-melted %>%
  split(.$variable) %>%
  map(~lmer(value ~ Genotype + (1|as.factor(Housing_ID)), data = .)) %>%
  map(summary) %>% 
  map(coef)

lmertest[["Serotonin"]]
lmertest[["Carnosine"]]

#Pull out pvalues from your test. The "10" is just the 10th position in the output.
pvalues<-map_dbl(lmertest,10)
p_ordered<-(sort(pvalues))
padj<-p.adjust(p_ordered, method = "fdr")
padj
padj_df<-as.data.frame(padj)

#Filter out insignificant taxa.
padj_df<-padj_df %>%
  filter(padj < 0.05)
padj_df<-rownames_to_column(padj_df, "names")

#Display results
padj_df

#Genotype: 12 & 18 mo males & females
# names        padj
#1 Carnosine 0.008905186

#Genotype: 18 mo males & females
#names       padj
#lysoPC.a.C18.1 0.02604629
#2      Carnosine 0.02604629

#Sex: 12 & 18 mo males & females
# names        padj
#1       SM.C20.2 0.004850558
#2    PC.ae.C42.0 0.009234041
#3    PC.aa.C34.3 0.016231419
#4    PC.aa.C42.6 0.016231419
#5            Lys 0.016231419
#6  SM..OH..C24.1 0.022236056
#7    PC.aa.C36.4 0.027860312
#8    PC.aa.C36.5 0.027860312
#9    PC.ae.C38.0 0.027860312
#10   PC.aa.C36.3 0.029567582
#11   PC.aa.C38.6 0.044784415

#Age: 12 & 18 mo males & females
# names        padjnames        padj
#1  Serotonin 0.001529308
#2   Spermine 0.011828063
#3 Spermidine 0.048890457
#4        Glu 0.048890457

#Genotype: 12 mo males & females
#names       padj
#0 rows> (or 0-length row.names)

#Genotype: 12 & 18 mo males
#names padj 
#0 rows> (or 0-length row.names)

#Genotype: 12 & 18 mo females
#      names       padj
#0 rows> (or 0-length row.names)

detach(var)

##################### Concentration Plots for Metabolites with Sex Diff ########
#clear objects from the workspace
rm(list=ls())

library("tidyverse")
library("readxl")
library("tibble")
library("ggpubr")
library("ggplot2")

#Read in OTU table
data_meta <- as.data.frame(read.table("./PlasmaMX_data.txt",
                                      fill = T,
                                      row.names = 1,
                                      check.names = F,
                                      quote=NULL, 
                                      sep="\t", 
                                      stringsAsFactors=F, 
                                      header= T))

#Order by age
data_meta$Cohort = factor(data_meta$Cohort, levels=c('12mo', '18mo'))

#=====mx Plots=======#
#testplot
MOI = "Lys"
data_meta[[MOI]]
data_meta[[MOI]] <- as.numeric(data_meta[[MOI]])
ggplot(data = data_meta) + 
  aes(x = Cohort, y = .data[[MOI]], fill = Sex, label = Mouse_ID) +
  geom_boxplot(outlier.shape = NA) +
  #scale_fill_manual(values=c('#d94801', '#2171b5')) +
  scale_fill_manual(values=c('blue', 'forestgreen')) +
  geom_point(size = 2, colour = "black", alpha = 1, 
             position = position_jitterdodge(jitter.width = 0.1),
             aes(fill=Sex, shape=Genotype)) + 
  #facet_wrap(Sex~., nrow = 1) +
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
  labs(x = "Age", y = "Plasma Concentration [nM]", fill = "Sex", 
       title = MOI) +
  stat_compare_means(method = "wilcox.test", label = "p.signif",
                     label.x = 1.5, vjust = 0.5, size = 4, hide.ns = F)

#ggsave(path = "./", filename=paste(MOI,".pdf",sep=""), 
#       device = "pdf", height = 1.75, width = 1.5)

#Create a list of selected analytes for plotting. Could also be uploaded from a file.
selects = as.data.frame(c("Carnosine", "lysoPC a C18:1", "Serotonin", "Spermine", "Spermidine", "Glu", 
                          "SM C20:2", "PC ae C42:0", "PC aa C34:3", "PC aa C42:6", 
                          "Lys", "SM (OH) C24:1", "PC aa C36:4","PC aa C36:5", "PC ae C38:0", 
                          "PC aa C36:3", "PC aa C38:6"))
colnames(selects) <- "moi"

#plot all selects using for loop and save as pdf 
for (i in 1:nrow(selects)) {
  MOI = selects[i,1]
  data_meta[[MOI]] <- as.numeric(data_meta[[MOI]])
  mxplot <- ggplot(data = data_meta) + 
    aes(x = Cohort, y = .data[[MOI]], fill = Sex, label = Mouse_ID) +
    geom_boxplot(outlier.shape = NA) +
    #scale_fill_manual(values=c('#d94801', '#2171b5')) +
    scale_fill_manual(values=c('blue', 'forestgreen')) +
    geom_point(size = 2, colour = "black", alpha = 1, 
               position = position_jitterdodge(jitter.width = 0.1),
               aes(fill=Sex, shape=Genotype)) + 
    #facet_wrap(Sex~., nrow = 1) +
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
    labs(x = "Age", y = "Plasma Concentration [nM]", fill = "Sex", 
         title = MOI)+
    stat_compare_means(method = "wilcox.test", label = "p.signif",
                       label.x = 1.5, vjust = 0.5, size = 4, hide.ns = F)
  
  mxplot
  ggsave(mxplot, path = "./", filename=paste(MOI,".pdf",sep=""), 
         device = "pdf", height = 1.75, width = 1.5)
}

##################### Concentration Plots for Metabolites with Age Diff ########
library("tidyverse")
library("readxl")
library("tibble")
library("ggpubr")
library("ggplot2")

#clear objects from the workspace
rm(list=ls())

#Read in OTU table
data_meta <- as.data.frame(read.table("./PlasmaMX_data.txt",
                                      fill = T,
                                      row.names = 1,
                                      check.names = F,
                                      quote=NULL, 
                                      sep="\t", 
                                      stringsAsFactors=F, 
                                      header= T))

#Order by age
data_meta$Cohort = factor(data_meta$Cohort, levels=c('12mo', '18mo'))

#testplot
MOI = "Glu"
data_meta[[MOI]]
data_meta[[MOI]] <- as.numeric(data_meta[[MOI]])
ggplot(data = data_meta) + 
  aes(x = Cohort, y = .data[[MOI]], fill = Genotype, label = Mouse_ID) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values=c('#d94801', '#2171b5')) +
  geom_point(size = 2, colour = "black", alpha = 1, 
             position = position_jitterdodge(jitter.width = 0.1),
             aes(fill=Genotype, shape=Sex)) + 
  #facet_wrap(Sex~., nrow = 1) +
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
  labs(x = "Age", y = "Plasma Concentration [uM]", fill = "Genotype", 
       title = MOI)
ggsave(path = "./", filename=paste(MOI,".pdf",sep=""), 
       device = "pdf", height = 1.75, width = 1.5)

#Create a list of selected analytes for plotting. Could also be uploaded from a file.
#selects = as.data.frame(c("PC ae C38:4", "Pro", "Lys", "Serotonin", "Taurine", "DOPA", "Phe", "C7-DC"))
#selects = as.data.frame(c("PC aa C36:6", "PC aa C38:0", "PC aa C38:6", "PC aa C40:1", 
#                          "PC aa C40:2", "PC aa C40:6", "PC ae C40:6", "lysoPC a C18:2", 
#                          "C3", "C16:1-OH"))
#selects = as.data.frame(c("PC ae C40:6", "PC aa C40:1", "PC aa C38:6", "PC aa C38:0", 
#                          "PC aa C36:6", "lysoPC a C18:2", "C3", "PC ae C36:4", "C10:2", 
#                          "C9", "PC ae C42:1", "PC aa C38:3", "C5", "ADMA", 
#                          "Asn", "PC aa C34:4", "C18:1-OH", "PC ae C34:0", 
#                          "C5-OH (C3-DC-M)", "PC aa C40:5", "PC aa C32:0", "C16:2", "C12:1", "C10:1"))
selects = as.data.frame(c("Carnosine", "lysoPC a C18:1", "Serotonin", "Spermine", "Spermidine", "Glu", 
                          "SM C20:2", "PC ae C42:0", "PC aa C34:3", "PC aa C42:6", 
                          "Lys", "SM (OH) C24:1", "PC aa C36:4","PC aa C36:5", "PC ae C38:0", 
                          "PC aa C36:3", "PC aa C38:6"))
colnames(selects) <- "moi"

#plot all selects using for loop and save to file. 
for (i in 1:nrow(selects)) {
  MOI = selects[i,1]
  data_meta[[MOI]] <- as.numeric(data_meta[[MOI]])
  mxplot <- ggplot(data = data_meta) + 
    aes(x = Cohort, y = .data[[MOI]], fill = Genotype, label = Mouse_ID, na.rm = TRUE) + 
    geom_boxplot(outlier.shape = NA) +
    scale_fill_manual(values=c('#d94801', '#2171b5')) +
    geom_point(size = 2, colour = "black", alpha = 1, 
               position = position_jitterdodge(jitter.width = 0.1),
               aes(fill=Genotype, shape=Sex)) + 
    #facet_wrap(Sex~., nrow = 1) +
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
    labs(x = "Age", y = "Plasma Concentration [uM]", fill = "Genotype", 
         title = MOI)+
    stat_compare_means(method = "wilcox.test", label = "p.signif",
                       label.x = 1.5, vjust = 0.5, size = 4, hide.ns = F)
  mxplot
  ggsave(mxplot, path = "./", filename=paste(MOI,".pdf",sep=""), 
         device = "pdf", height = 1.75, width = 1.5)
}

#plot all selects using for loop and save to file. 
# x-axis is geno
for (i in 1:nrow(selects)) {
  MOI = selects[i,1]
  data_meta[[MOI]] <- as.numeric(data_meta[[MOI]])
  mxplot <- ggplot(data = data_meta) + 
    aes(x = Genotype, y = .data[[MOI]], fill = Cohort, label = Mouse_ID, na.rm = TRUE) + 
    geom_boxplot(outlier.shape = NA) +
    #scale_fill_manual(values=c('#d94801', '#2171b5')) +
    geom_point(size = 2, colour = "black", alpha = 1, 
               position = position_jitterdodge(jitter.width = 0.1),
               aes(fill=Cohort, shape=Sex)) + 
    #facet_wrap(Sex~., nrow = 1) +
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
    labs(x = "", y = "Plasma Concentration [uM]", fill = "Cohort", 
         title = MOI)+
    stat_compare_means(method = "wilcox.test", label = "p.signif",
                       label.x = 1.5, vjust = 0.5, size = 4, hide.ns = F)
  mxplot
  ggsave(mxplot, path = "./", filename=paste(MOI,".pdf",sep=""), 
         device = "pdf", height = 1.75, width = 1.5)
}

##################### Concentration Plots for Metabolites with Genotype Diff #####
library("tidyverse")
library("readxl")
library("tibble")
library("ggpubr")
library("ggplot2")

#clear objects from the workspace
rm(list=ls())

#Read in OTU table
data_meta <- as.data.frame(read.table("./PlasmaMX_data.txt",
                                      fill = T,
                                      row.names = 1,
                                      check.names = F,
                                      quote=NULL, 
                                      sep="\t", 
                                      stringsAsFactors=F, 
                                      header= T))

#Order by age
data_meta$Cohort = factor(data_meta$Cohort, levels=c('12mo', '18mo'))

#testplot
MOI = "Serotonin"
data_meta[[MOI]]
data_meta[[MOI]] <- as.numeric(data_meta[[MOI]])
ggplot(data = data_meta) + 
  aes(x = Cohort, y = .data[[MOI]], fill = Genotype, label = Mouse_ID) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values=c('#d94801', '#2171b5')) +
  geom_point(size = 2, colour = "black", alpha = 1, 
             position = position_jitterdodge(jitter.width = 0.1),
             aes(fill=Genotype, shape=Sex)) + 
  #facet_wrap(Sex~., nrow = 1) +
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
  labs(x = "Age", y = "Plasma Concentration [nM]", fill = "Genotype", 
       title = MOI)
#ggsave(path = "./", filename=paste(MOI,".pdf",sep=""), 
#       device = "pdf", height = 1.75, width = 1.5)

#Create a list of selected analytes for plotting. Could also be uploaded from a file.
#selects = as.data.frame(c("PC ae C38:4", "Pro", "Lys", "Serotonin", "Taurine", "DOPA", "Phe", "C7-DC"))
#selects = as.data.frame(c("PC aa C36:6", "PC aa C38:0", "PC aa C38:6", "PC aa C40:1", 
#                          "PC aa C40:2", "PC aa C40:6", "PC ae C40:6", "lysoPC a C18:2", 
#                          "C3", "C16:1-OH"))
#selects = as.data.frame(c("PC ae C40:6", "PC aa C40:1", "PC aa C38:6", "PC aa C38:0", 
#                          "PC aa C36:6", "lysoPC a C18:2", "C3", "PC ae C36:4", "C10:2", 
#                          "C9", "PC ae C42:1", "PC aa C38:3", "C5", "ADMA", 
#                          "Asn", "PC aa C34:4", "C18:1-OH", "PC ae C34:0", 
#                          "C5-OH (C3-DC-M)", "PC aa C40:5", "PC aa C32:0", "C16:2", "C12:1", "C10:1"))
selects = as.data.frame(c("Carnosine", "lysoPC a C18:1", "Serotonin", "Spermine", "Spermidine", "Glu", 
                          "SM C20:2", "PC ae C42:0", "PC aa C34:3", "PC aa C42:6", 
                          "Lys", "SM (OH) C24:1", "PC aa C36:4","PC aa C36:5", "PC ae C38:0", 
                          "PC aa C36:3", "PC aa C38:6"))
colnames(selects) <- "moi"

#plot all selects using for loop and save to file. 
#"selects" is a data frame containing the list of analytes in column 1 

#ggplot(data = subset(data_meta_mod,!is.na(data_meta_mod[,]))) +

for (i in 1:nrow(selects)) {
  MOI = selects[i,1]
  data_meta[[MOI]] <- as.numeric(data_meta[[MOI]])
  mxplot <- ggplot(data = data_meta) + 
    aes(x = Cohort, y = .data[[MOI]], fill = Genotype, label = Mouse_ID, na.rm = TRUE) + 
    geom_boxplot(outlier.shape = NA) +
    scale_fill_manual(values=c('#d94801', '#2171b5')) +
    geom_point(size = 2, colour = "black", alpha = 1, 
               position = position_jitterdodge(jitter.width = 0.1),
               aes(fill=Genotype, shape=Sex)) + 
    #facet_wrap(Sex~., nrow = 1) +
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
    labs(x = "Age", y = "Plasma Concentration [nM]", fill = "Genotype", 
         title = MOI)+
    stat_compare_means(method = "wilcox.test", label = "p.signif",
                       label.x = 1.5, vjust = 0.5, size = 4, hide.ns = F)
  
  mxplot
  ggsave(mxplot, path = "./", filename=paste(MOI,".pdf",sep=""), 
         device = "pdf", height = 1.75, width = 1.5)
}
