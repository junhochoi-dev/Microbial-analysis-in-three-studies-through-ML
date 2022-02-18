######################################
##                                  ##    
##  Programmed by Junho Choi        ##
##  2022.01.24                      ##
##  CI_Study                        ##
##                                  ##
######################################

## loading packages
setRepositories(ind = 1:7)

library(data.table)
library(dplyr)
library(Rtsne)
library(cluster)
library(NbClust)
library(factoextra)
library(doParallel)
library(foreach)
library(lattice)
library(devtools)
library(dplyr)
library(rgl)

library(ggplot2)
library(ggbiplot)
library(ggpubr)

library(caret)
library(plyr)
library(party)
library(partykit)
library(mboost)

library(patchwork)
library(gridExtra)
library(doMC)

library(edgeR)
library(limma)

setwd("/home/bicjh/CI_final")
getwd()

##################################################################################################################################
# MGYS00005138
# PRJNA275918
# Infant nasopharyngeal bacteria microbiome

# MGYS00005184
# PRJEB14839
# Association of host genome with intestinal microbial composition in a large healthy cohort

# MGYS00005194
# PRJEB30316
# the upper respiratory tract microbiome of chronic rhinosinusitis patients

# 3 DISEASE : Healthy(건강) / Acute Respiratory Infection(급성 호흡기 전염) / Chronic rhinosinusitis(만성축농증)

##################################################################################################################################
emptyValue <- c("", " ", "not collected", "unspecified")

# disease study 1 : Acute Respiratory Infection Data
MGYS00005138_metadata <- fread(file = '/home/bicjh/CI_final/MGYS00005138/MGYS00005138_metadata.csv', na.strings = c(emptyValue, NA))
MGYS00005138_srarun_table <- fread(file = '/home/bicjh/CI_final/MGYS00005138/MGYS00005138_SraRunTable.txt', na.strings = c(emptyValue, NA))
MGYS00005138_otu_table <- fread('/home/bicjh/CI_final/MGYS00005138/MGYS00005138_SSU_OTUtable.csv')

for(i in 1:(ncol(MGYS00005138_otu_table)-5)){
  idx <- which(MGYS00005138_metadata$analysis_accession %in% colnames(MGYS00005138_otu_table)[i+5])
  colnames(MGYS00005138_otu_table)[i+5] <- MGYS00005138_metadata$run_accession[idx]
}

# disease study 2 : Healthy Data
MGYS00005184_metadata <- fread(file = '/home/bicjh/CI_final/MGYS00005184/MGYS00005184_metadata.csv', na.strings = c(emptyValue, NA))
MGYS00005184_srarun_table <- fread(file = '/home/bicjh/CI_final/MGYS00005184/MGYS00005184_SraRunTable.txt', na.strings = c(emptyValue, NA))
MGYS00005184_otu_table <- fread('/home/bicjh/CI_final/MGYS00005184/MGYS00005184_OTU_Table.tsv')

# disease study 3 : Chronic rhinosinusitis Data
MGYS00005194_metadata <- fread(file = '/home/bicjh/CI_final/MGYS00005194/MGYS00005194_metadata.csv', na.strings = c(emptyValue, NA))
MGYS00005194_srarun_table <- fread(file = '/home/bicjh/CI_final/MGYS00005194/MGYS00005194_SraRunTable.txt', na.strings = c(emptyValue, NA))
MGYS00005194_otu_table <- fread('/home/bicjh/CI_final/MGYS00005194/MGYS00005194_OTU_Table.tsv')

##################################################################################################################################
### Data preprocessing : similiar column

# MGYS00005138
MGYS00005138_metadata <- MGYS00005138_metadata[, -c('sample_collection date', 'sample_environment (biome)', 'sample_environment (feature)', 'sample_environment (material)')]
MGYS00005138_metadata <- MGYS00005138_metadata[, -c('sample_geographic location (latitude)', 'sample_geographic location (longitude)')]
MGYS00005138_metadata <- rename(MGYS00005138_metadata, 'sample_collection_date' = 'sample_collection-date')
MGYS00005138_metadata <- rename(MGYS00005138_metadata, 'sample_environment_feature' = 'sample_environment-feature')
MGYS00005138_metadata <- rename(MGYS00005138_metadata, 'sample_environment_material' = 'sample_environment-material')
MGYS00005138_srarun_table <- MGYS00005138_srarun_table[, -c('BioSample')]
MGYS00005138_srarun_table <- MGYS00005138_srarun_table[, -c('env_biome')]
MGYS00005138_srarun_table <- MGYS00005138_srarun_table[, -c('lat_lon', 'latitude_and_longitude')]
MGYS00005138_srarun_table <- rename(MGYS00005138_srarun_table, 'run_accession' = 'Run')

# MGYS00005184
MGYS00005184_metadata <- MGYS00005184_metadata[, -c('sample_collection date', 'sample_environment (material)', 'sample_host sex')]
MGYS00005184_metadata <- MGYS00005184_metadata[, -c('sample_geographic location (latitude)', 'sample_geographic location (longitude)')]
MGYS00005184_metadata <- rename(MGYS00005184_metadata, 'sample_collection_date' = 'sample_collection-date')
MGYS00005184_metadata <- rename(MGYS00005184_metadata, 'sample_environment_material' = 'sample_environment-material')
MGYS00005184_srarun_table <- MGYS00005184_srarun_table[, -c('BioSample')]
MGYS00005184_srarun_table <- MGYS00005184_srarun_table[, -c('ena_checklist')]
MGYS00005184_srarun_table <- MGYS00005184_srarun_table[, -c('lat_lon', 'latitude_and_longitude')]
MGYS00005184_srarun_table <- MGYS00005184_srarun_table[, -c('geographic_location_(latitude)', 'geographic_location_(longitude)')]
MGYS00005184_srarun_table <- MGYS00005184_srarun_table[, -c('Sample Name')]
MGYS00005184_srarun_table <- rename(MGYS00005184_srarun_table, "run_accession" = "Run")
# MGYS00005194
MGYS00005194_srarun_table <- MGYS00005194_srarun_table[, -c('BioSample')]
MGYS00005194_srarun_table <- MGYS00005194_srarun_table[, -c('ena_checklist')]
MGYS00005194_srarun_table <- MGYS00005194_srarun_table[, -c('Sample Name')]
MGYS00005194_srarun_table <- rename(MGYS00005194_srarun_table, "run_accession" = "Run")

##################################################################################################################################
#
MGYS00005138_data <- merge(MGYS00005138_metadata, MGYS00005138_srarun_table, by='run_accession', all = TRUE)
#
MGYS00005184_data <- merge(MGYS00005184_metadata, MGYS00005184_srarun_table, by='run_accession', all = TRUE)
#
MGYS00005194_data <- merge(MGYS00005194_metadata, MGYS00005194_srarun_table, by='run_accession', all = TRUE)
#
MGYS00005138_data[which(MGYS00005138_data$host_Age=="Infant")]$host_Age <- 1
MGYS00005138_data$host_Age <- as.integer(MGYS00005138_data$host_Age)
MGYS00005138_data$collection_date <- as.character(MGYS00005138_data$collection_date)
#
MGYS00005184_data$`Library Name` <- as.integer(MGYS00005184_data$`Library Name`)
MGYS00005184_data$collection_date <- as.character(MGYS00005184_data$collection_date)
#
MGYS00005138_data$Host_disease
MGYS00005184_data$Host_disease
MGYS00005194_data$Host_disease
MGYS00005184_data <- cbind(MGYS00005184_data, c(rep('Healthy', nrow(MGYS00005184_data))))
MGYS00005184_data <- rename(MGYS00005184_data, "Host_disease" = "V2")
MGYS00005194_data <- cbind(MGYS00005194_data, c(rep('Chronic rhinosinusitis', nrow(MGYS00005194_data))))
MGYS00005194_data <- rename(MGYS00005194_data, "Host_disease" = "V2")
MGYS00005138_data$Host_disease
MGYS00005184_data$Host_disease
MGYS00005194_data$Host_disease

##################################################################################################################################

MGYS00005138_Host_disease_null_data <- MGYS00005138_data[is.na(MGYS00005138_data$Host_disease)]
MGYS00005138_data <- MGYS00005138_data[ complete.cases(MGYS00005138_data[ , c("Host_disease")]), ]

##################################################################################################################################

metadata <- merge(MGYS00005138_data, MGYS00005184_data, by=intersect(names(MGYS00005138_data), names(MGYS00005184_data)), all=T)
metadata <- merge(metadata, MGYS00005194_data, by=intersect(names(metadata), names(MGYS00005194_data)), all=T)

##################################################################################################################################

# OTU Cleansing MGYS00005138
MGYS00005138_phylum <- MGYS00005138_otu_table[which(MGYS00005138_otu_table$V1 != 'p__' & MGYS00005138_otu_table$V1 != ''), -c(2,3,4,5)]
MGYS00005138_phylum <- aggregate(MGYS00005138_phylum[, -1], by=list(MGYS00005138_phylum$V1), FUN=sum)
MGYS00005138_class <- MGYS00005138_otu_table[which(MGYS00005138_otu_table$V2 != 'c__' & MGYS00005138_otu_table$V2 != ''), -c(1,3,4,5)]
MGYS00005138_class <- aggregate(MGYS00005138_class[, -1], by=list(MGYS00005138_class$V2), FUN=sum)
MGYS00005138_order <- MGYS00005138_otu_table[which(MGYS00005138_otu_table$V3 != 'o__' & MGYS00005138_otu_table$V3 != ''), -c(1,2,4,5)]
MGYS00005138_order <- aggregate(MGYS00005138_order[, -1], by=list(MGYS00005138_order$V3), FUN=sum)
MGYS00005138_family <- MGYS00005138_otu_table[which(MGYS00005138_otu_table$V4 != 'f__' & MGYS00005138_otu_table$V4 != ''), -c(1,2,3,5)]
MGYS00005138_family <- aggregate(MGYS00005138_family[, -1], by=list(MGYS00005138_family$V4), FUN=sum)
MGYS00005138_genus <- MGYS00005138_otu_table[which(MGYS00005138_otu_table$V5 != 'g__' & MGYS00005138_otu_table$V5 != ''), -c(1,2,3,4)]
MGYS00005138_genus <- aggregate(MGYS00005138_genus[, -1], by=list(MGYS00005138_genus$V5), FUN=sum)

# OTU Cleansing MGYS00005184
MGYS00005184_phylum <- MGYS00005184_otu_table[which(MGYS00005184_otu_table$V1 != 'p__' & MGYS00005184_otu_table$V1 != ''), -c(2,3,4,5)]
MGYS00005184_phylum <- aggregate(MGYS00005184_phylum[, -1], by=list(MGYS00005184_phylum$V1), FUN=sum)
MGYS00005184_class <- MGYS00005184_otu_table[which(MGYS00005184_otu_table$V2 != 'c__' & MGYS00005184_otu_table$V2 != ''), -c(1,3,4,5)]
MGYS00005184_class <- aggregate(MGYS00005184_class[, -1], by=list(MGYS00005184_class$V2), FUN=sum)
MGYS00005184_order <- MGYS00005184_otu_table[which(MGYS00005184_otu_table$V3 != 'o__' & MGYS00005184_otu_table$V3 != ''), -c(1,2,4,5)]
MGYS00005184_order <- aggregate(MGYS00005184_order[, -1], by=list(MGYS00005184_order$V3), FUN=sum)
MGYS00005184_family <- MGYS00005184_otu_table[which(MGYS00005184_otu_table$V4 != 'f__' & MGYS00005184_otu_table$V4 != ''), -c(1,2,3,5)]
MGYS00005184_family <- aggregate(MGYS00005184_family[, -1], by=list(MGYS00005184_family$V4), FUN=sum)
MGYS00005184_genus <- MGYS00005184_otu_table[which(MGYS00005184_otu_table$V5 != 'g__' & MGYS00005184_otu_table$V5 != ''), -c(1,2,3,4)]
MGYS00005184_genus <- aggregate(MGYS00005184_genus[, -1], by=list(MGYS00005184_genus$V5), FUN=sum)

# OTU Cleansing MGYS00005194
MGYS00005194_phylum <- MGYS00005194_otu_table[which(MGYS00005194_otu_table$V1 != 'p__' & MGYS00005194_otu_table$V1 != ''), -c(2,3,4,5)]
MGYS00005194_phylum <- aggregate(MGYS00005194_phylum[, -1], by=list(MGYS00005194_phylum$V1), FUN=sum)
MGYS00005194_class <- MGYS00005194_otu_table[which(MGYS00005194_otu_table$V2 != 'c__' & MGYS00005194_otu_table$V2 != ''), -c(1,3,4,5)]
MGYS00005194_class <- aggregate(MGYS00005194_class[, -1], by=list(MGYS00005194_class$V2), FUN=sum)
MGYS00005194_order <- MGYS00005194_otu_table[which(MGYS00005194_otu_table$V3 != 'o__' & MGYS00005194_otu_table$V3 != ''), -c(1,2,4,5)]
MGYS00005194_order <- aggregate(MGYS00005194_order[, -1], by=list(MGYS00005194_order$V3), FUN=sum)
MGYS00005194_family <- MGYS00005194_otu_table[which(MGYS00005194_otu_table$V4 != 'f__' & MGYS00005194_otu_table$V4 != ''), -c(1,2,3,5)]
MGYS00005194_family <- aggregate(MGYS00005194_family[, -1], by=list(MGYS00005194_family$V4), FUN=sum)
MGYS00005194_genus <- MGYS00005194_otu_table[which(MGYS00005194_otu_table$V5 != 'g__' & MGYS00005194_otu_table$V5 != ''), -c(1,2,3,4)]
MGYS00005194_genus <- aggregate(MGYS00005194_genus[, -1], by=list(MGYS00005194_genus$V5), FUN=sum)

# merge 5(phylum, class, order, family, genus)
data_phylum <- merge(merge(MGYS00005138_phylum,MGYS00005184_phylum, all = TRUE),MGYS00005194_phylum,all = TRUE)
data_class <- merge(merge(MGYS00005138_class,MGYS00005184_class, all = TRUE),MGYS00005194_class,all = TRUE)
data_order <- merge(merge(MGYS00005138_order,MGYS00005184_order, all = TRUE),MGYS00005194_order,all = TRUE)
data_family <- merge(merge(MGYS00005138_family,MGYS00005184_family, all = TRUE),MGYS00005194_family,all = TRUE)
data_genus <- merge(merge(MGYS00005138_genus,MGYS00005184_genus, all = TRUE),MGYS00005194_genus,all = TRUE)

data_phylum <- rename(data_phylum, 'phylum' = 'Group.1')
data_class <- rename(data_class, 'class' = 'Group.1')
data_order <- rename(data_order, 'order' = 'Group.1')
data_family <- rename(data_family, 'family' = 'Group.1')
data_genus <- rename(data_genus, 'genus' = 'Group.1')

data_phylum = setNames(data.frame(t(data_phylum[,-1])), data_phylum[,1])
data_class = setNames(data.frame(t(data_class[,-1])), data_class[,1])
data_order = setNames(data.frame(t(data_order[,-1])), data_order[,1])
data_family = setNames(data.frame(t(data_family[,-1])), data_family[,1])
data_genus = setNames(data.frame(t(data_genus[,-1])), data_genus[,1])

data_phylum[is.na(data_phylum)] <- 0
data_class[is.na(data_class)] <- 0
data_order[is.na(data_order)] <- 0
data_family[is.na(data_family)] <- 0
data_genus[is.na(data_genus)] <- 0

data_phylum <- cbind(rownames(data_phylum), data_phylum)
rownames(data_phylum) <- NULL
data_phylum <- rename(data_phylum, 'run_accession' = 'rownames(data_phylum)')
data_class <- cbind(rownames(data_class), data_class)
rownames(data_class) <- NULL
data_class <- rename(data_class, 'run_accession' = 'rownames(data_class)')
data_order <- cbind(rownames(data_order), data_order)
rownames(data_order) <- NULL
data_order <- rename(data_order, 'run_accession' = 'rownames(data_order)')
data_family <- cbind(rownames(data_family), data_family)
rownames(data_family) <- NULL
data_family <- rename(data_family, 'run_accession' = 'rownames(data_family)')
data_genus <- cbind(rownames(data_genus), data_genus)
rownames(data_genus) <- NULL
data_genus <- rename(data_genus, 'run_accession' = 'rownames(data_genus)')

# * 모든 feature의 데이터가 0인 샘플 제거
# * 분산이 0인 샘플 제거

##################################################################################################################################
### PCA & tSNE

# Host disease
# Study_arributes.bioproject

# Host sex
# Host Age

column_data <- data.frame(metadata$run_accession, metadata$Host_disease, metadata$study_attributes.bioproject, metadata$host_Age)
column_data <- rename(column_data, 'run_accession' = 'metadata.run_accession')
column_data <- rename(column_data, 'Host_disease' = 'metadata.Host_disease')
column_data <- rename(column_data, 'host_Age' = 'metadata.host_Age')
column_data <- rename(column_data, 'study_code' = 'metadata.study_attributes.bioproject')

##################################################################################################################################
# TMM Normalization

phylum <- column_data %>% inner_join(data_phylum)
class <- column_data %>% inner_join(data_class)
order <- column_data %>% inner_join(data_order)
family <- column_data %>% inner_join(data_family)
genus <- column_data %>% inner_join(data_genus)

tmp_column <- phylum[,c(1:4)]
tmp <- DGEList(t(as.matrix(phylum[,-c(1:4)])))
phylum <- calcNormFactors(tmp, method="TMM")
phylum <- cpm(phylum$counts, normalized.lib.sizes=TRUE, log=TRUE)
phylum <- t(phylum)
phylum <- cbind(tmp_column, phylum)

tmp_column <- class[,c(1:4)]
tmp <- DGEList(t(as.matrix(class[,-c(1:4)])))
class <- calcNormFactors(tmp, method="TMM")
class <- cpm(class$counts, normalized.lib.sizes=TRUE, log=TRUE)
class <- t(class)
class <- cbind(tmp_column, class)

tmp_column <- order[,c(1:4)]
tmp <- DGEList(t(as.matrix(order[,-c(1:4)])))
order <- calcNormFactors(tmp, method="TMM")
order <- cpm(order$counts, normalized.lib.sizes=TRUE, log=TRUE)
order <- t(order)
order <- cbind(tmp_column, order)

tmp_column <- family[,c(1:4)]
tmp <- DGEList(t(as.matrix(family[,-c(1:4)])))
family <- calcNormFactors(tmp, method="TMM")
family <- cpm(family$counts, normalized.lib.sizes=TRUE, log=TRUE)
family <- t(family)
family <- cbind(tmp_column, family)

tmp_column <- genus[,c(1:4)]
tmp <- DGEList(t(as.matrix(genus[,-c(1:4)])))
genus <- calcNormFactors(tmp, method="TMM")
genus <- cpm(genus$counts, normalized.lib.sizes=TRUE, log=TRUE)
genus <- t(genus)
genus <- cbind(tmp_column, genus)

set.seed(777)

pca_phylum <- prcomp(as.matrix(phylum[,-c(1,2,3,4)]))
pca_class <- prcomp(as.matrix(class[,-c(1,2,3,4)]))
pca_order <- prcomp(as.matrix(order[,-c(1,2,3,4)]))
pca_family <- prcomp(as.matrix(family[,-c(1,2,3,4)]))
pca_genus <- prcomp(as.matrix(genus[,-c(1,2,3,4)]))

par(mfrow=c(2, 3))

screeplot(pca_phylum, main = "PHYLUM PCA", col = "black", type = "lines", pch = 1, npcs = 10)
screeplot(pca_class, main = "CLASS PCA", col = "black", type = "lines", pch = 1, npcs = 10)
screeplot(pca_order, main = "ORDER PCA", col = "black", type = "lines", pch = 1, npcs = 10)
screeplot(pca_family, main = "FAMILY PCA", col = "black", type = "lines", pch = 1, npcs = 10)
screeplot(pca_genus, main = "GENUS PCA", col = "black", type = "lines", pch = 1, npcs = 10)

plot_pca_pylum <- autoplot(pca_phylum, data=phylum, colour="study_code", loadings=FALSE, loadings.colour = "black", scale = 0.5) +
  scale_colour_manual(values=c("forestgreen","red","blue")) +
  scale_fill_manual(values=c("forestgreen","red","blue")) +
  scale_shape_manual(values=c(25,22,23))

plot_pca_class <- autoplot(pca_class, data=class, colour="study_code", loadings=FALSE, loadings.colour = "black", scale = 0.5) +
  scale_colour_manual(values=c("forestgreen","red","blue")) +
  scale_fill_manual(values=c("forestgreen","red","blue")) +
  scale_shape_manual(values=c(25,22,23))

plot_pca_order <- autoplot(pca_order, data=order, colour="study_code", loadings=FALSE, loadings.colour = "black", scale = 0.5) +
  scale_colour_manual(values=c("forestgreen","red","blue")) +
  scale_fill_manual(values=c("forestgreen","red","blue")) +
  scale_shape_manual(values=c(25,22,23))

plot_pca_family <- autoplot(pca_family, data=family, colour="study_code", loadings=FALSE, loadings.colour = "black", scale = 0.5) +
  scale_colour_manual(values=c("forestgreen","red","blue")) +
  scale_fill_manual(values=c("forestgreen","red","blue")) +
  scale_shape_manual(values=c(25,22,23))

plot_pca_genus <- autoplot(pca_genus, data=genus, colour="study_code", loadings=FALSE, loadings.colour = "black", scale = 0.5) +
  scale_colour_manual(values=c("forestgreen","red","blue")) +
  scale_fill_manual(values=c("forestgreen","red","blue")) +
  scale_shape_manual(values=c(25,22,23))

grid.arrange(plot_pca_pylum, plot_pca_class, plot_pca_order, plot_pca_family, plot_pca_genus, nrow=2, ncol=3)

plot_pca_pylum <- autoplot(pca_phylum, data=phylum, colour="Host_disease", loadings=FALSE, loadings.colour = "black", scale = 0.5) +
  scale_colour_manual(values=c("forestgreen","red","blue")) +
  scale_fill_manual(values=c("forestgreen","red","blue")) +
  scale_shape_manual(values=c(25,22,23))

plot_pca_class <- autoplot(pca_class, data=class, colour="Host_disease", loadings=FALSE, loadings.colour = "black", scale = 0.5) +
  scale_colour_manual(values=c("forestgreen","red","blue")) +
  scale_fill_manual(values=c("forestgreen","red","blue")) +
  scale_shape_manual(values=c(25,22,23))

plot_pca_order <- autoplot(pca_order, data=order, colour="Host_disease", loadings=FALSE, loadings.colour = "black", scale = 0.5) +
  scale_colour_manual(values=c("forestgreen","red","blue")) +
  scale_fill_manual(values=c("forestgreen","red","blue")) +
  scale_shape_manual(values=c(25,22,23))

plot_pca_family <- autoplot(pca_family, data=family, colour="Host_disease", loadings=FALSE, loadings.colour = "black", scale = 0.5) +
  scale_colour_manual(values=c("forestgreen","red","blue")) +
  scale_fill_manual(values=c("forestgreen","red","blue")) +
  scale_shape_manual(values=c(25,22,23))

plot_pca_genus <- autoplot(pca_genus, data=genus, colour="Host_disease", loadings=FALSE, loadings.colour = "black", scale = 0.5) +
  scale_colour_manual(values=c("forestgreen","red","blue")) +
  scale_fill_manual(values=c("forestgreen","red","blue")) +
  scale_shape_manual(values=c(25,22,23))

grid.arrange(plot_pca_pylum, plot_pca_class, plot_pca_order, plot_pca_family, plot_pca_genus, nrow=2, ncol=3)

##################################################################################################################################

tsne_phylum <- Rtsne(as.matrix(phylum[,-c(1,2,3,4)]), PCA = FALSE, check_duplicates = FALSE, perplexity=30, theta=0.5, dims=2, num_threads = 30)
tsne_class <- Rtsne(as.matrix(class[,-c(1,2,3,4)]), PCA = FALSE, check_duplicates = FALSE, perplexity=30, theta=0.5, dims=2, num_threads = 30)
tsne_order <- Rtsne(as.matrix(order[,-c(1,2,3,4)]), PCA = FALSE, check_duplicates = FALSE, perplexity=30, theta=0.5, dims=2, num_threads = 30)
tsne_family <- Rtsne(as.matrix(family[,-c(1,2,3,4)]), PCA = FALSE, check_duplicates = FALSE, perplexity=30, theta=0.5, dims=2, num_threads = 30)
tsne_genus <- Rtsne(as.matrix(genus[,-c(1,2,3,4)]), PCA = FALSE, check_duplicates = FALSE, perplexity=30, theta=0.5, dims=2, num_threads = 30)

View(tsne_phylum$Y)

test_hd <- as.data.frame(phylum$Host_disease)
which(test_hd$`phylum$Host_disease`== 'Healthy')
unique(phylum$Host_disease)
test_hd[which(test_hd$`phylum$Host_disease`== 'Healthy'),] <- 1
test_hd[which(test_hd$`phylum$Host_disease`== 'Chronic rhinosinusitis'),] <- 2
test_hd[which(test_hd$`phylum$Host_disease`== 'Acute Respiratory Infection'),] <- 3

test <- cbind(tsne_phylum$Y, test_hd)

write.csv(test, '/home/bicjh/CI_final/tsne_genus.csv')

library(plotly)
plot_ly(data = as.data.frame(tsne_phylum$Y),x =  ~V1, y = ~V2, z = ~V3, color = ~genus$study_code, colors = c('#636EFA','#EF553B','#00CC96')) %>% 
  add_markers(size = 1)





par(mfrow=c(2, 3))

plot_t_p <- ggplot(as.data.frame(tsne_phylum$Y)) +
  geom_point(aes(x=tsne_phylum$Y[,1], y=tsne_phylum$Y[,2], color=phylum$Host_disease)) +
  ggtitle(label = "Phylum tSNE") + theme(legend.position = "none")

plot_t_c <- ggplot(as.data.frame(tsne_class$Y)) +
  geom_point(aes(x=tsne_class$Y[,1], y=tsne_class$Y[,2], color=class$Host_disease)) +
  ggtitle(label = "Class tSNE") + theme(legend.position = "none")

plot_t_o <- ggplot(as.data.frame(tsne_order$Y)) +
  geom_point(aes(x=tsne_order$Y[,1], y=tsne_order$Y[,2], color=order$Host_disease)) +
  ggtitle(label = "Order tSNE") + theme(legend.position = "none")

plot_t_f <- ggplot(as.data.frame(tsne_family$Y)) +
  geom_point(aes(x=tsne_family$Y[,1], y=tsne_family$Y[,2], color=family$Host_disease)) +
  ggtitle(label = "Family tSNE") + theme(legend.position = "none")

plot_t_g <- ggplot(as.data.frame(tsne_genus$Y)) +
  geom_point(aes(x=tsne_genus$Y[,1], y=tsne_genus$Y[,2], color=genus$Host_disease)) +
  ggtitle(label = "Genus tSNE") + theme(legend.position = "none")

ggarrange(plot_t_p, plot_t_c, plot_t_o, plot_t_f, plot_t_g, common.legend = TRUE)

plot(tsne_phylum$Y, col=as.factor(phylum$study_code), pch = 3, cex =0.5)
legend("bottomright", legend = levels(factor(phylum$study_code)), pch = 19, col = factor(levels(factor(phylum$study_code))), cex=0.3)
plot(tsne_class$Y, col=as.factor(class$study_code), pch = 3, cex =0.5)
legend("bottomright", legend = levels(factor(class$study_code)), pch = 19, col = factor(levels(factor(class$study_code))), cex=0.3)
plot(tsne_order$Y, col=as.factor(order$study_code), pch = 3, cex =0.5)
legend("bottomright", legend = levels(factor(order$study_code)), pch = 19, col = factor(levels(factor(order$study_code))), cex=0.3)
plot(tsne_family$Y, col=as.factor(family$study_code), pch = 3, cex =0.5)
legend("bottomright", legend = levels(factor(family$study_code)), pch = 19, col = factor(levels(factor(family$study_code))), cex=0.3)
plot(tsne_genus$Y, col=as.factor(genus$study_code), pch = 3, cex =0.5)
legend("bottomright", legend = levels(factor(ggg$study_code)), pch = 19, col = factor(levels(factor(ggg$study_code))), cex=0.3)


plot(tsne_phylum$Y, col=as.factor(phylum$host_Age), pch = 3, cex =0.5)
legend("bottomright", legend = levels(factor(phylum$host_Age)), pch = 19, col = factor(levels(factor(phylum$host_Age))), cex=0.3)
plot(tsne_class$Y, col=as.factor(class$host_Age), pch = 3, cex =0.5)
legend("bottomright", legend = levels(factor(class$host_Age)), pch = 19, col = factor(levels(factor(class$host_Age))), cex=0.3)
plot(tsne_order$Y, col=as.factor(order$host_Age), pch = 3, cex =0.5)
legend("bottomright", legend = levels(factor(order$host_Age)), pch = 19, col = factor(levels(factor(order$host_Age))), cex=0.3)
plot(tsne_family$Y, col=as.factor(family$host_Age), pch = 3, cex =0.5)
legend("bottomright", legend = levels(factor(family$host_Age)), pch = 19, col = factor(levels(factor(family$host_Age))), cex=0.3)
plot(tsne_genus$Y, col=as.factor(genus$host_Age), pch = 3, cex =0.5)
legend("bottomright", legend = levels(factor(genus$host_Age)), pch = 19, col = factor(levels(factor(genus$host_Age))), cex=0.3)

##################################################################################################################################

sil_phylum <- c(NA)
sil_class <- c(NA)
sil_order <- c(NA)
sil_family <- c(NA)
sil_genus <- c(NA)
#foreach(idx = 2:10,.packages = 'cluster') %dopar% {
for(idx in 2:10) {
  pam_fit_phylum <- pam(dist(data_phylum[,-c(1)]), diss=TRUE, k=idx)
  pam_fit_class <- pam(dist(data_class[,-c(1)]), diss=TRUE, k=idx)
  pam_fit_order <- pam(dist(data_order[,-c(1)]), diss=TRUE, k=idx)
  pam_fit_family <- pam(dist(data_family[,-c(1)]), diss=TRUE, k=idx)
  pam_fit_genus <- pam(dist(data_genus[,-c(1)]), diss=TRUE, k=idx)
  sil_phylum[idx] <- pam_fit_phylum$silinfo$avg.width
  sil_class[idx] <- pam_fit_class$silinfo$avg.width
  sil_order[idx] <- pam_fit_order$silinfo$avg.width
  sil_family[idx] <- pam_fit_family$silinfo$avg.width
  sil_genus[idx] <- pam_fit_genus$silinfo$avg.width
  print(idx)
}
par(mfrow=c(2, 3))

plot(1:10, sil_phylum, xlab = "Number of clusters", ylab = "Silhouette Width")
lines(1:10, sil_phylum)
plot(1:10, sil_class, xlab = "Number of clusters", ylab = "Silhouette Width")
lines(1:10, sil_class)
plot(1:10, sil_order, xlab = "Number of clusters", ylab = "Silhouette Width")
lines(1:10, sil_order)
plot(1:10, sil_family, xlab = "Number of clusters", ylab = "Silhouette Width")
lines(1:10, sil_family)
plot(1:10, sil_genus, xlab = "Number of clusters", ylab = "Silhouette Width")
lines(1:10, sil_genus)
View(sil_phylum)

plot_s_p <- ggplot(as.data.frame(sil_phylum)) +
  geom_point(aes(x=1:10, y=sil_phylum)) +
  geom_line(aes(x=1:10, y=sil_phylum)) +
  ggtitle(label = "Phylum Silhouette score")+
  scale_y_continuous(breaks = seq(-1, 1, 0.2),labels = scales::percent, limits = c(0, 1))

plot_s_c <- ggplot(as.data.frame(sil_class)) +
  geom_point(aes(x=1:10, y=sil_class)) +
  geom_line(aes(x=1:10, y=sil_class)) +
  ggtitle(label = "Class Silhouette score")+
  scale_y_continuous(breaks = seq(-1, 1, 0.2),labels = scales::percent, limits = c(0, 1))

plot_s_o <- ggplot(as.data.frame(sil_order)) +
  geom_point(aes(x=1:10, y=sil_order)) +
  geom_line(aes(x=1:10, y=sil_order)) +
  ggtitle(label = "Order Silhouette score")+
  scale_y_continuous(breaks = seq(-1, 1, 0.2),labels = scales::percent, limits = c(0, 1))

plot_s_f <- ggplot(as.data.frame(sil_family)) +
  geom_point(aes(x=1:10, y=sil_family)) +
  geom_line(aes(x=1:10, y=sil_family)) +
  ggtitle(label = "Family Silhouette score")+
  scale_y_continuous(breaks = seq(-1, 1, 0.2),labels = scales::percent, limits = c(0, 1))

plot_s_g <- ggplot(as.data.frame(sil_genus)) +
  geom_point(aes(x=1:10, y=sil_genus)) +
  geom_line(aes(x=1:10, y=sil_genus)) +
  ggtitle(label = "Genus Silhouette score")+
  scale_y_continuous(breaks = seq(-1, 1, 0.2),labels = scales::percent, limits = c(0, 1))

grid.arrange(plot_s_p, plot_s_c, plot_s_o, plot_s_f, plot_s_g, nrow=2, ncol=3)





#### 수정필요
hc <- hclust(dist(data_phylum), method='average')
cutree(hc,k=3)
silhouette(hc)
test <- pam(dist(phylum[,-c(1,2,3,4)]), k=3)
View(phylum[,-c(1,2,3,4)])
colnames(phylum)
table(test, )
plot(test)
plot(silhouette(cutree(hc,k=20),dist=dist(data_phylum[,-c(1,2,3,4)]),col=1:20))
fviz_nbclust(metadata, kmeans, method = "silhouette")

# plot(silhouette(cutree(hc,k=7),dist=dist(s_data),col=1:7))

plot(hc, k=3)
rect.hclust(hc,k=3,border="red")
plot(silhouette(cutree(hc,k=3),dist=metadata,col=1:3))

##################################################################################################################################
registerDoMC(cores = 40)
getDoParWorkers()

# K-Means Clustering
train_phylum <- phylum[,-c(1,3,4)]
train_class <- class[,-c(1,3,4)]
train_order <- order[,-c(1,3,4)]
train_family <- family[,-c(1,3,4)]
train_genus <- genus[,-c(1,3,4)]

clustering_phylum <- kmeans(train_phylum[,-1], center = 3, iter.max = 10000)
clustering_class <- kmeans(train_class[,-1], center = 3, iter.max = 10000)
clustering_order <- kmeans(train_order[,-1], center = 3, iter.max = 10000)
clustering_family <- kmeans(train_family[,-1], center = 3, iter.max = 10000)
clustering_genus <- kmeans(train_genus[,-1], center = 3, iter.max = 10000)

sil_phylum <- silhouette(clustering_phylum$cluster, dist(train_phylum[,-1]))
sil_class <- silhouette(clustering_class$cluster, dist(train_class[,-1]))
sil_order <- silhouette(clustering_order$cluster, dist(train_order[,-1]))
sil_family <- silhouette(clustering_family$cluster, dist(train_family[,-1]))
sil_genus <- silhouette(clustering_genus$cluster, dist(train_genus[,-1]))

plot_fc_p <- fviz_cluster(clustering_phylum, train_phylum[, -1], ellipse.type = "norm", geom = "point", main = "Phylum cluster map")
plot_fc_c <- fviz_cluster(clustering_class, train_class[, -1], ellipse.type = "norm", geom = "point", main = "Class cluster map")
plot_fc_o <- fviz_cluster(clustering_order, train_order[, -1], ellipse.type = "norm", geom = "point", main = "Order cluster map")
plot_fc_f <- fviz_cluster(clustering_family, train_family[, -1], ellipse.type = "norm", geom = "point", main = "Family cluster map")
plot_fc_g <- fviz_cluster(clustering_genus, train_genus[, -1], ellipse.type = "norm", geom = "point", main = "Genus cluster map")

grid.arrange(plot_fc_p, plot_fc_c, plot_fc_o, plot_fc_f, plot_fc_g, nrow=2, ncol=3)

# Kmeans clustering
# Hierarchical clustering
# PAM clustering

plot_fs_p <- fviz_silhouette(sil_phylum)
plot_fs_c <- fviz_silhouette(sil_class)
plot_fs_o <- fviz_silhouette(sil_order)
plot_fs_f <- fviz_silhouette(sil_family)
plot_fs_g <- fviz_silhouette(sil_genus)

?fviz_silhouette
grid.arrange(plot_fs_p, plot_fs_c, plot_fs_o, plot_fs_f, plot_fs_g, nrow=2, ncol=3)

##################################################################################################################################

registerDoMC(cores = 40)
getDoParWorkers()

train_phylum <- phylum[,-c(1,3,4)]
train_class <- class[,-c(1,3,4)]
train_order <- order[,-c(1,3,4)]
train_family <- family[,-c(1,3,4)]
train_genus <- genus[,-c(1,3,4)]

customMultiClassSummary <- function (data, lev = NULL, model = NULL) {
  if (!all(levels(data[, "pred"]) == levels(data[, "obs"]))) 
    stop("levels of observed and predicted data do not match")
  has_class_probs <- all(lev %in% colnames(data))
  if (has_class_probs) {
    lloss <- mnLogLoss(data = data, lev = lev, model = model)
    requireNamespaceQuietStop("pROC")
    requireNamespaceQuietStop("MLmetrics")
    prob_stats <- lapply(levels(data[, "pred"]), function(x) {
      obs <- ifelse(data[, "obs"] == x, 1, 0)
      prob <- data[, x]
      roc_auc <- try(pROC::roc(obs, data[, x], direction = "<", 
                               quiet = TRUE), silent = TRUE)
      roc_auc <- if (inherits(roc_auc, "try-error")) 
        NA
      else roc_auc$auc
      pr_auc <- try(MLmetrics::PRAUC(y_pred = data[, x], 
                                     y_true = obs), silent = TRUE)
      if (inherits(pr_auc, "try-error")) 
        pr_auc <- NA
      res <- c(ROC = roc_auc, AUC = pr_auc)
      return(res)
    })
    prob_stats <- do.call("rbind", prob_stats)
    prob_stats <- colMeans(prob_stats, na.rm = TRUE)
  }
  CM <- confusionMatrix(data[, "pred"], data[, "obs"], mode = "everything")
  if (length(levels(data[, "pred"])) == 2) {
    class_stats <- CM$byClass
  }
  else {
    class_stats <- colMeans(CM$byClass)
    names(class_stats) <- paste("Mean", names(class_stats))
  }
  overall_stats <- if (has_class_probs) 
    c(CM$overall, logLoss = as.numeric(lloss), AUC = unname(prob_stats["ROC"]), 
      prAUC = unname(prob_stats["AUC"]))
  else CM$overall
  stats <- c(overall_stats, class_stats)
  print('시행!')
  stats <- stats[!names(stats) %in% c("AccuracyNull", "AccuracyLower", 
                                      "AccuracyUpper", "AccuracyPValue", "McnemarPValue", 
                                      "Mean Prevalence", "Mean Detection Prevalence")]
  names(stats) <- gsub("[[:blank:]]+", "_", names(stats))
  stat_list <- c("Accuracy", "Kappa", "Mean_F1", "Mean_Sensitivity", 
                 "Mean_Specificity", "Mean_Pos_Pred_Value", "Mean_Neg_Pred_Value", 
                 "Mean_Precision", "Mean_Recall", "Mean_Detection_Rate", 
                 "Mean_Balanced_Accuracy")
  if (has_class_probs) 
    stat_list <- c("logLoss", "AUC", "prAUC", stat_list)
  if (length(levels(data[, "pred"])) == 2) 
    stat_list <- gsub("^Mean_", "", stat_list)
  stats <- stats[c(stat_list)]
  return(stats)
}

customSummary  <- function(data, lev = NULL, model = NULL){
  a1 <- defaultSummary(data, lev, model)
  b1 <- customMultiClassSummary(data, lev, model)
  # b1 <- multiClassSummary(data, lev, model)
  out <- c(a1,b1)
  out
  }

fitControl <- trainControl(method = "repeatedcv", number = 10, repeats = 3, savePredictions = TRUE, summaryFunction=customSummary)

model_phylum_randomforest <- train(Host_disease ~ ., data = train_phylum, method = 'rf', trControl = fitControl,verbose = F)
model_class_randomforest <- train(Host_disease ~ ., data = train_class, method = 'rf', trControl = fitControl,verbose = F)
model_order_randomforest <- train(Host_disease ~ ., data = train_order, method = 'rf', trControl = fitControl,verbose = F)
model_family_randomforest <- train(Host_disease ~ ., data = train_family, method = 'rf', trControl = fitControl,verbose = F)
model_genus_randomforest <- train(Host_disease ~ ., data = train_genus, method = 'rf', trControl = fitControl,verbose = F)

model_phylum_c50 <- train(Host_disease ~ ., data = train_phylum, method = 'C5.0', trControl = fitControl,verbose = F)
model_class_c50 <- train(Host_disease ~ ., data = train_class, method = 'C5.0', trControl = fitControl,verbose = F)
model_order_c50 <- train(Host_disease ~ ., data = train_order, method = 'C5.0', trControl = fitControl,verbose = F)
model_family_c50 <- train(Host_disease ~ ., data = train_family, method = 'C5.0', trControl = fitControl,verbose = F)
model_genus_c50 <- train(Host_disease ~ ., data = train_genus, method = 'C5.0', trControl = fitControl,verbose = F)

model_phylum_ranger <- train(Host_disease ~ ., data = train_phylum, method = 'ranger', trControl = fitControl,verbose = F)
model_class_ranger <- train(Host_disease ~ ., data = train_class, method = 'ranger', trControl = fitControl,verbose = F)
model_order_ranger <- train(Host_disease ~ ., data = train_order, method = 'ranger', trControl = fitControl,verbose = F)
model_family_ranger <- train(Host_disease ~ ., data = train_family, method = 'ranger', trControl = fitControl,verbose = F)
model_genus_ranger <- train(Host_disease ~ ., data = train_genus, method = 'ranger', trControl = fitControl,verbose = F)

model_phylum_Linda <- train(Host_disease ~ ., data = train_phylum, method = 'Linda', trControl = fitControl,verbose = F)
model_class_Linda <- train(Host_disease ~ ., data = train_class, method = 'Linda', trControl = fitControl,verbose = F)
model_order_Linda <- train(Host_disease ~ ., data = train_order, method = 'Linda', trControl = fitControl,verbose = F)
model_family_Linda <- train(Host_disease ~ ., data = train_family, method = 'Linda', trControl = fitControl,verbose = F)
model_genus_Linda <- train(Host_disease ~ ., data = train_genus, method = 'Linda', trControl = fitControl,verbose = F)

model_phylum_LogitBoost <- train(Host_disease ~ ., data = train_phylum, method = 'LogitBoost', trControl = fitControl,verbose = F)
model_class_LogitBoost <- train(Host_disease ~ ., data = train_class, method = 'LogitBoost', trControl = fitControl,verbose = F)
model_order_LogitBoost <- train(Host_disease ~ ., data = train_order, method = 'LogitBoost', trControl = fitControl,verbose = F)
model_family_LogitBoost <- train(Host_disease ~ ., data = train_family, method = 'LogitBoost', trControl = fitControl,verbose = F)
model_genus_LogitBoost <- train(Host_disease ~ ., data = train_genus, method = 'LogitBoost', trControl = fitControl,verbose = F)

model_phylum_svmRadial <- train(Host_disease ~ ., data = train_phylum, method = 'svmRadial', trControl = fitControl,verbose = F)
model_class_svmRadial <- train(Host_disease ~ ., data = train_class, method = 'svmRadial', trControl = fitControl,verbose = F)
model_order_svmRadial <- train(Host_disease ~ ., data = train_order, method = 'svmRadial', trControl = fitControl,verbose = F)
model_family_svmRadial <- train(Host_disease ~ ., data = train_family, method = 'svmRadial', trControl = fitControl,verbose = F)
model_genus_svmRadial <- train(Host_disease ~ ., data = train_genus, method = 'svmRadial', trControl = fitControl,verbose = F)

model_phylum_svmRadialWeights <- train(Host_disease ~ ., data = train_phylum, method = 'svmRadialWeights', trControl = fitControl,verbose = F)
model_class_svmRadialWeights <- train(Host_disease ~ ., data = train_class, method = 'svmRadialWeights', trControl = fitControl,verbose = F)
model_order_svmRadialWeights <- train(Host_disease ~ ., data = train_order, method = 'svmRadialWeights', trControl = fitControl,verbose = F)
model_family_svmRadialWeights <- train(Host_disease ~ ., data = train_family, method = 'svmRadialWeights', trControl = fitControl,verbose = F)
model_genus_svmRadialWeights <- train(Host_disease ~ ., data = train_genus, method = 'svmRadialWeights', trControl = fitControl,verbose = F)

model_phylum_LogitBoost <- train(Host_disease ~ ., data = train_phylum, method = 'LogitBoost', trControl = fitControl,verbose = F)
model_class_LogitBoost <- train(Host_disease ~ ., data = train_class, method = 'LogitBoost', trControl = fitControl,verbose = F)
model_order_LogitBoost <- train(Host_disease ~ ., data = train_order, method = 'LogitBoost', trControl = fitControl,verbose = F)
model_family_LogitBoost <- train(Host_disease ~ ., data = train_family, method = 'LogitBoost', trControl = fitControl,verbose = F)
model_genus_LogitBoost <- train(Host_disease ~ ., data = train_genus, method = 'LogitBoost', trControl = fitControl,verbose = F)

model_phylum_knn <- train(Host_disease ~ ., data = train_phylum, method = 'kknn', trControl = fitControl,verbose = F)
model_class_knn <- train(Host_disease ~ ., data = train_class, method = 'kknn', trControl = fitControl,verbose = F)
model_order_knn <- train(Host_disease ~ ., data = train_order, method = 'kknn', trControl = fitControl,verbose = F)
model_family_knn <- train(Host_disease ~ ., data = train_family, method = 'kknn', trControl = fitControl,verbose = F)
model_genus_knn <- train(Host_disease ~ ., data = train_genus, method = 'kknn', trControl = fitControl,verbose = F)
model_phylum_knn
##################################################################################################################################

accuracy_phylum <- c()
tmp <- data.frame(model_phylum_randomforest$resample[,'Accuracy'], rep('RandomForest', 30))
colnames(tmp) <- c('value', 'model')
accuracy_phylum <- rbind(accuracy_phylum, tmp)
tmp <- data.frame(model_phylum_c50$resample[,'Accuracy'], rep('C5.0', 30))
colnames(tmp) <- c('value', 'model')
accuracy_phylum <- rbind(accuracy_phylum, tmp)
tmp <- data.frame(model_phylum_ranger$resample[,'Accuracy'], rep('Ranger', 30))
colnames(tmp) <- c('value', 'model')
accuracy_phylum <- rbind(accuracy_phylum, tmp)
tmp <- data.frame(model_phylum_Linda$resample[,'Accuracy'], rep('Linda', 30))
colnames(tmp) <- c('value', 'model')
accuracy_phylum <- rbind(accuracy_phylum, tmp)
tmp <- data.frame(model_phylum_LogitBoost$resample[,'Accuracy'], rep('LogitBoost', 30))
colnames(tmp) <- c('value', 'model')
accuracy_phylum <- rbind(accuracy_phylum, tmp)
tmp <- data.frame(model_phylum_svmRadial$resample[,'Accuracy'], rep('svmRadial', 30))
colnames(tmp) <- c('value', 'model')
accuracy_phylum <- rbind(accuracy_phylum, tmp)
tmp <- data.frame(model_phylum_svmRadialWeights$resample[,'Accuracy'], rep('svmRadialWeights', 30))
colnames(tmp) <- c('value', 'model')
accuracy_phylum <- rbind(accuracy_phylum, tmp)
tmp <- data.frame(model_phylum_LogitBoost$resample[,'Accuracy'], rep('LogitBoost', 30))
colnames(tmp) <- c('value', 'model')
accuracy_phylum <- rbind(accuracy_phylum, tmp)
tmp <- data.frame(model_phylum_knn$resample[,'Accuracy'], rep('knn', 30))
colnames(tmp) <- c('value', 'model')
accuracy_phylum <- rbind(accuracy_phylum, tmp)

sensitivity_phylum <- c()
tmp <- data.frame(model_phylum_randomforest$resample[,'Mean_Sensitivity'], rep('RandomForest', 30))
colnames(tmp) <- c('value', 'model')
sensitivity_phylum <- rbind(sensitivity_phylum, tmp)
tmp <- data.frame(model_phylum_c50$resample[,'Mean_Sensitivity'], rep('C5.0', 30))
colnames(tmp) <- c('value', 'model')
sensitivity_phylum <- rbind(sensitivity_phylum, tmp)
tmp <- data.frame(model_phylum_ranger$resample[,'Mean_Sensitivity'], rep('Ranger', 30))
colnames(tmp) <- c('value', 'model')
sensitivity_phylum <- rbind(sensitivity_phylum, tmp)
tmp <- data.frame(model_phylum_Linda$resample[,'Mean_Sensitivity'], rep('Linda', 30))
colnames(tmp) <- c('value', 'model')
sensitivity_phylum <- rbind(sensitivity_phylum, tmp)
tmp <- data.frame(model_phylum_LogitBoost$resample[,'Mean_Sensitivity'], rep('LogitBoost', 30))
colnames(tmp) <- c('value', 'model')
sensitivity_phylum <- rbind(sensitivity_phylum, tmp)
tmp <- data.frame(model_phylum_svmRadial$resample[,'Mean_Sensitivity'], rep('svmRadial', 30))
colnames(tmp) <- c('value', 'model')
sensitivity_phylum <- rbind(sensitivity_phylum, tmp)
tmp <- data.frame(model_phylum_svmRadialWeights$resample[,'Mean_Sensitivity'], rep('svmRadialWeights', 30))
colnames(tmp) <- c('value', 'model')
sensitivity_phylum <- rbind(sensitivity_phylum, tmp)
tmp <- data.frame(model_phylum_LogitBoost$resample[,'Mean_Sensitivity'], rep('LogitBoost', 30))
colnames(tmp) <- c('value', 'model')
sensitivity_phylum <- rbind(sensitivity_phylum, tmp)
tmp <- data.frame(model_phylum_knn$resample[,'Mean_Sensitivity'], rep('knn', 30))
colnames(tmp) <- c('value', 'model')
sensitivity_phylum <- rbind(sensitivity_phylum, tmp)

specificity_phylum <- c()
tmp <- data.frame(model_phylum_randomforest$resample[,'Mean_Specificity'], rep('RandomForest', 30))
colnames(tmp) <- c('value', 'model')
specificity_phylum <- rbind(specificity_phylum, tmp)
tmp <- data.frame(model_phylum_c50$resample[,'Mean_Specificity'], rep('C5.0', 30))
colnames(tmp) <- c('value', 'model')
specificity_phylum <- rbind(specificity_phylum, tmp)
tmp <- data.frame(model_phylum_ranger$resample[,'Mean_Specificity'], rep('Ranger', 30))
colnames(tmp) <- c('value', 'model')
specificity_phylum <- rbind(specificity_phylum, tmp)
tmp <- data.frame(model_phylum_Linda$resample[,'Mean_Specificity'], rep('Linda', 30))
colnames(tmp) <- c('value', 'model')
specificity_phylum <- rbind(specificity_phylum, tmp)
tmp <- data.frame(model_phylum_LogitBoost$resample[,'Mean_Specificity'], rep('LogitBoost', 30))
colnames(tmp) <- c('value', 'model')
specificity_phylum <- rbind(specificity_phylum, tmp)
tmp <- data.frame(model_phylum_svmRadial$resample[,'Mean_Specificity'], rep('svmRadial', 30))
colnames(tmp) <- c('value', 'model')
specificity_phylum <- rbind(specificity_phylum, tmp)
tmp <- data.frame(model_phylum_svmRadialWeights$resample[,'Mean_Specificity'], rep('svmRadialWeights', 30))
colnames(tmp) <- c('value', 'model')
specificity_phylum <- rbind(specificity_phylum, tmp)
tmp <- data.frame(model_phylum_LogitBoost$resample[,'Mean_Specificity'], rep('LogitBoost', 30))
colnames(tmp) <- c('value', 'model')
specificity_phylum <- rbind(specificity_phylum, tmp)
tmp <- data.frame(model_phylum_knn$resample[,'Mean_Specificity'], rep('knn', 30))
colnames(tmp) <- c('value', 'model')
specificity_phylum <- rbind(specificity_phylum, tmp)

plot_ac_p <- ggplot(accuracy_phylum) +
  geom_boxplot(aes(x=model, y=value, fill=model)) +
  geom_jitter(aes(x=model, y=value), alpha=0.2) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1L)) +
  ggtitle("Phylum_Accuracy")

plot_ss_p <- ggplot(sensitivity_phylum) +
  geom_boxplot(aes(x=model, y=value, fill=model)) +
  geom_jitter(aes(x=model, y=value), alpha=0.2) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1L)) +
  ggtitle("Phylum_Sensitivity") +
  theme(legend.position = 'none')

plot_sp_c <- ggplot(specificity_phylum) +
  geom_boxplot(aes(x=model, y=value, fill=model)) +
  geom_jitter(aes(x=model, y=value), alpha=0.2) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1L)) +
  ggtitle("Phylum_Specificity") +
  theme(legend.position = 'none')

plot_ac_p | (plot_ss_p / plot_sp_c)

##################################################################################################################################

accuracy_class <- c()
tmp <- data.frame(model_class_randomforest$resample[,'Accuracy'], rep('RandomForest', 30))
colnames(tmp) <- c('value', 'model')
accuracy_class <- rbind(accuracy_class, tmp)
tmp <- data.frame(model_class_c50$resample[,'Accuracy'], rep('C5.0', 30))
colnames(tmp) <- c('value', 'model')
accuracy_class <- rbind(accuracy_class, tmp)
tmp <- data.frame(model_class_ranger$resample[,'Accuracy'], rep('Ranger', 30))
colnames(tmp) <- c('value', 'model')
accuracy_class <- rbind(accuracy_class, tmp)
tmp <- data.frame(model_class_Linda$resample[,'Accuracy'], rep('Linda', 30))
colnames(tmp) <- c('value', 'model')
accuracy_class <- rbind(accuracy_class, tmp)
tmp <- data.frame(model_class_LogitBoost$resample[,'Accuracy'], rep('LogitBoost', 30))
colnames(tmp) <- c('value', 'model')
accuracy_class <- rbind(accuracy_class, tmp)
tmp <- data.frame(model_class_svmRadial$resample[,'Accuracy'], rep('svmRadial', 30))
colnames(tmp) <- c('value', 'model')
accuracy_class <- rbind(accuracy_class, tmp)
tmp <- data.frame(model_class_svmRadialWeights$resample[,'Accuracy'], rep('svmRadialWeights', 30))
colnames(tmp) <- c('value', 'model')
accuracy_class <- rbind(accuracy_class, tmp)

sensitivity_class <- c()
tmp <- data.frame(model_class_randomforest$resample[,'Mean_Sensitivity'], rep('RandomForest', 30))
colnames(tmp) <- c('value', 'model')
sensitivity_class <- rbind(sensitivity_class, tmp)
tmp <- data.frame(model_class_c50$resample[,'Mean_Sensitivity'], rep('C5.0', 30))
colnames(tmp) <- c('value', 'model')
sensitivity_class <- rbind(sensitivity_class, tmp)
tmp <- data.frame(model_class_ranger$resample[,'Mean_Sensitivity'], rep('Ranger', 30))
colnames(tmp) <- c('value', 'model')
sensitivity_class <- rbind(sensitivity_class, tmp)
tmp <- data.frame(model_class_Linda$resample[,'Mean_Sensitivity'], rep('Linda', 30))
colnames(tmp) <- c('value', 'model')
sensitivity_class <- rbind(sensitivity_class, tmp)
tmp <- data.frame(model_class_LogitBoost$resample[,'Mean_Sensitivity'], rep('LogitBoost', 30))
colnames(tmp) <- c('value', 'model')
sensitivity_class <- rbind(sensitivity_class, tmp)
tmp <- data.frame(model_class_svmRadial$resample[,'Mean_Sensitivity'], rep('svmRadial', 30))
colnames(tmp) <- c('value', 'model')
sensitivity_class <- rbind(sensitivity_class, tmp)
tmp <- data.frame(model_class_svmRadialWeights$resample[,'Mean_Sensitivity'], rep('svmRadialWeights', 30))
colnames(tmp) <- c('value', 'model')
sensitivity_class <- rbind(sensitivity_class, tmp)

specificity_class <- c()
tmp <- data.frame(model_class_randomforest$resample[,'Mean_Specificity'], rep('RandomForest', 30))
colnames(tmp) <- c('value', 'model')
specificity_class <- rbind(specificity_class, tmp)
tmp <- data.frame(model_class_c50$resample[,'Mean_Specificity'], rep('C5.0', 30))
colnames(tmp) <- c('value', 'model')
specificity_class <- rbind(specificity_class, tmp)
tmp <- data.frame(model_class_ranger$resample[,'Mean_Specificity'], rep('Ranger', 30))
colnames(tmp) <- c('value', 'model')
specificity_class <- rbind(specificity_class, tmp)
tmp <- data.frame(model_class_Linda$resample[,'Mean_Specificity'], rep('Linda', 30))
colnames(tmp) <- c('value', 'model')
specificity_class <- rbind(specificity_class, tmp)
tmp <- data.frame(model_class_LogitBoost$resample[,'Mean_Specificity'], rep('LogitBoost', 30))
colnames(tmp) <- c('value', 'model')
specificity_class <- rbind(specificity_class, tmp)
tmp <- data.frame(model_class_svmRadial$resample[,'Mean_Specificity'], rep('svmRadial', 30))
colnames(tmp) <- c('value', 'model')
specificity_class <- rbind(specificity_class, tmp)
tmp <- data.frame(model_class_svmRadialWeights$resample[,'Mean_Specificity'], rep('svmRadialWeights', 30))
colnames(tmp) <- c('value', 'model')
specificity_class <- rbind(specificity_class, tmp)

plot_ac_c <- ggplot(accuracy_class) +
  geom_boxplot(aes(x=model, y=value, fill=model)) +
  geom_jitter(aes(x=model, y=value), alpha=0.2) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1L)) +
  ggtitle("Class_Accuracy")

plot_ss_c <- ggplot(sensitivity_class) +
  geom_boxplot(aes(x=model, y=value, fill=model)) +
  geom_jitter(aes(x=model, y=value), alpha=0.2) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1L)) +
  ggtitle("Class_Sensitivity") + 
  theme(legend.position = 'none')

plot_sp_c <- ggplot(specificity_class) +
  geom_boxplot(aes(x=model, y=value, fill=model)) +
  geom_jitter(aes(x=model, y=value), alpha=0.2) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1L)) +
  ggtitle("Class_Specificity") +
  theme(legend.position = 'none')

# grid.arrange(plot_ac_c,plot_ss_c,plot_sp_c, nrow=2, ncol=2)

plot_ac_c | (plot_ss_c / plot_sp_c)

##################################################################################################################################
### R
# ggsave()
# https://patchwork.data-imaginist.com/articles/guides/layout.html
##################################################################################################################################
### R STUDY

### NORMALIZATION
# https://davetang.org/muse/2011/01/24/normalisation-methods-for-dge-data/
# https://m.blog.naver.com/hanstar_1997/222099507315

### CLUSTERING
# https://rpubs.com/Evan_Jung/hierarchical_clustering