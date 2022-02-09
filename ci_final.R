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
library(factoextra)
library(lattice)
library(caret)
library(devtools)
library(dplyr)
library(ggplot2)
library(ggbiplot)
library(gridExtra)
library(doMC)

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
MGYS00005184_metadata <- MGYS00005184_metadata[, -c('sample_collection date', 'sample_environment (material)')]
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
MGYS00005138_data[which(MGYS00005138_data$host_Age=="Infant")]$host_Age <- 4
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
unique(metadata$Host_disease)

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

##################################################################################################################################
### PCA & tSNE

# Host disease
# Study_arributes.bioproject

# Host sex
# Host Age

tmp_data <- data.frame(metadata$run_accession, metadata$Host_disease, metadata$study_attributes.bioproject, metadata$host_Age)
tmp_data <- rename(tmp_data, 'run_accession' = 'metadata.run_accession')
tmp_data <- rename(tmp_data, 'Host_disease' = 'metadata.Host_disease')
tmp_data <- rename(tmp_data, 'host_Age' = 'metadata.host_Age')
tmp_data <- rename(tmp_data, 'study_code' = 'metadata.study_attributes.bioproject')

##################################################################################################################################

phylum <- tmp_data%>%
  inner_join(data_phylum)
class <- tmp_data %>%
  inner_join(data_class)
order <- tmp_data %>%
  inner_join(data_order)
family <- tmp_data %>%
  inner_join(data_family)
genus <- tmp_data %>%
  inner_join(data_genus)

set.seed(42)

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

plot_pca_pylum <- autoplot(pca_phylum, data=phylum, colour="host_Age", loadings=FALSE, loadings.colour = "black", scale = 0.5)+
  scale_colour_manual(values=c("forestgreen","red","blue")) +
  scale_fill_manual(values=c("forestgreen","red","blue")) +
  scale_shape_manual(values=c(25,22,23))

plot_pca_class <- autoplot(pca_class, data=class, colour="host_Age", loadings=FALSE, loadings.colour = "black", scale = 0.5)+
  scale_colour_manual(values=c("forestgreen","red","blue")) +
  scale_fill_manual(values=c("forestgreen","red","blue")) +
  scale_shape_manual(values=c(25,22,23))

plot_pca_order <- autoplot(pca_order, data=order, colour="host_Age", loadings=FALSE, loadings.colour = "black", scale = 0.5)+
  scale_colour_manual(values=c("forestgreen","red","blue")) +
  scale_fill_manual(values=c("forestgreen","red","blue")) +
  scale_shape_manual(values=c(25,22,23))

plot_pca_family <- autoplot(pca_family, data=family, colour="host_Age", loadings=FALSE, loadings.colour = "black", scale = 0.5)+
  scale_colour_manual(values=c("forestgreen","red","blue")) +
  scale_fill_manual(values=c("forestgreen","red","blue")) +
  scale_shape_manual(values=c(25,22,23))

plot_pca_genus <- autoplot(pca_genus, data=genus, colour="host_Age", loadings=FALSE, loadings.colour = "black", scale = 0.5)+
  scale_colour_manual(values=c("forestgreen","red","blue")) +
  scale_fill_manual(values=c("forestgreen","red","blue")) +
  scale_shape_manual(values=c(25,22,23))

grid.arrange(plot_pca_pylum, plot_pca_class, plot_pca_order, plot_pca_family, plot_pca_genus, nrow=2, ncol=3)

##################################################################################################################################

tsne_phylum <- Rtsne(as.matrix(phylum[,-c(1,2,3,4)]), PCA = FALSE, check_duplicates = FALSE, perplexity=30, theta=0.5, dims=2)
tsne_class <- Rtsne(as.matrix(class[,-c(1,2,3,4)]), PCA = FALSE, check_duplicates = FALSE, perplexity=30, theta=0.5, dims=2)
tsne_order <- Rtsne(as.matrix(order[,-c(1,2,3,4)]), PCA = FALSE, check_duplicates = FALSE, perplexity=30, theta=0.5, dims=2)
tsne_family <- Rtsne(as.matrix(family[,-c(1,2,3,4)]), PCA = FALSE, check_duplicates = FALSE, perplexity=30, theta=0.5, dims=2)
tsne_genus <- Rtsne(as.matrix(genus[,-c(1,2,3,4)]), PCA = FALSE, check_duplicates = FALSE, perplexity=30, theta=0.5, dims=2)

par(mfrow=c(2, 3))

plot(tsne_phylum$Y, col=as.factor(phylum$study_code), pch = 3, cex =0.5)
legend("bottomright", legend = levels(factor(phylum$study_code)), pch = 19, col = factor(levels(factor(phylum$study_code))), cex=0.5)
plot(tsne_class$Y, col=as.factor(class$study_code), pch = 3, cex =0.5)
legend("bottomright", legend = levels(factor(class$study_code)), pch = 19, col = factor(levels(factor(class$study_code))), cex=0.5)
plot(tsne_order$Y, col=as.factor(order$study_code), pch = 3, cex =0.5)
legend("bottomright", legend = levels(factor(order$study_code)), pch = 19, col = factor(levels(factor(order$study_code))), cex=0.5)
plot(tsne_family$Y, col=as.factor(family$study_code), pch = 3, cex =0.5)
legend("bottomright", legend = levels(factor(family$study_code)), pch = 19, col = factor(levels(factor(family$study_code))), cex=0.5)
plot(tsne_genus$Y, col=as.factor(genus$study_code), pch = 3, cex =0.5)
legend("bottomright", legend = levels(factor(ggg$study_code)), pch = 19, col = factor(levels(factor(ggg$study_code))), cex=0.5)

plot(tsne_phylum$Y, col=as.factor(phylum$Host_disease), pch = 3, cex =0.5)
legend("bottomright", legend = levels(factor(phylum$Host_disease)), pch = 19, col = factor(levels(factor(phylum$Host_disease))), cex=0.5)
plot(tsne_class$Y, col=as.factor(class$Host_disease), pch = 3, cex =0.5)
legend("bottomright", legend = levels(factor(class$Host_disease)), pch = 19, col = factor(levels(factor(class$Host_disease))), cex=0.5)
plot(tsne_order$Y, col=as.factor(order$Host_disease), pch = 3, cex =0.5)
legend("bottomright", legend = levels(factor(order$Host_disease)), pch = 19, col = factor(levels(factor(order$Host_disease))), cex=0.5)
plot(tsne_family$Y, col=as.factor(family$Host_disease), pch = 3, cex =0.5)
legend("bottomright", legend = levels(factor(family$Host_disease)), pch = 19, col = factor(levels(factor(family$Host_disease))), cex=0.5)
plot(tsne_genus$Y, col=as.factor(genus$Host_disease), pch = 3, cex =0.5)
legend("bottomright", legend = levels(factor(genus$Host_disease)), pch = 19, col = factor(levels(factor(ggg$Host_disease))), cex=0.5)

plot(tsne_phylum$Y, col=as.factor(phylum$host_Age), pch = 3, cex =0.5)
legend("bottomright", legend = levels(factor(phylum$host_Age)), pch = 19, col = factor(levels(factor(phylum$host_Age))), cex=0.5)
plot(tsne_class$Y, col=as.factor(class$host_Age), pch = 3, cex =0.5)
legend("bottomright", legend = levels(factor(class$host_Age)), pch = 19, col = factor(levels(factor(class$host_Age))), cex=0.5)
plot(tsne_order$Y, col=as.factor(order$host_Age), pch = 3, cex =0.5)
legend("bottomright", legend = levels(factor(order$host_Age)), pch = 19, col = factor(levels(factor(order$host_Age))), cex=0.5)
plot(tsne_family$Y, col=as.factor(family$host_Age), pch = 3, cex =0.5)
legend("bottomright", legend = levels(factor(family$host_Age)), pch = 19, col = factor(levels(factor(family$host_Age))), cex=0.5)
plot(tsne_genus$Y, col=as.factor(genus$host_Age), pch = 3, cex =0.5)
legend("bottomright", legend = levels(factor(genus$host_Age)), pch = 19, col = factor(levels(factor(ggg$host_Age))), cex=0.5)
##################################################################################################################################

# tsne_phylum_Host_disease <- as.data.frame(tsne_phylum_Host_disease$Y)

# ggplot(tsne_phylum_Host_disease, aes(x=V1, y=V2, fill=phylum$Host_disease)) +
#   geom_point(color=factor(levels(factor(family$Host_disease))))

##################################################################################################################################
# 5138 : 1021
# 5184 : 1557
# 5194 : 663

par(mfrow=c(1, 1))
glimpse(metadata)
hc <- hclust(d=dist(data_class), method='average')
plot(hc, main="HCLUSTERING RESULT")

sil_metadata <- c(NA)
sil_phylum <- c(NA)
sil_class <- c(NA)
sil_order <- c(NA)
sil_family <- c(NA)
sil_genus <- c(NA)
for(idx in 2:10){
  pam_fit <- pam(dist(metadata), diss=TRUE, k=idx)
  pam_fit_phylum <- pam(dist(data_phylum), diss=TRUE, k=idx)
  pam_fit_class <- pam(dist(data_class), diss=TRUE, k=idx)
  pam_fit_order <- pam(dist(data_order), diss=TRUE, k=idx)
  pam_fit_family <- pam(dist(data_family), diss=TRUE, k=idx)
  pam_fit_genus <- pam(dist(data_genus), diss=TRUE, k=idx)
  sil_metadata[idx] <- pam_fit$silinfo$avg.width
  sil_phylum[idx] <- pam_fit_phylum$silinfo$avg.width
  sil_class[idx] <- pam_fit_class$silinfo$avg.width
  sil_order[idx] <- pam_fit_order$silinfo$avg.width
  sil_family[idx] <- pam_fit_family$silinfo$avg.width
  sil_genus[idx] <- pam_fit_genus$silinfo$avg.width
  print(idx)
}

par(mfrow=c(2, 3))
plot(1:10, sil_metadata, xlab = "Number of clusters", ylab = "Silhouette Width")
lines(1:10, sil_metadata)
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

dummy <- dummyVars(" ~ .", data=metadata)
newdata <- data.frame(predict(dummy, newdata=metadata))

cutree(hc,k=3)
silhouette(hc)
fviz_nbclust(metadata, kmeans, method = "silhouette")

plot(hc, k=3)
rect.hclust(hc,k=3,border="red")
plot(silhouette(cutree(hc,k=20),dist=dist(metadata),col=1:20))

##################################################################################################################################
registerDoMC(cores = 30)
getDoParWorkers()

train_phylum <- phylum[,-c(1,3,4)]
train_class <- class[,-c(1,3,4)]
train_order <- order[,-c(1,3,4)]
train_family <- family[,-c(1,3,4)]
train_genus <- genus[,-c(1,3,4)]

MySummary  <- function(data, lev = NULL, model = NULL){
  a1 <- defaultSummary(data, lev, model)
  b1 <- multiClassSummary(data, lev, model) 
  out <- c(a1,b1)
  out}

fitControl <- trainControl(method = "repeatedcv", number = 10, repeats = 3, savePredictions = TRUE, summaryFunction=MySummary)

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

ggplot(accuracy_phylum) +
  geom_boxplot(aes(x=model, y=value, fill=model)) +
  geom_jitter(aes(x=model, y=value), alpha=0.2) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1L)) +
  ggtitle("Phylum_Accuracy")

ggplot(sensitivity_phylum) +
  geom_boxplot(aes(x=model, y=value, fill=model)) +
  geom_jitter(aes(x=model, y=value), alpha=0.2) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1L)) +
  ggtitle("Phylum_Sensitivity")

ggplot(specificity_phylum) +
  geom_boxplot(aes(x=model, y=value, fill=model)) +
  geom_jitter(aes(x=model, y=value), alpha=0.2) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1L)) +
  ggtitle("Phylum_Specificity")
##################################################################################################################################
# https://rpubs.com/Evan_Jung/hierarchical_clustering