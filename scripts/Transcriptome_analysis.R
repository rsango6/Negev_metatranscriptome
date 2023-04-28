library(DESeq2)
library(data.table)
library(xlsx)
library(dplyr)
library(stringr)
library(ggplot2)
library(ggrepel)
library(RColorBrewer)
library(pheatmap)
library(vegan)

options(scipen = 8)


All_Feature_counts <- read.csv("Negev_featureCounts.tsv", 
                               sep="\t", 
                               header=T, 
                               quote="",
                               check.names=F)
All_Feature_counts <- All_Feature_counts[,c(1,2,6,7:ncol(All_Feature_counts))]

Metadata = read.xlsx('Metadata.xlsx', sheetIndex = 1)

#Copying C_0 and E_0 to be C_0_dark and E_0_dark, respectively

All_Feature_counts$C_0_dark = All_Feature_counts$C_0
All_Feature_counts$E_0_dark = All_Feature_counts$E_0

colnames(All_Feature_counts)[4:ncol(All_Feature_counts)] = Metadata$X2_sample_name_short


#Fragments per kbp
FPK <- (All_Feature_counts[,4:ncol(All_Feature_counts)]) / (All_Feature_counts[,3]/1000)

#Transcripts per kbp per million
TPM <- sweep(FPK,2, colSums(FPK)/1000000, "/")

TPM$Geneid <- All_Feature_counts$Geneid

Annotation <- read.csv("Additional_data/Avdat_crust_all_annotations.tsv", 
                       header=T, 
                       sep="\t", 
                       quote="",
                       check.names=F)

write.table(merge(Annotation, TPM, by="Geneid", all.y=T), 
            file="Results/Negev_featureCounts_TPM.tsv", 
            sep="\t", 
            row.names=F, 
            quote=F)

#generate subtotals for columns 5 - 32. function to be applied is "sum"

# Transcripts_per_bin <- cube(merge(Annotation[,c("Bin","Geneid")], TPM, by="Geneid", all.y=T), 
#                             lapply(.SD, sum), 
#                             by="Bin", 
#                             .SDcol=c(3:26))

TPM_Annotation_Combined = merge(Annotation[,c("Bin","Geneid")], TPM, by="Geneid", all.y=T)

Transcripts_per_bin <- TPM_Annotation_Combined %>%
  select(-Geneid) %>%
  group_by(Bin) %>% 
  summarise(across(everything(), sum),
            .groups = 'drop')  %>%
  as.data.frame()

#generate subtotals for columns 5 - 32. function to be applied is "count non-zeroes"
# Transcribed_genes_per_bin <- cube(merge(Annotation[,c("Bin","Geneid")], TPM, by="Geneid", all.y=T), 
#                                         lapply(.SD, function(c)sum(c!=0)), 
#                                   by="Bin", 
#                                   .SDcol=c(3:26))

rowSums(Transcripts_per_bin[-1]) #removes first column from computation


write.csv(Transcripts_per_bin, file="Results/Negev_TPMs_per_bin.csv", row.names=F, quote=F)
#write.table(Transcribed_genes_per_bin, file="Avdat_transcribed_genes_per_bin.tsv", sep="\t", row.names=F, quote=F)


#---------------------------------------------------------------------------------------------------
#Generating TPM values normalized per bin
#----------------------------------------------------------------------------------------------------

#system("bash Subset_featureCounts_table.sh", wait=T)

Bin_list = unique(Transcripts_per_bin$Bin)

All_Feature_counts$Bin = sub("^[^_]*-", "", All_Feature_counts$Chr) #remove everything before a '-'
All_Feature_counts$Bin = sub("_[^_]+$", "", All_Feature_counts$Bin) #remove everything after a '_'


All_Feature_counts = All_Feature_counts %>%
  relocate(Bin, .after = Chr)

for(i in unique(All_Feature_counts$Bin)) {
  
  #Importing the feature Counts
  #countsfile <- file.path(paste("./Transcripts_per_bin/", i, "_featureCounts.tsv", sep=""))
  # FeatureCounts <- read.table(countsfile, 
  #                             sep="\t", 
  #                             header=F, 
  #                             quote="")
  
  FeatureCountsPerBin = All_Feature_counts %>%
    filter(Bin == i)
 
  colnames(FeatureCountsPerBin) <- colnames(All_Feature_counts)
  FPK <- (FeatureCountsPerBin[,5:ncol(FeatureCountsPerBin)]) / (FeatureCountsPerBin[,4]/1000)
  TPM <- sweep(FPK,2, colSums(FPK)/1000000, "/")
  TPM$Geneid <- FeatureCountsPerBin$Geneid
  TPM$Chr <- FeatureCountsPerBin$Chr
  TPM = TPM %>%
    relocate(Chr, .before = C_0)

  output <- file.path(paste("Results/Transcripts_per_bin/", i, "_featureCounts_TPM_PerBin.tsv", sep=""))
  
  write.table(
    merge(Annotation, TPM, by="Geneid", all.y=T),
    file=output,
    sep = "\t",
    row.names=F, 
    quote=F)
}


# hierarchical clustering on per-MAG-normalized-gene-expression
# on 1) only light dataset, on 2) dark dataset and 3) on everything combined

LightDf = Transcripts_per_bin[, !grepl("dark|Bin",colnames(Transcripts_per_bin))]

DarkDf = Transcripts_per_bin[, grepl("dark",colnames(Transcripts_per_bin))]

PerMAGlist = list(LightDf = LightDf, 
                  DarkDf = DarkDf, 
                  LightDarkCombined = Transcripts_per_bin[-1])


for (i in 1:length(PerMAGlist)) {
  
  dist <- vegdist(t(PerMAGlist[[i]]), method="jaccard")
  cluster <- hclust(dist, method="average")
  plot(cluster, hang="-1", ylab="Dissimilarity", main = NULL, xlab = names(PerMAGlist[i]))
}

for (i in 1:length(PerMAGlist)) {
  
  dist <- vegdist(t(PerMAGlist[[i]]), method="jaccard")
  
  MDS <- metaMDS(dist, engine = c("monoMDS"))
  
  MDSdf <- data.frame(MDS1=MDS$points[,1], 
                      MDS2=MDS$points[,2])
  
  MDSdf$Phase=Metadata$X7_phase[match(rownames(MDSdf), Metadata$X2_sample_name_short)]
  
  print(ggplot(MDSdf, aes(MDS1, MDS2, label = rownames(MDSdf))) +
          theme_bw() +
          geom_text_repel(aes(color = Phase)) +
          scale_color_manual(values=c("dry" = "#D95F02", 
                                     "hydrated" = "#1B9E77")) +
          labs(title = names(PerMAGlist[i]))
  )
}

for (i in 1:length(PerMAGlist)) {

  dist <- vegdist(t(PerMAGlist[[i]]), method="jaccard")
  sampleDistMatrix <- as.matrix(dist)
  rownames(sampleDistMatrix) <- colnames(PerMAGlist[[i]])
  colnames(sampleDistMatrix) <-  colnames(PerMAGlist[[i]])
  
  if (names(PerMAGlist)[i] == "DarkDf" | names(PerMAGlist)[i] == "LightDarkCombined") {
    PerMAGlist[[i]] =  PerMAGlist[[i]] %>%
      dplyr::relocate(C_0_dark, .before = H_0_dark) %>%
      dplyr::relocate(E_0_dark, .before = H_0_dark)
  }

  colors <- colorRampPalette(brewer.pal(8, "Greens"))(255)
  PhaseCols = Metadata$X7_phase[match(colnames(PerMAGlist[[i]]), Metadata$X2_sample_name_short)]
  anno <- data.frame(Timepoint=colnames(PerMAGlist[[i]]),
                      Phase = PhaseCols)

  rownames(anno) = colnames(sampleDistMatrix)

  annoCol <- colorRampPalette(brewer.pal(8, "YlOrRd"))(length(colnames(PerMAGlist[[i]])))
  names(annoCol) <- unique(anno$Timepoint)

  annoCol <- list(Timepoint = annoCol,
                Phase = c(dry = "#D95F02",
                          hydrated = "#1B9E77")
                )
  
  pheatmap(sampleDistMatrix, 
         clustering_distance_rows = dist,
         clustering_distance_cols = dist,
         col = colors,
         annotation_col = anno,
         annotation_colors = annoCol,
         main = paste(names(PerMAGlist[i]), "distances")
  )
  
}

#---------------------------------------------------------------------------------------------------------
#Differential expression analysis
#---------------------------------------------------------------------------------------------------------

#---------------------------------------------------------------------------------------------------------------------------------
#Analyzing differential expression between e.g. 0 min and 15 min timepoints of the rehydration experiment
#---------------------------------------------------------------------------------------------------------------------------------
#DeSeq2 compares one category vs another. Therefore a subsetting of the Metadata and Featurecounts tables is necessary.
#Importing metadata (subsetting to compare 0 min and 15 min timepoints)
rehydration_metadata_sub <- read.table("Metadata.csv", 
                                     sep="\t", 
                                     header=T, 
                                     row.names=1, 
                                     quote="")[c(1:6),]
#Bin_list_deseq.txt is the list of bins fulfilling the criteria for Deseq2 analysis. Otherwise, Deseq2 will abort.

for (i in read.table("Bin_list_deseq.txt")[,1]) {
  
  #Importing the feature Counts
  countsfile <- file.path(paste("./Transcripts_per_bin/", i, "_featureCounts.tsv", sep=""))
  
  FeatureCounts <- read.table(countsfile, 
                              row.names=1, 
                              sep="\t", 
                              header=F, 
                              quote="")[,c(4:9)]
  
  colnames(FeatureCounts) <- c("C_0","E_0","F_0","C_0.25","E_0.25","F_0.25")
  
  #The TPM table for the MAG will later be merged with the Deseq results, so that transcription profiles of significantly diff. expressed genes can be checked on the spot
  TPM_file <- file.path(paste("Transcripts_per_bin/", i, "_featureCounts_TPMperbin.tsv", sep=""))
  TPM_table <- read.csv(TPM_file, 
                        sep="\t", 
                        header=T, 
                        quote="")
  
  #Subsetting the count table in order to compare two conditions against each other. 
  #In this case 0 min and 15 min after rehydration
  
  FeatureCounts_sub <- FeatureCounts[,1:6]
  
  dds <- DESeqDataSetFromMatrix(countData = FeatureCounts_sub,
                                colData = rehydration_metadata_sub,
                                design= ~ Timepoint)
  
  dds <- DESeq(dds)
  #resultsNames(dds) # lists the coefficients
  res <- as.data.frame(results(dds, name="Timepoint_T00.25_vs_T00"))
  setDT(res, keep.rownames = "Geneid")
  res <- transform(res, Geneid = as.numeric(Geneid))
  
  Deseq_results <- merge(res, TPM_table, 
          by="Geneid", 
          all.x=T)
  
  output <- file.path(paste("./Diff_expression/", i, "_00.25_00.tsv", sep=""))
  
  write.table(Deseq_results, 
              file=output, 
              sep="\t",
              row.names=F,
              quote=F)
  
}


#---------------------------------------------------------------------------------------------------------------------------------
#Analyzing differential expression between all "dry" and all "hydrated" timepoints of the rehydration experiment
#---------------------------------------------------------------------------------------------------------------------------------

#Importing metadata (subsetting to compare all dry and all hydrated timepoints)
rehydration_metadata_sub <- read.table("Metadata.csv", 
                                     sep="\t", 
                                     header=T, 
                                     row.names=1, 
                                     quote="")[c(1:18,20:24),]

for(i in read.table("Bin_list_deseq.txt")[,1]) {
  
  #Importing the feature Counts
  countsfile <- file.path(paste("./Transcripts_per_bin/", i, "_featureCounts.tsv", sep=""))
  
  FeatureCounts <- read.table(countsfile, 
                              row.names=1, 
                              sep="\t", 
                              header=F, 
                                quote="")[,c(3:20,22:26)]
  
  colnames(FeatureCounts) <- c("C_0","E_0","F_0","C_0.25","E_0.25","F_0.25","C_0.5","E_0.5","F_0.5","C_3","E_3","F_3","C_6","E_6","F_6","C_12","E_12","F_12","E_39","F_39","C_55","E_55","F_55")
  
  
  TPM_file <- file.path(paste("Transcripts_per_bin/", i, "_featureCounts_TPMperbin.tsv", sep=""))
  TPM_table <- read.csv(TPM_file, 
                        sep="\t", 
                        header=T, 
                        quote="")
  
  
  #Subsetting the count table in order to compare two conditions against each other. 
  #In this case all hydrated vs all dry
  
  FeatureCounts_sub <- FeatureCounts[,1:23]
  
  dds <- DESeqDataSetFromMatrix(countData = FeatureCounts_sub,
                                colData = rehydration_metadata_sub,
                                design= ~ Phase1)
  
  dds <- DESeq(dds)
  #resultsNames(dds) # lists the coefficients
  res <- as.data.frame(results(dds, name="Phase1_hydrated_vs_dry"))
  setDT(res, keep.rownames = "Geneid")
  res <- transform(res, Geneid = as.numeric(Geneid))
  
  #Deseq_results <- merge(res, TPM_table, 
  #                       by="Geneid", 
  #                       all.x=T)
  
  output <- file.path(paste("./Diff_expression/", i, "_allhydrated_alldry.tsv", sep=""))
  
  write.table(res, 
              file=output, 
              sep="\t",
              row.names=F,
              quote=F)
  
}


#---------------------------------------------------------------------------------------------------------------------------------
#Analyzing differential expression between "dry" and "early hydrated" timepoints of the rehydration experiment
#---------------------------------------------------------------------------------------------------------------------------------

#Importing metadata (subsetting to compare eraly hydrated and dry timepoints)
rehydration_metadata_sub <- read.table("Metadata.csv", 
                                     sep="\t", 
                                     header=T, 
                                     row.names=1, 
                                     quote="")[c(1:9,20:24),]

for(i in read.table("Bin_list_deseq.txt")[,1]) {
  
  #Importing the feature Counts
  countsfile <- file.path(paste("./Transcripts_per_bin/", i, "_featureCounts.tsv", sep=""))
  
  FeatureCounts <- read.table(countsfile, 
                              row.names=1, 
                              sep="\t", 
                              header=F, 
                              quote="")[,c(3:20,22:26)]
  
  colnames(FeatureCounts) <- c("C_0","E_0","F_0","C_0.25","E_0.25","F_0.25","C_0.5","E_0.5","F_0.5","C_3","E_3","F_3","C_6","E_6","F_6","C_12","E_12","F_12","E_39","F_39","C_55","E_55","F_55")
  
  
  TPM_file <- file.path(paste("Transcripts_per_bin/", i, "_featureCounts_TPMperbin.tsv", sep=""))
  TPM_table <- read.csv(TPM_file, 
                        sep="\t", 
                        header=T, 
                        quote="")
  
  
  #Subsetting the count table in order to compare two conditions against each other. 
  #In this case dry and earlyhydrated
  
  FeatureCounts_sub <- FeatureCounts[,c(1:9,19:23)]
  
  
  dds <- DESeqDataSetFromMatrix(countData = FeatureCounts_sub,
                                colData = rehydration_metadata_sub,
                                design= ~ Phase2)
  
  dds <- DESeq(dds)
  #resultsNames(dds) # lists the coefficients
  res <- as.data.frame(results(dds, name="Phase2_earlyhydrated_vs_dry"))
  setDT(res, keep.rownames = "Geneid")
  res <- transform(res, Geneid = as.numeric(Geneid))
  
  #Deseq_results <- merge(res, TPM_table, 
  #                       by="Geneid", 
  #                       all.x=T)
  
  output <- file.path(paste("./Diff_expression/", i, "_earlyhydrated_dry.tsv", sep=""))
  
  write.table(res, 
              file=output, 
              sep="\t",
              row.names=F,
              quote=F)
  
}

#---------------------------------------------------------------------------------------------------------------------------------
#Analyzing differential expression between "early hydrated" and "hydrated" timepoints of the rehydration experiment
#---------------------------------------------------------------------------------------------------------------------------------

#Importing metadata (subsetting to compare early hydrated and hydrated timepoints)
rehydration_metadata_sub <- read.table("Metadata.csv", 
                                     sep="\t", 
                                     header=T, 
                                     row.names=1, 
                                     quote="")[c(4:9,10:18),]

for(i in read.table("Bin_list_deseq.txt")[,1]) {
  
  #Importing the feature Counts
  countsfile <- file.path(paste("./Transcripts_per_bin/", i, "_featureCounts.tsv", sep=""))
  
  FeatureCounts <- read.table(countsfile, 
                              row.names=1, 
                              sep="\t", 
                              header=F, 
                              quote="")[,c(3:20,22:26)]
  
  colnames(FeatureCounts) <- c("C_0","E_0","F_0","C_0.25","E_0.25","F_0.25","C_0.5","E_0.5","F_0.5","C_3","E_3","F_3","C_6","E_6","F_6","C_12","E_12","F_12","E_39","F_39","C_55","E_55","F_55")
  
  
  TPM_file <- file.path(paste("Transcripts_per_bin/", i, "_featureCounts_TPMperbin.tsv", sep=""))
  TPM_table <- read.csv(TPM_file, 
                        sep="\t", 
                        header=T, 
                        quote="")
  
  
  #Subsetting the count table in order to compare two conditions against each other. 
  #In this case hydrated and earlyhydrated
  
  FeatureCounts_sub <- FeatureCounts[,c(4:9,10:18)]
  
  
  dds <- DESeqDataSetFromMatrix(countData = FeatureCounts_sub,
                                colData = rehydration_metadata_sub,
                                design= ~ Phase2)
  
  dds <- DESeq(dds)
  #resultsNames(dds) # lists the coefficients
  res <- as.data.frame(results(dds, name="Phase2_hydrated_vs_earlyhydrated"))
  setDT(res, keep.rownames = "Geneid")
  res <- transform(res, Geneid = as.numeric(Geneid))
  
  #Deseq_results <- merge(res, TPM_table, 
  #                       by="Geneid", 
  #                       all.x=T)
  
  output <- file.path(paste("./Diff_expression/", i, "_hydrated_earlyhydrated.tsv", sep=""))
  
  write.table(res, 
              file=output, 
              sep="\t",
              row.names=F,
              quote=F)
  
}

#---------------------------------------------------------------------------------------------------------------------------------
#Analyzing differential expression between "hydrated" and "dry" timepoints of the rehydration experiment
#---------------------------------------------------------------------------------------------------------------------------------

#Importing metadata (subsetting to compare hydrated and dry timepoints)
rehydration_metadata_sub <- read.table("Metadata.csv", 
                                     sep="\t", 
                                     header=T, 
                                     row.names=1, 
                                     quote="")[c(1:3,10:18,20:24),]

for(i in read.table("Bin_list_deseq.txt")[,1]) {
  
  #Importing the feature Counts
  countsfile <- file.path(paste("./Transcripts_per_bin/", i, "_featureCounts.tsv", sep=""))
  
  FeatureCounts <- read.table(countsfile, 
                              row.names=1, 
                              sep="\t", 
                              header=F, 
                              quote="")[,c(3:20,22:26)]
  
  colnames(FeatureCounts) <- c("C_0","E_0","F_0","C_0.25","E_0.25","F_0.25","C_0.5","E_0.5","F_0.5","C_3","E_3","F_3","C_6","E_6","F_6","C_12","E_12","F_12","E_39","F_39","C_55","E_55","F_55")
  
  
  TPM_file <- file.path(paste("Transcripts_per_bin/", i, "_featureCounts_TPMperbin.tsv", sep=""))
  TPM_table <- read.csv(TPM_file, 
                        sep="\t", 
                        header=T, 
                        quote="")
  
  
  #Subsetting the count table in order to compare two conditions against each other. 
  #In this case hydrated and dry
  
  FeatureCounts_sub <- FeatureCounts[,c(1:3,10:23)]
  
  
  dds <- DESeqDataSetFromMatrix(countData = FeatureCounts_sub,
                                colData = rehydration_metadata_sub,
                                design= ~ Phase2)
  
  dds <- DESeq(dds)
  #resultsNames(dds) # lists the coefficients
  res <- as.data.frame(results(dds, name="Phase2_hydrated_vs_dry"))
  setDT(res, keep.rownames = "Geneid")
  res <- transform(res, Geneid = as.numeric(Geneid))
  
  #Deseq_results <- merge(res, TPM_table, 
  #                       by="Geneid", 
  #                       all.x=T)
  
  output <- file.path(paste("./Diff_expression/", i, "_hydrated_dry.tsv", sep=""))
  
  write.table(res, 
              file=output, 
              sep="\t",
              row.names=F,
              quote=F)
  
}

    

#---------------------------------------------------------------------------------------------------------------------------------
# Combining tables into one giant table
#---------------------------------------------------------------------------------------------------------------------------------
system("tail -n +2 -q Diff_expression/*_00.25_00.tsv > Diff_expression/Diff_expression_00.25_00.tsv", wait=T)
Log_pvalue_0_0.25 <- read.csv("Diff_expression/Diff_expression_00.25_00.tsv", 
                              header=T, 
                              sep="\t",
                              quote="",
                              check.names = F)[,c(1,3,7)]
colnames(Log_pvalue_0_0.25) <- c("Geneid", "Log2chnage_00_0.25", "padj_00_0.25")

system("tail -n +2 -q Diff_expression/*_00.5_00.25.tsv > Diff_expression/Diff_expression_00.5_00.25.tsv", wait=T)
Log_pvalue_0.25_0.5 <- read.csv("Diff_expression/Diff_expression_00.5_00.25.tsv", 
                                header=T, 
                                sep="\t",
                                quote="",
                                check.names = F)[,c(1,3,7)]
colnames(Log_pvalue_0.25_0.5) <- c("Geneid", "Log2chnage_0.25_0.5", "padj_0.25_0.5")

system("tail -n +2 -q Diff_expression/*_03_00.5.tsv > Diff_expression/Diff_expression_03_00.5.tsv", wait=T)

Log_pvalue_0.5_03 <- read.csv("Diff_expression/Diff_expression_03_00.5.tsv", 
                              header=T, 
                              sep="\t",
                              quote="",
                              check.names = F)[,c(1,3,7)]
colnames(Log_pvalue_0.5_03) <- c("Geneid", "Log2chnage_0.5_03", "padj_0.5_03")

system("tail -n +2 -q Diff_expression/*_06_03.tsv > Diff_expression/Diff_expression_06_03.tsv", wait=T)
Log_pvalue_03_06 <- read.csv("Diff_expression/Diff_expression_06_03.tsv", 
                             header=T, 
                             sep="\t",
                             quote="",
                             check.names = F)[,c(1,3,7)]
colnames(Log_pvalue_03_06) <- c("Geneid", "Log2chnage_03_06", "padj_03_06")

system("tail -n +2 -q Diff_expression/*_12_06.tsv > Diff_expression/Diff_expression_12_06.tsv", wait=T)
Log_pvalue_06_12 <- read.csv("Diff_expression/Diff_expression_12_06.tsv", 
                             header=T, 
                             sep="\t",
                             quote="",
                             check.names = F)[,c(1,3,7)]
colnames(Log_pvalue_06_12) <- c("Geneid", "Log2chnage_06_12", "padj_06_12")

system("tail -n +2 -q Diff_expression/*_39_12.tsv > Diff_expression/Diff_expression_39_12.tsv", wait=T)
Log_pvalue_12_39 <- read.csv("Diff_expression/Diff_expression_39_12.tsv", 
                             header=T, 
                             sep="\t",
                             quote="",
                             check.names = F)[,c(1,3,7)]
colnames(Log_pvalue_12_39) <- c("Geneid", "Log2chnage_12_39", "padj_12_39")

system("tail -n +2 -q Diff_expression/*_55_39.tsv > Diff_expression/Diff_expression_55_39.tsv", wait=T)
Log_pvalue_39_55 <- read.csv("Diff_expression/Diff_expression_55_39.tsv", 
                             header=T, 
                             sep="\t",
                             quote="",
                             check.names = F)[,c(1,3,7)]
colnames(Log_pvalue_39_55) <- c("Geneid", "Log2chnage_39_55", "padj_39_55")

system("tail -n +2 -q Diff_expression/*_allhydrated_vs_alldry.tsv > Diff_expression/Diff_expression_allhydrated_vs_alldry.tsv", wait=T)
Log_pvalue_allhydrated_vs_alldry <- read.csv("Diff_expression/Diff_expression_allhydrated_vs_alldry.tsv", 
                              header=T, 
                              sep="\t",
                              quote="",
                              check.names = F)[,c(1,3,7)]
colnames(Log_pvalue_allhydrated_vs_alldry) <- c("Geneid", "Log2chnage_allhydrated_vs_alldry", "padj_allhydrated_vs_alldry")

system("tail -n +2 -q Diff_expression/*_earlyhydrated_vs_dry.tsv > Diff_expression/Diff_expression_earlyhydrated_vs_dry.tsv", wait=T)
Log_pvalue_earlyhydrated_vs_dry <- read.csv("Diff_expression/Diff_expression_earlyhydrated_vs_dry.tsv", 
                                             header=T, 
                                             sep="\t",
                                             quote="",
                                             check.names = F)[,c(1,3,7)]
colnames(Log_pvalue_earlyhydrated_vs_dry) <- c("Geneid", "Log2chnage_earlyhydrated_vs_dry", "padj_earlyhydrated_vs_dry")

system("tail -n +2 -q Diff_expression/*_hydrated_vs_earlyhydrated.tsv > Diff_expression/Diff_expression_hydrated_vs_earlyhydrated.tsv", wait=T)
Log_pvalue_hydrated_vs_earlyhydrated <- read.csv("Diff_expression/Diff_expression_hydrated_vs_earlyhydrated.tsv", 
                                            header=T, 
                                            sep="\t",
                                            quote="",
                                            check.names = F)[,c(1,3,7)]
colnames(Log_pvalue_hydrated_vs_earlyhydrated) <- c("Geneid", "Log2chnage_hydrated_vs_earlyhydrated", "padj_hydrated_vs_earlyhydrated")

system("tail -n +2 -q Diff_expression/*_hydrated_vs_dry.tsv > Diff_expression/Diff_expression_hydrated_vs_dry.tsv", wait=T)
Log_pvalue_hydrated_vs_dry <- read.csv("Diff_expression/Diff_expression_hydrated_vs_dry.tsv", 
                                                 header=T, 
                                                 sep="\t",
                                                 quote="",
                                                 check.names = F)[,c(1,3,7)]
colnames(Log_pvalue_hydrated_vs_dry) <- c("Geneid", "Log2chnage_hydrated_vs_dry", "padj_hydrated_vs_dry")

Log_and_pvalues_all <- Reduce(function(x, y) 
  merge(x, y, by="Geneid", all=T), 
  list(Log_pvalue_allhydrated_vs_alldry, Log_pvalue_earlyhydrated_vs_dry, Log_pvalue_hydrated_vs_earlyhydrated, Log_pvalue_hydrated_vs_dry, Log_pvalue_0_0.25, Log_pvalue_0.25_0.5, Log_pvalue_0.5_03, Log_pvalue_03_06, Log_pvalue_06_12, Log_pvalue_12_39, Log_pvalue_39_55)
)

#calculating the smallest p_value. Important to quickly find genes that were significantly up- or down-regulated at lest at some transition.
#na.rm=T is important, because otherwise "NA" will always be the smallest value if it only appears once in a row.

Log_and_pvalues_all <- transform(Log_and_pvalues_all, 
                                 Min_pvalue = pmin(padj_allhydrated_vs_alldry, padj_earlyhydrated_vs_dry, padj_hydrated_vs_earlyhydrated, padj_hydrated_vs_dry, padj_00_0.25, padj_0.25_0.5, padj_0.5_03, padj_03_06, padj_06_12, padj_12_39, padj_39_55, 
                                                   na.rm = T))

write.table(Log_and_pvalues_all, file="Diff_expression/Differential_expression_Log2_and_padj.tsv", 
            sep="\t", 
            quote=F, 
            row.names=F)

All_TPM_annotated <- read.csv("Transcripts_per_bin/Annotated_TPM_perbin_all.tsv", 
                              header=T, 
                              sep="\t",
                              quote="",
                              check.names = F)

ALL_COMBINED <- merge(All_TPM_annotated, Log_and_pvalues_all, 
                      by="Geneid", 
                      all=T)

write.table(ALL_COMBINED, file="Avdat_per_bin_normalized_TPMs_DeSeq2.tsv", 
            sep="\t", 
            quote=F, 
            row.names=F)

#---------------------------------------------------------------------------------------------------------------------------------
# Calculating differentially expressed genes per MAG per timepoint
#---------------------------------------------------------------------------------------------------------------------------------

Avdat_pvalues <- data.table(ALL_COMBINED[,c(1,36,38,40,42,44,46,48,49)])

#replacing missing pvalues with a high number
Avdat_pvalues[is.na(Avdat_pvalues)] <- 100

Diff_expr_perMAG_per_time <- cube(Avdat_pvalues,
                                  lapply(.SD,
                                         function(c)sum(c<=0.05)), 
                                  by="MAG", 
                                  .SDcol=c(2:8))

write.table(Diff_expr_perMAG_per_time, file="Avdat_Diff_expr_perMAG_per_time.tsv", 
            sep="\t", 
            quote=F, 
            row.names=F)