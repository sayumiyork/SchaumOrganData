########################################
# Links
########################################

#Article
#https://pubmed.ncbi.nlm.nih.gov/32669715/ 

#Figshare
#https://figshare.com/projects/The_murine_transcriptome_reveals_global_aging_nodes_with_organ-specific_phase_and_amplitude/65126 

#GSE132040
#https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE132040

########################################
# Libraries
########################################

library(tidyverse)
library(DESeq2)

########################################
# Functions
########################################


formatDESeq2Results <- function( x ) {
  df <- as.data.frame(x)
  df <- data.frame(rownames(df), df)
  colnames(df) <- c("GeneID", colnames(df)[-1])
  rownames(df) <- c()
  return(df)
}

#Fly plotAcrossRegions example
plotAcrossRegions <- function( x ) {
  df <- data.frame(counts(midgut, normalized = TRUE)[x,], midgut_tsv$condition)
  colnames(df) <- c("counts", "region")
  df <- df[10:30,]
  df$region <- fct_relevel( factor( df$region ), "a1", "a2_3", "Cu", "LFCFe", "Fe", "p1", "p2_4" )
  ggplot( df ) +
    geom_bar( aes( region, counts ), stat="identity" ) +    
    ggtitle(x)+
    theme_bw() + 
    theme( axis.text.x = element_text( angle=90, vjust=0.5 ) )
}

#Mouse plotAcrossHostDonor example
plotAcrossHostDonor<- function( x ) {
  df <- data.frame(counts(brain, normalized=TRUE)[x,], brain_metadata$host_donor_tissue, row.names = row.names(brain_metadata))
  colnames(df) <- c("counts", "host_donor_tissue")
  df$host_donor_tissue <- fct_relevel( factor( df$host_donor_tissue), "TD_1_STR","TD_2_STR","TD_3_STR","TD_1_PFC", "TD_2_PFC", "TD_3_PFC")
  ggplot( df ) +
    geom_bar( aes( host_donor_tissue, counts ), stat="identity" ) +
    ggtitle(x)+
    theme_bw() +
    theme( axis.text.x = element_text( angle=90, vjust=0.5 ) )
}

plotAcrossAgeKidney<- function( x ) {
  df <- data.frame(counts(Schaum_DESeq2_Results_Kidney, normalized=TRUE)[x,], schaum_metadata_kidney$age, row.names = row.names(schaum_metadata_kidney))
  colnames(df) <- c("counts", "age")
  df$age <- fct_relevel( factor( df$age), "1","3","6","9","12", "15", "18", "21", "24", "27")
  ggplot( df ) +
    geom_bar( aes( age, counts ), stat="identity" ) +
    ggtitle(x)+
    xlab("Months") + 
    ylab("Normalized counts") +
    theme_bw() +
    theme( axis.text.x = element_text( angle=90, vjust=0.5 ) )
}


plotAcrossSexKidney<- function( x ) {
  df <- data.frame(counts(Schaum_DESeq2_Results_Kidney, normalized=TRUE)[x,], schaum_metadata_kidney$sex, row.names = row.names(schaum_metadata_kidney))
  colnames(df) <- c("counts", "sex")
  df$age <- fct_relevel( factor( df$sex), "m","f")
  ggplot( df ) +
    geom_bar( aes( sex, counts ), stat="identity" ) +
    ggtitle(x)+
    xlab("Sex") + 
    ylab("Normalized counts") +
    theme_bw() +
    theme( axis.text.x = element_text( angle=90, vjust=0.5 ) )
}


plotAcrossAgeSexKidney<- function( x ) {
  df <- data.frame(counts(Schaum_DESeq2_Results_Kidney, normalized=TRUE)[x,], schaum_metadata_kidney$age_sex, row.names = row.names(schaum_metadata_kidney))
  colnames(df) <- c("counts", "age_sex")
  df$age_sex <- fct_relevel( factor( df$age_sex), "1_m","3_m","6_m","9_m","12_m", "15_m", "18_m", "21_m", "24_m", "27_m", "1_f","3_f","6_f","9_f","12_f", "15_f", "18_f", "21_f")
  ggplot( df ) +
    geom_bar( aes( age_sex, counts ), stat="identity" ) +
    ggtitle(x)+
    xlab("Mice") + 
    ylab("Normalized counts") +
    theme_bw() +
    theme( axis.text.x = element_text( angle=90, vjust=0.5 ) )
}




########################################
# Approach 1: Recreate DESeq2 Object
########################################


# Load in raw counts table from figshare
schaum_counts <- read_tsv("Schaum_Raw_Counts.tsv")
schaum_counts <- column_to_rownames(schaum_counts, var = "...1")

# Load in metadata table from figshar
# Additional metadata columns can/have been added 
schaum_metadata <- read.delim("Schaum_Metadata.tsv", row.names = 1)

# Create DESeq2 object with comparision formula:
# To identify significant differential expression changes with age, we used the raw count 
# matrix as recommended for Deseq2’s standard analysis pipeline. Factors and dispersion
# estimates were calculated for each tissue separately. We conducted differential expression 
# analysis comparing samples from 3 months old mice to each consecutive time point, using 
# age and sex as covariates. P-values were adjusted for multiple testing and genes with an 
# adjusted p-value of below 0.05 were determined as statistically significant. In addition, we 
# ran similar analyses using 1 month or 6 months old mice as reference.



# Example subset code
# Select specific samples by name
keep_samples <- c("Sample1", "Sample3", "Sample5")
dds_subset <- dds[ , keep_samples]
# Subset to keep only 'Control' and 'TreatmentA'
dds_subset <- dds[ , dds$condition %in% c("Control", "TreatmentA")]


#Bone test
# Load in raw counts table from figshare
schaum_counts_bones <- read_tsv("Schaum_Raw_Counts_Bones.tsv")
schaum_counts_bones <- column_to_rownames(schaum_counts_bones, var = "...1")
# Load in metadata table from figshar
# Additional metadata columns can/have been added 
schaum_metadata_bones <- read.delim("Schaum_Metadata_Bones.tsv", row.names = 1)
dds_bones <- DESeqDataSetFromMatrix(countData = schaum_counts_bones,
                              colData = schaum_metadata_bones,
                              design = ~ age_categorical + sex)
Schaum_DESeq2_Results_Bones<-DESeq(dds_bones)
raw_results_bones <- results(Schaum_DESeq2_Results_Bones, contrast = c("age_categorical", "1_months", "24_months"))
formatted_results_bones <- formatDESeq2Results(raw_results_bones)
formatted_results_bones
sig_genes_bones_SY <- filter(formatted_results_bones, padj <= 0.01)
sig_genes_bones_SY 


########################################
# Approach 2: Extract information from original DESeq2 objects
########################################

# Load differential expression results from figshare

# list of tissues
ls(results_list_MACA_bulk)
# list of comparisons for Bone
ls(results_list_MACA_bulk$Bone)
# All results: comparison 3 vs 18 for Bone
results_list_MACA_bulk$Bone$`18_vs_3`$resall
# q<0.05 only: comparison 3 vs 18 for Bone
results_list_MACA_bulk$Bone$`18_vs_3`$ressig
# Extract DESeq2 object to new object
Bone_3_vs_18 <- results_list_MACA_bulk$Bone$`18_vs_3`$resall

# Load in raw counts table from figshare
schaum_counts_bones <- read_tsv("Schaum_Raw_Counts_Bones.tsv")
schaum_counts_bones <- column_to_rownames(schaum_counts_bones, var = "...1")
# Load in metadata table from figshar
# Additional metadata columns can/have been added 
schaum_metadata_bones <- read.delim("Schaum_Metadata_Bones.tsv", row.names = 1)
dds_bones <- DESeqDataSetFromMatrix(countData = schaum_counts_bones,
                                    colData = schaum_metadata_bones,
                                    design = ~ age_categorical + sex)
Schaum_DESeq2_Results_Bones<-DESeq(dds_bones)
#
raw_results_bones <- results(Schaum_DESeq2_Results_Bones, contrast = c("age_categorical", "12_months", "24_months"))
formatted_results_bones <- formatDESeq2Results(raw_results_bones)
formatted_results_bones
sig_genes_bones_SY <- filter(formatted_results_bones, padj <= 0.01)
sig_genes_bones_SY 
write.csv(sig_genes_bones_SY, "sig_genes_Bone_12v24_sig_SY.csv")

#Test_Code
Schaum_comp <- results_list_MACA_bulk$Bone$`12_vs_24`$resall
Schaum_comp_formatted <- formatDESeq2Results(Schaum_comp)
sig_genes_schaum <- filter(Schaum_comp_formatted, padj <= 0.01)
sig_genes_schaum
write.csv(sig_genes_schaum, "sig_genes_Bone_12v24_sig_Schaum.csv")


########################################
# Kidney
########################################
#Load in results from figshare
load("Schaum_Bulk_DESeq2_Results.bin")

# Load in raw counts table from figshare
schaum_counts_kidney <- read_tsv("Schaum_Raw_Counts_Kidney.tsv")
schaum_counts_kidney <- column_to_rownames(schaum_counts_kidney, var = "...1")
# Load in metadata table from figshar
# Additional metadata columns can/have been added 
schaum_metadata_kidney <- read.delim("Schaum_Metadata_Kidney.tsv", row.names = 1)
schaum_metadata_kidney$age <- as.character(schaum_metadata_kidney$age)
dds_kidney <- DESeqDataSetFromMatrix(countData = schaum_counts_kidney,
                                    colData = schaum_metadata_kidney,
                                    design = ~ age + sex)
Schaum_DESeq2_Results_Kidney<-DESeq(dds_kidney)
#

#
raw_results_kidney <- results(Schaum_DESeq2_Results_Kidney, contrast = c("age", "24", "27"))
formatted_results_kidney <- formatDESeq2Results(raw_results_kidney)
formatted_results_kidney
sig_genes_kidney_SY <- filter(formatted_results_kidney, padj <= 0.01)
sig_genes_kidney_SY 
write.csv(sig_genes_kidney_SY, "sig_genes_Kidney_24v27_sig_SY.csv")

Schaum_comp <- results_list_MACA_bulk$Kidney$`24_vs_27`$resall
Schaum_comp_formatted <- formatDESeq2Results(Schaum_comp)
sig_genes_schaum <- filter(Schaum_comp_formatted, padj <= 0.01)
sig_genes_schaum
write.csv(sig_genes_schaum, "sig_genes_Kidney_24v27_sig_Schaum.csv")

#

#Testing functions
library(clusterProfiler)

#Mouse clusterProfiler works out of the box, new functions added to the function area.
R1_vs_R2_clusters <- runClusterProfilerMouse(sig_genes_kidney_SY)
dim(R1_vs_R2_clusters)
dotplot(R1_vs_R2_clusters, showCategory=34, title="YOUR TITLE HERE", font.size=10, label_format = 50)



########################################
# Brain
########################################


# Load in raw counts table from figshare
schaum_counts_Brain <- read_tsv("Schaum_Raw_Counts_Brain.tsv")
schaum_counts_Brain <- column_to_rownames(schaum_counts_Brain, var = "...1")
# Load in metadata table from figshar
# Additional metadata columns can/have been added 
schaum_metadata_Brain <- read.delim("Schaum_Metadata_Brain.tsv", row.names = 1)
schaum_metadata_Brain$age <- as.character(schaum_metadata_Brain$age)
dds_Brain <- DESeqDataSetFromMatrix(countData = schaum_counts_Brain,
                                    colData = schaum_metadata_Brain,
                                    design = ~ age + sex)
Schaum_DESeq2_Results_Brain<-DESeq(dds_Brain)
#

#
raw_results_Brain <- results(Schaum_DESeq2_Results_Brain, contrast = c("age", "24", "27"))
formatted_results_Brain <- formatDESeq2Results(raw_results_Brain)
formatted_results_Brain
sig_genes_Brain_SY <- filter(formatted_results_Brain, padj <= 0.01)
sig_genes_Brain_SY 
write.csv(sig_genes_Brain_SY, "sig_genes_Brain_24v27_sig_SY.csv")

Schaum_comp <- results_list_MACA_bulk$Brain$`24_vs_27`$resall
Schaum_comp_formatted <- formatDESeq2Results(Schaum_comp)
sig_genes_schaum <- filter(Schaum_comp_formatted, padj <= 0.01)
sig_genes_schaum
write.csv(sig_genes_schaum, "sig_genes_Brain_24v27_sig_Schaum.csv")






########################################
# Heart
########################################


# Load in raw counts table from figshare
schaum_counts_Heart <- read_tsv("Schaum_Raw_Counts_Heart.tsv")
schaum_counts_Heart <- column_to_rownames(schaum_counts_Heart, var = "...1")
# Load in metadata table from figshare
# Additional metadata columns can/have been added 
schaum_metadata_Heart <- read.delim("Schaum_Metadata_Heart.tsv", row.names = 1)
schaum_metadata_Heart$age <- as.character(schaum_metadata_Heart$age)
dds_Heart <- DESeqDataSetFromMatrix(countData = schaum_counts_Heart,
                                    colData = schaum_metadata_Heart,
                                    design = ~ age + sex)
Schaum_DESeq2_Results_Heart<-DESeq(dds_Heart)
#

#
raw_results_Heart <- results(Schaum_DESeq2_Results_Heart, contrast = c("age", "24", "27"))
formatted_results_Heart <- formatDESeq2Results(raw_results_Heart)
formatted_results_Heart
sig_genes_Heart_SY <- filter(formatted_results_Heart, padj <= 0.01)
sig_genes_Heart_SY 
write.csv(sig_genes_Heart_SY, "sig_genes_Heart_24v27_sig_SY.csv")

Schaum_comp <- results_list_MACA_bulk$Heart$`24_vs_27`$resall
Schaum_comp_formatted <- formatDESeq2Results(Schaum_comp)
sig_genes_schaum <- filter(Schaum_comp_formatted, padj <= 0.01)
sig_genes_schaum
write.csv(sig_genes_schaum, "sig_genes_Heart_24v27_sig_Schaum.csv")
