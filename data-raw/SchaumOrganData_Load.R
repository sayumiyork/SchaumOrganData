############################################################
#Load libraries
############################################################
library("tidyverse")
library("DESeq2")

############################################################
# Processed files
############################################################

# kidney
kidney<-readRDS("data-raw/Schaum_Kidney.RDS") #DESeq2 object
kidney_metadata<-read.delim("data-raw/Schaum_Metadata_Kidney.tsv", row.names = 1)

# brain
brain<-readRDS("data-raw/Schaum_Brain.RDS") #DESeq2 object
brain_metadata<-read.delim("data-raw/Schaum_Metadata_Brain.tsv", row.names = 1)

# heart
heart<-readRDS("data-raw/Schaum_Heart.RDS") #DESeq2 object
heart_metadata<-read.delim("data-raw/Schaum_Metadata_Heart.tsv", row.names = 1)



############################################################
#usethis processed files

############################################################


#Kidney
usethis::use_data(kidney, overwrite = TRUE)
usethis::use_data(kidney_metadata, overwrite = TRUE)
#Brain
usethis::use_data(brain, overwrite = TRUE)
usethis::use_data(brain_metadata, overwrite = TRUE)
#Heart
usethis::use_data(heart, overwrite = TRUE)
usethis::use_data(heart_metadata, overwrite = TRUE)
