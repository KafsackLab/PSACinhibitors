# Load libraries
library("MSnbase")
library("MSstats")
library("limma")
library("EnhancedVolcano")
library("DEP")
library("tidyverse")

#Set directory
setwd("~/proteomics/")

# Set per-arm missing value percent cutoff, minimum number of peptides, generic identifiers
filter_percent_cutoff <- 0.5
peptide_count_cutoff <- 2
data_identifier <- "LFQ"

#Volcano Plot Labeling Cutoffs
PValue_cutoff = 1e-2
logFC_cutoff = 1

#logic setting up data structure
  gene_column <- "Gene"
  protein_column <- "Protein.ID"
  metadata_columns <- 11
  filename <- "Proteomics.C1&C2_microMap_peptideData.tsv"

# Load data, make unique name column, assign NA as required (MSfragger reports 0)
data <- read.delim(filename)
data_unique <- make_unique(data, gene_column, protein_column)
rownames(data_unique) <- data_unique[,"name"]
names_column <- grep(data_identifier, colnames(data_unique))
data_unique[,names_column][data_unique[,names_column] == 0] <- NA

#sample key
sampleKey<-read.delim("Proteomics.sample_key.tab")
colnames(sampleKey)<-c("sampleID","date","bRep","tRep","treat","ug.IP","uL.tryp")
sampleKey$bRep<-LETTERS[sampleKey$bRep]
sampleKey$tRep<-letters[sampleKey$tRep]
sampleKey$catalyst<-substr(sampleKey$treat,1,13)
sampleKey$competitor<-rep(c(T,F),each=3)
sampleKey$dataCol<-grep("MaxLFQ",colnames(data_unique),value=F)
sampleKey$dataColName<-grep("MaxLFQ",colnames(data_unique),value=T)

sampleKey

################################################################################################
### Uncomment the desired set of 3 lines below, then run the rest of the code
################################################################################################

###################################
# Compound 1 analog TDI-015668-Ir #
###################################
title= "TDI-015668-Ir vs TDI-015668-Ir + excess TDI-015668 \n all three Replicates combined"
name="TDI-015668-Ir_rep1-3"
neg_column_names <- sampleKey$dataCol[sampleKey$catalyst=="TDI-015668-Ir" & sampleKey$competitor==F]
pos_column_names <- sampleKey$dataCol[sampleKey$catalyst=="TDI-015668-Ir" & sampleKey$competitor==T]

# title= "TDI-015668-Ir vs TDI-015668-Ir + excess TDI-015668 \n Biological Replicate 1"
# name="TDI-015668-Ir_rep1"
# neg_column_names <- sampleKey$dataCol[sampleKey$catalyst=="TDI-015668-Ir" & sampleKey$competitor==F & sampleKey$bRep=="A"]
# pos_column_names <- sampleKey$dataCol[sampleKey$catalyst=="TDI-015668-Ir" & sampleKey$competitor==T & sampleKey$bRep=="A"]

# title= "TDI-015668-Ir vs TDI-015668-Ir + excess TDI-015668 \n Biological Replicate 2"
# name="TDI-015668-Ir_rep2"
# neg_column_names <- sampleKey$dataCol[sampleKey$catalyst=="TDI-015668-Ir" & sampleKey$competitor==F & sampleKey$bRep=="B"]
# pos_column_names <- sampleKey$dataCol[sampleKey$catalyst=="TDI-015668-Ir" & sampleKey$competitor==T & sampleKey$bRep=="B"]

# title= "TDI-015668-Ir vs TDI-015668-Ir + excess TDI-015668 \n Biological Replicate 3"
# name="TDI-015668-Ir_rep3"
# neg_column_names <- sampleKey$dataCol[sampleKey$catalyst=="TDI-015668-Ir" & sampleKey$competitor==F & sampleKey$bRep=="C"]
# pos_column_names <- sampleKey$dataCol[sampleKey$catalyst=="TDI-015668-Ir" & sampleKey$competitor==T & sampleKey$bRep=="C"]

###################################
# Compound 2 analog TDI-015836-Ir #
###################################

# title= "TDI-015836-Ir vs TDI-015836-Ir + excess TDI-015836 \n all three Replicates combined"
# name="TDI-015668-Ir_rep1-3"
# neg_column_names <- sampleKey$dataCol[sampleKey$catalyst=="TDI-015836-Ir" & sampleKey$competitor==F]
# pos_column_names <- sampleKey$dataCol[sampleKey$catalyst=="TDI-015836-Ir" & sampleKey$competitor==T]

# title= "TDI-015836-Ir vs TDI-015836-Ir + excess TDI-015836 \n Biological Replicate 1"
# name="TDI-015836-Ir_rep1"
# neg_column_names <- sampleKey$dataCol[sampleKey$catalyst=="TDI-015836-Ir" & sampleKey$competitor==F & sampleKey$bRep=="A"]
# pos_column_names <- sampleKey$dataCol[sampleKey$catalyst=="TDI-015836-Ir" & sampleKey$competitor==T & sampleKey$bRep=="A"]

# title= "TDI-015836-Ir vs TDI-015836-Ir + excess TDI-015836 \n Biological Replicate 2"
# name="TDI-015836-Ir_rep2"
# neg_column_names <- sampleKey$dataCol[sampleKey$catalyst=="TDI-015836-Ir" & sampleKey$competitor==F & sampleKey$bRep=="B"]
# pos_column_names <- sampleKey$dataCol[sampleKey$catalyst=="TDI-015836-Ir" & sampleKey$competitor==T & sampleKey$bRep=="B"]

# title= "TDI-015836-Ir vs TDI-015836-Ir + excess TDI-015836 \n Biological Replicate 3"
# name="TDI-015836-Ir_rep3"
# neg_column_names <- sampleKey$dataCol[sampleKey$catalyst=="TDI-015836-Ir" & sampleKey$competitor==F & sampleKey$bRep=="C"]
# pos_column_names <- sampleKey$dataCol[sampleKey$catalyst=="TDI-015836-Ir" & sampleKey$competitor==T & sampleKey$bRep=="C"]

##################################
# Run the rest of the code below #
##################################

# Rename arm-containing columns
colnames(data_unique)[pos_column_names]  <- "high_a"
colnames(data_unique)[neg_column_names]  <- "low_a"

# Define imputation cases for each protein, filter proteins based on missing value cutoff. Filtration of missing values is on a one-arm basis for separate arms,
# so that data is retained if one arm passes the cutoff and "all-negative" arms are retained. For uniform analysis missing values are filtered on an all-arm basis.
data_unique$NA_percent_neg <- rowSums(apply(is.na(data_unique[,neg_column_names]),2,as.numeric))/ncol(data_unique[,neg_column_names])
data_unique$NA_percent_pos <- rowSums(apply(is.na(data_unique[,pos_column_names]),2,as.numeric))/ncol(data_unique[,pos_column_names])

data_unique$fully_imputed <- as.numeric(data_unique[,"NA_percent_pos"] == 1 | data_unique[,"NA_percent_neg"] == 1)
data_unique$partially_imputed <- as.numeric(!data_unique[,"NA_percent_pos"] %in% c(0,1) | !data_unique[,"NA_percent_neg"] %in% c(0,1))

data_unique$not_imputed <- as.numeric(data_unique[,"NA_percent_pos"] == 0 & data_unique[,"NA_percent_neg"] == 0)
data_unique_filtered <- subset(data_unique, NA_percent_neg < filter_percent_cutoff | NA_percent_pos < filter_percent_cutoff)
data_unique_filtered <- subset(data_unique_filtered, Combined.Total.Peptides > peptide_count_cutoff)

#For imputation and normalization of separate experimental arms:
# Convert into treatment-group unique msnset objects for downstream processing
  msnset_neg <- readMSnSet2(data_unique_filtered, neg_column_names, "name")
  msnset_pos <- readMSnSet2(data_unique_filtered, pos_column_names, "name")
  
  a<- readMSnSet2(data_unique_filtered, neg_column_names, "name")
  b<- readMSnSet2(data_unique_filtered, neg_column_names, "name")
  
  # Normalize in each treatment group separately
  msnset_neg_norm <- normalise(msnset_neg, "vsn")
  msnset_pos_norm <- normalise(msnset_pos, "vsn")
  
  # Combine and generate matrix of unimputed data
  combined_preprocessing_postnorm <- MSnbase::combine(msnset_neg_norm, msnset_pos_norm)
  data_combined_preprocessing_postnorm_names <- combined_preprocessing_postnorm@featureData@data[1:metadata_columns]
  data_combined_preprocessing_postnorm_exprs <- combined_preprocessing_postnorm@assayData[["exprs"]]
  data_preprocessed_postnorm <- cbind(data_combined_preprocessing_postnorm_names, data_combined_preprocessing_postnorm_exprs)
  data_preprocessed_postnorm_unique <- make_unique(data_preprocessed_postnorm, gene_column, protein_column)
  write.csv(data_preprocessed_postnorm_unique, file = paste0(name,"_report_no_imputation.csv"))
  
  # Impute in each treatment group separately, combine into one msnset object, generate matrix of imputed data
  msnset_neg_norm_imputed <- MSnbase::impute(msnset_neg_norm, "MinProb")
  msnset_pos_norm_imputed <- MSnbase::impute(msnset_pos_norm, "MinProb")
  combined_preprocessing <- MSnbase::combine(msnset_neg_norm_imputed,
                                             msnset_pos_norm_imputed)

# Combine and generate matrix of imputed data
data_combined_preprocessing_names <- combined_preprocessing@featureData@data[1:metadata_columns]
data_combined_preprocessing_exprs <- combined_preprocessing@assayData[["exprs"]]
data_preprocessed <- cbind(data_combined_preprocessing_names, data_combined_preprocessing_exprs)
data_preprocessed_unique <- make_unique(data_preprocessed, gene_column, protein_column)
write.csv(data_preprocessed_unique, file =  paste0(name,"_report_imputed.csv"))

# Find data in imputed matrix, convert into suitable format for limma, create limma design matrix
data_limma <- data.matrix(data_preprocessed_unique[c(grep("low_a", colnames(data_preprocessed_unique)),
                                                     grep("high_a", colnames(data_preprocessed_unique)))])
colnames <- 1:sum(length(pos_column_names), length(neg_column_names))
colnames(data_limma) <- 1:sum(length(pos_column_names), length(neg_column_names))
design <- cbind(grp1=1,grp2=c(rep(0, length(neg_column_names)),rep(1, length(pos_column_names))))

# Run limma with Benjamimi Hochberg moderated T for P-Adjust
limma_output <- eBayes(lmFit(data_limma, design))

# Extract all DEP data from limma output, combine with descriptive columns and generate export matrix
# For DIA data incorporate a column for the total peptides used for each protein.
limma_output_toptable <- topTable(limma_output, coef=2, n = 10000)
limma_output_toptable$neg_log_p <- -log10(limma_output_toptable[,"P.Value"])
limma_output_toptable <- cbind(row.names(limma_output_toptable), limma_output_toptable)
colnames(limma_output_toptable)[1] <- "name"
full_report <- merge(limma_output_toptable,
                     data_unique_filtered[grep("imputed", colnames(data_unique_filtered))],
                     by = 0)
full_report <- full_report[,-1]
rownames(full_report) <- full_report[,"name"]
if(data_type == "DIA") {
  full_report <- merge(full_report,
                       data_unique_filtered[grep("Combined.Total.Peptides", colnames(data_unique_filtered))],
                       by = 0)
  full_report <- full_report[,-1]
  colnames(full_report)[1] <- "name"
  rownames(full_report) <- full_report[,"name"]
}
full_report <- merge(full_report,
                     data_preprocessed_unique[c(grep("low_a", colnames(data_preprocessed_unique)),grep("high_a", colnames(data_preprocessed_unique)))],
                     by = 0)
full_report <- full_report[,-1]
colnames(full_report)[1] <- "name"
rownames(full_report) <- full_report[,"name"]
full_report <- merge(full_report,
                     data_unique_filtered[,1:metadata_columns],
                     by = 0)
full_report <- full_report[,-1]
rownames(full_report) <- full_report[,"name"]
colnames(full_report)[1] <- "name"
full_report <- full_report[order(full_report$P.Value),]


rownames(full_report) <- full_report[,"Protein"]

write.csv(full_report, file = paste0(name,"_report_imputed.csv"))

full_report$logFC<- full_report$logFC
full_report<-full_report[grep("PF3D7",full_report$name),]
# Define volcano plot custom colors, labels, and cutoffs
keyvals <- ifelse(
  full_report$fully_imputed == 1 & full_report$logFC > logFC_cutoff & full_report$P.Value < PValue_cutoff, "firebrick1",
  ifelse(full_report$partially_imputed == 1 & full_report$logFC > logFC_cutoff & full_report$P.Value < PValue_cutoff, "firebrick2",
         ifelse(full_report$not_imputed == 1 & full_report$logFC > logFC_cutoff & full_report$P.Value < PValue_cutoff,  "firebrick3",
                "black")))
keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == 'firebrick1'] <- 'fully imputed'
names(keyvals)[keyvals == 'firebrick2'] <- 'partially imputed'
names(keyvals)[keyvals == 'firebrick3'] <- 'not imputed'  

# Plot volcano plot
plot_max_logfc <- max(full_report$logFC) + 0.5
plot_max_logp <- max(full_report$neg_log_p)
plot_min_logfc <- min(full_report$logFC) -0.5

EnhancedVolcano(full_report,
                lab = row.names(full_report),
                x = 'logFC',
                y = 'P.Value',
                #xlim = c(-2,2),
                ylim = c(0,plot_max_logp),
                col=c('black', 'black', 'black', 'red3'),
                FCcutoff = logFC_cutoff,
                pCutoff = PValue_cutoff,
                labSize = 3,
                colAlpha = 1,
                drawConnectors = TRUE,
                widthConnectors = 0.5,
                arrowheads = FALSE,
                gridlines.major = FALSE,
                gridlines.minor = FALSE,
                subtitle = NULL,
                legendPosition = "none",
                title = title,
                cutoffLineType = "blank",
                colCustom = keyvals,
                selectLab = rownames(full_report)[which(names(keyvals) %in% c('fully imputed', 'partially imputed', 'not imputed'))]
                )
