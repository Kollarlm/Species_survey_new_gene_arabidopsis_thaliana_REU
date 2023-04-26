# Set working directory and load necessary packages
library(DESeq2)
library(dplyr)
library(stringr)
library(ggplot2)
############################################# RUNNING DESEQ2 #############################################
##Preparing files for DESeq 
# ff <- list.files( path = "./data/DEG_L1/", pattern = "*L1_star_2pass", full.names = TRUE )
# counts.files <- lapply( ff, read.table, skip = 4 )
# counts <- as.data.frame( sapply( counts.files, function(x) x[ , 4 ] ) )
# ff <- gsub( "[.]ReadsPerGene[.]out[.]tab", "", ff )
# ff <- gsub( "[.]/counts/", "", ff )
# colnames(counts) <- ff
# row.names(counts) <- counts.files[[1]]$V1
# colnames(counts) <- str_replace_all(colnames(counts),"./all-counts-V2/","")
# colnames(counts) <- str_replace_all(colnames(counts),"ReadsPerGene.out.tab","")
colData <- read.table(file="./data/DEG_L1/expression-design.tsv", header = T)
colData <- arrange(colData, desc(Genotype))
data <- read.table("./data/DEG_L1/L1_star_2pass.tsv", header = T, row.names = 1)
#counts <- counts[,c(7,8,9,10,11,12,1,2,3,4,5,6)]

##Making DESeq Data Set
dds = DESeqDataSetFromMatrix(countData = data, colData = colData, design = ~Genotype + Condition + Genotype:Condition)

#Make sure that dds is properly formatted
as.data.frame(colData(dds))

##Run DESeq
dds <- DESeq(dds)
res <- results(dds)

############################### EXTRACTING AND FILTERING SPECIFIC SAMPLES ################################
## Double check result names
resultsNames(dds)

## Filter differentially expressed genes between L1_coastal and L1_inland
L1IL1C <- results(dds, contrast=c("Condition", "Coastal", "Inland"))
L1 <- as.data.frame(L1IL1C)
L11<-L1[(L1$log2FoldChange>1|L1$log2FoldChange< -1)& !is.na(L1$log2FoldChange),]
L12 <- L11[(L11$padj<0.05)&!is.na(L11$padj),]
L12$gene_model <- row.names(L12)

## Filter differentially expressed genes between S1_coastal and S1_inland
S1IS1C <- results(dds, list(c("Condition_Inland_vs_Coastal","GenotypeS1.ConditionInland")))
S1 <- as.data.frame(S1IS1C)
S11<-S1[(S1$log2FoldChange>1|S1$log2FoldChange< -1)& !is.na(S1$log2FoldChange),]
S12 <- S11[(S11$padj<0.05)&!is.na(S11$padj),]
S12$gene_model <- row.names(S12)

# ## What genes are different between field sites samples AND between genotypes (i.e. are the genes that are different between sites the same across genotypes)?
field_siteres <- results(dds, name="GenotypeS1.ConditionInland")
field_site <- as.data.frame(field_siteres)
field_site1<-field_site[(field_site$log2FoldChange>1|field_site$log2FoldChange< -1)& !is.na(field_site$log2FoldChange),]
field_site2 <- field_site1[(field_site1$padj<0.05)&!is.na(field_site1$padj),]
field_site2$gene_model <- row.names(field_site2)

## What genes are different between genotypes?
genores <- results(dds, name="Genotype_S1_vs_L1")
geno <- as.data.frame(genores)
geno1<-geno[(geno$log2FoldChange>1|geno$log2FoldChange< -1)& !is.na(geno$log2FoldChange),]
geno2 <- geno1[(geno1$padj<0.05)&!is.na(geno1$padj),]
geno2$gene_model <- row.names(geno2)

## Pulling out intersting genes

# Genes we found interesting because they are ...
  # 1. orthologs with genes found in GOuld et al 2017
  # 2. near breakpoint regions
  # 3. intersting outliers

int_genes <- read.csv("./data/Intersting_genes_for_expression.csv", header = T)

#### Shared between our interesting outliers and the DEGs between field sites between genotypes ####
DE_outliers_int <- subset(int_genes,   int_genes$gene_model %in% field_site2$gene_model)

DE_outliers_int2 <- DE_outliers_int %>% 
  inner_join(field_site2, by ="gene_model", multiple = "all") # Adding L1 outlier data to dataset of syntenic only genes
colnames(DE_outliers_int2) <- c("gene_model", "chromosome", "location", "function", "reason", "high_pi", "low_pi", "gstat", "fst", "baseMean",                 
                                "log2FoldChange","lfcSE","stat" ,"pvalue" , "padj" )

## Shared between our interesting outliers and DE in specific genotype

L1_outliers_int <- subset(int_genes,  int_genes$gene_model %in% L12$gene_model)

L1_outliers_int2 <- L1_outliers_int %>% 
  inner_join(L12, by ="gene_model", multiple = "all") # Adding L1 outlier data to dataset of syntenic only genes
colnames(L1_outliers_int2) <- c("gene_model", "chromosome", "location", "function", "reason", "high_pi", "low_pi", "gstat", "fst", "baseMean",                 
                                "log2FoldChange","lfcSE","stat" ,"pvalue" , "padj" )

### Checking for DE in all outlier genes ###

# DE between genotypes and between field sites
L1_outliers <- read.csv("/Users/lesliekollar/Desktop/mimulus-gould2018-reanalysis/data/poolseq_output/L1_aligned/genes/L1_all_data_unique_collapsed.csv")
#L1_outliers <-read.csv("/Users/lesliekollar/Desktop/mimulus-gould2018-reanalysis/data/poolseq_output/L1_aligned/genes/L1_all_data_unique_collapsed.csv")

DE_outliers_L1 <- subset(L1_outliers,   L1_outliers$gene_model %in% field_site2$gene_model)

DE_outliers_L12 <- DE_outliers_L1 %>% 
  inner_join(field_site2, by ="gene_model", multiple = "all") # Adding L1 outlier data to dataset of syntenic only genes
colnames(DE_outliers_L12) <- c("gene_model", "chromosome", "location", "function", "reason", "high_pi", "low_pi", "gstat", "fst", "baseMean",                 
                               "log2FoldChange","lfcSE","stat" ,"pvalue" , "padj" )

# Shared between all outliers and DE in specific genotype

L1_outliers_all <- subset(L1_outliers,  L1_outliers$gene_model %in% L12$gene_model)

L1_outliers_all2 <- L1_outliers_all %>% 
  inner_join(L12, by ="gene_model", multiple = "all") # Adding L1 outlier data to dataset of syntenic only genes
colnames(L1_outliers_all2) <- c("gene_model", "chromosome", "location", "function", "reason", "high_pi", "low_pi", "gstat", "fst", "baseMean",                 
                                "log2FoldChange","lfcSE","stat" ,"pvalue" , "padj" )

# Making plots

gene <- plotCounts(dds,"MgL1_08g16120",intgroup=c("Genotype","Condition"),returnData = TRUE)
gene$genotype_condition <- paste(colData$Genotype,colData$Condition,sep="_")
ggplot(gene,aes(x=genotype_condition,y=count,color=genotype_condition)) + geom_boxplot(fill=NA) + geom_jitter() +  theme_bw() + ggtitle("MgL1_08g16120")

gene <- plotCounts(dds,"MgL1_14g10910",intgroup=c("Genotype","Condition"),returnData = TRUE)
gene$genotype_condition <- paste(colData$Genotype,colData$Condition,sep="_")
ggplot(gene,aes(x=genotype_condition,y=count,color=genotype_condition)) + geom_boxplot(fill=NA) + geom_jitter() +  theme_bw() + ggtitle("MgL1_14g10910")

gene <- plotCounts(dds,"MgL1_05g11230",intgroup=c("Genotype","Condition"),returnData = TRUE)
gene$genotype_condition <- paste(colData$Genotype,colData$Condition,sep="_")
ggplot(gene,aes(x=genotype_condition,y=count,color=genotype_condition)) + geom_boxplot(fill=NA) + geom_jitter() +  theme_bw() + ggtitle("MgL1_05g11230")

gene <- plotCounts(dds,"MgL1_05g15660",intgroup=c("Genotype","Condition"),returnData = TRUE)
gene$genotype_condition <- paste(colData$Genotype,colData$Condition,sep="_")
ggplot(gene,aes(x=genotype_condition,y=count,color=genotype_condition)) + geom_boxplot(fill=NA) + geom_jitter() +  theme_bw() + ggtitle("MgL1_05g15660")

gene <- plotCounts(dds,"MgL1_08g06810",intgroup=c("Genotype","Condition"),returnData = TRUE)
gene$genotype_condition <- paste(colData$Genotype,colData$Condition,sep="_")
ggplot(gene,aes(x=genotype_condition,y=count,color=genotype_condition)) + geom_boxplot(fill=NA) + geom_jitter() +  theme_bw() + ggtitle("MgL1_08g06810")

gene <- plotCounts(dds,"MgL1_08g23640",intgroup=c("Genotype","Condition"),returnData = TRUE)
gene$genotype_condition <- paste(colData$Genotype,colData$Condition,sep="_")
ggplot(gene,aes(x=genotype_condition,y=count,color=genotype_condition)) + geom_boxplot(fill=NA) + geom_jitter() +  theme_bw() + ggtitle("MgL1_08g23640")

gene <- plotCounts(dds,"MgL1_08g01750",intgroup=c("Genotype","Condition"),returnData = TRUE)
gene$genotype_condition <- paste(colData$Genotype,colData$Condition,sep="_")
ggplot(gene,aes(x=genotype_condition,y=count,color=genotype_condition)) + geom_boxplot(fill=NA) + geom_jitter() +  theme_bw() + ggtitle("MgL1_08g01750")


# Myb anthocyanin

# amMyb2
# This one has been duplicated.
# Myb1, Myb2, and 3 were tandem duplicates (happened in mg specifically) but evolved seperately. (L1 and S1 both have)
# Common ancester did not have them.
# Expressed higher in L1 than S1 in all locations
# Based on Billies its only happening in the leaf.

gene <- plotCounts(dds,"MgL1_08g07000",intgroup=c("Genotype","Condition"),returnData = TRUE)
gene$genotype_condition <- paste(colData$Genotype,colData$Condition,sep="_")
ggplot(gene,aes(x=genotype_condition,y=count,color=genotype_condition)) + geom_boxplot(fill=NA) + geom_jitter() +  theme_bw() + ggtitle("MgL1_08g07000")

