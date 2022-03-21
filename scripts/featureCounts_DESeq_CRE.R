# Bioconductor version 3.14 (BiocManager 1.30.16), R 4.1.2 (2021-11-01)
# install.packages("vegan")
library( "DESeq2" )
library(ggplot2)
library(readxl)
library(vegan) # This is vegan 2.5-7

setwd("/Users/shuaishuai/OneDrive - Northwestern University/Probiotic_Competition")
# read in feature count output file
data<-read.table("./data/featureCounts-CRE.txt",
                 header = TRUE,skip=1,row.names = 1,
                 check.names = FALSE)
# read in metadata (sample information)
metadata <- read_excel("./data/metadata_r2_CRE.xlsx")

# exclude feature info
countdata<-data[,7:ncol(data)]
# chop off "_hisat-CRE.bam" in column names
colnames(countdata) <- sub("_hisat-CRE.bam","", colnames(countdata))
# add pseudocount 1
#countdata_pseudo <- countdata + 1

# chop off "_S12_R1_001" in sample names
metadata$sample <- substring(metadata$sample,1, nchar(metadata$sample)-11 )
# match column names
metadata<-metadata[match(colnames(countdata),metadata$sample),]
metadata$Condition <- factor(metadata$Condition)
metadata$BG <- factor(metadata$BG)
metadata$Time[which(metadata$Time=="A")]<-"3 Hours"
metadata$Time[which(metadata$Time=="B")]<-"24 Hours"
metadata$Time[which(metadata$Time=="C")]<-"72 Hours"
metadata$Time = factor(metadata$Time, 
                         levels=c('3 Hours','24 Hours','72 Hours'))

########## DESeq2
dds<-DESeqDataSetFromMatrix(countData = countdata_pseudo,
                            colData = metadata,
                            design=~Condition+Time+BG)
##### filter out zero count features (# of samples = 61)
keep <- rowSums(counts(dds)) >= 61
dds <- dds[keep,]

##### estimate size factor (get rid of error in estimateSizeFactorsForMatrix(counts(object), locfunc = locfunc,)
#dds <- estimateSizeFactors(dds,type="iterate") # "iterate" takes long time. Ended up using pseudocount

###################### Differential expression
##### second batch data, including samples with only BG and no CRE
dds <- DESeq(dds,
              test = "Wald",
              fitType = "parametric", 
              sfType = "ratio")
### Results
# res <- results(dds)
res <- results(dds, contrast=c("BG","VEG","N"))

###################### Differential expression subsamples###################################
######### dss1: Only keep samples with CRE231, use countdata with pseudocount
metadata1<-metadata[which(metadata$CRE=="Y"),]
countdata1<-countdata[,match(metadata1$sample,colnames(countdata))]

dds1<-DESeqDataSetFromMatrix(countData = countdata1,
                             colData = metadata1,
                             design=~Condition+Time+BG)
dds1 <- DESeq(dds1,
             test = "Wald",
             fitType = "parametric", 
             sfType = "ratio")

##### Results
#### define reference level
res <- results(dds1, contrast=c("Time","72 Hours","24 Hours"))
res
summary(res)
sum(res$padj <= 0.05 & abs(res$log2FoldChange) >= 1, na.rm=TRUE)
resFilt <- res[which(res$padj <= 0.05 & abs(res$log2FoldChange) >= 1), ]
write.csv(resFilt, file="data/DESeq2_72h_vs_24h_filtered_new.csv")

##### Results
#### define reference level
res <- results(dds1, contrast=c("BG","AF","VEG"))
res
summary(res)
sum(res$padj <= 0.05 & abs(res$log2FoldChange) >= 1, na.rm=TRUE)
resFilt <- res[which(res$padj <= 0.05 & abs(res$log2FoldChange) >= 1), ]
write.csv(resFilt, file="data/DESeq2_CAF_vs_CVEG_filtered_new.csv")
####################################################################################
########## dss2: Only keep samples with APPC+CRE231, use countdata
metadata2<-metadata[which(metadata$BG=="BG"),]
countdata2<-countdata[,match(metadata2$sample,colnames(countdata))]

dds2<-DESeqDataSetFromMatrix(countData = countdata2,
                             colData = metadata2,
                             design=~Condition+Time)
dds2 <- DESeq(dds2,
              test = "Wald",
              fitType = "parametric", 
              sfType = "ratio")

##### Results
### Wet as baseline
res <- results(dds2, contrast=c("Condition","Ambient","Wet"))
res
summary(res)
sum(res$padj <= 0.05 & abs(res$log2FoldChange) >= 1, na.rm=TRUE)
resFilt <- res[which(res$padj <= 0.05 & abs(res$log2FoldChange) >= 1), ]
write.csv(resFilt, file="data/DESeq2_APPCa_vs_APPCw_filtered_new.csv")
##############################################################################################
########### dds3: Only keep samples with no cleaning CRE231, use countdata
metadata3<-metadata[which(metadata$BG=="N"),]
countdata3<-countdata[,match(metadata3$sample,colnames(countdata))]

dds3<-DESeqDataSetFromMatrix(countData = countdata3,
                             colData = metadata3,
                             design=~Condition+Time)
dds3 <- DESeq(dds3,
              test = "Wald",
              fitType = "parametric", 
              sfType = "ratio")

######### Results
##### Wet as baseline
res <- results(dds3, contrast=c("Condition","Ambient","Wet"))
res
summary(res)
sum(res$padj <= 0.05 & abs(res$log2FoldChange) >= 1, na.rm=TRUE)
resFilt <- res[which(res$padj <= 0.05 & abs(res$log2FoldChange) >= 1), ]
write.csv(resFilt, file="data/DESeq2_CREa_vs_CREw_filtered_new.csv")

###############################################################################################
##################### Variance stabilization and PCA ##########################################
###############################################################################################
vsd <- vst(dds1, blind=FALSE)
pcaData <- plotPCA(vsd, intgroup=c("Condition","Time","BG"), 
                   returnData=TRUE)
## percentage variance *100
percentVar <- round(100 * attr(pcaData, "percentVar"))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
######## PCA color coded with different Cleaning Scenario
P1<-ggplot(pcaData, 
           aes(PC1, PC2, shape=Condition,
               color=factor(BG, 
                            levels = c("N","AF","BG","VEG"),
                            labels = c("CRE231", 
                                       "CRE231 + APC", 
                                       "CRE231 + APPC",
                                       "CRE231 + VB")))) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  theme_bw()+
  #stat_ellipse(type = "norm", linetype = 2,size=1)+
  theme(axis.text.x = element_text(size=10,face="bold"),
        axis.title.x = element_text(face="bold"),
        axis.title.y = element_text(face="bold"),
        legend.text = element_text(size=10),
        legend.title = element_text(face="bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))+
  scale_color_brewer(palette="RdYlBu",name="Cleaning Scenario")+
  scale_shape_discrete(name="Temperature/Humidity Condition",
                       labels=c("25°C and 20% RH",
                                "37°C and 90% RH"))
P1

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
######## PCA color coded with different Time
P2<-ggplot(pcaData, 
           aes(PC1, PC2, shape=Condition,
               color=factor(Time, 
                            levels = c("3 Hours","24 Hours","72 Hours"),
                            labels = c("3 Hours",
                                       "24 Hours",
                                       "72 Hours")))) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  theme_bw()+
  #stat_ellipse(type = "norm", linetype = 2,size=1)+
  theme(axis.text.x = element_text(size=10,face="bold"),
        axis.title.x = element_text(face="bold"),
        axis.title.y = element_text(face="bold"),
        legend.text = element_text(size=10),
        legend.title = element_text(face="bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))+
  scale_color_brewer(palette="Dark2",name="Cleaning Time")+
  scale_shape_discrete(name="Temperature/Humidity Condition",
                       labels=c("25°C and 20% RH",
                                "37°C and 90% RH"))
P2


#############################################################################################
############################ PREMANOVA ######################################################
##### Data preparation (dds1: CRE data)
# rlog transformation 
rld <- rlog(dds1,blind = FALSE)
plotPCA(rld, intgroup = c("BG", "Condition"))

# variance stabilization
vsd <- vst(dds1, blind=FALSE)

##### Perform PERMANOVA
permanova <- adonis(t(counts(dds1,normalized=TRUE)) ~ Condition+BG+Time,
                    data = metadata1, permutations=9999, method = "euclidean")

permanova <- adonis(assay(vsd)%>%t() ~ Condition+BG+Time,
                    data = metadata1, permutations=9999, method = "euclidean")

permanova <- adonis(assay(rld)%>%t() ~ Condition+BG+Time,
                    data = metadata1, permutations=9999, method = "euclidean")

permanova$aov.tab

#############################################################################################


