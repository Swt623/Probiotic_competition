library( "DESeq2" )
library(ggplot2)
library(readxl)
library(vegan) # This is vegan 2.5-7

setwd("/Users/shuaishuai/OneDrive - Northwestern University/Probiotic_Competition")
# read in feature count output file
data_BG<-read.table("./data/featureCounts_BG_r1_r2.txt",
                    header = TRUE,row.names = 1,
                    check.names = FALSE)
# read in metadata (sample information)
metadata_BG <- read_excel("./data/metadata_for_BG.xlsx")

##### filter out zero count features
keep <- rowSums(data_BG) > 0
data_BG <- data_BG[keep,]

# add pseudocount 1
######## use psedocounts b/c encounter this error even after using 
## > keep <- rowSums(counts(dds)) >=1
## > dds <- dds[keep,]
######## to get rid of all zero entries
## estimating size factors
## Error in estimateSizeFactorsForMatrix(counts(object), locfunc = locfunc,  : 
##  every gene contains at least one zero, cannot compute log geometric means
#countdata_pseudo <- data_BG + 1

# match column names
#metadata_BG<-metadata_BG[match(colnames(countdata_pseudo),metadata_BG$sample),]
metadata_BG<-metadata_BG[match(colnames(countdata),metadata_BG$sample),]
metadata_BG$Condition <- factor(metadata_BG$Condition)
metadata_BG$BG <- factor(metadata_BG$BG)
metadata_BG$Batch <- factor(metadata_BG$Batch)
metadata_BG$Pathogen <- factor(metadata_BG$Pathogen)
metadata_BG$Sp <- factor(metadata_BG$Sp)
metadata_BG$Time[which(metadata_BG$Time=="A")]<-"3 Hours"
metadata_BG$Time[which(metadata_BG$Time=="B")]<-"24 Hours"
metadata_BG$Time[which(metadata_BG$Time=="C")]<-"72 Hours"
metadata_BG$Time = factor(metadata_BG$Time, 
                          levels=c('3 Hours','24 Hours','72 Hours'))

###################### Only keep data for vegetative BG ######################
###################### Condition: Ambient for all samples
metadata3<-metadata_BG[which(metadata_BG$BG=="VEG"),]
#countdata3<-countdata_pseudo[,match(metadata3$sample,colnames(countdata_pseudo))]
countdata3<-countdata_pseudo[,match(metadata3$sample,colnames(countdata))]

###################### DEG for VB+ABBL, VB+CRE vs. VB ######################
dds3<-DESeqDataSetFromMatrix(countData = countdata3,
                             colData = metadata3,
                             design=~Time+Pathogen+Batch)

"dds3<-DESeqDataSetFromMatrix(countData = countdata3,
                             colData = metadata3,
                             design=~Time+Sp)"

###### design = ~ Time + Pathogen + Batch is ok, but design = ~ Time + Sp + Batch is not full-rank
dds3 <- DESeq(dds3,
              test = "Wald",
              fitType = "parametric", 
              sfType = "ratio")

######### Results
res <- results(dds3, contrast=c("Sp","CRE231","N"))
res
summary(res)
sum(res$padj <= 0.05 & abs(res$log2FoldChange) >= 1, na.rm=TRUE)
resFilt <- res[which(res$padj <= 0.05 & abs(res$log2FoldChange) >= 1), ]

###################### Variant stablization ######################
vsd <- vst(dds3, blind=FALSE)

###################### PCA ######################
pcaData <- plotPCA(vsd, intgroup=c("Time","Sp"), 
                   returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar")) ## percentage variance *100

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
######## PCA color coded with different Cleaning Scenario
P1<-ggplot(pcaData, 
           aes(PC1, PC2, shape=Time,
               color=factor(Sp, 
                            levels = c("ABBL18","CRE231", "N")))) +
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
  scale_color_brewer(palette="Dark2",name="Pathogen Species")+
  scale_shape_discrete(name="Time")
P1

########################################################################################
########################## KEEG pathway ###############################################
#### Read in data files
library(readr)
library(dplyr)
library(tidyr)
metabolism <- readxl::read_xlsx("./data/round1/KEGG.xlsx",sheet = "Metabolism")
genetic<-readxl::read_xlsx("./data//round1/KEGG.xlsx",sheet = "Genetic")
env<-readxl::read_xlsx("./data//round1/KEGG.xlsx",sheet = "Environmental")
cellular<-readxl::read_xlsx("./data/round1/KEGG.xlsx",sheet = "Cellular")
unclassified<-readxl::read_xlsx("./data/round1/KEGG.xlsx",sheet = "Sheet1")
brite<-readxl::read_xlsx("./data/round1/KEGG.xlsx",sheet = "Brite")
others<-readxl::read_xlsx("./data/round1/KEGG.xlsx",sheet = "Sheet1")
kegg<-rbind(metabolism,genetic)%>%
  rbind(env)%>%
  rbind(cellular)%>%
  rbind(unclassified)%>%
  rbind(brite)%>%
  rbind(others)

#### prokka geneid vs ko number
#taxa_B_o <-readxl::read_xlsx("./data/round1/blastkoala.xlsx", sheet = "BG")%>% drop_na()  # id start with BJKFFFGJ
taxa_B <-readxl::read_xlsx("./data/round1/blastkoala.xlsx", sheet = "BG_WO_rRNA")%>% drop_na() # ID start with GHANMGCD
#taxa <- left_join(taxa_B,taxa_B_o)

#### prokka annotation tsv file
name_B<-read_delim("./data/round1/PROKKA_06022021.tsv", 
                 "\t", escape_double = FALSE, trim_ws = TRUE) # id start with GHANMGCD

#### construct results table, filter with padj <= 0.05, abs(log2FoldChange) >= 1
res_B <- results(dds3, contrast = c("Time","24 Hours",'3 Hours'))%>%
  as.data.frame()%>%
  rownames_to_column(var="rowname")%>%
  left_join(taxa_B,by=c("rowname"="ID"))%>%
  left_join(kegg,by=c("KEGG"="K"))
res_B <- res_B[which(abs(res_B$log2FoldChange)>=1 &
             res_B$padj<=0.05),]

res2<-results(dds1, contrast=c("Time","72 Hours","3 Hours"))%>%
  as.data.frame()%>%
  rownames_to_column(var="rowname")%>%
  left_join(taxa_B,by=c("rowname"="ID"))%>%
  left_join(kegg,by=c("KEGG"="K"))
res2 <- res2[which(abs(res2$log2FoldChange)>=1 &
                     res2$padj<=0.05),]

res3<-results(dds1, contrast=c("Time","72 Hours","24 Hours"))%>%
  as.data.frame()%>%
  rownames_to_column(var="rowname")%>%
  left_join(taxa_B,by=c("rowname"="ID"))%>%
  left_join(kegg,by=c("KEGG"="K"))
res3 <- res3[which(abs(res3$log2FoldChange)>=1 &
                     res3$padj<=0.05),]

#### add index to the three results dataframe
res1$Index<-"1"   # for 24h vs 3h
res2$Index<-"2"   # for 72h vs 3h
res3$Index<-"3"   # for 72h vs 24h
resb<-rbind(res1,res2,res3) # construct resb dataframe containing all 3 contrast level results

write_csv(resb, "./data/KEGG/BG_DEGs_time.csv")

#### Heatmap of log2FoldChange for sporulation genes
ggplot(resb[which(resb$Name%like%"sporulation"),],
       aes(x=Index,y=sub("\\;.*", "",Name),fill=log2FoldChange))+
  #geom_tile(aes(fill=log2FoldChange),color="white")+
  geom_tile()+
  theme_bw()+
  labs(x="VB Samples",
       y="Sporulation of Bacillus")


#### vst read counts
df<-kegg[which(kegg$Name%like%"sporulation"),][,1:2]%>%
  unique()%>%
  left_join(taxa_B,by=c("K"="KEGG"))%>%
  na.omit()%>%
  left_join(assay(vsd)%>%as.data.frame()%>%rownames_to_column(),by=c("ID"="rowname"))

pheatmap(df[,4:30]%>%na.omit()%>%as.matrix())

##################
library(tidyr)
library(stringr)
library(radiant.data)
library(data.table)

#resb<-left_join(resb,name,by=c("rowname"="locus_tag"))

# vs_name <- name[grep("chelin",name$product),]
vs<-name_B[grep("bactin",name_B$product),]

vsd_d <- as.data.frame(assay(vsd)) # convert from matrix to dataframe
vsd_d <- rownames_to_column(vsd_d) # convert rownames to a column in the dataframe for easier operation

vsd_id_conversion = left_join(vsd_d, taxa_B, by =c( "rowname" = "ID")) # correct ids 

#vs<-left_join(vs,vsd_id_conversion,by=c("locus_tag"="id"))
#v <- distinct(vs,gene, .keep_all = TRUE) # only keep distinct genes in the table
#vs<-v[which(v$`AVEG24-1`!="NA"),]
vs <- left_join(vs, vsd_d, by=c("locus_tag"="rowname"))
vs <- vs[which(vs$`AVEG24-1`!="NA"),]
vs$Siderophore<-word(vs$product,1)
vs$Siderophore[which(vs$Siderophore%like%"Ferri")]<-"Bacillibactin"
vs$Siderophore[which(vs$Siderophore%like%"Apo")]<-"Petrobactin"
vs$Siderophore[which(vs$Siderophore%like%"Petrobactin")]<-"Petrobactin"
vs$Siderophore[which(vs$Siderophore%like%"Aerobactin")]<-"Aerobactin"
vs$Siderophore[which(vs$product%like%"Ferric enterobactin transport")]<-"Enterobactin"
vs$Siderophore[which(vs$gene%like%"fat")]<-"Anguibactin"
colnames(vs) <- sub("CVB","CVEG",colnames(vs))

vsl<-gather(vs, Sample,Value,8:34)

ggplot(vsl,aes(x =Sample, y = gene))+
  geom_tile(aes(fill=Value),color="white")+
  facet_grid(Siderophore~.,space="free",scale="free")+
  scale_fill_viridis(name="Variance Stablized Count")+
  theme_bw()+
  theme(strip.text.y = element_text(face="bold",angle=0),
        panel.grid.major = element_blank(),
        axis.title = element_text(face="bold"),
        legend.title = element_text(face="bold"),
        axis.text.y = element_text(face="bold"),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text.x = element_text(face="bold",angle=45,hjust=1,vjust=1),
        strip.text.x = element_text(face="bold"),
        legend.text = element_text(face="bold"),
        axis.ticks = element_blank(),
        text = element_text(size=12),
        strip.background = element_rect(fill=NA),
        plot.margin = margin(0.1,0.1,0.1,0.1,"cm"))+
  labs(x="VB Samples",
       y="Siderophore and Iron Acquisition of Bacillus")


