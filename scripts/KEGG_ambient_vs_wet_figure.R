#### Figures for KEGG enrichment pathways
######### Created by Jinglin Hu for ABBL18 results
######### Notes added by Weitao Shuai; Edited for CRE231 results
rm(list= ls())
library(readr)
library(ggpubr)
library(data.table)
library(pheatmap)
library(clusterProfiler)
library(tidyr)
library(forcats)
library(readxl)
library(DESeq2)
library(corrplot)
library(viridis)
library(clusterProfiler)
library(gplots)
library(phyloseq)
library(dplyr)
library(tibble)
library(reshape2)
library(tximport)
library(tximportData)
library(ggplot2)
library(ggh4x)
library(RColorBrewer)
library(nord)
library(viridis)
library(vegan)
library(stringr)
library(circlize)
library(cowplot)


setwd("/Users/shuaishuai/OneDrive - Northwestern University/Probiotic_Competition/data/")
set.seed(1)

############################################################################################################################################
#### Read in data files for KEGG pathway
taxa_C<-read.delim("./KEGG/CRE231_geneid_ko.txt", header = FALSE, col.names = c('ID',"KEGG")) # taxa_C: prokka geneid vs KEGG ko for CRE 231

metabolism <- readxl::read_xlsx("./round1/KEGG.xlsx",sheet = "Metabolism")
genetic<-readxl::read_xlsx("./round1/KEGG.xlsx",sheet = "Genetic")
env<-readxl::read_xlsx("./round1/KEGG.xlsx",sheet = "Environmental")
cellular<-readxl::read_xlsx("./round1/KEGG.xlsx",sheet = "Cellular")
unclassified<-readxl::read_xlsx("./round1/KEGG.xlsx",sheet = "Sheet1")
brite<-readxl::read_xlsx("./round1/KEGG.xlsx",sheet = "Brite")
others<-readxl::read_xlsx("./round1/KEGG.xlsx",sheet = "Sheet1")
name<-read_delim("./KEGG/CRE231_scaffolds.tsv", 
                 "\t", escape_double = FALSE, trim_ws = TRUE)
rm(kegg)
kegg<-rbind(metabolism,genetic)%>%
  rbind(env)%>%
  rbind(cellular)%>%
  rbind(unclassified)%>%
  rbind(brite)%>%
  rbind(others)

#### Read in data files for DESeq2
metadata <- read_excel("metadata_r2_CRE.xlsx")

# read in feature count output file
data<-read.table("featureCounts-CRE.txt",
                 header = TRUE,skip=1,row.names = 1,
                 check.names = FALSE)
# exclude feature info
countdata<-data[,7:ncol(data)]
# chop off "_hisat-CRE.bam" in column names
colnames(countdata) <- sub("_hisat-CRE.bam","", colnames(countdata))
# chop off "_S12_R1_001" in sample names
metadata$sample <- substring(metadata$sample,1, nchar(metadata$sample)-11 )
# match column names
metadata<-metadata[match(colnames(countdata),metadata$sample),]

#### make metadata factors
metadata$Condition <- factor(metadata$Condition)
metadata$BG <- factor(metadata$BG)
metadata$Time[which(metadata$Time=="A")]<-"3 Hours"
metadata$Time[which(metadata$Time=="B")]<-"24 Hours"
metadata$Time[which(metadata$Time=="C")]<-"72 Hours"
metadata$Time = factor(metadata$Time, 
                       levels=c('3 Hours','24 Hours','72 Hours'))

# add pseudocount 1
#countdata_pseudo <- countdata + 1

############################################################################################################################################
#### No bacillus (N, only CRE231) Wet vs Ambient  CRE231 wet vs CRE231 ambient
metadata2<-metadata[which(metadata$BG=="N"),]

countdata2<-countdata[,match(metadata2$sample,colnames(countdata))]
dds2<-DESeqDataSetFromMatrix(countData = countdata2,
                             colData = metadata2,
                             design=~Condition)

dds2$Condition <- relevel(dds2$Condition, ref = "Ambient")

dds2 <- DESeq(dds2,
              test = "Wald",
              fitType = "parametric", 
              sfType = "ratio")

res4<-results(dds2)%>%
  as.data.frame()%>%
  rownames_to_column(var="rowname")%>%
  left_join(taxa_C,by=c("rowname"="ID"))%>%
  left_join(kegg,by=c("KEGG"="K"))%>%
  left_join(name,by=c("rowname"="locus_tag"))

res4$Index<-"CRE231"

res4_LFC1_padj005 <- res4[which(abs(res4$log2FoldChange)>=1 &
                                  res4$padj<=0.05),]

unique(res4_LFC1_padj005$rowname)%>%length() # 1868

df4 <- res4[,c(8,3)] %>%unique%>%
  mutate(rank = rank(log2FoldChange,  ties.method = "random"))
df4<-df4[which(df4$log2FoldChange!="NA"),]
df4<-df4[order(-df4$rank),]

l4<-df4$log2FoldChange %>%as.numeric()
names(l4)<-df4$KEGG

l4<-res4$log2FoldChange %>%as.numeric()
names(l4)<-res4$KEGG
l4<-sort(l4, decreasing = TRUE)
set.seed(225)
ke4 <- gseKEGG(geneList     = l4,
               organism     = 'ko',
               pvalueCutoff = 1,
               seed=TRUE,
               verbose      = TRUE)

er4<-ke4@result%>%as.data.frame()

er4$Index<-"CRE231 Wet vs Ambient"

############################################################################################################################################
##### CRE231+APPC vs CRE231 under Wet condition
metadata10<-metadata[which(metadata$Condition!="Ambient"),]

countdata10<-countdata[,match(metadata10$sample,colnames(countdata))]
dds12<-DESeqDataSetFromMatrix(countData = countdata10,
                              colData = metadata10,
                              design=~BG+Time)
dds12 <- DESeq(dds12,
               test = "Wald",
               fitType = "parametric", 
               sfType = "ratio")

res39<-results(dds12, contrast=c("BG","BG","N"))%>% 
  as.data.frame()%>%
  rownames_to_column(var="rowname")%>%
  left_join(taxa_C,by=c("rowname"="ID"))%>%
  left_join(kegg,by=c("KEGG"="K"))%>%
  left_join(name,by=c("rowname"="locus_tag"))

res39_LFC1_padj005 <- res39[which(abs(res39$log2FoldChange)>=1 &
                                    res39$padj<=0.05),]
unique(res39_LFC1_padj005$rowname)%>%length() # 21
############################################################################################################################################
#### KEGG pathway enrichment
l39<-res39$log2FoldChange %>%as.numeric()
names(l39)<-res39$KEGG

df39 <- res39[,c(8,3)] %>%unique%>%
  mutate(rank = rank(log2FoldChange,  ties.method = "random"))
df39<-df39[which(df39$log2FoldChange!="NA"),]
df39<-df39[order(-df39$rank),]

l39<-df39$log2FoldChange %>%as.numeric()
names(l39)<-df39$KEGG
l39<-sort(l39, decreasing = TRUE)
set.seed(225)
ke39 <- gseKEGG(geneList     = l39,
                organism     = 'ko',
                pvalueCutoff = 1,
                seed=TRUE,
                verbose      = TRUE,
                eps=0)

er39<-ke39@result%>%as.data.frame()

er4$Index<-"aCRE231"
er39$Index<-"APPC"
er_th<-rbind(er4,er39)

erthp<-er_th[which(er_th$p.adjust<=0.05),]
erth_count<-count(erthp,erthp$ID)
er_th<-left_join(erth_count,er_th,by=c("erthp$ID"="ID")) # for results with p.adjust <=0.05
#### add significance level
er_th$sig<-""
er_th$sig[which(er_th$p.adjust<=0.05)]<-"*"
er_th$sig[which(er_th$p.adjust<=0.01)]<-"**"
er_th$sig[which(er_th$p.adjust<=0.001)]<-"***"
#write.csv(er_th,"./KEGG/wet_vs_ambient_KEGG_normalized_enrichment_score.csv")

#### Figure for KEGG pathway enrichment
p5_A<-ggplot(er_th,aes(x=factor(Index,labels = c("CRE231 Wet vs CRE231 Ambient",
                                                "APPC Wet vs CRE231 Wet"))))+
  geom_tile(aes(y=Description,fill=NES),
            color="white")+
  geom_text(aes(y=Description,label=sig),color="black")+
  scale_fill_distiller(palette = "BrBG",
                       limits=c(-3,3),
                       breaks=c(-3,-1.5,0,1.5,3),
                       name="Normalized Enrichment Score")+
  facet_grid(n~.,scale="free_y",space="free")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        axis.title = element_text(face="bold"),
        legend.title = element_text(face="bold"),
        axis.text.y = element_text(face="bold"),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text.x = element_text(angle=45,
                                   vjust=1,hjust=1,face="bold"),
        strip.text.x = element_text(face="bold"),
        strip.text.y=element_blank(),
        legend.text = element_text(face="bold"),
        axis.ticks = element_blank(),
        text = element_text(size=12),
        plot.margin = margin(0.1,0.1,0.1,0.1,"cm"))+
  labs(x="",
       y="")
p5_A

############################################################################################################################################
#### sulfur metabolism
res4$Index<-"aCRE231"
res39$Index<-"APPC"
res45<-rbind(res4,res39)

res45$sig<-""
res45$sig[which(res45$padj<=0.05)]<-"*"
res45$sig[which(res45$padj<=0.01)]<-"**"
res45$sig[which(res45$padj<=0.001)]<-"***"

sul45<-rbind(res4[which(res4$C%like%"Sulfur"&
                          res4$log2FoldChange>=1),][1],
             res39[which(res39$C%like%"Sulfur"&
                           res39$log2FoldChange>=1),][1])%>%unique()
sul45s<-left_join(sul45,res45,by=c("rowname"="rowname"))%>%
  left_join(name,by=c("rowname"="locus_tag"))
sul45s$sig<-""
sul45s$sig[which(sul45s$padj<=0.05)]<-"*"
sul45s$sig[which(sul45s$padj<=0.01)]<-"**"
sul45s$sig[which(sul45s$padj<=0.001)]<-"***"

sul45s <- sul45s[which(sul45s$product.x!="hypothetical protein"),]
#write.csv(sul45s,"./KEGG/wet_vs_ambient_Sulfur_DESeq_log2foldChange.csv")
#sul45s <- read.csv("C:/Users/jingl/OneDrive/Desktop/CHRISAL/sulfur45.csv", row.names=NULL)

p5_C<-ggplot(sul45s[which(sul45s$gene.x!="NA"),],aes(x=factor(Index,labels = c("CRE231 Wet vs CRE231 Ambient",
                                                                              "APPC Wet vs CRE231 Wet"))))+
  geom_tile(aes(y=gene.x,fill=log2FoldChange),
            color="white")+
  geom_text(aes(y=gene.x,label=sig))+
  scale_fill_distiller(palette = "Spectral",limits=c(-5.5,5.5),
                       name="log2 Fold Change")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        axis.title = element_text(face="bold"),
        legend.title = element_text(face="bold"),
        axis.text.y = element_text(face="bold"),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text.x = element_text(angle=45,vjust=1,hjust=1,
                                   face="bold"),
        text=element_text(size=12),
        strip.text.x = element_text(face="bold"),
        strip.text.y=element_text(angle=0,face="bold"),
        legend.text = element_text(face="bold"),
        axis.ticks = element_blank(),
        strip.background = element_rect(color = "black",fill=NA),
        plot.margin = margin(0.1,0.1,0.1,0.1,unit = "cm"))+
  xlab("")+ylab("Sulfur Metabolism")

p5_C

############################################################################################################################################
#### Type IV Secretion System proteins
write.csv(res45[which(res45$Name%like%"type IV"),][,-c(10,11,12,13)]%>%unique(),"./KEGG/wet_vs_ambient_TypeIVsecretion.csv")
t45<-read.csv("./KEGG/wet_vs_ambient_TypeIVsecretion.csv", row.names=NULL)

t45$sig<-""
t45$sig[which(t45$padj<0.05)]<-"*"
t45$sig[which(t45$padj<0.01)]<-"**"
t45$sig[which(t45$padj<0.001)]<-"***"

t45<-t45%>%separate(Name,c("kn","kl"),',')
t45 <- t45[which(t45$product!="hypothetical protein"),]

unique(t45$rowname)%>%length # 9
unique(t45$kn)%>%length # 7

#### Figure for Type IV secretion system genes
p5_x<-ggplot(t45,aes(x=factor(Index,labels = c("CRE231 Wet vs CRE231 Ambient",
                                              "APPC Wet vs CRE231 Wet"))))+
  geom_tile(aes(y=kn,fill=log2FoldChange),color="white")+
  geom_text(aes(y=kn,label=sig))+
  scale_fill_distiller(palette = "Spectral",limits=c(-5.5,5.5),name="log2 Fold Change")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        axis.title = element_text(face="bold"),
        legend.title = element_text(face="bold"),
        axis.text.y = element_text(face="bold"),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text.x = element_text(angle=45,vjust=1,
                                   hjust=1,face="bold"),
        text=element_text(size=12),
        strip.text.x = element_text(face="bold"),
        strip.text.y=element_text(angle=0,face="bold"),
        legend.text = element_text(face="bold"),
        axis.ticks = element_blank(),
        strip.background = element_rect(color = "black",fill=NA),
        plot.margin = margin(0.1,0.1,0.1,0.1,unit = "cm"))+
  xlab("")+ylab("Type IV Secretion System")
p5_x

############################################################################################################################################
#### Biofilm formation genes
b45<-res45[which(res45$C%like%"form"),]%>%unique()
e45<-rbind(b45[which(b45$padj<=0.05 & b45$Index == "aCRE231" & abs(b45$log2FoldChange)>=1),][1],
           b45[which(b45$padj<=0.05 & b45$Index == "APPC" & abs(b45$log2FoldChange)>=1),][1])%>%unique()

b45<-left_join(e45,b45,by=c("rowname"="rowname"))
b45 <- b45[which(b45$product!="hypothetical protein"),]

b45$sig<-""
b45$sig[which(b45$padj<0.05)]<-"*"
b45$sig[which(b45$padj<0.01)]<-"**"
b45$sig[which(b45$padj<0.001)]<-"***"
#b45<-b45%>%separate(Name,c("kn","kl"),',')

unique(b45$rowname)%>%length # 20
unique(b45$gene)%>%length # 20
write.csv(b45,"./KEGG/wet_vs_ambient_Biofilm.csv")
#### Figure for Biofilm formation genes
p5_y<-ggplot(b45,aes(x=factor(Index,labels = c("CRE231 Wet vs CRE231 Ambient",
                                               "APPC Wet vs CRE231 Wet"))))+
  geom_tile(aes(y=gene,fill=log2FoldChange),color="white")+
  geom_text(aes(y=gene,label=sig))+
  scale_fill_distiller(palette = "Spectral",limits=c(-5.5,5.5),name="log2 Fold Change")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        axis.title = element_text(face="bold"),
        legend.title = element_text(face="bold"),
        axis.text.y = element_text(face="bold"),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text.x = element_text(angle=45,vjust=1,
                                   hjust=1,face="bold"),
        text=element_text(size=12),
        strip.text.x = element_text(face="bold"),
        strip.text.y=element_text(angle=0,face="bold"),
        legend.text = element_text(face="bold"),
        axis.ticks = element_blank(),
        strip.background = element_rect(color = "black",fill=NA),
        plot.margin = margin(0.1,0.1,0.1,0.1,unit = "cm"))+
  xlab("")+ylab("Biofilm Formation")
p5_y

############################################################################################################################################
#### Combining figures
p5l<-plot_grid( get_legend(p5_A),get_legend(p5_y),nrow = 2,align="v",
                rel_heights = c(1,1))

plot_grid(p5_A+rremove("legend"),NULL,p5_y+rremove("legend"),NULL,p5_C+rremove("legend"),NULL,p5l,
          nrow=1,rel_widths  = c(3,0.25,1,0.25,1.15,0.25,2.5),
          axis = "bt",align="vh",
          labels = c("D",NA,"E",NA,"F"))


#### Type 1 pili
pili45<-res45[which(res45$Name%like%"Type IV"),] # no type IV pili
pili45<-res45[which(res45$Name%like%"pili"),] 
pili45<-pili45[which(pili45$product!="hypothetical protein"),] 

pili45$sig<-""
pili45$sig[which(pili45$padj<0.05)]<-"*"
pili45$sig[which(pili45$padj<0.01)]<-"**"
pili45$sig[which(pili45$padj<0.001)]<-"***"

