########HISAT2
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
setwd("C:/Users/jingl/OneDrive/Desktop/CHRISAL//")
set.seed(1)


taxa_A<-readxl::read_xlsx("blastkoala.xlsx",sheet="ABBL18_WO_rRNA")
metadata <- read_excel("metadata.xlsx")
mapping_rate<-read_excel("RNA_Summary.xlsx",sheet="hisat2")


metabolism <- readxl::read_xlsx("C:/Users/jingl/OneDrive/Desktop/CHRISAL/KEGG.xlsx",sheet = "Metabolism")
genetic<-readxl::read_xlsx("C:/Users/jingl/OneDrive/Desktop/CHRISAL/KEGG.xlsx",sheet = "Genetic")
env<-readxl::read_xlsx("C:/Users/jingl/OneDrive/Desktop/CHRISAL/KEGG.xlsx",sheet = "Environmental")
cellular<-readxl::read_xlsx("C:/Users/jingl/OneDrive/Desktop/CHRISAL/KEGG.xlsx",sheet = "Cellular")
unclassified<-readxl::read_xlsx("C:/Users/jingl/OneDrive/Desktop/CHRISAL/KEGG.xlsx",sheet = "Sheet1")
brite<-readxl::read_xlsx("C:/Users/jingl/OneDrive/Desktop/CHRISAL/KEGG.xlsx",sheet = "Brite")
others<-readxl::read_xlsx("C:/Users/jingl/OneDrive/Desktop/CHRISAL/KEGG.xlsx",sheet = "Sheet1")
rm(kegg)
kegg<-rbind(metabolism,genetic)%>%
  rbind(env)%>%
  rbind(cellular)%>%
  rbind(unclassified)%>%
  rbind(brite)%>%
  rbind(others)
#Mapping rate
########
#Mapping Rate 
mr<-left_join(mapping_rate[,1:4],metadata[,c(1,2,4,5)],by=c("Sample"="sample"))%>%
  melt()

ggplot(mr[which(mr$BG!='NA'),],
       aes(x=factor(variable,
                    labels=c("Acinetobacter","Bacillus","Sum")),
           y=value,color=factor(BG,levels=c("N","AF","BG","VEG"),
                                labels=c("ABBL18","ABBL18 + APC",
                                          "ABBL18 + APPC", "ABBL18 + VB"))))+
  geom_boxplot()+
  geom_jitter()+
  facet_grid(.~factor(BG,levels=c("N","AF","BG","VEG"),
                      labels=c("ABBL18","ABBL18 + APC",
                               "ABBL18 + APPC", "ABBL18 + VB")))+
  ylab("Mapping Rate (%)")+
  theme_bw()+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1),
        axis.title.x = element_text(face="bold"),
        axis.title.y = element_text(face="bold"),
        strip.text.x = element_text(face="bold"),
        strip.background = element_rect(fill=NA),
        legend.position = "right",
        legend.title = element_text(face="bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        aspect.ratio = 1)+
  scale_color_brewer(palette="RdYlBu",name="Cleaning Scenario")+
  xlab("Reference Genome")

#########
### Import Count Data
data<-read.table("./hisat/final_featureCounts.txt",
                 header = TRUE,skip=1,row.names = 1,
                 check.names = FALSE)
countdata<-data[,6:ncol(data)]
rownames(metadata)<-metadata$sample
countdata<-countdata[,-19]
metadata<-metadata[match(colnames(countdata),metadata$sample),]
dds<-DESeqDataSetFromMatrix(countData = countdata,
                            colData = metadata,
                            design=~Condition+Time+BG)

dds <- estimateSizeFactors(dds,type="ratio")

unnormalized_counts<-counts(dds,normalize=FALSE)
normalized_counts<-counts(dds,normalize=TRUE)

colSums(normalized_counts != 0)
summary(normalized_counts)

metadata$Time[which(metadata$Time=="A")]<-"3 Hours"
metadata$Time[which(metadata$Time=="B")]<-"24 Hours"
metadata$Time[which(metadata$Time=="C")]<-"72 Hours"
metadata$Time_f = factor(metadata$Time, 
                         levels=c('3 Hours','24 Hours','72 Hours'))

#### QC Plots ####
nc<-melt(normalized_counts)%>%
  left_join(metadata,by=c("Var2"="sample"))
unc<-melt(unnormalized_counts)%>%
  left_join(metadata,by=c("Var2"="sample"))


ggplot(nc)+
  geom_boxplot(aes(x=Var2,y=log2(value+1),
                   color= factor(BG, labels = c("ABBL18 + Abiotic Filtrate", 
                                                "ABBL18 + BactoGreen", 
                                                "ABBL18",
                                                "ABBL18 + Vegetative BactoGreen"))))+
  theme_bw()+
  theme(axis.text.x = element_text(angle=90,vjust = 0.5,size=10,face="bold"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(face="bold"),
        strip.text.x = element_text(size=12,face = "bold"),
        strip.background = element_blank(),
        legend.text = element_text(size=10),
        legend.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))+
  ylab("log2(Normalized Counts)")+
  scale_color_brewer(palette="RdYlBu")+
  facet_nested(.~Condition+Time_f, scales="free",space="free")


ggplot(unc)+
  geom_boxplot(aes(x=Var2,y=log2(value+1),
                   color= factor(BG, labels = c("ABBL18 + Abiotic Filtrate", 
                                                "ABBL18 + BactoGreen", 
                                                "ABBL18",
                                                "ABBL18 + Vegetative BactoGreen"))))+
  theme_bw()+
  theme(axis.text.x = element_text(angle=90,vjust = 0.5,size=10,face="bold"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(face="bold"),
        strip.text.x = element_text(size=12,face = "bold"),
        strip.background = element_blank(),
        legend.text = element_text(size=10),
        legend.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))+
  ylab("log2(Unnormalized Counts)")+
  scale_color_brewer(palette="RdYlBu")+
  facet_nested(.~Condition+Time_f, scales="free",space="free")

############
##### Size Factor ####
x<-sizeFactors(dds)%>%as.data.frame()

x1<-colSums(counts(dds))%>%as.data.frame()
x<-merge(x,x1,by=0)%>%
  left_join(metadata,by=c("Row.names"="sample"))

ggplot(x)+
  geom_point(aes(x=x$..x,y=x$..y,color= factor(BG, labels = c("ABBL18 + Abiotic Filtrate", 
                                                                       "ABBL18 + BactoGreen", 
                                                                       "ABBL18",
                                                                       "ABBL18 + Vegetative BactoGreen"))))+
  theme_bw()+
  theme(axis.text.x = element_text(size=10,face="bold"),
        axis.title.x = element_text(face="bold"),
        axis.title.y = element_text(face="bold"),
        legend.text = element_text(size=10),
        legend.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))+
  ylab("Library Size")+xlab("Size Factor")+
  scale_color_brewer(palette="RdYlBu")


#### Variant stablization and PCoA ####
vsd <- vst(dds, blind=FALSE)
vsdc<-vst(ddsColl, blind=FALSE)

plotPCA(vsd, intgroup=c("BG"))
plotPCA(vsd, intgroup=c("Time"))

pcaData <- plotPCA(vsd, intgroup=c("Condition","Time","BG"), 
                   returnData=TRUE)

percentVar <- round(100 * attr(pcaData, "percentVar"))

P1<-ggplot(pcaData, 
       aes(PC1, PC2, shape=Condition,
           color=factor(BG, 
                        levels = c("N","AF","BG","VEG"),
                        labels = c("ABBL18", 
                                       "ABBL18 + APC", 
                                       "ABBL18 + APPC",
                                       "ABBL18 + VB")))) +
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

###########

meanSdPlot(assay(vsd))


dds <- DESeq(dds,
             test = "Wald",
             fitType = "parametric", 
             sfType = "ratio")

plotDispEsts(dds)

metadata1<-metadata[which(metadata$Condition=="Ambient"),]
countdata1<-countdata[,match(metadata1$sample,colnames(countdata))]
dds1<-DESeqDataSetFromMatrix(countData = countdata1,
                            colData = metadata1,
                            design=~Time+BG)
dds1 <- DESeq(dds1,
              test = "Wald",
              fitType = "parametric", 
              sfType = "ratio")

res1<-results(dds1, contrast=c("BG","VEG","N"))%>%
  as.data.frame()%>%
  rownames_to_column(var="rowname")%>%
  left_join(taxa_A,by=c("rowname"="ID"))%>%
  left_join(kegg,by=c("KEGG"="K"))

sum(res1$padj < 0.05, na.rm=TRUE)
sum(!is.na(res1$padj))

res1_LFC1_padj005 <- res1[which(abs(res1$log2FoldChange)>=1 &
                                  res1$padj<=0.05),]
unique(res1_LFC1_padj005$rowname)%>%length()

res1$Index<-"VB"       

df1 <- res1[,c(8,3)] %>%unique%>%
  mutate(rank = rank(log2FoldChange,  ties.method = "random"))
df1<-df1[which(df1$log2FoldChange!="NA"),]
df1<-df1[order(-df1$rank),]

l1<-df1$log2FoldChange %>%as.numeric()
names(l1)<-df1$KEGG


l1<-res1$log2FoldChange %>%as.numeric()
names(l1)<-res1$KEGG
l1<-sort(l1, decreasing = TRUE)
ke1 <- gseKEGG(geneList     = l1,
               organism     = 'ko',
               pvalueCutoff = 1,
               verbose      = TRUE,
               eps=0)
clusterProfiler::dotplot(ke1,showCategory=500)
er1<-ke1@result%>%as.data.frame()
er1$count<-sapply(strsplit(er1$core_enrichment,'/'), data.table::uniqueN)
er1$Percent<-er1$count/er1$setSize

er1$E[er1$NES<0]<-"Negative"
er1$E[er1$NES>0]<-"Positive"
er1$Index<-"VB"
ggplot(er1, aes(x = Percent, y = fct_reorder(Description, Percent))) + 
  geom_point(aes(size = count, color = p.adjust)) +
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank())+
  facet_grid(E~.,space="free",scale="free")+
  scale_colour_gradient(limits=c(0, 0.05), low="red",high="blue") +
  ylab(NULL)

ribosome1<-res1[which(res1$C%like%"Ribosome"&
                        res1$padj<=0.05&
                        res1$log2FoldChange%>%abs()>=1),][,c(1,3,7)]%>%unique

# BG vs N
res2<-results(dds1, contrast=c("BG","BG","N"))%>%
  as.data.frame()%>%
  rownames_to_column(var="rowname")%>%
  left_join(taxa_A,by=c("rowname"="ID"))%>%
  left_join(kegg,by=c("KEGG"="K"))

sum(res2$padj < 0.05, na.rm=TRUE)

res2_LFC1_padj005 <- res2[which(abs(res2$log2FoldChange)>=1 &
                                  res2$padj<=0.05),]
unique(res2_LFC1_padj005$rowname)%>%length()
res2$Index<-"APPC"
df2 <- res2[,c(8,3)] %>%unique%>%
  mutate(rank = rank(log2FoldChange,  ties.method = "random"))
df2<-df2[which(df2$log2FoldChange!="NA"),]
df2<-df2[order(-df2$rank),]
l2<-df2$log2FoldChange %>%as.numeric()
names(l2)<-df2$KEGG

l2<-res2$log2FoldChange %>%as.numeric()
names(l2)<-res2$KEGG
l2<-sort(l2, decreasing = TRUE)
ke2 <- gseKEGG(geneList     = l2,
               organism     = 'ko',
               pvalueCutoff = 1,
               verbose      = TRUE,
               eps=0)
clusterProfiler::dotplot(ke2,showCategory=500)
er2<-ke2@result%>%as.data.frame()
er2$count<-sapply(strsplit(er2$core_enrichment,'/'), data.table::uniqueN)
er2$Percent<-er2$count/er2$setSize

er2$E[er2$NES<0]<-"Negative"
er2$E[er2$NES>0]<-"Positive"
er2$Index<-"APPC"

ggplot(er2, aes(x = Percent, y = fct_reorder(Description, Percent))) + 
  geom_point(aes(size = count, color = p.adjust)) +
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank())+
  facet_grid(E~.,space="free",scale="free")+
  scale_colour_gradient(limits=c(0, 0.05), low="red",high="blue") +
  ylab(NULL)

l123<-res123[which(abs(res123$log2FoldChange)>3),][,c(1,3)]%>%
  unique()%>%
  as.data.frame()%>%
  left_join(name,by=c("rowname"="locus_tag"))
library(readr)
name<-read_delim("PROKKA_06012021.tsv", 
                 "\t", escape_double = FALSE, trim_ws = TRUE)
res123<-rbind(res1,res2,res3)
res1$Index<-"VB"
res3$Index<-"AF"
res123<-left_join(res123,name,by=c("rowname"="locus_tag"))

ggplot(res123[which(abs(res123$log2FoldChange)>3),])+
  geom_tile(aes(x=Index,y=rowname,fill=log2FoldChange))+
  scale_fill_viridis()

# AF vs N
res3<-results(dds1, contrast=c("BG","AF","N"))%>%
  as.data.frame()%>%
  rownames_to_column(var="rowname")%>%
  left_join(taxa_A,by=c("rowname"="ID"))%>%
  left_join(kegg,by=c("KEGG"="K"))

sum(res3$padj < 0.05, na.rm=TRUE)

res3_LFC1_padj005 <- res3[which(abs(res3$log2FoldChange)>=1 &
                                  res3$padj<=0.05),]

x3<-res3_LFC1_padj005[which(res3_LFC1_padj005$A!="Metabolism"),]%>%
  left_join(name,by=c("rowname"="locus_tag"))

write.csv(x3,"luisa_question_AF.csv")

unique(res3_LFC1_padj005$rowname)%>%length()
res3$Index<-"APC"
res123<-rbind(res1,res2,res3)

df3 <- res3[,c(8,3)] %>%unique%>%
  mutate(rank = rank(log2FoldChange,  ties.method = "random"))
df3<-df3[which(df3$log2FoldChange!="NA"),]
df3<-df3[order(-df3$rank),]
l3<-df3$log2FoldChange %>%as.numeric()
names(l3)<-df3$KEGG

l3<-res3$log2FoldChange %>%as.numeric()
names(l3)<-res3$KEGG
l3<-sort(l3, decreasing = TRUE)
ke3 <- gseKEGG(geneList     = l3,
               organism     = 'ko',
               pvalueCutoff = 1,
               verbose      = TRUE,
               eps=0)
clusterProfiler::dotplot(ke3,showCategory=500)
er3<-ke3@result%>%as.data.frame()
er3$count<-sapply(strsplit(er3$core_enrichment,'/'), data.table::uniqueN)
er3$Percent<-er3$count/er3$setSize

er3$E[er3$NES<0]<-"Negative"
er3$E[er3$NES>0]<-"Positive"

er3$Index<-"APC"
er2$Index<-"APPC"
er_bio3<-rbind(er1,er2,er3)
er_bio3l<-er_bio3$Description%>%unique()

res123p<-res123[which(res123$C=="Ribosome"|
                        res123$C=="Oxidative phosphorylation"|
                        res123$C=="Geraniol degradation" |
                        res123$C%like%"Fatty acid" |
                        res123$C=="Propanoate metabolism"|
                        res123$C%like%"Peptidoglycan"|
                        res123$C== "Valine, leucine and isoleucine degradation"|
                        res123$C%like%"secondary metabolites" |
                        res123$C== "Purine metabolism"|
                        res123$C== "Histidine metabolism" |
                        res123$C== "Carbon metabolism"   |
                        res123$C=="Aminoacyl-tRNA"   |
                        res123$C== "Sulfur metabolism"   |
                        res123$C== "Citrate cycle (TCA cycle)"|           
                        res123$C=="Phenylalanine metabolism"|
                        res123$C=="Biosynthesis of amino acids"),]

res123p$C%>%unique()%>%length()
res123pp<-res123p[which(res123p$log2FoldChange>0),]
pc<-dcast(res123pp, Index ~ C, fun.aggregate = length)%>%t()%>%as.data.frame()
pc$E<-"Positive"
res123pn<-res123p[which(res123p$log2FoldChange<0),]
nc<-dcast(res123pn, Index ~ C, fun.aggregate = length)%>%t()%>%as.data.frame()
nc<-mutate_all(nc, function(x) as.numeric(as.character(x)))
nc[sapply(nc, is.numeric)] <- nc[sapply(nc, is.numeric)] * -1
nc$E<-"Negative"


c<-rbind(pc[-1,]%>%rownames_to_column(),nc[-1,]%>%rownames_to_column())

cl<-gather(c, Index, Value, V1:V3, factor_key=TRUE)

ggplot(cl)+
  geom_bar(aes(x=factor(rowname),y=as.numeric(Value),fill=E),stat = "identity")+
  facet_grid(.~Index,space="free",scale="free")+coord_flip()



er_bio3$Percent1<-er_bio3$Percent*-1
er_bio3$Percent1[which(er_bio3$NES>0)]<-er_bio3$Percent


er_bio3<-left_join(er_bio3,er_count,by=c("ID"="er_bio3$ID"))

ggplot(er_bio3)+
  geom_bar(aes(x=Description,y=as.numeric(Percent1),
               fill=Index),
           stat="identity",position = position_dodge(preserve = "single"))+
  theme_bw()+coord_flip()+
  theme(panel.grid.major = element_blank(),
        axis.title = element_text(face="bold"),
        legend.title = element_text(face="bold"),
        axis.text.y = element_text(face="bold"),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text.x = element_text(size=8),
        strip.text.x = element_text(face="bold"),
        strip.text.y=element_blank())+
  geom_hline(yintercept=0, linetype="dotted")+
  facet_grid(.~Index,space="free",scales = "free",drop = TRUE,
             labeller = label_wrap_gen(1))+
  labs(x="",y="Proportion of Genes Detected within Corresponding Gene Sets")+
  scale_fill_manual(values=c("#D7191C" ,"#FDAE61" ,"#2C7BB6"),
                    name="Surface Cleaning Scenarios")+
  scale_y_continuous(limits = c(-1,1),
                     breaks=c(-1,-0.5,0,0.5,1))

e123<-rbind(er1[which(er1$p.adjust<=0.05),][1],
            er2[which(er2$p.adjust<=0.05),][1],
            er3[which(er3$p.adjust<=0.05),][1])%>%unique()

er123s<-left_join(e123,er_bio3,by=c("ID"="ID"))%>%
  left_join(er_count,by=c("ID"="er_bio3s$ID"))
er_bio3s<-rbind(er1[which(er1$p.adjust<=0.05),],
            er2[which(er2$p.adjust<=0.05),],
            er3[which(er3$p.adjust<=0.05),])
er_count<-dplyr::count(er_bio3s,er_bio3s$ID)

er123s$sig<-""
er123s$sig[which(er123s$p.adjust<=0.05)]<-"*"
er123s$sig[which(er123s$p.adjust<=0.01)]<-"**"
er123s$sig[which(er123s$p.adjust<=0.001)]<-"***"

er123s$NES[which(er123s$NES>3)]<-as.numeric(3)

p41<-ggplot(er123s,aes(x=factor(Index,labels = c("ABBL18 + APC",
                                                 "ABBL18 + APPC",
                                                 "ABBL18 +VB"))))+
  geom_tile(aes(y=Description,fill=log10(`p.adjust`)),
            color="white")+
  geom_text(aes(y=Description,label=sig),color="white")+
  scale_fill_distiller(palette = "Purples",na.value = "#3F007D",
                       limits=c(-5,0),
                       breaks=c(0,-1.3,-2,-3,-4,-5),
                       labels=c("1","0.05","0.01","0.001","0.0001","<0.00001"),
                       name="Adjusted p-value")+
  facet_grid(n.y~.,scale="free_y",space="free")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        axis.title = element_text(face="bold"),
        legend.title = element_text(face="bold"),
        axis.text.y = element_text(face="bold"),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text.x = element_text(face="bold",angle=45,hjust=1,vjust=1),
        strip.text.x = element_text(face="bold"),
        strip.text.y=element_blank(),
        legend.text = element_text(face="bold"),
        axis.ticks = element_blank(),
        text = element_text(size=10),
        plot.margin = margin(0.1,0.1,0.1,0.1,"cm"))+
  labs(x="",
       y="KEGG Pathway")


p41

p42<-ggplot(er123s,aes(x=factor(Index,labels = c("ABBL18 + APC",
                                            "ABBL18 + APPC",
                                            "ABBL18 +VB"))))+
  geom_tile(aes(y=Description,fill=as.numeric(NES)),
            color="white")+
  geom_text(aes(y=Description,label=sig))+
  scale_fill_distiller(palette = "BrBG",
                       limits=c(-3,3),
                       breaks=c(-3,-1.5,0,1.5,3),
                       name="Normalized Enrichment Score")+
  facet_grid(er123s$n~.,scale="free",space="free")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        axis.title = element_text(face="bold"),
        legend.title = element_text(face="bold"),
        axis.text.y = element_text(face="bold"),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text.x = element_text(face="bold",angle=45,hjust=1,vjust=1),
        strip.text.x = element_text(face="bold"),
        strip.text.y=element_blank(),
        text = element_text(size=12),
        legend.text = element_text(face="bold"),
        axis.ticks = element_blank())+
  labs(x="",
       y="")
p42

p42<-ggplot(er123s)+
  geom_bar(aes(x=Description,y=as.numeric(NES),
               fill=factor(Index,labels = c("ABBL18 + APC",
                                            "ABBL18 + APPC",
                                            "ABBL18 +VB"))),
           stat="identity",position = position_dodge(preserve = "single"))+
  theme_bw()+coord_flip()+
  theme(panel.grid.major = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_text(face="bold"),
        legend.title = element_text(face="bold"),
        legend.text = element_text(face="bold"),
        axis.text.y = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        strip.background = element_rect(colour = "black", fill = NA),
        axis.text.x = element_text(size=8),
        strip.text.x = element_text(face="bold"),
        strip.text.y=element_blank())+
  geom_hline(yintercept=0, linetype="dotted")+
  facet_grid(`n.x`~factor(Index,labels = c("ABBL18 + APC",
                                       "ABBL18 + APPC",
                                       "ABBL18 + VB")),space="free",scales = "free",drop = TRUE)+
  labs(x="",y="Normal Enrichment Score")+
  scale_fill_manual(values=c( "#FDAE61" ,"#ABD9E9","#2C7BB6"),
                    name="Cleaning Scenario")+
  scale_y_continuous(limits = c(-3,3))
p42

max(er123s$NES)
library(cowplot)
p4l<-plot_grid(get_legend(p42),get_legend(p44),NULL,nrow=3,
               align="hv",
               rel_heights = c(1,1,1))
p4l


plot_grid(p42+rremove("legend"),pr,p4l,ncol=3,
          axis = "bt",align="v",rel_widths = c(2.5,2.25,1.5),
          labels = c("A",NULL,NULL))

# BG vs AF
res26<-results(dds1, contrast=c("BG","BG","AF"))%>%
  as.data.frame()%>%
  rownames_to_column(var="rowname")%>%
  left_join(taxa_A,by=c("rowname"="ID"))%>%
  left_join(kegg,by=c("KEGG"="K"))

sum(res26$padj < 0.05, na.rm=TRUE)

res26_LFC1_padj005 <- res26[which(abs(res26$log2FoldChange)>=1 &
                                    res26$padj<=0.05),]

unique(res26_LFC1_padj005$rowname)%>%length
l26<-res26$log2FoldChange %>%as.numeric()
names(l26)<-res26$KEGG
l26<-sort(l26, decreasing = TRUE)
ke26<- gseKEGG(geneList     = l26,
               organism     = 'ko',
               pvalueCutoff = 0.05,
               verbose      = TRUE,
               eps=0)
clusterProfiler::dotplot(ke26,showCategory=500)
er26<-ke26@result%>%as.data.frame()
er26$count<-sapply(strsplit(er26$core_enrichment,'/'), data.table::uniqueN)
er26$Percent<-er26$count/er26$setSize

er26$E[er26$NES<0]<-"Negative"
er26$E[er26$NES>0]<-"Positive"
ggplot(er26, aes(x = Percent, y = fct_reorder(Description, Percent))) + 
  geom_point(aes(size = count, color = p.adjust)) +
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank())+
  facet_grid(E~.,space="free",scale="free")+
  scale_colour_gradient(limits=c(0, 0.05), low="red",high="blue") +
  ylab(NULL)


# BG vs VEG
res27<-results(dds1, contrast=c("BG","BG","VEG"))%>%
  as.data.frame()%>%
  rownames_to_column(var="rowname")%>%
  left_join(taxa_A,by=c("rowname"="ID"))%>%
  left_join(kegg,by=c("KEGG"="K"))

sum(res27$padj < 0.05, na.rm=TRUE)

res27_LFC1_padj005 <- res27[which(abs(res27$log2FoldChange)>=1 &
                                    res27$padj<=0.05),]
unique(res27_LFC1_padj005$rowname)%>%length()

### VEG vs AF
res28<-results(dds1, contrast=c("BG","VEG","AF"))%>%
  as.data.frame()%>%
  rownames_to_column(var="rowname")%>%
  left_join(taxa_A,by=c("rowname"="ID"))%>%
  left_join(kegg,by=c("KEGG"="K"))

sum(res28$padj < 0.05, na.rm=TRUE)

res28_LFC1_padj005 <- res28[which(abs(res28$log2FoldChange)>=1 &
                                    res28$padj<=0.05),]

unique(res28_LFC1_padj005$rowname)%>%length
#

#### N Ambient vs Wet
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
  left_join(taxa_A,by=c("rowname"="ID"))%>%
  left_join(kegg,by=c("KEGG"="K"))%>%
  left_join(name,by=c("rowname"="locus_tag"))


sum(res4$padj < 0.05, na.rm=TRUE)
sum(!is.na(res4$padj))

res4$Index<-"ABBL18"

res4_LFC1_padj005 <- res4[which(abs(res4$log2FoldChange)>=1 &
                                  res4$padj<=0.05),]
unique(res4_LFC1_padj005$rowname)%>%length

cf45<-count(res4_LFC1_padj005,res4_LFC1_padj005$C)%>%
  left_join(count(res5_LFC1_padj005,res5_LFC1_padj005$C),
            by=c("res4_LFC1_padj005$C"="res5_LFC1_padj005$C"))


df4 <- res4[,c(8,3)] %>%unique%>%
  mutate(rank = rank(log2FoldChange,  ties.method = "random"))
df4<-df4[which(df4$log2FoldChange!="NA"),]
df4<-df4[order(-df4$rank),]

l4<-df4$log2FoldChange %>%as.numeric()
names(l4)<-df4$KEGG

ggplot(res4[which(res4$C%like%"Ascorbate and aldarate metabolism"),])+
  geom_tile(aes(x=Index,y=rowname,fill=log2FoldChange))+
  scale_fill_distiller(palette = "Spectral",limits=c(-5.5,5.5))

ke4 <- gseKEGG(geneList     = l4,
               organism     = 'ko',
               pvalueCutoff = 1,
               verbose      = TRUE)

clusterProfiler::dotplot(ke4,showCategory=500)
er4<-ke4@result%>%as.data.frame()
er4$count<-sapply(strsplit(er4$core_enrichment,'/'), data.table::uniqueN)
er4$Percent<-er4$count/er4$setSize

er4$E[er4$NES<0]<-"Negative"
er4$E[er4$NES>0]<-"Positive"
er4$Index<-"ABBL18 Wet vs Ambient"

res4_LFC1_padj005p<-res4_LFC1_padj005[which(res4_LFC1_padj005$log2FoldChange>0),]
unique(res4_LFC1_padj005p$rowname)%>%length()

res4_LFC1_padj005n<-res4_LFC1_padj005[which(res4_LFC1_padj005$log2FoldChange<0),]
unique(res4_LFC1_padj005n$rowname)%>%length()

cf4n<-count(res4_LFC1_padj005n,res4_LFC1_padj005n$C)
cf4c<-full_join(cf4n,cf4p,by=c("res4_LFC1_padj005n$C"="res4_LFC1_padj005p$C"))
#### BG Ambient vs Wet
metadata3<-metadata[which(metadata$BG=="BG"),]

countdata3<-countdata[,match(metadata3$sample,colnames(countdata))]
dds3<-DESeqDataSetFromMatrix(countData = countdata3,
                             colData = metadata3,
                             design=~Condition)

dds3 <- DESeq(dds3,
              test = "Wald",
              fitType = "parametric", 
              sfType = "ratio")

res5<-results(dds3, contrast=c("Condition","Wet","Ambient"))%>%
  as.data.frame()%>%
  rownames_to_column(var="rowname")%>%
  left_join(taxa_A,by=c("rowname"="ID"))%>%
  left_join(kegg,by=c("KEGG"="K"))

sum(res5$padj < 0.05, na.rm=TRUE)
sum(!is.na(res5$padj))
res5$Index<-"ABBL18 + APPC"
res45<-rbind(res4,res5)

res5_LFC1_padj005 <- res5[which(abs(res5$log2FoldChange)>=1 &
                                  res5$padj<=0.05),]
unique(res5_LFC1_padj005$rowname)%>%length

l5<-res5$log2FoldChange %>%as.numeric()
names(l5)<-res5$KEGG
l5<-sort(l5, decreasing = TRUE)
ke5 <- gseKEGG(geneList     = l5,
               organism     = 'ko',
               pvalueCutoff = 1,
               verbose      = TRUE,
               eps=0)
clusterProfiler::dotplot(ke5,showCategory=500)
er5<-ke5@result%>%as.data.frame()
er5$count<-sapply(strsplit(er5$core_enrichment,'/'), data.table::uniqueN)
er5$Percent<-er5$count/er5$setSize

er5$E[er5$NES<0]<-"Negative"
er5$E[er5$NES>0]<-"Positive"

er5$Index<-"ABBL18 + APPC"

res5_LFC1_padj005p<-res5_LFC1_padj005[which(res5_LFC1_padj005$log2FoldChange>0),]
cf5p<-count(res5_LFC1_padj005p,res5_LFC1_padj005p$C)

res5_LFC1_padj005n<-res5_LFC1_padj005[which(res5_LFC1_padj005$log2FoldChange<0),]
cf5n<-count(res5_LFC1_padj005n,res5_LFC1_padj005n$C)
cf5c<-full_join(cf5n,cf5p,by=c("res5_LFC1_padj005n$C"="res5_LFC1_padj005p$C"))


unique(res5_LFC1_padj005p$rowname)%>%length()

unique(res5_LFC1_padj005n$rowname)%>%length()

er4$Index<-"ABBL"
er39$Index<-"APPC"
er_th<-rbind(er4,er39)

erthp<-er_th[which(er_th$p.adjust<0.05),]

erth_count<-count(erthp,erthp$ID)
er_th<-left_join(erth_count,er_th,by=c("erthp$ID"="ID"))

er_th$sig<-""
er_th$sig[which(er_th$p.adjust<=0.05)]<-"*"
er_th$sig[which(er_th$p.adjust<=0.01)]<-"**"
er_th$sig[which(er_th$p.adjust<=0.001)]<-"***"

p51<-ggplot(er_th,aes(x=factor(Index,labels = c("ABBL18 Wet vs ABBL18 Ambient",
                                                "APPC Wet vs ABBL18 Wet"))))+
  geom_tile(aes(y=Description,fill=NES),
            color="white")+
  geom_text(aes(y=Description,label=sig),color="black")+
  scale_fill_distiller(palette = "BrBG",
                       limits=c(-3,3),
                       breaks=c(-3,-1.5,0,1.5,3),
                       name="Normalized Enrichment Score (NES)")+
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
p51

res45<-rbind(res4,res39)

res45<-left_join(res45,name,by=c("rowname"="locus_tag"))

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
#write.csv(sul45s,"sulfur45.csv")
sul45s <- read.csv("C:/Users/jingl/OneDrive/Desktop/CHRISAL/sulfur45.csv", row.names=NULL)


p53<-ggplot(sul45s,aes(x=factor(Index,labels = c("ABBL18 Wet vs ABBL18 Ambient",
                                            "APPC Wet vs ABBL18 Wet"))))+
  geom_tile(aes(y=gene,fill=log2FoldChange),
            color="white")+
  geom_text(aes(y=gene,label=sig))+
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
  


ggplot(res45[which(res45$Name%like%"type VI"),])+
  geom_tile(aes(x=Index,y=rowname,fill=log2FoldChange))+
  geom_text(aes(x=Index,y=rowname,label=sig))+
  scale_fill_distiller(palette = "Spectral",limits=c(-5.5,5.5))+
  theme()

#write.csv(res45[which(res45$Name%like%"type IV"),][,-c(10,11,12,13)]%>%
            unique(),"t4_5.csv")
t45<-read.csv("C:/Users/jingl/OneDrive/Desktop/CHRISAL/t4_5.csv", row.names=NULL)

t45$sig<-""
t45$sig[which(t45$padj<0.05)]<-"*"
t45$sig[which(t45$padj<0.01)]<-"**"
t45$sig[which(t45$padj<0.001)]<-"***"

unique(t45$rowname)%>%length
unique(t45$gene)%>%length

p52<-ggplot(t45,aes(x=factor(Index,labels = c("ABBL18 Wet vs ABBL18 Ambient",
                                              "APPC Wet vs ABBL18 Wet"))))+
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
  xlab("")+ylab("Type IV Pili")
p52

s45<-res4[which(res4$C%like%"Sulfur"&abs(res4$log2FoldChange)>=2),]%>%
  rbind(res39[which(res39$C%like%"Sulfur"&abs(res39$log2FoldChange)>=2),])
s45s<-left_join(s45$rowname%>%unique()%>%as.data.frame(),
                res45,by=c("."="rowname"))

#write.csv(s45s[,-c(10,11,12,13)]%>%unique(),"s45.csv")
s45s<-read.csv("C:/Users/jingl/OneDrive/Desktop/CHRISAL/s45.csv", row.names=NULL)


p53<-ggplot(s45s,aes(x=factor(Index,labels = c("ABBL18 Wet vs Ambient",
                                               "Wet APPC vs ABBL18"))))+
  geom_tile(aes(y=gene,fill=log2FoldChange),color="white")+
  geom_text(aes(y=gene,label=sig))+
  scale_fill_distiller(palette = "Spectral",limits=c(-5.5,5.5))+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        axis.title = element_text(face="bold"),
        legend.title = element_text(face="bold"),
        axis.text.y = element_text(face="bold"),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text.x = element_text(angle=45,vjust=1,hjust=1,face="bold"),
        text=element_text(size=12),
        strip.text.x = element_text(face="bold"),
        strip.text.y=element_text(angle=0,face="bold"),
        legend.text = element_text(face="bold"),
        axis.ticks = element_blank(),
        strip.background = element_rect(color = "black",fill=NA),
        plot.margin = margin(0.1,0.1,0,0.1,unit = "cm"))+
  xlab("")+ylab("Type VI Secretion System")

p53

ggplot(res45[which(res45$product%like%"bactin"),])+
  geom_tile(aes(x=Index,y=rowname,fill=log2FoldChange))+
  geom_text(aes(x=Index,y=rowname,label=sig))+
  scale_fill_distiller(palette = "Spectral",limits=c(-5.5,5.5))+
  theme()

ggplot(res45[which(res45$B%like%"secondary"),])+
  geom_tile(aes(x=Index,y=rowname,fill=log2FoldChange))+
  geom_text(aes(x=Index,y=rowname,label=sig))+
  scale_fill_distiller(palette = "Spectral",limits=c(-5.5,5.5))+
  facet_grid(C~.,space="free",scale="free")+
  theme_bw()

ggplot(res45[which(res45$C%like%"Toluene"),])+
  geom_tile(aes(x=Index,y=rowname,fill=log2FoldChange))+
  geom_text(aes(x=Index,y=rowname,label=sig))+
  scale_fill_distiller(palette = "Spectral",limits=c(-5.5,5.5))+
  theme()



card45s<-card45[which(card45$log2FoldChange%>%abs()>=2),][,2]%>%unique%>%
  left_join(card45,by=c("Contig"="Contig"))

card45s$sig<-""
card45s$sig[which(card45s$padj<=0.05)]<-"*"
card45s$sig[which(card45s$padj<=0.01)]<-"**"
card45s$sig[which(card45s$padj<=0.001)]<-"***"
card45s$`Resistance Mechanism`[which(card45s$`Resistance Mechanism`=="antibiotic efflux; reduced permeability to antibiotic")]<-"antibiotic efflux"

card45s$Best_Hit_ARO[which(card45s$Best_Hit_ARO%like%"Mycobacterium")]<-"ethA"
card45s$Best_Hit_ARO[which(card45s$Best_Hit_ARO%like%"fabI")]<-"fabI"
card45s$Best_Hit_ARO[which(card45s$Best_Hit_ARO%like%"AmvA")]<-"AmvA"
card45s$Best_Hit_ARO[which(card45s$Best_Hit_ARO%like%"AbaQ")]<-"AbaQ"
ggplot(card45s,
       aes(x=Index,y=Best_Hit_ARO,fill=log2FoldChange))+
  geom_tile(color="white")+
  facet_grid(`Resistance Mechanism`~.,scale="free",space="free",
             labeller = label_wrap_gen())+
  geom_text(aes(label=sig))+
  scale_fill_distiller(palette = "Spectral",limits=c(-5.5,5.5),
                       name="log2 Fold Change")+theme_bw()+
  theme(panel.grid.major = element_blank(),
        axis.title = element_text(face="bold"),
        legend.title = element_text(face="bold"),
        axis.text.y = element_text(face="bold"),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text.x = element_blank(),
        text=element_text(size=12),
        strip.text.x = element_text(face="bold"),
        strip.text.y=element_text(angle=0,face="bold"),
        legend.text = element_text(face="bold"),
        axis.ticks = element_blank(),
        strip.background = element_rect(color = "black",fill=NA),
        plot.margin = margin(0.5,0.1,0,0.1,unit = "cm"))+
  xlab("")+ylab("Antimicrobial Resistance Gene")


p5l<-plot_grid( get_legend(p51),get_legend(p52),nrow = 2,align="v",
                rel_heights = c(1,1))
p5l

plot_grid(p51+rremove("legend"),NULL,p52+rremove("legend"),NULL,p53+rremove("legend"),NULL,p5l,
                nrow=1,rel_widths  = c(3.35,0.25,1,0.25,1.15,0.25,2.5),
          axis = "bt",align="vh",
                labels = c("A",NA,"B",NA,"C"))



### Time ####

res18<-results(dds, contrast=c("Time","24 Hours","3 Hours"))%>%
  as.data.frame()

sum(res18$padj < 0.05, na.rm=TRUE)

res18_LFC1_padj005 <- res18[which(abs(res18$log2FoldChange)>=1 &
                                    res18$padj<=0.05),]

res19<-results(dds, contrast=c("Time","72 Hours","3 Hours"))%>%
  as.data.frame()

sum(res19$padj < 0.05, na.rm=TRUE)

res19_LFC1_padj005 <- res19[which(abs(res19$log2FoldChange)>=1 &
                                    res19$padj<=0.05),]

res20<-results(dds, contrast=c("Time","72 Hours","24 Hours"))%>%
  as.data.frame()

sum(res20$padj < 0.05, na.rm=TRUE)

res20_LFC1_padj005 <- res20[which(abs(res20$log2FoldChange)>=1 &
                                    res20$padj<=0.05),]
### VEG vs N
res21<-results(dds, contrast=c("BG","VEG","N"))%>%
  as.data.frame()%>%
  rownames_to_column(var="rowname")%>%
  left_join(taxa_A,by=c("rowname"="ID"))%>%
  left_join(kegg,by=c("KEGG"="K"))
res21<-results(dds, contrast=c("BG","VEG","N"))

sum(res21$padj < 0.05, na.rm=TRUE)

res21_LFC1_padj005 <- res21[which(abs(res21$log2FoldChange)>=1 &
                                    res21$padj<=0.05),]

unique(res21_LFC1_padj005$rowname)%>%length

### BG vs N
res22<-results(dds, contrast=c("BG","BG","N"))%>%
  as.data.frame()%>%
  rownames_to_column(var="rowname")%>%
  left_join(taxa_A,by=c("rowname"="ID"))%>%
  left_join(kegg,by=c("KEGG"="K"))

sum(res22$padj < 0.05, na.rm=TRUE)

res22_LFC1_padj005 <- res22[which(abs(res22$log2FoldChange)>=1 &
                                    res22$padj<=0.05),]

### AF vs N
res23<-results(dds, contrast=c("BG","AF","N"))%>%
  as.data.frame()%>%
  rownames_to_column(var="rowname")%>%
  left_join(taxa_A,by=c("rowname"="ID"))%>%
  left_join(kegg,by=c("KEGG"="K"))

sum(res23$padj < 0.05, na.rm=TRUE)

res23_LFC1_padj005 <- res23[which(abs(res23$log2FoldChange)>=1 &
                                    res23$padj<=0.05),]

unique(res21_LFC1_padj005$rowname)%>%length()
unique(res22_LFC1_padj005$rowname)%>%length()
unique(res23_LFC1_padj005$rowname)%>%length()
### BG vs AF
res24<-results(dds, contrast=c("BG","BG","AF"))%>%
  as.data.frame()

sum(res24$padj < 0.05, na.rm=TRUE)

res24_LFC1_padj005 <- res24[which(abs(res24$log2FoldChange)>=1 &
                                    res24$padj<=0.05),]
### BG vs VEG
res25<-results(dds, contrast=c("BG","BG","VEG"))%>%
  as.data.frame()

sum(res25$padj < 0.05, na.rm=TRUE)

res25_LFC1_padj005 <- res25[which(abs(res25$log2FoldChange)>=1 &
                                    res25$padj<=0.05),]


# VEG vs AF
res29<-results(dds, contrast=c("BG","VEG","AF"))%>%
  as.data.frame()

sum(res29$padj < 0.05, na.rm=TRUE)

res29_LFC1_padj005 <- res29[which(abs(res29$log2FoldChange)>=1 &
                                    res29$padj<=0.05),]



###



#### ABBL Time ####
#Ambient 
metadata6<-metadata[which(metadata$BG=="N"&metadata$Condition=="Ambient"),]


countdata6<-countdata[,match(metadata6$sample,colnames(countdata))]
dds8<-DESeqDataSetFromMatrix(countData = countdata6,
                             colData = metadata6,
                             design=~Time)
dds8 <- DESeq(dds8,
              test = "Wald",
              fitType = "parametric", 
              sfType = "ratio")

res30<-results(dds8, contrast=c("Time","24 Hours","3 Hours"))%>%
  as.data.frame()%>%
  rownames_to_column(var="rowname")%>%
  left_join(taxa_A,by=c("rowname"="ID"))%>%
  left_join(kegg,by=c("KEGG"="K"))

sum(res30$padj < 0.05, na.rm=TRUE)

res30_LFC1_padj005 <- res30[which(abs(res30$log2FoldChange)>=1 &
                                    res30$padj<=0.05),]
l30<-res30$log2FoldChange %>%as.numeric()
names(l30)<-res30$KEGG
l30<-sort(l30, decreasing = TRUE)
ke30 <- gseKEGG(geneList     = l30,
                organism     = 'ko',
                pvalueCutoff = 0.05,
                verbose      = TRUE,
                eps=0)
clusterProfiler::dotplot(ke30,showCategory=500)
er30<-ke30@result%>%as.data.frame()
er30$count<-sapply(strsplit(er30$core_enrichment,'/'), data.table::uniqueN)
er30$Percent<-er30$count/er30$setSize

er30$E[er30$NES<0]<-"Negative"
er30$E[er30$NES>0]<-"Positive"


###
res31<-results(dds8, contrast=c("Time","72 Hours","3 Hours"))%>%
  as.data.frame()%>%
  rownames_to_column(var="rowname")%>%
  left_join(taxa_A,by=c("rowname"="ID"))%>%
  left_join(kegg,by=c("KEGG"="K"))

sum(res31$padj < 0.05, na.rm=TRUE)

res31_LFC1_padj005 <- res31[which(abs(res31$log2FoldChange)>=1 &
                                    res31$padj<=0.05),]

l31<-res31$log2FoldChange %>%as.numeric()
names(l31)<-res31$KEGG
l31<-sort(l31, decreasing = TRUE)
ke31 <- gseKEGG(geneList     = l31,
                organism     = 'ko',
                pvalueCutoff = 0.05,
                verbose      = TRUE,
                eps=0)
clusterProfiler::dotplot(ke31,showCategory=500)
er31<-ke31@result%>%as.data.frame()
er31<-ke31@result%>%as.data.frame()
er31$count<-sapply(strsplit(er31$core_enrichment,'/'), data.table::uniqueN)
er31$Percent<-er31$count/er31$setSize

er31$E[er31$NES<0]<-"Negative"
er31$E[er31$NES>0]<-"Positive"


# 
res32<-results(dds8, contrast=c("Time","72 Hours","24 Hours"))%>%
  as.data.frame()%>%
  rownames_to_column(var="rowname")%>%
  left_join(taxa_A,by=c("rowname"="ID"))%>%
  left_join(kegg,by=c("KEGG"="K"))


res32_LFC1_padj005 <- res32[which(abs(res32$log2FoldChange)>=1 &
                                    res32$padj<=0.05),]


unique(res30_LFC1_padj005$rowname)%>%length
unique(res31_LFC1_padj005$rowname)%>%length
unique(res32_LFC1_padj005$rowname)%>%length





# Wet
metadata7<-metadata[which(metadata$BG=="N"&metadata$Condition=="Wet"),]


countdata7<-countdata[,match(metadata7$sample,colnames(countdata))]
dds9<-DESeqDataSetFromMatrix(countData = countdata7,
                             colData = metadata7,
                             design=~Time)

dds9 <- DESeq(dds9,
              test = "Wald",
              fitType = "parametric", 
              sfType = "ratio")

res33<-results(dds9, contrast=c("Time","24 Hours","3 Hours"))%>%
  as.data.frame()%>%
  rownames_to_column(var="rowname")%>%
  left_join(taxa_A,by=c("rowname"="ID"))%>%
  left_join(kegg,by=c("KEGG"="K"))

res33_LFC1_padj005 <- res33[which(abs(res33$log2FoldChange)>=1 &
                                    res33$padj<=0.05),]

l33<-res33$log2FoldChange %>%as.numeric()
names(l33)<-res33$KEGG
l33<-sort(l33, decreasing = TRUE)
ke33 <- gseKEGG(geneList     = l33,
                organism     = 'ko',
                pvalueCutoff = 0.05,
                verbose      = TRUE,
                eps=0)
clusterProfiler::dotplot(ke33,showCategory=500)
er33<-ke33@result%>%as.data.frame()
er33$count<-sapply(strsplit(er33$core_enrichment,'/'), data.table::uniqueN)
er33$Percent<-er33$count/er33$setSize

er33$E[er33$NES<0]<-"Negative"
er33$E[er33$NES>0]<-"Positive"

###
res34<-results(dds9, contrast=c("Time","72 Hours","3 Hours"))%>%
  as.data.frame()%>%
  rownames_to_column(var="rowname")%>%
  left_join(taxa_A,by=c("rowname"="ID"))%>%
  left_join(kegg,by=c("KEGG"="K"))

sum(res34$padj < 0.05, na.rm=TRUE)

res34_LFC1_padj005 <- res34[which(abs(res34$log2FoldChange)>=1 &
                                    res34$padj<=0.05),]
res34_LFC1_padj005$rowname%>%unique()%>%length()
l34<-res34$log2FoldChange %>%as.numeric()
names(l34)<-res34$KEGG
l34<-sort(l34, decreasing = TRUE)
ke34 <- gseKEGG(geneList     = l34,
                organism     = 'ko',
                pvalueCutoff = 0.05,
                verbose      = TRUE,
                eps=0)
clusterProfiler::dotplot(ke34,showCategory=500)
er34<-ke34@result%>%as.data.frame()
er34$count<-sapply(strsplit(er34$core_enrichment,'/'), data.table::uniqueN)
er34$Percent<-er34$count/er34$setSize

er34$E[er34$NES<0]<-"Negative"
er34$E[er34$NES>0]<-"Positive"

res35<-results(dds9, contrast=c("Time","72 Hours","24 Hours"))%>%
  as.data.frame()%>%
  rownames_to_column(var="rowname")%>%
  left_join(taxa_A,by=c("rowname"="ID"))%>%
  left_join(kegg,by=c("KEGG"="K"))


res35_LFC1_padj005 <- res35[which(abs(res35$log2FoldChange)>=1 &
                                    res35$padj<=0.05),]
res35_LFC1_padj005$rowname%>%unique()%>%length()
l35<-res35$log2FoldChange %>%as.numeric()
names(l35)<-res35$KEGG
l35<-sort(l35, decreasing = TRUE)
ke35 <- gseKEGG(geneList     = l35,
                organism     = 'ko',
                pvalueCutoff = 0.05,
                verbose      = TRUE,
                eps=0)
clusterProfiler::dotplot(ke35,showCategory=500)
er35<-ke35@result%>%as.data.frame()

er35$count<-sapply(strsplit(er35$core_enrichment,'/'), data.table::uniqueN)
er35$Percent<-er35$count/er35$setSize
unique(res33_LFC1_padj005$rowname)%>%length
unique(res34_LFC1_padj005$rowname)%>%length
unique(res35_LFC1_padj005$rowname)%>%length


###

## BG TIME 
metadata8<-metadata[which(metadata$BG=="BG"&metadata$Condition=="Ambient"),]


countdata8<-countdata[,match(metadata8$sample,colnames(countdata))]
dds10<-DESeqDataSetFromMatrix(countData = countdata8,
                             colData = metadata8,
                             design=~Time)
dds10 <- DESeq(dds10,
              test = "Wald",
              fitType = "parametric", 
              sfType = "ratio")

res36<-results(dds10, contrast=c("Time","24 Hours","3 Hours"))%>%
  as.data.frame()%>%
  rownames_to_column(var="rowname")%>%
  left_join(taxa_A,by=c("rowname"="ID"))%>%
  left_join(kegg,by=c("KEGG"="K"))

res36_LFC1_padj005 <- res36[which(abs(res36$log2FoldChange)>=1 &
                                  res36$padj<=0.05),]

####
res37<-results(dds10, contrast=c("Time","72 Hours","3 Hours"))%>%
  as.data.frame()%>%
  rownames_to_column(var="rowname")%>%
  left_join(taxa_A,by=c("rowname"="ID"))%>%
  left_join(kegg,by=c("KEGG"="K"))

res37_LFC1_padj005 <- res37[which(abs(res37$log2FoldChange)>=1 &
                                    res37$padj<=0.05),]
#######
res38<-results(dds10, contrast=c("Time","72 Hours","24 Hours"))%>%
  as.data.frame()%>%
  rownames_to_column(var="rowname")%>%
  left_join(taxa_A,by=c("rowname"="ID"))%>%
  left_join(kegg,by=c("KEGG"="K"))


res38_LFC1_padj005 <- res38[which(abs(res38$log2FoldChange)>=1 &
                                    res38$padj<=0.05),]

unique(res36_LFC1_padj005$rowname)%>%length
unique(res37_LFC1_padj005$rowname)%>%length
unique(res38_LFC1_padj005$rowname)%>%length

## BG TIME 
metadata9<-metadata[which(metadata$BG=="BG"&metadata$Condition=="Wet"),]


countdata9<-countdata[,match(metadata9$sample,colnames(countdata))]
dds11<-DESeqDataSetFromMatrix(countData = countdata9,
                              colData = metadata9,
                              design=~Time)
dds11 <- DESeq(dds11,
               test = "Wald",
               fitType = "parametric", 
               sfType = "ratio")

res39<-results(dds11, contrast=c("Time","24 Hours","3 Hours"))%>%
  as.data.frame()%>%
  rownames_to_column(var="rowname")%>%
  left_join(taxa_A,by=c("rowname"="ID"))%>%
  left_join(kegg,by=c("KEGG"="K"))

res39_LFC1_padj005 <- res39[which(abs(res39$log2FoldChange)>=1 &
                                    res39$padj<=0.05),]

####
res40<-results(dds11, contrast=c("Time","72 Hours","3 Hours"))%>%
  as.data.frame()%>%
  rownames_to_column(var="rowname")%>%
  left_join(taxa_A,by=c("rowname"="ID"))%>%
  left_join(kegg,by=c("KEGG"="K"))

res40_LFC1_padj005 <- res40[which(abs(res40$log2FoldChange)>=1 &
                                    res40$padj<=0.05),]
#######
res41<-results(dds11, contrast=c("Time","72 Hours","24 Hours"))%>%
  as.data.frame()%>%
  rownames_to_column(var="rowname")%>%
  left_join(taxa_A,by=c("rowname"="ID"))%>%
  left_join(kegg,by=c("KEGG"="K"))


res41_LFC1_padj005 <- res41[which(abs(res41$log2FoldChange)>=1 &
                                    res41$padj<=0.05),]

unique(res39_LFC1_padj005$rowname)%>%length
unique(res40_LFC1_padj005$rowname)%>%length
unique(res41_LFC1_padj005$rowname)%>%length

bar_df <- read_excel("Bar_plot.xlsx",sheet = "hisat2")
P2<-ggplot(bar_df)+
  geom_bar(aes(x=fct_reorder(Comparison,bar_df$`# DE`),
               y=bar_df$`# DE`,
               fill=bar_df$`Main variable`),
           stat = "identity")+
  facet_grid(.~`Main variable`,space="free",scale="free",
             labeller = label_wrap_gen(20))+
  scale_fill_nord(palette = "aurora")+
  theme_bw()+
  theme(legend.title = element_blank(),
        axis.text.x.bottom = element_text(face="bold"),
        legend.position = "none",
        axis.text.x = element_text(angle=45,vjust=1,hjust=1),
        axis.title = element_text(face="bold"),
        strip.text.x = element_blank(),
        strip.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  xlab("Pairwise Comparison")+ylim(0,1500)+
  ylab("Number of Differentially\nExpressed Genes")

plot_grid(P1,NULL,P2,ncol=3,labels = c("A","B",NA),
          rel_widths = c(2,0.1,1.25),align="hv")

bar_df_time <- read_excel("Bar_plot.xlsx",sheet = "hisat2_2")
ggplot(bar_df_time)+
  geom_bar(aes(x=fct_reorder(Comparison,bar_df_time$`# DE`),
               y=bar_df_time$`# DE`,
               fill=bar_df_time$`Temperature/Humidity`),
           stat = "identity")+
  scale_fill_nord(name="Temperature and Humidity",
                  labels=c("25°C + Ambient Humidity",
                           "37°C + Elevated Humidity"),
                  palette = "aurora")+
  theme_bw()+
  theme(legend.title = element_text(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  xlab("Pairwise Comparison")+
  ylab("Number of Differentially Expressed Genes")+
  coord_flip()+
  facet_grid(bar_df_time$`Temperature/Humidity`~.,
             space="free",scale="free")
####
sampleDists <- dist(t(assay(vsd)))
sampleDists <- dist(t(assay(vsd)[c(4,5,6,13,14,18,22,23,24,40,41,42),
                                 c(4,5,6,13,14,18,22,23,24,40,41,42)]))
sampleDistMatrix <- as.matrix(sampleDists)
#rownames(sampleDistMatrix) <- paste(vsd$Condition, vsd$BG, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Spectral")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors,labels_row=paste0(rownames(sampleDistMatrix)), 
         labels_col=rownames(sampleDistMatrix))

pheatmap(sampleDistMatrix[c(4,5,6,13,14,18,22,23,24,40,41,42),
                          c(4,5,6,13,14,18,22,23,24,40,41,42)],
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

sampleDists <- dist(t(assay(vsdc)))
sampleDistMatrix <- as.matrix(sampleDists)
#rownames(sampleDistMatrix) <- paste(vsd$Condition, vsd$BG, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

corrs <- cor(assay(vsd)[c(4,5,6,13,14,18,22,23,24,40,41,42),
                        c(4,5,6,13,14,18,22,23,24,40,41,42)], method="spearman")

mat<-assay(vsd)
colnames(mat)<-rep(c("A","AAF","ABG","ABG","AVEG","A"),each=9)
mat <- mat[, order(colnames(mat))]
corrs <- cor(mat,method="pearson")
corrf<- cor(mat,method="pearson")
min(corrs)
diag(corrs) <- NA
corr.dists <- as.dist(1 - corrs)

colors <- colorRampPalette(brewer.pal(9, "Spectral"))(520)
pheatmap(corrs, 
         clustering_distance_rows=corr.dists,
         clustering_distance_cols=corr.dists,
         col=colors)
library(ComplexHeatmap)
rownames(corrs)<-colnames(corrs)
corrs[upper.tri(corrs)] <- NA

Heatmap(corrs,cluster_rows = FALSE,cluster_columns  = FALSE,
        show_row_names  = FALSE,show_column_names =  FALSE,
        na_col="white",
        col=colorRamp2(seq(min(corrf), max(corrf), length = 10), 
                       rev(viridis(10))),
        row_title_rot = 0,border = TRUE,
        column_split = rep(c("ABBL18","ABBL18","APC",
                          "APC","APPC",
                          "VB"),each=9),
        row_split = rep(c("ABBL18","ABBL18","APC",
                          "APC","APPC",
                          "VB"),each=9),
        heatmap_legend_param = list(
          title = "Pearson Correlation", 
          at = seq(min(corrf), max(corrf), length = 3)%>%round(2)))

lgd<-Legend(col_fun = colorRamp2(seq(min(corrf), max(corrf), length = 10), 
                            rev(viridis(10))), 
       title = "Pearson Correlation", 
       at = seq(min(corrf), max(corrf), length = 3)%>%round(2))

vsd1<-vst(dds1, blind=FALSE)
mat1<-assay(vsd1)
corrs1 <- cor(mat1,method="pearson")
corrs1[upper.tri(corrs1)] <- NA
corrf1<- cor(mat1,method="pearson")
Heatmap(corrs1,cluster_rows = FALSE,cluster_columns  = FALSE,
        show_row_names  = FALSE,show_column_names =  FALSE,
        na_col="white",
        col=colorRamp2(seq(min(corrf1), max(corrf1), length = 10), 
                       rev(viridis(10))),
        row_title_rot = 0,border = TRUE,
        column_split = rep(c("ABBL18","APC",
                             "APPC",
                             "VB"),each=9),
        row_split = rep(c("ABBL18","APC",
                          "APPC",
                          "VB"),each=9),
        heatmap_legend_param = list(
          title = "Pearson Correlation", 
          at = seq(min(corrf1), max(corrf1), length = 3)%>%round(2)))
########

rm(df2)
df2<-rbind(res1_LFC1_padj005$rowname%>%as.data.frame(),
           res2_LFC1_padj005$rowname%>%as.data.frame())%>%
  rbind(res3_LFC1_padj005$rowname%>%as.data.frame())%>%
  unique()

df2<-left_join(df2, res1[,c(1,3,7)],by=c('.'="rowname"))
colnames(df2)[2:3]<-c("log2FoldChange_VEG","padj_VEG")
df2<-left_join(df2,res2[,c(1,3,7)],by=c('.'="rowname"))
colnames(df2)[4:5]<-c("log2FoldChange_BG","padj_BG")
df2<-left_join(df2,res3[,c(1,3,7)],by=c('.'="rowname"))%>%
  unique()
colnames(df2)[6:7]<-c("log2FoldChange_AF","padj_AF")

df3<-df2[which(df2$padj_VEG<=0.05 &
                 df2$padj_BG<=0.05&
                 df2$padj_AF<=0.05),]
df3<-df2
dfm3<-as.matrix(sapply(df3[,c(2,4,6)], as.numeric))  
rownames(dfm3)<-df3$.
colnames(dfm3)
colnames(dfm3)<-c("APC","APPC","VB")
max(dfm3)
min(dfm3)
col_fun1 = colorRamp2(c(5,0,-5), c("#D7191C", "white", "#2B83BA"))

circos.clear()
dfm4<-as.matrix(sapply(df3[,c(3,5,7)], as.numeric))  
rownames(dfm4)<-df3$.
col_fun2 = colorRamp2(c( 0, 0.05), c("#5E4FA2", "white"))

circos.clear()
circos.par(gap.after = c(4,4,4,4,4,4,25))
circos.heatmap(dfm3,col=col_fun1,dend.side = "outside",
               dend.track.height = 0.25,split=7)

circos.track(track.index = get.current.track.index(), 
             panel.fun = function(x, y) {
               if(CELL_META$sector.numeric.index == 7) { # the last sector
                 cn = colnames(dfm3)
                 n = length(cn)
                 circos.text(rep(CELL_META$cell.xlim[2], n) + convert_x(1, "mm"), 
                             1:n - 0.5, cn, 
                             cex = 1, adj = c(0, 0.5), facing = "inside")
               }
             }, bg.border = NA)
circos.heatmap(dfm4,col=col_fun2)
lgd1 = Legend(title = "log2FC\n", col_fun = col_fun1,
              border = "gray",at=c(-5,-2.5,0,2.5,5),
              legend_height = unit(2.5, "cm"),
              grid_width = unit(0.7, "cm"),
              labels_gp = gpar(fontsize=14),
              title_gp = gpar(fontsize = 16, fontface = "bold"))
lgd2 = Legend(title = "BH-adjusted p-value\n", col_fun = col_fun2,
              border = "gray",
              at = c(0,0.001, 0.01, 0.05),
              legend_height = unit(2.5, "cm"),
              grid_width = unit(0.7, "cm"),
              labels_gp = gpar(fontsize=14),
              title_gp = gpar(fontsize = 16, fontface = "bold"))
h = dev.size()[1]
lgd_list=packLegend(lgd1,lgd2,max_height = unit(50, "in"),
                    row_gap = unit(1, "cm"))
draw(lgd_list, x = unit(7, "mm"), y = unit(20, "mm"), just = c("left", "bottom"))

#### clusters ####
#Cluster 1 VEG+ BG/AF-
cluster1<-df3[which(df3$log2FoldChange_VEG>0&
                      df3$log2FoldChange_BG<0&
                      df3$log2FoldChange_AF<0),]%>%
  left_join(taxa_A,by=c("."="ID"))%>%
  left_join(kegg,by=c("KEGG"="K"))

#Cluster 2 VEG- BG/AF+
cluster2<-df3[which(df3$log2FoldChange_VEG<0&
                      df3$log2FoldChange_BG>0&
                      df3$log2FoldChange_AF>0),]%>%
  left_join(taxa_A,by=c("."="ID"))%>%
  left_join(kegg,by=c("KEGG"="K"))

#Cluster 3 all +
cluster3<-df3[which(df3$log2FoldChange_VEG>0&
                      df3$log2FoldChange_BG>0&
                      df3$log2FoldChange_AF>0),]%>%
  left_join(taxa_A,by=c("."="ID"))%>%
  left_join(kegg,by=c("KEGG"="K"))

#Cluster 4 all -
cluster4<-df3[which(df3$log2FoldChange_VEG<0&
                      df3$log2FoldChange_BG<0&
                      df3$log2FoldChange_AF<0),]%>%
  left_join(taxa_A,by=c("."="ID"))%>%
  left_join(kegg,by=c("KEGG"="K"))

cf1<-count(cluster1, cluster1$C)%>%
  full_join(count(cluster2,cluster2$C),
            by=c("cluster1$C"="cluster2$C"))%>%
  full_join(count(cluster3,cluster3$C),
            by=c("cluster1$C"="cluster3$C"))%>%
  full_join(count(cluster4,cluster4$C),
            by=c("cluster1$C"="cluster4$C"))

colnames(cf1)<-c("C","1","2","3","4")
cf1<-gather(cf1,Cluster,Value,2:5)%>%
  left_join(kegg[,3:5]%>%unique(),by=c("C"="C"))

cf1$Value[which(cf1$Value=="NA")]<-0
cf2<-cf1[which(cf1$A=="Environmental Information Processing"|
                 cf1$A=="Cellular Processes"|cf1$A=="Genetic Information Processing"),]
cf3<-cf1[which(cf1$A=="Metabolism"),]

ggplot(cf2)+
  geom_bar(aes(x=C,y=Value,fill=Cluster),
           stat="identity",
           position = "stack")+coord_flip()+
  theme_bw()+
  facet_grid(A~.,scale="free",space="free")+
  theme(strip.text.y = element_text(angle=0))

ggplot(cf3)+
  geom_bar(aes(x=C,y=Value,fill=Cluster),
           stat = "identity",
           position = "stack")+coord_flip()+
  facet_grid(B~.,space="free",scale="free")+
  theme_bw()+
  theme(strip.text.y = element_text(angle=0))

er_th<-rbind(er4,er5)
er_th$Comparison<-"ABBL18"
er_th$Comparison[2:8]<-"ABBL18 w/ APPC"


er_th$Percent1<-er_th$Percent
er_th$Percent1[which(er_th$NES<0)]<-er_th$Percent*-1

ggplot(er_th[which(er_th$p.adjust<0.05),])+
  geom_bar(aes(x=Description,y=as.numeric(Percent1),fill=Comparison),
           stat="identity",position = position_dodge(preserve = "single"))+
  theme_bw()+coord_flip()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text.x = element_text(angle=30,hjust=1,vjust=1),
        strip.text.y = element_text(angle=0))+
  facet_grid(Comparison~.,scale="free",space="free")+ylim(-1,1)

er30$Comparison<-"24"
er31$Comparison<-"72"
er33$Comparison<-"24"
er34$Comparison<-"72"

er30$Comparison1<-"Ambient"
er31$Comparison1<-"Ambient"
er33$Comparison1<-"Wet"
er34$Comparison1<-"Wet"

er_time<-rbind(er30,er31,er33,er34)

er_time$Percent1<-er_time$Percent*-1
er_time$Percent1[which(er_time$NES>0)]<-er_time$Percent
er_tc<-count(er_time,er_time$ID)
er_time<-left_join(er_time,er_tc,by=c("ID"="er_time$ID"))

ggplot(er_time)+
  geom_bar(aes(x=Description,y=as.numeric(Percent1),
               fill=factor(er_time$Comparison, 
                           labels=c('24 vs 3 hours','72 vs 3 hours'))),
           stat="identity",position = position_dodge(preserve = "single"))+
  theme_bw()+coord_flip()+
  theme(panel.grid.major = element_blank(),
        axis.title = element_text(face="bold"),
        legend.title = element_text(face="bold"),
        axis.text.y = element_text(face="bold"),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text.x = element_text(size=8),
        strip.text.x = element_text(face="bold"),
        strip.text.y=element_blank())+
  geom_hline(yintercept=0, linetype="dotted")+
  facet_grid(n~Comparison1,space="free",scales = "free",drop = TRUE)+
  labs(x="",y="Gene Ratio")+
  scale_fill_nord(palette = "silver_mine",name="Time Comparison")+ 
  labs(x="",y="Proportion of Genes Detected within Corresponding Gene Sets")+
  scale_y_continuous(limits = c(-1,1),
                     breaks=c(-1,-0.5,0,0.5,1))

library(readr)
card <- read_delim("2021-06-17_16_53_23.774_PROKKA_06012021.part_001.ffn.txt", 
                                                                    "\t", escape_double = FALSE, trim_ws = TRUE)

card123<-left_join(card,res123,by=c("Contig"="rowname"))%>%unique
card45<-left_join(card,res45,c("Contig"="rowname"))%>%unique

card123s$Best_Hit_ARO[which(card123s$Best_Hit_ARO=="Acinetobacter baumannii AbaQ")]<-"AbaQ"
card123s$Best_Hit_ARO[which(card123s$Best_Hit_ARO=="Acinetobacter baumannii AbaF")]<-"AbaF"
card123s$Best_Hit_ARO[which(card123s$Best_Hit_ARO%like%"ethA")]<-"ethA"


p43<-ggplot(card123s,
       aes(x=factor(Index,labels = c("ABBL18 + APC",
                                     "ABBL18 + APPC",
                                     "ABBL18 + VB")),y=Best_Hit_ARO,fill=log2FoldChange))+
  geom_tile(color="white")+
  facet_grid(`Resistance Mechanism`~.,scale="free",space="free",
             labeller = label_wrap_gen())+
  geom_text(aes(label=sig))+
  scale_fill_distiller(palette = "Spectral",limits=c(-5.5,5.5),
                       name="log2 Fold Change")+theme_bw()+
  theme(panel.grid.major = element_blank(),
        axis.title = element_text(face="bold"),
        legend.title = element_text(face="bold"),
        axis.text.y = element_text(face="bold"),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text.x = element_blank(),
        text=element_text(size=12),
        strip.text.x = element_text(face="bold"),
        strip.text.y=element_text(angle=0,face="bold"),
        legend.text = element_text(face="bold"),
        axis.ticks = element_blank(),
        strip.background = element_rect(color = "black",fill=NA),
        plot.margin = margin(0.5,0.1,0,0.1,unit = "cm"))+
  xlab("")+ylab("Antibiotic Resistance")
p43

card123s<-card123[which(card123$log2FoldChange%>%abs()>=2),][,2]%>%unique%>%
  left_join(card123,by=c("Contig"="Contig"))

card123s$sig<-""
card123s$sig[which(card123s$padj<=0.05)]<-"*"
card123s$sig[which(card123s$padj<=0.01)]<-"**"
card123s$sig[which(card123s$padj<=0.001)]<-"***"

max(card123$log2FoldChange)
min(card123$log2FoldChange)

sulfur<-res123[which(res123$C%like%"Sulfur"),]
sulfur$sig<-""
sulfur$sig[which(sulfur$padj<=0.05)]<-"*"
sulfur$sig[which(sulfur$padj<=0.01)]<-"**"
sulfur$sig[which(sulfur$padj<=0.001)]<-"***"

sulfurs<-sulfur[which(sulfur$log2FoldChange%>%abs()>=2),][1]%>%unique()%>%
  left_join(sulfur,by=c("rowname"="rowname"))
sulfurs[which(sulfurs$Name%like%"ssuE"),][17]<-"ssuE"
sulfurs[which(sulfurs$Name%like%"cysE"),][17]<-"cysE"

sulfurs<-sulfurs%>%separate(Name,c("kn","kl"),';')


p44<-ggplot(sulfurs[which(sulfurs$kn!="ETHE1"),],
       aes(x=factor(Index,labels = c("ABBL18 + APC",
                                     "ABBL18 + APPC","ABBL18 + VB")),
                    y=kn))+
  geom_tile(aes(fill=log2FoldChange),color="white")+
  geom_text(aes(label=sig))+
  scale_fill_distiller(palette = "Spectral",limits=c(-5.5,5.5),
                       name="log2 Fold Change")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        text = element_text(size=12),
        axis.title.x = element_blank(),
        axis.title.y = element_text(face="bold"),
        legend.title = element_text(face="bold"),
        axis.text.y = element_text(face="bold"),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text.x = element_text(face="bold",angle=45,hjust=1,vjust=1),
        strip.text.x = element_text(face="bold"),
        strip.text.y=element_text(angle=0,face="bold"),
        legend.text = element_text(face="bold"),
        axis.ticks = element_blank(),
        plot.margin = margin(0,0,0.5,0,"cm"),
        strip.background = element_rect(color = "black",fill=NA))+
  ylab("Sulfur Metabolism")
p44

min(sulfur$log2FoldChange)

pr<-plot_grid(p43+rremove("legend"),p44+rremove("legend"),
              nrow=2,rel_heights  = c(1.25,1),
              axis = "tblr",align="v",
              labels = c("B","C"))
pr

library(readxl)
Gates_Run1 <- read_excel("Gates_Run1.xlsx", 
                         sheet = "Sheet2")
rna<-left_join(Gates_Run1[,1:2],metadata,by=c("Sample"="sample"))
ggplot(rna[which(rna$BG!="NA"),])+
  geom_boxplot(aes(x=BG,y=Concentration,
                   color=factor(Time,
                                levels = c("3 Hours", "24 Hours","72 Hours"))))+
  facet_grid(.~Condition,space="free",scale="free")

ggplot(rna[which(rna$BG!="NA"),])+
  geom_bar(aes(x=factor(BG,labels = c("APC","APPC","ABBL18","VB")),y=Concentration,
                  fill=factor(Time,
                                levels = c("3 Hours", "24 Hours","72 Hours"))),
           stat="identity",position = "dodge")+
  facet_grid(.~Condition,space="free",scale="free")+
  theme_bw()+
  theme(axis.title.x = element_blank(),
        legend.title = element_blank())

P00010<- res1[res1$C %like% "Glycolysis", ]%>%as.data.frame()
pheatmap(P00010[,-1],cluster_cols = FALSE)
ggplot(res123[which(res123$C %like% "Glycolysis"),])+
  geom_tile(aes(x=Index,y=Name,fill=log2FoldChange))+
  scale_fill_distiller(palette = "Spectral")


res123m<-spread(res123[,c(1,3,14)]%>%unique(), Index, log2FoldChange)
res123mm<-left_join(res123m,res123[,c(1,9:12)],by=c("rowname"="rowname"))%>%
  na.omit%>%unique()%>%
  left_join(name[,c(1,4)],by=c("rowname"="locus_tag"))
pheatmap(res123mm[which(res123mm$C%like% "Biofilm"),][,2:4]%>%na.omit(),
         labels_row =res123mm[which(res123mm$C%like% "Biofilm"),][,9]%>%na.omit() )

pheatmap(res123mm[which(res123mm$C%like% "Glycolysis"),][,2:4]%>%na.omit(),
         labels_row =res123mm[which(res123mm$C%like% "Glycolysis"),][,9] )
pheatmap(res123mm[which(res123mm$C%like% "Citrate"),][,2:4]%>%na.omit(),
         labels_row =res123mm[which(res123mm$C%like% "Citrate"),][,9] )
pheatmap(res123mm[which(res123mm$C%like% "Fatty acid degradation"),][,2:4]%>%na.omit(),
         labels_row =res123mm[which(res123mm$C%like% "Fatty acid degradation"),][,9] )
pheatmap(res123mm[which(res123mm$C%like% "Oxidative"),][,2:4]%>%na.omit(),
         labels_row =res123mm[which(res123mm$C%like% "Oxidative"),][,9] )
pheatmap(res123mm[which(res123mm$C%like% "Sulfur"),][,2:4]%>%na.omit(),
         labels_row =res123mm[which(res123mm$C%like% "Sulfur"),][,9] ,
         breaks = seq(-5.5,5.5,length.out=101))
res123mm[which(res123mm$C%like% "Sulfur"),][,2:4]%>%min

pheatmap(res123mm[which(res123mm$C%like% "Aminoacyl"),][,2:4]%>%na.omit(),
         labels_row =res123mm[which(res123mm$C%like% "Aminoacyl"),][,9] )
pheatmap(res123mm[which(res123mm$C%like% "xenobiotic"),][,2:4]%>%na.omit(),
         labels_row =res123mm[which(res123mm$C%like% "xenobiotic"),][,9] )

pheatmap(res123mm[which(res123mm$C%like% "siderophore"),][,2:4]%>%na.omit(),
         labels_row =res123mm[which(res123mm$C%like% "siderophore"),][,9] )

#pheatmap(res123mm[which(res123mm$C%like% "Quorum"),][,2:4]%>%na.omit(),
         labels_row =res123mm[which(res123mm$C%like% "Quorum"),][,9] )
pheatmap(res123mm[which(res123mm$C%like% "Ribosome"),][,2:4]%>%na.omit(),
         labels_row =res123mm[which(res123mm$C%like% "Ribosome"),][,9] )
pheatmap(res123mm[which(res123mm$C%like% "RNA degr"),][,2:4]%>%na.omit(),
         labels_row =res123mm[which(res123mm$C%like% "RNA degr"),][,9] )
pheatmap(res123mm[which(res123mm$C%like% "RNA polymerase"),][,2:4]%>%na.omit(),
         labels_row =res123mm[which(res123mm$C%like% "RNA polymerase"),][,9] )
pheatmap(res123mm[which(res123mm$C%like% "Protein ex"),][,2:4]%>%na.omit(),
         labels_row =res123mm[which(res123mm$C%like% "Protein ex"),][,9] )
pheatmap(res123mm[which(res123mm$C%like% "secre"),][,2:4]%>%na.omit(),
         labels_row =res123mm[which(res123mm$C%like% "secre"),][,9] ,
         limits=c(-5.5,5.5))

biof<-name[which(name$gene%like%"pga"),][c(1,7)]%>%
  unique()%>%
  left_join(res123,
            by=c("locus_tag"="rowname"))%>%unique

biof<-res123[which(res123$Name%like%"pga"),]%>%unique
biof$sig<-""
biof$sig[which(biof$padj<0.05)]<-"*"
biof$sig[which(biof$padj<0.01)]<-"**"
biof$sig[which(biof$padj<0.001)]<-"***"

biof$gene[which(biof$gene=='pgaB')]<-"pgaB_1"
biof$gene[which(is.na(biof$gene)==TRUE & biof$Name%like%"pgaB")]<-"pgaB_2"
biof$gene[which(is.na(biof$gene)==TRUE & biof$Name%like%"pgaD")]<-"pgaD"


biof <- biof[order(biof$gene),]

ps2<-ggplot(biof,aes(x=factor(Index,labels = c("ABBL18 + APC",
                                          "ABBL18 + APPC",
                                          "ABBL18 + VB"))))+
  geom_tile(aes(y=gene,fill=log2FoldChange),color="white")+
  scale_fill_distiller(palette = "Spectral",limits=c(-5.5,5.5))+
  geom_text(aes(y=gene,label=sig))+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        axis.title = element_text(face="bold"),
        legend.title = element_text(face="bold"),
        axis.text.y = element_text(face="bold"),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text.x = element_text(face="bold",angle=45,hjust=1,vjust=1),
        strip.text.x = element_text(face="bold"),
        strip.text.y=element_text(angle = 0,face="bold"),
        text = element_text(size=12),
        legend.position = "none",
        legend.text = element_text(face="bold"),
        axis.ticks = element_blank(),
        strip.background = element_rect(fill=NA,color="black"))+
  labs(x="",
       y="??-1,6-N-acetyl-D-glucosamine\npolymer (PGA) and\nBiofilm Formation")

ps2

type6<-res123[which(res123$Name%like%"type VI"),]%>%unique
write.csv(type6[,-c(10,11,12,13)]%>%unique(),"type6.csv")
type6 <- read.csv("C:/Users/jingl/OneDrive/Desktop/CHRISAL/type6.csv", row.names=NULL)

type6$sig<-""
type6$sig[which(type6$padj<0.05)]<-"*"
type6$sig[which(type6$padj<0.01)]<-"**"
type6$sig[which(type6$padj<0.001)]<-"***"

ps3<-ggplot(type6,aes(x=factor(Index,labels = c("ABBL18 + APC",
                                           "ABBL18 + APPC",
                                           "ABBL18 + VB"))))+
  geom_tile(aes(y=gene,fill=log2FoldChange),color="white")+
  scale_fill_distiller(palette = "Spectral",limits=c(-5.5,5.5)) +
  geom_text(aes(y=gene,label=sig))+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        axis.title = element_text(face="bold"),
        legend.title = element_text(face="bold"),
        axis.text.y = element_text(face="bold"),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text.x = element_text(face="bold",angle=45,hjust=1,vjust=1),
        strip.text.x = element_text(face="bold"),
        strip.text.y=element_text(angle = 0,face="bold"),
        text = element_text(size=12),
        legend.position = "none",
        legend.text = element_text(face="bold"),
        axis.ticks = element_blank(),
        strip.background = element_rect(fill=NA,color="black"))+
  labs(x="",
       y="Type VI Secretion System (T6SS)")

ps3

t4<-res123[which(res123$Name%like%"type IV"),]%>%unique
#write.csv(t4[,-c(10,11,12,13)]%>%unique(),"t4.csv")
t4 <- read.csv("C:/Users/jingl/OneDrive/Desktop/CHRISAL/t4.csv", row.names=NULL)

t4$sig<-""
t4$sig[which(t4$padj<0.05)]<-"*"
t4$sig[which(t4$padj<0.01)]<-"**"
t4$sig[which(t4$padj<0.001)]<-"***"

ps4<-ggplot(t4,aes(x=factor(Index,labels = c("ABBL18 + APC",
                                           "ABBL18 + APPC",
                                           "ABBL18 + VB"))))+
  geom_tile(aes(y=gene,fill=log2FoldChange),color="white")+
  scale_fill_distiller(palette = "Spectral",limits=c(-5.5,5.5))+
  geom_text(aes(y=gene,label=sig))+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        axis.title = element_text(face="bold"),
        legend.title = element_text(face="bold"),
        axis.text.y = element_text(face="bold"),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text.x = element_text(angle = 45,face="bold",hjust = 1,vjust = 1),
        strip.text.x = element_text(face="bold"),
        strip.text.y=element_text(angle = 0,face="bold"),
        text = element_text(size=12),
        legend.position = "none",
        legend.text = element_text(face="bold"),
        axis.ticks = element_blank(),
        strip.background = element_rect(fill=NA,color="black"))+
  labs(x="",
       y="Type IV Pili")



sdph<-name[which(name$product%like%"bactin"),][c(1,7)]%>%
  unique()%>%
  left_join(res123,by=c("locus_tag"="rowname"))%>%unique


sdph<-sdph[,c(1,2,4,8,15,18,21)]%>%unique()

sdph$sig<-""
sdph$sig[which(sdph$padj<0.05)]<-"*"
sdph$sig[which(sdph$padj<0.01)]<-"**"
sdph$sig[which(sdph$padj<0.001)]<-"***"

sdph$n<-'NA'
sdph$n[which(sdph$product.y%like%"enterobactin")]<-"Enterobactin"
sdph$n[which(sdph$product.y%like%"anguibactin")]<-"Anguibactin"
sdph$n[which(sdph$product.y%like%"Vibriobactin")]<-"Vibriobactin"
sdph$n[which(sdph$product.y%like%"Petrobactin")]<-"Petrobactin"
sdph$n[which(sdph$product.y%like%"aerobactin")]<-"Aerobactin"

ps1<-ggplot(sdph,aes(x=factor(Index,labels = c("ABBL18 + APC",
                                          "ABBL18 + APPC",
                                          "ABBL18 + VB"))))+
  geom_tile(aes(y=gene ,fill=log2FoldChange),color="white")+
  scale_fill_distiller(palette="Spectral",limits=c(-5.5,5.5),
                       name="log2 Fold Change")+
  geom_text(aes(y=gene,label=sig))+
  theme_bw()+
  facet_grid(n~.,space="free",scale="free")+
  theme(panel.grid.major = element_blank(),
        axis.title = element_text(face="bold"),
        legend.title = element_text(face="bold"),
        axis.text.y = element_text(face="bold"),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text.x = element_text(angle=45,hjust = 1,vjust = 1,face="bold"),
        strip.text.x = element_text(face="bold"),
        strip.text.y=element_text(angle = 0,face="bold"),
        text = element_text(size=12),
        legend.text = element_text(face="bold"),
        axis.ticks = element_blank(),
        strip.background = element_rect(fill=NA,color="black"))+
  labs(x="",
       y="Siderophore and\nIron Acquisition")

max(sdph$log2FoldChange)
min(sdph$log2FoldChange)

permanova <- adonis(t(counts(dds,normalized=TRUE)) ~ Condition+BG+Time,
                    data = metadata, permutations=9999, method = "euclidean")

permanova <- adonis(assay(vsd)%>%t() ~ Condition+BG+Time,
                    data = metadata, permutations=9999, method = "euclidean")

permanova <- adonis(assay(rld)%>%t() ~ Condition+BG+Time,
                    data = metadata, permutations=9999, method = "euclidean")

permanova$aov.tab
##

ps2n<-plot_grid(NULL,ps2,get_legend(ps1),nrow=3,rel_heights = c(0.1,1,0.5))
plot_grid(ps1+rremove("legend"),ps3,ps4,ps2n,ncol=4,
          rel_widths = c(1.5,1.5,1,1.35))



####
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
  left_join(taxa_A,by=c("rowname"="ID"))%>%
  left_join(kegg,by=c("KEGG"="K"))

res39_LFC1_padj005 <- res39[which(abs(res39$log2FoldChange)>=1 &
                                    res39$padj<=0.05),]
l39<-res39[,c(8,3)]%>%unique()
l39s<-l39$log2FoldChange %>%as.numeric()
names(l39s)<-l39$KEGG
l39s<-sort(l39s, decreasing = TRUE)
ke39 <- gseKEGG(geneList     = l39s,
                organism     = 'ko',
                pvalueCutoff = 1,
                verbose      = TRUE,
                eps=0)

er39<-ke39@result%>%as.data.frame()

res39$Index<-"APPC vs N"
res45<-rbind(res4,res39)

ggplot(res45[which(res45$C%like%"Secretion syste"),])+
  geom_tile(aes(y=rowname,fill=log2FoldChange,x=Index))+
  scale_fill_distiller(palette="Spectral",limits=c(-5.5,5.5))
ggplot(res45[which(res45$C%like%"Sulfur"),])+
  geom_tile(aes(y=rowname,fill=log2FoldChange,x=Index))+
  scale_fill_distiller(palette="Spectral",limits=c(-5.5,5.5))


ggplot(res123[which(res123$B%like%"repair"),])+
  geom_tile(aes(y=rowname,fill=log2FoldChange,x=Index))+
  scale_fill_distiller(palette="Spectral",limits=c(-5.5,5.5))+
  facet_grid(C~.,space="free",scale="free")+
  theme(strip.text.y = element_text(angle=0))
