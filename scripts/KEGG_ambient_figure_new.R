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
# countdata_pseudo <- countdata + 1


############################################################################################################################################
#### Vegetative BG (VEG) vs no BG(N)
metadata1<-metadata[which(metadata$Condition=="Ambient"),]
metadata1<-metadata[which(metadata$CRE=="Y"),]
#countdata1<-countdata_pseudo[,match(metadata1$sample,colnames(countdata_pseudo))]
countdata1<-countdata[,match(metadata1$sample,colnames(countdata))] # pseudo counts added because 
##"Error in estimateSizeFactorsForMatrix(counts(object), locfunc = locfunc,  : 
## every gene contains at least one zero, cannot compute log geometric means"
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
  left_join(taxa_C,by=c("rowname"="ID"))%>%
  left_join(kegg,by=c("KEGG"="K"))

res1_LFC1_padj005 <- res1[which(abs(res1$log2FoldChange)>=1 &
                                  res1$padj<=0.05),]
unique(res1_LFC1_padj005$rowname)%>%length() # 337

res1$Index<-"VB"    # add index "VB" for results for VEG vs N
df1 <- res1[,c(8,3)] %>%unique%>%
  mutate(rank = rank(log2FoldChange,  ties.method = "random"))
df1<-df1[which(df1$log2FoldChange!="NA"),]
df1<-df1[order(-df1$rank),]

l1<-df1$log2FoldChange %>%as.numeric()
names(l1)<-df1$KEGG

l1<-res1$log2FoldChange %>%as.numeric()
names(l1)<-res1$KEGG
l1<-sort(l1, decreasing = TRUE)

set.seed(225)
ke1 <- gseKEGG(geneList     = l1,
               organism     = 'ko',
               pvalueCutoff = 1,
               verbose      = TRUE,
               eps=0,
               seed = TRUE)
# clusterProfiler::dotplot(ke1,showCategory=500)
er1<-ke1@result%>%as.data.frame()
er1$Index<-"VB"  # add index "VB" for results for VEG vs N

############################################################################################################################################
#### Bacillus in APPC (BG) vs no BG (N): res2, df2, ke2, er2
res2<-results(dds1, contrast=c("BG","BG","N"))%>%
  as.data.frame()%>%
  rownames_to_column(var="rowname")%>%
  left_join(taxa_C,by=c("rowname"="ID"))%>%
  left_join(kegg,by=c("KEGG"="K"))

res2_LFC1_padj005 <- res2[which(abs(res2$log2FoldChange)>=1 &
                                  res2$padj<=0.05),]
unique(res2_LFC1_padj005$rowname)%>%length() # 127

res2$Index<-"APPC"

df2 <- res2[,c(8,3)] %>%unique%>%
  mutate(rank = rank(log2FoldChange,  ties.method = "random"))
df2<-df2[which(df2$log2FoldChange!="NA"),]
df2<-df2[order(-df2$rank),]
l2<-df2$log2FoldChange %>%as.numeric()
names(l2)<-df2$KEGG

l2<-res2$log2FoldChange %>%as.numeric()
names(l2)<-res2$KEGG
set.seed(225)
l2<-sort(l2, decreasing = TRUE)
ke2 <- gseKEGG(geneList     = l2,
               organism     = 'ko',
               pvalueCutoff = 1,
               seed = TRUE,
               verbose      = TRUE,
               eps=0)
#clusterProfiler::dotplot(ke2,showCategory=500)
er2<-ke2@result%>%as.data.frame()
er2$Index<-"APPC"

############################################################################################################################################
#### Abiotic filtration (AF) vs no BG (N) 
res3<-results(dds1, contrast=c("BG","AF","N"))%>%
  as.data.frame()%>%
  rownames_to_column(var="rowname")%>%
  left_join(taxa_C,by=c("rowname"="ID"))%>%
  left_join(kegg,by=c("KEGG"="K"))

res3_LFC1_padj005 <- res3[which(abs(res3$log2FoldChange)>=1 &
                                  res3$padj<=0.05),]
unique(res3_LFC1_padj005$rowname)%>%length() # 465
res3$Index<-"APC"

df3 <- res3[,c(8,3)] %>%unique%>%
  mutate(rank = rank(log2FoldChange,  ties.method = "random"))
df3<-df3[which(df3$log2FoldChange!="NA"),]
df3<-df3[order(-df3$rank),]
l3<-df3$log2FoldChange %>%as.numeric()
names(l3)<-df3$KEGG

l3<-res3$log2FoldChange %>%as.numeric()
names(l3)<-res3$KEGG
l3<-sort(l3, decreasing = TRUE)
set.seed(225)
ke3 <- gseKEGG(geneList     = l3,
               organism     = 'ko',
               pvalueCutoff = 1,
               verbose      = TRUE,
               seed = TRUE,
               eps=0)
#clusterProfiler::dotplot(ke3,showCategory=500)
er3<-ke3@result%>%as.data.frame()
er3$Index<-"APC"

############################################################################################################################################
################ Figure 4_A: KEGG Normalized Enrichment Score
#### Read in gene annotation from prokka
name<-read_delim("./KEGG/CRE231_scaffolds.tsv", 
                 "\t", escape_double = FALSE, trim_ws = TRUE)
#### res123: DESeq2 results from 3 contrasts, attach gene information from prokka
res123<-rbind(res1,res2,res3)
res123<-left_join(res123,name,by=c("rowname"="locus_tag"))

#### gseKEGG results from 3 contrasts
er_bio3<-rbind(er1,er2,er3)
er_bio3l<-er_bio3$Description%>%unique() # list of unique kegg pathways (91)

#### ko number of p.adjust <= 0.05
e123<-rbind(er1[which(er1$p.adjust<=0.05),][1],
            er2[which(er2$p.adjust<=0.05),][1],
            er3[which(er3$p.adjust<=0.05),][1])%>%unique()
#### gseKEGG results from 3 contrasts when p.adjust <= 0.05
er_bio3s<-rbind(er1[which(er1$p.adjust<=0.05),],
                er2[which(er2$p.adjust<=0.05),],
                er3[which(er3$p.adjust<=0.05),])
#### Count how many times out of 3 did the ko has p.adjust <=0.05
er_count<-dplyr::count(er_bio3s,er_bio3s$ID)
#### er123s: summary of kegg pathway enrichment results with enrichment n >= 1 out of 3 contrasts
er123s<-left_join(e123,er_bio3,by=c("ID"="ID"))%>%
  left_join(er_count,by=c("ID"="er_bio3s$ID"))

#er_bio3<-left_join(er_bio3,er_count,by=c("ID"="er_bio3$ID"))

#### Add significant level information
er123s$sig<-""
er123s$sig[which(er123s$p.adjust<=0.05)]<-"*"
er123s$sig[which(er123s$p.adjust<=0.01)]<-"**"
er123s$sig[which(er123s$p.adjust<=0.001)]<-"***"

# er123s$NES[which(er123s$NES>3)]<-as.numeric(3) # this result does not have NES > 3

write.csv(er123s,"./KEGG/KEGG_normalized_enrichment_score_new.csv")

p4_A<-ggplot(er123s,aes(x=factor(Index,labels = c("CRE231 + APC",
                                                 "CRE231 + APPC",
                                                 "CRE231 +VB"))))+
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
       y="KEGG Pathway Enrichment")

p4_A

############################################################################################################################################
################ Figure 4_? not sure what this is
l<-read_excel("kegg summary.xlsx")%>%  # do not have this file
  left_join(er123s,by=c("Name"="Description"))

p42<-ggplot(l,aes(x="KEGG Pathway"))+
  geom_tile(aes(y=Name,fill=as.factor(L1)),
            color="white")+
  scale_fill_nord("aurora",name="KEGG Pathway")+
  facet_grid(n~.,scale="free",space="free")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        axis.title = element_text(face="bold"),
        legend.title = element_text(face="bold"),
        axis.text.y = element_blank(),
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

############################################################################################################################################
################ Figure 4_B: CARD results
#### Based on the R processing code, CARD was run with prokka annotated output ffn file, where in CARD output Contig name = prokka geneid
card <- read_delim("./KEGG/CRE231_CARD.txt", 
                   "\t", escape_double = FALSE, trim_ws = TRUE)
card$Contig <- substring(card$Contig,1, nchar(card$Contig)-2 ) # Chop off "_1" at the end of Contig

card123<-left_join(card,res123,by=c("Contig"="rowname"))%>%unique
card123s<-card123[which(card123$log2FoldChange%>%abs()>=1),][,2]%>%unique%>%  # find genes with log2foldchange >= 1
  left_join(card123,by=c("Contig"="Contig"))

#### shorten long names
card123s$Best_Hit_ARO[which(card123s$Best_Hit_ARO=="Klebsiella pneumoniae KpnE")]<-"KpnE"
card123s$Best_Hit_ARO[which(card123s$Best_Hit_ARO=="Klebsiella pneumoniae KpnF")]<-"KpnF"
card123s$Best_Hit_ARO[which(card123s$Best_Hit_ARO=="Trimethoprim-resistant dihydrofolate reductase DfrA42")]<-"DfrA42"
card123s$Best_Hit_ARO[which(card123s$Best_Hit_ARO=="Escherichia coli mdfA")]<-"mdfA"
card123s$Best_Hit_ARO[which(card123s$Best_Hit_ARO=="Escherichia coli UhpT with mutation conferring resistance to fosfomycin")]<-"UhpT"
card123s$Best_Hit_ARO[which(card123s$Best_Hit_ARO=="Escherichia coli marR mutant conferring antibiotic resistance")]<-"marR"

card123s$`Resistance Mechanism`[which(card123s$`Resistance Mechanism`=="antibiotic efflux; reduced permeability to antibiotic")]<-"efflux; reduced permeability"
card123s$`Resistance Mechanism`[which(card123s$`Resistance Mechanism`=="antibiotic target alteration; antibiotic efflux")]<-"target alteration; efflux"


#### add significance level
card123s$sig<-""
card123s$sig[which(card123s$padj<=0.05)]<-"*"
card123s$sig[which(card123s$padj<=0.01)]<-"**"
card123s$sig[which(card123s$padj<=0.001)]<-"***"
write.csv(card123s,"./KEGG/CARD_DESeq_log2foldChange_new.csv")

p4_B<-ggplot(card123s,
            aes(x=factor(Index,labels = c("CRE231 + APC",
                                          "CRE231 + APPC",
                                          "CRE231 + VB")),y=Best_Hit_ARO,fill=log2FoldChange))+
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
p4_B

############################################################################################################################################
################ Figure 4_C: Sulfur metabolism

sulfur<-res123[which(res123$C%like%"Sulfur"&
                       res123$product!="hypothetical protein"),]
sulfur$sig<-""
sulfur$sig[which(sulfur$padj<=0.05)]<-"*"
sulfur$sig[which(sulfur$padj<=0.01)]<-"**"
sulfur$sig[which(sulfur$padj<=0.001)]<-"***"

sulfurs<-sulfur[which(sulfur$log2FoldChange%>%abs()>=1),][1]%>%unique()%>%
  left_join(sulfur,by=c("rowname"="rowname"))
#sulfurs[which(sulfurs$Name%like%"ssuE"),][17]<-"ssuE"
#sulfurs[which(sulfurs$Name%like%"cysE"),][17]<-"cysE"

sulfurs<-sulfurs%>%separate(Name,c("kn","kl"),';')
write.csv(sulfurs,"./KEGG/Sulfur_DESeq_log2foldChange_new.csv")

#p44<-ggplot(sulfurs[which(sulfurs$kn!="ETHE1"),],
p4_C<-ggplot(sulfurs,
            aes(x=factor(Index,labels = c("CRE231 + APC",
                                          "CRE231 + APPC",
                                          "CRE231 + VB")),
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
p4_C

###########
library("cowplot")
library("ggpubr")
pl<-plot_grid(get_legend(p4_A),get_legend(p4_B),nrow=2,
               align="hv",
               rel_heights = c(1.25,1))
pl

pr<-plot_grid(p4_B+rremove("legend"),p4_C+rremove("legend"),
              nrow=2,rel_heights  = c(1,1),
              axis = "tblr",align="v",
              labels = c("E","F"),
              label_x = c(0,0), label_y = c(1,1),
              hjust = c(-0.1,-0.1), vjust = c(1,-1))
pr

plot_grid(p4_A+rremove("legend"),pr,pl,ncol=3,
          axis = "bt",align="v",rel_widths = c(3,2.5,2),
          labels = c("D",NULL,NULL))


############################################################################################################################################
################ Figure S4: Siderophore
sdph<-name[which(name$product%like%"bactin"),][c(1,7)]%>%
  unique()%>%
  left_join(res123,by=c("locus_tag"="rowname"))%>%unique
sdph<-sdph[,c(1,2,4,8,15,18,21)]%>%unique()

sdph$sig<-""
sdph$sig[which(sdph$padj<0.05)]<-"*"
sdph$sig[which(sdph$padj<0.01)]<-"**"
sdph$sig[which(sdph$padj<0.001)]<-"***"

sdph$n[which(sdph$product.y%like%"enterobactin"|
               sdph$product.y%like%"Enterobactin")]<-"Enterobactin"
sdph$n[which(sdph$product.y%like%"anguibactin"|
               sdph$product.y%like%"Anguibactin")]<-"Anguibactin"

sdph$n[which(sdph$product.y%like%"Vibriobactin")]<-"Vibriobactin"
sdph$n[which(sdph$product.y%like%"Petrobactin"|
               sdph$product.y%like%"petrobactin")]<-"Petrobactin"
sdph$n[which(sdph$product.y%like%"aerobactin"|
               sdph$product.y%like%"Aerobactin")]<-"Aerobactin"
write.csv(sdph,"./KEGG/CRE_siderophore_log2foldChange_new.csv")
ps1<-ggplot(sdph,aes(x=factor(Index,labels = c("CRE231 + APC",
                                               "CRE231 + APPC",
                                               "CRE231 + VB"))))+
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
ps1

############################################################################################################################################
################ Figure S4: Type VI secretion system
t6<-res123[which(res123$Name%like%"VI secretion"),]

t6$sig<-""
t6$sig[which(t6$padj<0.05)]<-"*"
t6$sig[which(t6$padj<0.01)]<-"**"
t6$sig[which(t6$padj<0.001)]<-"***"

t6<-t6[which(t6$product != "hypothetical protein"), ]
write.csv(t6,"./KEGG/CRE_TypeVIsecretion_log2foldChange_new.csv")
ps2<-ggplot(t6,aes(x=factor(Index,labels = c("CRE231 + APC",
                                               "CRE231 + APPC",
                                               "CRE231 + VB"))))+
  geom_tile(aes(y=gene ,fill=log2FoldChange),color="white")+
  scale_fill_distiller(palette="Spectral",limits=c(-5.5,5.5),
                       name="log2 Fold Change")+
  geom_text(aes(y=gene,label=sig))+
  theme_bw()+
  #facet_grid(n~.,space="free",scale="free")+
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
       y="Type VI Secretion System (T6SS)")
ps2

############################################################################################################################################
################ Figure S4: Type 1 pili system
pili<-res123[which(res123$Name%like%"Type IV"),] # no type IV pili
pili<-res123[which(res123$Name%like%"pili"),] 
pili<-pili[which(pili$product!="hypothetical protein"),] 

pili$sig<-""
pili$sig[which(pili$padj<0.05)]<-"*"
pili$sig[which(pili$padj<0.01)]<-"**"
pili$sig[which(pili$padj<0.001)]<-"***"
write.csv(pili,"./KEGG/CRE_pili_log2foldChange_new.csv")
ps3<-ggplot(pili,aes(x=factor(Index,labels = c("CRE231 + APC",
                                             "CRE231 + APPC",
                                             "CRE231 + VB"))))+
  geom_tile(aes(y=gene ,fill=log2FoldChange),color="white")+
  scale_fill_distiller(palette="Spectral",limits=c(-5.5,5.5),
                       name="log2 Fold Change")+
  geom_text(aes(y=gene,label=sig))+
  theme_bw()+
  #facet_grid(n~.,space="free",scale="free")+
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
       y="Type 1 Pili")
ps3

############################################################################################################################################
################ Figure S4: poly-beta-1,6-N-acetyl-D-glucosamine and biofilm formation
biof <- res123[which(res123$C%like%"Biofilm form"|res123$Name%like%"poly-beta-1,6-N-acetyl-D-glucosamine"),]%>%unique()
biof_l <- rbind(biof[which(biof$padj<=0.05 & biof$Index == "APC" & abs(biof$log2FoldChange)>=1),][1],
              biof[which(biof$padj<=0.05 & biof$Index == "APPC" & abs(biof$log2FoldChange)>=1),][1],
              biof[which(biof$padj<=0.05 & biof$Index == "VB" & abs(biof$log2FoldChange)>=1),][1])%>%unique()
biof<-left_join(biof_l, biof, by=c("rowname"="rowname"))
biof <- biof[which(biof$product!="hypothetical protein"),] 

biof$sig<-""
biof$sig[which(biof$padj<0.05)]<-"*"
biof$sig[which(biof$padj<0.01)]<-"**"
biof$sig[which(biof$padj<0.001)]<-"***"
write.csv(biof, "./KEGG/CRE_biofilm_log2foldChange_new.csv")
ps4<-ggplot(biof,aes(x=factor(Index,labels = c("CRE231 + APC",
                                               "CRE231 + APPC",
                                               "CRE231 + VB"))))+
  geom_tile(aes(y=gene ,fill=log2FoldChange),color="white")+
  scale_fill_distiller(palette="Spectral",limits=c(-5.5,5.5),
                       name="log2 Fold Change")+
  geom_text(aes(y=gene,label=sig))+
  theme_bw()+
  #facet_grid(n~.,space="free",scale="free")+
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
       y="beta-1,6-N-acetyl-D-glucosamine\npolymer (PGA) and\nBiofilm Formation")
ps4

#### combine 4 figures for supplemental 
pls<-plot_grid(ps3+rremove("legend"), get_legend(ps1),
               nrow=2, 
               labels = c("D", NULL),
               label_x = c(0,0), label_y = c(1,0),
               hjust = c(-0.1,-0.1), vjust = c(1,1))

pls

prs<-plot_grid(ps1+rremove("legend"),ps2+rremove("legend"),ps4+rremove("legend"),
              ncol=3,rel_widths  = c(1.5,1,1.25),
              axis = "tblr",align="h",
              labels = c("A","B","C"),
              label_x = c(0,0,0), label_y = c(1,1,1),
              hjust = c(-0.1,-0.1,-0.1), vjust = c(1,1,1))
prs

plot_grid(prs,pls,ncol=2,
          axis = "bt",align="v",rel_widths = c(3.75,1))
