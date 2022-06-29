library(stringr)
library(edgeR)
library(DESeq2)
library(clusterProfiler)

#This script does 2 thing:
#calculates optimal comparisons for all tumor types against all 
#normal tissues and all combinations of normal tissues
#Also (about line 187) calculated DEGs for tumor types compared to optimal normal tissues
#and generates lists of all DEGs as well as just protein coding DEGs

########  Load the prepared data ##########################
baseDir<-"C:\\Users\\idaho\\Dropbox\\"
baseDir<-"C:\\Users\\Eric Shelden\\Dropbox\\"
load(str_c(baseDir,"Shelden lab\\Canine SARC analysis\\Input Data\\all_data.rda"))
mDat<-all_data

#Use this Final approach on 7-6-21#
omit<-unlist(str_split("SRR9870570,SRR9870558,SRR9870562,SRR9870560,SRR9870564,SRR9870566,SRR9870572,SRR9870548",","))
length(omit) #8
omit<-match(omit,colnames(mDat))
mDat<-mDat[,-omit]
colnames(mDat)
ncSTS<-16
dim(mDat) #41818 x 56
sTypes = c(rep("cSTS",ncSTS), rep("cVSMS", 9), rep("cSK",6), rep("cCT",18), rep("cNS", 7))
length(sTypes) #56
#retain genes where cpm is greater than 10 for at least 2 tumor and 4 normal tissues
mDat<-mDat[rowSums(cpm(mDat[,1:16]) > 10) > 2,]
mDat<-mDat[rowSums(cpm(mDat[,17:56]) > 10) > 4,] #changed to 4 from 2 on 4-28-21
dim(mDat) #11455 x 56
############################## done loading and preparing ######################


mTypes<-names(table(sTypes)) #list is alphabetical
omit<-match("cSTS",mTypes)
otherTypes<-mTypes[-omit]

#function returns list of lists of all possible combinations of samples
mComb<-function(x){
  out<-list()
  count <- 1
  for (i in 1:length(x)){
    mMat<-combn(x,i)
    for (j in 1:(ncol(mMat))){
      out[[count]]<-paste(mMat[,j])
      count <- count+1
    }
  }
  return(out)
}

all_combinations<-mComb(otherTypes)
length(all_combinations) #15

#get canine STS data for first set of data
metaDataFile<-"C:\\Users\\idaho\\Dropbox\\Shelden lab\\Canine SARC analysis\\Input Data\\Subject Data\\STS Data for study-ES-reassigned.csv"
metaDataFile<-"C:\\Users\\idaho\\Dropbox\\Shelden lab\\Canine SARC analysis\\Input Data\\Subject Data\\STS Data for study-ES-reassigned.csv"
#metaDataFile<-"C:\\Users\\Eric Shelden\\Dropbox\\Shelden lab\\Canine SARC analysis\\Input Data\\Subject Data\\STS Data for study-ES-reassigned.csv"
metaD<-read.csv(metaDataFile)
dim(metaD) #16 x 43
table(metaD$Histology.group) #unassigned, 5,5,6/ 4,6,6 after reassingment

############ Start here for tumor type specific comparisons ##################

cSTS<-mDat[,sTypes %in% "cSTS"] #create dataframe with just canine tumor data
head(cSTS)

#get just FS samples
FS.Samples<-metaD[metaD$Histology.group == "FS",]$Sample.ID
cSTS<-cSTS[,FS.Samples]
head(cSTS)
#or
#get just PNST samples
FS.Samples<-metaD[metaD$Histology.group == "PNST",]$Sample.ID
cSTS<-cSTS[,FS.Samples]
head(cSTS)
#or
#get just PWT samples
FS.Samples<-metaD[metaD$Histology.group == "PWT",]$Sample.ID
cSTS<-cSTS[,FS.Samples]
head(cSTS)

#setup dataframe needed by clusterProfiler for enricher
#sarcoma tumor association gene list downloaded from here:
#https://hpo.jax.org/app/browse/term/HP:0030448
library(readr)
sarcomaGenes <- read_csv("C:/Users/idaho/Dropbox/Shelden lab/Canine SARC analysis/Input Data/Gene Sets/genes_for_HP_0030448.csv")
sarcomaGenes<-sarcomaGenes$GENE_SYMBOL
disease2gene<-data.frame(disease = character(nrow(all_data)), gene = rownames(all_data))
head(disease2gene)
disease2gene$gene<-str_remove(disease2gene$gene, "gene-")
disease2gene$disease<-rep("other",nrow(all_data))
mSarcomaGenes<-sarcomaGenes[sarcomaGenes %in% disease2gene$gene]
disease2gene[match(mSarcomaGenes, disease2gene$gene),]$disease <- "sarcoma"

#run the test. Find significantly different genes using DESEQ

#create matrix to record results
mRes<-matrix(0,nrow = length(all_combinations), ncol = 2)
#records how many DEGs and significance of GSEA for soft tissue sarcoma genes
resMatrix<-matrix(0,nrow = nrow(mDat), ncol = 2*nrow(mRes))

for (i in 1:length(all_combinations)){
  mSearch <- paste0(all_combinations[[i]], collapse = "|")
  print(mSearch) #samples selected
  tDat<-mDat[,grepl(mSearch,sTypes)] #pull data for current selection of samples
  ddDat<-cbind(cSTS,tDat)
  print(colnames(ddDat))
  ddTypes<-c(rep("cSTS",ncol(cSTS)), rep("NT",ncol(tDat))) #types are our tumors vs everything else not tumor (NT)
  #run deseq
  ddSDat<-data.frame(samples = colnames(ddDat), types = ddTypes)
  dds<-DESeqDataSetFromMatrix(ddDat,ddSDat, design = ~types)
  dds<-DESeq(dds)
  res<-results(dds, contrast = c("types", "cSTS", "NT"))
  index<-(2*i)-1
  resMatrix[,index]<-res$padj
  resMatrix[,index+1]<-res$log2FoldChange
  res_sig<-subset(res, padj< .05)
  res_sig<-subset(res_sig, abs(log2FoldChange) >= 1)
  #res_sig<-subset(res_sig, (log2FoldChange) >= 1) #look at just upregulated genes..
  mRes[i,1]<-nrow(res_sig) #save the number of genes where padj < .05 and abs(log2FoldChange >= 1)
  print(nrow(res_sig))
  #this is where we call cluterprofiler
  deg<-rownames(res_sig)
  deg<-str_remove(deg, "gene-")
  x = enricher(deg, TERM2GENE=disease2gene)
  mRes[i,2]<-x@result$p.adjust #save padj value from enrichr function)
  print(x@result$p.adjust)
}

mRes<-as.data.frame(mRes)
mRes<-cbind(paste0(all_combinations), mRes) #add labels to mRes
colnames(mRes)<-c("group","nDEGs","padj")
mRes<-mRes[order(mRes$padj),]
View(mRes) #using the data without reassignment does not reproduce the padj values in the supp data.
           #using the reassigned samples IDs does generate the same values as in supp data

setwd("C:\\Users\\idaho\\Dropbox\\Shelden lab\\Canine SARC analysis\\Results")
write.csv(mRes, file = "FS_DEGs_vs_padj_comparisons_7-12-21.csv")

plot(x = mRes$nDEGs, y = -log10(mRes$padj))
plot(x = seq(1,15,1), y = -log10(mRes$padj))

library(ggplot2)
library(ggrepel)

#make a nice plot:
pData<-mRes
pData$padj <- -log10(pData$padj)
p<-ggplot(pData)
p <- p+geom_point(aes(x = seq(1,15,1), y = padj), size = 5)
p <- p+    theme_bw() +
  theme (
    plot.title = element_text(size=16,face="bold", hjust = 0.5),
    axis.text.y = element_text(size=16,face="bold"),
    axis.title.y = element_text(size=16,face="bold"),
    axis.text.x = element_text(size=16,face="bold"),
    axis.title.x = element_text(size=16,face="bold"),
    axis.line = element_line(size = 1.5),
    axis.ticks = element_line(size = 1.5),
    legend.text = element_text(face = "bold", size = 16),
    legend.title = element_text(size=16, face="bold"),
    panel.border = element_rect(colour = "black", fill=NA, size=2),
    plot.margin = unit(c(1,1,1,1), "cm")
  )
  p <- p + geom_text_repel(aes(x = seq(1,15,1), y = padj, label = group), force = 20, max.overlaps = 20) +
    ylab("-log10(padj)") +
    xlab("")
print(p)


mNames<-character(15)
for (i in 1:15){
  mNames[i]<-paste0(all_combinations[[i]], collapse = "+")
}
allColNames<-character(30)
allColNames[seq(1,29,2)]<-paste0(mNames,"_padj")
allColNames[seq(2,30,2)]<-paste0(mNames,"_l2FC")
colnames(resMatrix)<-allColNames
View(resMatrix)
rownames(resMatrix)<-rownames(mDat)
resMatrix<-as.data.frame(resMatrix)


#7/12/2021 Run DESEQ for each tumor type against the optimum normal data set

#set up things...

#get canine STS data for first set of data
load("C:\\Users\\idaho\\Dropbox\\Shelden lab\\Canine SARC analysis\\Input Data\\geneInfo.rda")
#name of protein coding genes
protein_coding_genes<-mGeneInfo[mGeneInfo$gene_biotype == "protein_coding",]$external_gene_name
#make names the same format as our genes.
protein_coding_genes<-paste0("gene-",protein_coding_genes)

sTypes = c(rep("cSTS",ncSTS), rep("cVSMS", 9), rep("cSK",6), rep("cCT",18), rep("cNS", 7))

metaDataFile<-"C:\\Users\\idaho\\Dropbox\\Shelden lab\\Canine SARC analysis\\Input Data\\Subject Data\\STS Data for study-ES.csv"
#metaDataFile<-"C:\\Users\\Eric Shelden\\Dropbox\\Shelden lab\\Canine SARC analysis\\Input Data\\Subject Data\\STS Data for study-ES.csv"

metaDataFile<-"C:\\Users\\idaho\\Dropbox\\Shelden lab\\Canine SARC analysis\\Input Data\\Subject Data\\STS Data for study-ES.csv"

metaDataFile<-"C:\\Users\\idaho\\Dropbox\\Shelden lab\\Canine SARC analysis\\Input Data\\Subject Data\\STS Data for study-ES_reassigned.csv"
#STS Data for study-ES_reassigned
metaD<-read.csv(metaDataFile)
dim(metaD) #16 x 43
table(metaD$Histology.group) #4,6,6

#Run analysis below

cSTS<-mDat[,sTypes %in% "cSTS"] #create dataframe with just canine tumor data
head(cSTS)

tDat<-mDat[,(grepl("cSK",sTypes))]
head(tDat)
sDat<-c(rep("cSTS",ncol(cSTS)), rep("NT",ncol(tDat)))
dds<-DESeqDataSetFromMatrix(cbind(cSTS,tDat),data.frame(IDs<-c(colnames(cSTS),colnames(tDat)), type = sDat), design = ~ type)
dds<-DESeq(dds)
res<-results(dds, contrast = c("type","cSTS","NT"))
res_sig<-subset(res, padj < .05 & abs(log2FoldChange) > 1)
dim(res_sig) #4387 x 6
res_sig<-as.data.frame(res_sig)
setwd("C:\\Users\\idaho\\Dropbox\\Shelden lab\\Canine SARC analysis\\Results")
write.csv(res_sig, file = "DEGs_all_tumors_vs_cSKIN.csv")
pc_res_sig<-res_sig[rownames(res_sig) %in% protein_coding_genes,]
dim(pc_res_sig)
write.csv(res_sig, file = "PC_DEGs_all_tumors_vs_cSKIN.csv")
#View(res_sig)
#plot(as.numeric(cbind(cSTS,tDat)["gene-RTL1",]))

#get just FS samples vs cSK
cSTS<-mDat[,sTypes %in% "cSTS"] #create dataframe with just canine tumor data
head(cSTS)
FS.Samples<-metaD[metaD$Histology.group == "FS",]$Sample.ID
FS.Samples
cSTS<-cSTS[,FS.Samples]
head(cSTS)
tDat<-mDat[,(grepl("cSK",sTypes))]
head(tDat)
sDat<-c(rep("cSTS",ncol(cSTS)), rep("NT",ncol(tDat)))
dds<-DESeqDataSetFromMatrix(cbind(cSTS,tDat),data.frame(IDs<-c(colnames(cSTS),colnames(tDat)), type = sDat), design = ~ type)
dds<-DESeq(dds)
res<-results(dds, contrast = c("type", "cSTS", "NT"))
res_sig<-subset(res, padj < .05 & abs(log2FoldChange) > 1)
dim(res_sig) #4534 x 6
setwd("C:\\Users\\idaho\\Dropbox\\Shelden lab\\Canine SARC analysis\\Results")
write.csv(res_sig, file = "DEGs_FS_tumors_vs_cSKIN.csv")
pc_res_sig<-res_sig[rownames(res_sig) %in% protein_coding_genes,]
dim(pc_res_sig) #3802 x 6!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
write.csv(pc_res_sig, file = "PC_DEGs_FS_tumors_vs_cSKIN.csv")


#get just PNST samples vs cSK
cSTS<-mDat[,sTypes %in% "cSTS"] #create dataframe with just canine tumor data
head(cSTS)
FS.Samples<-metaD[metaD$Histology.group == "PNST",]$Sample.ID
cSTS<-cSTS[,FS.Samples]
head(cSTS)
tDat<-mDat[,(grepl("cSK",sTypes))]
head(tDat)
sDat<-c(rep("cSTS",ncol(cSTS)), rep("NT",ncol(tDat)))
dds<-DESeqDataSetFromMatrix(cbind(cSTS,tDat),data.frame(IDs<-c(colnames(cSTS),colnames(tDat)), type = sDat), design = ~ type)
dds<-DESeq(dds)
res<-results(dds, contrast = c("type", "cSTS", "NT"))
res_sig<-subset(res, padj < .05 & abs(log2FoldChange) > 1)
dim(res_sig) #4505 x 6
setwd("C:\\Users\\idaho\\Dropbox\\Shelden lab\\Canine SARC analysis\\Results")
write.csv(res_sig, file = "DEGs_PNST_tumors_vs_cSKIN.csv")
pc_res_sig<-res_sig[rownames(res_sig) %in% protein_coding_genes,]
dim(pc_res_sig) #3747 x 6
write.csv(pc_res_sig, file = "PC_DEGs_PNST_tumors_vs_cSKIN.csv")

#get just PWT samples vs vSM
cSTS<-mDat[,sTypes %in% "cSTS"] #create dataframe with just canine tumor data
head(cSTS)
FS.Samples<-metaD[metaD$Histology.group == "PWT",]$Sample.ID
cSTS<-cSTS[,FS.Samples]
head(cSTS)
tDat<-mDat[,(grepl("cVSMS",sTypes))]
head(tDat)
sDat<-c(rep("cSTS",ncol(cSTS)), rep("NT",ncol(tDat)))
dds<-DESeqDataSetFromMatrix(cbind(cSTS,tDat),data.frame(IDs<-c(colnames(cSTS),colnames(tDat)), type = sDat), design = ~ type)
dds<-DESeq(dds)
res<-results(dds, contrast = c("type", "cSTS", "NT"))
res_sig<-subset(res, padj < .05 & abs(log2FoldChange) > 1)
dim(res_sig) #5318 x 6
setwd("C:\\Users\\idaho\\Dropbox\\Shelden lab\\Canine SARC analysis\\Results")
write.csv(res_sig, file = "DEGs_PWT_tumors_vs_cVSMS.csv")
dim(pc_res_sig)
pc_res_sig<-res_sig[rownames(res_sig) %in% protein_coding_genes,]
write.csv(pc_res_sig, file = "PC_DEGs_PWT_tumors_vs_cVSMS.csv")