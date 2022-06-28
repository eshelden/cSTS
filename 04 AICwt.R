
library(stringr)
library(ggplot2)
library(readr)
library(SummarizedExperiment)
library(DESeq2)
library(edgeR)
library(biomaRt)
library(gplots)
library(AICcmodavg) #250

#returns a character string containing the first letter of each word.
mFun<-function(x){
  parts<-unlist(str_split(x, " "))
  res<-character(length(parts))
  for (i in 1:length(parts)){
    res[i]<-substr(parts[i],1,1)
  }
  return(paste0(res, collapse = ""))
}


PlotAIC<-function(genes, samples){
  cS.names<-str_remove(genes, "gene-")
  common.genes<-intersect(cS.names,rownames(final.human.data))
  length(common.genes) #for FS, we get 212 genes of the top 250 that we can find human data for
  mData<-as.data.frame(matrix(0,nrow = length(common.genes), ncol = nCols))
  dim (mData) #for FS, 212 x 8, we have 7 TCGA tumor types and the canine (FS, PNST, or PWT data)
  colnames(mData)<-c("cS",codes)
  rownames(mData)<-common.genes
  #Average normalized expression of cs sample gene expression values
  cS.expression <- apply(mCPM[genes,samples],1,mean)
  cS.expression <- unlist(cS.expression)
  names(cS.expression) <- str_remove(names(cS.expression),"gene-")
  mData$cS <- cS.expression[common.genes]
  #average normalized expression of the human tumor type expression values
  for (i in 1:length(TCGA.types)){
    #get data for the tumor type listed.
    Hum.data <- final.human.data[,str_detect(SARC_data@colData$primary_diagnosis, names(TCGA.types)[i])]
    #print(dim(Hum.data))
    Hum.averages <- apply(Hum.data,1,mean)
    Hum.averages <- unlist(Hum.averages)
    mData[,i+1]<- as.numeric(Hum.averages[common.genes])
  }
  
  lm1<-lm(cS ~ Dl, data = mData)
  lm2<-lm(cS ~ F, data = mData)
  lm3<-lm(cS ~ LN, data = mData)
  lm4<-lm(cS ~ Mfh, data = mData)
  lm5<-lm(cS ~ Mpnst, data = mData)
  lm6<-lm(cS ~ Sssc, data = mData)
  lm7<-lm(cS ~ Us, data = mData)
  
  models<-list(lm1,lm2,lm3,lm4,lm5,lm6,lm7)
  model.names<-c("Dl","F","LN","Mfh","Mpnst","Sssc","Us")
  
  my.aictab<-aictab(cand.set = models, modnames = model.names)
  
  
  TCGA.names<-data.frame(name = names(TCGA.types), code = unlist(codes))
  
  mDF<-as.data.frame(my.aictab)
  TCGA.names$height <- mDF[match(TCGA.names$code,mDF$Modnames),]$AICcWt
  
  pDat<-TCGA.names
  colnames(pDat)[3]<-"height"
  pDat$code<-c("DL","FS","LN","MFH","MPNST","SSSC","US")
  
  p<-ggplot(pDat, aes(x=code, y = height)) +
    geom_bar(stat = "identity") + 
    theme_bw() +
    ylim(0,1) +
    theme(
      plot.title = element_text(size=16,face="bold", hjust = 0.5),
      axis.text.y = element_text(size=16,face="bold"),
      #axis.title.y = element_text(size=16,face="bold"),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x = element_text(size=16,face="bold", angle = 35, vjust = 0.7),
      #axis.title.x = element_text(size=16,face="bold"),
      axis.line = element_line(size = 1.5),
      axis.ticks = element_line(size = 1.5),
      legend.text = element_text(face = "bold", size = 16),
      legend.title = element_text(size=16, face="bold")
    )
  print(p)
}

############### If you've completed processing the data, just complete the following, otherwise
############### Complete the processing steps beneath the end of this code section starting on
############### line 155.

setwd("C:\\Users\\idaho\\Dropbox\\Shelden lab\\Canine SARC analysis\\Results\\Tumor type specific gene expression\\Eric")
setwd("C:\\Users\\Eric Shelden\\Dropbox\\Shelden lab\\Canine SARC analysis\\Results\\Tumor type specific gene expression\\Eric")
setwd("C:\\Users\\eshel\\Dropbox\\Shelden lab\\Canine SARC analysis\\Results\\Tumor type specific gene expression\\Eric")

mDir<-getwd()

PNST.F<-read.csv(paste0(mDir,"\\PNST_specific_genes_final.csv"))
FS.F<-read.csv(paste0(mDir,"\\FS_specific_genes_final.csv"))
PWT.F<-read.csv(paste0(mDir,"\\PWT_specific_genes_final.csv"))
dim(PNST.F) #291
dim(FS.F) #2090
dim(PWT.F) #875

FS.F <- FS.F[order(FS.F$padj, decreasing = FALSE),]
PNST.F <- PNST.F[order(PNST.F$padj, decreasing = FALSE),]
PWT.F <- PWT.F[order(PWT.F$padj, decreasing = FALSE),]


load("TCGA_SARC_RLOG_counts.rda")
hNormCounts<-assay(sNormCounts, normalized = TRUE)
rm(sNormCounts)

load("TCGA_data.rda")
FS.samples <- rownames(metaD[metaD$Histology.group == "FS",])
PNST.samples <- rownames(metaD[metaD$Histology.group == "PNST",])
PWT.samples <- rownames(metaD[metaD$Histology.group == "PWT",])

load("mHGeneInfo.rda") #name of all the human genes
dim(mHGeneInfo) #56716 x 5

load("C:\\Users\\idaho\\Dropbox\\Shelden lab\\SARC_analysis\\SARC-data.RData")
load("C:\\Users\\eshel\\Dropbox\\Shelden lab\\SARC_analysis\\SARC-data.RData")
load("C:\\Users\\Eric Shelden\\Dropbox\\Shelden lab\\SARC_analysis\\SARC-data.RData")

mHGeneInfo<-mHGeneInfo[!is.na(mHGeneInfo$external_gene_name),]
rownames(mHGeneInfo)<-make.unique(mHGeneInfo$external_gene_name)
dim(hNormCounts) #56716 x 265
dim(mHGeneInfo) #56404 x 5
length(intersect(rownames(hNormCounts), mHGeneInfo$ensembl_gene_id)) #56404 
final.human.data <- hNormCounts[mHGeneInfo$ensembl_gene_id,] #normalized counts of the human expression data
rownames(final.human.data)<-mHGeneInfo[match(rownames(final.human.data),mHGeneInfo$ensembl_gene_id),]$external_gene_name


TCGA.types<-table(SARC_data@colData$primary_diagnosis)
TCGA.types<-TCGA.types[TCGA.types > 3]
codes<-lapply(names(TCGA.types), mFun)
nCols<-length(TCGA.types)+1

#build data frame for data used by AIC linear models
#nCols
cS.genes<-FS.F[1:250,]$X #NAMES of top 250 expressed FS genes
mSamples<-FS.samples #four FS samples
PlotAIC(cS.genes,mSamples)

cS.genes<-PNST.F[1:250,]$X
mSamples<-PNST.samples
PlotAIC(cS.genes,mSamples)

cS.genes<-PWT.F[1:250,]$X
mSamples<-PWT.samples
PlotAIC(cS.genes,mSamples)

############################################################################################
############################################################################################
#####code below looks at genes that define tumor types....################
setwd("C:\\Users\\idaho\\Dropbox\\Shelden lab\\Canine SARC analysis\\Results\\Tumor type specific gene expression\\Eric")
setwd("C:\\Users\\Eric Shelden\\Dropbox\\Shelden lab\\Canine SARC analysis\\Results\\Tumor type specific gene expression\\Eric")
setwd("C:\\Users\\eshel\\Dropbox\\Shelden lab\\Canine SARC analysis\\Results\\Tumor type specific gene expression\\Eric")

mDir<-getwd()

PNST.F<-read.csv(paste0(mDir,"\\PNST_specific_genes_final.csv"))
FS.F<-read.csv(paste0(mDir,"\\FS_specific_genes_final.csv"))
PWT.F<-read.csv(paste0(mDir,"\\PWT_specific_genes_final.csv"))
dim(PNST.F) #291
dim(FS.F) #2090
dim(PWT.F) #875

FS.F <- FS.F[order(FS.F$padj, decreasing = FALSE),]
PNST.F <- PNST.F[order(PNST.F$padj, decreasing = FALSE),]
PWT.F <- PWT.F[order(PWT.F$padj, decreasing = FALSE),]

#normalized counts.... already done, just load the file on 122 ##############

load("C:\\Users\\idaho\\Dropbox\\Shelden lab\\SARC_analysis\\SARC-data.RData")
load("C:\\Users\\eshel\\Dropbox\\Shelden lab\\SARC_analysis\\SARC-data.RData")
load("C:\\Users\\Eric Shelden\\Dropbox\\Shelden lab\\SARC_analysis\\SARC-data.RData")

SARC.counts<-assay(SARC_data)
m.SARC.counts<-aggregate(SARC.counts, list(rownames(SARC.counts)), max)
rownames(m.SARC.counts)<-m.SARC.counts[,1]
m.SARC.counts<-m.SARC.counts[,-1]

dim(m.SARC.counts) #56716 x 265
dds<-DESeqDataSetFromMatrix(countData = m.SARC.counts, colData = SARC_data@colData, design = ~1)
sNormCounts<-rlog(dds) #takes a long time and a lot of memmory.
save(sNormCounts, file = "TCGA_SARC_RLOG_counts.rda")

load("TCGA_SARC_RLOG_counts.rda")
hNormCounts<-assay(sNormCounts, normalized = TRUE)
rm(sNormCounts)

#we need normalized canine counts too, do this ###############
## OR load the data on line 165 ####
baseDir<-"C:\\Users\\idaho\\Dropbox\\"
baseDir<-"C:\\Users\\Eric Shelden\\Dropbox\\"
baseDir<-"C:\\Users\\eshel\\Dropbox\\"
load(str_c(baseDir,"Shelden lab\\Canine SARC analysis\\Input Data\\all_data.rda"))
mDat<-all_data
colnames(mDat)
mDat<-mDat[,1:16] #just the scrolls
colnames(mDat)

metaDataFile<-"C:\\Users\\Eric Shelden\\Dropbox\\Shelden lab\\Canine SARC analysis\\Input Data\\Subject Data\\STS Data for study-ES_reassigned.csv"
metaDataFile<-"C:\\Users\\idaho\\Dropbox\\Shelden lab\\Canine SARC analysis\\Input Data\\Subject Data\\STS Data for study-ES_reassigned.csv"
metaDataFile<-"C:\\Users\\eshel\\Dropbox\\Shelden lab\\Canine SARC analysis\\Input Data\\Subject Data\\STS Data for study-ES_reassigned.csv"

metaD<-read.csv(metaDataFile)
dim(metaD)
metaD<-metaD[1:16,1:43]
#View(metaD)
rownames(metaD)<-metaD$Sample.ID
#reorder the metaD data
metaD<-metaD[colnames(mDat),]
all.equal(colnames(mDat),metaD$Sample.ID) #TRUE
table(metaD$Histology.group)
#re-assigned: FS: 4, PNST: 6, PWT: 6

#get normalized canine tumor data for plotting:
######### get normalized data for plotting ################
dds<-DESeqDataSetFromMatrix(mDat, metaD, design = ~1)
mCPM<-rlog(dds)
mCPM<-assay(mCPM, normalized = TRUE)

#which samples belong to which tumor type group?
FS.samples <- rownames(metaD[metaD$Histology.group == "FS",])
PNST.samples <- rownames(metaD[metaD$Histology.group == "PNST",])
PWT.samples <- rownames(metaD[metaD$Histology.group == "PWT",])

save(metaD, mCPM, file = "TCGA_data.rda")

load("TCGA_data.rda")

FS.samples <- rownames(metaD[metaD$Histology.group == "FS",])
PNST.samples <- rownames(metaD[metaD$Histology.group == "PNST",])
PWT.samples <- rownames(metaD[metaD$Histology.group == "PWT",])
###############################################################################

#deal with gene names, if already done, just load the file on line 192
ensembl <- useMart("ensembl")
ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)
#mAtt<-attributes(ensembl)
#View(mAtt$attributes)
genemap <- getBM( attributes = c("ensembl_gene_id", "entrezgene_id", "external_gene_name", "description", "transcript_length"),
                  values = rownames(hNormCounts),
                  mart = ensembl )
dim(genemap)

#idx <- match( rownames(SARC.counts), genemap$ensembl_gene_id )


idx <- match( rownames(hNormCounts), genemap$ensembl_gene_id )
length(idx)
mHGeneInfo<-genemap[idx,]
dim(mHGeneInfo) #56716 x 5
length(unique(mHGeneInfo$external_gene_name))
save(mHGeneInfo, file = "mHGeneInfo.rda")

load("mHGeneInfo.rda") #name of all the human genes
dim(mHGeneInfo) #56716 x 5

# ############################################################################

mHGeneInfo<-mHGeneInfo[!is.na(mHGeneInfo$external_gene_name),]
rownames(mHGeneInfo)<-make.unique(mHGeneInfo$external_gene_name)
dim(hNormCounts) #56716 x 265
dim(mHGeneInfo) #56404 x 5
length(intersect(rownames(hNormCounts), mHGeneInfo$ensembl_gene_id)) #56404 
final.human.data <- hNormCounts[mHGeneInfo$ensembl_gene_id,] #normalized counts of the human expression data
rownames(final.human.data)<-mHGeneInfo[match(rownames(final.human.data),mHGeneInfo$ensembl_gene_id),]$external_gene_name

###########################################################################################
load("C:\\Users\\idaho\\Dropbox\\Shelden lab\\SARC_analysis\\SARC-data.RData")
load("C:\\Users\\eshel\\Dropbox\\Shelden lab\\SARC_analysis\\SARC-data.RData")
load("C:\\Users\\Eric Shelden\\Dropbox\\Shelden lab\\SARC_analysis\\SARC-data.RData")

TCGA.types<-table(SARC_data@colData$primary_diagnosis)
TCGA.types<-TCGA.types[TCGA.types > 3]
codes<-lapply(names(TCGA.types), mFun)
nCols<-length(TCGA.types)+1

#build data frame for data used by AIC linear models
#nCols
cS.genes<-FS.F[1:250,]$X #NAMES of top 250 expressed FS genes
mSamples<-FS.samples #four FS samples
PlotAIC(cS.genes,mSamples)

cS.genes<-PNST.F[1:250,]$X
mSamples<-PNST.samples
PlotAIC(cS.genes,mSamples)

cS.genes<-PWT.F[1:250,]$X
mSamples<-PWT.samples
PlotAIC(cS.genes,mSamples)

