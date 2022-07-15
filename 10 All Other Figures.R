library(FactoMineR)
library(factoextra)
library(DESeq2)
library(SummarizedExperiment)
library(ggrepel)
library(ggplot2)
library(stringr)
library(edgeR)

########################### Figure S1 PCA of replicates for first scrolls ####################################
load("C:\\Users\\idaho\\Dropbox\\Shelden lab\\Canine SARC analysis\\raw_reads\\CanFam-Newest\\FirstScrolls\\first_scrolls_replicates.rda")
tail(mDat)
mDat<-mDat[1:(nrow(mDat)-5),]
dds<-DESeqDataSetFromMatrix(mDat,colData = as.data.frame(colnames(mDat)), design = ~ 1)
mCPM<-rlog(dds)
mCPM<-assay(mCPM, normalized = TRUE)
mVar<-apply(mCPM,1,var)
mCPM<-mCPM[order(mVar, decreasing = TRUE),]
head(mCPM)
mPCA<-PCA(t(mCPM[1:500,]))
plot.PCA(mPCA)
pDat<-mPCA$ind$coord
pIDs<-unlist(str_split("10,10,10,11,11,11,12,12,12,2,2,2,3,3,3,4,4,4,5,5,5,7,7,7,9,9,9",","))

p<-ggplot(data = as.data.frame(pDat),aes(x = Dim.1, y = Dim.2))
p+geom_point(col = as.factor(pIDs), size = 3) +
  geom_text_repel(label = rownames(pDat), force = 5, size = 6) +
  theme_bw()

################# PCA clustering of original tumor data, Figure SF2A ################################
#PCA of top 1000 most variable genes

load("C:\\Users\\idaho\\Dropbox\\Shelden lab\\Canine SARC analysis\\Input Data\\all_data.rda")
mDat<-all_data
colnames(mDat)
mDat<-mDat[,1:16] #analyze just the cSTS
colnames(mDat)

metaDataFile<-"C:\\Users\\idaho\\Dropbox\\Shelden lab\\Canine SARC analysis\\Input Data\\Subject Data\\STS Data for study-original.csv"
metaD<-read.csv(metaDataFile)
dim(metaD)
metaD<-metaD[1:16,1:43]
rownames(metaD)<-metaD$Sample.ID
metaD<-metaD[colnames(mDat),]
all.equal(colnames(mDat),metaD$Sample.ID) #TRUE
table(metaD$Histology.group)
#initially, FS, PNST, PWT = 5, 5, 6

dds<-DESeqDataSetFromMatrix(mDat,colData = as.data.frame(colnames(mDat)), design = ~ 1)
mCPM<-rlog(dds)
mCPM<-assay(mCPM, normalized = TRUE)
mL2<-mCPM
mVar<-apply(mL2,1,var)
mL2<-mL2[order(mVar, decreasing = TRUE),]
sTypes<-as.factor(metaD$Histology.group)
MyPCAPlot(mL2,computeLogs = FALSE,text = TRUE,ngenes = 1000,sTypes = sTypes) #MyPCAPlot is at the end of this file

############## kmeans clustering of 500 most variable genes, figure SF2B ##################
#this requires the data loaded and processed for figure SF2A, above.
nGenes<-500
data2Cluster<-mL2[1:nGenes,]
res<-PCA(t(mL2[1:nGenes,])) #calculate PCA using factominer
res.km <- kmeans(scale(t(data2Cluster)), 3, nstart = 100)
ind.coord <- as.data.frame(get_pca_ind(res)$coord)
ind.coord$cluster <- factor(res.km$cluster)
# eigenvalue <- round(get_eigenvalue(res), 1)
# eigenvalue<-as.data.frame(eigenvalue)
# variance.percent <- eigenvalue$variance.percent
# 
# type<-factor(metaD$Histology.group, levels = c ("FS","PNST","PWT"))

cluster<-factor(res.km$cluster)

p<-ggplot(data = ind.coord) +
  geom_point(aes(x = Dim.1, y = Dim.2, color = cluster), size = 5)
p



MyPCAPlot<-function(mDat, computeLogs = TRUE, text = FALSE, ngenes = 2000, sTypes){
  #mDat - data to plot, if computeLogs is false is is expected that
  #    you've already computed them, otherwise the function uses
  #    edger's cpm function to normalize the data and computes log2
  #computeLogs - see above
  #text - label the points or not
  #ngenes - number of most variable genes to use for PCA calculations
  #sTypes - factors for point colors
  if(computeLogs){
    mCPM<-cpm(mDat)
    mL2<-log2(mCPM+.01)
    mVar<-apply(mL2,1,var)
    mL2<-mL2[order(mVar, decreasing = TRUE),]
  }
  else {
    mL2<-mDat
  }
  
  
  topVGenes<-mL2[1:ngenes,]
  
  res<-PCA(t(topVGenes)) #calculate PCA using factominer
  mVD<-res$eig[,1]/sum(res$eig[,1])
  barplot(mVD*100)
  print(sum(mVD[1:5])) #83% of variability in first 4 PCA components
  
  #Plot some PCA plots with various sample data
  
  pDat<-as.data.frame(res$ind$coord)
  pDat$name<-colnames(topVGenes)
  
  tTypes<-as.factor(sTypes)
  
  p<-ggplot(data = pDat)
  p<-p+geom_point(aes(x = Dim.1, y = Dim.2, color = tTypes), size = 5)
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
      panel.border = element_rect(colour = "black", fill=NA, size=2)
    )
  if (text)
    p <- p + geom_text_repel(aes(x = Dim.1, y = Dim.2, label = colnames(mDat)), force = 20, max.overlaps = 20)
  print(p)
  
  p<-ggplot(data = pDat)
  p <- p+geom_point(aes(x = Dim.3, y = Dim.4, color = tTypes), size = 5)
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
      panel.border = element_rect(colour = "black", fill=NA, size=2)
    )
  if (text)
    p <- p + geom_text_repel(aes(x = Dim.3, y = Dim.4, label = colnames(mDat)), force = 20, max.overlaps = 20)
  print(p)
  
}

##### PCA of all samples #######################
########  Load the prepared data ##########################
baseDir<-"C:\\Users\\idaho\\Dropbox\\"
baseDir<-"C:\\Users\\Eric Shelden\\Dropbox\\"
load(str_c(baseDir,"Shelden lab\\Canine SARC analysis\\Input Data\\all_data.rda"))
mDat<-all_data
colnames(mDat)

#retain genes where cpm is greater than 10 for at least 2 tumor and 4 normal tissues
mDat<-mDat[rowSums(cpm(mDat[,1:15]) > 10) > 2,]
mDat<-mDat[rowSums(cpm(mDat[,16:55]) > 10) > 4,] #changed to 4 from 2 on 4-28-21
sTypes = c(rep("cSTS",16), rep("cVSMS", 9), rep("cSK",6), rep("cCT",18), rep("cNS", 15))
length(sTypes)
dim(mDat)

MyPCAPlot(mDat, text = TRUE, sTypes = sTypes) #function defined above in this file.
MyPCAPlot(mDat, text = FALSE, sTypes = sTypes)
