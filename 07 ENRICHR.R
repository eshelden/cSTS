library(readr)
library(enrichR)
library(stringr)
library(lsa)
library(dendextend)
library(SummarizedExperiment)
library(DESeq2)
library(gplots)
library(dendextend)
library(ggplot2)

#require("13 Functional Analysis Utils 5.R")


#Figure 3A
setwd("C:\\Users\\idaho\\Dropbox\\Shelden lab\\Canine SARC analysis\\Results\\Tumor type specific gene expression\\Eric")
setwd("C:\\Users\\Eric Shelden\\Dropbox\\Shelden lab\\Canine SARC analysis\\Results\\Tumor type specific gene expression\\Eric")
mDir<-getwd()
mDB<-c("ARCHS4_Tissues","ARCHS4_Cell-lines")

PNST.F<-read.csv(paste0(mDir,"\\PNST_specific_genes_final.csv"))
FS.F<-read.csv(paste0(mDir,"\\FS_specific_genes_final.csv"))
PWT.F<-read.csv(paste0(mDir,"\\PWT_specific_genes_final.csv"))
dim(PNST.F) #291
dim(FS.F) #2090
dim(PWT.F) #875

names<-FS.F$X
names<-str_remove(names,"gene-")
names<-names[!grepl("^LOC",names)]
length(names)#1637
FS.res<-enrichr(names,mDB)
FS.res <- ProcessEnrichrRes(FS.res,FS.F, nRes = 25) #formats the data for inclusion in supplemental data table
write.csv(FS.res, file = "FS_tissue_and_cell_enrichments.csv")
temp<-(FS.res[1:15,])
write.clip(temp) #export to paste into supplemental data table

names<-PNST.F$X
names<-str_remove(names,"gene-")
names<-names[!grepl("^LOC",names)]
length(names)#210
PNST.res<-enrichr(names,mDB)
PNST.res <- ProcessEnrichrRes(PNST.res,PNST.F, nRes = 25)
write.csv(PNST.res, file = "PNST_tissue_and_cell_enrichments.csv")
temp<-cbind(temp,PNST.res[1:15,])
write.clip(temp)

names<-PWT.F$X
names<-str_remove(names,"gene-")
names<-names[!grepl("^LOC",names)]
length(names)#540
PWT.res<-enrichr(names,mDB)
PWT.res <- ProcessEnrichrRes(PWT.res,PWT.F, nRes = 25)
write.csv(PWT.res, file = "PWT_tissue_and_cell_enrichments.csv")
temp<-cbind(temp,PWT.res[1:15,])
write.clip(temp)

#this is where the plotting for Figure 3A gets set up and done.
#just the top 15 most significant enrichments
FS.res<-FS.res[1:15,]
PNST.res<-PNST.res[1:15,]
PWT.res<-PWT.res[1:15,]

#compute negative log 10 of padj
FS.res$padj<- -log10(FS.res$Adjusted.P.value)
PNST.res$padj<- -log10(PNST.res$Adjusted.P.value)
PWT.res$padj<- -log10(PWT.res$Adjusted.P.value)

#find values needed to get all the plots on the same axis values for significance and number of genes
maxPadj<-max(c(FS.res$padj, PNST.res$padj,PWT.res$padj))
temp<-c(FS.res[,2],PNST.res[,2],PWT.res[,2])
temp<-unlist(lapply(temp,getSize))
maxCount<-max(as.numeric(temp))

#round scales to nearest 5
maxPadj<-ceiling(maxPadj/5)*5
maxCount<-ceiling(maxCount/5)*5

pDat<-data.frame(name = factor(FS.res[,1], levels = rev(FS.res[,1])), padj = FS.res$padj, count = as.numeric(lapply(FS.res[,2], getSize)))
PrintDotPlot(pDat, colorLimits = c(1.3,maxPadj), sizeLimits = c(1,maxCount), maxSize = 12)

pDat<-data.frame(name = factor(PNST.res[,1], levels = rev(PNST.res[,1])), padj = PNST.res$padj, count = as.numeric(lapply(PNST.res[,2], getSize)))
PrintDotPlot(pDat, colorLimits = c(1.3,maxPadj), sizeLimits = c(1,maxCount), maxSize = 12)

pDat<-data.frame(name = factor(PWT.res[,1], levels = rev(PWT.res[,1])), padj = PWT.res$padj, count = as.numeric(lapply(PWT.res[,2], getSize)))
PrintDotPlot(pDat, colorLimits = c(1.3,maxPadj), sizeLimits = c(1,maxCount), maxSize = 12)

####################### Data for Figures 4 and 5 #############################################
setwd("C:\\Users\\idaho\\Dropbox\\Shelden lab\\Canine SARC analysis\\Results\\Normal Tissue Comparisons by subtype\\Eric")
setwd("C:\\Users\\Eric Shelden\\Dropbox\\Shelden lab\\Canine SARC analysis\\Results\\Normal Tissue Comparisons by subtype\\Eric")
mDir<-getwd()

#load expression data
PC_DEGs_FS_tumors_vs_cSKIN <- read_csv("PC_DEGs_FS_tumors_vs_cSKIN.csv")
dim(PC_DEGs_FS_tumors_vs_cSKIN) #3802

PC_DEGs_PNST_tumors_vs_cSKIN <- read_csv("PC_DEGs_PNST_tumors_vs_cSKIN.csv")
dim(PC_DEGs_PNST_tumors_vs_cSKIN) #3747

PC_DEGs_PWT_tumors_vs_cVSMS <- read_csv("PC_DEGs_PWT_tumors_vs_cVSMS.csv")
dim(PC_DEGs_PWT_tumors_vs_cVSMS) #4442

colnames(PC_DEGs_FS_tumors_vs_cSKIN)[1]<-"name"
colnames(PC_DEGs_PNST_tumors_vs_cSKIN)[1]<-"name"
colnames(PC_DEGs_PWT_tumors_vs_cVSMS)[1]<-"name"

#find upregulated genes
FS_up<-PC_DEGs_FS_tumors_vs_cSKIN[PC_DEGs_FS_tumors_vs_cSKIN$log2FoldChange > 0,]
dim(FS_up)  #2144
2144/3802*100 #56.4% is the fraction of upregulated genes in total genes

PWT_up<-PC_DEGs_PWT_tumors_vs_cVSMS[PC_DEGs_PWT_tumors_vs_cVSMS$log2FoldChange > 0,]
dim(PWT_up) #2528
2528/4442*100 #56.9%

PNST_up<-PC_DEGs_PNST_tumors_vs_cSKIN[PC_DEGs_PNST_tumors_vs_cSKIN$log2FoldChange > 0,]
dim(PNST_up) #1903
1903/3747*100 #50.78%

#find down regulated genes
FS_dn<-PC_DEGs_FS_tumors_vs_cSKIN[PC_DEGs_FS_tumors_vs_cSKIN$log2FoldChange < 0,]
dim(FS_dn) #1658

PWT_dn<-PC_DEGs_PWT_tumors_vs_cVSMS[PC_DEGs_PWT_tumors_vs_cVSMS$log2FoldChange < 0,]
dim(PWT_dn) #1914

PNST_dn<-PC_DEGs_PNST_tumors_vs_cSKIN[PC_DEGs_PNST_tumors_vs_cSKIN$log2FoldChange < 0,]
dim(PNST_dn) #1844

save(FS_up, file = "FS_UP_DEGS.rda")
save(PNST_up, file = "PNST_UP_DEGS.rda")
save(PWT_up, file = "PWT_UP_DEGS.rda")
save(FS_dn, file = "FS_DN_DEGS.rda")
save(PNST_dn, file = "PNST_DN_DEGS.rda")
save(PWT_dn, file = "PWT_DN_DEGS.rda")


### Generate data for transcription factor enrichments for genes upregulated in comparison #####
### to normal tissues counterparts for all three tumor types.                              #####
##  This is the data for Figure 4                                                          #####
################################################################################################

setwd("C:\\Users\\idaho\\Dropbox\\Shelden lab\\Canine SARC analysis\\Results\\Normal Tissue Comparisons by subtype\\Eric")
setwd("C:\\Users\\Eric Shelden\\Dropbox\\Shelden lab\\Canine SARC analysis\\Results\\Normal Tissue Comparisons by subtype\\Eric")
mDir<-getwd()

#load files below on line 195 if you've already done this analysis

nRes<-200

#This is an analysis of the upregulated genes....
names<-FS_up$name
names<-str_remove(names,"gene-")
names<-names[!grepl("^LOC",names)]
length(names)#2144
FS.res.all.from.up<-Get.ENRICHR.Results(names,nRes, doPlot = F, minSig = .1, type = "FS", fromFile = F)
save(FS.res.all.from.up, file = "FS_res_all_from_up.rda")

names<-PNST_up$name
names<-str_remove(names,"gene-")
names<-names[!grepl("^LOC",names)]
length(names)#1903
PNST.res.all.from.up<-Get.ENRICHR.Results(names,nRes, doPlot = F, minSig = .1, type = "PNST", fromFile = F)
save(PNST.res.all.from.up, file = "PNST_res_all_from_up.rda")

names<-PWT_up$name
names<-str_remove(names,"gene-")
names<-names[!grepl("^LOC",names)]
length(names)#2528
PWT.res.all.from.up<-Get.ENRICHR.Results(names,nRes, doPlot = F, minSig = .1, type = "PWT", fromFile = F)
save(PWT.res.all.from.up, file = "PWT_res_all_from_up.rda")

#this is an analysis of the down regulated genes
names<-FS_dn$name
names<-str_remove(names,"gene-")
names<-names[!grepl("^LOC",names)]
length(names)#1658
FS.res.all.from.dn<-Get.ENRICHR.Results(names,nRes, doPlot = F, minSig = .1, type = "FS", fromFile = F)
save(FS.res.all.from.dn, file = "FS_res_all_from_dn.rda")

names<-PNST_dn$name
names<-str_remove(names,"gene-")
names<-names[!grepl("^LOC",names)]
length(names)#1844
PNST.res.all.from.dn<-Get.ENRICHR.Results(names,nRes, doPlot = F, minSig = .1, type = "PNST", fromFile = F)
save(PNST.res.all.from.dn, file = "PNST_res_all_from_dn.rda")

names<-PWT_dn$name
names<-str_remove(names,"gene-")
names<-names[!grepl("^LOC",names)]
length(names)#1914
PWT.res.all.from.dn<-Get.ENRICHR.Results(names,nRes, doPlot = F, minSig = .1, type = "PWT", fromFile = F)
save(PWT.res.all.from.dn, file = "PWT_res_all_from_dn.rda")

#if you've already completed the above, then load the data directly
load(file = "FS_res_all_from_up.rda")
load(file = "PNST_res_all_from_up.rda")
load(file = "PWT_res_all_from_up.rda")
load(file = "FS_res_all_from_dn.rda")
load(file = "PNST_res_all_from_dn.rda")
load(file = "PWT_res_all_from_dn.rda")

#################################################################################
#Create tables for transcription factors common to all tumor types
#paste the data into an excel spreadsheet
#Generate Figure 4 using FinalHeatMaps.R

transcriptionFactors<-1

mRes<-HandleCommonGenes(FS.res.all.from.up, PNST.res.all.from.up, PWT.res.all.from.up,transcriptionFactors,nRes = 12, nTerms = 50, trimNames = F, genesReturn = "both")
write.clip(mRes)

##############   analysis of DOWN-regulated genes ###########################################

mRes<-HandleCommonGenes(FS.res.all.from.dn, PNST.res.all.from.dn, PWT.res.all.from.dn,transcriptionFactors,nRes = 13, nTerms = 55, trimNames = F, genesReturn = "down")
write.clip(mRes)

################################ Figure 5 data ########################################
#plot the things that are unique to each tumor type.. 10-28-21
#need these for plots of things specific to individual tumor types
#This is Figure 5

setwd("C:\\Users\\idaho\\Dropbox\\Shelden lab\\Canine SARC analysis\\Results\\Normal Tissue Comparisons by subtype\\Eric")
setwd("C:\\Users\\Eric Shelden\\Dropbox\\Shelden lab\\Canine SARC analysis\\Results\\Normal Tissue Comparisons by subtype\\Eric")
mDir<-getwd()

#skip below to load data if already processed

load(file = "FS_UP_DEGS.rda")
load(file = "PNST_UP_DEGS.rda")
load(file = "PWT_UP_DEGS.rda")

#Query enrichr with names of genes that are upregulated and unique to each
#tumor type

names<- setdiff(FS_up$name,union(PNST_up$name,PWT_up$name))
names<-str_remove(names,"gene-")
names<-names[!grepl("^LOC",names)]
length(names) #622
FS.res<-Get.ENRICHR.Results(t=names,nRes = 200, doPlot = F,minSig = .1, fromFile = F)

names<- setdiff(PNST_up$name,union(FS_up$name,PWT_up$name))
names<-str_remove(names,"gene-")
names<-names[!grepl("^LOC",names)]
length(names) #260
PNST.res<-Get.ENRICHR.Results(t=names,nRes = 200, doPlot = F,minSig = .1, fromFile = F)

names<- setdiff(PWT_up$name,union(FS_up$name,PNST_up$name))
names<-str_remove(names,"gene-")
names<-names[!grepl("^LOC",names)]
length(names) #1472
PWT.res<-Get.ENRICHR.Results(t=names,nRes = 200, doPlot = F,minSig = .1, fromFile = F)

save(FS.res, file = "FS_specific_genes_enrichr_results.rda")
save(PNST.res, file = "PNST_specific_genes_enrichr_results.rda")
save(PWT.res, file = "PWT_specific_genes_enrichr_results.rda")

load(file = "FS_specific_genes_enrichr_results.rda")
load(file = "PNST_specific_genes_enrichr_results.rda")
load(file = "PWT_specific_genes_enrichr_results.rda")

#load expression level data
PC_DEGs_FS_tumors_vs_cSKIN <- read_csv("PC_DEGs_FS_tumors_vs_cSKIN.csv")
dim(PC_DEGs_FS_tumors_vs_cSKIN) #3802
PC_DEGs_PNST_tumors_vs_cSKIN <- read_csv("PC_DEGs_PNST_tumors_vs_cSKIN.csv")
dim(PC_DEGs_PNST_tumors_vs_cSKIN) #3747
PC_DEGs_PWT_tumors_vs_cVSMS <- read_csv("PC_DEGs_PWT_tumors_vs_cVSMS.csv")
dim(PC_DEGs_PWT_tumors_vs_cVSMS)#4442
colnames(PC_DEGs_FS_tumors_vs_cSKIN)[1]<-"name"
colnames(PC_DEGs_PNST_tumors_vs_cSKIN)[1]<-"name"
colnames(PC_DEGs_PWT_tumors_vs_cVSMS)[1]<-"name"

transcriptionFactors<-1

#tumor drivers
tRes<-HandleDotPlots(index = transcriptionFactors,type = "drivers",nRes = 15) #formats data, creates plots
FinishDotPlotRes(tRes) #(function in Functional Analysis Utils, just writes the table to the clipboard



########################  down regulated genes specific to individual tumors #####################
load(file = "FS_DN_DEGS.rda")
load(file = "PNST_DN_DEGS.rda")
load(file = "PWT_DN_DEGS.rda")

names<- setdiff(FS_dn$name,union(PNST_dn$name,PWT_dn$name))
names<-str_remove(names,"gene-")
names<-names[!grepl("^LOC",names)]
length(names) #466
FS.res.dn<-Get.ENRICHR.Results(t=names,nRes = 200, doPlot = F,minSig = .1, fromFile = F)

names<- setdiff(PNST_dn$name,union(FS_dn$name,PWT_dn$name))
names<-str_remove(names,"gene-")
names<-names[!grepl("^LOC",names)]
length(names) #448
PNST.res.dn<-Get.ENRICHR.Results(t=names,nRes = 200, doPlot = F,minSig = .1, fromFile = F)

names<- setdiff(PWT_dn$name,union(FS_dn$name,PNST_dn$name))
names<-str_remove(names,"gene-")
names<-names[!grepl("^LOC",names)]
length(names) #1294
PWT.res.dn<-Get.ENRICHR.Results(t=names,nRes = 200, doPlot = F,minSig = .1, fromFile = F)

save(FS.res.dn, file = "FS_specific_dn_genes_enrichr_results.rda")
save(PNST.res.dn, file = "PNST_specific_dn_genes_enrichr_results.rda")
save(PWT.res.dn, file = "PWT_specific_dn_genes_enrichr_results.rda")

load(file = "FS_specific_dn_genes_enrichr_results.rda")
load(file = "PNST_specific_dn_genes_enrichr_results.rda")
load(file = "PWT_specific_dn_genes_enrichr_results.rda")

#HandleDotPlots assumes the three data sets are called FS.res, PNST.res and PWT.res
FS.res<-FS.res.dn
PNST.res<-PNST.res.dn
PWT.res<-PWT.res.dn

transcriptionFactors<-1

#tumor drivers
tRes<-HandleDotPlots(index = transcriptionFactors,type = "drivers",nRes = 15)
FinishDotPlotRes(tRes)

sessionInfo()
