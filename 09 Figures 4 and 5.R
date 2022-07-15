library(stringr)
library(dplyr)
library(ggplot2)

#You'll need to get the TF ENRICHR data from the supplemental data tables and 
#produce a csv file from them for input below

## Figure 4 ##
setwd("C:\\Users\\idaho\\Dropbox\\Shelden lab\\Canine SARC analysis\\Results")

mDat<-read.csv(file = "C:\\Users\\idaho\\Dropbox\\Shelden lab\\Canine SARC analysis\\Results\\Normal Tissue Comparisons by subtype\\Eric\\TFVals.csv")

#work on data from upregulated genes
upGenes<-mDat[mDat$GeneD == "up",]
View(upGenes)
#get columns I need
upGenes<-upGenes[,c(2,4,8,9,11,15,16,18,22)]

#isolate each tumor type
upGenesFS<-upGenes[,c(1,2,3)]
upGenesPNST<-upGenes[,c(4,5,6)]
upGenesPWT<-upGenes[,c(7,8,9)]
#order by significance
upGenesFS<-upGenesFS[order(upGenesFS$Adjusted.P.value),]
upGenesPNST<-upGenesPNST[order(upGenesPNST$Adjusted.P.value),]
upGenesPWT<-upGenesPWT[order(upGenesPWT$Adjusted.P.value),]

#
View(upGenesFS)
#

#create ranking
upGenesFS$rank<-1:12
upGenesPNST$rank<-1:12
upGenesPWT$rank<-1:12

#set all implied gene effects to up as default
upGenesFS$Direction<-"up"
upGenesPNST$Direction<-"up"
upGenesPWT$Direction<-"up"

#switch the direction of effects as appropriate
mFirstTerms<-unlist(str_split("KO ,KD ,INHIBITION ,SIRNA, SHRNA", ","))
mSwitch<-ifelse(str_detect(upGenesFS$original.term, paste0(mFirstTerms,collapse = "|")) & str_detect(upGenesFS$original.term, "UP$"), T, F)
upGenesFS$Direction[mSwitch]<-"down"                
mFirstTerms<-"OE "
mSwitch<-ifelse(str_detect(upGenesFS$original.term, mFirstTerms) & str_detect(upGenesFS$original.term, "DOWN$|DN$"), T, F)
upGenesFS$Direction[mSwitch]<-"down"

mFirstTerms<-unlist(str_split("KO ,KD ,INHIBITION ,SIRNA, SHRNA", ","))
mSwitch<-ifelse(str_detect(upGenesPNST$original.term, paste0(mFirstTerms,collapse = "|")) & str_detect(upGenesPNST$original.term, "UP$"), T, F)
upGenesPNST$Direction[mSwitch]<-"down"                
mFirstTerms<-"OE "
mSwitch<-ifelse(str_detect(upGenesPNST$original.term, mFirstTerms) & str_detect(upGenesPNST$original.term, "DOWN$|DN$"), T, F)
upGenesPNST$Direction[mSwitch]<-"down"

mFirstTerms<-unlist(str_split("KO ,KD ,INHIBITION ,SIRNA, SHRNA", ","))
mSwitch<-ifelse(str_detect(upGenesPWT$original.term, paste0(mFirstTerms,collapse = "|")) & str_detect(upGenesPWT$original.term, "UP$"), T, F)
upGenesPWT$Direction[mSwitch]<-"down"                
mFirstTerms<-"OE "
mSwitch<-ifelse(str_detect(upGenesPWT$original.term, mFirstTerms) & str_detect(upGenesPWT$original.term, "DOWN$|DN$"), T, F)
upGenesPWT$Direction[mSwitch]<-"down"

#order by Term aphabetical order
upGenesFS<-upGenesFS[order(upGenesFS$Term),]
upGenesPNST<-upGenesPNST[order(upGenesPNST$Term.1),]
upGenesPWT<-upGenesPWT[order(upGenesPWT$original.term.2),]
#check order, should be TRUE
all.equal(upGenesFS$Term,upGenesPNST$Term.1)
all.equal(upGenesFS$Term,upGenesPWT$Term.2)
rankDF<-data.frame(term<-upGenesFS$Term, mean.rank<-(upGenesFS$rank+upGenesPNST$rank+upGenesPWT$rank)/3); colnames(rankDF)<-c("term","mean.rank")
rankDF<-rankDF[order(rankDF$mean.rank),]

#add type
upGenesFS$type<-"FS"
upGenesPNST$type<-"PNST"
upGenesPWT$type<-"PWT"

colnames(upGenesFS)<-c("term","adjp","original.term","rank","Regulation","type")
colnames(upGenesPNST)<-c("term","adjp","original.term","rank","Regulation","type")
colnames(upGenesPWT)<-c("term","adjp","original.term","rank","Regulation","type")

#plot the data as bars.
upGenesFS$pAdj <- -log10(upGenesFS$adjp)
upGenesPNST$pAdj <- -log10(upGenesPNST$adjp)
upGenesPWT$pAdj <- -log10(upGenesPWT$adjp)

upGenesFS<-upGenesFS[order(upGenesFS$pAdj, decreasing = F),]
upGenesPNST<-upGenesPNST[order(upGenesPNST$pAdj, decreasing = F),]
upGenesPWT<-upGenesPWT[order(upGenesPWT$pAdj, decreasing = F),]

PlotBars<-function(pDat){
  pDat$term<-factor(pDat$term, levels = pDat$term)
  p<-ggplot(pDat, aes(x=term, y=pAdj, fill = Regulation)) +
    geom_bar(stat='identity') +
    scale_fill_manual(values = c("orange","blue")) +
    theme_classic()+
    theme (
      axis.text.y = element_text(size=16,face="bold", colour = "black"),
      axis.title.y = element_blank(),
      axis.text.x = element_text(size=16,face="bold", color = "black"),
      axis.title.x = element_blank(),
      axis.line = element_blank(),
      axis.ticks = element_blank(),
      legend.text = element_text(face = "bold", size = 16),
      legend.title = element_text(size=16, face="bold"),
      panel.border = element_blank()
    ) +
    ylim(0,maxPAdj) +
    coord_flip()
  print(p)
}

#work on data from downregulated genes
myGenes<-mDat[mDat$GeneD == "down",]
#get columns I need
myGenes<-myGenes[,c(2,4,8,9,11,15,16,18,22)]

#isolate each tumor type
myGenesFS<-myGenes[,c(1,2,3)]
myGenesPNST<-myGenes[,c(4,5,6)]
myGenesPWT<-myGenes[,c(7,8,9)]

#order by significance
myGenesFS<-myGenesFS[order(myGenesFS$Adjusted.P.value),]
myGenesPNST<-myGenesPNST[order(myGenesPNST$Adjusted.P.value.1),]
myGenesPWT<-myGenesPWT[order(myGenesPWT$Adjusted.P.value.2),]
#create ranking
myGenesFS$rank<-1:13
myGenesPNST$rank<-1:13
myGenesPWT$rank<-1:13
#set all implied gene effects to down as default
myGenesFS$Direction<-"down"
myGenesPNST$Direction<-"down"
myGenesPWT$Direction<-"down"

#switch the direction of effects as appropriate
mFirstTerms<-unlist(str_split("KO ,KD ,INHIBITION ,SIRNA, SHRNA, DEPLETION", ","))
mSwitch<-ifelse(str_detect(myGenesFS$original.term, paste0(mFirstTerms,collapse = "|")) & str_detect(myGenesFS$original.term, "UP$"), T, F)
myGenesFS$Direction[mSwitch]<-"up"                
mFirstTerms<-" OE "
mSwitch<-ifelse(str_detect(myGenesFS$original.term, mFirstTerms) & str_detect(myGenesFS$original.term, "DOWN$|DN$"), T, F)
myGenesFS$Direction[mSwitch]<-"up"

mFirstTerms<-unlist(str_split("KO ,KD ,INHIBITION ,SIRNA, SHRNA, DEPLETION", ","))
mSwitch<-ifelse(str_detect(myGenesPNST$original.term, paste0(mFirstTerms,collapse = "|")) & str_detect(myGenesPNST$original.term, "UP$"), T, F)
myGenesPNST$Direction[mSwitch]<-"up"                
mFirstTerms<-" OE "
mSwitch<-ifelse(str_detect(myGenesPNST$original.term, mFirstTerms) & str_detect(myGenesPNST$original.term, "DOWN$|DN$"), T, F)
myGenesPNST$Direction[mSwitch]<-"up"

mFirstTerms<-unlist(str_split("KO ,KD ,INHIBITION ,SIRNA, SHRNA, DEPLETION", ","))
mSwitch<-ifelse(str_detect(myGenesPWT$original.term, paste0(mFirstTerms,collapse = "|")) & str_detect(myGenesPWT$original.term, "UP$"), T, F)
myGenesPWT$Direction[mSwitch]<-"up"                
mFirstTerms<-" OE "
mSwitch<-ifelse(str_detect(myGenesPWT$original.term, mFirstTerms) & str_detect(myGenesPWT$original.term, "DOWN$|DN$"), T, F)
myGenesPWT$Direction[mSwitch]<-"up"

#order by Term aphabetical order
myGenesFS<-myGenesFS[order(myGenesFS$Term, decreasing = F),]
myGenesPNST<-myGenesPNST[order(myGenesPNST$Term.1, decreasing = F),]
myGenesPWT<-myGenesPWT[order(myGenesPWT$original.term.2, decreasing = F),]
#check order, should be TRUE
all.equal(myGenesFS$Term,myGenesPNST$Term.1)
all.equal(myGenesFS$Term,myGenesPWT$Term.2)

#add type
myGenesFS$type<-"FS"
myGenesPNST$type<-"PNST"
myGenesPWT$type<-"PWT"

colnames(myGenesFS)<-c("term","adjp","original.term","rank","Regulation","type")
colnames(myGenesPNST)<-c("term","adjp","original.term","rank","Regulation","type")
colnames(myGenesPWT)<-c("term","adjp","original.term","rank","Regulation","type")

myGenesFS$pAdj <- -log10(myGenesFS$adjp)
myGenesPNST$pAdj <- -log10(myGenesPNST$adjp)
myGenesPWT$pAdj <- -log10(myGenesPWT$adjp)

myGenesFS<-myGenesFS[order(myGenesFS$pAdj, decreasing = F),]
myGenesPNST<-myGenesPNST[order(myGenesPNST$pAdj, decreasing = F),]
myGenesPWT<-myGenesPWT[order(myGenesPWT$pAdj, decreasing = F),]

maxPAdj<-max(c(upGenesFS$pAdj, upGenesPNST$pAdj, upGenesPWT$pAdj))
maxPAdj<-max(c(myGenesFS$pAdj, myGenesPNST$pAdj, myGenesPWT$pAdj, maxPAdj))


PlotBars(upGenesFS)
PlotBars(upGenesPNST)
PlotBars(upGenesPWT)
PlotBars(myGenesFS)
PlotBars(myGenesPNST)
PlotBars(myGenesPWT)

# 
# ###################### dotplots Figure 5 ####################
# 
HandleUpGenes<-function(terms){
  res<-rep("up",length(terms))
  terms<-toupper(terms)
  mFirstTerms<-unlist(str_split("KO ,KD ,INHIBITION ,SIRNA, SHRNA, DEPLETION, LOF, INACTIVATION, DEFICIENCY", ","))
  mSwitch<-ifelse(str_detect(terms, paste0(mFirstTerms,collapse = "|")) & str_detect(terms, "UP$"), T, F)
  res[mSwitch]<-"down"
  mFirstTerms<-" OE | UNK | ACTIVATION | GOF | KNOCKIN"
  mSwitch<-ifelse(str_detect(terms, mFirstTerms) & str_detect(terms, "DOWN$|DN$"), T, F)
  res[mSwitch]<-"down"
  return(res)
}

HandleDownGenes<-function(terms){
  res<-rep("down",length(terms))
  terms<-toupper(terms)
  mFirstTerms<-unlist(str_split("KO ,KD ,INHIBITION ,SIRNA, SHRNA, DEPLETION, LOF, INACTIVATION, DEFICIENCY", ","))
  mSwitch<-ifelse(str_detect(terms, paste0(mFirstTerms,collapse = "|")) & str_detect(terms, "UP$"), T, F)
  res[mSwitch]<-"up"
  mFirstTerms<-" OE | UNK | ACTIVATION | GOF | KNOCKIN"
  mSwitch<-ifelse(str_detect(terms, mFirstTerms) & str_detect(terms, "DOWN$|DN$"), T, F)
  res[mSwitch]<-"up"
  return(res)
}

##############################################
setwd("C:\\Users\\idaho\\Dropbox\\Shelden lab\\Canine SARC analysis\\Results\\Normal Tissue Comparisons by subtype\\Eric")
mDir<-getwd()

mFun<-function(x){
  parts<-unlist(str_split(x," "))
  return(parts[1])
}

getSize<-function(x){
  parts<-unlist(str_split(x,"/"))
  return(as.numeric(parts[1]))
}

#data files produced by ENRICHR scripts

load(file = "FS_res_all_from_up.rda")
load(file = "PNST_res_all_from_up.rda")
load(file = "PWT_res_all_from_up.rda")

FS.up.res<-FS.res.all.from.up
PNST.up.res<-PNST.res.all.from.up
PWT.up.res<-PWT.res.all.from.up

load(file = "FS_res_all_from_dn.rda")
load(file = "PNST_res_all_from_dn.rda")
load(file = "PWT_res_all_from_dn.rda")

FS.dn.res<-FS.res.all.from.dn
PNST.dn.res<-PNST.res.all.from.dn
PWT.dn.res<-PWT.res.all.from.dn

FSTFs<-FS.up.res$TFs
#View(FSTFs)
FSTFsUp<-FSTFs[,c(11,10,4)]
FSTFsUp$nGenes<-unlist(lapply(FSTFs$Overlap,getSize))
head(FSTFsUp)
FSTFsUp$Regulation<-HandleUpGenes(FSTFsUp$original.term)
rownames(FSTFsUp)<-NULL

PNSTTFs<-PNST.up.res$TFs
#View(PNSTTPNST)
PNSTTFsUp<-PNSTTFs[,c(11,10,4)]
PNSTTFsUp$nGenes<-unlist(lapply(PNSTTFs$Overlap,getSize))
head(PNSTTFsUp)
PNSTTFsUp$Regulation<-HandleUpGenes(PNSTTFs$original.term)
rownames(PNSTTFsUp)<-NULL

PWTTFs<-PWT.up.res$TFs
#View(PWTTPWT)
PWTTFsUp<-PWTTFs[,c(11,10,4)]
PWTTFsUp$nGenes<-unlist(lapply(PWTTFs$Overlap,getSize))
head(PWTTFsUp)
PWTTFsUp$Regulation<-HandleUpGenes(PWTTFs$original.term)
rownames(PWTTFsUp)<-NULL

#### now down ##

FSTFs<-FS.dn.res$TFs
#View(FSTFs)
FSTFsDn<-FSTFs[,c(11,10,4)]
FSTFsDn$nGenes<-unlist(lapply(FSTFs$Overlap,getSize))
head(FSTFsDn)
FSTFsDn$Regulation<-HandleDownGenes(FSTFsDn$original.term)
rownames(FSTFsDn)<-NULL

PNSTTFs<-PNST.dn.res$TFs
#View(PNSTTPNST)
PNSTTFsDn<-PNSTTFs[,c(11,10,4)]
PNSTTFsDn$nGenes<-unlist(lapply(PNSTTFs$Overlap,getSize))
head(PNSTTFsDn)
PNSTTFsDn$Regulation<-HandleDownGenes(PNSTTFs$original.term)
rownames(PNSTTFsDn)<-NULL

PWTTFs<-PWT.dn.res$TFs
#View(PWTTPWT)
PWTTFsDn<-PWTTFs[,c(11,10,4)]
PWTTFsDn$nGenes<-unlist(lapply(PWTTFs$Overlap,getSize))
head(PWTTFsDn)
PWTTFsDn$Regulation<-HandleDownGenes(PWTTFs$original.term)
rownames(PWTTFsDn)<-NULL

FS<-rbind(FSTFsUp, FSTFsDn)
PNST<-rbind(PNSTTFsUp, PNSTTFsDn)
PWT<-rbind(PWTTFsUp, PWTTFsDn)

FS<-FS[order(FS$Adjusted.P.value),]
PNST<-PNST[order(PNST$Adjusted.P.value),]
PWT<-PWT[order(PWT$Adjusted.P.value),]
head(FS)
head(PNST)
head(PWT)

nTerms<-75
FSUnique<-setdiff(FS$raw.term[1:nTerms],union(PNST$raw.term[1:nTerms],PWT$raw.term[1:nTerms]))
FSUnique

nTerms<-50
PNSTUnique<-setdiff(PNST$raw.term[1:nTerms],union(FS$raw.term[1:nTerms],PWT$raw.term[1:nTerms]))
PNSTUnique


nTerms<-15
PWTUnique<-setdiff(PWT$raw.term[1:nTerms],union(FS$raw.term[1:nTerms],PNST$raw.term[1:nTerms]))
PWTUnique

pDat<-FS[match(FSUnique,FS$raw.term),]
PlotDots(pDat)
pDat<-PNST[match(PNSTUnique,PNST$raw.term),]
PlotDots(pDat)
pDat<-PWT[match(PWTUnique,PWT$raw.term),]
PlotDots(pDat)

PlotDots<-function(pDat){
  colnames(pDat)<-c("term","original.term","padj","ngenes","regulation")
  pDat[pDat$ngenes>200,]$ngenes<-195
  pDat$term<-factor(pDat$term, levels = rev(pDat$term))
  p <- ggplot(pDat, aes(x = -log10(padj), y = (term), size = ngenes, color = regulation)) + 
    geom_point() +
    #scale_size_area(max_size = 5) +
    scale_size(breaks = c(50,100,150,200), limits = c(1,200), range = c(2,10)) +
    scale_colour_manual(values = c("orange","blue")) +
    labs(x = "-log10(padj)", y = element_blank(),color = "regulation", size = "Genes") +
    theme_bw() +
    scale_x_continuous(breaks = c(65,70,75,80,85)) +
    theme (
      axis.text.y = element_text(size=16,face="bold", colour = "black"),
      axis.title.y = element_text(size=16,face="bold", colour = "black"),
      axis.text.x = element_text(size=16,face="bold", color = "black"),
      axis.title.x = element_text(size=16,face="bold", color = "black"),
      axis.line = element_blank(),
      axis.ticks = element_blank(),
      legend.text = element_text(face = "bold", size = 16),
      legend.title = element_text(size=16, face="bold")
      #panel.border = element_blank()
    )
  print(p)
}
