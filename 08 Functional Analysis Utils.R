# #write a csv file with the data that makes up cDat...
write.clip<-function(x,...){
  write.table(x,"clipboard-256",sep="\t",row.names=F,col.names=colnames(x),...)
}
# 
# mFun<-function(x){
#   parts<-unlist(str_split(x," "))
#   return(parts[1])
# }
# 
getSize<-function(x){
  parts<-unlist(str_split(x,"/"))
  return(as.numeric(parts[1]))
}
# 
#
#this takes data from ENRICHR and formats the data and some derivative values
#into a dataframe for export into an excel file
ProcessEnrichrRes<-function(myRes,expData, nRes){
  myRes<-do.call("rbind", myRes)
  myRes<-as.data.frame(myRes)
  myRes<-myRes[order(myRes[,'Adjusted.P.value'], decreasing = F),]
  myRes<-myRes[1:nRes,]
  myRes$padj<- -log10(myRes[,'Adjusted.P.value'])
  myRes$matched.genes<-unlist(lapply(myRes$Overlap,str_remove,pattern = "/.*$"))
  myRes$total.genes<-unlist(lapply(myRes$Overlap,str_remove,pattern = "^.*/"))
  myRes$nUpGenes <- 0
  myRes$upGenes <- ""
  myRes$nDnGenes <- 0
  myRes$dnGenes <- ""

  expData<-expData[order(expData$log2FoldChange, decreasing = T),] #sort expression data by l2fChange

  for (i in 1:nrow(myRes)){
    genes<-myRes[i,]$Genes
    genes<-unlist(str_split(genes, ";"))
    genes<-paste0("gene-",genes)

    myExpData<-expData[expData$X %in% genes,]
    upGenes<-myExpData[myExpData$log2FoldChange > 0,]$X
    upGenes<-str_remove(upGenes,"gene-")
    myRes[i,]$nUpGenes<-length(upGenes)
    myRes[i,]$upGenes<-paste(upGenes,collapse = ";")
    dnGenes<-myExpData[myExpData$log2FoldChange < 0,]$X
    dnGenes<-str_remove(dnGenes,"gene-")
    dnGenes<-rev(dnGenes)#want the most highly downregulated genes first...
    myRes[i,]$nDnGenes<-length(dnGenes)
    myRes[i,]$dnGenes<-paste(dnGenes,collapse = ";")
  }
  myRes$db<-str_remove(rownames(myRes), "\\..*$")
  return(myRes[,c('Term','matched.genes','total.genes','P.value',
                  'Adjusted.P.value','padj','nUpGenes','upGenes','nDnGenes','dnGenes','db')])
}

PrintDotPlot<-function(df,sizeLimits = c(1,100), colorLimits = c(1,100), maxSize = 10){
  colnames(df)<-c("col","padj", "count")
  # create graph
  p <- ggplot(df)+
    geom_point(aes(x = padj, y = col, size = count, color = padj)) +
    scale_color_gradient(low = "steelblue", high = "red", limits = colorLimits) +
    xlab("-log10*padj") +
    ylab("") +
    scale_size_area(limits = sizeLimits, max_size = maxSize)

  p<-p +
    theme_bw() +
    theme (
      plot.title = element_text(size=16,face="bold", hjust = 0.5),
      axis.text.y = element_text(size=16,face="bold"),
      axis.title.y = element_text(size=16,face="bold"),
      axis.text.x = element_text(size=16,face="bold"),
      axis.title.x = element_text(size=16,face="bold"),
      axis.line = element_line(size = 1.5),
      axis.ticks = element_line(size = 1.5),
      legend.text = element_text(face = "bold", size = 16),
      legend.title = element_text(size=16, face="bold")
    ) +
    theme(plot.margin=unit(c(1,1,1,1),"cm"))
  print(p)
}

# #formats data for export to supplemental data tables
ProcessHeatMapRes<-function(myRes,expData, nRes, genesReturn = "both"){
  myRes<-as.data.frame(myRes)
  myRes<-myRes[order(myRes[,'Adjusted.P.value'], decreasing = F),]
  myRes<-myRes[1:nRes,]
  myRes$matched.genes<-unlist(lapply(myRes$Overlap,str_remove,pattern = "/.*$"))
  myRes$total.genes<-unlist(lapply(myRes$Overlap,str_remove,pattern = "^.*/"))
  myRes$nUpGenes <- 0
  myRes$upGenes <- ""
  myRes$nDnGenes <- 0
  myRes$dnGenes <- ""

  expData<-expData[order(expData$log2FoldChange, decreasing = T),]

  for (i in 1:nrow(myRes)){
    genes<-myRes[i,]$Genes
    genes<-unlist(str_split(genes, ";"))
    genes<-paste0("gene-",genes)

    myExpData<-expData[expData$name %in% genes,]
    upGenes<-myExpData[myExpData$log2FoldChange > 0,]$name
    upGenes<-str_remove(upGenes,"gene-")
    myRes[i,]$nUpGenes<-length(upGenes)
    myRes[i,]$upGenes<-paste(upGenes,collapse = ";")
    dnGenes<-myExpData[myExpData$log2FoldChange < 0,]$name
    dnGenes<-str_remove(dnGenes,"gene-")
    dnGenes<-rev(dnGenes)#want the most highly downregulated genes first...
    myRes[i,]$nDnGenes<-length(dnGenes)
    myRes[i,]$dnGenes<-paste(dnGenes,collapse = ";")
  }
  myRes$db<-str_remove(rownames(myRes), "\\..*$")
  if (genesReturn == "both"){
    return(myRes[,c('raw.term','matched.genes','total.genes','P.value',
                    'Adjusted.P.value','padj','nUpGenes','upGenes',
                    'nDnGenes','dnGenes','db','original.term')])
  }
  if (genesReturn == "up"){
    return(myRes[,c('raw.term','matched.genes','total.genes','P.value',
                    'Adjusted.P.value','padj','nUpGenes','upGenes',
                    'db','original.term')])
  }
  if (genesReturn == "down"){
    return(myRes[,c('raw.term','matched.genes','total.genes','P.value',
                    'Adjusted.P.value','padj','nDnGenes','dnGenes',
                    'db','original.term')])
  }
}

ProcessDotPlotRes<-function(myRes,expData, nRes){
  #myRes<-do.call("rbind", myRes)
  myRes<-as.data.frame(myRes)
  myRes<-myRes[order(myRes[,4], decreasing = F),]
  myRes<-myRes[1:nRes,]
  #myRes$padj<- -log10(myRes[,4]) already done
  myRes$matched.genes<-unlist(lapply(myRes$Overlap,str_remove,pattern = "/.*$"))
  myRes$total.genes<-unlist(lapply(myRes$Overlap,str_remove,pattern = "^.*/"))
  myRes$nUpGenes <- 0
  myRes$upGenes <- ""
  myRes$nDnGenes <- 0
  myRes$dnGenes <- ""

  expData<-expData[order(expData$log2FoldChange, decreasing = T),]

  for (i in 1:nrow(myRes)){
    genes<-myRes[i,]$Genes
    genes<-unlist(str_split(genes, ";"))
    genes<-paste0("gene-",genes)

    myExpData<-expData[expData$name %in% genes,]
    upGenes<-myExpData[myExpData$log2FoldChange > 0,]$name
    upGenes<-str_remove(upGenes,"gene-")
    myRes[i,]$nUpGenes<-length(upGenes)
    myRes[i,]$upGenes<-paste(upGenes,collapse = ";")
    dnGenes<-myExpData[myExpData$log2FoldChange < 0,]$name
    dnGenes<-str_remove(dnGenes,"gene-")
    dnGenes<-rev(dnGenes)#want the most highly downregulated genes first...
    myRes[i,]$nDnGenes<-length(dnGenes)
    myRes[i,]$dnGenes<-paste(dnGenes,collapse = ";")
  }
  myRes$db<-str_remove(rownames(myRes), "\\..*$")
  return(myRes)
}

PrintPlot<-function(df,limits = NA){
  colnames(df)<-c("col","val")
  # create graph
  if(length(limits) == 1) {
    p<-ggplot(data = df, aes(x = col, y = val)) +
      geom_bar(stat = 'identity', aes(fill = val)) +
      scale_fill_gradientn(colors = c("steelblue", "red")) +
      coord_flip() +
      ylab("-log10*padj") +
      xlab("") +
      labs(fill = "")
  }
  else {
    p<-ggplot(data = df, aes(x = col, y = val)) +
      geom_bar(stat = 'identity', aes(fill = val)) +
      scale_fill_gradientn(colors = c("steelblue", "red"), limits = limits) +
      coord_flip() +
      ylab("-log10*padj") +
      xlab("") +
      labs(fill = "")
  }


  p<-p +
    theme_bw() +
    theme (
      plot.title = element_text(size=16,face="bold", hjust = 0.5),
      axis.text.y = element_text(size=16,face="bold"),
      axis.title.y = element_text(size=16,face="bold"),
      axis.text.x = element_text(size=16,face="bold"),
      axis.title.x = element_text(size=16,face="bold"),
      axis.line = element_line(size = 1.5),
      axis.ticks = element_line(size = 1.5),
      legend.text = element_text(face = "bold", size = 16),
      legend.title = element_text(size=16, face="bold")
    ) +
    theme(plot.margin=unit(c(1,1,1,1),"cm"))
  print(p)
}

# GetTFs<-function(x){
#   parts<-unlist(str_split(x,"-"))
#   return(toupper(parts[1]))
# }
# 
# Get.Most.Significant.TFs<-function(mData){
#   mData$TFID<-unlist(lapply(mData$Term,GetTFs)) #removes stuff after the dash and makes it all upper case
#   mNames<-names(table(mData$TFID))
#   index<-match(mNames,mData$TFID)
#   return(mData[index,])
# }
# 
# ################### Transcription factor filter ################
Filter.TFs<-function(tData,lib = NA){
  terms<-tData[,1]
  if(!is.na(lib)){
    if(str_detect(lib,"ChEA_2016")){
      terms<-Filter.ChEA_2016(terms)
    }
    if(str_detect(lib,"ENCODE_and_ChEA_Consensus_TFs_from_ChIP")){
      terms<-Filter.ENCODE_and_ChEA_Consensus_TFs_from_ChIP(terms)
    }
    if(str_detect(lib,"ENCODE_TF_ChIP")){
      terms<-Filter.ENCODE_TF_ChIP(terms)
    }
    if(str_detect(lib,"TF_Perturbations_Followed_by_Expression")){
      terms<-Filter.TF_Perturbations_Followed_by_Expression(terms)
    }
    if(str_detect(lib,"LOF_Expression_from_GEO")){
      terms<-Filter.LOF_Expression_from_GEO(terms)
    }
    if(str_detect(lib,"Transcription_Factor_PPIs")){
      terms<-Filter.Transcription_Factor_PPIs(terms)
    }
    if(str_detect(lib,"TRANSFAC_and_JASPAR_PWMs")){
      terms<-Filter.TRANSFAC_and_JASPAR_PWMs(terms)
    }
    if(str_detect(lib,"TRRUST_Transcription_Factors_2019")){
      terms<-Filter.TRRUST_Transcription_Factors_2019(terms)
    }
  }
  tData[,1]<-terms
  return(tData)
}
# #################################################################################

Filter.ChEA_2016<-function(x){
  return(str_to_upper(unlist(lapply(x,str_remove,pattern=" .*$"))))
}
Filter.ENCODE_and_ChEA_Consensus_TFs_from_ChIP<-function(x){
  return(str_to_upper(unlist(lapply(x,str_remove,pattern=" .*$"))))
}
Filter.ENCODE_TF_ChIP<-function(x){
  return(str_to_upper(unlist(lapply(x,str_remove,pattern=" .*$"))))
}
Filter.TF_Perturbations_Followed_by_Expression<-function(x){
  return(str_to_upper(unlist(lapply(x,str_remove,pattern=" .*$"))))
}
Filter.LOF_Expression_from_GEO<-function(x){
  return(str_to_upper(unlist(lapply(x,str_remove,pattern=" .*$"))))
}
Filter.Transcription_Factor_PPIs<-function(x){
  return(str_to_upper(x))
}
Filter.TRANSFAC_and_JASPAR_PWMs<-function(x){
  return(str_to_upper(unlist(lapply(x,str_remove,pattern="[ :-].*$"))))
}
Filter.TRRUST_Transcription_Factors_2019<-function(x){
  return(str_to_upper(unlist(lapply(x,str_remove,pattern=" .*$"))))
}

# ####################################################################
Get.ENRICHR.Results<-function(t,nRes,doPlot = TRUE, minSig = .05, type = NA, fromFile = F){
  RetValues<-list()
  mDBs<-"ChEA_2016,ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X,ENCODE_TF_ChIP-seq_2015,TF_Perturbations_Followed_by_Expression,TF-LOF_Expression_from_GEO,Transcription_Factor_PPIs,TRANSFAC_and_JASPAR_PWMs,TRRUST_Transcription_Factors_2019"
  mDBs<-unlist(str_split(mDBs,","))

  if (fromFile & !is.na(type)){
    mFile = paste0(type,"-Enrichr_TFs.rda")
    load(mFile)
  }
  else mRes <- enrichr(t, mDBs) #call to ENRICHR website

  if(!is.na(type) & !fromFile){
    save(mRes, file = paste0(type,"-Enrichr_TFs.rda"))
  }
  mData<-FilterResults(mRes,"TF",nRes, minSig)
  nRows<-min(nRes,nrow(mData))
  if(doPlot){
    pDat<-data.frame(name = factor(mData$Term[1:nRows], levels = rev(mData$Term[1:nRows])),
                     padj = mData$padj[1:nRows])
    colnames(pDat)<-c("name","padj")
    PrintPlot(pDat[1:50,])
  }
  RetValues[[1]]<-mData[1:nRows,]
  names(RetValues)<-c("TFs")
  return(RetValues)
}

# ###################### Filter Results ##########################################
FilterResults<-function(mRes, type, nRes, minSig){
  mDBs<-"ChEA_2016,ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X,ENCODE_TF_ChIP-seq_2015,TF_Perturbations_Followed_by_Expression,TF-LOF_Expression_from_GEO,Transcription_Factor_PPIs,TRANSFAC_and_JASPAR_PWMs,TRRUST_Transcription_Factors_2019"
  mDBs<-unlist(str_split(mDBs,","))
  allDBs<-data.frame(libraryName = mDBs, ID = seq(1:length(mDBs)))

  for (i in 1:length(mRes)){ #keep only things that are significant as defined by minSig
    temp<-mRes[[i]]
    temp<-temp[temp[,4] < minSig,]
    mRes[[i]]<-temp
  }

  for (i in 1:length(mRes)){ #make the terms comparable
    temp<-mRes[[i]]
    if (nrow(temp) >= 1){
      temp$original.term<-temp$Term #add the original term returned by enrichr
      temp<-Filter.TFs(temp,names(mRes)[i])
    }
    mRes[[i]]<-temp
  }

  #mRes is a list, I want to annotate the terms with the index of the database..
  for(i in 1:length(mRes)){
    mRes[[i]]$raw.term<-mRes[[i]]$Term
    index <- allDBs[match(names(mRes)[i],allDBs$libraryName),]$ID
    suffix <- paste0("-(",index,")", collapse = "")
    if (nrow(mRes[[i]]>0)){
      mRes[[i]]$Term <- paste0(mRes[[i]]$Term, suffix)
    }
  }

  mData<-do.call("rbind", mRes)
  mData$lib<-rownames(mData)
  mData<-mData[order(mData$Adjusted.P.value, decreasing = FALSE),]
  mData$padj <- -log10(mData$Adjusted.P.value)
  mNames<-names(table(mData$raw.term))
  mIndex<-match(mNames,mData$raw.term)
  mData<-mData[mIndex,]
  mData<-mData[order(mData$Adjusted.P.value, decreasing = FALSE),]
  nRows<-min(nRes,nrow(mData))
  return(mData)
}

# ########################################################################################
# 
# #this does not generate figures, it just formats the data for export to the supplemental data table
HandleCommonGenes<-function(FS.res.all, PNST.res.all, PWT.res.all,index,nRes,nTerms, trimNames = F, genesReturn = "both"){
  #index sets which result is used (TFs, kinases, drugs, pathways)
  FS.res<-FS.res.all[[index]]
  PNST.res<-PNST.res.all[[index]]
  PWT.res<-PWT.res.all[[index]]

  #sort by padj, smallest values first
  FS.res<-FS.res[order(FS.res[,'Adjusted.P.value'], decreasing = F),]
  PNST.res<-PNST.res[order(PNST.res[,'Adjusted.P.value'], decreasing = F),]
  PWT.res<-PWT.res[order(PWT.res[,'Adjusted.P.value'], decreasing = F),]

  #number of results to consider (or all, if nrow is less than nterms)
  FS.res<-FS.res[1:(min(nTerms,nrow(FS.res))),]
  PNST.res<-PNST.res[1:(min(nTerms,nrow(PNST.res))),]
  PWT.res<-PWT.res[1:(min(nTerms,nrow(PWT.res))),]

  #get the list of terms to combine
  mTerm<-"raw.term"
  all.terms<-unique(c(FS.res[,mTerm], PNST.res[,mTerm],PWT.res[,mTerm]))

  #matrix to store values for plotting, will get 'NA' if the term doesn't exist for the tumor type.
  mDat<-matrix(ncol = 3, nrow = length(all.terms))
  rownames(mDat)<-all.terms
  colnames(mDat)<-c("FS","PNST","PWT")

  padj<-"Adjusted.P.value"

  mDat[,1]<- -log10(FS.res[match(all.terms,FS.res[,mTerm]),padj])
  mDat[,2]<- -log10(PNST.res[match(all.terms,PNST.res[,mTerm]),padj])
  mDat[,3]<- -log10(PWT.res[match(all.terms,PWT.res[,mTerm]),padj])

  #keep terms that are present in all three tumor types
  mFun<-function(x, thresh = 2){
    return(sum(is.na(x)) < thresh)
  }
  keep<-apply(mDat,1,mFun, thresh = 1)


    cDat<-mDat[keep,]
    cDat[is.na(cDat)]<-0
    pDat<-cDat

    if (trimNames){
      #handle long term names
      for(i in 1:nrow(pDat)){
        mWords<-unlist(str_split(rownames(pDat)[i], " "))
        if (length(mWords) > 3){
          rownames(pDat)[i]<-paste(paste(mWords[1:3],collapse = " "),"-\n        ",
                                   paste(mWords[4:length(mWords)],collapse = " "), sep = "")
          }
        }
      } #if trimnames


  #format the data to the clipboard to put in supplemental data....
  if (sum(keep) < 2){
    mTerms<-names(keep[keep])
  } else {
    mTerms<-rownames(cDat)
  }

  mCData<-FS.res[FS.res$raw.term %in% mTerms,]
  output<-ProcessHeatMapRes(mCData,PC_DEGs_FS_tumors_vs_cSKIN,nRes,genesReturn = genesReturn)
  colnames(output)[1]<-"term"
  output$spacer<-" "; #add a column of blank spaces so the long names don't show up in the csv file.

  mCData<-PNST.res[PNST.res$raw.term %in% mTerms,]
  temp<-ProcessHeatMapRes(mCData,PC_DEGs_PNST_tumors_vs_cSKIN,nRes,genesReturn = genesReturn)
  colnames(temp)[1]<-"term"
  temp$spacer<-" "
  output<-cbind(output,temp)

  mCData<-PWT.res[PWT.res$raw.term %in% mTerms,]
  temp<-ProcessHeatMapRes(mCData,PC_DEGs_PWT_tumors_vs_cVSMS,nRes,genesReturn = genesReturn)
  colnames(temp)[1]<-"term"
  output<-cbind(output,temp)
  return(output)

}

HandleDotPlots<-function(index, type, nRes){
  FS.unique.data<-FS.res[[index]]
  PNST.unique.data<-PNST.res[[index]]
  PWT.unique.data<-PWT.res[[index]]


  if (nrow(FS.unique.data) > nRes)
    FS.unique.data<-FS.unique.data[1:nRes,]
  if (nrow(PNST.unique.data) > nRes)
    PNST.unique.data<-PNST.unique.data[1:nRes,]
  if (nrow(PWT.unique.data) > nRes)
    PWT.unique.data<-PWT.unique.data[1:nRes,]

  maxPadj<-max(c(FS.unique.data$padj, PNST.unique.data$padj,PWT.unique.data$padj))
  temp<-c(FS.unique.data[,2],PNST.unique.data[,2],PWT.unique.data[,2])
  temp<-unlist(lapply(temp,getSize))
  maxCount<-max(as.numeric(temp))
  #round scales to nearest 5
  maxPadj<-ceiling(maxPadj/5)*5
  maxCount<-ceiling(maxCount/5)*5

  mFSData<-FS.unique.data
  pDat<-data.frame(name = factor(mFSData$raw.term, levels = rev(mFSData$raw.term)),
                   padj = mFSData$padj, count = as.numeric(lapply(mFSData[,2], getSize)))
  PrintDotPlot(pDat, colorLimits = c(1.3,maxPadj), sizeLimits = c(1,maxCount), maxSize = 12)
  mTerms<-as.character(pDat$name)
  mCData<-FS.unique.data[FS.unique.data$raw.term %in% mTerms,]

  if(type == "drivers"){ #next calls are going to sort genes by expression data
    cDat<-ProcessDotPlotRes(mCData,PC_DEGs_FS_tumors_vs_cSKIN,nRes = nRes)
  } else {
    cDat<-ProcessDotPlotRes(mCData,FS.exp,nRes = 15)
  }


  mPNSTData<-PNST.unique.data
  pDat<-data.frame(name = factor(mPNSTData$raw.term, levels = rev(mPNSTData$raw.term)),
                   padj = mPNSTData$padj, count = as.numeric(lapply(mPNSTData[,2], getSize)))
  PrintDotPlot(pDat, colorLimits = c(1.3,maxPadj), sizeLimits = c(1,maxCount), maxSize = 12)
  mTerms<-as.character(pDat$name)
  mCData<-PNST.unique.data[PNST.unique.data$raw.term %in% mTerms,]

  if(type == "drivers"){
    temp<-ProcessDotPlotRes(mCData,PC_DEGs_PNST_tumors_vs_cSKIN,nRes = nRes)
  } else {
    temp<-ProcessDotPlotRes(mCData,PNST.exp,nRes = nRes)
  }
  cDat<-cbind(cDat,temp)


  mPWTData<-PWT.unique.data
  pDat<-data.frame(name = factor(mPWTData$raw.term, levels = rev(mPWTData$raw.term)),
                   padj = mPWTData$padj, count = as.numeric(lapply(mPWTData[,2], getSize)))
  PrintDotPlot(pDat, colorLimits = c(1.3,maxPadj), sizeLimits = c(1,maxCount), maxSize = 12)
  mTerms<-as.character(pDat$name)
  mCData<-PWT.unique.data[PWT.unique.data$raw.term %in% mTerms,]
  # colnames(mCData)[10]<-"Term"
  # colnames(mCData)[11]<-"db"
  # mCData<-mCData[,c(10,2,3,4,9,11)]
  if(type == "drivers"){
    temp<-ProcessDotPlotRes(mCData,PC_DEGs_PWT_tumors_vs_cVSMS,nRes = nRes)
  } else {
    temp<-ProcessDotPlotRes(mCData,PWT.exp,nRes = nRes)
  }
  cDat<-cbind(cDat,temp)

  return(cDat)
}

FinishDotPlotRes<-function(tRes){
  colnames(tRes)<-make.unique(colnames(tRes))
  tRes<-tRes[,c("raw.term","P.value","Adjusted.P.value","nUpGenes","Genes","db","original.term",
                "raw.term.1","P.value.1","Adjusted.P.value.1","nUpGenes.1","Genes.1","db.1","original.term.1",
                "raw.term.2","P.value.2","Adjusted.P.value.2","nUpGenes.2","Genes.2","db.2","original.term.2")]
  write.clip(tRes)
}
