library(stringr)
library(readr)
library(edgeR)
library(DESeq2)

#change this for your own computer:
baseDir<-"C:\\Users\\idaho\\Dropbox\\"
baseDir<-"C:\\Users\\Eric Shelden\\Dropbox\\"

#################### Prepare the data ####################################
#First set of scrolls
setwd(str_c(baseDir,"Shelden lab\\Canine SARC analysis\\raw_reads\\CanFam-Newest\\FirstScrolls"))
mDir<-getwd()
list.files(mDir)
mFiles<-list.files(mDir)
mFiles<-mFiles[str_detect(mFiles,".tsv")]
mFiles

mFun<-function(x){
  read_delim(x, "\t", escape_double = FALSE, 
                      col_names = FALSE, trim_ws = TRUE)
}

mDat<-lapply(mFiles,mFun)
mDat<-as.data.frame(mDat)
dim(mDat) #41823 x 54
rownames(mDat)<-mDat[,1]
mDat<-mDat[,seq(2,ncol(mDat),2)]

cNames<-unlist(str_split("2c,3c,4c,5c,7c,9c,10c,11c,12c,2b,2a,3b,3a,4b,4a,5b,5a,7b,7a,9b,9a,10b,10a,11b,11a,12b,12a",","))
colnames(mDat)<-cNames
#order data columns by name
mDat<-mDat[,order(colnames(mDat))]
dim(mDat) #41823 x 27 or 3x9 samples

#add the reads from all the runs together
colnames(mDat)

#use DESEQ collapse function to combine technical replicates for First scrolls.
#deseq needs some column data in order to create a valid deseq object
sData<-data.frame(sID<-colnames(mDat))
colnames(mDat)
save(mDat, file = "first_scrolls_replicates.rda")
dds<-DESeqDataSetFromMatrix(mDat, sData, design = ~1)

#collapse function needs factor of samples to group together, I call this reps
reps<-colnames(mDat)
reps<-str_remove_all(reps,"[a-z]")
reps<-as.factor(reps)
reps
ddsCol<-collapseReplicates(dds,reps)
colnames(ddsCol)
dds<-ddsCol
rm(ddsCol)
cNames<-unlist(str_split("S1-10,S1-11,S1-12,S1-2,S1-3,S1-4,S1-5,S1-7,S1-9",","))
mSums<-counts(dds)
rm(dds)
colnames(mSums)<-cNames
mSums<-as.data.frame(mSums)
rownames(mSums)<-rownames(mDat)

#reorder mSums (we don't want them alphabetical)
mSums<-mSums[,c(4,5,6,7,8,9,1,2,3)]
colnames(mSums) #"S1-2"  "S1-3"  "S1-4"  "S1-5"  "S1-7"  "S1-9"  "S1-10" "S1-11" "S1-12"
dim(mSums) #41823 x 9
tail(mSums)

#get rid of the 5 lines of garbage at the end of the htseq-count file
S1Data<-mSums[1:(nrow(mSums)-5),]
save(S1Data, file = str_c(baseDir,"Shelden lab\\Canine SARC analysis\\Input Data\\S1DataFinal.rda"))

########################
#second set of scrolls
setwd(str_c(baseDir,"Shelden lab\\Canine SARC analysis\\raw_reads\\CanFam-Newest\\Second scrolls"))
mDir<-getwd()
mFiles<-list.files(mDir)
mFiles<-mFiles[str_detect(mFiles,".tsv")]
mFiles
mDat<-lapply(mFiles, mFun)
mDat<-as.data.frame(mDat)
rownames(mDat)<-mDat[,1]
#View(mDat)
dim(mDat) #41823 x 14
mDat<-mDat[,seq(2,ncol(mDat),2)]
cNames<-unlist(str_split("S2-3,S2-4,S2-6,S2-7,S2-8,S2-9,S2-11",","))
colnames(mDat)<-cNames
S2Data<-mDat
#remove last 5 lines added by HTSEqCount
S2Data<-S2Data[1:(nrow(S2Data)-5),]
tail(S2Data)
all.equal(rownames(S1Data),rownames(S2Data)) #TRUE
save(S2Data,file = str_c(baseDir,"Shelden lab\\Canine SARC analysis\\Input Data\\S2DataFinal.rda"))

########################
#smooth muscle cell data Vascular_smooth_muscle
mDir<-str_c(baseDir,"Shelden lab\\Canine SARC analysis\\raw_reads\\CanFam-Newest\\Vascular_smooth_muscle")
setwd(mDir)
mFiles<-list.files(mDir)
mFiles<-mFiles[str_detect(mFiles,".tsv")]
mFiles
mDat<-lapply(mFiles, mFun)
mDat<-as.data.frame(mDat)
rownames(mDat)<-mDat[,1]
mDat<-mDat[,seq(2,ncol(mDat),2)]

GetName<-function(x){
  return(unlist(str_split(x,"\\."))[1])
}
cNames<-unlist(lapply(mFiles,GetName))
cNames
colnames(mDat)<-cNames
vSMDat<-mDat
#remove last 5 lines added by HTSEqCount
vSMDat<-vSMDat[1:(nrow(vSMDat)-5),]
tail(vSMDat)
all.equal(rownames(S1Data),rownames(vSMDat)) #TRUE
save(vSMDat,file = str_c(baseDir,"Shelden lab\\Canine SARC analysis\\Input Data\\vSMDat.rda"))
########################
#skin
mDir<-str_c(baseDir,"Shelden lab\\Canine SARC analysis\\raw_reads\\CanFam-Newest\\Canine Skin")
setwd(mDir)
mFiles<-list.files(mDir)
mFiles<-mFiles[str_detect(mFiles,".tsv")]
mFiles
mDat<-lapply(mFiles, mFun)
mDat<-as.data.frame(mDat)
rownames(mDat)<-mDat[,1]
mDat<-mDat[,seq(2,ncol(mDat),2)]

GetName<-function(x){
  return(unlist(str_split(x,"Aligned"))[1])
}
cNames<-unlist(lapply(mFiles,GetName))
cNames
colnames(mDat)<-cNames
cSKINDat<-mDat
#remove last 5 lines
cSKINDat<-cSKINDat[1:(nrow(cSKINDat)-5),]
all.equal(rownames(S1Data),rownames(cSKINDat)) #TRUE
save(cSKINDat,file = str_c(baseDir,"Shelden lab\\Canine SARC analysis\\Input Data\\cSKINDat.rda"))
##################
#soft tissues
mDir<-str_c(baseDir,mDir<-"Shelden lab\\Canine SARC analysis\\raw_reads\\CanFam-Newest\\Capsule and Tendon")
setwd(mDir)
mFiles<-list.files(mDir)
mFiles<-mFiles[str_detect(mFiles,".tsv")]
mDat<-lapply(mFiles, mFun)
mDat<-as.data.frame(mDat)
rownames(mDat)<-mDat[,1]
mDat<-mDat[,seq(2,ncol(mDat),2)]
GetName<-function(x){
  return(unlist(str_split(x,"\\.tsv"))[1])
}
cNames<-unlist(lapply(mFiles,GetName))
cNames
colnames(mDat)<-cNames
controlSTDat<-mDat
#delete last 5 lines
controlSTDat<-controlSTDat[1:(nrow(controlSTDat)-5),]
all.equal(rownames(S1Data),rownames(controlSTDat)) #TRUE
save(controlSTDat, file = str_c(baseDir,"Shelden lab\\Canine SARC analysis\\Input Data\\controlSTDat.rda"))
########################
#normal stroma
mDir<-str_c(baseDir,"Shelden lab\\Canine SARC analysis\\raw_reads\\CanFam-Newest\\normal_stroma")
setwd(mDir)
mFiles<-list.files(mDir)
mFiles<-mFiles[str_detect(mFiles,".tsv")]
mDat<-lapply(mFiles, mFun)
mDat<-as.data.frame(mDat)
rownames(mDat)<-mDat[,1]
mDat<-mDat[,seq(2,ncol(mDat),2)]
GetName<-function(x){
  return(unlist(str_split(x,"\\.tsv"))[1])
}
cNames<-unlist(lapply(mFiles,GetName))
cNames
colnames(mDat)<-cNames
normalStroma<-mDat
normalStroma<-normalStroma[1:(nrow(normalStroma)-5),]
dim(normalStroma)
all.equal(rownames(S1Data),rownames(normalStroma)) #TRUE
save(normalStroma,file = str_c(baseDir,"Shelden lab\\Canine SARC analysis\\Input Data\\normalStroma.rda"))

############## Combine data for S1 and S2 files
setwd(str_c(baseDir,"Shelden lab\\Canine SARC analysis\\raw_reads\\CanFam-Newest"))
sTypes<-unlist(str_split("PNST,PWT,FS",","))
SampleData <- read_csv("SampleData_2.csv")

SampleData<-as.data.frame(SampleData)
rownames(SampleData)<-SampleData$ID
dim(SampleData) #23 x 11
SampleData<-SampleData[,3:10] #choose columns we need
colnames(SampleData)<-unlist(str_split("Gender,Type,Batch,Grade,Location,Size,Age,Recur",","))
SampleData
#discard frozen sample data
SampleData<-SampleData[str_detect(rownames(SampleData),"S"),]

S1_S2_data<-cbind(S1Data,S2Data)
head(S1_S2_data)

all.equal(rownames(SampleData),colnames(S1_S2_data)) #TRUE
save(S1_S2_data, file = "S1_S2_dataFinal.rda")

################################################################################################

setwd("C:/Users/idaho/Dropbox/Shelden lab/Canine SARC analysis/raw_reads/CanFam-Newest")
setwd("C:/Users/Eric Shelden/Dropbox/Shelden lab/Canine SARC analysis/raw_reads/CanFam-Newest")

load("S1_S2_dataFinal.rda")

#combine all data into one data frame and save it

mDat<-cbind(S1_S2_data,vSMDat,cSKINDat,controlSTDat,normalStroma)

sTypes<-c(rep("cSTS",16), rep("cVSMS", 9), rep("cSK",6), rep("cCT",18), rep("nStr",15))
all_data<-mDat
save(all_data, file = str_c(baseDir,"Shelden lab\\Canine SARC analysis\\Input Data\\all_data_final.rda"))
save(sTypes, file = str_c(baseDir,"Shelden lab\\Canine SARC analysis\\Input Data\\sTypes.rda"))

###################### Done preparing the data #######################################