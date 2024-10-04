#!/usr/bin/env Rscript


folder<-"/home/samman/Documents/MyGitHub/ChickpeaOmicsR/data/normalizing/Count/Normalized"
#read all files in the folder
files<-list.files(folder)
read_data<-function(file){
  data<-read.table(file,header=TRUE,sep=",",row.names=1)
  #remove the first 2 columns
  data<-data[,-c(1,2)]
  return(data)
}
#read all files
data<-lapply(files,function(x) read_data(paste(folder,x,sep="/")))
#get all rownames
allrownames<-lapply(data,function(x) rownames(x))
allrownames<-unlist(allrownames)
allrownames<-unique(allrownames)
allcolnams<-lapply(data,function(x) colnames(x))
allcolnams<-unlist(allcolnams)

numberOfColumns<-lapply(data,function(x) ncol(x))
numberOfColumns<-unlist(numberOfColumns)
numberOfColumns<-sum(numberOfColumns)
# create a matrix with all rownames
mergeData<-matrix(NA,nrow=length(allrownames),ncol=numberOfColumns)
rownames(mergeData)<-allrownames
colnames(mergeData)<-allcolnams

#fill the matrix
for(i in 1:length(data)){
  for(j in 1:ncol(data[[i]])){
    mergeData[rownames(data[[i]]),colnames(data[[i]])[j]]<-data[[i]][,j]
  }
}
#write the merged data
write.table(mergeData,file=paste(folder,"/merged.csv",sep="/"),sep=",",quote=FALSE)
 


