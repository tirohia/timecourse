# this file loads data from 3 files. One containing read counts that map to the host genome,
# one gene per row, one colum per sample. One containing read counts for reads map to the
# pathogen and one file that assigns column numbers to condition and time. An example file
# called lookupsheet is included, showing the assignments for a timecourse with 3 replicates
# for each time point, 10 time points, two treatments. If anyone can think of a better way to
# do this, it would be much appreciated.

#source("https://bioconductor.org/biocLite.R")
#biocLite()
#biocLite("DESeq2")

library(DESeq2)

setwd("/home/ben/googleDrive/code/timecourse/R")

design<-read.csv("lookupsheet.csv")
columnOrder<-as.vector(t(design[,3:5]))
columnNames<-c(paste("control",1:30,sep=""),paste("inoculated",1:30,sep=""))
times<-c(0,1.5,3,6,12,24,48,72,96,120)
pData<-data.frame(treatment=c(rep("control",each=30),rep("inoculated",each=30)),time=(rep(times,each=3)))
rownames(pData)<-columnNames
metadata <- data.frame(labelDescription= c("Sample treatment","timepoint"),row.names=c("treatment", "time"))
phenoData <- new("AnnotatedDataFrame", data=pData, varMetadata=metadata)

experimentData <- new("MIAME",
                      name="Ben Curran",
                      lab="PFR",
                      contact="b.curran@auckland.ac.nz",
                      title="Kiwifruit-Psa time course experiment",
                      abstract="An examination of the initial stages of pathogen incursion",
                      url="www.pfr.co.nz"
)


rawCounts<-read.csv("cornellGenomeModels-FM-ordered.csv",row.names = 1)
colnames(rawCounts)<-columnNames
hostEset <-new("ExpressionSet",exprs= as.matrix(rawCounts),phenoData=phenoData,experimentData=experimentData)
rawCounts<-read.csv("psaGenomeModels-FM-ordered.csv",row.names = 1)
colnames(rawCounts)<-columnNames
pathogenEset <-new("ExpressionSet",exprs= as.matrix(rawCounts),phenoData=phenoData,experimentData=experimentData)

hostEset
pathogenEset


