# this file loads data from 3 files. One containing read counts that map to the host genome,
# one gene per row, one colum per sample. One containing read counts for reads map to the
# pathogen and one file that assigns column numbers to condition and time. An example file
# called lookupsheet is included, showing the assignments for a timecourse with 3 replicates
# for each time point, 10 time points, two treatments. If anyone can think of a better way to
# do this, it would be much appreciated.

# 4 contructs are created, two Expression Sets, for use with limma/edgeR etc.
# And two summarized experiments, for use with DESeq2 etc. One each for host and pathogen.

source("https://bioconductor.org/biocLite.R")
library(Biobase)

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

library("SummarizedExperiment")
library("rtracklayer")

rawCounts<-read.csv("cornellGenomeModels-FM-ordered.csv",row.names = 1)
rawCounts<-as.matrix(rawCounts) #has to be a matrix. One of those wonderful things that every single manual/tutorial etc fails to mention.

times<-c(0,1.5,3,6,12,24,48,72,96,120)
pdata<-data.frame(treatment=c(rep("control",each=30),rep("inoculated",each=30)),time=(rep(times,each=3)),replicate=(rep(c("1","2","3"),times=20)))
names(pdata) <- c("treatment","time","replicate")

pdataclean <- data.frame(treatment=ifelse(grepl("control",pdata$treatment),"con","ino"),
                         time=sub("time  \\(min\\): (.*)","\\1",pdata$time),
                         replicate=paste0("r",sub("replicate: (.*)","\\1",pdata$replicate)),
                         row.names=rownames(pdata))
pdataclean$treatment <- relevel(pdataclean$treatment, "con")
pdataclean$time <- factor(pdataclean$time, levels=unique(pdataclean$time))
pdataclean$id <- paste(pdataclean$treatment,pdataclean$time,pdataclean$replicate,sep="_")
colData<-DataFrame(pdataclean) #The DataFrame declaration is important for some reason. No idea why. Again, not mentioned in manuals/tutorials.


gffRangedData<-import.gff3("Kiwifruit_pseudomolecule.gff3")
gffRange<-gffRangedData[gffRangedData$type=="gene",]
gffRange<-gffRange[-1,]
gffRanges<-as(gffRange, "GRanges")


hostSE <- SummarizedExperiment(assays=list(counts=rawCounts), rowRanges=gffRanges,colData=colData)

rawCounts<-read.csv("psaGenomeModels-FM-ordered.csv",row.names = 1)
rawCounts<-as.matrix(rawCounts)

gffRangedData<-import.gff3("Nz13v.gff")
gffRange<-gffRangedData[gffRangedData$type=="gene",]
gffRanges<-as(gffRange, "GRanges")
gffRanges<-gffRanges[!gffRanges$locus_tag=="IYO_00005",]

pathogenSE <- SummarizedExperiment(assays=list(counts=rawCounts), rowRanges=gffRanges,colData=colData)




hostSE
hostEset
pathogenSE
pathogenEset
