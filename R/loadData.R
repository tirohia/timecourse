# this file loads data from 3 files. One containing read counts that map to the host genome,
# one gene per row, one colum per sample. One containing read counts for reads map to the
# pathogen and one file that assigns column numbers to condition and time. An example file
# called lookupsheet is included, showing the assignments for a timecourse with 3 replicates
# for each time point, 10 time points, two treatments. If anyone can think of a better way to
# do this, it would be much appreciated.

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
rawCounts<-as.matrix(rawCounts)

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
colData<-DataFrame(pdataclean)


gffRangedData<-import.gff3("Kiwifruit_pseudomolecule.gff3")
gffRange<-gffRangedData[gffRangedData$type=="gene",]
gffRange<-gffRange[-1,]
gffRanges<-as(gffRange, "GRanges")


SummarizedExperiment(assays=list(counts=rawCounts), colData=colData)

#                     rowRanges=gffRanges )
#                     metadata=list(metadata))
counts
colData

xrawCounts<-rawCounts
dim(rawCounts)
dim(colData)
dim(counts)
counts<-as.matrix(counts)
class(counts)
rownames(counts)<-NULL
counts<-rawCounts[1:200,1:6]
colData<-colData[1:6,]
dim(counts)
length(gffRanges)








colData <- DataFrame(Treatment=rep(c("ChIP", "Input"), 3), row.names=LETTERS[1:6],time=c(1:6))
colData
rowRanges <- GRanges(rep(c("chr1", "chr2"), c(50, 150)),
                     IRanges(floor(runif(200, 1e5, 1e6)), width=100),
                     strand=sample(c("+", "-"), 200, TRUE),
                     feature_id=sprintf("ID%03d", 1:200))


rowRanges
nrows <- 200
ncols <- 6
counts1 <- matrix(runif(nrows * ncols, 1, 1e4), nrows)
dim(counts)
colnames(counts)<-NULL
dim(counts1)
class(counts)
head(counts)
head(counts1)
rownames(rawCounts)<-NULL
colnames(rawCounts)<-1:60
counts<-as.matrix(rawCounts[1:200,1:20])
head(counts)
se<-SummarizedExperiment(assays=list(counts=counts),colData=colData[1:20,])
                     rowRanges=rowRanges,


coldata <- DataFrame(pdataclean)

colData(se)
head(assay(se))
rowRanges(se)

dim(myGranges)
class(gffRangedData)
gfffile<-gffRangedData[gffRangedData$type=="gene",]
gfffile
coldata

length(rownames(rawCounts))
tail(myGranges)

biocLite("rtracklayer")
library("rtracklayer")

gse
gse<-fission
colnames(rawCounts)<-pdataclean$id

library("SummarizedExperiment")
coldata <- DataFrame(pdataclean)



genes<-rownames(rawCounts)


hostSE
pathogenSE
