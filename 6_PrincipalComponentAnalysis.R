## Load libraries

library(mixOmics)
library(phyloseq)
library(vegan)
library(ape)
library(ggplot2)


## Read in species-level table from MetaPhlAn3 and Sample metadata file
OTUIN<-otu_table(read.table("PCA_SpeciesTable.txt",header=TRUE,row.names=1),taxa_are_rows = TRUE)
SAMPDATA<-sample_data(read.table("PCA_MetaData.txt",header=TRUE,row.names=1,sep="\t"))

psobject<-phyloseq(OTUIN,SAMPDATA)

PRobject<-t(otu_table(psobject))+1

pca.result<-pca(PRobject,logratio='CLR')

coords_pca<-as.data.frame(pca.result$variates)

FinalTable<-merge(coords_pca,as.data.frame(SAMPDATA),by="row.names")

## Plotting
plot1<-ggplot(data=FinalTable,aes(X.PC1,X.PC2))+geom_point(alpha=0.0001)+theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+xlim(-4,6)+ylim(-3,4)
plot1<-plot1+geom_point(data=subset(FinalTable,Population %in% c("Abusir")), shape=17, colour="burlywood3", size=3)
plot1<-plot1+geom_point(data=subset(FinalTable,Population %in% c("American")), shape=21, fill="lightskyblue", size=3)
plot1<-plot1+geom_point(data=subset(FinalTable,Population %in% c("Californian")), shape=23, fill="deepskyblue1", size=3)
plot1<-plot1+geom_point(data=subset(FinalTable,Population %in% c("Dalheim")), shape=2, colour="gray10", size=3)
plot1<-plot1+geom_point(data=subset(FinalTable,Population %in% c("Edo")), shape=15, colour="goldenrod1", size=3)
plot1<-plot1+geom_point(data=subset(FinalTable,Population %in% c("ElSidron")), shape=22, fill="gray50", size=3)
plot1<-plot1+geom_point(data=subset(FinalTable,Population %in% c("HMP")), shape=18, colour="deepskyblue3", size=3)
plot1<-plot1+geom_point(data=subset(FinalTable,Population %in% c("Jomon")), shape=22, fill="goldenrod1", size=3)
plot1<-plot1+geom_point(data=subset(FinalTable,Population %in% c("Kilteasheen")), shape=2, colour="coral2", size=2)
plot1<-plot1+geom_point(data=subset(FinalTable,Population %in% c("Maya")), shape=1, colour="red", size=2)
plot1<-plot1+geom_point(data=subset(FinalTable,Population %in% c("Mexican")), shape=17, colour="red", size=2)
plot1<-plot1+geom_point(data=subset(FinalTable,Population %in% c("Nuragic")), shape=23, fill="yellow2", size=2)
plot1<-plot1+geom_point(data=subset(FinalTable,Population %in% c("Oneota")), shape=24, fill="deepskyblue1", size=3)
plot1<-plot1+geom_point(data=subset(FinalTable,Population %in% c("PasoDelIndio")), shape=5, colour="red", size=2)
plot1<-plot1+geom_point(data=subset(FinalTable,Population %in% c("Polyoak")), shape=15, colour="burlywood3", size=3)
plot1<-plot1+geom_point(data=subset(FinalTable,Population %in% c("Radcliffe")), shape=5, colour="lightskyblue", size=2)
plot1<-plot1+geom_point(data=subset(FinalTable,Population %in% c("Spanish")), shape=18, colour="grey30", size=3)
plot1<-plot1+geom_point(data=subset(FinalTable,Population %in% c("Wichita")), shape=21, colour="blue4", size=3, stroke=1.5)

plot1
