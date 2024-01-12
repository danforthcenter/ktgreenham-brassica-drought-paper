# Script for analysis of VIS image analysis and water data
library(ggplot2)
library(lubridate)
library(MASS)
library(plyr)
library(car)
library(RColorBrewer)
library(Hmisc)
library(gridExtra)
library(fields)
library(gtools)
library(gplots)
library(scales)
library(reshape)
library(moments)
library(nlme)
library(grid)
library(lme4)
library(gganimate)
library(plotly)
library(data.table)
library(magick)


setwd("/Users/mgehan/Documents/kt-greenham-brassica/brassica-2023/")

############################################
# Brassica green color over time line graphs
############################################
vis.data = read.table(file="VIS_SV_color_results_without_yellow_leaves.csv",sep=",",header=TRUE, stringsAsFactors = FALSE)

vis.data.parents<-vis.data[!grepl("x",vis.data$Genotype.ID),]
vis.data.parents<-vis.data.parents[!grepl("R500-43",vis.data.parents$Genotype.ID),]

uniquegenotypes=unique(vis.data.parents$Genotype.ID)

vis.data.subset<-data.frame(vis.data.parents$plantbarcode,vis.data.parents$days,
                            vis.data.parents$Treatment.2,
                            vis.data.parents$Genotype.ID, 
                            vis.data.parents$frame,
                            vis.data.parents$percent_scenesed)

colnames(vis.data.subset)<-c("barcode","day","treatment","genotype","frame","percent_yellow")

vis.data.subset$percent_green<-1-vis.data.subset$percent_yellow

vis.data.subset$treatment[grep("100", vis.data.subset$treatment)] <- 'control'      	    
vis.data.subset$treatment[grep("20", vis.data.subset$treatment)] <- 'drought'      	    

green.avg=ddply(vis.data.subset, c("genotype","treatment","day"), summarise, N=length(percent_green), mean=mean(percent_green), sd=sd(percent_green),se=sd/sqrt(N))

yellow.avg=ddply(vis.data.subset, c("genotype","treatment","day"), summarise, N=length(percent_yellow), mean=mean(percent_yellow), sd=sd(percent_yellow),se=sd/sqrt(N))

pdf(file="brassica_percent_yellow.pdf",height=7,width=7,useDingbats=FALSE)
ggplot(yellow.avg,aes(x=day,y=mean, group= treatment, colour = factor(treatment))) +
  geom_smooth()+
  scale_color_manual(labels = c("WW", "WL"),values=c("chartreuse4", "darkgoldenrod2"))+
  labs(x= "Days", y="Mean Percentage Yellow Area", col="Treatment")+
  theme_bw()+
  theme(legend.text = element_text(size = 14),legend.title = element_text(size = 14), axis.title = element_text(size=14))+
  facet_wrap( ~ genotype, ncol = 4)+
  theme(axis.text=element_text(size=12), strip.text.x = element_text(size = 12))
dev.off()

############################################
# Brassica green color over time line graphs
############################################

vis.data = read.table(file="VIS_SV_color_results_without_yellow_leaves.csv",sep=",",header=TRUE, stringsAsFactors = FALSE)

vis.data.parents<-vis.data[!grepl("x",vis.data$Genotype.ID),]
vis.data.parents<-vis.data.parents[!grepl("R500-43",vis.data.parents$Genotype.ID),]

uniquegenotypes=unique(vis.data.parents$Genotype.ID)

vis.data.subset<-data.frame(vis.data.parents$plantbarcode,vis.data.parents$days,
                            vis.data.parents$Treatment.2,
                            vis.data.parents$Genotype.ID, 
                            vis.data.parents$frame,
                            vis.data.parents$hue_circular_mean)

colnames(vis.data.subset)<-c("barcode","day","treatment","genotype","frame","hue_circular_mean")

vis.data.subset$treatment[grep("100", vis.data.subset$treatment)] <- 'control'      	    
vis.data.subset$treatment[grep("20", vis.data.subset$treatment)] <- 'drought'      	    

hue.avg=ddply(vis.data.subset, c("genotype","treatment","day"), summarise, N=length(hue_circular_mean), mean=mean(hue_circular_mean), sd=sd(hue_circular_mean),se=sd/sqrt(N))

pdf(file="brassica_hue.pdf",height=12,width=12,useDingbats=FALSE)
ggplot(hue.avg,aes(x=day,y=mean, group= treatment, colour = factor(treatment))) +
  geom_smooth()+
  scale_color_manual(labels = c("WW", "WL"),values=c("chartreuse4", "darkgoldenrod2"))+
  labs(x= "Days", y="Hue Circular Mean", col="Treatment")+
  theme_bw()+
  theme(legend.text = element_text(size = 14),legend.title = element_text(size = 14), axis.title = element_text(size=14))+
  facet_wrap( ~ genotype, ncol = 4)+
  theme(axis.text=element_text(size=12), strip.text.x = element_text(size = 12))
dev.off()

############################################
# Brassica green color over time heatmaps
############################################

hue.avg.subset<-data.frame(hue.avg$genotype,hue.avg$treatment, hue.avg$day,hue.avg$mean)
colnames(hue.avg.subset)<-c("genotype","treatment","day","hue.avg")

hue.avg.control<-hue.avg.subset[(hue.avg.subset$treatment=='control'),]
hue.avg.drought<-hue.avg.subset[(hue.avg.subset$treatment=='drought'),]

control.cast<-cast(hue.avg.control,genotype~day,value ='hue.avg')
drought.cast<-cast(hue.avg.drought,genotype~day,value ='hue.avg')

#get order of genotypes from one of the datasets
distance    = dist(control.cast)
cluster     = hclust(distance, method="ward.D")
dendrogram  = as.dendrogram(cluster)

color.palette = colorRampPalette(c("lightyellow","lightgreen","darkgreen"),space="rgb")
pdf(file="control.hue.pdf",width = 7,height = 7,pointsize = 8,useDingbats = FALSE)
control.hue<- heatmap.2(as.matrix(control.cast),
                  #Rowv=dendrogram,
                  Colv=FALSE,
                  Rowv=FALSE,
                  #dendrogram='row',
                  scale='none',
                  col=color.palette,
                  trace='none',
                  symbreaks=FALSE,
                  sepwidth=c(0.01,0.01),
                  colsep=1:ncol(control.cast),
                  rowsep=1:nrow(control.cast),
                  sepcolor="darkgray",
                  cex.main=0.75,
                  keysize=1.4,
                  breaks = seq(60,120, length.out = 120),
                  xlab="Days",
                  ylab="Genotype",
                  main="Well-Watered Hue Circular Mean")
dev.off()

pdf(file="drought.hue.pdf",width = 7,height = 7,pointsize = 8,useDingbats = FALSE)
control.hue<- heatmap.2(as.matrix(drought.cast),
                       # Rowv=dendrogram,
                        Colv=FALSE,
                        Rowv=FALSE,
                        #dendrogram='row',
                        scale='none',
                        col=color.palette,
                        trace='none',
                        symbreaks=FALSE,
                        sepwidth=c(0.01,0.01),
                        colsep=1:ncol(control.cast),
                        rowsep=1:nrow(control.cast),
                        sepcolor="darkgray",
                        cex.main=0.75,
                        keysize=1.4,
                        breaks = seq(60,120, length.out = 120),
                        xlab="Days",
                        ylab="Genotype",
                        main="Water-Limited Hue Circular Mean")
dev.off()




