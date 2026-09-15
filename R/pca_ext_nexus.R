#Load libraries
library(readxl)
library(ggplot2)

#Read in metadata
metadata <- read_xlsx("/Users/andrewblack/Documents/Research/GROUSE/sarek_nexus_new_plus_shotguns/nexus_extension_heterozygosity.xlsx")                             

#Read in covariation matrix
cov<-as.matrix(read.table("~/final.cov"))

#Extract and calculate eplained variation
axes<-eigen(cov)
head(axes$values/sum(axes$values)*100)
#[1] 22.6180474  2.9970624  2.1391689  0.6050842  0.5466266
#[6]  0.4565761


#Bind vectors with metadata and plot
PC1_3<-as.data.frame(axes$vectors[,1:3])
x<-cbind(PC1_3,metadata)
 #By species and DPS
ggplot(data=x, aes(y=V2, x=V1))+geom_point(size=7,color="black",aes(shape=metadata$DPS,fill=metadata$SPECIES))+ theme_classic() + xlab("PC1 (22.6%)") +ylab("PC2 (3.0%)")+geom_hline(yintercept=0,linetype="dashed")+geom_vline(xintercept =0,linetype="dashed")+scale_fill_manual("Species", values=c("goldenrod","brown","black","grey"))+scale_shape_manual("DPS", values=c(25,21,21,21))+ theme(legend.position = "top")

