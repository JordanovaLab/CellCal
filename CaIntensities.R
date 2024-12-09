library(dplyr)
library(ggplot2)
library(magrittr)
library(gridExtra)
library(ggbeeswarm)
library(data.table)
library(tidyverse)
library(tidyr)
library(writexl)
library(poolr)
library(broom)
library(readxl)

Cyto <- read.table("Intensities.txt", sep = "\t", header = T)
Cyto <- Cyto [-c(6:8)]

#Add Genotype and Treatment
Cyto$Geno<-factor(substr(as.character(Cyto$ImageName),1,7))
Cyto$Treatment<-factor(substr(as.character(Cyto$ImageName),9,11))
Cyto$Well<-factor(substr(as.character(Cyto$ImageName),13,14))

#Adding Identifiers per measurement

Cyto<- mutate(Cyto, Identifier=paste(Geno, ImageName, ROI_Number, Well, Treatment))
Cyto<- mutate(Cyto, GenoTreat = paste(Geno, Treatment))
Cyto<- mutate(Cyto, ImageName2=ImageName)

#Table Before Calculations and after defining groups

write.table(Cyto,"Cyto_intensities.csv",sep="\t",col.names=TRUE,quote=F,row.names=F)

#Cyto <- read.table("Cyto_intensities.csv", sep = "\t", header = T)

#Checking ROI types
Cyto$ROI.type

#SUBSTRACTING BACKGROUND

#Selection of background ROIs
Cyto_FBGtemp <-Cyto[which(Cyto$ROI.type==0),]

#Selection and removal fo background
Cyto_FBG<-subset(Cyto_FBGtemp, select= c('Geno', "ImageName2", "Frame", "Mean_Intensity")) 
Cyto_FBG<-rename(Cyto_FBG, BG_Int_ROI = Mean_Intensity)
Cyto=left_join(Cyto, Cyto_FBG) 
rm(Cyto_FBG, Cyto_FBGtemp)
Cyto<-mutate(Cyto, Int_BG = Mean_Intensity-BG_Int_ROI)


#CALCULATING F/F0
F0_temp<-Cyto[which(Cyto$ROI.type==4),] 

F0Int<- subset(F0_temp, select=c('Geno', "ImageName2", "Frame", "Identifier", "Mean_Intensity"))
F0Int <-summarise(group_by(Cyto, Identifier), mean(Int_BG[5]))
names(F0Int)[2]<-"F0"
Cyto=left_join(Cyto, F0Int)

rm(F0Int,F0_temp)

Cyto<-Cyto%>%
  filter (F0>0)

Cyto<-mutate(Cyto, F_F0=Int_BG/F0)



#Making average
NBGtemp<-Cyto[which(Cyto$ROI.type==4),] 

Cyto_NBG<- subset(NBGtemp, select=c('Geno', "ImageName2", "Frame", "Identifier", "Mean_Intensity", "F_F0","GenoTreat"))
Cyto_NBG<-mutate (Cyto_NBG, GenoTreatFrame = paste (GenoTreat, Frame))
MF_F0 <-summarise(group_by(Cyto_NBG, GenoTreatFrame), mean(F_F0), sd(F_F0))

Cyto_NBG=left_join(Cyto_NBG,MF_F0)
names(Cyto_NBG)[9]<-"Average"
names(Cyto_NBG)[10]<-"SDev"

rm(NBGtemp)

rm(MF_F0)

#Table after calculations
write.table(Cyto_NBG,"Cyto_NBG.csv",sep="\t",col.names=TRUE,quote=F,row.names=F)

#Selecting strong responders
F_F0limit<-2
react<-summarise(group_by(Cyto_NBG, Identifier), max(F_F0))
names(react)[2]<-"max_F_F0"
react$react<-"0"
react$react[react$max_F_F0>F_F0limit]<-1
Cyto_NBG=left_join(Cyto_NBG, react)
react$Geno<-factor(substr(as.character(react$Identifier),1,7))
react$Treatment<-factor(substr(as.character(react$Identifier),17,19))


#Removing negative values
Positives<-summarise(group_by(Cyto_NBG,Identifier),min(F_F0))
names(Positives)[2]<-"minF_F0"
Positives$Positives<-"0"
Positives$Positives[Positives$minF_F0>0]<-1
Cyto_NBG=left_join(Cyto_NBG,Positives)
rm(Positives)




#Make table with Area under the curve
AUC_F_F0=summarise(group_by(Cyto_NBG, Identifier,GenoTreat),
                   sum(diff(Frame) * (head(F_F0,-1)+tail(F_F0,-1)))/2)
names(AUC_F_F0)[3]<-'AUC_F_F0'

AUC_F_F0$Geno<-factor(substr(as.character(AUC_F_F0$Identifier),1,7))
AUC_F_F0$Treatment<-factor(substr(as.character(AUC_F_F0$Identifier),17,18))


write.table(AUC_F_F0,"AUC_NBG.csv",sep="\t",col.names=TRUE,quote=F,row.names=F)
write.table(Cyto_NBG,"Cyto_allCalc_NBG.csv",sep="\t",col.names=TRUE,quote=F,row.names=F)
write.table(react,"React_Cyto_all_merged_NBG.csv",sep="\t",col.names=TRUE,quote=F,row.names=F)
#write.table(Cells,"Cell_count_ATP.csv",sep="\t",col.names=TRUE,quote=F,row.names=F)
