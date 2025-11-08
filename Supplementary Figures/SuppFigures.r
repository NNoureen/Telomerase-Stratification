##############  libraries required
library('ggplot2')
library('gplots')
library('dplyr')
library('ggpubr')
stringsAsFactors=FALSE
library(circlize)
library(stringr)
library(Hmisc)
library(forcats)
library(data.table)
library(grid)
library(ComplexHeatmap)
library(tidyverse)
library(broom)
library(stats)

##setwd('D:/TTUHSC/Project_Extend_Classification/Analysis_WithTumorCasesOnly/Manuscript/Version_102024/Version_102024/Computational_StructuralBiology/Revisions/AllFigures_andData/AnalysisData')
################################################################################################
############################ Supplementary FIGURES #################################

### Data provided in the github

###################################################################################################################################################
#############################################
#############################################               Supplementary Figure 1 
#############################################
###################################################################################################################################################


############ Supp Fig 1a : Female Ratios

Fig = read.table('SuppFig1a.txt',sep='\t',head=T)

labes1 = Fig$Cancer[which(Fig$significance == "Sig")]

Fig <- Fig %>%
  mutate(
    Size = case_when(
      Extend_Low < 200 ~ 1,
      Extend_Low >= 200 & Extend_Low < 300 ~ 2,
      Extend_Low >= 300 ~ 3,
      TRUE ~ NA_real_  # In case LowExtendFrac is NA
    )
  )


Fig$Size = as.factor(Fig$Size)


pdf('SuppFig1A.pdf')

P1 = ggscatter(Fig, x="Female_ratio_High", y="Female_ratio_Low",size ="Size",stroke = 0.6,color = "significance",palette = c('Sig'='red','NS'='skyblue'),alpha=0.6,label = "Cancer",label.select = labes1,repel=TRUE,font.label = c(7,"plain","black"))+theme_classic()+theme(axis.text.x=element_text(angle=0,size=7,vjust=0,hjust=1),legend.key.size=unit(0.05,"cm"),legend.position="top")+geom_abline(slope=1,intercept=0,linetype="dashed")


grid.newpage()
# Create layout : nrow = 2, ncol =2
pushViewport(viewport(layout = grid.layout(nrow =10, ncol = 10)))
# A helper function to define a region on the layout
define_region <- function(row, col){
  viewport(layout.pos.row = row, layout.pos.col = col)
} 

print(ggpar(P1,font.xtickslab =c("black",10),font.ytickslab =c("black",10),font.x = 10,font.y=10,font.legend=10,xlab="High EXTEND Female Ratio",ylab="Low EXTEND Female Ratio"),vp = define_region(row = 1:4, col = 1:4))


dev.off()


############ Supp Fig 1b : Older Adult Ratios

Fig = read.table('SuppFig1b.txt',sep='\t',head=T)

labes1 = Fig$Cancer[which(Fig$significance == "Sig")]

Fig <- Fig %>%
  mutate(
    Size = case_when(
      Extend_High < 100 ~ 1,
      Extend_High >= 100 & Extend_High < 200 ~ 2,
      Extend_High >= 200 ~ 3,
      TRUE ~ NA_real_  # In case HighExtendFrac is NA
    )
  )


Fig$Size = as.factor(Fig$Size)



pdf('SuppFig1b.pdf')

P1 = ggscatter(Fig, x="OA_ratio_High", y="OA_ratio_Low",size ="Size",stroke = 0.6,color = "significance",palette = c('Sig'='red','NS'='skyblue'),alpha=0.6,label = "Cancer",label.select = labes1,repel=TRUE,font.label = c(7,"plain","black"))+theme_classic()+theme(axis.text.x=element_text(angle=0,size=7,vjust=0,hjust=1),legend.key.size=unit(0.05,"cm"),legend.position="top")+geom_abline(slope=1,intercept=0,linetype="dashed")


grid.newpage()
# Create layout : nrow = 2, ncol =2
pushViewport(viewport(layout = grid.layout(nrow =10, ncol = 10)))
# A helper function to define a region on the layout
define_region <- function(row, col){
  viewport(layout.pos.row = row, layout.pos.col = col)
} 

print(ggpar(P1,font.xtickslab =c("black",10),font.ytickslab =c("black",10),font.x = 10,font.y=10,font.legend=10,xlab="High EXTEND OA Ratio",ylab="Low EXTEND OA Ratio"),vp = define_region(row = 1:4, col = 1:4))


dev.off()


############ Supp Fig 1c  : Tumor Stage Ratios 

Fig = read.table('SuppFig1c.txt',sep='\t',head=T)

labes1 = Fig$Cancer[which(Fig$significance == "Sig")]

Fig <- Fig %>%
  mutate(
    Size = case_when(
      Extend_Low < 100 ~ 1,
      Extend_Low >= 100 & Extend_Low < 200 ~ 2,
      Extend_Low >= 200 ~ 3,
      TRUE ~ NA_real_  # In case HighExtendFrac is NA
    )
  )


Fig$Size = as.factor(Fig$Size)


pdf('SuppFig1c.pdf')

P1 = ggscatter(Fig, x="LowStage_ratio_High", y="LowStage_ratio_Low",size ="Size",stroke = 0.6,color = "significance",palette = c('Sig'='red','NS'='skyblue'),alpha=0.6,label = "Cancer",label.select = labes1,repel=TRUE,font.label = c(7,"plain","black"))+theme_classic()+theme(axis.text.x=element_text(angle=0,size=7,vjust=0,hjust=1),legend.key.size=unit(0.05,"cm"),legend.position="top")+geom_abline(slope=1,intercept=0,linetype="dashed")


grid.newpage()
# Create layout : nrow = 2, ncol =2
pushViewport(viewport(layout = grid.layout(nrow =10, ncol = 10)))
# A helper function to define a region on the layout
define_region <- function(row, col){
  viewport(layout.pos.row = row, layout.pos.col = col)
} 

print(ggpar(P1,font.xtickslab =c("black",10),font.ytickslab =c("black",10),font.x = 10,font.y=10,font.legend=10,xlab="High EXTEND Low Stage Ratio",ylab="Low EXTEND Low Stage Ratio"),vp = define_region(row = 1:4, col = 1:4))


dev.off()


############ Supp Fig 1d : Hazard Ratio plot for survival analysis

HR_Results = read.table('SuppFig1d.txt',sep='\t',head=TRUE)


HR_Results$Significance = NULL
HR_Results$Significance = ifelse(HR_Results$p.value<= 0.054,"Sig",'NS')


cols2 <- c("NS" = "azure4",  "Sig" = "red")
HR_Results$HR = as.numeric(HR_Results$HR)

HR_Results$HR = log(HR_Results$HR)
HR_Results$LowerInterval = as.numeric(HR_Results$LowerInterval)
HR_Results$UpperInterval = as.numeric(HR_Results$UpperInterval)

HR_Results$LowerInterval = log(HR_Results$LowerInterval)
HR_Results$UpperInterval = log(HR_Results$UpperInterval)

### Sigificant cases (*) are marked manually on top of the figure
pdf('SuppFig1d.pdf')

  plot2 <- ggplot(data=HR_Results, aes(x=fct_reorder(Cancer,HR,.desc=F), y=HR,ymin=LowerInterval, ymax=UpperInterval),color=Significance) +
  geom_point(size=1,position=position_dodge(1),) +scale_fill_manual(values = cols2)+scale_color_manual(values = cols2)+
  geom_errorbar(width=0.2,size=0.1,position=position_dodge(1))+
        geom_hline(yintercept=0, lty=2,color="red",size=0.2) + xlab("Cancer") + ylab("Hazard Ratio") +
        theme_classic()+  # use a white background
		theme(axis.text.x=element_text(angle=45,size=7,vjust=1,hjust=1),legend.key.size= unit(0.3,"cm"),legend.key.width = unit(0.3,"cm"),legend.position = "bottom",axis.text.y = element_blank(), axis.title.y = element_blank())
  
grid.newpage()
# Create layout : nrow = 2, ncol =2
pushViewport(viewport(layout = grid.layout(nrow = 11, ncol = 11)))
# A helper function to define a region on the layout
define_region <- function(row, col){
  viewport(layout.pos.row = row, layout.pos.col = col)
} 
print(ggpar(plot2,font.xtickslab =c("black",7),font.ytickslab =c("black",10),font.x = 10,font.y=0),vp = define_region(row = 1:3, col = 1:7))  

dev.off()


###################################################################################################################################################
#############################################
#############################################               Supplementary Figure 2 : Whole genome doubling Ratios
#############################################
###################################################################################################################################################


AS_Scores_Pancan = read.table('SuppFig2.txt',sep='\t',head=T)


ExtendClassFracs = read.table('Fig1A_SignificanceInterval.txt',sep='\t',head=T)
ExtendClassFracs = ExtendClassFracs[,c(1:3,12:14)]
colnames(ExtendClassFracs) = c('Cancer','HighExtendSamples','LowExtendSamples','TotalSamples','HighExtendFrac','LowExtendFrac')

AS_Scores_Pancan = merge(AS_Scores_Pancan,ExtendClassFracs,by='Cancer')

AS_Scores_Pancan$significance = ifelse(AS_Scores_Pancan$p_value<0.02,'Sig','NS')

labes1 = AS_Scores_Pancan$Cancer[which(AS_Scores_Pancan$significance == "Sig")]


AS_Scores_Pancan$HighEXTEND_WGDFrac = AS_Scores_Pancan$doubled_High/AS_Scores_Pancan$HighExtendSamples
AS_Scores_Pancan$LowEXTEND_WGDFrac = AS_Scores_Pancan$doubled_Low/AS_Scores_Pancan$LowExtendSamples



AS_Scores_Pancan <- AS_Scores_Pancan %>%
  mutate(
    Size = case_when(
      HighExtendFrac < 0.25 ~ 1,
      HighExtendFrac >= 0.25 & HighExtendFrac < 0.5 ~ 2,
      HighExtendFrac >= 0.5 ~ 3,
      TRUE ~ NA_real_  # In case LowExtendFrac is NA
    )
  )


AS_Scores_Pancan$Size = as.factor(AS_Scores_Pancan$Size)


pdf('SuppFig2.pdf')

P1 = ggscatter(AS_Scores_Pancan, x="HighEXTEND_WGDFrac", y="LowEXTEND_WGDFrac",size ="Size",alpha=0.6,stroke = 0.6,color = "significance",palette = c('Sig'='red','NS'='skyblue'),label="Cancer",label.select= labes1,repel=T,font.label=c("black",5))+theme_classic()+theme(axis.text.x=element_text(angle=0,size=7,vjust=1,hjust=1,color="black"),legend.key.size=unit(0.05,"cm"),legend.position="top")+geom_abline(slope=1,intercept=0)


grid.newpage()
# Create layout : nrow = 2, ncol =2
pushViewport(viewport(layout = grid.layout(nrow =10, ncol = 10)))
# A helper function to define a region on the layout
define_region <- function(row, col){
  viewport(layout.pos.row = row, layout.pos.col = col)
} 

print(ggpar(P1,font.xtickslab =c("black",10),font.ytickslab =c("black",10),font.x = 10,font.y=10,font.legend=9,xlab="High EXTEND WGD Fraction",ylab="Low EXTEND WGD Fraction",ylim=c(0,0.9),xlim=c(0,0.9)),vp = define_region(row = 1:4, col = 1:4))

dev.off()


###################################################################################################################################################
#############################################
#############################################               Supplementary Figure 3 : LGG comparisons for telomerase high and low groups
#############################################
###################################################################################################################################################

combLGG = read.table('SuppFig3.txt',sep='\t',head=T)

my_comparisons1 = list(c("ATRX/DAXXalt","TERTexpr"))

pdf('SuppFig3.pdf')

P1 = ggplot(combLGG, aes(x=TMM , y=CNA_frac_altered,fill=TMM,color=TMM))+ geom_boxplot(size=0.1,outlier.size = 0.0,outlier.shape=NA,alpha=0.5)+geom_point(position = position_jitter(width = .15), size = 0.7) +scale_fill_manual(values=c("ATRX/DAXXalt"="darkseagreen2","TERTexpr"="lightpink2"))+scale_color_manual(values=c("ATRX/DAXXalt"="darkseagreen4","TERTexpr"="deeppink4"))+theme_classic()+theme(axis.text.x=element_text(angle=0,size=7,vjust=1,hjust=0.5),legend.position = "none",legend.key.size=unit(0.4,"cm"))+stat_compare_means(comparisons = my_comparisons1, method= "t.test",method.args = list(var.equal = TRUE,alternative ="two.sided"),size=3)

P2 = ggplot(combLGG, aes(x=TMM , y=LOH_frac_altered,fill=TMM,color=TMM))+ geom_boxplot(size=0.1,outlier.size = 0.0,outlier.shape=NA,alpha=0.5)+geom_point(position = position_jitter(width = .15), size = 0.7) +scale_fill_manual(values=c("ATRX/DAXXalt"="darkseagreen2","TERTexpr"="lightpink2"))+scale_color_manual(values=c("ATRX/DAXXalt"="darkseagreen4","TERTexpr"="deeppink4"))+theme_classic()+theme(axis.text.x=element_text(angle=0,size=7,vjust=1,hjust=0.5),legend.position = "none",legend.key.size=unit(0.4,"cm"))+stat_compare_means(comparisons = my_comparisons1, method= "t.test",method.args = list(alternative ="two.sided"),size=3)

P3 = ggplot(combLGG, aes(x=TMM , y=TLratio,fill=TMM,color=TMM))+ geom_boxplot(size=0.1,outlier.size = 0.0,outlier.shape=NA,alpha=0.5)+geom_point(position = position_jitter(width = .15), size = 0.7) +scale_fill_manual(values=c("ATRX/DAXXalt"="darkseagreen2","TERTexpr"="lightpink2"))+scale_color_manual(values=c("ATRX/DAXXalt"="darkseagreen4","TERTexpr"="deeppink4"))+theme_classic()+theme(axis.text.x=element_text(angle=0,size=7,vjust=1,hjust=0.5),legend.position = "none",legend.key.size=unit(0.4,"cm"))+stat_compare_means(comparisons = my_comparisons1, method= "t.test",method.args = list(var.equal = TRUE,alternative ="two.sided"),size=3)

P3 = ggplot(combLGG, aes(x=TMM , y=EXTENDScores,fill=TMM,color=TMM))+ geom_boxplot(size=0.1,outlier.size = 0.0,outlier.shape=NA,alpha=0.5)+geom_point(position = position_jitter(width = .15), size = 0.7) +scale_fill_manual(values=c("ATRX/DAXXalt"="darkseagreen2","TERTexpr"="lightpink2"))+scale_color_manual(values=c("ATRX/DAXXalt"="darkseagreen4","TERTexpr"="deeppink4"))+theme_classic()+theme(axis.text.x=element_text(angle=0,size=7,vjust=1,hjust=0.5),legend.position = "none",legend.key.size=unit(0.4,"cm"))+stat_compare_means(comparisons = my_comparisons1, method= "t.test",method.args = list(var.equal = TRUE,alternative ="two.sided"),size=3)

  
grid.newpage()
# Create layout : nrow = 2, ncol =2
pushViewport(viewport(layout = grid.layout(nrow = 4, ncol = 4)))
# A helper function to define a region on the layout
define_region <- function(row, col){
  viewport(layout.pos.row = row, layout.pos.col = col)
} 
print(ggpar(P1,font.xtickslab =8,font.ytickslab =8,font.x = 9,font.y=9,font.legend=8,ylab='CNA fractions altered'),vp = define_region(row = 1, col = 1))  
print(ggpar(P2,font.xtickslab =8,font.ytickslab =8,font.x = 9,font.y=9,font.legend=8,ylab='LOH fractions altered'),vp = define_region(row = 1, col = 2))  
print(ggpar(P3,font.xtickslab =8,font.ytickslab =8,font.x = 9,font.y=9,font.legend=8,ylab='TL Ratio'),vp = define_region(row = 2, col = 1))  
print(ggpar(P4,font.xtickslab =8,font.ytickslab =8,font.x = 9,font.y=9,font.legend=8,ylab='EXTEND Scores'),vp = define_region(row = 2, col = 2))  

dev.off()



###################################################################################################################################################
#############################################
############################################# Supplementary Figure 4 : Pathway Enrichment Analysis for low and high telomerase activity groups
#############################################
###################################################################################################################################################

Data = read.table('SuppFig4.txt',sep='\t',head=T)


Data$logFDR = -log10(Data$FDR)

HighEnr = Data[which(Data$EXTENDclass == 'High'),]
LowEnr = Data[which(Data$EXTENDclass == 'Low'),]


pdf('SuppFig4.pdf')
P1 = ggplot(HighEnr, aes(x=fct_reorder(GeneSetName,logFDR,.desc =FALSE), y=logFDR)) + geom_bar(stat = "identity",width=0.7,size=0.2, fill="pink",color="darkred")+theme_classic()+theme(legend.key.size=unit(0.1,"cm"))+coord_flip()
  
P2 = ggplot(LowEnr, aes(x=fct_reorder(GeneSetName,logFDR,.desc =FALSE), y=logFDR)) +  geom_bar(stat = "identity",width=0.7,size=0.2, fill="skyblue",color="darkblue")+theme_classic()+theme(legend.key.size=unit(0.1,"cm"))+coord_flip()


grid.newpage()
# Create layout : nrow = 2, ncol =2
pushViewport(viewport(layout = grid.layout(nrow = 7, ncol = 9)))
# A helper function to define a region on the layout
define_region <- function(row, col){
  viewport(layout.pos.row = row, layout.pos.col = col)
} 
print(ggpar(P1,font.xtickslab =8,font.ytickslab =7,font.x = 10,font.y=10,font.title=6,font.legend=6,xlab="Pathways",ylab="-log10(Pval)"),vp = define_region(row = 1:3, col = 1:6))  
print(ggpar(P2,font.xtickslab =8,font.ytickslab =7,font.x = 10,font.y=10,font.title=6,font.legend=6,xlab="Pathways",ylab="-log10(Pval)"),vp = define_region(row = 4:6, col = 1:5))  

dev.off()





###################################################################################################################################################
#############################################
#############################################               Supplementary Figure 5 : Overall Fusion frequencies across TCGA data
#############################################
###################################################################################################################################################


########################################## Supp Fig 5A

TCGA_TM_Fusions = read.table('SuppFig5A.txt',sep='\t',head=T)

### Significance levels added on top of barplot manually from the Figure data

pdf('SuppFig5A.pdf')
P1 = ggplot(TCGA_TM_Fusions, aes(x=Cancer , y=FusionFreq,fill=EXTENDclass)) + geom_bar(stat = "identity",size=0.3,width=0.8, position=position_dodge())+scale_fill_manual(values=c("High"="red","Low"="blue"))+theme_classic()+theme(axis.text.x=element_text(angle=35,size=7,vjust=1,hjust=1),legend.position = "top",legend.key.size=unit(0.4,"cm"))#+coord_flip()


  
grid.newpage()
# Create layout : nrow = 2, ncol =2
pushViewport(viewport(layout = grid.layout(nrow = 10, ncol = 10)))
# A helper function to define a region on the layout
define_region <- function(row, col){
  viewport(layout.pos.row = row, layout.pos.col = col)
} 
print(ggpar(P1,font.xtickslab =8,font.ytickslab =8,font.x = 9,font.y=9,font.legend=8,ylab ="Fusion Frequency"),vp = define_region(row = 1:3, col = 1:10))  

dev.off()

############################################### Supp Fig 5B


comball = read.table('SuppFig5B.txt',sep='\t',head=T)
pdf('SuppFig5B.pdf')
P1 = ggplot(RecFusionsFreqTotal, aes(x=Cancer , y=RecFusionsNo,fill=factor(EXTENDclass,levels = c("Low","High")))) + geom_bar(stat = "identity",size=0.1,width=0.7)+geom_text(aes(label=RecFusionsNo), position=position_dodge(width=0.9),hjust=0.35, vjust=-0.5,size=2)+scale_fill_manual(values=c("High"="hotpink","Low"="skyblue"))+theme_classic()+theme(axis.text.x=element_text(angle=35,size=7,vjust=1,hjust=1),legend.position = "top",legend.key.size=unit(0.4,"cm"))#+coord_flip()


  
grid.newpage()
# Create layout : nrow = 2, ncol =2
pushViewport(viewport(layout = grid.layout(nrow = 10, ncol = 11)))
# A helper function to define a region on the layout
define_region <- function(row, col){
  viewport(layout.pos.row = row, layout.pos.col = col)
} 
print(ggpar(P1,font.xtickslab =8,font.ytickslab =8,font.x = 9,font.y=9,font.legend=8,ylab ="Fusion Frequency"),vp = define_region(row = 1:3, col = 1:7))  

dev.off()

###################################################################################################################################################
#############################################
#############################################  Supplementary Figure 6 : Senescence and telomerase scores correlations in CCLE and Ackerman's Data
#############################################
###################################################################################################################################################

####### Supp 6A

subdata = read.table('SuppFig6A.txt',sep='\t',head=T)

pdf('SuppFig6A.pdf',onefile=FALSE)


P4 = ggscatter(subdata, x = "EXTENDScores", y = "SenescenceScores",color = "Mechanism",palette = c("TERT_Rearrangement"="maroon","TERT_High"="pink","MYCN_Amp"="red"),size = 2,add = "reg.line",  # Add regressin line
add.params = list(color = "black", fill = "grey",size=0.5), # Customize reg. line
conf.int = TRUE # Add confidence interval
)+stat_cor(method="spearman",size=3)



grid.newpage()
# Create layout : nrow = 3, ncol = 3
pushViewport(viewport(layout = grid.layout(nrow = 3, ncol = 3)))
# A helper function to define a region on the layout
define_region <- function(row, col){
  viewport(layout.pos.row = row, layout.pos.col = col)
} 

print(ggpar(P4,font.xtickslab =8,font.ytickslab =8,font.x = 8,font.y=8,xlab="EXTEND Scores",ylab ="Senescence Signature"),vp = define_region(row = 2, col = 2))

dev.off()



######## Supp 6B

combdata = read.table('SuppFig6B.txt',sep='\t',head=T)
combdata$colorTiss = paste(combdata$EXTENDclass,combdata$Tissues,sep='_')

source("https://gist.githubusercontent.com/benmarwick/2a1bb0133ff568cbe28d/raw/fb53bd97121f7f9ce947837ef1a4c65a73bffb3f/geom_flat_violin.R")

my_comparisons1 = list(c("High","Low"))



colpa = c("Low_A_Ganglia"="blue","Low_BiliaryTract"="blue","Low_Bone"="blue","Low_Breast"="blue","Low_Cervix"="blue","Low_CNS"="blue","Low_Endometrium"="blue","Low_Fibroblast" ="skyblue","Low_Haematpoetic&Lymphoid"="blue","Low_Kidney"="blue","Low_LargeIntestine"="blue","Low_Liver"="blue","Low_Lung"="blue","Low_Oesophagus"="blue","Low_Ovary"="blue","Low_Pancreas"="blue","Low_Pleura"="blue","Low_Prostate"="blue","Low_SalivaryGland"="blue","Low_Skin"="blue","Low_SmallIntestine"="blue","Low_SoftTissue"="blue","Low_Stomach"="blue","Low_Thyroid"="blue","Low_UpperAeroDigestiveTract"="blue","Low_UrinaryTract"="blue","High_A_Ganglia"="red","High_BiliaryTract"="red","High_Bone"="red","High_Breast"="red","High_Cervix"="red","High_CNS"="red","High_Endometrium"="red","High_Fibroblast" ="red","High_Haematpoetic&Lymphoid"="pink","High_Kidney"="red","High_LargeIntestine"="red","High_Liver"="red","High_Lung"="red","High_Oesophagus"="red","High_Ovary"="red","High_Pancreas"="red","High_Pleura"="red","High_Prostate"="red","High_SalivaryGland"="red","High_Skin"="red","High_SmallIntestine"="red","High_SoftTissue"="red","High_Stomach"="red","High_Thyroid"="red","High_UpperAeroDigestiveTract"="red","High_UrinaryTract"="red")


pdf('SuppFig6B.pdf')
  
P2 = ggscatter(combdata, y="SenescenceScores",x="EXTENDScores",size =1,color = "colorTiss",palette = colpa, add = "reg.line",
add.params = list(color = "blue", fill = "grey",size=0.2),conf.int = TRUE)+ stat_cor(method = "spearman",size=3)  +theme_classic()+theme(legend.key.size= unit(0.4,"cm"),legend.key.width = unit(0.2,"cm"),legend.position = "none")
   
 
  ##$#coord_flip()+
grid.newpage()
# Create layout : nrow = 2, ncol =2
pushViewport(viewport(layout = grid.layout(nrow = 8, ncol = 6)))
# A helper function to define a region on the layout
define_region <- function(row, col){
  viewport(layout.pos.row = row, layout.pos.col = col)
} 

print(ggpar(P2,font.xtickslab =c("black",10),font.ytickslab =c("black",10),font.x = 10,font.y=10,ylab="Senescence Scores",xlab="EXTEND Class"),vp = define_region(row = 1:2, col = 1:2))
dev.off()

###################################################################################################################################################
#############################################
#############################################               Supplementary Figure 7: Single Cell GBM and HNSC scores comparisons
#############################################
###################################################################################################################################################

############## Supp.Fig.7A-B 

comball1 = read.table('SuppFig7A_C.txt',sep='\t',head=T)
comball2 = read.table('SuppFig7B_D.txt',sep='\t',head=T)


cols2 = c("NonCycling"="skyblue","G2_M"="pink","G1_S" = "hotpink")
cols3 = c("NonCycling"="skyblue","G2_M"="red","G1_S" = "maroon")

my_comparisons2 = list(c("High","Low"))

my_comparisons1 = list(c("G1_S","G2_M"),c("G1_S","NonCycling"),c("G2_M","NonCycling"))

pdf('SuppFig7.pdf')
  
P1 = ggplot(data = comball1, 
         aes(x = EXTENDclass, y = Senescence, fill = "dodgerblue4")) +
  geom_flat_violin(position = position_nudge(x = .2, y = 0),fill="dodgerblue4",color="darkblue",alpha=0.5) +
   #geom_point(aes(y = Senescence),color="blue",fill = "white",shape = 21,stroke = .1,  
            # position = position_jitter(width = .15), size = .2)+  
  geom_boxplot(width = .15, outlier.shape = NA,fill="white",color="black",size=0.2) +
  guides(fill = FALSE) +
  guides(color = FALSE) +
  scale_color_brewer(palette = "Spectral") +
  scale_fill_brewer(palette = "Spectral") +  
  theme_classic()+theme(axis.text.x=element_text(size=8,vjust=1,hjust=1),legend.position="none")+stat_compare_means(comparisons = my_comparisons2, method= "t.test",size=3)# coord_flip() + # flip or not

  
P3 = ggplot(data = comball2, 
         aes(x = EXTENDclass, y = Senescence, fill = "dodgerblue4")) +
  geom_flat_violin(position = position_nudge(x = .2, y = 0),fill="dodgerblue4",color="darkblue",alpha=0.5) +
   #geom_point(aes(y = Senescence),color="blue",fill = "white",shape = 21,stroke = .1,  
            # position = position_jitter(width = .15), size = .2)+  
  geom_boxplot(width = .15, outlier.shape = NA,fill="white",color="black",size=0.2) +
  guides(fill = FALSE) +
  guides(color = FALSE) +
  scale_color_brewer(palette = "Spectral") +
  scale_fill_brewer(palette = "Spectral") +  
  theme_classic()+theme(axis.text.x=element_text(size=8,vjust=1,hjust=1),legend.position="none")+stat_compare_means(comparisons = my_comparisons2, method= "t.test",size=3)# coord_flip() + # flip or not




P2 = ggplot(data = comball1, 
         aes(x = Phase, y = Senescence, fill = Phase,color=Phase)) +
  geom_flat_violin(position = position_nudge(x = .2, y = 0)) +
  geom_point(aes(y = Senescence),
             position = position_jitter(width = .07), size = .1) +
  geom_boxplot(size=0.5, width = .15, outlier.shape = NA,fill="white",color="black") +
  #expand_limits(x = 5.25) +
  guides(fill = FALSE) +
  guides(color = FALSE) +
  scale_color_manual(values = cols3) +
  scale_fill_manual(values = cols3) +
  theme_classic()+theme(axis.text.x=element_text(size=8,vjust=1,hjust=1),legend.position="none")+stat_compare_means(comparisons = my_comparisons1, method= "t.test",size=2.5)# coord_flip() + # flip or not
 

P4 = ggplot(data = comball2, 
         aes(x = Phase, y = Senescence, fill = Phase,color=Phase)) +
  geom_flat_violin(position = position_nudge(x = .2, y = 0)) +
  geom_point(aes(y = Senescence),
             position = position_jitter(width = .07), size = .1) +
  geom_boxplot(size=0.5, width = .15, outlier.shape = NA,fill="white",color="black") +
  #expand_limits(x = 5.25) +
  guides(fill = FALSE) +
  guides(color = FALSE) +
  scale_color_manual(values = cols3) +
  scale_fill_manual(values = cols3) +
  theme_classic()+theme(axis.text.x=element_text(size=8,vjust=1,hjust=1),legend.position="none")+stat_compare_means(comparisons = my_comparisons1, method= "t.test",size=2.5)# coord_flip() + # flip or not
 


 
  ##$#coord_flip()+
grid.newpage()
# Create layout : nrow = 2, ncol =2
pushViewport(viewport(layout = grid.layout(nrow = 7, ncol = 6)))
# A helper function to define a region on the layout
define_region <- function(row, col){
  viewport(layout.pos.row = row, layout.pos.col = col)
} 


print(ggpar(P1,font.xtickslab =c("black",10),font.ytickslab =c("black",10),font.x = 10,font.y=10,ylab="Senescence Scores",xlab="EXTEND Class"),vp = define_region(row = 1:2, col = 1:2))

print(ggpar(P2,font.xtickslab =c("black",10),font.ytickslab =c("black",10),font.x = 10,font.y=10,ylab="Senescence Scores",xlab="EXTEND Class"),vp = define_region(row = 1:2, col = 3:4))

print(ggpar(P3,font.xtickslab =c("black",10),font.ytickslab =c("black",10),font.x = 10,font.y=10),vp = define_region(row = 3:4, col = 1:2))

print(ggpar(P4,font.xtickslab =c("black",10),font.ytickslab =c("black",10),font.x = 10,font.y=10),vp = define_region(row = 3:4, col = 3:4))

dev.off()



###################################################################################################################################################
#############################################
#############################################               Supplementary Figure 8 : Immune cell distributions 
#############################################
###################################################################################################################################################


dataImm = read.table('cibersort_EXTENDClass_TCGA_052023.txt',sep='\t',head=T)


colours = c('Low'='blue','High'='red')

dataImm = dataImm[order(dataImm$Cancer),]

can = unique(dataImm$Cancer)
my_comparisons1 = list(c("Low","High"))

subM0_Can = dataImm[which(dataImm$Cancer %in% c('ACC','BRCA','GBM','LUAD','OV','STAD')),]

my_comparisons1 = list(c("Low","High"))

pdf('Supp.Fig8a_Macrophage_M0.pdf')

subACC = dataImm[which(dataImm$Cancer == 'ACC'),]
subBRCA = dataImm[which(dataImm$Cancer == 'BRCA'),]
subGBM = dataImm[which(dataImm$Cancer == 'GBM'),]
subLUAD = dataImm[which(dataImm$Cancer == 'LUAD'),]
subOV = dataImm[which(dataImm$Cancer == 'OV'),]
subSTAD = dataImm[which(dataImm$Cancer == 'STAD'),]

P1 = ggplot(subACC, aes(x=EXTENDclass , y=Macrophages.M0,fill=EXTENDclass,color=EXTENDclass))+ geom_boxplot(size=0.1,outlier.size = 0.0,outlier.shape=NA,alpha=0.5,width=0.5)+geom_point(position = position_jitter(width = .1), size = 0.3) +scale_fill_manual(values=c("Low"="skyblue","High"="pink"))+scale_color_manual(values=c("Low"="blue","High"="darkred"))+theme_classic()+theme(axis.text.x=element_text(angle=0,size=7,vjust=1,hjust=0.5),legend.position = "none",legend.key.size=unit(0.4,"cm"))+stat_compare_means(comparisons = my_comparisons1, method= "t.test",method.args = list(var.equal = TRUE,alternative ="two.sided"),size=3)

P2 = ggplot(subBRCA, aes(x=EXTENDclass , y=Macrophages.M0,fill=EXTENDclass,color=EXTENDclass))+ geom_boxplot(size=0.1,outlier.size = 0.0,outlier.shape=NA,alpha=0.5,width=0.5)+geom_point(position = position_jitter(width = .1), size = 0.3) +scale_fill_manual(values=c("Low"="skyblue","High"="pink"))+scale_color_manual(values=c("Low"="blue","High"="darkred"))+theme_classic()+theme(axis.text.x=element_text(angle=0,size=7,vjust=1,hjust=0.5),legend.position = "none",legend.key.size=unit(0.4,"cm"))+stat_compare_means(comparisons = my_comparisons1, method= "t.test",method.args = list(var.equal = TRUE,alternative ="two.sided"),size=3)

P3 = ggplot(subGBM, aes(x=EXTENDclass , y=Macrophages.M0,fill=EXTENDclass,color=EXTENDclass))+ geom_boxplot(size=0.1,outlier.size = 0.0,outlier.shape=NA,alpha=0.5,width=0.5)+geom_point(position = position_jitter(width = .1), size = 0.3) +scale_fill_manual(values=c("Low"="skyblue","High"="pink"))+scale_color_manual(values=c("Low"="blue","High"="darkred"))+theme_classic()+theme(axis.text.x=element_text(angle=0,size=7,vjust=1,hjust=0.5),legend.position = "none",legend.key.size=unit(0.4,"cm"))+stat_compare_means(comparisons = my_comparisons1, method= "t.test",method.args = list(var.equal = TRUE,alternative ="two.sided"),size=3)

P4 = ggplot(subLUAD, aes(x=EXTENDclass , y=Macrophages.M0,fill=EXTENDclass,color=EXTENDclass))+ geom_boxplot(size=0.1,outlier.size = 0.0,outlier.shape=NA,alpha=0.5,width=0.5)+geom_point(position = position_jitter(width = .1), size = 0.3) +scale_fill_manual(values=c("Low"="skyblue","High"="pink"))+scale_color_manual(values=c("Low"="blue","High"="darkred"))+theme_classic()+theme(axis.text.x=element_text(angle=0,size=7,vjust=1,hjust=0.5),legend.position = "none",legend.key.size=unit(0.4,"cm"))+stat_compare_means(comparisons = my_comparisons1, method= "t.test",method.args = list(var.equal = TRUE,alternative ="two.sided"),size=3)

P5 = ggplot(subOV, aes(x=EXTENDclass , y=Macrophages.M0,fill=EXTENDclass,color=EXTENDclass))+ geom_boxplot(size=0.1,outlier.size = 0.0,outlier.shape=NA,alpha=0.5,width=0.5)+geom_point(position = position_jitter(width = .1), size = 0.3) +scale_fill_manual(values=c("Low"="skyblue","High"="pink"))+scale_color_manual(values=c("Low"="blue","High"="darkred"))+theme_classic()+theme(axis.text.x=element_text(angle=0,size=7,vjust=1,hjust=0.5),legend.position = "none",legend.key.size=unit(0.4,"cm"))+stat_compare_means(comparisons = my_comparisons1, method= "t.test",method.args = list(var.equal = TRUE,alternative ="two.sided"),size=3)

P6 = ggplot(subSTAD, aes(x=EXTENDclass , y=Macrophages.M0,fill=EXTENDclass,color=EXTENDclass))+ geom_boxplot(size=0.1,outlier.size = 0.0,outlier.shape=NA,alpha=0.5,width=0.5)+geom_point(position = position_jitter(width = .1), size = 0.3) +scale_fill_manual(values=c("Low"="skyblue","High"="pink"))+scale_color_manual(values=c("Low"="blue","High"="darkred"))+theme_classic()+theme(axis.text.x=element_text(angle=0,size=7,vjust=1,hjust=0.5),legend.position = "none",legend.key.size=unit(0.4,"cm"))+stat_compare_means(comparisons = my_comparisons1, method= "t.test",method.args = list(var.equal = TRUE,alternative ="two.sided"),size=3)


grid.newpage()
# Create layout : nrow = 2, ncol =2
pushViewport(viewport(layout = grid.layout(nrow =10, ncol = 8)))
# A helper function to define a region on the layout
define_region <- function(row, col){
  viewport(layout.pos.row = row, layout.pos.col = col)
} 

print(ggpar(P1,font.xtickslab =c("black",9),font.ytickslab =c("black",10),font.x = 9,font.y=9,font.legend=7,ylab="M0 Macrophages",xlab='ACC'),vp = define_region(row = 1:2, col = 1:2))
print(ggpar(P2,font.xtickslab =c("black",9),font.ytickslab =c("black",10),font.x = 9,font.y=9,font.legend=7,ylab="M0 Macrophages",xlab='BRCA'),vp = define_region(row = 1:2, col = 3:4))
print(ggpar(P3,font.xtickslab =c("black",9),font.ytickslab =c("black",10),font.x = 9,font.y=9,font.legend=7,ylab="M0 Macrophages",xlab='GBM'),vp = define_region(row = 1:2, col = 5:6))
print(ggpar(P4,font.xtickslab =c("black",9),font.ytickslab =c("black",10),font.x = 9,font.y=9,font.legend=7,ylab="M0 Macrophages",xlab='LUAD'),vp = define_region(row = 3:4, col = 1:2))
print(ggpar(P5,font.xtickslab =c("black",9),font.ytickslab =c("black",10),font.x = 9,font.y=9,font.legend=7,ylab="M0 Macrophages",xlab='OV'),vp = define_region(row = 3:4, col = 3:4))
print(ggpar(P6,font.xtickslab =c("black",9),font.ytickslab =c("black",10),font.x = 9,font.y=9,font.legend=7,ylab="M0 Macrophages",xlab='STAD'),vp = define_region(row = 3:4, col = 5:6))


dev.off()



my_comparisons1 = list(c("Low","High"))

pdf('SuppFig8b_Macrophage_M1.pdf')

subBRCA = dataImm[which(dataImm$Cancer == 'BRCA'),]
subKIRP = dataImm[which(dataImm$Cancer == 'KIRP'),]
subLUAD = dataImm[which(dataImm$Cancer == 'LUAD'),]
subPCPG = dataImm[which(dataImm$Cancer == 'PCPG'),]
subSTAD = dataImm[which(dataImm$Cancer == 'STAD'),]
subTGCT = dataImm[which(dataImm$Cancer == 'TGCT'),]
subTHYM = dataImm[which(dataImm$Cancer == 'THYM'),]


P1 = ggplot(subBRCA, aes(x=EXTENDclass , y=Macrophages.M1,fill=EXTENDclass,color=EXTENDclass))+ geom_boxplot(size=0.1,outlier.size = 0.0,outlier.shape=NA,alpha=0.5,width=0.5)+geom_point(position = position_jitter(width = .1), size = 0.3) +scale_fill_manual(values=c("Low"="skyblue","High"="pink"))+scale_color_manual(values=c("Low"="blue","High"="darkred"))+theme_classic()+theme(axis.text.x=element_text(angle=0,size=7,vjust=1,hjust=0.5),legend.position = "none",legend.key.size=unit(0.4,"cm"))+stat_compare_means(comparisons = my_comparisons1, method= "t.test",method.args = list(var.equal = TRUE,alternative ="two.sided"),size=3)

P2 = ggplot(subKIRP, aes(x=EXTENDclass , y=Macrophages.M1,fill=EXTENDclass,color=EXTENDclass))+ geom_boxplot(size=0.1,outlier.size = 0.0,outlier.shape=NA,alpha=0.5,width=0.5)+geom_point(position = position_jitter(width = .1), size = 0.3) +scale_fill_manual(values=c("Low"="skyblue","High"="pink"))+scale_color_manual(values=c("Low"="blue","High"="darkred"))+theme_classic()+theme(axis.text.x=element_text(angle=0,size=7,vjust=1,hjust=0.5),legend.position = "none",legend.key.size=unit(0.4,"cm"))+stat_compare_means(comparisons = my_comparisons1, method= "t.test",method.args = list(var.equal = TRUE,alternative ="two.sided"),size=3)

P3 = ggplot(subLUAD, aes(x=EXTENDclass , y=Macrophages.M1,fill=EXTENDclass,color=EXTENDclass))+ geom_boxplot(size=0.1,outlier.size = 0.0,outlier.shape=NA,alpha=0.5,width=0.5)+geom_point(position = position_jitter(width = .1), size = 0.3) +scale_fill_manual(values=c("Low"="skyblue","High"="pink"))+scale_color_manual(values=c("Low"="blue","High"="darkred"))+theme_classic()+theme(axis.text.x=element_text(angle=0,size=7,vjust=1,hjust=0.5),legend.position = "none",legend.key.size=unit(0.4,"cm"))+stat_compare_means(comparisons = my_comparisons1, method= "t.test",method.args = list(var.equal = TRUE,alternative ="two.sided"),size=3)

P4 = ggplot(subPCPG, aes(x=EXTENDclass , y=Macrophages.M1,fill=EXTENDclass,color=EXTENDclass))+ geom_boxplot(size=0.1,outlier.size = 0.0,outlier.shape=NA,alpha=0.5,width=0.5)+geom_point(position = position_jitter(width = .1), size = 0.3) +scale_fill_manual(values=c("Low"="skyblue","High"="pink"))+scale_color_manual(values=c("Low"="blue","High"="darkred"))+theme_classic()+theme(axis.text.x=element_text(angle=0,size=7,vjust=1,hjust=0.5),legend.position = "none",legend.key.size=unit(0.4,"cm"))+stat_compare_means(comparisons = my_comparisons1, method= "t.test",method.args = list(var.equal = TRUE,alternative ="two.sided"),size=3)

P5 = ggplot(subSTAD, aes(x=EXTENDclass , y=Macrophages.M1,fill=EXTENDclass,color=EXTENDclass))+ geom_boxplot(size=0.1,outlier.size = 0.0,outlier.shape=NA,alpha=0.5,width=0.5)+geom_point(position = position_jitter(width = .1), size = 0.3) +scale_fill_manual(values=c("Low"="skyblue","High"="pink"))+scale_color_manual(values=c("Low"="blue","High"="darkred"))+theme_classic()+theme(axis.text.x=element_text(angle=0,size=7,vjust=1,hjust=0.5),legend.position = "none",legend.key.size=unit(0.4,"cm"))+stat_compare_means(comparisons = my_comparisons1, method= "t.test",method.args = list(var.equal = TRUE,alternative ="two.sided"),size=3)

P6 = ggplot(subTGCT, aes(x=EXTENDclass , y=Macrophages.M1,fill=EXTENDclass,color=EXTENDclass))+ geom_boxplot(size=0.1,outlier.size = 0.0,outlier.shape=NA,alpha=0.5,width=0.5)+geom_point(position = position_jitter(width = .1), size = 0.3) +scale_fill_manual(values=c("Low"="skyblue","High"="pink"))+scale_color_manual(values=c("Low"="blue","High"="darkred"))+theme_classic()+theme(axis.text.x=element_text(angle=0,size=7,vjust=1,hjust=0.5),legend.position = "none",legend.key.size=unit(0.4,"cm"))+stat_compare_means(comparisons = my_comparisons1, method= "t.test",method.args = list(var.equal = TRUE,alternative ="two.sided"),size=3)

P7 = ggplot(subTHYM, aes(x=EXTENDclass , y=Macrophages.M1,fill=EXTENDclass,color=EXTENDclass))+ geom_boxplot(size=0.1,outlier.size = 0.0,outlier.shape=NA,alpha=0.5,width=0.5)+geom_point(position = position_jitter(width = .1), size = 0.3) +scale_fill_manual(values=c("Low"="skyblue","High"="pink"))+scale_color_manual(values=c("Low"="blue","High"="darkred"))+theme_classic()+theme(axis.text.x=element_text(angle=0,size=7,vjust=1,hjust=0.5),legend.position = "none",legend.key.size=unit(0.4,"cm"))+stat_compare_means(comparisons = my_comparisons1, method= "t.test",method.args = list(var.equal = TRUE,alternative ="two.sided"),size=3)

grid.newpage()
# Create layout : nrow = 2, ncol =2
pushViewport(viewport(layout = grid.layout(nrow =10, ncol = 8)))
# A helper function to define a region on the layout
define_region <- function(row, col){
  viewport(layout.pos.row = row, layout.pos.col = col)
} 

print(ggpar(P1,font.xtickslab =c("black",9),font.ytickslab =c("black",10),font.x = 9,font.y=9,font.legend=7,ylab="M1 Macrophages",xlab='BRCA'),vp = define_region(row = 1:2, col = 1:2))
print(ggpar(P2,font.xtickslab =c("black",9),font.ytickslab =c("black",10),font.x = 9,font.y=9,font.legend=7,ylab="M1 Macrophages",xlab='KIRP'),vp = define_region(row = 1:2, col = 3:4))
print(ggpar(P3,font.xtickslab =c("black",9),font.ytickslab =c("black",10),font.x = 9,font.y=9,font.legend=7,ylab="M1 Macrophages",xlab='LUAD'),vp = define_region(row = 1:2, col = 5:6))
print(ggpar(P4,font.xtickslab =c("black",9),font.ytickslab =c("black",10),font.x = 9,font.y=9,font.legend=7,ylab="M1 Macrophages",xlab='PCPG'),vp = define_region(row = 3:4, col = 1:2))
print(ggpar(P5,font.xtickslab =c("black",9),font.ytickslab =c("black",10),font.x = 9,font.y=9,font.legend=7,ylab="M1 Macrophages",xlab='STAD'),vp = define_region(row = 3:4, col = 3:4))
print(ggpar(P6,font.xtickslab =c("black",9),font.ytickslab =c("black",10),font.x = 9,font.y=9,font.legend=7,ylab="M1 Macrophages",xlab='TGCT'),vp = define_region(row = 3:4, col = 5:6))
print(ggpar(P7,font.xtickslab =c("black",9),font.ytickslab =c("black",10),font.x = 9,font.y=9,font.legend=7,ylab="M1 Macrophages",xlab='THYM'),vp = define_region(row = 3:4, col = 7:8))

dev.off()




###################################################################################################################################################
#############################################
#############################################               Supplementary Figure 9 : Immune cell abundances 
#############################################
###################################################################################################################################################


dataImm = read.table('cibersort_EXTENDClass_TCGA_052023.txt',sep='\t',head=T)


colours = c('Low'='blue','High'='red')

dataImm = dataImm[order(dataImm$Cancer),]

can = unique(dataImm$Cancer)
my_comparisons1 = list(c("Low","High"))


pdf('SuppFig9_Macrophage_M2.pdf')

subCOAD = dataImm[which(dataImm$Cancer == 'COAD'),]
subKIRP = dataImm[which(dataImm$Cancer == 'KIRP'),]
subLUAD = dataImm[which(dataImm$Cancer == 'LUAD'),]
subOV = dataImm[which(dataImm$Cancer == 'OV'),]
subREAD = dataImm[which(dataImm$Cancer == 'READ'),]
subTGCT = dataImm[which(dataImm$Cancer == 'TGCT'),]
subTHYM = dataImm[which(dataImm$Cancer == 'THYM'),]


P1 = ggplot(subCOAD, aes(x=EXTENDclass , y=Macrophages.M2,fill=EXTENDclass,color=EXTENDclass))+ geom_boxplot(size=0.1,outlier.size = 0.0,outlier.shape=NA,alpha=0.5,width=0.5)+geom_point(position = position_jitter(width = .1), size = 0.3) +scale_fill_manual(values=c("Low"="skyblue","High"="pink"))+scale_color_manual(values=c("Low"="blue","High"="darkred"))+theme_classic()+theme(axis.text.x=element_text(angle=0,size=7,vjust=1,hjust=0.5),legend.position = "none",legend.key.size=unit(0.4,"cm"))+stat_compare_means(comparisons = my_comparisons1, method= "t.test",method.args = list(var.equal = TRUE,alternative ="two.sided"),size=3)

P2 = ggplot(subKIRP, aes(x=EXTENDclass , y=Macrophages.M2,fill=EXTENDclass,color=EXTENDclass))+ geom_boxplot(size=0.1,outlier.size = 0.0,outlier.shape=NA,alpha=0.5,width=0.5)+geom_point(position = position_jitter(width = .1), size = 0.3) +scale_fill_manual(values=c("Low"="skyblue","High"="pink"))+scale_color_manual(values=c("Low"="blue","High"="darkred"))+theme_classic()+theme(axis.text.x=element_text(angle=0,size=7,vjust=1,hjust=0.5),legend.position = "none",legend.key.size=unit(0.4,"cm"))+stat_compare_means(comparisons = my_comparisons1, method= "t.test",method.args = list(var.equal = TRUE,alternative ="two.sided"),size=3)

P3 = ggplot(subLUAD, aes(x=EXTENDclass , y=Macrophages.M2,fill=EXTENDclass,color=EXTENDclass))+ geom_boxplot(size=0.1,outlier.size = 0.0,outlier.shape=NA,alpha=0.5,width=0.5)+geom_point(position = position_jitter(width = .1), size = 0.3) +scale_fill_manual(values=c("Low"="skyblue","High"="pink"))+scale_color_manual(values=c("Low"="blue","High"="darkred"))+theme_classic()+theme(axis.text.x=element_text(angle=0,size=7,vjust=1,hjust=0.5),legend.position = "none",legend.key.size=unit(0.4,"cm"))+stat_compare_means(comparisons = my_comparisons1, method= "t.test",method.args = list(var.equal = TRUE,alternative ="two.sided"),size=3)

P4 = ggplot(subOV, aes(x=EXTENDclass , y=Macrophages.M2,fill=EXTENDclass,color=EXTENDclass))+ geom_boxplot(size=0.1,outlier.size = 0.0,outlier.shape=NA,alpha=0.5,width=0.5)+geom_point(position = position_jitter(width = .1), size = 0.3) +scale_fill_manual(values=c("Low"="skyblue","High"="pink"))+scale_color_manual(values=c("Low"="blue","High"="darkred"))+theme_classic()+theme(axis.text.x=element_text(angle=0,size=7,vjust=1,hjust=0.5),legend.position = "none",legend.key.size=unit(0.4,"cm"))+stat_compare_means(comparisons = my_comparisons1, method= "t.test",method.args = list(var.equal = TRUE,alternative ="two.sided"),size=3)

P5 = ggplot(subREAD, aes(x=EXTENDclass , y=Macrophages.M2,fill=EXTENDclass,color=EXTENDclass))+ geom_boxplot(size=0.1,outlier.size = 0.0,outlier.shape=NA,alpha=0.5,width=0.5)+geom_point(position = position_jitter(width = .1), size = 0.3) +scale_fill_manual(values=c("Low"="skyblue","High"="pink"))+scale_color_manual(values=c("Low"="blue","High"="darkred"))+theme_classic()+theme(axis.text.x=element_text(angle=0,size=7,vjust=1,hjust=0.5),legend.position = "none",legend.key.size=unit(0.4,"cm"))+stat_compare_means(comparisons = my_comparisons1, method= "t.test",method.args = list(var.equal = TRUE,alternative ="two.sided"),size=3)

P6 = ggplot(subTGCT, aes(x=EXTENDclass , y=Macrophages.M2,fill=EXTENDclass,color=EXTENDclass))+ geom_boxplot(size=0.1,outlier.size = 0.0,outlier.shape=NA,alpha=0.5,width=0.5)+geom_point(position = position_jitter(width = .1), size = 0.3) +scale_fill_manual(values=c("Low"="skyblue","High"="pink"))+scale_color_manual(values=c("Low"="blue","High"="darkred"))+theme_classic()+theme(axis.text.x=element_text(angle=0,size=7,vjust=1,hjust=0.5),legend.position = "none",legend.key.size=unit(0.4,"cm"))+stat_compare_means(comparisons = my_comparisons1, method= "t.test",method.args = list(var.equal = TRUE,alternative ="two.sided"),size=3)

P7 = ggplot(subTHYM, aes(x=EXTENDclass , y=Macrophages.M2,fill=EXTENDclass,color=EXTENDclass))+ geom_boxplot(size=0.1,outlier.size = 0.0,outlier.shape=NA,alpha=0.5,width=0.5)+geom_point(position = position_jitter(width = .1), size = 0.3) +scale_fill_manual(values=c("Low"="skyblue","High"="pink"))+scale_color_manual(values=c("Low"="blue","High"="darkred"))+theme_classic()+theme(axis.text.x=element_text(angle=0,size=7,vjust=1,hjust=0.5),legend.position = "none",legend.key.size=unit(0.4,"cm"))+stat_compare_means(comparisons = my_comparisons1, method= "t.test",method.args = list(var.equal = TRUE,alternative ="two.sided"),size=3)

grid.newpage()
# Create layout : nrow = 2, ncol =2
pushViewport(viewport(layout = grid.layout(nrow =10, ncol = 8)))
# A helper function to define a region on the layout
define_region <- function(row, col){
  viewport(layout.pos.row = row, layout.pos.col = col)
} 

print(ggpar(P1,font.xtickslab =c("black",9),font.ytickslab =c("black",10),font.x = 9,font.y=9,font.legend=7,ylab="M2 Macrophages",xlab='COAD'),vp = define_region(row = 1:2, col = 1:2))
print(ggpar(P2,font.xtickslab =c("black",9),font.ytickslab =c("black",10),font.x = 9,font.y=9,font.legend=7,ylab="M2 Macrophages",xlab='KIRP'),vp = define_region(row = 1:2, col = 3:4))
print(ggpar(P3,font.xtickslab =c("black",9),font.ytickslab =c("black",10),font.x = 9,font.y=9,font.legend=7,ylab="M2 Macrophages",xlab='LUAD'),vp = define_region(row = 1:2, col = 5:6))
print(ggpar(P4,font.xtickslab =c("black",9),font.ytickslab =c("black",10),font.x = 9,font.y=9,font.legend=7,ylab="M2 Macrophages",xlab='OV'),vp = define_region(row = 3:4, col = 1:2))
print(ggpar(P5,font.xtickslab =c("black",9),font.ytickslab =c("black",10),font.x = 9,font.y=9,font.legend=7,ylab="M2 Macrophages",xlab='READ'),vp = define_region(row = 3:4, col = 3:4))
print(ggpar(P6,font.xtickslab =c("black",9),font.ytickslab =c("black",10),font.x = 9,font.y=9,font.legend=7,ylab="M2 Macrophages",xlab='TGCT'),vp = define_region(row = 3:4, col = 5:6))
print(ggpar(P7,font.xtickslab =c("black",9),font.ytickslab =c("black",10),font.x = 9,font.y=9,font.legend=7,ylab="M2 Macrophages",xlab='THYM'),vp = define_region(row = 3:4, col = 7:8))

dev.off()



my_comparisons1 = list(c("Low","High"))

pdf('SuppFig9_Monocytes.pdf')

subBRCA = dataImm[which(dataImm$Cancer == 'BRCA'),]
subGBM = dataImm[which(dataImm$Cancer == 'GBM'),]
subLAML = dataImm[which(dataImm$Cancer == 'LAML'),]
subPAAD = dataImm[which(dataImm$Cancer == 'PAAD'),]



P1 = ggplot(subBRCA, aes(x=EXTENDclass , y=Monocytes,fill=EXTENDclass,color=EXTENDclass))+ geom_boxplot(size=0.1,outlier.size = 0.0,outlier.shape=NA,alpha=0.5,width=0.5)+geom_point(position = position_jitter(width = .1), size = 0.3) +scale_fill_manual(values=c("Low"="skyblue","High"="pink"))+scale_color_manual(values=c("Low"="blue","High"="darkred"))+theme_classic()+theme(axis.text.x=element_text(angle=0,size=7,vjust=1,hjust=0.5),legend.position = "none",legend.key.size=unit(0.4,"cm"))+stat_compare_means(comparisons = my_comparisons1, method= "t.test",method.args = list(var.equal = TRUE,alternative ="two.sided"),size=3)

P2 = ggplot(subGBM, aes(x=EXTENDclass , y=Monocytes,fill=EXTENDclass,color=EXTENDclass))+ geom_boxplot(size=0.1,outlier.size = 0.0,outlier.shape=NA,alpha=0.5,width=0.5)+geom_point(position = position_jitter(width = .1), size = 0.3) +scale_fill_manual(values=c("Low"="skyblue","High"="pink"))+scale_color_manual(values=c("Low"="blue","High"="darkred"))+theme_classic()+theme(axis.text.x=element_text(angle=0,size=7,vjust=1,hjust=0.5),legend.position = "none",legend.key.size=unit(0.4,"cm"))+stat_compare_means(comparisons = my_comparisons1, method= "t.test",method.args = list(var.equal = TRUE,alternative ="two.sided"),size=3)

P3 = ggplot(subLAML, aes(x=EXTENDclass , y=Monocytes,fill=EXTENDclass,color=EXTENDclass))+ geom_boxplot(size=0.1,outlier.size = 0.0,outlier.shape=NA,alpha=0.5,width=0.5)+geom_point(position = position_jitter(width = .1), size = 0.3) +scale_fill_manual(values=c("Low"="skyblue","High"="pink"))+scale_color_manual(values=c("Low"="blue","High"="darkred"))+theme_classic()+theme(axis.text.x=element_text(angle=0,size=7,vjust=1,hjust=0.5),legend.position = "none",legend.key.size=unit(0.4,"cm"))+stat_compare_means(comparisons = my_comparisons1, method= "t.test",method.args = list(var.equal = TRUE,alternative ="two.sided"),size=3)

P4 = ggplot(subPAAD, aes(x=EXTENDclass , y=Monocytes,fill=EXTENDclass,color=EXTENDclass))+ geom_boxplot(size=0.1,outlier.size = 0.0,outlier.shape=NA,alpha=0.5,width=0.5)+geom_point(position = position_jitter(width = .1), size = 0.3) +scale_fill_manual(values=c("Low"="skyblue","High"="pink"))+scale_color_manual(values=c("Low"="blue","High"="darkred"))+theme_classic()+theme(axis.text.x=element_text(angle=0,size=7,vjust=1,hjust=0.5),legend.position = "none",legend.key.size=unit(0.4,"cm"))+stat_compare_means(comparisons = my_comparisons1, method= "t.test",method.args = list(var.equal = TRUE,alternative ="two.sided"),size=3)


grid.newpage()
# Create layout : nrow = 2, ncol =2
pushViewport(viewport(layout = grid.layout(nrow =10, ncol = 8)))
# A helper function to define a region on the layout
define_region <- function(row, col){
  viewport(layout.pos.row = row, layout.pos.col = col)
} 

print(ggpar(P1,font.xtickslab =c("black",9),font.ytickslab =c("black",10),font.x = 9,font.y=9,font.legend=7,ylab="Monocytes",xlab='BRCA'),vp = define_region(row = 1:2, col = 1:2))
print(ggpar(P2,font.xtickslab =c("black",9),font.ytickslab =c("black",10),font.x = 9,font.y=9,font.legend=7,ylab="Monocytes",xlab='GBM'),vp = define_region(row = 1:2, col = 3:4))
print(ggpar(P3,font.xtickslab =c("black",9),font.ytickslab =c("black",10),font.x = 9,font.y=9,font.legend=7,ylab="Monocytes",xlab='LAML'),vp = define_region(row = 1:2, col = 5:6))
print(ggpar(P4,font.xtickslab =c("black",9),font.ytickslab =c("black",10),font.x = 9,font.y=9,font.legend=7,ylab="Monocytes",xlab='PAAD'),vp = define_region(row = 3:4, col = 1:2))
dev.off()



###################################################################################################################################################
#############################################
#############################################               Supplementary Figure 10: Immune cell populations in Single cell data sets
#############################################
###################################################################################################################################################
comball = read.table('SuppFigure10.txt',sep='\t',head=T)


comball$curated_anno <- factor(comball$curated_anno, levels = c("Mono/Macro", "DC","NK","B","Treg","CD8T","CD4Tconv","Tprolif"))
tiss = unique(comball$Tissue)

pdf('SuppFig10.pdf')   #### Plotted separately and then combined in illustrator
for(i in 1:length(tiss)){

subsetdata = comball[which(comball$Tissue %in% tiss[i]),]

P2 = ggplot(subsetdata, aes(x=curated_anno , y=EXTENDscores))+geom_violin(size=0.1,draw_quantiles=c(0.5),alpha=0.5,fill='red',color='darkred')+theme_classic()+theme(axis.text.x=element_text(angle=35,size=7,vjust=1,hjust=1),legend.position = "none",legend.key.size=unit(0.4,"cm"))+stat_compare_means(comparisons = my_comparisons2, method= "t.test",method.args = list(var.equal = TRUE,alternative ="two.sided"),size=3)

  
grid.newpage()
# Create layout : nrow = 2, ncol =2
pushViewport(viewport(layout = grid.layout(nrow = 6, ncol = 7)))
# A helper function to define a region on the layout
define_region <- function(row, col){
  viewport(layout.pos.row = row, layout.pos.col = col)
} 
print(ggpar(P2,font.xtickslab =8,font.ytickslab =8,font.x = 9,font.y=9,font.legend=8,ylab='EXTEND Scores',xlab=tiss[i]),vp = define_region(row = 2:3, col = 2:6))  

}
dev.off()

