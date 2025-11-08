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

##setwd('D:/TTUHSC/Project_Extend_Classification/Analysis_WithTumorCasesOnly/Manuscript/Version_102024/Version_102024/Computational_StructuralBiology/Revisions/AllFigures_andData/Revisions_FinalVersion/MainFigures_Data_Code')
################################################################################################
############################ MAIN FIGURES #################################


###			Data for the main figures is uploaded in Github repository 

###################################################################################################################################################
#############################################
#############################################                   Figure 1 : Stratification and Genome Instability
#############################################
###################################################################################################################################################



############# 1A ###################################

comball = read.table('Fig1A.txt',sep='\t',head=T)


# Step 1: Pivot to wide format for each cancer type
comball_wide <- comball %>%
  select(Cancer, EXTENDclass, Frequency) %>%
  pivot_wider(names_from = EXTENDclass, values_from = Frequency)

# Step 2: Total High/Low across all cancers
totals <- comball %>%
  group_by(EXTENDclass) %>%
  summarise(Total = sum(Frequency)) %>%
  pivot_wider(names_from = EXTENDclass, values_from = Total)

total_high <- totals$High
total_low  <- totals$Low

# Step 3: For each cancer type, build 2x2 table: [Cancer vs All Others]
results <- comball_wide %>%
  rowwise() %>%
  mutate(
    high_cancer = High,
    low_cancer = Low,
    high_others = total_high - High,
    low_others = total_low - Low,
    
    # Build 2x2 matrix
    test = list(fisher.test(matrix(c(high_cancer, low_cancer, high_others, low_others), nrow = 2))),
    
    odds_ratio = test$estimate,
    ci_low = test$conf.int[1],
    ci_high = test$conf.int[2],
    p_value = test$p.value
  ) %>%
  ungroup()

# Step 4: FDR correction
results <- results %>%
  mutate(
    fdr_p_value = p.adjust(p_value, method = "fdr"),
    significant = fdr_p_value < 0.05
  )

# Step 5: Final output
final_results <- results %>%
  select(Cancer, high_cancer, low_cancer, high_others, low_others,
         odds_ratio, ci_low, ci_high, p_value, fdr_p_value, significant)

# Save results
final_results = data.frame(final_results)

final_results = read.table('Fig1A_SignificanceInterval.txt',sep='\t',head=T)
final_results$TotalSamples = final_results$high_cancer+final_results$low_cancer
final_results$HighExtendFrac = final_results$high_cancer/final_results$TotalSamples
final_results$LowExtendFrac = final_results$low_cancer/final_results$TotalSamples

write.table(final_results,'Fig1A_SignificanceInterval.txt',sep='\t',quote=FALSE,row.names=FALSE)


## Figure 1A; The significant values asterick are added on top of the bar plot manually for the significant cases

pdf('Fig1A.pdf')
P1 = ggplot(comball,aes(x=Cancer, y=Fraction,fill=EXTENDclass,color=EXTENDclass))+
  geom_bar(stat = "identity",width=0.6,position="dodge",linewidth=0.1)+scale_fill_manual(values=c("High"="darkred","Low"="pink"))+scale_color_manual(values=c("High"="pink","Low"="darkred"))+scale_x_discrete(limits=c("COAD","GBM","HNSC","LAML","LUSC","READ","STAD","TGCT","THYM","UCEC","UVM","ACC","BLCA","BRCA","CESC","CHOL","DLBC","ESCA","KICH","KIRC","KIRP","LGG","LIHC","LUAD","MESO","OV","PAAD","PCPG","PRAD","SARC","SKCM","THCA","UCS"))+theme_classic()+theme(axis.text.x=element_text(angle=35,size=7,vjust=1,hjust=1,color="black"),legend.key.size=unit(0.1,"cm"),legend.position="top")#+coord_flip()

grid.newpage()
# Create layout : nrow = 2, ncol =2
pushViewport(viewport(layout = grid.layout(nrow =9, ncol = 1)))
# A helper function to define a region on the layout
define_region <- function(row, col){
  viewport(layout.pos.row = row, layout.pos.col = col)
} 

print(ggpar(P1,font.xtickslab =c("black",9),font.ytickslab =c("black",10),font.x = 9,font.y=9,font.legend=7),vp = define_region(row = 1:3, col = 1))

dev.off()


############# 1B ###################################

SigTest = read.table('Fig1B.txt',sep='\t',head=T)
ExtendClassFracs = read.table('Fig1A_SignificanceInterval.txt',sep='\t',head=T)
ExtendClassFracs = ExtendClassFracs[,c(1,12:14)]

SigTest = merge(SigTest,ExtendClassFracs,by='Cancer')


SigTest <- SigTest %>%
  mutate(
    Size = case_when(
      LowExtendFrac < 0.3 ~ 1,
      LowExtendFrac >= 0.3 & LowExtendFrac < 0.6 ~ 2,
      LowExtendFrac >= 0.6 ~ 3,
      TRUE ~ NA_real_  # In case LowExtendFrac is NA
    )
  )



SigTest$Size = as.factor(SigTest$Size)

pdf('Fig1B.pdf')


P2 = ggscatter(SigTest, x="meanTL_High", y="meanTL_Low",size ="Size",shape = 1,stroke = 0.8,color = "Significance",palette = c('NS'='azure3','Sig'='red'),label = "Cancer",label.select = c("GBM","TGCT","KIRP","BRCA","THYM","LGG"),repel=TRUE,font.label = c(6,"plain","black"))+theme_classic()+theme(axis.text.x=element_text(angle=0,size=7,vjust=1,hjust=1),legend.key.size=unit(0.1,"cm"),legend.position="top")+geom_abline(slope=1,intercept=0)


grid.newpage()
# Create layout : nrow = 2, ncol =2
pushViewport(viewport(layout = grid.layout(nrow =11, ncol = 11)))
# A helper function to define a region on the layout
define_region <- function(row, col){
  viewport(layout.pos.row = row, layout.pos.col = col)
} 

print(ggpar(P2,font.xtickslab =c("black",10),font.ytickslab =c("black",10),font.x = 10,font.y=10,ylim=c(-1,0.6),font.legend=10,xlab="High EXTEND TL Ratio",ylab="Low EXTEND TL Ratio"),vp = define_region(row = 1:4, col = 1:4))

dev.off()



############# 1C ###################################
AS_Scores_Pancan = read.table('Fig1C.txt',sep='\t',head=T)

ExtendClassFracs = read.table('Fig1A_SignificanceInterval.txt',sep='\t',head=T)
ExtendClassFracs = ExtendClassFracs[,c(1,12:14)]


AS_Scores_Pancan = merge(AS_Scores_Pancan,ExtendClassFracs,by='Cancer')

labes1 = AS_Scores_Pancan$Cancer[which(AS_Scores_Pancan$significance == "Sig")]

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

pdf('Fig1C.pdf')

P1 = ggscatter(AS_Scores_Pancan, x="mean_high", y="mean_low",size ="Size",shape = 1,stroke = 0.8,color = "significance",palette = c('NS'='azure3','Sig'='red'),label="Cancer",label.select= labes1,repel=T,font.label=c("black",5))+theme_classic()+theme(axis.text.x=element_text(angle=0,size=7,vjust=1,hjust=1,color="black"),legend.key.size=unit(0.05,"cm"),legend.position="top")+geom_abline(slope=1,intercept=0)


grid.newpage()
# Create layout : nrow = 2, ncol =2
pushViewport(viewport(layout = grid.layout(nrow =11, ncol = 11)))
# A helper function to define a region on the layout
define_region <- function(row, col){
  viewport(layout.pos.row = row, layout.pos.col = col)
} 

print(ggpar(P1,font.xtickslab =c("black",10),font.ytickslab =c("black",10),font.x = 10,font.y=10,font.legend=9,xlab="High EXTEND CNA",ylab="Low EXTEND CNA",ylim=c(0,0.85)),vp = define_region(row = 1:4, col = 1:4))

dev.off()

############# 1D ###################################


AS_Scores_Pancan = read.table('Fig1D.txt',sep='\t',head=T)


ExtendClassFracs = read.table('Fig1A_SignificanceInterval.txt',sep='\t',head=T)
ExtendClassFracs = ExtendClassFracs[,c(1,12:14)]


AS_Scores_Pancan = merge(AS_Scores_Pancan,ExtendClassFracs,by='Cancer')


labes1 = AS_Scores_Pancan$Cancer[which(AS_Scores_Pancan$significance == "Sig")]


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

pdf('Fig1D.pdf')

P1 = ggscatter(AS_Scores_Pancan, x="mean_high", y="mean_low",size ="Size",shape = 1,stroke = 0.8,color = "significance",palette = c('NS'='azure3','Sig'='red'),label="Cancer",label.select= labes1,repel=T,font.label=c("black",5))+theme_classic()+theme(axis.text.x=element_text(angle=0,size=7,vjust=1,hjust=1,color="black"),legend.key.size=unit(0.05,"cm"),legend.position="top")+geom_abline(slope=1,intercept=0)


grid.newpage()
# Create layout : nrow = 2, ncol =2
pushViewport(viewport(layout = grid.layout(nrow =11, ncol = 11)))
# A helper function to define a region on the layout
define_region <- function(row, col){
  viewport(layout.pos.row = row, layout.pos.col = col)
} 

print(ggpar(P1,font.xtickslab =c("black",10),font.ytickslab =c("black",10),font.x = 10,font.y=10,font.legend=9,xlab="High EXTEND LOH",ylab="Low EXTEND LOH",ylim=c(0,0.5)),vp = define_region(row = 1:4, col = 1:4))

dev.off()


###################################################################################################################################################
#############################################
#############################################                   Figure 2 : Tumor Mutation Burden
#############################################
###################################################################################################################################################

############################################# 2A #############################################


samples <- read.table('Fig2A.txt',sep='\t',head=T)

samples$mutLoad_nonsilentLog = log(samples$mutLoad_nonsilent)

samples = samples[order(samples$Cancer),]
can_types <- unique(samples$Cancer)

Extend_high  <- samples[which(samples$EXTENDclass=='High'),]
Extend_low  <- samples[which(samples$EXTENDclass=='Low'),]


pdf('Fig2A.pdf')
par(xaxs='i',mar=c(25,5,4,2))
plot.new()
plot.window(xlim=c(0.3,length(can_types)+1.0),ylim=c(-4,6))

for (i in 1:length(can_types)){

  j <- i-0.2
  k <- i+0.2
  
  can_high <- Extend_high[which(Extend_high$Cancer==can_types[i]),]
  N=dim(can_high)[1]
  
  idx=seq(j-0.2,j+0.2,length=N)
  a=can_high$mutLoad_nonsilentLog
  points(idx,a,pch=20,type='p',cex=0.3, col='red')
  mv=median(a)
  lines(x=c(j-0.15,j+0.15),y=c(mv,mv),lwd=2,col='black')
  
  can_low <- Extend_low[which(Extend_low$Cancer==can_types[i]),]
  N=dim(can_low)[1]
  
  idx=seq(k-0.2,k+0.2,length=N)
  a=can_low$mutLoad_nonsilentLog
  points(idx,a,pch=20,type='p',cex=0.3, col='blue')
  mv=median(a)
  lines(x=c(k-0.15,k+0.15),y=c(mv,mv),lwd=2,col='black')

}

par(xpd=FALSE)
abline(h=10,lwd=1,col='blue')
axis(1,labels=NA,tick=F)
axis(2,at=seq(-4,6,by=1),las=2, cex.axis=0.6)
mtext('TMB', side=2,line=3)
box()
text(1:length(can_types), par("usr")[3]-0.0, 
     srt = 45, pos=2,offset=-0.1, xpd = TRUE,
     labels=can_types, cex = 0.6)
	 
dev.off()


####################### Fig 2B

Fig = read.table('Fig2B.txt',sep='\t',head=T)

Fig$Size = as.factor(Fig$Size)

pdf('Fig2B.pdf')


P3 = ggscatter(Fig, x = "Gene", y = "Cancer",color = "ClassColor",palette = c('Low'='blue2','High'='deeppink'),size = "Size")+theme_linedraw()+theme(axis.text.x=element_text(angle= 50,size=8,vjust=1,hjust=1),legend.key.height=unit(0.1,"cm"),legend.position="bottom")

grid.newpage()
# Create layout : nrow = 2, ncol =2
pushViewport(viewport(layout = grid.layout(nrow = 11, ncol = 9)))
# A helper function to define a region on the layout
define_region <- function(row, col){
  viewport(layout.pos.row = row, layout.pos.col = col)
}
print(ggpar(P3,font.xtickslab =6,font.ytickslab =6,font.x = 7,font.y=7,font.title=6,font.legend=6,),vp = define_region(row = 1:7, col = 2:5))  

dev.off()

###################### Fig 2C

### The figure is plotted in 2 sections and then combined in adobe illustrator

# ---- Load and split data ----
FigC <- read.table("Fig2C.txt", sep = "\t", header = TRUE)


pdf('Fig.2C.pdf')
P1 = ggplot(FigC, aes(x=GeneOrder, y=HighPercent,fill=Cancer)) + 
  geom_bar(stat = "identity",width=0.7,color='black',size=0.1)+scale_fill_manual(values=c("BLCA"="red","BRCA"="blue","HNSC"="green4","LGG"="khaki","LIHC"="limegreen","LUAD"="darkblue","GBM"="orange","PAAD"="pink","PRAD"="darkolivegreen3","SARC"="purple","SKCM"="bisque1","STAD"="brown4","TGCT"="grey","THCA"="darkgoldenrod1","THYM"="blueviolet","UCS"="hotpink","UCEC"="cornflowerblue"))+theme_classic()+theme(legend.key.size=unit(0.1,"cm"))+coord_flip()


P2 = ggplot(FigC, aes(x=GeneOrder, y=LowPercent,fill=Cancer)) + 
  geom_bar(stat = "identity",width=0.7,color='black',size=0.1)+scale_fill_manual(values=c("BLCA"="red","BRCA"="blue","HNSC"="green4","LGG"="khaki","LIHC"="limegreen","LUAD"="darkblue","GBM"="orange","PAAD"="pink","PRAD"="darkolivegreen3","SARC"="purple","SKCM"="bisque1","STAD"="brown4","TGCT"="grey","THCA"="darkgoldenrod1","THYM"="blueviolet","UCS"="hotpink","UCEC"="cornflowerblue"))+theme_classic()+theme(legend.key.size=unit(0.1,"cm"))+coord_flip()


grid.newpage()
# Create layout : nrow = 2, ncol =2
pushViewport(viewport(layout = grid.layout(nrow = 7, ncol = 6)))
# A helper function to define a region on the layout
define_region <- function(row, col){
  viewport(layout.pos.row = row, layout.pos.col = col)
} 
print(ggpar(P1,font.xtickslab =6,font.ytickslab =6,font.x = 7,font.y=7,font.title=6,font.legend=6),vp = define_region(row = 1:4, col = 1:2))  
print(ggpar(P2,font.xtickslab =6,font.ytickslab =6,font.x = 7,font.y=7,font.title=6,font.legend=6,),vp = define_region(row = 1:4, col = 3:4))  

dev.off()  #### High and Low EXTEND Mutation Panels of the Figure are combined in Adobe to finalize the plots


################### Fig 2D

FigD = read.table('Fig2D.txt',sep='\t',head=T)

FigD = FigD[which(FigD$Status == "Mutant"),]

### Significance is added manually on top of the plot from this data, and is only for the cases where both Extend groups are available in mutant category

pdf('Fig2D.pdf')
P1 = ggplot(FigD, aes(x=Cancer, y=RatioClass,fill=EXTENDclass,color=EXTENDclass)) + 
  geom_bar(stat = "identity", position=position_dodge(),width=0.7,size=0.1)+scale_fill_manual(values=c("High"="red","Low"="blue"))+scale_color_manual(values=c("High"="red","Low"="blue"))+theme_classic()+theme(axis.text.x=element_text(angle= 50,size=8,vjust=1,hjust=1),legend.key.height=unit(0.1,"cm"),legend.position="bottom")
##+coord_flip()



grid.newpage()
# Create layout : nrow = 2, ncol =2
pushViewport(viewport(layout = grid.layout(nrow = 7, ncol = 8)))
# A helper function to define a region on the layout
define_region <- function(row, col){
  viewport(layout.pos.row = row, layout.pos.col = col)
} 
print(ggpar(P1,font.xtickslab =8,font.ytickslab =8,font.x = 8,font.y=8,font.title=6,font.legend=6),vp = define_region(row = 1:2, col = 1:4))  


dev.off()

####### Fig 2E 

FigE = read.table('Fig2E.txt',sep='\t',head=T)

FigE = FigE[which(FigE$Status == "Mutant"),]
FigE = FigE[which(FigE$Cancer %in% c("GBM","LGG","SARC")),]

pdf('Fig2E.pdf')
P1 = ggplot(FigE, aes(x=Cancer, y=RatioClass,fill=EXTENDclass,color=EXTENDclass)) + 
  geom_bar(stat = "identity", position=position_dodge(),width=0.8,size=0.1)+scale_fill_manual(values=c("High"="red","Low"="blue"))+scale_color_manual(values=c("High"="red","Low"="blue"))+theme_classic()+theme(axis.text.x=element_text(angle= 50,size=8,vjust=1,hjust=1),legend.key.height=unit(0.1,"cm"),legend.position="bottom")
##+coord_flip()



grid.newpage()
# Create layout : nrow = 2, ncol =2
pushViewport(viewport(layout = grid.layout(nrow = 7, ncol = 8)))
# A helper function to define a region on the layout
define_region <- function(row, col){
  viewport(layout.pos.row = row, layout.pos.col = col)
} 
print(ggpar(P1,font.xtickslab =8,font.ytickslab =8,font.x = 8,font.y=8,font.title=6,font.legend=6),vp = define_region(row = 1:2, col = 1:2))  


dev.off()


###################################################################################################################################################
#############################################
#############################################                   Figure 3 : Fusion Analysis
#############################################
###################################################################################################################################################

##### Fusion Analysis


################## Figure 3A

SigLevel = read.table('Fig3A.txt',sep='\t',head=T)

SigLevel$Size = as.factor(SigLevel$Size)

SigLevel$Low_Ratio2 = SigLevel$Low_Ratio/100
SigLevel$High_Ratio2 = SigLevel$High_Ratio/100
SigLevel$High_Ratio2 = round(SigLevel$High_Ratio2,2)

SigLevel$Low_Ratio2 = round(SigLevel$Low_Ratio2,2)

SigLevel= SigLevel[order(SigLevel$Cancer),]
SigLevel$Order = paste(SigLevel$RepOrd,SigLevel$Fusion,sep='.')

pdf('Fig3A.pdf')
P1 = ggplot(SigLevel, aes(x=Order, y=Low_Ratio2,fill=Cancer)) + 
  geom_bar(stat = "identity",width=0.8,color='black',size=0.1)+scale_fill_manual(values=c("ACC"="palevioletred4",'BLCA'='skyblue',"CESC"="aquamarine",'CHOL'='blue',"COAD"="blanchedalmond","GBM"="red","KIRP"="plum2","LAML"="grey","LGG"="green","LIHC"="black","OV"="darkolivegreen","PAAD"="orange","PCPG"="cyan", "PRAD"="gold4", "SARC"="purple", "STAD"="darkolivegreen1","TGCT"="magenta3", "THCA"="lightcoral","UCEC"="deeppink","ESCA"="cornflowerblue","THYM"="aquamarine4","LUSC"="green3"))+theme_classic()+theme(legend.key.size=unit(0.1,"cm"),legend.position = 'top')+coord_flip()


P2 = ggplot(SigLevel, aes(x=Order, y=High_Ratio2,fill=Cancer)) + 
  geom_bar(stat = "identity",width=0.8,color='black',size=0.1)+scale_fill_manual(values=c("ACC"="palevioletred4",'BLCA'='skyblue',"CESC"="aquamarine",'CHOL'='blue',"COAD"="blanchedalmond","GBM"="red","KIRP"="plum2","LAML"="grey","LGG"="green","LIHC"="black","OV"="darkolivegreen","PAAD"="orange","PCPG"="cyan", "PRAD"="gold4", "SARC"="purple", "STAD"="darkolivegreen1","TGCT"="magenta3", "THCA"="lightcoral","UCEC"="deeppink","ESCA"="cornflowerblue","THYM"="aquamarine4","LUSC"="green3"))+theme_classic()+theme(legend.key.size=unit(0.1,"cm"),legend.position = 'top')+coord_flip()

  
grid.newpage()
# Create layout : nrow = 2, ncol =2
pushViewport(viewport(layout = grid.layout(nrow = 8, ncol = 6)))
# A helper function to define a region on the layout
define_region <- function(row, col){
  viewport(layout.pos.row = row, layout.pos.col = col)
} 
print(ggpar(P1,font.xtickslab =7,font.ytickslab =6,font.x = 7,font.y=7,font.title=6,font.legend=6),vp = define_region(row = 1:6, col = 1:3))  
print(ggpar(P2,font.xtickslab =7,font.ytickslab =6,font.x = 7,font.y=7,font.title=6,font.legend=6,),vp = define_region(row = 1:6, col =4:6))  

dev.off()



################## Figure 3B


dataFusions = read.table('Fig3B.txt',sep='\t',head=T)

pdf('Fig3B.pdf')


P2 = ggscatter(dataFusions, x="Cancer", y="DrugTarget",size ="Size",color = "Label",facet.by='EXTENDClass',palette = c('OffLabel'='salmon','OnLabel'='mediumvioletred'))+theme_linedraw()+theme(axis.text.x=element_text(angle=90,size=7,vjust=1,hjust=1),legend.key.size=unit(0.1,"cm"),legend.position="top")

grid.newpage()
# Create layout : nrow = 2, ncol =2
pushViewport(viewport(layout = grid.layout(nrow =5, ncol = 6)))
# A helper function to define a region on the layout
define_region <- function(row, col){
  viewport(layout.pos.row = row, layout.pos.col = col)
} 

print(ggpar(P2,font.xtickslab =6,font.ytickslab =6,font.x = 8,font.y=8,font.legend=9),vp = define_region(row = 1:3, col = 1:6))

dev.off()





###################################################################################################################################################
#############################################
#############################################                   Figure 4 : Senescence
#############################################
###################################################################################################################################################

##### 4A

SenSig = read.table('Fig4A.txt',sep='\t',head=T)


labsel = SenSig$Cancer[which(SenSig$Sig == "Sig")]


SenSig$Size = NULL

for(i in 1:nrow(SenSig)){
if(SenSig$LowExtendFrac[i] < 0.5){
SenSig$Size[i]=1}else if(SenSig$LowExtendFrac[i] >= 0.5 & SenSig$LowExtendFrac[i] < 0.6){
SenSig$Size[i]=2}else if(SenSig$LowExtendFrac[i] >= 0.6){
 SenSig$Size[i]=3}
 }


SenSig$Size = as.factor(SenSig$Size)
pdf('Fig4A.pdf')


P1 = ggscatter(SenSig, y="MeanLow", x="MeanHigh",size ="Size",color = "Sig",palette = c('Sig'='red','NS'='blue'),alpha=0.6,,label = "Cancer",font.label = c(6,"plain","black"),label.select = labsel,repel=TRUE)+theme_classic()+theme(axis.text.x=element_text(angle=0,size=7,vjust=1,hjust=1),legend.key.size=unit(0.1,"cm"),legend.position="top")+geom_abline(slope=1,intercept=0)##+geom_hline(yintercept=0.45)+geom_vline(xintercept =0.4)



grid.newpage()
# Create layout : nrow = 2, ncol =2
pushViewport(viewport(layout = grid.layout(nrow =5, ncol = 5)))
# A helper function to define a region on the layout
define_region <- function(row, col){
  viewport(layout.pos.row = row, layout.pos.col = col)
} 

print(ggpar(P1,font.xtickslab =c("black",10),font.ytickslab =c("black",10),font.x = 10,font.y=10,font.legend=9,xlab="High EXTEND Senescence",ylab="Low EXTEND Senescence",xlim=c(0,0.6),ylim=c(0,0.6)),vp = define_region(row = 1:2, col = 1:2))


dev.off()


##### 4B CCLE Senescence Scores 

Fig4B = read.table('Fig4B.txt',sep='\t',head=T)


source("https://gist.githubusercontent.com/benmarwick/2a1bb0133ff568cbe28d/raw/fb53bd97121f7f9ce947837ef1a4c65a73bffb3f/geom_flat_violin.R")

my_comparisons1 = list(c("High","Low"))


pdf('Fig4B.pdf')
  
P1 = ggplot(data = Fig4B, 
         aes(x = EXTENDclass, y = SenescenceGenes, fill = "dodgerblue4")) +
  geom_flat_violin(position = position_nudge(x = .2, y = 0),fill="dodgerblue4",color="darkblue",alpha=0.5) +
   #geom_point(aes(y = Senescence),color="blue",fill = "white",shape = 21,stroke = .1,  
            # position = position_jitter(width = .15), size = .2)+  
  geom_boxplot(width = .15, outlier.shape = NA,fill="white",color="black",size=0.2) +
  guides(fill = FALSE) +
  guides(color = FALSE) +
  scale_color_brewer(palette = "Spectral") +
  scale_fill_brewer(palette = "Spectral") +  
  theme_classic()+theme(axis.text.x=element_text(size=8,vjust=1,hjust=1),legend.position="top")+stat_compare_means(comparisons = my_comparisons1, method= "t.test",size=3)# coord_flip() + # flip or not

  ##$#coord_flip()+
grid.newpage()
# Create layout : nrow = 2, ncol =2
pushViewport(viewport(layout = grid.layout(nrow = 7, ncol = 6)))
# A helper function to define a region on the layout
define_region <- function(row, col){
  viewport(layout.pos.row = row, layout.pos.col = col)
} 


print(ggpar(P1,font.xtickslab =c("black",10),font.ytickslab =c("black",10),font.x = 10,font.y=10,ylab="Senescence Scores",xlab="EXTEND Class"),vp = define_region(row = 1:2, col = 1:2))

dev.off()


###### Fig 4C , 4D and 4G


FigC = readRDS('Fig4C_Panel.RDS')

pdf('Fig4C-G.pdf')

P1= FeaturePlot(FigC, features = "EXTEND",min.cutoff = 0.2, max.cutoff = 0.5)+ scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))+theme_classic()+theme(axis.text.x=element_text(angle= 0,size=8,vjust=1,hjust=1),legend.key.height=unit(0.3,"cm"),legend.key.width=unit(0.5,"cm"),legend.position="top")

P2 = FeaturePlot(FigC, features = "SenescenceScores",min.cutoff = 0.2, max.cutoff = 0.5)+ scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))+theme_classic()+theme(axis.text.x=element_text(angle= 0,size=8,vjust=1,hjust=1),legend.key.height=unit(0.3,"cm"),legend.key.width=unit(0.5,"cm"),legend.position="top")

P3 = DimPlot(FigC, reduction = 'umap',group.by='Phase',cols = c('NonCycling'='#FF68A1','G2_M'='#00A9FF','G1_S'='#00BFC4'))+theme_classic()+theme(axis.text.x=element_text(angle= 0,size=8,vjust=1,hjust=1),legend.key.height=unit(0.5,"cm"),legend.key.width=unit(0.5,"cm"),legend.position="top")

grid.newpage()
# Create layout : nrow = 2, ncol =2
pushViewport(viewport(layout = grid.layout(nrow = 5, ncol = 3)))
# A helper function to define a region on the layout
define_region <- function(row, col){
  viewport(layout.pos.row = row, layout.pos.col = col)
} 

print(ggpar(P1,font.xtickslab =8,font.ytickslab =8,font.x = 8,font.y=8, font.legend=8),vp = define_region(row = 1:2, col = 1))
print(ggpar(P2,font.xtickslab =8,font.ytickslab =8,font.x = 8,font.y=8, font.legend=8),vp = define_region(row = 1:2, col = 2))
print(ggpar(P3,font.xtickslab =8,font.ytickslab =8,font.x = 8,font.y=8, font.legend=8),vp = define_region(row = 1:2, col = 3))

dev.off()


###### Fig 4E, 4F and 4H


mydata = readRDS('Fig4E_Panel.RDS')


pdf('Fig4E_Panel.pdf')
P1= FeaturePlot(mydata, features = "EXTEND",min.cutoff = 0, max.cutoff = 0.57,pt.size=0.5)+ scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))+theme_classic()+theme(axis.text.x=element_text(angle= 0,size=8,vjust=1,hjust=1),legend.key.height=unit(0.3,"cm"),legend.key.width=unit(0.5,"cm"),legend.position="top")

P2 = FeaturePlot(mydata, features = "SenescenceScores",min.cutoff = 0.1, max.cutoff = 0.57,pt.size=0.5)+ scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))+theme_classic()+theme(axis.text.x=element_text(angle= 0,size=8,vjust=1,hjust=1),legend.key.height=unit(0.3,"cm"),legend.key.width=unit(0.5,"cm"),legend.position="top")


P3 = DimPlot(mydata, reduction = 'umap',group.by='Phase',cols = c('NonCycling'='#FF68A1','G2_M'='#00A9FF','G1_S'='#00BFC4'),pt.size=0.5)+theme_classic()+theme(axis.text.x=element_text(angle= 0,size=8,vjust=1,hjust=1),legend.key.height=unit(0.5,"cm"),legend.key.width=unit(0.5,"cm"),legend.position="top")

grid.newpage()
# Create layout : nrow = 2, ncol =2
pushViewport(viewport(layout = grid.layout(nrow = 5, ncol = 3)))
# A helper function to define a region on the layout
define_region <- function(row, col){
  viewport(layout.pos.row = row, layout.pos.col = col)
} 

print(ggpar(P1,font.xtickslab =8,font.ytickslab =8,font.x = 8,font.y=8, font.legend=8),vp = define_region(row = 1:2, col = 1))
print(ggpar(P2,font.xtickslab =8,font.ytickslab =8,font.x = 8,font.y=8, font.legend=8),vp = define_region(row = 1:2, col = 2))
print(ggpar(P3,font.xtickslab =8,font.ytickslab =8,font.x = 8,font.y=8, font.legend=8),vp = define_region(row = 1:2, col = 3))

dev.off()

###################################################################################################################################################
#############################################
#############################################                   Figure 5 : Spatial transcriptomics Data
#############################################
###################################################################################################################################################


################### Spatial Figure 5

#### 

combMeta = read.table('Fig5_Lung.txt',sep='\t',head=T)
combMeta2 = read.table('Fig5_Breast.txt',sep='\t',head=T)

pdf('Fig5.pdf')

P1 = ggplot(combMeta, aes(x = PosX, y =PosY ,color=EXTENDScores))+geom_point(size=0.4)+scale_color_gradientn(colors = c("yellow","yellow","yellow", "red","red", "red"), limits = c(0, 1),breaks = c(0,0.5,1))+theme_classic()+theme(axis.text.x=element_text(size=8,vjust=1,hjust=1),legend.key.height=unit(0.2,"cm"),legend.position="top")

P2 = ggplot(combMeta, aes(x = PosX, y =PosY ,color=Senescence_AUCell))+geom_point(size=0.4)+scale_color_gradientn(colors = c("yellow","yellow","red", "red"), limits = c(0, 0.2),breaks = c(0,0.1,0.2))+theme_classic()+theme(axis.text.x=element_text(size=8,vjust=1,hjust=1),legend.key.height=unit(0.2,"cm"),legend.position="top")

P3 =ggscatter(combMeta, y="PosY", x="PosX",size =0.4,color = "Phase",palette=c("G1"="yellow","S"="red","G2M"="hotpink"),alpha=0.9)+theme_classic()+theme(legend.key.size= unit(0.4,"cm"),legend.key.width = unit(0.2,"cm"),legend.position = "top")



P4 = ggplot(combMeta2, aes(x = PosX, y =PosY ,color=EXTENDScores))+geom_point(size=0.4)+scale_color_gradientn(colors = c("yellow","yellow","yellow","red", "red","red", "red"), limits = c(0, 1),breaks = c(0,0.5,1))+theme_classic()+theme(axis.text.x=element_text(size=8,vjust=1,hjust=1),legend.key.height=unit(0.2,"cm"),legend.position="top")



P5 = ggplot(combMeta2, aes(x = PosX, y =PosY ,color=Senescence_AUCell))+geom_point(size=0.4)+scale_color_gradientn(colors = c("yellow","yellow","red", "red"), limits = c(0, 0.2),breaks = c(0,0.1,0.2))+theme_classic()+theme(axis.text.x=element_text(size=8,vjust=1,hjust=1),legend.key.height=unit(0.2,"cm"),legend.position="top")



P6 =ggscatter(combMeta2, y="PosY", x="PosX",size =0.4,color = "Phase",palette=c("G1"="yellow","S"="red","G2M"="hotpink"),alpha=0.9)+theme_classic()+theme(legend.key.size= unit(0.4,"cm"),legend.key.width = unit(0.2,"cm"),legend.position = "top")



grid.newpage()
# Create layout : nrow = 2, ncol =2
pushViewport(viewport(layout = grid.layout(nrow = 3, ncol =3)))
# A helper function to define a region on the layout
define_region <- function(row, col){
  viewport(layout.pos.row = row, layout.pos.col = col)
} 



print(ggpar(P1,font.xtickslab =8,font.ytickslab =8,font.x = 8,font.y=8, font.legend=8),vp = define_region(row = 1, col = 1))
print(ggpar(P2,font.xtickslab =8,font.ytickslab =8,font.x = 8,font.y=8, font.legend=8),vp = define_region(row = 1, col = 2))
print(ggpar(P3,font.xtickslab =8,font.ytickslab =7,font.x = 8,font.y=8, font.legend=8),vp = define_region(row = 1, col = 3))

print(ggpar(P4,font.xtickslab =8,font.ytickslab =8,font.x = 8,font.y=8, font.legend=8),vp = define_region(row = 2, col = 1))
print(ggpar(P5,font.xtickslab =8,font.ytickslab =8,font.x = 8,font.y=8, font.legend=8),vp = define_region(row = 2, col = 2))
print(ggpar(P6,font.xtickslab =8,font.ytickslab =8,font.x = 8,font.y=8, font.legend=8),vp = define_region(row = 2, col = 3))


dev.off()


###################################################################################################################################################
#############################################
#############################################                   Figure 6 : Aging and Senescence
#############################################
###################################################################################################################################################


############### Figure 6A Sen TCGA

library(dplyr)
library(tidyr)


comball = read.table('Fig6A.txt',sep='\t',head=T)

comball = comball[order(comball$Cancer),]

comball = comball[-which(comball$Sig=='NS'),]

comball$Size = NA

for(i in 1:nrow(comball)){
if(comball$OA_Total[i]<=100){
comball$Size[i]=1}else if(comball$OA_Total[i]>100 & comball$OA_Total[i]<=200){
comball$Size[i]=2}else if(comball$OA_Total[i]>200 & comball$OA_Total[i]<=300){
comball$Size[i]=3}else if(comball$OA_Total[i]>300){
comball$Size[i]=4}}



pdf('Fig6A.pdf')
   
P3 = ggplot(comball, aes(y = Cancer, x = EXTENDclass,color=Mean))+geom_point(size=comball$Size,shape=19,stroke=0.7)+scale_color_gradient2(midpoint=0.3, low="aquamarine",mid="blue",high="red")+theme_linedraw()+theme(axis.text.x=element_text(angle= 40,size=8,vjust=1,hjust=1),legend.key.height=unit(0.3,"cm"),legend.key.width=unit(0.4,"cm"),legend.position="bottom")




grid.newpage()
# Create layout : nrow = 2, ncol =2
pushViewport(viewport(layout = grid.layout(nrow = 10, ncol = 13)))
# A helper function to define a region on the layout
define_region <- function(row, col){
  viewport(layout.pos.row = row, layout.pos.col = col)
} 

print(ggpar(P3,font.xtickslab =8,font.ytickslab =8,font.x = 8,font.y=8, font.legend=8),vp = define_region(row = 3:8, col = 1:2))

dev.off()


############ Figure 6B to 6E   GTEX violin

combdata2 = read.table('Fig6B_E.txt',sep='\t',head=T)

LungOnly = combdata2[which(combdata2$Tissue == "Lung"),]
EsophagusOnly = combdata2[which(combdata2$Tissue == "Esophagus"),]
SkinOnly = combdata2[which(combdata2$Tissue == "Skin"),]
PituitaryOnly = combdata2[which(combdata2$Tissue == "Brain"),]


pdf('Fig6B_E.pdf',onefile=FALSE)

cols2 <- c("Low" = "skyblue", "High" = "darkblue")

P1 =  ggviolin(LungOnly,x="AgeGroup",y="SenescenceGenes",size =  0.4, width = 0.9,draw_quantiles=c(0.5),outlier.shape=NA,color = "white",fill="EXTENDclass",order = c("AYA","Adult","OA"))+scale_fill_manual(values = cols2, name = "EXTENDclass")+theme_classic()+theme(axis.text.x=element_text(angle=0,size=7,vjust=0,hjust=0.5),legend.position = "top",legend.key.size= unit(0.3,"cm"),legend.key.width = unit(0.3,"cm"),legend.title = element_text(size=7),legend.text =element_text(size=6))+stat_compare_means(aes(group = EXTENDclass,label = paste0("p = ", ..p.format..)),method = "t.test",size=2.5)
  


P2 =  ggviolin(EsophagusOnly,x="AgeGroup",y="SenescenceGenes",size =  0.4, width = 0.9,draw_quantiles=c(0.5),outlier.shape=NA,color = "white",fill="EXTENDclass",order = c("AYA","Adult","OA"))+scale_fill_manual(values = cols2, name = "EXTENDclass")+theme_classic()+theme(axis.text.x=element_text(angle=0,size=7,vjust=0,hjust=0.5),legend.position = "top",legend.key.size= unit(0.3,"cm"),legend.key.width = unit(0.3,"cm"),legend.title = element_text(size=7),legend.text =element_text(size=6))+stat_compare_means(aes(group = EXTENDclass,label = paste0("p = ", ..p.format..)),method = "t.test",size=2.5)
  

P3 =  ggviolin(SkinOnly,x="AgeGroup",y="SenescenceGenes",size =  0.4, width = 0.9,draw_quantiles=c(0.5),outlier.shape=NA,color = "white",fill="EXTENDclass",order = c("AYA","Adult","OA"))+scale_fill_manual(values = cols2, name = "EXTENDclass")+theme_classic()+theme(axis.text.x=element_text(angle=0,size=7,vjust=0,hjust=0.5),legend.position = "top",legend.key.size= unit(0.3,"cm"),legend.key.width = unit(0.3,"cm"),legend.title = element_text(size=7),legend.text =element_text(size=6))+stat_compare_means(aes(group = EXTENDclass,label = paste0("p = ", ..p.format..)),method = "t.test",size=2.5)
  
P4 =  ggviolin(PituitaryOnly,x="AgeGroup",y="SenescenceGenes",size =  0.4, width = 0.9,draw_quantiles=c(0.5),outlier.shape=NA,color = "white",fill="EXTENDclass",order = c("AYA","Adult","OA"))+scale_fill_manual(values = cols2, name = "EXTENDclass")+theme_classic()+theme(axis.text.x=element_text(angle=0,size=7,vjust=0,hjust=0.5),legend.position = "top",legend.key.size= unit(0.3,"cm"),legend.key.width = unit(0.3,"cm"),legend.title = element_text(size=7),legend.text =element_text(size=6))+stat_compare_means(aes(group = EXTENDclass,label = paste0("p = ", ..p.format..)),method = "t.test",size=2.5)
  
grid.newpage()
# Create layout : nrow = 3, ncol = 3
pushViewport(viewport(layout = grid.layout(nrow =12, ncol = 4)))
# A helper function to define a region on the layout
define_region <- function(row, col){
  viewport(layout.pos.row = row, layout.pos.col = col)
} 

print(ggpar(P1,font.xtickslab =c("black",0),font.ytickslab =c("black",10),font.x = 10,font.y=10,font.legend=10,xlab="Lung Tissue",ylab="Senescence Score"),vp = define_region(row = 1:3, col = 1))
print(ggpar(P2,font.xtickslab =c("black",0),font.ytickslab =c("black",10),font.x = 10,font.y=10,font.legend=10,xlab="Esophagus Tissue",ylab="Senescence Score"),vp = define_region(row = 1:3, col = 2))
print(ggpar(P3,font.xtickslab =c("black",0),font.ytickslab =c("black",10),font.x = 10,font.y=10,font.legend=10,xlab="Skin Tissue",ylab="Senescence Score"),vp = define_region(row = 4:6, col = 1))
print(ggpar(P4,font.xtickslab =c("black",0),font.ytickslab =c("black",10),font.x = 10,font.y=10,font.legend=10,xlab="Brain Tissue",ylab="Senescence Score"),vp = define_region(row = 4:6, col = 2))

dev.off()



############# Figure 6F and G Human Development


HumanLiver = read.table('Fig6F.txt',sep='\t',head=T)
HumanHeart = read.table('Fig6G.txt',sep='\t',head=T)


for(i in 1:nrow(HumanLiver)){
  if(HumanLiver$Group[i]=='S1'){
    HumanLiver$Groups[i]='Fetal'}else if(HumanLiver$Group[i] %in% c('S2','S3','S4','S5','S6')){
    HumanLiver$Groups[i]='Young'}else if(HumanLiver$Group[i] %in% c('S7','S8','S9')){
    HumanLiver$Groups[i]='Adults'
    }
    }



for(i in 1:nrow(HumanHeart)){
  if(HumanHeart$Group[i]=='S1'){
    HumanHeart$Groups[i]='Fetal'}else if(HumanHeart$Group[i] %in% c('S2','S3','S4','S5','S6')){
    HumanHeart$Groups[i]='Young'}else if(HumanHeart$Group[i] %in% c('S7')){
    HumanHeart$Groups[i]='Adults'
    }
    }

pdf('Fig6F_G.pdf')


P1 = ggscatter(HumanLiver, y="EXTENDScores", x="SenescenceGenes",size =2.5,shape=21,stroke =0.8,color = "Groups",palette = c("Fetal"="red","Young"="blue","Adults"="green4"), add = "reg.line", 
add.params = list(color = "blue", fill = "grey",size=0.2),conf.int = TRUE)+theme_classic()+theme(axis.text.x=element_text(angle=0,size=7,vjust=1,hjust=1),legend.key.size=unit(0.1,"cm"),legend.position="top")+stat_cor(method="spearman",size=3)


P2 = ggscatter(HumanHeart, y="EXTENDScores", x="SenescenceGenes",size =2.5,shape=21,stroke =0.8,color = "Groups",palette = c("Fetal"="red","Young"="blue","Adults"="green4"), add = "reg.line", 
add.params = list(color = "blue", fill = "grey",size=0.2),conf.int = TRUE)+theme_classic()+theme(axis.text.x=element_text(angle=0,size=7,vjust=1,hjust=1),legend.key.size=unit(0.1,"cm"),legend.position="top")+stat_cor(method="spearman",size=3)


grid.newpage()
# Create layout : nrow = 2, ncol =2
pushViewport(viewport(layout = grid.layout(nrow =5, ncol = 5)))
# A helper function to define a region on the layout
define_region <- function(row, col){
  viewport(layout.pos.row = row, layout.pos.col = col)
} 

print(ggpar(P1,font.xtickslab =c("black",10),font.ytickslab =c("black",10),font.x = 10,font.y=10,font.legend=9),vp = define_region(row = 1:2, col = 1:2))

print(ggpar(P2,font.xtickslab =c("black",10),font.ytickslab =c("black",10),font.x = 10,font.y=10,font.legend=9),vp = define_region(row = 1:2, col = 3:4))


dev.off()


###################################################################################################################################################
#############################################
#############################################                   Figure 7 : Immune Cell Types
#############################################
###################################################################################################################################################

comball = read.table('Fig7A.txt',sep='\t',head=T)
comball = comball[order(comball$Ratio),]
comball$Size = NULL

for(i in 1:nrow(comball)){
if(comball$Ratio[i]<=0.2){
comball$Size[i]=1}else if(comball$Ratio[i]>0.2 & comball$Ratio[i]<=0.4){
comball$Size[i]=2}else if(comball$Ratio[i]>0.4 & comball$Ratio[i]<=0.6){
comball$Size[i]=3}else if(comball$Ratio[i]>0.6){
comball$Size[i]=4}}



comball$Size = as.factor(comball$Size)
comball = comball[order(comball$ImmuneCT),]

comball$Class = paste(comball$EXTENDclass,comball$Significance,sep="_")

subcomb = comball[which(comball$Significance == "Sig"),]

pdf('Fig7A.pdf')

P2 = ggscatter(subcomb, x="ImmuneCT", y="Cancer",size ="Size",shape=21,stroke = 1,alpha=0.9,color = "Class",palette=c("High_Sig"="red","High_NS"="skyblue","Low_Sig"="blue","Low_NS"="pink"))+theme_linedraw()+theme(axis.text.x=element_text(angle=90,size=7,vjust=1,hjust=1,color="black"),legend.key.size=unit(0.05,"cm"),legend.position="top")


grid.newpage()
# Create layout : nrow = 2, ncol =2
pushViewport(viewport(layout = grid.layout(nrow =12, ncol = 13)))
# A helper function to define a region on the layout
define_region <- function(row, col){
  viewport(layout.pos.row = row, layout.pos.col = col)
} 

print(ggpar(P2,font.xtickslab =c("black",9),font.ytickslab =c("black",10),font.x = 10,font.y=10,font.legend=7),vp = define_region(row = 1:8, col = 1:4))

dev.off()



###################################### Figure 7B 


SenSig = read.table('Fig7B.txt',sep='\t',head=T)
colnames(SenSig)[2] = 'ImmuneCT'
SenSig = SenSig[complete.cases(SenSig),]

SenSig$log10FDR = -log10(SenSig$FDR)

SenSig = SenSig[order(SenSig$EffectSize),]
SenSig$Size = NULL

for(i in 1:nrow(SenSig)){
if(SenSig$EffectSize[i] <= -1){
SenSig$Size[i]=3}else if(SenSig$EffectSize[i] > -1 & SenSig$EffectSize[i] <= -0.5){
 SenSig$Size[i]=2}else if(SenSig$EffectSize[i] > -0.5 & SenSig$EffectSize[i] < 0.5){
 SenSig$Size[i]=1}else if(SenSig$EffectSize[i] >= 0.5  & SenSig$EffectSize[i] < 1){
 SenSig$Size[i]=2}else if(SenSig$EffectSize[i] >= 1 ){
 SenSig$Size[i]=3}
 }


SenSig = SenSig[order(SenSig$log10FDR),]


SenSig$Size = as.factor(SenSig$Size)

SenSigsub = SenSig[which(SenSig$Sig == "Sig"),]
SenSigsub = SenSigsub[order(SenSigsub$Cancer),]


pdf('Fig7B.pdf')

P2= ggscatter(SenSigsub, x="ImmuneCT", y="Cancer",size ="Size",shape=21,alpha=0.8,stroke=0.8,color="HighClass",palette = c('High'='red','Low'='blue'))+theme_linedraw()+theme(axis.text.x=element_text(angle=35,size=7,vjust=1,hjust=1),legend.key.size=unit(0.1,"cm"),legend.position="top")+stat_cor(method='spearman',size=3)

##+scale_shape_manual(values=c(7,8))
grid.newpage()
# Create layout : nrow = 2, ncol =2
pushViewport(viewport(layout = grid.layout(nrow =10, ncol = 7)))
# A helper function to define a region on the layout
define_region <- function(row, col){
  viewport(layout.pos.row = row, layout.pos.col = col)
} 

print(ggpar(P2,font.xtickslab =c("black",7),font.ytickslab =c("black",9),font.x = 10,font.y=10,font.legend=9),vp = define_region(row = 1:8, col = 1:4))

dev.off()



###################################################################################################################################################
#############################################
#############################################                   Figure 8 : ROS, MAPK
#############################################
###################################################################################################################################################
############################################# Figure  8A-B


SenSig1 = read.table('Fig8A_ROS.txt',sep='\t',head=T)
SenSig2 = read.table('Fig8B_MAPK.txt',sep='\t',head=T)

pdf('Fig8A_B.pdf')


P1 = ggscatter(SenSig1, y="MeanLow", x="MeanHigh",size ="Size",color = "Sig",palette = c('Sig'='red','NS'='blue'),alpha=0.6,,label = "Cancer",font.label = c(6,"plain","black"),label.select = labsel,repel=TRUE)+theme_classic()+theme(axis.text.x=element_text(angle=0,size=7,vjust=1,hjust=1),legend.key.size=unit(0.1,"cm"),legend.position="top")+geom_abline(slope=1,intercept=0)##+geom_hline(yintercept=0.45)+geom_vline(xintercept =0.4)


P2 = ggscatter(SenSig2, y="MeanLow", x="MeanHigh",size ="Size",color = "Sig",palette = c('Sig'='red','NS'='blue'),alpha=0.6,,label = "Cancer",font.label = c(6,"plain","black"),label.select = labsel,repel=TRUE)+theme_classic()+theme(axis.text.x=element_text(angle=0,size=7,vjust=1,hjust=1),legend.key.size=unit(0.1,"cm"),legend.position="top")+geom_abline(slope=1,intercept=0)##+geom_hline(yintercept=0.45)+geom_vline(xintercept =0.4)



grid.newpage()
# Create layout : nrow = 2, ncol =2
pushViewport(viewport(layout = grid.layout(nrow =5, ncol = 5)))
# A helper function to define a region on the layout
define_region <- function(row, col){
  viewport(layout.pos.row = row, layout.pos.col = col)
} 

print(ggpar(P1,font.xtickslab =c("black",10),font.ytickslab =c("black",10),font.x = 10,font.y=10,font.legend=9,xlab="High EXTEND ROS",ylab="Low EXTEND ROS",xlim=c(1.2,1.8),ylim=c(1.2,1.8)),vp = define_region(row = 1:2, col = 1:2))
##
print(ggpar(P2,font.xtickslab =c("black",10),font.ytickslab =c("black",10),font.x = 10,font.y=10,font.legend=9,xlab="High EXTEND MAPK Signalling",ylab="Low EXTEND MAPK Signalling",xlim=c(1,1.35),ylim=c(1,1.35)),vp = define_region(row = 1:2, col = 3:4))
dev.off()


############################################# Figure  8C


comballDF = read.table('Fig8C.txt',sep='\t',head=T)
comballDF$Size = NULL



for(i in 1:nrow(comballDF)){
if(comballDF$Corr[i] <= 0.55){
comballDF$Size[i]=1}else if(comballDF$Corr[i] > 0.55 & comballDF$Corr[i] <= 0.65){
comballDF$Size[i]=2}else if(comballDF$Corr[i] > 0.65){
 comballDF$Size[i]=3}
 }

comballDF$Size = as.factor(comballDF$Size)
#
pdf('Fig8C.pdf')


P1 = ggscatter(comballDF,y="Cancer", x="Type",shape=21,,size ="Size",stroke=0.8,color = "Sig",palette = c('Sig'='red','NS'='skyblue3'))+theme_linedraw()+theme(axis.text.x=element_text(angle=0,size=7,vjust=1,hjust=1),legend.key.size=unit(0.1,"cm"),legend.position="right")##+geom_abline(slope=1,intercept=0)##+geom_hline(yintercept=0.45)+geom_vline(xintercept =0.4)


grid.newpage()
# Create layout : nrow = 2, ncol =2
pushViewport(viewport(layout = grid.layout(nrow =5, ncol = 11)))
# A helper function to define a region on the layout
define_region <- function(row, col){
  viewport(layout.pos.row = row, layout.pos.col = col)
} 

print(ggpar(P1,font.xtickslab =c("black",10),font.ytickslab =c("black",10),font.x = 10,font.y=10,font.legend=9),vp = define_region(row = 1:5, col = 1:3))


dev.off()



############################################# Figure  8D

comballC3 = read.table('Fig8D.txt',sep='\t',head=T)

library(ggsankey)

pdf('Fig8D.pdf')
ggplot(df, aes(x = x, 
               next_x = next_x, 
               node = node, 
               next_node = next_node,
               fill = factor(node),
               label = node)) +
  geom_sankey(flow.alpha = 0.5, node.color = 1) +
  geom_sankey_label(size = 3.5, color = 1, fill = "white") +
  scale_fill_viridis_d(option = "A", alpha = 0.95) +
  theme_sankey(base_size = 16)
dev.off()
