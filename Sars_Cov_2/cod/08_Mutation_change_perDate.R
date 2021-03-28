library(ggplot2)
library(lubridate)
library(dplyr)
library(plotly)
library(data.table)
rm(list=ls(all=TRUE))
#df = read.csv("../../Sars_Cov_2/data_obtained/ideal_table.csv")
df = read.csv('D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/data/norm_data_modernized_nextstrain.csv')
#nstrain = read.csv("../../Sars_Cov_2/data/norm_data_modernized_nextstrain.csv")
ann = read.csv("D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/data_obtained/ideal_table.csv")
nstrain = subset(df, RefNuc != '-' & AltNuc != '-')

for (i in 1:nrow(nstrain)){
  temporary_subset = subset(ann, Pos == nstrain$Pos[i] & AltNuc == nstrain$AltNuc[i])
  if (length(temporary_subset$RefNuc) != 0){
    nstrain$GenName[i] = temporary_subset$GenName
    nstrain$GenType[i] = temporary_subset$GenType
    nstrain$AaSub[i] = temporary_subset$AaSub
    nstrain$NeighL[i] = temporary_subset$NeighL
    nstrain$NeighR[i] = temporary_subset$NeighR
  }
}
nstrain <- nstrain %>%
    arrange(Time)

nstrain = nstrain[nstrain$GenType == 'translated',]
ann = ann[ann$GenType == 'translated',]
nstrain = nstrain[nstrain$AaSub == 'S',]
ann = ann[ann$AaSub == 'S',]

ann$MutExp <- 1
ann$NucSubst = paste(ann$RefNuc,ann$AltNuc,sep='>')
Exp12Comp = aggregate(ann$MutExp, by = list(ann$NucSubst), FUN = sum);
names(Exp12Comp) = c('NucSubst','For_all_dates')

MutSpec12Comp = c('A>C','A>G','A>T','C>A','C>G','C>T','G>A','G>C','G>T','T>A','T>C','T>G')
MutSpec12Comp = data.frame(MutSpec12Comp)
names(MutSpec12Comp)=c('NucSubst')
t1 = MutSpec12Comp
for (i in seq(501, length(as.Date(nstrain$Time)), 500)){
  if ((i+500) <= length(nstrain$Time)){
    time_set = nstrain[i:(i+500),]
    vremia = nstrain$Time[i]
    For_Mut_Spec = merge(t1,Exp12Comp, by = 'NucSubst',all.x = TRUE)
    
    time_set$MutObs = 1
    time_set$NucSubst = paste(time_set$RefNuc,time_set$AltNuc,sep='>')
    Obs12Comp = aggregate(time_set$MutObs, by = list(time_set$NucSubst), FUN = sum);
    names(Obs12Comp) = c('NucSubst', vremia)
    
    For_Mut_Spec = merge(Exp12Comp,Obs12Comp, by = 'NucSubst',all.x = TRUE)
    For_Mut_Spec[is.na(For_Mut_Spec)] <- 0
    For_Mut_Spec$ObsToExp = For_Mut_Spec[[vremia]]/For_Mut_Spec$For_all_dates
    For_Mut_Spec$ObsToExp[is.na(For_Mut_Spec$ObsToExp)] <- 0
    For_Mut_Spec[[vremia]] = For_Mut_Spec$ObsToExp/sum(For_Mut_Spec$ObsToExp)
    Mut_NucSub = cbind(For_Mut_Spec$NucSubst, For_Mut_Spec[[vremia]]);
    colnames(Mut_NucSub) = c('NucSubst', vremia)
    MutSpec12Comp = merge(MutSpec12Comp, Mut_NucSub, by = 'NucSubst',all.x = TRUE)
    MutSpec12Comp[is.na(MutSpec12Comp)] <- 0
  }
}

write.csv(MutSpec12Comp, "D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/data_obtained/Single_Mut_PerDate/08_MutSpec12Comp_AllMut_AllDate.csv")

molten <- melt(setDT(MutSpec12Comp), id.vars = c("NucSubst"))
fig = ggplot(molten, aes(x = variable, y = as.numeric(value), group = NucSubst, colour = NucSubst)) +
  theme(axis.text.x = element_text(angle = 90))+
  xlab("Dates") +
  ylab("Share Of Mutations")+
  ggtitle('Mutation change over time')+
  geom_line() + 
  geom_point() +
  scale_y_continuous(limits=c(0, 1))
ggsave('D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/R_figure/Single_Mut_PerDate/08_MutChange_12spec.svg', plot = fig, device = 'svg',  limitsize = F)

ann$MutExp = 1
ann$FromWithNeigh = paste(ann$NeighL,ann$RefNuc,ann$NeighR, sep='')
ann$ToWithNeigh = paste(ann$NeighL,ann$AltNuc,ann$NeighR, sep='')
ann$MutSpecWithNeigh = paste(ann$FromWithNeigh, ann$ToWithNeigh, sep = '>')
Exp192Comp = aggregate(ann$MutExp, by = list(ann$MutSpecWithNeigh), FUN = sum);
names(Exp192Comp) = c('NucSubst','For_all_dates')

MutSpec192Comp = Exp192Comp$NucSubst
MutSpec192Comp = data.frame(MutSpec192Comp)
names(MutSpec192Comp)=c('NucSubst')
t2 = MutSpec192Comp

for (i in seq(501, length(as.Date(nstrain$Time)), 500)){
  if ((i+500) <= length(nstrain$Time)){
    time_set = nstrain[i:(i+500),]
    vremia = nstrain$Time[i]
    For_Mut_Spec = merge(t2,Exp192Comp, by = 'NucSubst',all.x = TRUE)
    
    time_set$MutObs = 1
    
    time_set$FromWithNeigh = paste(time_set$NeighL,time_set$RefNuc,time_set$NeighR, sep='')
    time_set$ToWithNeigh = paste(time_set$NeighL,time_set$AltNuc,time_set$NeighR, sep='')
    time_set$MutSpecWithNeigh = paste(time_set$FromWithNeigh, time_set$ToWithNeigh, sep = '>')
    Obs192Comp = aggregate(time_set$MutObs, by = list(time_set$MutSpecWithNeigh), FUN = sum);
    names(Obs192Comp) = c('NucSubst', vremia)
    
    For_Mut_Spec = merge(Exp192Comp,Obs192Comp, by = 'NucSubst',all.x = TRUE)
    For_Mut_Spec[is.na(For_Mut_Spec)] <- 0
    For_Mut_Spec$ObsToExp = For_Mut_Spec[[vremia]]/For_Mut_Spec$For_all_dates
    For_Mut_Spec$ObsToExp[is.na(For_Mut_Spec$ObsToExp)] <- 0
    For_Mut_Spec[[vremia]] = For_Mut_Spec$ObsToExp/sum(For_Mut_Spec$ObsToExp)
    Mut_NucSub = cbind(For_Mut_Spec$NucSubst, For_Mut_Spec[[vremia]]);
    colnames(Mut_NucSub) = c('NucSubst', vremia)
    MutSpec192Comp = merge(MutSpec192Comp, Mut_NucSub, by = 'NucSubst',all.x = TRUE)
    MutSpec192Comp[is.na(MutSpec192Comp)] <- 0
  }
}

write.csv(MutSpec192Comp, "D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/data_obtained/Single_Mut_PerDate/08_MutSpec192Comp_AllMut_AllDate.csv")
for (i in MutSpec192Comp$NucSubst){
  new_i = i
  substr(new_i, 4, 4) <- "_"
  date_df = subset(MutSpec192Comp, NucSubst == i, row.names = F)
  date_df2 <- date_df[,-1]
  rownames(date_df2) <- date_df[,1]
  new_date = as.data.frame(t(date_df2))
  fig = ggplot(data = new_date, aes(x = rownames(new_date), y = as.numeric(new_date[,1]), group = 1)) + 
    theme(axis.text.x = element_text(angle = 90))+
    xlab("Dates") +
    ylab("Share Of Mutations")+
    ggtitle(sprintf('Mutation change over time %s', i))+
    geom_line(linetype = "dashed") + 
    geom_point() +
    scale_y_continuous(limits=c(0, 0.1))
  ggsave(sprintf('D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/R_figure/Single_Mut_PerDate/08_MutChange_192spec_%s.svg', new_i), plot = fig, device = 'svg',  limitsize = F)
}