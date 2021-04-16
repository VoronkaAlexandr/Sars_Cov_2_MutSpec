library("heatmaply")
library(ggplot2)
library(lubridate)
library(dplyr)
library(plotly)
library(data.table)

df = read.csv("D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/data_obtained/Single_Mut_PerDate/08_MutSpec12Comp_AllMut_AllDate.csv", check.names=FALSE, row.names = 2)
df_2 = read.csv("D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/data_obtained/Single_Mut_PerDate/08_MutSpec192Comp_AllMut_AllDate.csv", check.names=FALSE, row.names = 2)
df[[1]] <- NULL
df_2[[1]] <- NULL
log_df = log(df+1, base = 2)
log_df_2 = log(df_2+1, base = 2)
heatmaply(
  as.matrix(log_df),
  Rowv = NULL,
  Colv = NULL,
  xlab = "Dates",
  ylab = 'Mutations',
  main = 'Mutation Change 12 components per Date',
  show_dendrogram = c(FALSE, FALSE),
  plot_method = "plotly",
  file = "D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/R_figure/HeatMap/09_Mut_change_12_comp_Date.html"
)
heatmaply(
  as.matrix(log_df_2),
  Rowv = NULL,
  Colv = NULL,
  xlab = "Dates",
  ylab = 'Mutations',
  main = 'Mutation Change 192 components per Date',
  show_dendrogram = c(FALSE, FALSE),
  plot_method = "plotly",
  cexRow = 0.5,
  file = "D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/R_figure/HeatMap/09_Mut_change_192_comp_Date.html"
)

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

nstrain = nstrain[nstrain$GenType == 'translated',]
ann = ann[ann$GenType == 'translated',]
nstrain = nstrain[nstrain$AaSub == 'S',]
ann = ann[ann$AaSub == 'S',]

MutSpec12Comp_heatmap = c('A>C','A>G','A>T','C>A','C>G','C>T','G>A','G>C','G>T','T>A','T>C','T>G')
MutSpec12Comp_heatmap = data.frame(MutSpec12Comp_heatmap)
names(MutSpec12Comp_heatmap)=c('NucSubst')
t1 = MutSpec12Comp_heatmap

for (gen in unique(nstrain$GenName)){
  gen_df = nstrain[nstrain$GenName == gen,]
  gen_ann = ann[ann$GenName == gen,]
  
  gen_ann$MutExp = 1
  gen_ann$NucSubst = paste(gen_ann$RefNuc,gen_ann$AltNuc,sep='>')
  Exp12Comp = aggregate(gen_ann$MutExp, by = list(gen_ann$NucSubst), FUN = sum);
  names(Exp12Comp) = c('NucSubst','ExpFr')
  
  gen_df$MutObs = 1
  gen_df$NucSubst = paste(gen_df$RefNuc,gen_df$AltNuc,sep='>')
  Obs12Comp = aggregate(gen_df$MutObs, by = list(gen_df$NucSubst), FUN = sum);
  names(Obs12Comp) = c('NucSubst','ObsFr')
  
  For_Mut_Spec = merge(t1,Exp12Comp, by = 'NucSubst',all.x = TRUE)
  For_Mut_Spec = merge(For_Mut_Spec,Obs12Comp, by = 'NucSubst',all.x = TRUE)
  For_Mut_Spec[is.na(For_Mut_Spec)] <- 0
  For_Mut_Spec$ObsToExp = For_Mut_Spec$ObsFr/For_Mut_Spec$ExpFr
  For_Mut_Spec$ObsToExp[is.na(For_Mut_Spec$ObsToExp)] <- 0
  For_Mut_Spec[[gen]] = For_Mut_Spec$ObsToExp/sum(For_Mut_Spec$ObsToExp)
  Mut_NucSub = cbind(For_Mut_Spec$NucSubst, For_Mut_Spec[[gen]]);
  colnames(Mut_NucSub) = c('NucSubst', gen)
  MutSpec12Comp_heatmap = merge(MutSpec12Comp_heatmap, Mut_NucSub, by = 'NucSubst',all.x = TRUE)
  MutSpec12Comp_heatmap[is.na(MutSpec12Comp_heatmap)] <- 0
}
write.csv(MutSpec12Comp_heatmap, "D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/data_obtained/Single_Mut_PerDate/09_MutSpec12Comp_AllMut_Genes.csv")

ann$MutExp = 1
ann$FromWithNeigh = paste(ann$NeighL,ann$RefNuc,ann$NeighR, sep='')
ann$ToWithNeigh = paste(ann$NeighL,ann$AltNuc,ann$NeighR, sep='')
ann$MutSpecWithNeigh = paste(ann$FromWithNeigh, ann$ToWithNeigh, sep = '>')
Exp192Comp = aggregate(ann$MutExp, by = list(ann$MutSpecWithNeigh), FUN = sum);
names(Exp192Comp) = c('NucSubst','For_all_dates')

MutSpec192Comp_heatmap = Exp192Comp$NucSubst
MutSpec192Comp_heatmap = data.frame(MutSpec192Comp_heatmap)
names(MutSpec192Comp_heatmap)=c('NucSubst')
t2 = MutSpec192Comp_heatmap

for (gen in unique(nstrain$GenName)){
  gen_df = nstrain[nstrain$GenName == gen,]
  gen_ann = ann[ann$GenName == gen,]
  gen_ann$MutExp = 1
  gen_ann$FromWithNeigh = paste(gen_ann$NeighL,gen_ann$RefNuc,gen_ann$NeighR, sep='')
  gen_ann$ToWithNeigh = paste(gen_ann$NeighL,gen_ann$AltNuc,gen_ann$NeighR, sep='')
  gen_ann$MutSpecWithNeigh = paste(gen_ann$FromWithNeigh, gen_ann$ToWithNeigh, sep = '>')
  Exp192Comp = aggregate(gen_ann$MutExp, by = list(gen_ann$MutSpecWithNeigh), FUN = sum);
  names(Exp192Comp) = c('NucSubst','ExpFr')
 
  gen_df$MutObs = 1
  gen_df$FromWithNeigh = paste(gen_df$NeighL,gen_df$RefNuc,gen_df$NeighR, sep='')
  gen_df$ToWithNeigh = paste(gen_df$NeighL,gen_df$AltNuc,gen_df$NeighR, sep='')
  gen_df$MutSpecWithNeigh = paste(gen_df$FromWithNeigh, gen_df$ToWithNeigh, sep = '>')
  Obs192Comp = aggregate(gen_df$MutObs, by = list(gen_df$MutSpecWithNeigh), FUN = sum);
  names(Obs192Comp) = c('NucSubst','ObsFr')
  
  MutSpec192Comp = merge(t2,Exp192Comp, by = 'NucSubst', all.x = TRUE)
  MutSpec192Comp = merge(MutSpec192Comp,Obs192Comp, by = 'NucSubst', all.x = TRUE)
  MutSpec192Comp[is.na(MutSpec192Comp)] <- 0
  MutSpec192Comp$ObsToExp = MutSpec192Comp$ObsFr/MutSpec192Comp$ExpFr
  MutSpec192Comp$ObsToExp[is.na(MutSpec192Comp$ObsToExp)] <- 0
  MutSpec192Comp$MutSpec = MutSpec192Comp$ObsToExp/sum(MutSpec192Comp$ObsToExp)
  
  Mut_NucSub = cbind(MutSpec192Comp$NucSubst, MutSpec192Comp$MutSpec);
  colnames(Mut_NucSub) = c('NucSubst', gen)
  MutSpec192Comp_heatmap = merge(MutSpec192Comp_heatmap, Mut_NucSub, by = 'NucSubst',all.x = TRUE)
  MutSpec192Comp_heatmap[is.na(MutSpec192Comp_heatmap)] <- 0
}

write.csv(MutSpec192Comp_heatmap, "D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/data_obtained/Single_Mut_PerDate/09_MutSpec192Comp_AllMut_Genes.csv")

gen_map_12 = as.data.frame(MutSpec12Comp_heatmap, row.names = MutSpec12Comp_heatmap[[1]])
gen_map_12[[1]] <- NULL
for (i in colnames(gen_map_12)){
  gen_map_12[[i]] = as.numeric(gen_map_12[[i]])
}
log_gen_map_12 = log(gen_map_12+1, base = 3)
log_gen_map_12 = log_gen_map_12[c(1,2,8,11,5,3,10,7,6,4,9)]
heatmaply(
  as.matrix(log_gen_map_12),
  Rowv = NULL,
  Colv = NULL,
  xlab = "Genes",
  ylab = 'Mutations',
  main = 'Mutation Change 12 components per Gen',
  show_dendrogram = c(FALSE, FALSE),
  plot_method = "plotly",
  file = "D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/R_figure/HeatMap/09_Mut_change_12_comp_Gen.html"
)

gen_map_192 = as.data.frame(MutSpec192Comp_heatmap, row.names = MutSpec192Comp_heatmap[[1]])
gen_map_192[[1]] <- NULL
for (i in colnames(gen_map_192)){
  gen_map_192[[i]] = as.numeric(gen_map_192[[i]])
}
log_gen_map_192 = log(gen_map_192+1, base = 3)
log_gen_map_192 = log_gen_map_192[c(1,2,8,11,5,3,10,7,6,4,9)]
heatmaply(
  as.matrix(log_gen_map_192),
  Rowv = NULL,
  Colv = NULL,
  xlab = "Dates",
  ylab = 'Mutations',
  main = 'Mutation Change 192 components per Gen',
  show_dendrogram = c(FALSE, FALSE),
  plot_method = "plotly",
  cexRow = 0.5,
  file = "D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/R_figure/HeatMap/09_Mut_change_192_comp_Gen.html"
)