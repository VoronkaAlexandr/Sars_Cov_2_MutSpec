library("heatmaply")
library(ggplot2)
library(lubridate)
library(dplyr)
library(plotly)
library(data.table)

#Old data
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
  fontsize_row = 15,
  fontsize_col = 15,
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
  fontsize_row = 15,
  fontsize_col = 15,
  xlab = "Dates",
  ylab = 'Mutations',
  main = 'Mutation Change 192 components per Date',
  show_dendrogram = c(FALSE, FALSE),
  plot_method = "plotly",
  cexRow = 0.5,
  file = "D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/R_figure/HeatMap/09_Mut_change_192_comp_Date.html"
)

rm(list=ls(all=TRUE))

#New Data
gis = read.csv('D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/new_data/data/U_mutations.csv')
#nstrain = read.csv("../../Sars_Cov_2/data/norm_data_modernized_nextstrain.csv")
ann = read.csv("D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/new_data/data/U_ideal_table.csv")

gis = gis[gis$GenType == 'translated',]
ann = ann[ann$GenType == 'translated',]
gis = gis[gis$AaSub == 'S',]
ann = ann[ann$AaSub == 'S',]

MutSpec12Comp_heatmap = c('A>C','A>G','A>U','C>A','C>G','C>U','G>A','G>C','G>U','U>A','U>C','U>G')
MutSpec12Comp_heatmap = data.frame(MutSpec12Comp_heatmap)
names(MutSpec12Comp_heatmap)=c('NucSubst')
t1 = MutSpec12Comp_heatmap

#Drawing Syn MutSpec per Gen

for (gen in unique(gis$GenName)){
  gen_df = gis[gis$GenName == gen,]
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
  if (gen == 'ORF6'){
    For_Mut_Spec[For_Mut_Spec$ExpFr == 0,]$ExpFr = sum(gen_ann[gen_ann$RefNuc=='G',]$MutExp)
  }
  For_Mut_Spec[[gen]] = For_Mut_Spec$ObsFr/For_Mut_Spec$ExpFr
  For_Mut_Spec[[gen]][is.na(For_Mut_Spec[[gen]])] <- 0
  #For_Mut_Spec[[gen]] = For_Mut_Spec$ObsToExp/sum(For_Mut_Spec$ObsToExp)
  Mut_NucSub = cbind(For_Mut_Spec$NucSubst, For_Mut_Spec[[gen]]);
  colnames(Mut_NucSub) = c('NucSubst', gen)
  MutSpec12Comp_heatmap = merge(MutSpec12Comp_heatmap, Mut_NucSub, by = 'NucSubst',all.x = TRUE)
  MutSpec12Comp_heatmap[is.na(MutSpec12Comp_heatmap)] <- 0
}
write.csv(MutSpec12Comp_heatmap, "D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/new_data/data_obtained/09_MutSpec12Comp_AllMut_Genes.csv")



# ann$MutExp = 1
# ann$FromWithNeigh = paste(ann$NeighL,ann$RefNuc,ann$NeighR, sep='')
# ann$ToWithNeigh = paste(ann$NeighL,ann$AltNuc,ann$NeighR, sep='')
# ann$MutSpecWithNeigh = paste(ann$FromWithNeigh, ann$ToWithNeigh, sep = '>')
# Exp192Comp = aggregate(ann$MutExp, by = list(ann$MutSpecWithNeigh), FUN = sum);
# names(Exp192Comp) = c('NucSubst','For_all_dates')
# 
# MutSpec192Comp_heatmap = Exp192Comp$NucSubst
# MutSpec192Comp_heatmap = data.frame(MutSpec192Comp_heatmap)
# names(MutSpec192Comp_heatmap)=c('NucSubst')
# t2 = MutSpec192Comp_heatmap
# 
# for (gen in unique(nstrain$GenName)){
#   gen_df = nstrain[nstrain$GenName == gen,]
#   gen_ann = ann[ann$GenName == gen,]
#   gen_ann$MutExp = 1
#   gen_ann$FromWithNeigh = paste(gen_ann$NeighL,gen_ann$RefNuc,gen_ann$NeighR, sep='')
#   gen_ann$ToWithNeigh = paste(gen_ann$NeighL,gen_ann$AltNuc,gen_ann$NeighR, sep='')
#   gen_ann$MutSpecWithNeigh = paste(gen_ann$FromWithNeigh, gen_ann$ToWithNeigh, sep = '>')
#   Exp192Comp = aggregate(gen_ann$MutExp, by = list(gen_ann$MutSpecWithNeigh), FUN = sum);
#   names(Exp192Comp) = c('NucSubst','ExpFr')
#   
#   gen_df$MutObs = 1
#   gen_df$FromWithNeigh = paste(gen_df$NeighL,gen_df$RefNuc,gen_df$NeighR, sep='')
#   gen_df$ToWithNeigh = paste(gen_df$NeighL,gen_df$AltNuc,gen_df$NeighR, sep='')
#   gen_df$MutSpecWithNeigh = paste(gen_df$FromWithNeigh, gen_df$ToWithNeigh, sep = '>')
#   Obs192Comp = aggregate(gen_df$MutObs, by = list(gen_df$MutSpecWithNeigh), FUN = sum);
#   names(Obs192Comp) = c('NucSubst','ObsFr')
#   
#   MutSpec192Comp = merge(t2,Exp192Comp, by = 'NucSubst', all.x = TRUE)
#   MutSpec192Comp = merge(MutSpec192Comp,Obs192Comp, by = 'NucSubst', all.x = TRUE)
#   MutSpec192Comp[is.na(MutSpec192Comp)] <- 0
#   MutSpec192Comp$ObsToExp = MutSpec192Comp$ObsFr/MutSpec192Comp$ExpFr
#   MutSpec192Comp$ObsToExp[is.na(MutSpec192Comp$ObsToExp)] <- 0
#   MutSpec192Comp$MutSpec = MutSpec192Comp$ObsToExp/sum(MutSpec192Comp$ObsToExp)
#   
#   Mut_NucSub = cbind(MutSpec192Comp$NucSubst, MutSpec192Comp$MutSpec);
#   colnames(Mut_NucSub) = c('NucSubst', gen)
#   MutSpec192Comp_heatmap = merge(MutSpec192Comp_heatmap, Mut_NucSub, by = 'NucSubst',all.x = TRUE)
#   MutSpec192Comp_heatmap[is.na(MutSpec192Comp_heatmap)] <- 0
# }
# 
# write.csv(MutSpec192Comp_heatmap, "D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/data_obtained/Single_Mut_PerDate/09_MutSpec192Comp_AllMut_Genes.csv")

gen_map_12 = as.data.frame(MutSpec12Comp_heatmap, row.names = MutSpec12Comp_heatmap[[1]])
gen_map_12[[1]] <- NULL
for (i in colnames(gen_map_12)){
  gen_map_12[[i]] = as.numeric(gen_map_12[[i]])
}
#log_gen_map_12 = log(gen_map_12+1, base = 3)
#log_gen_map_12 = log_gen_map_12[c(1,3,9,8,5,11,4,10,7,2,6)]
heatmaply(
  as.matrix(gen_map_12),
  Rowv = NULL,
  Colv = NULL,
  fontsize_row = 16,
  fontsize_col = 16,
  xlab = "Genes",
  ylab = 'Mutations',
  main = 'Twelve component mutation spectrum for each gen of SARS_CoV-2',
  show_dendrogram = c(FALSE, FALSE),
  plot_method = "plotly",
  file = "D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/new_data/figures/HeatMaps/09_Mut_change_12_comp_Gen.html"
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
  fontsize_row = 15,
  fontsize_col = 15,
  xlab = "Dates",
  ylab = 'Mutations',
  main = 'Mutation Change 192 components per Gen',
  show_dendrogram = c(FALSE, FALSE),
  plot_method = "plotly",
  cexRow = 0.5,
  file = "D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/R_figure/HeatMap/09_Mut_change_192_comp_Gen.html"
)

gis = read.csv('D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/new_data/data/uppdate_mutations.csv')
gis = gis[gis$RefNuc != '-'|gis$AltNuc!='-',]
orf_df = gis[gis$GenName=='ORF1ab',]
orf_df = orf_df[orf_df$AaSub=='S',]

orf_df$nsp = ''
orf_df[orf_df$Pos>=266 & orf_df$Pos<=805,]$nsp = 'nsp1'
orf_df[orf_df$Pos>=806 & orf_df$Pos<=2719,]$nsp = 'nsp2'
orf_df[orf_df$Pos>=2720 & orf_df$Pos<=8554,]$nsp = 'nsp3'
orf_df[orf_df$Pos>=8555 & orf_df$Pos<=10054,]$nsp = 'nsp4'
orf_df[orf_df$Pos>=10055 & orf_df$Pos<=10972,]$nsp = '_3C-like proteinase'
orf_df[orf_df$Pos>=10973 & orf_df$Pos<=11842,]$nsp = 'nsp6'
orf_df[orf_df$Pos>=11843 & orf_df$Pos<=12091,]$nsp = 'nsp7'
orf_df[orf_df$Pos>=12092 & orf_df$Pos<=12685,]$nsp = 'nsp8'
orf_df[orf_df$Pos>=12686 & orf_df$Pos<=13024,]$nsp = 'nsp9'
orf_df[orf_df$Pos>=13025 & orf_df$Pos<=13441,]$nsp = 'nsp10'
orf_df[orf_df$Pos>=13442 & orf_df$Pos<=16236,]$nsp = '_RNA-dependent_RNA_polymerase'
orf_df[orf_df$Pos>=16237 & orf_df$Pos<=18039,]$nsp = 'helicase'
orf_df[orf_df$Pos>=18040 & orf_df$Pos<=19620,]$nsp = "_3'-to-5'_exonuclease"
orf_df[orf_df$Pos>=19621 & orf_df$Pos<=20658,]$nsp = "endoRNAse"
orf_df[orf_df$Pos>=20659 & orf_df$Pos<=21555,]$nsp = "_2'-O-ribose_methyltransferase"
orf_df$MutObs = 1

ideal_table = read.csv("D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/data_obtained/ideal_table.csv")
ann = ideal_table[ideal_table$GenName=='ORF1ab',]
ann = ann[ann$AaSub=='S',]
ann$nsp = ''
names(ann)[3] = 'pos'
ann[ann$pos>=266 & ann$pos<=805,]$nsp = 'nsp1'
ann[ann$pos>=806 & ann$pos<=2719,]$nsp = 'nsp2'
ann[ann$pos>=2720 & ann$pos<=8554,]$nsp = 'nsp3'
ann[ann$pos>=8555 & ann$pos<=10054,]$nsp = 'nsp4'
ann[ann$pos>=10055 & ann$pos<=10972,]$nsp = '_3C-like proteinase'
ann[ann$pos>=10973 & ann$pos<=11842,]$nsp = 'nsp6'
ann[ann$pos>=11843 & ann$pos<=12091,]$nsp = 'nsp7'
ann[ann$pos>=12092 & ann$pos<=12685,]$nsp = 'nsp8'
ann[ann$pos>=12686 & ann$pos<=13024,]$nsp = 'nsp9'
ann[ann$pos>=13025 & ann$pos<=13441,]$nsp = 'nsp10'
ann[ann$pos>=13442 & ann$pos<=16236,]$nsp = '_RNA-dependent_RNA_polymerase'
ann[ann$pos>=16237 & ann$pos<=18039,]$nsp = 'helicase'
ann[ann$pos>=18040 & ann$pos<=19620,]$nsp = "_3'-to-5'_exonuclease"
ann[ann$pos>=19621 & ann$pos<=20658,]$nsp = "endoRNAse"
ann[ann$pos>=20659 & ann$pos<=21555,]$nsp = "_2'-O-ribose_methyltransferase"

nsp_csv = data.frame(unique(orf_df$nsp))
names(nsp_csv) = 'nsp'
nsp_csv$len = 0
nsp_csv$MutObs = 0
nsp_csv$MutObs_S = 0
nsp_csv$MutObs_NS = 0
for (nsp in unique(orf_df$nsp)){
  nsp_csv[nsp_csv$nsp == nsp,]$len = length(ann[ann$nsp==nsp,]$pos / 3)
  nsp_csv[nsp_csv$nsp == nsp,]$MutObs = sum(orf_df[orf_df$nsp == nsp,]$MutObs)
  nsp_csv[nsp_csv$nsp == nsp,]$MutObs_S = sum(orf_df[orf_df$nsp == nsp & orf_df$AaRef==orf_df$AaAlt,]$MutObs)
  nsp_csv[nsp_csv$nsp == nsp,]$MutObs_NS = sum(orf_df[orf_df$nsp == nsp & orf_df$AaRef!=orf_df$AaAlt,]$MutObs)
}
nsp_csv$Kn_Ks = nsp_csv$MutObs_NS / nsp_csv$MutObs_S

MutSpec12Comp_heatmap = c('A>C','A>G','A>T','C>A','C>G','C>T','G>A','G>C','G>T','T>A','T>C','T>G')
MutSpec12Comp_heatmap = data.frame(MutSpec12Comp_heatmap)
names(MutSpec12Comp_heatmap)=c('NucSubst')
t1 = MutSpec12Comp_heatmap
for (gen in unique(orf_df$nsp)){
  gen_df = orf_df[orf_df$nsp == gen,]
  gen_ann = ann[ann$nsp == gen,]
  
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
  if (gen == 'ORF6'){
    For_Mut_Spec[For_Mut_Spec$ExpFr == 0,]$ExpFr = sum(gen_ann[gen_ann$RefNuc=='G',]$MutExp)
  }
  For_Mut_Spec[[gen]]  = For_Mut_Spec$ObsFr/For_Mut_Spec$ExpFr
  For_Mut_Spec[[gen]][is.na(For_Mut_Spec[[gen]] )] <- 0
  #For_Mut_Spec[[gen]] = For_Mut_Spec$ObsToExp/sum(For_Mut_Spec$ObsToExp)
  Mut_NucSub = cbind(For_Mut_Spec$NucSubst, For_Mut_Spec[[gen]]);
  colnames(Mut_NucSub) = c('NucSubst', gen)
  MutSpec12Comp_heatmap = merge(MutSpec12Comp_heatmap, Mut_NucSub, by = 'NucSubst',all.x = TRUE)
  MutSpec12Comp_heatmap[is.na(MutSpec12Comp_heatmap)] <- 0
}
gen_map_12 = as.data.frame(MutSpec12Comp_heatmap, row.names = MutSpec12Comp_heatmap[[1]])
gen_map_12[[1]] <- NULL
#gen_map_12 = gen_map_12[c(15,8,9,1,2,3,12,4,10,13,5,6,7,11,14)]
for (i in colnames(gen_map_12)){
  gen_map_12[[i]] = as.numeric(gen_map_12[[i]])
}
heatmaply(
  as.matrix(gen_map_12),
  Rowv = NULL,
  Colv = NULL,
  fontsize_row = 16,
  fontsize_col = 16,
  xlab = "Genes",
  ylab = 'Mutations',
  main = 'Twelve components mutation spectrum of ORF1ab gen',
  show_dendrogram = c(FALSE, FALSE),
  plot_method = "plotly",
  file = "D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/new_data/figures/HeatMaps/09_Mut_change_12_comp_orf1ab.html"
)

write.csv(nsp_csv, 'D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/new_data/data/09.nsp_len_mut.csv')
orf_df=orf_df[orf_df$parent_aa == orf_df$child_aa,]
ann = ann[ann$AaSub == 'S',]

mut_list = c('A>C','A>G','A>T','C>A','C>G','C>T','G>A','G>C','G>T','T>A','T>C','T>G')
mut_list = data.frame(mut_list)
names(mut_list)=c('NucSubst')
ann$MutExp <- 1
ann$NucSubst = paste(ann$RefNuc,ann$AltNuc,sep='>')
Exp12Comp = aggregate(ann$MutExp, by = list(ann$NucSubst), FUN = sum);
names(Exp12Comp) = c('NucSubst','ExpFr')

orf_df$Mutobs = 1
orf_df$NucSubst = paste(orf_df$parent_nucl,orf_df$child_nucl,sep='>')
Obs12Comp = aggregate(orf_df$MutObs, by = list(orf_df$NucSubst), FUN = sum);
names(Obs12Comp) = c('NucSubst','ObsFr')

MutSpec12Comp = merge(mut_list,Exp12Comp, by = 'NucSubst', all.x = TRUE)
MutSpec12Comp = merge(MutSpec12Comp,Obs12Comp, by = 'NucSubst',all.x = TRUE)
MutSpec12Comp[is.na(MutSpec12Comp)] <- 0
MutSpec12Comp$ObsToExp = MutSpec12Comp$ObsFr/MutSpec12Comp$ExpFr
MutSpec12Comp$ObsToExp[is.na(MutSpec12Comp$ObsToExp)] <- 0
MutSpec12Comp$MutSpec = MutSpec12Comp$ObsToExp/sum(MutSpec12Comp$ObsToExp)

barplot(MutSpec12Comp$MutSpec, names =  MutSpec12Comp$NucSubst, main = 'MutSpec 12 components for ORF1ab', cex.names = 0.7)