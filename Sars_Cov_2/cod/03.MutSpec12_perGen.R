rm(list=ls(all=TRUE))

ann = read.csv("D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/data_obtained/ideal_table.csv")

ann = ann[ann$GenType == 'translated',]
ann = ann[ann$RefAa != '*',]
ann = ann[ann$RefAa == ann$AltAa,]
ann$NucSubst = paste(ann$RefNuc,ann$AltNuc,sep='>')
ann$MutExp <-1
ann[is.na(ann$MutObsNC),]$MutObsNC <- 0

for (gen in unique(ann$GenName)){
  gen_ann = ann[ann$GenName == gen,]
  Exp12Comp = aggregate(gen_ann$MutExp, by = list(gen_ann$NucSubst), FUN = sum);
  names(Exp12Comp)=c('NucSubst','ExpFr')
  png(sprintf('D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/R_figure/03.Expec_MutSpec12_%s.png', gen))
  barplot(Exp12Comp$ExpFr, names =  Exp12Comp$NucSubst, main = sprintf('Expected Covid Mut spec 12 components in %s', gen), cex.names = 0.7)
  dev.off()
}

for (gen in unique(ann$GenName)){
  gen_ann = ann[ann$GenName == gen,]
  Obs12Comp = aggregate(gen_ann$MutObsNC, by = list(gen_ann$NucSubst), FUN = sum);
  names(Obs12Comp)=c('NucSubst','ObsFr')
  png(sprintf('D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/R_figure/03.Observ_MutSpec12_%s.png', gen))
  barplot(Obs12Comp$ObsFr, names =  Obs12Comp$NucSubst, main = sprintf('Observed Covid Mut spec 12 components in %s', gen), cex.names = 0.7)
  dev.off()
}