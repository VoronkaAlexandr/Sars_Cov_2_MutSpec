library(ggplot2)
library(plotly)
rm(list=ls(all=TRUE))

ann = read.csv("D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/data_obtained/ideal_table.csv")
ann$MutExp <-1
ann[is.na(ann$MutObsNC),]$MutObsNC <- 0
ann = ann[ann$GenType == 'translated',]
ann = ann[ann$RefAa != '*',]
ann = ann[ann$RefAa == ann$AltAa,]

ann$FromWithNeigh = paste(ann$NeighL,ann$RefNuc,ann$NeighR, sep='')
ann$ToWithNeigh = paste(ann$NeighL,ann$AltNuc,ann$NeighR, sep='')
ann$MutSpecWithNeigh = paste(ann$FromWithNeigh, ann$ToWithNeigh, sep = '>')

Exp192Comp = aggregate(ann$MutExp, by = list(ann$MutSpecWithNeigh), FUN = sum);
names(Exp192Comp)=c('NucSubst','ExpFr')
pdf("D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/my_plot.")
fig <- plot_ly(Exp192Comp, x = ~ExpFr, y = ~NucSubst, type = 'bar', orientation = 'h', name = 'MutSpec',
               marker = list(color = 'rgb(55, 83, 109)')) %>% layout(xaxis = list(
                   size = 1,
                   color = 'rgb(107, 107, 107)'))
fig
dev.off()
for (gen in unique(ann$GenName)){
  gen_ann = ann[ann$GenName == gen,]
  Exp12Comp = aggregate(gen_ann$MutExp, by = list(gen_ann$MutSpecWithNeigh), FUN = sum);
  names(Exp12Comp)=c('NucSubst','ExpFr')
  png(sprintf('D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/R_figure/04.Expec_MutSpec192_%s.png', gen))
  barplot(Exp12Comp$ExpFr, names =  Exp12Comp$NucSubst, main = sprintf('Expected Covid Mut spec 192 components in %s', gen), cex.names = 0.7, las=2, space=1)
  dev.off()
  print(lengths(Exp12Comp))
  print(gen)
}

for (gen in unique(ann$GenName)){
  gen_ann = ann[ann$GenName == gen,]
  Obs12Comp = aggregate(gen_ann$MutObsNC, by = list(gen_ann$MutSpecWithNeigh), FUN = sum);
  names(Obs12Comp)=c('NucSubst','ObsFr')
  png(sprintf('D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/R_figure/04.Observ_MutSpec192_%s.png', gen))
  barplot(Obs12Comp$ObsFr, names =  Obs12Comp$NucSubst, main = sprintf('Observed Covid Mut spec 192 components in %s', gen), cex.names = 0.7,las=2, space=1)
  dev.off()
  print(lengths(Obs12Comp))
}