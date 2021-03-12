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

for (gen in unique(ann$GenName)){
  gen_ann = ann[ann$GenName == gen,]
  Exp12Comp = aggregate(gen_ann$MutExp, by = list(gen_ann$MutSpecWithNeigh), FUN = sum);
  names(Exp12Comp)=c('NucSubst','ExpFr')
  fig <- ggplot(Exp12Comp,aes(x = NucSubst, y = ExpFr, ymin = 0, fill = NucSubst))+
    geom_bar(stat = "identity", width = 0.6)+ 
    coord_flip() + 
    xlab("Mutations") +
    ylab("Number Of Mutations")+
    theme(legend.position="none")+
    theme(axis.text.y = element_text(size=3.2))+
    ggtitle(sprintf('Expected Covid Mut spec 192 components in %s', gen))+
    theme(axis.text.y = element_text(hjust = 3.4))+
    geom_text(aes(label=NucSubst), size = 1)
  fig %>% ggplotly
  ggsave(sprintf('D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/R_figure/04.Expec_MutSpec192_%s.svg', gen), plot = fig, device = 'svg', limitsize = F)
}

for (gen in unique(ann$GenName)){
  gen_ann = ann[ann$GenName == gen,]
  Obs12Comp = aggregate(gen_ann$MutObsNC, by = list(gen_ann$MutSpecWithNeigh), FUN = sum);
  names(Obs12Comp)=c('NucSubst','ObsFr')
  fig <- ggplot(Obs12Comp,aes(x = NucSubst, y = ObsFr, ymin = 0, fill = NucSubst))+
    geom_bar(stat = "identity", width = 0.6)+ 
    coord_flip() + 
    xlab("Mutations") +
    ylab("Number Of Mutations")+
    theme(legend.position="none")+
    theme(axis.text.y = element_text(size=3.2))+
    ggtitle(sprintf('Observed Covid Mut spec 192 components in %s', gen))+
    theme(axis.text.y = element_text(hjust = 3.4))+
    geom_text(aes(label=NucSubst), size = 1)
  fig %>% ggplotly
  ggsave(sprintf('D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/R_figure/04.Observ_MutSpec192_%s.svg', gen), plot = fig, device = 'svg',  limitsize = F)
}