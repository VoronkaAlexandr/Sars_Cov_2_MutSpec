library(ggplot2)
library(plotly)
library(seqinr)

df = read.csv('D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/data/norm_data_modernized_nextstrain.csv')
ref = read.fasta('D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/data/Covid_ref.fasta', forceDNAtolower = FALSE)
ann = read.csv("D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/data_obtained/ideal_table.csv")
ann = ann[ann$GenType == 'translated',]
ann = ann[ann$RefAa != '*',]
ann = ann[ann$RefAa == ann$AltAa,]
ann$NucSubst = paste(ann$RefNuc,ann$AltNuc,sep='>')
ann$MutExp <-1
ann[is.na(ann$MutObsNC),]$MutObsNC <- 0
only_mut = subset(df, RefNuc != '-' & AltNuc != '-' & NucInCodon != 'None')
for (i in (1 : nrow(only_mut))){
  if (only_mut$NucInCodon[i] == '1'){
    only_mut$RefAa[i] = translate(s2c(paste(df$RefNuc[i], ref$NC_045512.2[df$Pos[i] + 1],ref$NC_045512.2[df$Pos[i] + 2], sep = '')))
    only_mut$AltAa[i] = translate(s2c(paste(df$AltNuc[i], ref$NC_045512.2[df$Pos[i] + 1],ref$NC_045512.2[df$Pos[i] + 2], sep = '')))
  }
  else if(only_mut$NucInCodon[i] == '2'){
    only_mut$RefAa[i] = translate(s2c(paste(ref$NC_045512.2[df$Pos[i] -1], df$RefNuc[i], ref$NC_045512.2[df$Pos[i] + 1], sep = '')))
    only_mut$AltAa[i] = translate(s2c(paste(ref$NC_045512.2[df$Pos[i] -1], df$AltNuc[i], ref$NC_045512.2[df$Pos[i] + 1], sep = '')))
  }
  else if(only_mut$NucInCodon[i] == '3'){
    only_mut$RefAa[i] = translate(s2c(paste(ref$NC_045512.2[df$Pos[i] -2], ref$NC_045512.2[df$Pos[i] -1], df$RefNuc[i], sep = '')))
    only_mut$AltAa[i] = translate(s2c(paste(ref$NC_045512.2[df$Pos[i] -2], ref$NC_045512.2[df$Pos[i] -1], df$AltNuc[i], sep = '')))
  }
}
only_mut = only_mut[only_mut$RefAa == only_mut$AltAa,]
only_mut$NucSubst = paste(only_mut$RefNuc,only_mut$AltNuc,sep='>')

ann$FromWithNeigh = paste(ann$NeighL,ann$RefNuc,ann$NeighR, sep='')
ann$ToWithNeigh = paste(ann$NeighL,ann$AltNuc,ann$NeighR, sep='')
ann$MutSpecWithNeigh = paste(ann$FromWithNeigh, ann$ToWithNeigh, sep = '>')
Exp192Comp = aggregate(ann$MutExp, by = list(ann$MutSpecWithNeigh), FUN = sum);
names(Exp192Comp)=c('NucSubst','ExpFr')

only_mut$RefNeigh = paste(ref$NC_045512.2[only_mut$Pos - 1], only_mut$RefNuc, ref$NC_045512.2[only_mut$Pos + 1], sep = '')
only_mut$AltNeigh = paste(ref$NC_045512.2[only_mut$Pos - 1], only_mut$AltNuc, ref$NC_045512.2[only_mut$Pos + 1], sep = '')
only_mut$MutWithNeigh = paste(only_mut$RefNeigh, only_mut$AltNeigh, sep = '>')


for (x in unique(only_mut$Year)){
  year_set = subset(only_mut, Year == x)
  for (i in unique(year_set$Month)){
    temporary = subset(year_set, Month == i)
    Month192Comp = aggregate(temporary$MutCount, by = list(temporary$MutWithNeigh), FUN = sum);
    names(Month192Comp)=c('NucSubst','Month_Mut_Count')
    
    fig <- ggplot(Month192Comp,aes(x = NucSubst, y = Month_Mut_Count, ymin = 0, fill = NucSubst))+
      geom_bar(stat = "identity", width = 0.6)+ 
      #coord_flip() + 
      xlab("Mutations") +
      ylab("Number Of Mutations")+
      theme(legend.position="none")+
      theme(axis.text.x = element_text(size=3.2, angle = 90))+
      ggtitle(sprintf('Observed Covid Mut spec 192 components at %s', paste(x,i,sep='-')))+
      #theme(axis.text.x = element_text(hjust = 3.4))+
      geom_text(aes(label=NucSubst), size = 1, angle = 90)
    fig %>% ggplotly
    ggsave(sprintf('D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/R_figure/Mut_Per_Date/06.Obs_Month_MutSpec192_%s.svg', paste(x,i,sep='-')), plot = fig, device = 'svg', limitsize = F)
    
    write.csv(x = Month192Comp, sprintf('D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/data_obtained/MutPerMonth/06.Obs_Month_MutSpec192_%s.csv', paste(x,i,sep='-')))
  }
}

for (x in unique(only_mut$Year)){
  year_set = subset(only_mut, Year == x)
  for (i in unique(year_set$Month)){
    temporary = subset(year_set, Month == i)
    Month192Comp = aggregate(temporary$MutCount, by = list(temporary$MutWithNeigh), FUN = sum);
    names(Month192Comp)=c('NucSubst','Month_Mut_Count')
    MutSpec192Comp = merge(Exp192Comp,Month192Comp, by = 'NucSubst')
    MutSpec192Comp$ObsToExp = MutSpec192Comp$Month_Mut_Count/MutSpec192Comp$ExpFr
    MutSpec192Comp$MutSpec = MutSpec192Comp$ObsToExp/sum(MutSpec192Comp$ObsToExp)
    
    fig <- ggplot(MutSpec192Comp,aes(x = NucSubst, y = MutSpec, ymin = 0, fill = NucSubst))+
      geom_bar(stat = "identity", width = 0.6)+ 
      #coord_flip() + 
      xlab("Mutations") +
      ylab("Number Of Mutations")+
      theme(legend.position="none")+
      theme(axis.text.x = element_text(size=3.2, angle = 90))+
      ggtitle(sprintf('Covid Mut spec 192 components at %s', paste(x,i,sep='-')))+
      #theme(axis.text.x = element_text(hjust = 3.4))+
      geom_text(aes(label=NucSubst), size = 1, angle = 90)
    fig %>% ggplotly
    ggsave(sprintf('D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/R_figure/Mut_Per_Date/06.Month_MutSpec192_%s.svg', paste(x,i,sep='-')), plot = fig, device = 'svg', limitsize = F)
    
    write.csv(x = Month192Comp, sprintf('D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/data_obtained/MutPerMonth/06.Month_MutSpec192_%s.csv', paste(x,i,sep='-')))
  }
}