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

Exp12Comp = aggregate(ann$MutExp, by = list(ann$NucSubst), FUN = sum);
names(Exp12Comp)=c('NucSubst','ExpFr')

for (x in unique(only_mut$Year)){
  year_set = subset(only_mut, Year == x)
  for (i in unique(year_set$Month)){
    temporary = subset(year_set, Month == i)
    Month12Comp = aggregate(temporary$MutCount, by = list(temporary$NucSubst), FUN = sum);
    names(Month12Comp)=c('NucSubst','Month_Mut_Count')
    png(sprintf('D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/R_figure/Mut_Per_Date/05.Obs_Month_MutSpec12_%s.png', paste(x,i,sep='-')))
    barplot(Month12Comp$Month_Mut_Count, names =  Month12Comp$NucSubst, main = sprintf('Observed Covid Mut spec 12 components at %s', paste(x,i,sep='-')), cex.names = 0.7)
    dev.off()
    write.csv(x = Month12Comp, sprintf('D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/data_obtained/MutPerMonth/05.Obs_Month_MutSpec12_%s.csv', paste(x,i,sep='-')))
  }
}

for (x in unique(only_mut$Year)){
  year_set = subset(only_mut, Year == x)
  for (i in unique(year_set$Month)){
    temporary = subset(year_set, Month == i)
    Month12Comp = aggregate(temporary$MutCount, by = list(temporary$NucSubst), FUN = sum);
    names(Month12Comp)=c('NucSubst','Month_Mut_Count')
    MutSpec12Comp = merge(Exp12Comp,Month12Comp, by = 'NucSubst')
    MutSpec12Comp$ObsToExp = MutSpec12Comp$Month_Mut_Count/MutSpec12Comp$ExpFr
    MutSpec12Comp$MutSpec = MutSpec12Comp$ObsToExp/sum(MutSpec12Comp$ObsToExp)
    print(sum(MutSpec12Comp$MutSpec))
    png(sprintf('D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/R_figure/Mut_Per_Date/05.Month_MutSpec12_%s.png', paste(x,i,sep='-')))
    barplot(MutSpec12Comp$MutSpec, names =  Month12Comp$NucSubst, main = sprintf('Covid Mut spec 12 components at %s', paste(x,i,sep='-')), cex.names = 0.7)
    dev.off()
    write.csv(x = MutSpec12Comp, sprintf('D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/data_obtained/MutPerMonth/05.Month_MutSpec12_%s.csv', paste(x,i,sep='-')))
  }
}
