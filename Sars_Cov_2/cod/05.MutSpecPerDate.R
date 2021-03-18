df = read.csv('D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/data/norm_data_modernized_nextstrain.csv')

only_mut = subset(df, RefNuc != '-' & AltNuc != '-')
only_mut$NucSubst = paste(only_mut$RefNuc,only_mut$AltNuc,sep='>')


for (x in unique(only_mut$Year)){
  year_set = subset(only_mut, Year == x)
  for (i in unique(year_set$Month)){
    temporary = subset(year_set, Month == i)
    Month12Comp = aggregate(temporary$MutCount, by = list(temporary$NucSubst), FUN = sum);
    names(Month12Comp)=c('NucSubst','Month_Mut_Count')
    png(sprintf('D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/R_figure/Mut_Per_Date/05.Month_MutSpec12_%s.png', paste(x,i,sep='-')))
    barplot(Month12Comp$Month_Mut_Count, names =  Month12Comp$NucSubst, main = sprintf('Covid Mut spec 12 components at %s', paste(x,i,sep='-')), cex.names = 0.7)
    dev.off()
  }
}

table(df$Continent)