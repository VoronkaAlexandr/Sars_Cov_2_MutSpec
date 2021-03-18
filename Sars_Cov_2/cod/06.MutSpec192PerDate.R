library(ggplot2)
library(plotly)
library(seqinr)

df = read.csv('D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/data/norm_data_modernized_nextstrain.csv')
ref = read.fasta('D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/data/Covid_ref.fasta', forceDNAtolower = FALSE)
df = subset(df, RefNuc != '-' & AltNuc != '-')

df$RefNeigh = paste(ref$NC_045512.2[df$Pos - 1], df$RefNuc, ref$NC_045512.2[df$Pos + 1], sep = '')
df$AltNeigh = paste(ref$NC_045512.2[df$Pos - 1], df$AltNuc, ref$NC_045512.2[df$Pos + 1], sep = '')
df$MutWithNeigh = paste(df$RefNeigh, df$AltNeigh, sep = '>')

for (x in unique(df$Year)){
  year_set = subset(df, Year == x)
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
      ggtitle(sprintf('Covid Mut spec 192 components at %s', paste(x,i,sep='-')))+
      #theme(axis.text.x = element_text(hjust = 3.4))+
      geom_text(aes(label=NucSubst), size = 1, angle = 90)
    fig %>% ggplotly
    ggsave(sprintf('D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/R_figure/Mut_Per_Date/06.Month_MutSpec192_%s.svg', paste(x,i,sep='-')), plot = fig, device = 'svg', limitsize = F)
    
    write.csv(x = Month192Comp, sprintf('D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/data_obtained/MutPerMonth/06.Month_MutSpec192_%s.csv', paste(x,i,sep='-')))
  }
}