library(plotly)
library(data.table)
library(seqinr)
library(ggplot2)

gis = read.csv('D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/data/gisaid_mutations_annotation.csv')

gis = gis[gis$parent_nucl!='-' & gis$child_nucl!='-',]

gis_aa_change = gis[gis$GenType == 'translated',]
gis_aa_change$AaSub = ifelse(gis_aa_change$child_aa == gis_aa_change$parent_aa, 'S', 'NS')
gis_aa_change = gis_aa_change[gis_aa_change$AaSub == 'NS',]

ideal_table = read.csv("D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/data_obtained/ideal_table.csv")
ideal_table = ideal_table[ideal_table$GenType == 'translated',]

gis_aa_change[gis_aa_change$child_codon == 'CTT' | gis_aa_change$child_codon == 'CTA' | gis_aa_change$child_codon == 'CTG' | gis_aa_change$child_codon == 'CTC',]$child_aa = 'L_CT'
gis_aa_change[gis_aa_change$child_codon == 'TTA' | gis_aa_change$child_codon == 'TTG',]$child_aa = 'L_TT'

gis_aa_change[gis_aa_change$child_codon == 'TCT' | gis_aa_change$child_codon == 'TCA' | gis_aa_change$child_codon == 'TCG' | gis_aa_change$child_codon == 'TCC',]$child_aa = 'S_TC'
gis_aa_change[gis_aa_change$child_codon == 'AGT' | gis_aa_change$child_codon == 'AGC',]$child_aa = 'S_AG'

gis_aa_change[gis_aa_change$child_codon == 'CGC' | gis_aa_change$child_codon == 'CGA' | gis_aa_change$child_codon == 'CGT' | gis_aa_change$child_codon == 'CGG',]$child_aa = 'R_CG'
gis_aa_change[gis_aa_change$child_codon == 'AGG' | gis_aa_change$child_codon == 'AGA',]$child_aa = 'R_AG'

gis_aa_change[gis_aa_change$child_codon == 'TAG' | gis_aa_change$child_codon == 'TAA',]$child_aa = '*_TA'
gis_aa_change[gis_aa_change$child_codon == 'TGA',]$child_aa = '*_TG'

##

gis_aa_change[gis_aa_change$parent_codon == 'CTT' | gis_aa_change$parent_codon == 'CTA' | gis_aa_change$parent_codon == 'CTG' | gis_aa_change$parent_codon == 'CTC',]$parent_aa = 'L_CT'
gis_aa_change[gis_aa_change$parent_codon == 'TTA' | gis_aa_change$parent_codon == 'TTG',]$parent_aa = 'L_TT'

gis_aa_change[gis_aa_change$parent_codon == 'TCT' | gis_aa_change$parent_codon == 'TCA' | gis_aa_change$parent_codon == 'TCG' | gis_aa_change$parent_codon == 'TCC',]$parent_aa = 'S_TC'
gis_aa_change[gis_aa_change$parent_codon == 'AGT' | gis_aa_change$parent_codon == 'AGC',]$parent_aa = 'S_AG'

gis_aa_change[gis_aa_change$parent_codon == 'CGC' | gis_aa_change$parent_codon == 'CGA' | gis_aa_change$parent_codon == 'CGT' | gis_aa_change$parent_codon == 'CGG',]$parent_aa = 'R_CG'
gis_aa_change[gis_aa_change$parent_codon == 'AGG' | gis_aa_change$parent_codon == 'AGA',]$parent_aa = 'R_AG'

gis_aa_change[gis_aa_change$parent_codon == 'TAG' | gis_aa_change$parent_codon == 'TAA',]$parent_aa = '*_TA'
gis_aa_change[gis_aa_change$parent_codon == 'TGA',]$parent_aa = '*_TG'

gis_aa_change$AaSub = ifelse(gis_aa_change$child_aa == gis_aa_change$parent_aa, 'S', 'NS')
gis_aa_change = gis_aa_change[gis_aa_change$AaSub == 'NS',]

ideal_table[ideal_table$RefCodon == 'CTT' | ideal_table$RefCodon == 'CTA' | ideal_table$RefCodon == 'CTG' | ideal_table$RefCodon == 'CTC',]$RefAa = 'L_CT'
ideal_table[ideal_table$RefCodon == 'TTA' | ideal_table$RefCodon == 'TTG',]$RefAa = 'L_TT'

ideal_table[ideal_table$RefCodon == 'TCT' | ideal_table$RefCodon == 'TCA' | ideal_table$RefCodon == 'TCG' | ideal_table$RefCodon == 'TCC',]$RefAa = 'S_TC'
ideal_table[ideal_table$RefCodon == 'AGT' | ideal_table$RefCodon == 'AGC',]$RefAa = 'S_AG'

ideal_table[ideal_table$RefCodon == 'CGC' | ideal_table$RefCodon == 'CGA' | ideal_table$RefCodon == 'CGT' | ideal_table$RefCodon == 'CGG',]$RefAa = 'R_CG'
ideal_table[ideal_table$RefCodon == 'AGG' | ideal_table$RefCodon == 'AGA',]$RefAa = 'R_AG'

ideal_table[ideal_table$RefCodon == 'TAG' | ideal_table$RefCodon == 'TAA',]$RefAa = '*_TA'

ideal_table$ExpAa = 1
Aa_Count = aggregate(ideal_table$ExpAa, by = list(ideal_table$RefAa), FUN = sum)
Aa_Count$x = Aa_Count$x / 3
names(Aa_Count) = c('Aa', 'Ref_Obs')

gis_aa_change$AaMutCount = 1
gis_aa_change = gis_aa_change[gis_aa_change$child_aa != '',]
gis_aa_change$AaMut = paste(gis_aa_change$child_aa, gis_aa_change$parent_aa, sep = '>')


Aa_Mut = aggregate(gis_aa_change$AaMutCount, by = list(gis_aa_change$AaMut), FUN = sum)
names(Aa_Mut) = c('Aa_Mutation', 'Count')
mutations = c('P>L_CT', 'P>S_TC', 'S_TC>F', 'S_TC>L_TT', 'A>S_TC', 'A>V', 'T>I', 'T>M', 'L_CT>L_TT', 'L_CT>F', 'V>L_TT', 'V>F', 'R_CG>L_CT', 'R_CG>C', 'C>F', '*_TG>L_TT', 'W>L_TT', 'G>V', 'G>W', 'G>*_TG', 'S_AG>I', 'R_AG>M', 'R_AG>I', 'H>Y', 'Q>*_TA', 'E>*_TA', 'D>Y')

Aa_Mut$Group = 0
Aa_Mut[Aa_Mut$Aa_Mutation %in% mutations,]$Group = 1

Aa_Mut$Aa_parent = ''
for (i in 1:nrow(Aa_Mut)){
  Aa_Mut$Aa_parent[i] = strsplit(Aa_Mut$Aa_Mutation[i], split=">")[[1]][1]
}
Aa_Mut$Ref_Count = 0
for (i in 1:nrow(Aa_Mut)){
  if (Aa_Mut$Aa_parent[i] %in% Aa_Count$Aa){
  Aa_Mut$Ref_Count[i] = Aa_Count[Aa_Count$Aa == Aa_Mut$Aa_parent[i],]$Ref_Obs
  }else{
    Aa_Mut$Ref_Count[i] = 0
  }
}
Aa_Mut$Speed_Mut = Aa_Mut$Count / Aa_Mut$Ref_Count
Aa_Mut = Aa_Mut[Aa_Mut$Ref_Count!= 0,]

  
fig <- ggplot(Aa_Mut,aes(x = Aa_Mutation, y = Speed_Mut, ymin = 0, fill = Group))+
  geom_bar(stat = "identity", width = 0.6)+ 
  #coord_flip() + 
  xlab("Amino Acid Mutations") +
  ylab("Speed of Mutation")+
  theme(legend.position="none")+
  theme(axis.text.x = element_text(size=3, angle = 90))+
  ggtitle('Speed of Mutation amino acids in SARS-CoV-2 Gisaid')+
  #theme(axis.text.y = element_text(hjust = 3.4))+
  geom_text(aes(label=Aa_Mutation), size = 1, angle = 90)
fig %>% ggplotly
ggsave('D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/R_function/gisaid_mutspec/figures/17.Speed_aa_mut.svg', plot = fig, device = 'svg',  limitsize = F)

only_interest = Aa_Mut[Aa_Mut$Group == 1,]

fig <- ggplot(only_interest,aes(x = Aa_Mutation, y = Speed_Mut, ymin = 0, fill = Group))+
  geom_bar(stat = "identity", width = 0.6)+ 
  #coord_flip() + 
  xlab("Amino Acid Mutations") +
  ylab("Speed of Mutation")+
  theme(legend.position="none")+
  theme(axis.text.x = element_text(size=5, angle = 90))+
  ggtitle('Speed of Mutation amino acids in SARS-CoV-2 Gisaid')+
  #theme(axis.text.y = element_text(hjust = 3.4))+
  geom_text(aes(label=Aa_Mutation), size = 3, angle = 90)
fig %>% ggplotly
ggsave('D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/R_function/gisaid_mutspec/figures/17.Speed_aa_mut_only_interest.svg', plot = fig, device = 'svg',  limitsize = F)