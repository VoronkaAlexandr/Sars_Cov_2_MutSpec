library(ggplot2)
library(lubridate)
library(dplyr)
library(plotly)
library(data.table)
# I should count number of nucleotidef, that are different from Uhani reference
reference = read.csv("D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/data_obtained/ideal_table.csv")
nstrain = read.csv('D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/data/norm_data_modernized_nextstrain.csv')
df = subset(nstrain, RefNuc != '-' & AltNuc != '-')

different_mut = 0
all_mut = 0

for (i in (1:nrow(df))){
  all_mut = all_mut + 1
  pos_set = subset(reference, Pos == df[i, 'Pos'])
  if (df[i, "RefNuc"] != pos_set[1, 'RefNuc']){
    different_mut = different_mut + 1
  }
}

procent_of_dif_mut = different_mut/all_mut*100

write(procent_of_dif_mut, 'D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/data_obtained/10.Procent_of_dif_mut.txt')

# Now i shuld count nucleotide usage among translateble zones

for_nuc_count = subset(reference, GenType == 'translated')
for_nuc_count = for_nuc_count[for_nuc_count$NucInCodon == 3,]
for_nuc_count$NucExp = 1

nuc_usage = aggregate(for_nuc_count$NucExp, by = list(for_nuc_count$RefNuc), FUN = sum);
names(nuc_usage) = c('Nucleotides', 'Number_of_nucleotides')
nuc_usage$Number_of_nucleotides = nuc_usage$Number_of_nucleotides/3
write.csv(nuc_usage, 'D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/data_obtained/10.Ref_Nuc_usage_nuc_in_codon3.csv')

# now drawing figure for C>T and G>T per date

for (i in 1:nrow(df)){
  temporary_subset = subset(reference, Pos == df$Pos[i] & AltNuc == df$AltNuc[i])
  if (length(temporary_subset$RefNuc) != 0){
    df$GenName[i] = temporary_subset$GenName
    df$GenType[i] = temporary_subset$GenType
    df$AaSub[i] = temporary_subset$AaSub
    df$NeighL[i] = temporary_subset$NeighL
    df$NeighR[i] = temporary_subset$NeighR
  }
}
tranlated_df = df[df$GenType == 'translated',]
for (i in 1:nrow(tranlated_df)){
  pos_set = subset(reference, Pos == tranlated_df[i, 'Pos'])
  if (tranlated_df[i, "RefNuc"] != pos_set[1, 'RefNuc']){
    tranlated_df$Equally = 0
  }else{
    tranlated_df$Equally = 1
  }
}
tranlated_df = tranlated_df[tranlated_df$Equally == 1,]
tranlated_df = tranlated_df[tranlated_df$AaSub == 'S',]

MutSpec_CG_to_T = c('C>T','G>T')
MutSpec_CG_to_T = data.frame(MutSpec_CG_to_T)
names(MutSpec_CG_to_T)=c('NucSubst')

date_df <- tranlated_df %>%
  arrange(Time)
for (i in seq(501, length(as.Date(date_df$Time)), 500)){
  if ((i+500) <= length(date_df$Time)){
    time_set = date_df[i:(i+500),]
    vremia = date_df$Time[i]
    time_set$MutObs = 1
    time_set$NucSubst = paste(time_set$RefNuc,time_set$AltNuc,sep='>')
    time_set = time_set[time_set$NucSubst == 'C>T' | time_set$NucSubst == 'G>T',]
    Obs12Comp = aggregate(time_set$MutObs, by = list(time_set$NucSubst), FUN = sum);
    names(Obs12Comp) = c('NucSubst', vremia)
    
    MutSpec_CG_to_T = merge(MutSpec_CG_to_T, Obs12Comp, by = 'NucSubst',all.x = TRUE)
    MutSpec_CG_to_T[is.na(MutSpec_CG_to_T)] <- 0
  }
}
  
write.csv(x = MutSpec_CG_to_T, 'D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/data_obtained/10.Date_Mut_CG_to_T_100.csv')

molten <- melt(setDT(MutSpec_CG_to_T), id.vars = c("NucSubst"))
fig = ggplot(molten, aes(x = variable, y = as.numeric(value), group = NucSubst, colour = NucSubst)) +
  theme(axis.text.x = element_text(angle = 90))+
  xlab("Dates") +
  ylab("Number of Mutation")+
  ggtitle('Mutation change over time C to T and G to T')+
  geom_line() + 
  geom_point()
ggsave('D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/R_figure/Single_Mut_PerDate/10_MutChange_CG_to_T_100.svg', plot = fig, device = 'svg',  limitsize = F)

# count Kn/Ks per 500 mutations, sorting for date

date_df

for (i in 1:nrow(date_df)){
  set_s_ns = subset(reference, Pos == date_df[i, 'Pos'])
  if (date_df[i,"RefNuc"] == set_s_ns[1, 'RefNuc']){
    date_df[i, "AaSub"] = set_s_ns[1, 'AaSub']
  }else{
    date_df[i, "AaSub"] = 0
  }
}
only_s_ns_df = date_df[date_df$AaSub != 0,]

date_sub = c('S','NS')
date_sub = data.frame(date_sub)
names(date_sub)=c('S_or_NS')

for (i in seq(501, length(as.Date(only_s_ns_df$Time)), 500)){
  if ((i+500) <= length(only_s_ns_df$Time)){
    time_set = only_s_ns_df[i:(i+500),]
    vremia = only_s_ns_df$Time[i]
    time_set$MutObs = 1
    
    Aa_sub = aggregate(time_set$MutObs, by = list(time_set$AaSub), FUN = sum);
    names(Aa_sub) = c('S_or_NS', vremia)
    
    date_sub = merge(date_sub, Aa_sub, by = 'S_or_NS',all.x = TRUE)
    date_sub[is.na(date_sub)] <- 0
  }
}

date_sub[3,(2:ncol(date_sub))] = date_sub[1,(2:ncol(date_sub))] / date_sub[2,(2:ncol(date_sub))]
date_sub[3, 'S_or_NS'] = 'N/S'

write.csv(x = date_sub, 'D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/data_obtained/10.Kn_Ks_500mut.csv')

molten <- melt(setDT(date_sub), id.vars = c("S_or_NS"))
kn_ks = molten[molten$S_or_NS == 'N/S',]
fig = ggplot(kn_ks, aes(x = variable, y = as.numeric(value), group = S_or_NS, colour = S_or_NS)) +
  theme(axis.text.x = element_text(angle = 90))+
  xlab("Dates") +
  ylab("Number of Mutation")+
  ggtitle('Kn/Ks for 500 mutations')+
  geom_line() + 
  geom_point()
ggsave('D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/R_figure/Single_Mut_PerDate/10_Kn_Ks_500mut.svg', plot = fig, device = 'svg',  limitsize = F)
n_s = molten[molten$S_or_NS != 'N/S',]
fig = ggplot(n_s, aes(x = variable, y = as.numeric(value), group = S_or_NS, colour = S_or_NS)) +
  theme(axis.text.x = element_text(angle = 90))+
  xlab("Dates") +
  ylab("Number of Mutation")+
  ggtitle('Number of Kn and Ks for 500 mutations')+
  geom_line() + 
  geom_point()
ggsave('D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/R_figure/Single_Mut_PerDate/10_Kn_and_Ks_500mut.svg', plot = fig, device = 'svg',  limitsize = F)