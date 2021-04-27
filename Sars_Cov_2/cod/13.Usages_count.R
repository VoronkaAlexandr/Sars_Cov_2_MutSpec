library(seqinr)

reference = read.csv("D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/data_obtained/ideal_table.csv")
ref_fasta = read.fasta("D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/data/Covid_ref.fasta",forceDNAtolower = FALSE)
thousand = read.fasta('D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/data/first_last_1000_mulal.fasta',forceDNAtolower = FALSE)

# First, I will count codon usage and AA usage for referecnce of SARS-CoV-2

uhani_codon_usage = data.frame(Codons = character(), Aacids = character(), Number_of_codons = numeric())

first_nuc_pos = reference[reference$NucInCodon == 1,]$Pos
first_nuc_pos = unique(first_nuc_pos)

for (nuc_number in first_nuc_pos){
  codon = paste(ref_fasta[[1]][nuc_number], ref_fasta[[1]][nuc_number+1], ref_fasta[[1]][nuc_number+2], sep = '')
  aa = translate(s2c(codon))
  if (codon %in% uhani_codon_usage$Codons){
    uhani_codon_usage[uhani_codon_usage$Codons == codon,]$Number_of_codons = as.numeric(uhani_codon_usage[uhani_codon_usage$Codons == codon,]$Number_of_codons) + 1
  }else{
    uhani_codon_usage[nrow(uhani_codon_usage) + 1,] = c(codon, aa, 1)
  }
}

uhani_codon_usage[uhani_codon_usage$Codons == 'CTT' | uhani_codon_usage$Codons == 'CTA' | uhani_codon_usage$Codons == 'CTG' | uhani_codon_usage$Codons == 'CTC',]$Aacids = 'L_CT'
uhani_codon_usage[uhani_codon_usage$Codons == 'TTA' | uhani_codon_usage$Codons == 'TTG',]$Aacids = 'L_TT'

uhani_codon_usage[uhani_codon_usage$Codons == 'TCT' | uhani_codon_usage$Codons == 'TCA' | uhani_codon_usage$Codons == 'TCG' | uhani_codon_usage$Codons == 'TCC',]$Aacids = 'S_TC'
uhani_codon_usage[uhani_codon_usage$Codons == 'AGT' | uhani_codon_usage$Codons == 'AGC',]$Aacids = 'S_AG'

uhani_codon_usage[uhani_codon_usage$Codons == 'CGC' | uhani_codon_usage$Codons == 'CGA' | uhani_codon_usage$Codons == 'CGT' | uhani_codon_usage$Codons == 'CGG',]$Aacids = 'R_CG'
uhani_codon_usage[uhani_codon_usage$Codons == 'AGG' | uhani_codon_usage$Codons == 'AGA',]$Aacids = 'R_AG'

uhani_codon_usage[uhani_codon_usage$Codons == 'TAG' | uhani_codon_usage$Codons == 'TAA',]$Aacids = '*_TA'

uhani_codon_usage$normed_number = as.numeric(uhani_codon_usage$Number_of_codons) / sum(as.numeric(uhani_codon_usage$Number_of_codons))

uhani_codon_usage[nrow(uhani_codon_usage)+1,] = c("TGA", "*", "0", "0")
uhani_codon_usage[uhani_codon_usage$Codons == 'TGA',]$Aacids = '*_TG'

uhani_aa_usage = aggregate(as.numeric(uhani_codon_usage$normed_number), by = list(uhani_codon_usage$Aacids), FUN = sum);
names(uhani_aa_usage) = c('Aacids', 'Normed_number')

write.csv(uhani_codon_usage, 'D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/data_obtained/Usage/13.Codon_usage_NC_045512.2.csv', row.names = F)
write.csv(uhani_aa_usage, 'D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/data_obtained/Usage/13.AA_usage_NC_045512.2.csv', row.names = F)

# Secondly, I will count nucleotide usage, codon usage, aa usage for first 1000 seq and last 1000 seq to see evolution

first_thousand = read.csv('D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/data_obtained/12.Codon_Usage_begin.csv')
second_thousand = read.csv('D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/data_obtained/12.Codon_Usage_end.csv')
first_thousand[[1]] <- NULL
second_thousand[[1]] <- NULL
first_thousand[[1]] <- NULL
second_thousand[[1]] <- NULL


first_thousand$Aa = ''
for (i in first_thousand$Codons){
  first_thousand[first_thousand$Codons == i,]$Aa = translate(s2c(i))
}
second_thousand$Aa = ''
for (i in second_thousand$Codons){
  second_thousand[second_thousand$Codons == i,]$Aa = translate(s2c(i))
}

first_thousand[first_thousand$Codons == 'CTT' | first_thousand$Codons == 'CTA' | first_thousand$Codons == 'CTG' | first_thousand$Codons == 'CTC',]$Aa = 'L_CT'
first_thousand[first_thousand$Codons == 'TTA' | first_thousand$Codons == 'TTG',]$Aa = 'L_TT'

first_thousand[first_thousand$Codons == 'TCT' | first_thousand$Codons == 'TCA' | first_thousand$Codons == 'TCG' | first_thousand$Codons == 'TCC',]$Aa = 'S_TC'
first_thousand[first_thousand$Codons == 'AGT' | first_thousand$Codons == 'AGC',]$Aa = 'S_AG'

first_thousand[first_thousand$Codons == 'CGC' | first_thousand$Codons == 'CGA' | first_thousand$Codons == 'CGT' | first_thousand$Codons == 'CGG',]$Aa = 'R_CG'
first_thousand[first_thousand$Codons == 'AGG' | first_thousand$Codons == 'AGA',]$Aa = 'R_AG'

first_thousand[first_thousand$Codons == 'TAG' | first_thousand$Codons == 'TAA',]$Aa = '*_TA'
first_thousand[first_thousand$Codons == 'TGA',]$Aa = '*_TG'

##
second_thousand[second_thousand$Codons == 'CTT' | second_thousand$Codons == 'CTA' | second_thousand$Codons == 'CTG' | second_thousand$Codons == 'CTC',]$Aa = 'L_CT'
second_thousand[second_thousand$Codons == 'TTA' | second_thousand$Codons == 'TTG',]$Aa = 'L_TT'

second_thousand[second_thousand$Codons == 'TCT' | second_thousand$Codons == 'TCA' | second_thousand$Codons == 'TCG' | second_thousand$Codons == 'TCC',]$Aa = 'S_TC'
second_thousand[second_thousand$Codons == 'AGT' | second_thousand$Codons == 'AGC',]$Aa = 'S_AG'

second_thousand[second_thousand$Codons == 'CGC' | second_thousand$Codons == 'CGA' | second_thousand$Codons == 'CGT' | second_thousand$Codons == 'CGG',]$Aa = 'R_CG'
second_thousand[second_thousand$Codons == 'AGG' | second_thousand$Codons == 'AGA',]$Aa = 'R_AG'

second_thousand[second_thousand$Codons == 'TAG' | second_thousand$Codons == 'TAA',]$Aa = '*_TA'
second_thousand[second_thousand$Codons == 'TGA',]$Aa = '*_TG'


first_thousand$count = 1
aa_summ = aggregate(first_thousand$count, by = list(first_thousand$Aa), FUN = sum)
names(aa_summ) = c('aa', 'virojdenosti')
cetirejdi_virojdenie = aa_summ[aa_summ$virojdenosti >= 4,]$aa
cetirejdi_virojdenie_codons = first_thousand[first_thousand$Aa %in% cetirejdi_virojdenie,]
cetirejdi_virojdenie_codons = cetirejdi_virojdenie_codons[cetirejdi_virojdenie_codons$Codons != "AGG" & cetirejdi_virojdenie_codons$Codons != "AGA" & cetirejdi_virojdenie_codons$Codons != 'AGC' & cetirejdi_virojdenie_codons$Codons != "TTA" & cetirejdi_virojdenie_codons$Codons != 'AGT' & cetirejdi_virojdenie_codons$Codons != 'TTC',]

cetirejdi_virojdenie_codons$Nuc = substring(cetirejdi_virojdenie_codons$Codons, 3)
start_nuc_usage = aggregate(cetirejdi_virojdenie_codons$Number_of_codons, by=list(cetirejdi_virojdenie_codons$Nuc), FUN = sum)
names(start_nuc_usage) = c('Nucleotides', 'count_in_thousand_seq')
start_nuc_usage$mean_count = start_nuc_usage$count_in_thousand_seq/1000
start_nuc_usage$normed = start_nuc_usage$mean_count / sum(start_nuc_usage$mean_count)

second_thousand$count = 1
aa_summ = aggregate(second_thousand$count, by = list(second_thousand$Aa), FUN = sum)
names(aa_summ) = c('aa', 'virojdenosti')
cetirejdi_virojdenie = aa_summ[aa_summ$virojdenosti >= 4,]$aa
cetirejdi_virojdenie_codons = second_thousand[second_thousand$Aa %in% cetirejdi_virojdenie,]
cetirejdi_virojdenie_codons = cetirejdi_virojdenie_codons[cetirejdi_virojdenie_codons$Codons != "AGG" & cetirejdi_virojdenie_codons$Codons != "AGA" & cetirejdi_virojdenie_codons$Codons != 'AGC' & cetirejdi_virojdenie_codons$Codons != "TTA" & cetirejdi_virojdenie_codons$Codons != 'AGT' & cetirejdi_virojdenie_codons$Codons != 'TTC',]

cetirejdi_virojdenie_codons$Nuc = substring(cetirejdi_virojdenie_codons$Codons, 3)
end_nuc_usage = aggregate(cetirejdi_virojdenie_codons$Number_of_codons, by=list(cetirejdi_virojdenie_codons$Nuc), FUN = sum)
names(end_nuc_usage) = c('Nucleotides', 'count_in_thousand_seq')
end_nuc_usage$mean_count = end_nuc_usage$count_in_thousand_seq/1000
end_nuc_usage$normed = end_nuc_usage$mean_count / sum(end_nuc_usage$mean_count)

write.csv(start_nuc_usage, 'D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/data_obtained/Usage/13.Nuc_usage_start1000.csv', row.names = F)
write.csv(end_nuc_usage, 'D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/data_obtained/Usage/13.Nuc_usage_end1000.csv', row.names = F)
#codon usage
colnames(first_thousand)[5] = 'Normed'
first_thousand$Normed = as.numeric(first_thousand$Number_of_codons) / sum(as.numeric(first_thousand$Number_of_codons))
first_codon_usage = data.frame(first_thousand$Codons,  first_thousand$Normed)
names(first_codon_usage) = c('Codons', 'Normed_mean_count')

colnames(second_thousand)[5] = 'Normed'
second_thousand$Normed = as.numeric(second_thousand$Number_of_codons) / sum(as.numeric(second_thousand$Number_of_codons))
last_codon_usage = data.frame(second_thousand$Codons,  second_thousand$Normed)
names(last_codon_usage) = c('Codons', 'Normed_mean_count')

write.csv(first_codon_usage, 'D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/data_obtained/Usage/13.Codon_usage_start1000.csv', row.names = F)
write.csv(last_codon_usage, 'D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/data_obtained/Usage/13.Codon_usage_end1000.csv', row.names = F)
# aa usage
first_aa_usage = aggregate(first_thousand$Normed, by=list(first_thousand$Aa), FUN = sum)
names(first_aa_usage) = c('Aacids', 'Normed_mean_count')

last_aa_usage = aggregate(second_thousand$Normed, by=list(second_thousand$Aa), FUN = sum)
names(last_aa_usage) = c('Aacids', 'Normed_mean_count')

write.csv(first_aa_usage, 'D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/data_obtained/Usage/13.AA_usage_start1000.csv', row.names = F)
write.csv(last_aa_usage, 'D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/data_obtained/Usage/13.AA_usage_end1000.csv', row.names = F)
