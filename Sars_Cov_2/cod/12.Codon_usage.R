library(seqinr)

covid_seq = read.fasta('D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/data/first_last_1000_mulal.fasta',forceDNAtolower = FALSE)
reference = read.csv("D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/data_obtained/ideal_table.csv")


codon_usage_begin = data.frame(Codons = character(), Number_of_codons = numeric())
codon_usage_end = data.frame(Codons = character(), Number_of_codons = numeric())

first_nuc_pos = reference[reference$NucInCodon == 1,]$Pos
first_nuc_pos = unique(first_nuc_pos)

for (seq_number in (1 : length(names(covid_seq)))){
  if (seq_number != 1){
    if (seq_number<1002){
      for (nuc_number in first_nuc_pos){
        if(covid_seq[[seq_number]][nuc_number] != '-' & covid_seq[[seq_number]][nuc_number+1] != '-' & covid_seq[[seq_number]][nuc_number+2] != '-' ){
          codon = paste(covid_seq[[seq_number]][nuc_number], covid_seq[[seq_number]][nuc_number+1], covid_seq[[seq_number]][nuc_number+2], sep = '')
          if (codon %in% codon_usage_begin$Codons){
            codon_usage_begin[codon_usage_begin$Codons == codon,]$Number_of_codons = as.numeric(codon_usage_begin[codon_usage_begin$Codons == codon,]$Number_of_codons) + 1
          }else{
            codon_usage_begin[nrow(codon_usage_begin) + 1,] = c(codon, 1)
          }
        }
      }
    }else{
      for (nuc_number in first_nuc_pos){
        if(covid_seq[[seq_number]][nuc_number] != '-' & covid_seq[[seq_number]][nuc_number+1] != '-' & covid_seq[[seq_number]][nuc_number+2] != '-'){
          codon = paste(covid_seq[[seq_number]][nuc_number], covid_seq[[seq_number]][nuc_number+1], covid_seq[[seq_number]][nuc_number+2], sep = '')
          if (codon %in% codon_usage_end$Codons){
            codon_usage_end[codon_usage_end$Codons == codon,]$Number_of_codons = as.numeric(codon_usage_end[codon_usage_end$Codons == codon,]$Number_of_codons) + 1
          }else{
            codon_usage_end[nrow(codon_usage_end) + 1,] = c(codon, 1)
          }
        }
      }
    }
  }
}

write.csv(codon_usage_begin, 'D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/data_obtained/12.Codon_Usage_begin.csv')
write.csv(codon_usage_end, 'D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/data_obtained/12.Codon_Usage_end.csv')

codon_usage_begin = read.csv('D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/data_obtained/12.Codon_Usage_begin.csv')
codon_usage_end = read.csv('D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/data_obtained/12.Codon_Usage_end.csv')
codon_usage_begin$Mean_codon_number = codon_usage_begin$Number_of_codons / 1000
codon_usage_end$Mean_codon_number = codon_usage_end$Number_of_codons / 1000
write.csv(codon_usage_begin, 'D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/data_obtained/12.Codon_Usage_begin.csv')
write.csv(codon_usage_end, 'D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/data_obtained/12.Codon_Usage_end.csv')
