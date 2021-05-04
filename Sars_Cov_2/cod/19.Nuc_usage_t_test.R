library(seqinr)
library(ggpubr)

thousand = read.fasta('D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/data/first_last_1000_mulal.fasta',forceDNAtolower = FALSE)
ideal_table = read.csv("D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/data_obtained/ideal_table.csv")

ideal_table[ideal_table$RefCodon == 'CTT' | ideal_table$RefCodon == 'CTA' | ideal_table$RefCodon == 'CTG' | ideal_table$RefCodon == 'CTC',]$RefAa = 'L_CT'
ideal_table[ideal_table$RefCodon == 'TTA' | ideal_table$RefCodon == 'TTG',]$RefAa = 'L_TT'

ideal_table[ideal_table$RefCodon == 'TCT' | ideal_table$RefCodon == 'TCA' | ideal_table$RefCodon == 'TCG' | ideal_table$RefCodon == 'TCC',]$RefAa = 'S_TC'
ideal_table[ideal_table$RefCodon == 'AGT' | ideal_table$RefCodon == 'AGC',]$RefAa = 'S_AG'

ideal_table[ideal_table$RefCodon == 'CGC' | ideal_table$RefCodon == 'CGA' | ideal_table$RefCodon == 'CGT' | ideal_table$RefCodon == 'CGG',]$RefAa = 'R_CG'
ideal_table[ideal_table$RefCodon == 'AGG' | ideal_table$RefCodon == 'AGA',]$RefAa = 'R_AG'

ideal_table[ideal_table$RefCodon == 'TAG' | ideal_table$RefCodon == 'TAA',]$RefAa = '*_TA'
###
ideal_table[ideal_table$AltCodon == 'CTT' | ideal_table$AltCodon == 'CTA' | ideal_table$AltCodon == 'CTG' | ideal_table$AltCodon == 'CTC',]$AltAa = 'L_CT'
ideal_table[ideal_table$AltCodon == 'TTA' | ideal_table$AltCodon == 'TTG',]$AltAa = 'L_TT'

ideal_table[ideal_table$AltCodon == 'TCT' | ideal_table$AltCodon == 'TCA' | ideal_table$AltCodon == 'TCG' | ideal_table$AltCodon == 'TCC',]$AltAa = 'S_TC'
ideal_table[ideal_table$AltCodon == 'AGT' | ideal_table$AltCodon == 'AGC',]$AltAa = 'S_AG'

ideal_table[ideal_table$AltCodon == 'CGC' | ideal_table$AltCodon == 'CGA' | ideal_table$AltCodon == 'CGT' | ideal_table$AltCodon == 'CGG',]$AltAa = 'R_CG'
ideal_table[ideal_table$AltCodon == 'AGG' | ideal_table$AltCodon == 'AGA',]$AltAa = 'R_AG'

ideal_table[ideal_table$AltCodon == 'TAG' | ideal_table$AltCodon == 'TAA',]$AltAa = '*_TA'
ideal_table[ideal_table$AltCodon == 'TGA',]$AltAa = '*_TG'
###
ideal_table$AaSub = ifelse(ideal_table$RefAa == ideal_table$AltAa, 'S', 'NS')



nuc_start = data.frame(matrix(vector(), 1000, length(unique(ideal_table$RefNuc)),
                       dimnames=list(c(), unique(ideal_table$RefNuc))),
                stringsAsFactors=F)
nuc_end = data.frame(matrix(vector(), 1000, length(unique(ideal_table$RefNuc)),
                            dimnames=list(c(), unique(ideal_table$RefNuc))),
                     stringsAsFactors=F)
syn_nuc_start = data.frame(matrix(vector(), 1000, length(unique(ideal_table$RefNuc)),
                                  dimnames=list(c(), unique(ideal_table$RefNuc))),
                           stringsAsFactors=F)
syn_nuc_end = data.frame(matrix(vector(), 1000, length(unique(ideal_table$RefNuc)),
                                dimnames=list(c(), unique(ideal_table$RefNuc))),
                         stringsAsFactors=F)
ff_nuc_start = data.frame(matrix(vector(), 1000, length(unique(ideal_table$RefNuc)),
                                  dimnames=list(c(), unique(ideal_table$RefNuc))),
                           stringsAsFactors=F)
ff_nuc_end = data.frame(matrix(vector(), 1000, length(unique(ideal_table$RefNuc)),
                                dimnames=list(c(), unique(ideal_table$RefNuc))),
                         stringsAsFactors=F)
codon_start = data.frame(matrix(vector(), 1000, length(unique(ideal_table[ideal_table$GenType == 'translated',]$AltCodon)),
                                dimnames=list(c(), unique(ideal_table[ideal_table$GenType == 'translated',]$AltCodon))),
                         stringsAsFactors=F)
codon_end = data.frame(matrix(vector(), 1000, length(unique(ideal_table[ideal_table$GenType == 'translated',]$AltCodon)),
                              dimnames=list(c(), unique(ideal_table[ideal_table$GenType == 'translated',]$AltCodon))),
                       stringsAsFactors=F)
first_nuc_pos = ideal_table[ideal_table$NucInCodon == 1,]$Pos
first_nuc_pos = unique(first_nuc_pos)
syn_nuc_pos = ideal_table[ideal_table$NucInCodon == 3 & ideal_table$AaSub == 'S',]$Pos
syn_nuc_pos = unique(syn_nuc_pos)
four_fold_nuc_pos = ideal_table[ideal_table$NucInCodon == 3 & ideal_table$AaSub == 'S' & (ideal_table$RefAa=='P' | ideal_table$RefAa=='S_TC'|ideal_table$RefAa=='A'|ideal_table$RefAa=='T'|ideal_table$RefAa=='L_CT'|ideal_table$RefAa=='V'|ideal_table$RefAa=='R_CG'|ideal_table$RefAa=='G'),]$Pos
four_fold_nuc_pos = unique(four_fold_nuc_pos)

for (seq_number in (1 : length(names(thousand)))){
  print(seq_number)
  if (seq_number != 1){
    if (seq_number<1002){
      for (x in 1 : length(table(thousand[[seq_number]]))){
        nucleotide = names(table(thousand[[seq_number]])[x])
        nuc_count = table(thousand[[seq_number]])[[x]]
        if (nucleotide!='-'){
          nuc_start[,nucleotide][length(nuc_start[!is.na(nuc_start[,nucleotide]),][,nucleotide])+1] = nuc_count
        }
      }
      syn_seq = ''
      for (syn_nuc_number in syn_nuc_pos){
        syn_seq = append(syn_seq, thousand[[seq_number]][syn_nuc_number], after = length(syn_seq))
      }
      syn_seq = syn_seq[syn_seq!='']
      for (syn_nuc in 1 : length(table(syn_seq))){
        if (names(table(syn_seq)[syn_nuc]) != '-'){
          nuc_syn = names(table(syn_seq)[syn_nuc])
          nuc_syn_count = table(syn_seq)[[syn_nuc]]
          syn_nuc_start[,nuc_syn][length(syn_nuc_start[!is.na(syn_nuc_start[,nuc_syn]),][,nuc_syn])+1] = nuc_syn_count
        }
      }
      four_fold_seq = ''
      for (ff_nuc_number in four_fold_nuc_pos){
        four_fold_seq = append(four_fold_seq, thousand[[seq_number]][ff_nuc_number], after = length(four_fold_seq))
      }
      four_fold_seq = four_fold_seq[four_fold_seq!='']
      for (ff_nuc in 1 : length(table(four_fold_seq))){
        if (names(table(four_fold_seq)[ff_nuc]) != '-'){
          nuc_ff = names(table(four_fold_seq)[ff_nuc])
          nuc_ff_count = table(four_fold_seq)[[ff_nuc]]
          ff_nuc_start[,nuc_ff][length(ff_nuc_start[!is.na(ff_nuc_start[,nuc_ff]),][,nuc_ff])+1] = nuc_ff_count
        }
      }
      
      codon_string = ''
      for (nuc_number in first_nuc_pos){
        if(thousand[[seq_number]][nuc_number] != '-' & thousand[[seq_number]][nuc_number+1] != '-' & thousand[[seq_number]][nuc_number+2] != '-' ){
          codon_string = append(codon_string, paste(thousand[[seq_number]][nuc_number], thousand[[seq_number]][nuc_number+1], thousand[[seq_number]][nuc_number+2], sep = ''), after = length(codon_string))
        }
      }
      codon_string = codon_string[codon_string!='']
      for (y in 1 : length(table(codon_string))){
        codon = names(table(codon_string)[y])
        codon_count = table(codon_string)[[y]]
        codon_start[,codon][length(codon_start[!is.na(codon_start[,codon]),][,codon])+1] = codon_count
    }
    }else{
      for (x in 1 : length(table(thousand[[seq_number]]))){
        nucleotide = names(table(thousand[[seq_number]])[x])
        nuc_count = table(thousand[[seq_number]])[[x]]
        if (nucleotide!='-'){
          nuc_end[,nucleotide][length(nuc_end[!is.na(nuc_end[,nucleotide]),][,nucleotide])+1] = nuc_count
        }
      }
      syn_seq = ''
      for (syn_nuc_number in syn_nuc_pos){
        syn_seq = append(syn_seq, thousand[[seq_number]][syn_nuc_number], after = length(syn_seq))
      }
      syn_seq = syn_seq[syn_seq!='']
      for (syn_nuc in 1 : length(table(syn_seq))){
        if (names(table(syn_seq)[syn_nuc]) != '-'){
          nuc_syn = names(table(syn_seq)[syn_nuc])
          nuc_syn_count = table(syn_seq)[[syn_nuc]]
          syn_nuc_end[,nuc_syn][length(syn_nuc_end[!is.na(syn_nuc_end[,nuc_syn]),][,nuc_syn])+1] = nuc_syn_count
        }
      }
      four_fold_seq = ''
      for (ff_nuc_number in four_fold_nuc_pos){
        four_fold_seq = append(four_fold_seq, thousand[[seq_number]][ff_nuc_number], after = length(four_fold_seq))
      }
      four_fold_seq = four_fold_seq[four_fold_seq!='']
      for (ff_nuc in 1 : length(table(four_fold_seq))){
        if (names(table(four_fold_seq)[ff_nuc]) != '-'){
          nuc_ff = names(table(four_fold_seq)[ff_nuc])
          nuc_ff_count = table(four_fold_seq)[[ff_nuc]]
          ff_nuc_end[,nuc_ff][length(ff_nuc_end[!is.na(ff_nuc_end[,nuc_ff]),][,nuc_ff])+1] = nuc_ff_count
        }
      }
      codon_string = ''
      for (nuc_number in first_nuc_pos){
        if(thousand[[seq_number]][nuc_number] != '-' & thousand[[seq_number]][nuc_number+1] != '-' & thousand[[seq_number]][nuc_number+2] != '-' ){
          codon_string = append(codon_string, paste(thousand[[seq_number]][nuc_number], thousand[[seq_number]][nuc_number+1], thousand[[seq_number]][nuc_number+2], sep = ''), after = length(codon_string))
        }
      }
      codon_string = codon_string[codon_string!='']
      for (y in 1 : length(table(codon_string))){
        codon = names(table(codon_string)[y])
        codon_count = table(codon_string)[[y]]
        codon_end[,codon][length(codon_end[!is.na(codon_end[,codon]),][,codon])+1] = codon_count
      }
    }
  }
}

write.csv(nuc_start, 'D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/data_obtained/19.Nuc_start_1000.csv')
write.csv(nuc_end, 'D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/data_obtained/19.Nuc_end_1000.csv')
write.csv(codon_start, 'D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/data_obtained/19.Codons_start_1000.csv')
write.csv(codon_end, 'D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/data_obtained/19.Codons_end_1000.csv')
write.csv(syn_nuc_start, 'D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/data_obtained/19.Syn_nuc_start_1000.csv')
write.csv(syn_nuc_end, 'D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/data_obtained/19.Syn_nuc_end_1000.csv')
write.csv(ff_nuc_start, 'D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/data_obtained/19.FF_start_1000.csv')
write.csv(ff_nuc_end, 'D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/data_obtained/19.FF_end_1000.csv')
nucleotides = c('A', 'T', 'G', 'C')
nuc_usage = data.frame(nucleotides)
nuc_usage$Normed_all_start = 0
nuc_usage$Normed_all_end = 0
nuc_usage$Normed_syn_start = 0
nuc_usage$Normed_syn_end = 0
nuc_usage$Normed_ff_start = 0
nuc_usage$Normed_ff_end = 0
nuc_start$A_normed = nuc_start$A / (nuc_start[,1] + nuc_start[,2] + nuc_start[,3] + nuc_start[,4])
nuc_start$T_normed = nuc_start[,2] / (nuc_start[,1] + nuc_start[,2] + nuc_start[,3] + nuc_start[,4])
nuc_start$G_normed = nuc_start$G / (nuc_start[,1] + nuc_start[,2] + nuc_start[,3] + nuc_start[,4])
nuc_start$C_normed = nuc_start$C / (nuc_start[,1] + nuc_start[,2] + nuc_start[,3] + nuc_start[,4])
A_all_start_mean = t.test(nuc_start[,5])$estimate[[1]]
conf_a_all_start = A_all_start_mean - t.test(nuc_start[,5])$conf[[1]]
T_all_start_mean = t.test(nuc_start[,6])$estimate[[1]]
conf_t_all_start = T_all_start_mean - t.test(nuc_start[,6])$conf[[1]]
G_all_start_mean = t.test(nuc_start[,7])$estimate[[1]]
conf_g_all_start = G_all_start_mean - t.test(nuc_start[,7])$conf[[1]]
C_all_start_mean = t.test(nuc_start[,8])$estimate[[1]]
conf_c_all_start = C_all_start_mean - t.test(nuc_start[,8])$conf[[1]]
nuc_usage[nuc_usage$nucleotides == 'A',]$Normed_all_start = paste(round(A_all_start_mean,4), round(conf_a_all_start,8), sep = '+-')
nuc_usage[nuc_usage$nucleotides == 'T',]$Normed_all_start = paste(round(T_all_start_mean,4), round(conf_t_all_start,8), sep = '+-')
nuc_usage[nuc_usage$nucleotides == 'G',]$Normed_all_start = paste(round(G_all_start_mean,4), round(conf_g_all_start,8), sep = '+-')
nuc_usage[nuc_usage$nucleotides == 'C',]$Normed_all_start = paste(round(C_all_start_mean,4), round(conf_c_all_start,8), sep = '+-')
#nuc usage end
nuc_end$A_normed = nuc_end$A / (nuc_end[,1] + nuc_end[,2] + nuc_end[,3] + nuc_end[,4])
nuc_end$T_normed = nuc_end[,2] / (nuc_end[,1] + nuc_end[,2] + nuc_end[,3] + nuc_end[,4])
nuc_end$G_normed = nuc_end$G / (nuc_end[,1] + nuc_end[,2] + nuc_end[,3] + nuc_end[,4])
nuc_end$C_normed = nuc_end$C / (nuc_end[,1] + nuc_end[,2] + nuc_end[,3] + nuc_end[,4])
A_all_start_mean = t.test(nuc_end[,5])$estimate[[1]]
conf_a_all_start = A_all_start_mean - t.test(nuc_end[,5])$conf[[1]]
T_all_start_mean = t.test(nuc_end[,6])$estimate[[1]]
conf_t_all_start = T_all_start_mean - t.test(nuc_end[,6])$conf[[1]]
G_all_start_mean = t.test(nuc_end[,7])$estimate[[1]]
conf_g_all_start = G_all_start_mean - t.test(nuc_end[,7])$conf[[1]]
C_all_start_mean = t.test(nuc_end[,8])$estimate[[1]]
conf_c_all_start = C_all_start_mean - t.test(nuc_end[,8])$conf[[1]]
nuc_usage[nuc_usage$nucleotides == 'A',]$Normed_all_end = paste(round(A_all_start_mean,4), round(conf_a_all_start,8), sep = '+-')
nuc_usage[nuc_usage$nucleotides == 'T',]$Normed_all_end = paste(round(T_all_start_mean,4), round(conf_t_all_start,8), sep = '+-')
nuc_usage[nuc_usage$nucleotides == 'G',]$Normed_all_end = paste(round(G_all_start_mean,4), round(conf_g_all_start,8), sep = '+-')
nuc_usage[nuc_usage$nucleotides == 'C',]$Normed_all_end = paste(round(C_all_start_mean,4), round(conf_c_all_start,8), sep = '+-')
#nuc syn start
syn_nuc_start$A_normed = syn_nuc_start$A / (syn_nuc_start[,1] + syn_nuc_start[,2] + syn_nuc_start[,3] + syn_nuc_start[,4])
syn_nuc_start$T_normed = syn_nuc_start[,2] / (syn_nuc_start[,1] + syn_nuc_start[,2] + syn_nuc_start[,3] + syn_nuc_start[,4])
syn_nuc_start$G_normed = syn_nuc_start$G / (syn_nuc_start[,1] + syn_nuc_start[,2] + syn_nuc_start[,3] + syn_nuc_start[,4])
syn_nuc_start$C_normed = syn_nuc_start$C / (syn_nuc_start[,1] + syn_nuc_start[,2] + syn_nuc_start[,3] + syn_nuc_start[,4])
A_all_start_mean = t.test(syn_nuc_start[,5])$estimate[[1]]
conf_a_all_start = A_all_start_mean - t.test(syn_nuc_start[,5])$conf[[1]]
T_all_start_mean = t.test(syn_nuc_start[,6])$estimate[[1]]
conf_t_all_start = T_all_start_mean - t.test(syn_nuc_start[,6])$conf[[1]]
G_all_start_mean = t.test(syn_nuc_start[,7])$estimate[[1]]
conf_g_all_start = G_all_start_mean - t.test(syn_nuc_start[,7])$conf[[1]]
C_all_start_mean = t.test(syn_nuc_start[,8])$estimate[[1]]
conf_c_all_start = C_all_start_mean - t.test(syn_nuc_start[,8])$conf[[1]]
nuc_usage[nuc_usage$nucleotides == 'A',]$Normed_syn_start = paste(round(A_all_start_mean,4), round(conf_a_all_start,8), sep = '+-')
nuc_usage[nuc_usage$nucleotides == 'T',]$Normed_syn_start = paste(round(T_all_start_mean,4), round(conf_t_all_start,8), sep = '+-')
nuc_usage[nuc_usage$nucleotides == 'G',]$Normed_syn_start = paste(round(G_all_start_mean,4), round(conf_g_all_start,8), sep = '+-')
nuc_usage[nuc_usage$nucleotides == 'C',]$Normed_syn_start = paste(round(C_all_start_mean,4), round(conf_c_all_start,8), sep = '+-')
#nuc syn end
syn_nuc_end$A_normed = syn_nuc_end$A / (syn_nuc_end[,1] + syn_nuc_end[,2] + syn_nuc_end[,3] + syn_nuc_end[,4])
syn_nuc_end$T_normed = syn_nuc_end[,2] / (syn_nuc_end[,1] + syn_nuc_end[,2] + syn_nuc_end[,3] + syn_nuc_end[,4])
syn_nuc_end$G_normed = syn_nuc_end$G / (syn_nuc_end[,1] + syn_nuc_end[,2] + syn_nuc_end[,3] + syn_nuc_end[,4])
syn_nuc_end$C_normed = syn_nuc_end$C / (syn_nuc_end[,1] + syn_nuc_end[,2] + syn_nuc_end[,3] + syn_nuc_end[,4])
A_all_start_mean = t.test(syn_nuc_end[,5])$estimate[[1]]
conf_a_all_start = A_all_start_mean - t.test(syn_nuc_end[,5])$conf[[1]]
T_all_start_mean = t.test(syn_nuc_end[,6])$estimate[[1]]
conf_t_all_start = T_all_start_mean - t.test(syn_nuc_end[,6])$conf[[1]]
G_all_start_mean = t.test(syn_nuc_end[,7])$estimate[[1]]
conf_g_all_start = G_all_start_mean - t.test(syn_nuc_end[,7])$conf[[1]]
C_all_start_mean = t.test(syn_nuc_end[,8])$estimate[[1]]
conf_c_all_start = C_all_start_mean - t.test(syn_nuc_end[,8])$conf[[1]]
nuc_usage[nuc_usage$nucleotides == 'A',]$Normed_syn_end = paste(round(A_all_start_mean,4), round(conf_a_all_start,8), sep = '+-')
nuc_usage[nuc_usage$nucleotides == 'T',]$Normed_syn_end = paste(round(T_all_start_mean,4), round(conf_t_all_start,8), sep = '+-')
nuc_usage[nuc_usage$nucleotides == 'G',]$Normed_syn_end = paste(round(G_all_start_mean,4), round(conf_g_all_start,8), sep = '+-')
nuc_usage[nuc_usage$nucleotides == 'C',]$Normed_syn_end = paste(round(C_all_start_mean,4), round(conf_c_all_start,8), sep = '+-')
#nuc ff start
ff_nuc_start$A_normed = ff_nuc_start$A / (ff_nuc_start[,1] + ff_nuc_start[,2] + ff_nuc_start[,3] + ff_nuc_start[,4])
ff_nuc_start$T_normed = ff_nuc_start[,2] / (ff_nuc_start[,1] + ff_nuc_start[,2] + ff_nuc_start[,3] + ff_nuc_start[,4])
ff_nuc_start$G_normed = ff_nuc_start$G / (ff_nuc_start[,1] + ff_nuc_start[,2] + ff_nuc_start[,3] + ff_nuc_start[,4])
ff_nuc_start$C_normed = ff_nuc_start$C / (ff_nuc_start[,1] + ff_nuc_start[,2] + ff_nuc_start[,3] + ff_nuc_start[,4])
A_all_start_mean = t.test(ff_nuc_start[,5])$estimate[[1]]
conf_a_all_start = A_all_start_mean - t.test(ff_nuc_start[,5])$conf[[1]]
T_all_start_mean = t.test(ff_nuc_start[,6])$estimate[[1]]
conf_t_all_start = T_all_start_mean - t.test(ff_nuc_start[,6])$conf[[1]]
G_all_start_mean = t.test(ff_nuc_start[,7])$estimate[[1]]
conf_g_all_start = G_all_start_mean - t.test(ff_nuc_start[,7])$conf[[1]]
C_all_start_mean = t.test(ff_nuc_start[,8])$estimate[[1]]
conf_c_all_start = C_all_start_mean - t.test(ff_nuc_start[,8])$conf[[1]]
nuc_usage[nuc_usage$nucleotides == 'A',]$Normed_ff_start = paste(round(A_all_start_mean,4), round(conf_a_all_start,8), sep = '+-')
nuc_usage[nuc_usage$nucleotides == 'T',]$Normed_ff_start = paste(round(T_all_start_mean,4), round(conf_t_all_start,8), sep = '+-')
nuc_usage[nuc_usage$nucleotides == 'G',]$Normed_ff_start = paste(round(G_all_start_mean,4), round(conf_g_all_start,8), sep = '+-')
nuc_usage[nuc_usage$nucleotides == 'C',]$Normed_ff_start = paste(round(C_all_start_mean,4), round(conf_c_all_start,8), sep = '+-')
#nuc ff end
ff_nuc_end$A_normed = ff_nuc_end$A / (ff_nuc_end[,1] + ff_nuc_end[,2] + ff_nuc_end[,3] + ff_nuc_end[,4])
ff_nuc_end$T_normed = ff_nuc_end[,2] / (ff_nuc_end[,1] + ff_nuc_end[,2] + ff_nuc_end[,3] + ff_nuc_end[,4])
ff_nuc_end$G_normed = ff_nuc_end$G / (ff_nuc_end[,1] + ff_nuc_end[,2] + ff_nuc_end[,3] + ff_nuc_end[,4])
ff_nuc_end$C_normed = ff_nuc_end$C / (ff_nuc_end[,1] + ff_nuc_end[,2] + ff_nuc_end[,3] + ff_nuc_end[,4])
A_all_start_mean = t.test(ff_nuc_end[,5])$estimate[[1]]
conf_a_all_start = A_all_start_mean - t.test(ff_nuc_end[,5])$conf[[1]]
T_all_start_mean = t.test(ff_nuc_end[,6])$estimate[[1]]
conf_t_all_start = T_all_start_mean - t.test(ff_nuc_end[,6])$conf[[1]]
G_all_start_mean = t.test(ff_nuc_end[,7])$estimate[[1]]
conf_g_all_start = G_all_start_mean - t.test(ff_nuc_end[,7])$conf[[1]]
C_all_start_mean = t.test(ff_nuc_end[,8])$estimate[[1]]
conf_c_all_start = C_all_start_mean - t.test(ff_nuc_end[,8])$conf[[1]]
nuc_usage[nuc_usage$nucleotides == 'A',]$Normed_ff_end = paste(round(A_all_start_mean,4), round(conf_a_all_start,8), sep = '+-')
nuc_usage[nuc_usage$nucleotides == 'T',]$Normed_ff_end = paste(round(T_all_start_mean,4), round(conf_t_all_start,8), sep = '+-')
nuc_usage[nuc_usage$nucleotides == 'G',]$Normed_ff_end = paste(round(G_all_start_mean,4), round(conf_g_all_start,8), sep = '+-')
nuc_usage[nuc_usage$nucleotides == 'C',]$Normed_ff_end = paste(round(C_all_start_mean,4), round(conf_c_all_start,8), sep = '+-')

t.test(nuc_start$A_normed, nuc_end$A_normed)$p.value

png('D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/R_figure/boxplot/19.All_A.png')
boxplot(nuc_start$A_normed, nuc_end$A_normed, names = c('A_normed_start', 'A_normed_end'), main = sprintf('Nucleotide "A" for full genome, p.value = %s', t.test(nuc_start$A_normed, nuc_end$A_normed, paired=T)$p.value))
dev.off()
png('D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/R_figure/boxplot/19.All_T.png')
boxplot(nuc_start$T_normed, nuc_end$T_normed, names = c('T_normed_start', 'T_normed_end'), main = sprintf('Nucleotide "T" for full genome, p.value = %s', t.test(nuc_start$T_normed, nuc_end$T_normed, paired=T)$p.value))
dev.off()
png('D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/R_figure/boxplot/19.All_G.png')
boxplot(nuc_start$G_normed, nuc_end$G_normed, names = c('G_normed_start', 'G_normed_end'), main = sprintf('Nucleotide "G" for full genome, p.value = %s', t.test(nuc_start$G_normed, nuc_end$G_normed, paired=T)$p.value))
dev.off()
png('D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/R_figure/boxplot/19.All_C.png')
boxplot(nuc_start$C_normed, nuc_end$C_normed, names = c('C_normed_start', 'C_normed_end'), main = sprintf('Nucleotide "C" for full genome, p.value = %s', t.test(nuc_start$C_normed, nuc_end$C_normed, paired=T)$p.value))
dev.off()
# Syn
png('D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/R_figure/boxplot/19.Syn_A.png')
boxplot(syn_nuc_start$A_normed, syn_nuc_end$A_normed, names = c('A_normed_start', 'A_normed_end'), main = sprintf('Nucleotide "A" for syn position, p.value = %s', t.test(syn_nuc_start$A_normed, syn_nuc_end$A_normed, paired=T)$p.value))
dev.off()
png('D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/R_figure/boxplot/19.Syn_T.png')
boxplot(syn_nuc_start$T_normed, syn_nuc_end$T_normed, names = c('T_normed_start', 'T_normed_end'), main = sprintf('Nucleotide "T" for syn position, p.value = %s', t.test(syn_nuc_start$T_normed, syn_nuc_end$T_normed, paired=T)$p.value))
dev.off()
png('D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/R_figure/boxplot/19.Syn_G.png')
boxplot(syn_nuc_start$G_normed, syn_nuc_end$G_normed, names = c('G_normed_start', 'G_normed_end'), main = sprintf('Nucleotide "G" for syn position, p.value = %s', t.test(syn_nuc_start$G_normed, syn_nuc_end$G_normed, paired=T)$p.value))
dev.off()
png('D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/R_figure/boxplot/19.Syn_C.png')
boxplot(syn_nuc_start$C_normed, syn_nuc_end$C_normed, names = c('C_normed_start', 'C_normed_end'), main = sprintf('Nucleotide "C" for syn position, p.value = %s', t.test(syn_nuc_start$C_normed, syn_nuc_end$C_normed, paired=T)$p.value))
dev.off()
# Four Fold
png('D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/R_figure/boxplot/19.FF_A.png')
boxplot(ff_nuc_start$A_normed, ff_nuc_end$A_normed, names = c('A_normed_start', 'A_normed_end'), main = sprintf('Nucleotide "A" for ff position, p.value = %s', t.test(ff_nuc_start$A_normed, ff_nuc_end$A_normed, paired=T)$p.value))
dev.off()
png('D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/R_figure/boxplot/19.FF_T.png')
boxplot(ff_nuc_start$T_normed, ff_nuc_end$T_normed, names = c('T_normed_start', 'T_normed_end'), main = sprintf('Nucleotide "T" for ff position, p.value = %s', t.test(ff_nuc_start$T_normed, ff_nuc_end$T_normed, paired=T)$p.value))
dev.off()
png('D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/R_figure/boxplot/19.FF_G.png')
boxplot(ff_nuc_start$G_normed, ff_nuc_end$G_normed, names = c('G_normed_start', 'G_normed_end'), main = sprintf('Nucleotide "G" for ff position, p.value = %s', t.test(ff_nuc_start$G_normed, ff_nuc_end$G_normed, paired=T)$p.value))
dev.off()
png('D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/R_figure/boxplot/19.FF_C.png')
boxplot(ff_nuc_start$C_normed, ff_nuc_end$C_normed, names = c('C_normed_start', 'C_normed_end'), main = sprintf('Nucleotide "C" for ff position, p.value = %s', t.test(ff_nuc_start$C_normed, ff_nuc_end$C_normed, paired=T)$p.value))
dev.off()