library(seqinr)
library(ggpubr)
library(ggplot2)
library(plotly)

start = read.fasta('D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/new_data/data/splits/mulal_split_1.fasta',forceDNAtolower = FALSE)
inter = read.fasta('D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/new_data/data/splits/mulal_split_2.fasta',forceDNAtolower = FALSE)
end = read.fasta('D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/new_data/data/splits/mulal_split_3.fasta',forceDNAtolower = FALSE)



ideal_table = read.csv("D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/data_obtained/ideal_table.csv")
mut_df = read.csv("D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/new_data/data/U_mutations.csv")

###
ideal_table$AaSub = ifelse(ideal_table$RefAa == ideal_table$AltAa, 'S', 'NS')
ideal_table[ideal_table$GenType != 'translated',]$AaSub = ''
mut_df$AaSub = ifelse(mut_df$RefAa == mut_df$AltAa, 'S', 'NS')
mut_df[mut_df$GenType != 'translated',]$AaSub = ''

nuc = c('A', 'T', 'G', 'C')
nuc_start = data.frame(matrix(vector(), length(start), length(nuc),
                       dimnames=list(c(), nuc)),
                stringsAsFactors=F)
nuc_inter = data.frame(matrix(vector(), length(inter), length(nuc),
                              dimnames=list(c(), nuc)),
                       stringsAsFactors=F)
nuc_end = data.frame(matrix(vector(), length(end), length(nuc),
                            dimnames=list(c(), nuc)),
                     stringsAsFactors=F)

ff_nuc_start = data.frame(matrix(vector(), length(start), length(nuc),
                                 dimnames=list(c(), nuc)),
                          stringsAsFactors=F)
ff_nuc_inter = data.frame(matrix(vector(), length(inter), length(nuc),
                                 dimnames=list(c(), nuc)),
                          stringsAsFactors=F)
ff_nuc_end = data.frame(matrix(vector(), length(end), length(nuc),
                               dimnames=list(c(), nuc)),
                        stringsAsFactors=F)
codon_start = data.frame(matrix(vector(), length(start), length(unique(ideal_table[ideal_table$GenType == 'translated',]$AltCodon)),
                                dimnames=list(c(), unique(ideal_table[ideal_table$GenType == 'translated',]$AltCodon))),
                         stringsAsFactors=F)
codon_inter = data.frame(matrix(vector(), length(inter), length(unique(ideal_table[ideal_table$GenType == 'translated',]$AltCodon)),
                              dimnames=list(c(), unique(ideal_table[ideal_table$GenType == 'translated',]$AltCodon))),
                       stringsAsFactors=F)
codon_end = data.frame(matrix(vector(), length(end), length(unique(ideal_table[ideal_table$GenType == 'translated',]$AltCodon)),
                              dimnames=list(c(), unique(ideal_table[ideal_table$GenType == 'translated',]$AltCodon))),
                       stringsAsFactors=F)

first_nuc_pos = ideal_table[ideal_table$NucInCodon == 1,]$Pos
first_nuc_pos = unique(first_nuc_pos)
ff_aa = c('S', 'P', 'S_UC', 'A', 'T', 'L_CU', 'V', 'R_CG')
four_fold_nuc_pos = ideal_table[ideal_table$NucInCodon == 3 & ideal_table$RefAa %in% ff_aa,]$Pos
four_fold_nuc_pos = unique(four_fold_nuc_pos)
length(table(start[[2004]]))
for (seq_number in (1 : length(names(start)))){
  print(seq_number)
  for (x in 1 : length(table(start[[seq_number]]))){
    nucleotide = names(table(start[[seq_number]])[x])
    nuc_count = table(start[[seq_number]])[[x]]
    if (nucleotide!='-'){
      nuc_start[,nucleotide][seq_number] = nuc_count
    }
    four_fold_seq = ''
    for (ff_nuc_number in four_fold_nuc_pos){
      four_fold_seq = append(four_fold_seq, start[[seq_number]][ff_nuc_number], after = length(four_fold_seq))
    }
    four_fold_seq = four_fold_seq[four_fold_seq!='']
    for (ff_nuc in 1 : length(table(four_fold_seq))){
      if (names(table(four_fold_seq)[ff_nuc]) != '-'){
        nuc_ff = names(table(four_fold_seq)[ff_nuc])
        nuc_ff_count = table(four_fold_seq)[[ff_nuc]]
        ff_nuc_start[,nuc_ff][seq_number] = nuc_ff_count
      }
    }
    codon_string = ''
    for (nuc_number in first_nuc_pos){
      if(start[[seq_number]][nuc_number] != '-' & start[[seq_number]][nuc_number+1] != '-' & start[[seq_number]][nuc_number+2] != '-' ){
        codon_string = append(codon_string, paste(start[[seq_number]][nuc_number], start[[seq_number]][nuc_number+1], start[[seq_number]][nuc_number+2], sep = ''), after = length(codon_string))
      }
    }
    codon_string = codon_string[codon_string!='']
    for (y in 1 : length(table(codon_string))){
      codon = names(table(codon_string)[y])
      codon_count = table(codon_string)[[y]]
      codon_start[,codon][seq_number] = codon_count
    }
  }
}
write.csv(nuc_start, 'D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/new_data/data_obtained/19.Nuc_start.csv')
write.csv(ff_nuc_start, 'D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/new_data/data_obtained/19.FF_Nuc_Start.csv')
write.csv(codon_start, 'D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/new_data/data_obtained/19.Codon_Start.csv')
nuc_start
#Inter
for (seq_number in (1 : length(names(inter)))){
  print(seq_number)
  for (x in 1 : length(table(inter[[seq_number]]))){
    nucleotide = names(table(inter[[seq_number]])[x])
    nuc_count = table(inter[[seq_number]])[[x]]
    if (nucleotide!='-'){
      nuc_inter[,nucleotide][seq_number] = nuc_count
    }
    four_fold_seq = ''
    for (ff_nuc_number in four_fold_nuc_pos){
      four_fold_seq = append(four_fold_seq, inter[[seq_number]][ff_nuc_number], after = length(four_fold_seq))
    }
    four_fold_seq = four_fold_seq[four_fold_seq!='']
    for (ff_nuc in 1 : length(table(four_fold_seq))){
      if (names(table(four_fold_seq)[ff_nuc]) != '-'){
        nuc_ff = names(table(four_fold_seq)[ff_nuc])
        nuc_ff_count = table(four_fold_seq)[[ff_nuc]]
        ff_nuc_inter[,nuc_ff][seq_number] = nuc_ff_count
      }
    }
    codon_string = ''
    for (nuc_number in first_nuc_pos){
      if(inter[[seq_number]][nuc_number] != '-' & inter[[seq_number]][nuc_number+1] != '-' & inter[[seq_number]][nuc_number+2] != '-' ){
        codon_string = append(codon_string, paste(inter[[seq_number]][nuc_number], inter[[seq_number]][nuc_number+1], inter[[seq_number]][nuc_number+2], sep = ''), after = length(codon_string))
      }
    }
    codon_string = codon_string[codon_string!='']
    for (y in 1 : length(table(codon_string))){
      codon = names(table(codon_string)[y])
      codon_count = table(codon_string)[[y]]
      codon_inter[,codon][seq_number] = codon_count
    }
  }
}
write.csv(nuc_inter, 'D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/new_data/data_obtained/19.Nuc_inter.csv')
write.csv(ff_nuc_inter, 'D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/new_data/data_obtained/19.FF_Nuc_inter.csv')
write.csv(codon_inter, 'D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/new_data/data_obtained/19.Codon_inter.csv')

#End
for (seq_number in (1 : length(names(end)))){
  print(seq_number)
  for (x in 1 : length(table(end[[seq_number]]))){
    nucleotide = names(table(end[[seq_number]])[x])
    nuc_count = table(end[[seq_number]])[[x]]
    if (nucleotide!='-'){
      nuc_end[,nucleotide][seq_number] = nuc_count
    }
    four_fold_seq = ''
    for (ff_nuc_number in four_fold_nuc_pos){
      four_fold_seq = append(four_fold_seq, end[[seq_number]][ff_nuc_number], after = length(four_fold_seq))
    }
    four_fold_seq = four_fold_seq[four_fold_seq!='']
    for (ff_nuc in 1 : length(table(four_fold_seq))){
      if (names(table(four_fold_seq)[ff_nuc]) != '-'){
        nuc_ff = names(table(four_fold_seq)[ff_nuc])
        nuc_ff_count = table(four_fold_seq)[[ff_nuc]]
        ff_nuc_end[,nuc_ff][seq_number] = nuc_ff_count
      }
    }
    codon_string = ''
    for (nuc_number in first_nuc_pos){
      if(end[[seq_number]][nuc_number] != '-' & end[[seq_number]][nuc_number+1] != '-' & end[[seq_number]][nuc_number+2] != '-' ){
        codon_string = append(codon_string, paste(end[[seq_number]][nuc_number], end[[seq_number]][nuc_number+1], end[[seq_number]][nuc_number+2], sep = ''), after = length(codon_string))
      }
    }
    codon_string = codon_string[codon_string!='']
    for (y in 1 : length(table(codon_string))){
      codon = names(table(codon_string)[y])
      codon_count = table(codon_string)[[y]]
      codon_end[,codon][seq_number] = codon_count
    }
  }
}
write.csv(nuc_end, 'D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/new_data/data_obtained/19.Nuc_end.csv')
write.csv(ff_nuc_end, 'D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/new_data/data_obtained/19.FF_Nuc_end.csv')
write.csv(codon_end, 'D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/new_data/data_obtained/19.Codon_end.csv')


# for (seq_number in (1 : length(names(thousand)))){
#   print(seq_number)
#   if (seq_number != 1){
#     if (seq_number<1002){
#       for (x in 1 : length(table(thousand[[seq_number]]))){
#         nucleotide = names(table(thousand[[seq_number]])[x])
#         nuc_count = table(thousand[[seq_number]])[[x]]
#         if (nucleotide!='-'){
#           nuc_start[,nucleotide][length(nuc_start[!is.na(nuc_start[,nucleotide]),][,nucleotide])+1] = nuc_count
#         }
#       }
#       syn_seq = ''
#       for (syn_nuc_number in syn_nuc_pos){
#         syn_seq = append(syn_seq, thousand[[seq_number]][syn_nuc_number], after = length(syn_seq))
#       }
#       syn_seq = syn_seq[syn_seq!='']
#       for (syn_nuc in 1 : length(table(syn_seq))){
#         if (names(table(syn_seq)[syn_nuc]) != '-'){
#           nuc_syn = names(table(syn_seq)[syn_nuc])
#           nuc_syn_count = table(syn_seq)[[syn_nuc]]
#           syn_nuc_start[,nuc_syn][length(syn_nuc_start[!is.na(syn_nuc_start[,nuc_syn]),][,nuc_syn])+1] = nuc_syn_count
#         }
#       }
#       four_fold_seq = ''
#       for (ff_nuc_number in four_fold_nuc_pos){
#         four_fold_seq = append(four_fold_seq, thousand[[seq_number]][ff_nuc_number], after = length(four_fold_seq))
#       }
#       four_fold_seq = four_fold_seq[four_fold_seq!='']
#       for (ff_nuc in 1 : length(table(four_fold_seq))){
#         if (names(table(four_fold_seq)[ff_nuc]) != '-'){
#           nuc_ff = names(table(four_fold_seq)[ff_nuc])
#           nuc_ff_count = table(four_fold_seq)[[ff_nuc]]
#           ff_nuc_start[,nuc_ff][length(ff_nuc_start[!is.na(ff_nuc_start[,nuc_ff]),][,nuc_ff])+1] = nuc_ff_count
#         }
#       }
#       
#       codon_string = ''
#       for (nuc_number in first_nuc_pos){
#         if(thousand[[seq_number]][nuc_number] != '-' & thousand[[seq_number]][nuc_number+1] != '-' & thousand[[seq_number]][nuc_number+2] != '-' ){
#           codon_string = append(codon_string, paste(thousand[[seq_number]][nuc_number], thousand[[seq_number]][nuc_number+1], thousand[[seq_number]][nuc_number+2], sep = ''), after = length(codon_string))
#         }
#       }
#       codon_string = codon_string[codon_string!='']
#       for (y in 1 : length(table(codon_string))){
#         codon = names(table(codon_string)[y])
#         codon_count = table(codon_string)[[y]]
#         codon_start[,codon][length(codon_start[!is.na(codon_start[,codon]),][,codon])+1] = codon_count
#     }
#     }else{
#       for (x in 1 : length(table(thousand[[seq_number]]))){
#         nucleotide = names(table(thousand[[seq_number]])[x])
#         nuc_count = table(thousand[[seq_number]])[[x]]
#         if (nucleotide!='-'){
#           nuc_end[,nucleotide][length(nuc_end[!is.na(nuc_end[,nucleotide]),][,nucleotide])+1] = nuc_count
#         }
#       }
#       syn_seq = ''
#       for (syn_nuc_number in syn_nuc_pos){
#         syn_seq = append(syn_seq, thousand[[seq_number]][syn_nuc_number], after = length(syn_seq))
#       }
#       syn_seq = syn_seq[syn_seq!='']
#       for (syn_nuc in 1 : length(table(syn_seq))){
#         if (names(table(syn_seq)[syn_nuc]) != '-'){
#           nuc_syn = names(table(syn_seq)[syn_nuc])
#           nuc_syn_count = table(syn_seq)[[syn_nuc]]
#           syn_nuc_end[,nuc_syn][length(syn_nuc_end[!is.na(syn_nuc_end[,nuc_syn]),][,nuc_syn])+1] = nuc_syn_count
#         }
#       }
#       four_fold_seq = ''
#       for (ff_nuc_number in four_fold_nuc_pos){
#         four_fold_seq = append(four_fold_seq, thousand[[seq_number]][ff_nuc_number], after = length(four_fold_seq))
#       }
#       four_fold_seq = four_fold_seq[four_fold_seq!='']
#       for (ff_nuc in 1 : length(table(four_fold_seq))){
#         if (names(table(four_fold_seq)[ff_nuc]) != '-'){
#           nuc_ff = names(table(four_fold_seq)[ff_nuc])
#           nuc_ff_count = table(four_fold_seq)[[ff_nuc]]
#           ff_nuc_end[,nuc_ff][length(ff_nuc_end[!is.na(ff_nuc_end[,nuc_ff]),][,nuc_ff])+1] = nuc_ff_count
#         }
#       }
#       codon_string = ''
#       for (nuc_number in first_nuc_pos){
#         if(thousand[[seq_number]][nuc_number] != '-' & thousand[[seq_number]][nuc_number+1] != '-' & thousand[[seq_number]][nuc_number+2] != '-' ){
#           codon_string = append(codon_string, paste(thousand[[seq_number]][nuc_number], thousand[[seq_number]][nuc_number+1], thousand[[seq_number]][nuc_number+2], sep = ''), after = length(codon_string))
#         }
#       }
#       codon_string = codon_string[codon_string!='']
#       for (y in 1 : length(table(codon_string))){
#         codon = names(table(codon_string)[y])
#         codon_count = table(codon_string)[[y]]
#         codon_end[,codon][length(codon_end[!is.na(codon_end[,codon]),][,codon])+1] = codon_count
#       }
#     }
#   }
# }

write.csv(nuc_start, 'D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/data_obtained/19.Nuc_start_1000.csv')
write.csv(nuc_end, 'D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/data_obtained/19.Nuc_end_1000.csv')
write.csv(codon_start, 'D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/data_obtained/19.Codons_start_1000.csv')
write.csv(codon_end, 'D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/data_obtained/19.Codons_end_1000.csv')
write.csv(syn_nuc_start, 'D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/data_obtained/19.Syn_nuc_start_1000.csv')
write.csv(syn_nuc_end, 'D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/data_obtained/19.Syn_nuc_end_1000.csv')
write.csv(ff_nuc_start, 'D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/data_obtained/19.FF_start_1000.csv')
write.csv(ff_nuc_end, 'D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/data_obtained/19.FF_end_1000.csv')

nuc_start = read.csv('D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/data_obtained/19.Nuc_start_1000.csv')
nuc_end = read.csv('D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/data_obtained/19.Nuc_end_1000.csv')
syn_nuc_start = read.csv('D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/data_obtained/19.Syn_nuc_start_1000.csv')
syn_nuc_end = read.csv('D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/data_obtained/19.Syn_nuc_end_1000.csv')
ff_nuc_start = read.csv('D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/data_obtained/19.FF_start_1000.csv')
ff_nuc_end = read.csv('D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/data_obtained/19.FF_end_1000.csv')
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
write.csv(nuc_usage, 'D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/data_obtained/19.Nuc_Usage.csv')
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



#pdf('D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/R_figure/boxplot/19.All.pdf')
boxplot(nuc_start$A_normed, nuc_end$A_normed, names = c('A_normed_start', 'A_normed_end'), main = 'Nucleotide "A" for full genome')
boxplot(nuc_start$T_normed, nuc_end$T_normed, names = c('T_normed_start', 'T_normed_end'), add = T)
boxplot(nuc_start$G_normed, nuc_end$G_normed, names = c('G_normed_start', 'G_normed_end'), add = T)
boxplot(nuc_start$C_normed, nuc_end$C_normed, names = c('C_normed_start', 'C_normed_end'), add = T)
boxplot(syn_nuc_start$A_normed, syn_nuc_end$A_normed, names = c('A_normed_start', 'A_normed_end'), add = T)
boxplot(syn_nuc_start$T_normed, syn_nuc_end$T_normed, names = c('T_normed_start', 'T_normed_end'), add = T)
boxplot(syn_nuc_start$G_normed, syn_nuc_end$G_normed, names = c('G_normed_start', 'G_normed_end'), add = T)
boxplot(syn_nuc_start$C_normed, syn_nuc_end$C_normed, names = c('C_normed_start', 'C_normed_end'), add = T)
boxplot(ff_nuc_start$A_normed, ff_nuc_end$A_normed, names = c('A_normed_start', 'A_normed_end'), add = T)
boxplot(ff_nuc_start$T_normed, ff_nuc_end$T_normed, names = c('T_normed_start', 'T_normed_end'), add = T)
boxplot(ff_nuc_start$G_normed, ff_nuc_end$G_normed, names = c('G_normed_start', 'G_normed_end'), add = T)
boxplot(ff_nuc_start$C_normed, ff_nuc_end$C_normed, names = c('C_normed_start', 'C_normed_end'), add = T)
#dev.off()
p2 <- ggplot(data, aes(x=variety, y=note, fill=treatment)) + 
  geom_boxplot() +
  facet_wrap(~variety, scale="free")

count = ''
nucleotide = ''
group = ''
for(i in 1: nrow(nuc_start)){
  count = append(count, nuc_start[i,]$A_normed, after = length(count))
  count = append(count, nuc_start[i,]$T_normed, after = length(count))
  count = append(count, nuc_start[i,]$G_normed, after = length(count))
  count = append(count, nuc_start[i,]$C_normed, after = length(count))
  
  nucleotide = append(nucleotide, 'A_all', after = length(nucleotide))
  nucleotide = append(nucleotide, 'T_all', after = length(nucleotide))
  nucleotide = append(nucleotide, 'G_all', after = length(nucleotide))
  nucleotide = append(nucleotide, 'C_all', after = length(nucleotide))
  
  group = append(group, 'start', after = length(group))
  group = append(group, 'start', after = length(group))
  group = append(group, 'start', after = length(group))
  group = append(group, 'start', after = length(group))
  
  count = append(count, nuc_end[i,]$A_normed, after = length(count))
  count = append(count, nuc_end[i,]$T_normed, after = length(count))
  count = append(count, nuc_end[i,]$G_normed, after = length(count))
  count = append(count, nuc_end[i,]$C_normed, after = length(count))
  
  nucleotide = append(nucleotide, 'A_all', after = length(nucleotide))
  nucleotide = append(nucleotide, 'T_all', after = length(nucleotide))
  nucleotide = append(nucleotide, 'G_all', after = length(nucleotide))
  nucleotide = append(nucleotide, 'C_all', after = length(nucleotide))
  
  group = append(group, 'end', after = length(group))
  group = append(group, 'end', after = length(group))
  group = append(group, 'end', after = length(group))
  group = append(group, 'end', after = length(group))
  
  count = append(count, syn_nuc_start[i,]$A_normed, after = length(count))
  count = append(count, syn_nuc_start[i,]$T_normed, after = length(count))
  count = append(count, syn_nuc_start[i,]$G_normed, after = length(count))
  count = append(count, syn_nuc_start[i,]$C_normed, after = length(count))
  
  nucleotide = append(nucleotide, 'A_3_pos', after = length(nucleotide))
  nucleotide = append(nucleotide, 'T_3_pos', after = length(nucleotide))
  nucleotide = append(nucleotide, 'G_3_pos', after = length(nucleotide))
  nucleotide = append(nucleotide, 'C_3_pos', after = length(nucleotide))
  
  group = append(group, 'start', after = length(group))
  group = append(group, 'start', after = length(group))
  group = append(group, 'start', after = length(group))
  group = append(group, 'start', after = length(group))
  
  count = append(count, syn_nuc_end[i,]$A_normed, after = length(count))
  count = append(count, syn_nuc_end[i,]$T_normed, after = length(count))
  count = append(count, syn_nuc_end[i,]$G_normed, after = length(count))
  count = append(count, syn_nuc_end[i,]$C_normed, after = length(count))
  
  nucleotide = append(nucleotide, 'A_3_pos', after = length(nucleotide))
  nucleotide = append(nucleotide, 'T_3_pos', after = length(nucleotide))
  nucleotide = append(nucleotide, 'G_3_pos', after = length(nucleotide))
  nucleotide = append(nucleotide, 'C_3_pos', after = length(nucleotide))
  
  group = append(group, 'end', after = length(group))
  group = append(group, 'end', after = length(group))
  group = append(group, 'end', after = length(group))
  group = append(group, 'end', after = length(group))
  
  count = append(count, ff_nuc_start[i,]$A_normed, after = length(count))
  count = append(count, ff_nuc_start[i,]$T_normed, after = length(count))
  count = append(count, ff_nuc_start[i,]$G_normed, after = length(count))
  count = append(count, ff_nuc_start[i,]$C_normed, after = length(count))
  
  nucleotide = append(nucleotide, 'A_ff', after = length(nucleotide))
  nucleotide = append(nucleotide, 'T_ff', after = length(nucleotide))
  nucleotide = append(nucleotide, 'G_ff', after = length(nucleotide))
  nucleotide = append(nucleotide, 'C_ff', after = length(nucleotide))
  
  group = append(group, 'start', after = length(group))
  group = append(group, 'start', after = length(group))
  group = append(group, 'start', after = length(group))
  group = append(group, 'start', after = length(group))
  
  count = append(count, ff_nuc_end[i,]$A_normed, after = length(count))
  count = append(count, ff_nuc_end[i,]$T_normed, after = length(count))
  count = append(count, ff_nuc_end[i,]$G_normed, after = length(count))
  count = append(count, ff_nuc_end[i,]$C_normed, after = length(count))
  
  nucleotide = append(nucleotide, 'A_ff', after = length(nucleotide))
  nucleotide = append(nucleotide, 'T_ff', after = length(nucleotide))
  nucleotide = append(nucleotide, 'G_ff', after = length(nucleotide))
  nucleotide = append(nucleotide, 'C_ff', after = length(nucleotide))
  
  group = append(group, 'end', after = length(group))
  group = append(group, 'end', after = length(group))
  group = append(group, 'end', after = length(group))
  group = append(group, 'end', after = length(group))
}

count = count[count!='']
nucleotide = nucleotide[nucleotide!='']
group = group[group!='']

nuc_usage_count = data.frame(count,nucleotide,group)
names(nuc_usage_count) = c("Nuc_Count","Nuc","Group")
write.csv(nuc_usage_count, 'D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/R_function/gisaid_mutspec/data/19.Nuc_Usage.csv')
scaleFUN <- function(x) sprintf("%.2f", x)
fig = ggplot(nuc_usage_count, aes(x = Nuc, y = Nuc_Count, fill = Group)) + geom_boxplot() + 
  facet_grid(~Nuc, scale="free") + scale_y_continuous(labels=scaleFUN)
fig %>% ggplotly
ggsave('D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/R_function/gisaid_mutspec/figures/19.Nuc_Usage.svg', plot = fig, device = 'svg',  limitsize = F)