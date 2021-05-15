library(plotly)
library(data.table)
library("heatmaply")
library(seqinr)

thousand = read.fasta('D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/data/first_last_1000_mulal.fasta',forceDNAtolower = FALSE)
ideal_table = read.csv("D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/data_obtained/ideal_table.csv")

orf1_pos = ideal_table[ideal_table$GenName == 'ORF1ab' & ideal_table$NucInCodon==1,]$Pos
orf1_pos = unique(orf1_pos)
no_ofr1_pos = ideal_table[ideal_table$GenName != 'ORF1ab' & ideal_table$GenType == 'translated' & ideal_table$NucInCodon==1,]$Pos
no_ofr1_pos = unique(no_ofr1_pos)

codon_orf1_start = data.frame(matrix(vector(), 1000, length(unique(ideal_table[ideal_table$GenType == 'translated',]$AltCodon)),
                                dimnames=list(c(), unique(ideal_table[ideal_table$GenType == 'translated',]$AltCodon))),
                         stringsAsFactors=F)
codon_orf1_end = data.frame(matrix(vector(), 1000, length(unique(ideal_table[ideal_table$GenType == 'translated',]$AltCodon)),
                              dimnames=list(c(), unique(ideal_table[ideal_table$GenType == 'translated',]$AltCodon))),
                       stringsAsFactors=F)
codon_no_orf1_start = data.frame(matrix(vector(), 1000, length(unique(ideal_table[ideal_table$GenType == 'translated',]$AltCodon)),
                                     dimnames=list(c(), unique(ideal_table[ideal_table$GenType == 'translated',]$AltCodon))),
                              stringsAsFactors=F)
codon_no_orf1_end = data.frame(matrix(vector(), 1000, length(unique(ideal_table[ideal_table$GenType == 'translated',]$AltCodon)),
                                   dimnames=list(c(), unique(ideal_table[ideal_table$GenType == 'translated',]$AltCodon))),
                            stringsAsFactors=F)

for (seq_number in (1 : length(names(thousand)))){
  print(seq_number)
  if (seq_number != 1){
    if (seq_number<1002){
      codon_string = ''
      for (orf_pos in orf1_pos){
        if(thousand[[seq_number]][orf_pos] != '-' & thousand[[seq_number]][orf_pos+1] != '-' & thousand[[seq_number]][orf_pos+2] != '-' ){
          codon_string = append(codon_string, paste(thousand[[seq_number]][orf_pos], thousand[[seq_number]][orf_pos+1], thousand[[seq_number]][orf_pos+2], sep = ''), after = length(codon_string))
        }
      }
      codon_string = codon_string[codon_string!='']
      for (y in 1 : length(table(codon_string))){
        codon = names(table(codon_string)[y])
        codon_count = table(codon_string)[[y]]
        codon_orf1_start[,codon][length(codon_orf1_start[!is.na(codon_orf1_start[,codon]),][,codon])+1] = codon_count
      }
      codon_string = ''
      for (no_orf_pos in no_ofr1_pos){
        if(thousand[[seq_number]][no_orf_pos] != '-' & thousand[[seq_number]][no_orf_pos+1] != '-' & thousand[[seq_number]][no_orf_pos+2] != '-' ){
          codon_string = append(codon_string, paste(thousand[[seq_number]][no_orf_pos], thousand[[seq_number]][no_orf_pos+1], thousand[[seq_number]][no_orf_pos+2], sep = ''), after = length(codon_string))
        }
      }
      codon_string = codon_string[codon_string!='']
      for (y in 1 : length(table(codon_string))){
        codon = names(table(codon_string)[y])
        codon_count = table(codon_string)[[y]]
        codon_no_orf1_start[,codon][length(codon_no_orf1_start[!is.na(codon_no_orf1_start[,codon]),][,codon])+1] = codon_count
      }
    }else{
      codon_string = ''
      for (orf_pos in orf1_pos){
        if(thousand[[seq_number]][orf_pos] != '-' & thousand[[seq_number]][orf_pos+1] != '-' & thousand[[seq_number]][orf_pos+2] != '-' ){
          codon_string = append(codon_string, paste(thousand[[seq_number]][orf_pos], thousand[[seq_number]][orf_pos+1], thousand[[seq_number]][orf_pos+2], sep = ''), after = length(codon_string))
        }
      }
      codon_string = codon_string[codon_string!='']
      for (y in 1 : length(table(codon_string))){
        codon = names(table(codon_string)[y])
        codon_count = table(codon_string)[[y]]
        codon_orf1_end[,codon][length(codon_orf1_end[!is.na(codon_orf1_end[,codon]),][,codon])+1] = codon_count
      }
      codon_string = ''
      for (no_orf_pos in no_ofr1_pos){
        if(thousand[[seq_number]][no_orf_pos] != '-' & thousand[[seq_number]][no_orf_pos+1] != '-' & thousand[[seq_number]][no_orf_pos+2] != '-' ){
          codon_string = append(codon_string, paste(thousand[[seq_number]][no_orf_pos], thousand[[seq_number]][no_orf_pos+1], thousand[[seq_number]][no_orf_pos+2], sep = ''), after = length(codon_string))
        }
      }
      codon_string = codon_string[codon_string!='']
      for (y in 1 : length(table(codon_string))){
        codon = names(table(codon_string)[y])
        codon_count = table(codon_string)[[y]]
        codon_no_orf1_end[,codon][length(codon_no_orf1_end[!is.na(codon_no_orf1_end[,codon]),][,codon])+1] = codon_count
      }
    }
  }
}
write.csv(codon_orf1_start, 'D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/data_obtained/20.Codon_orf1ab_start_1000.csv')
write.csv(codon_no_orf1_start, 'D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/data_obtained/20.Codon_no_orf1ab_start_1000.csv')
write.csv(codon_orf1_end, 'D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/data_obtained/20.Codon_orf1ab_end_1000.csv')
write.csv(codon_no_orf1_end, 'D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/data_obtained/20.Codon_no_orf1ab_end_1000.csv')

#orf start
codon_orf1_start[is.na(codon_orf1_start)] = 0
codons = names(codon_orf1_start)
first_thousand_orf1ab = data.frame(codons)
first_thousand_orf1ab$codon_summ = 0
for (codon in first_thousand_orf1ab$codons){
  first_thousand_orf1ab[first_thousand_orf1ab$codons == codon,]$codon_summ = sum(codon_orf1_start[,codon])
}
first_thousand_orf1ab$Aa = ''
for (i in first_thousand_orf1ab$codons){
  first_thousand_orf1ab[first_thousand_orf1ab$codons == i,]$Aa = translate(s2c(i))
}
first_thousand_orf1ab[first_thousand_orf1ab$codons == 'CTT' | first_thousand_orf1ab$codons == 'CTA' | first_thousand_orf1ab$codons == 'CTG' | first_thousand_orf1ab$codons == 'CTC',]$Aa = 'L_CT'
first_thousand_orf1ab[first_thousand_orf1ab$codons == 'TTA' | first_thousand_orf1ab$codons == 'TTG',]$Aa = 'L_TT'

first_thousand_orf1ab[first_thousand_orf1ab$codons == 'TCT' | first_thousand_orf1ab$codons == 'TCA' | first_thousand_orf1ab$codons == 'TCG' | first_thousand_orf1ab$codons == 'TCC',]$Aa = 'S_TC'
first_thousand_orf1ab[first_thousand_orf1ab$codons == 'AGT' | first_thousand_orf1ab$codons == 'AGC',]$Aa = 'S_AG'

first_thousand_orf1ab[first_thousand_orf1ab$codons == 'CGC' | first_thousand_orf1ab$codons == 'CGA' | first_thousand_orf1ab$codons == 'CGT' | first_thousand_orf1ab$codons == 'CGG',]$Aa = 'R_CG'
first_thousand_orf1ab[first_thousand_orf1ab$codons == 'AGG' | first_thousand_orf1ab$codons == 'AGA',]$Aa = 'R_AG'

first_thousand_orf1ab[first_thousand_orf1ab$codons == 'TAG' | first_thousand_orf1ab$codons == 'TAA',]$Aa = '*_TA'
first_thousand_orf1ab[first_thousand_orf1ab$codons == 'TGA',]$Aa = '*_TG'
first_thousand_orf1ab$codon_one_seq = first_thousand_orf1ab$codon_summ / 1000
first_thousand_orf1ab$normed_codon = first_thousand_orf1ab$codon_one_seq / sum(first_thousand_orf1ab$codon_one_seq)

first_thousand_orf1ab$first = substring(first_thousand_orf1ab$codons, 1, 1)
first_thousand_orf1ab$second = substring(first_thousand_orf1ab$codons, 2, 2)
first_thousand_orf1ab$three = substring(first_thousand_orf1ab$codons, 3, 3)
first_thousand_orf1ab$first = factor(first_thousand_orf1ab$first, levels = c("C","T","G","A"))
first_thousand_orf1ab$second = factor(first_thousand_orf1ab$second, levels = c("C","T","G","A"))
first_thousand_orf1ab$three = factor(first_thousand_orf1ab$three, levels = c("C","T","G","A"))

df_first_thousand_orf1ab <-first_thousand_orf1ab[order(first_thousand_orf1ab$second, first_thousand_orf1ab$first, first_thousand_orf1ab$three, first_thousand_orf1ab$Aa),]
rownames(df_first_thousand_orf1ab) = 1: nrow(df_first_thousand_orf1ab)

df_first_thousand_orf1ab$normed_codons_in_aa = 0
df_first_thousand_orf1ab$normed_aa = 0
for (i in unique(df_first_thousand_orf1ab$Aa)){
  df_first_thousand_orf1ab[df_first_thousand_orf1ab$Aa == i,]$normed_codons_in_aa = df_first_thousand_orf1ab[df_first_thousand_orf1ab$Aa == i,]$codon_one_seq / sum(df_first_thousand_orf1ab[df_first_thousand_orf1ab$Aa == i,]$codon_one_seq)
  df_first_thousand_orf1ab[df_first_thousand_orf1ab$Aa == i,]$normed_aa = sum(df_first_thousand_orf1ab[df_first_thousand_orf1ab$Aa == i,]$codon_one_seq) / sum(df_first_thousand_orf1ab$codon_one_seq)
}
#no orf start
codon_no_orf1_start[is.na(codon_no_orf1_start)] = 0
codons = names(codon_no_orf1_start)
first_thousand_no_orf1ab = data.frame(codons)
first_thousand_no_orf1ab$codon_summ = 0
for (codon in first_thousand_orf1ab$codons){
  first_thousand_no_orf1ab[first_thousand_no_orf1ab$codons == codon,]$codon_summ = sum(codon_no_orf1_start[,codon])
}
first_thousand_no_orf1ab$Aa = ''
for (i in first_thousand_no_orf1ab$codons){
  first_thousand_no_orf1ab[first_thousand_no_orf1ab$codons == i,]$Aa = translate(s2c(i))
}
first_thousand_no_orf1ab[first_thousand_no_orf1ab$codons == 'CTT' | first_thousand_no_orf1ab$codons == 'CTA' | first_thousand_no_orf1ab$codons == 'CTG' | first_thousand_no_orf1ab$codons == 'CTC',]$Aa = 'L_CT'
first_thousand_no_orf1ab[first_thousand_no_orf1ab$codons == 'TTA' | first_thousand_no_orf1ab$codons == 'TTG',]$Aa = 'L_TT'

first_thousand_no_orf1ab[first_thousand_no_orf1ab$codons == 'TCT' | first_thousand_no_orf1ab$codons == 'TCA' | first_thousand_no_orf1ab$codons == 'TCG' | first_thousand_no_orf1ab$codons == 'TCC',]$Aa = 'S_TC'
first_thousand_no_orf1ab[first_thousand_no_orf1ab$codons == 'AGT' | first_thousand_no_orf1ab$codons == 'AGC',]$Aa = 'S_AG'

first_thousand_no_orf1ab[first_thousand_no_orf1ab$codons == 'CGC' | first_thousand_no_orf1ab$codons == 'CGA' | first_thousand_no_orf1ab$codons == 'CGT' | first_thousand_no_orf1ab$codons == 'CGG',]$Aa = 'R_CG'
first_thousand_no_orf1ab[first_thousand_no_orf1ab$codons == 'AGG' | first_thousand_no_orf1ab$codons == 'AGA',]$Aa = 'R_AG'

first_thousand_no_orf1ab[first_thousand_no_orf1ab$codons == 'TAG' | first_thousand_no_orf1ab$codons == 'TAA',]$Aa = '*_TA'
first_thousand_no_orf1ab[first_thousand_no_orf1ab$codons == 'TGA',]$Aa = '*_TG'
first_thousand_no_orf1ab$codon_one_seq = first_thousand_no_orf1ab$codon_summ / 1000
first_thousand_no_orf1ab$normed_codon = first_thousand_no_orf1ab$codon_one_seq / sum(first_thousand_no_orf1ab$codon_one_seq)

first_thousand_no_orf1ab$first = substring(first_thousand_no_orf1ab$codons, 1, 1)
first_thousand_no_orf1ab$second = substring(first_thousand_no_orf1ab$codons, 2, 2)
first_thousand_no_orf1ab$three = substring(first_thousand_no_orf1ab$codons, 3, 3)
first_thousand_no_orf1ab$first = factor(first_thousand_no_orf1ab$first, levels = c("C","T","G","A"))
first_thousand_no_orf1ab$second = factor(first_thousand_no_orf1ab$second, levels = c("C","T","G","A"))
first_thousand_no_orf1ab$three = factor(first_thousand_no_orf1ab$three, levels = c("C","T","G","A"))

df_first_thousand_no_orf1ab <-first_thousand_no_orf1ab[order(first_thousand_no_orf1ab$second, first_thousand_no_orf1ab$first, first_thousand_no_orf1ab$three, first_thousand_no_orf1ab$Aa),]
rownames(df_first_thousand_no_orf1ab) = 1: nrow(df_first_thousand_no_orf1ab)

df_first_thousand_no_orf1ab$normed_codons_in_aa = 0
df_first_thousand_no_orf1ab$normed_aa = 0
for (i in unique(df_first_thousand_no_orf1ab$Aa)){
  df_first_thousand_no_orf1ab[df_first_thousand_no_orf1ab$Aa == i,]$normed_codons_in_aa = df_first_thousand_no_orf1ab[df_first_thousand_no_orf1ab$Aa == i,]$codon_one_seq / sum(df_first_thousand_no_orf1ab[df_first_thousand_no_orf1ab$Aa == i,]$codon_one_seq)
  df_first_thousand_no_orf1ab[df_first_thousand_no_orf1ab$Aa == i,]$normed_aa = sum(df_first_thousand_no_orf1ab[df_first_thousand_no_orf1ab$Aa == i,]$codon_one_seq) / sum(df_first_thousand_no_orf1ab$codon_one_seq)
}
#orf1ab end
codon_orf1_end[is.na(codon_orf1_end)] = 0
codons = names(codon_orf1_end)
last_thousand_orf1ab = data.frame(codons)
last_thousand_orf1ab$codon_summ = 0
for (codon in last_thousand_orf1ab$codons){
  last_thousand_orf1ab[last_thousand_orf1ab$codons == codon,]$codon_summ = sum(codon_orf1_end[,codon])
}
last_thousand_orf1ab$Aa = ''
for (i in last_thousand_orf1ab$codons){
  last_thousand_orf1ab[last_thousand_orf1ab$codons == i,]$Aa = translate(s2c(i))
}
last_thousand_orf1ab[last_thousand_orf1ab$codons == 'CTT' | last_thousand_orf1ab$codons == 'CTA' | last_thousand_orf1ab$codons == 'CTG' | last_thousand_orf1ab$codons == 'CTC',]$Aa = 'L_CT'
last_thousand_orf1ab[last_thousand_orf1ab$codons == 'TTA' | last_thousand_orf1ab$codons == 'TTG',]$Aa = 'L_TT'

last_thousand_orf1ab[last_thousand_orf1ab$codons == 'TCT' | last_thousand_orf1ab$codons == 'TCA' | last_thousand_orf1ab$codons == 'TCG' | last_thousand_orf1ab$codons == 'TCC',]$Aa = 'S_TC'
last_thousand_orf1ab[last_thousand_orf1ab$codons == 'AGT' | last_thousand_orf1ab$codons == 'AGC',]$Aa = 'S_AG'

last_thousand_orf1ab[last_thousand_orf1ab$codons == 'CGC' | last_thousand_orf1ab$codons == 'CGA' | last_thousand_orf1ab$codons == 'CGT' | last_thousand_orf1ab$codons == 'CGG',]$Aa = 'R_CG'
last_thousand_orf1ab[last_thousand_orf1ab$codons == 'AGG' | last_thousand_orf1ab$codons == 'AGA',]$Aa = 'R_AG'

last_thousand_orf1ab[last_thousand_orf1ab$codons == 'TAG' | last_thousand_orf1ab$codons == 'TAA',]$Aa = '*_TA'
last_thousand_orf1ab[last_thousand_orf1ab$codons == 'TGA',]$Aa = '*_TG'
last_thousand_orf1ab$codon_one_seq = last_thousand_orf1ab$codon_summ / 1000
last_thousand_orf1ab$normed_codon = last_thousand_orf1ab$codon_one_seq / sum(last_thousand_orf1ab$codon_one_seq)

last_thousand_orf1ab$first = substring(last_thousand_orf1ab$codons, 1, 1)
last_thousand_orf1ab$second = substring(last_thousand_orf1ab$codons, 2, 2)
last_thousand_orf1ab$three = substring(last_thousand_orf1ab$codons, 3, 3)
last_thousand_orf1ab$first = factor(last_thousand_orf1ab$first, levels = c("C","T","G","A"))
last_thousand_orf1ab$second = factor(last_thousand_orf1ab$second, levels = c("C","T","G","A"))
last_thousand_orf1ab$three = factor(last_thousand_orf1ab$three, levels = c("C","T","G","A"))

df_last_thousand_orf1ab <-last_thousand_orf1ab[order(last_thousand_orf1ab$second, last_thousand_orf1ab$first, last_thousand_orf1ab$three, last_thousand_orf1ab$Aa),]
rownames(df_last_thousand_orf1ab) = 1: nrow(df_last_thousand_orf1ab)

df_last_thousand_orf1ab$normed_codons_in_aa = 0
df_last_thousand_orf1ab$normed_aa = 0
for (i in unique(df_last_thousand_orf1ab$Aa)){
  df_last_thousand_orf1ab[df_last_thousand_orf1ab$Aa == i,]$normed_codons_in_aa = df_last_thousand_orf1ab[df_last_thousand_orf1ab$Aa == i,]$codon_one_seq / sum(df_last_thousand_orf1ab[df_last_thousand_orf1ab$Aa == i,]$codon_one_seq)
  df_last_thousand_orf1ab[df_last_thousand_orf1ab$Aa == i,]$normed_aa = sum(df_last_thousand_orf1ab[df_last_thousand_orf1ab$Aa == i,]$codon_one_seq) / sum(df_last_thousand_orf1ab$codon_one_seq)
}
#no orf1ab end
codon_no_orf1_end[is.na(codon_no_orf1_end)] = 0
codons = names(codon_no_orf1_end)
last_thousand_no_orf1ab = data.frame(codons)
last_thousand_no_orf1ab$codon_summ = 0
for (codon in first_thousand_orf1ab$codons){
  last_thousand_no_orf1ab[last_thousand_no_orf1ab$codons == codon,]$codon_summ = sum(codon_no_orf1_end[,codon])
}
last_thousand_no_orf1ab$Aa = ''
for (i in last_thousand_no_orf1ab$codons){
  last_thousand_no_orf1ab[last_thousand_no_orf1ab$codons == i,]$Aa = translate(s2c(i))
}
last_thousand_no_orf1ab[last_thousand_no_orf1ab$codons == 'CTT' | last_thousand_no_orf1ab$codons == 'CTA' | last_thousand_no_orf1ab$codons == 'CTG' | last_thousand_no_orf1ab$codons == 'CTC',]$Aa = 'L_CT'
last_thousand_no_orf1ab[last_thousand_no_orf1ab$codons == 'TTA' | last_thousand_no_orf1ab$codons == 'TTG',]$Aa = 'L_TT'

last_thousand_no_orf1ab[last_thousand_no_orf1ab$codons == 'TCT' | last_thousand_no_orf1ab$codons == 'TCA' | last_thousand_no_orf1ab$codons == 'TCG' | last_thousand_no_orf1ab$codons == 'TCC',]$Aa = 'S_TC'
last_thousand_no_orf1ab[last_thousand_no_orf1ab$codons == 'AGT' | last_thousand_no_orf1ab$codons == 'AGC',]$Aa = 'S_AG'

last_thousand_no_orf1ab[last_thousand_no_orf1ab$codons == 'CGC' | last_thousand_no_orf1ab$codons == 'CGA' | last_thousand_no_orf1ab$codons == 'CGT' | last_thousand_no_orf1ab$codons == 'CGG',]$Aa = 'R_CG'
last_thousand_no_orf1ab[last_thousand_no_orf1ab$codons == 'AGG' | last_thousand_no_orf1ab$codons == 'AGA',]$Aa = 'R_AG'

last_thousand_no_orf1ab[last_thousand_no_orf1ab$codons == 'TAG' | last_thousand_no_orf1ab$codons == 'TAA',]$Aa = '*_TA'
last_thousand_no_orf1ab[last_thousand_no_orf1ab$codons == 'TGA',]$Aa = '*_TG'
last_thousand_no_orf1ab$codon_one_seq = last_thousand_no_orf1ab$codon_summ / 1000
last_thousand_no_orf1ab$normed_codon = last_thousand_no_orf1ab$codon_one_seq / sum(last_thousand_no_orf1ab$codon_one_seq)

last_thousand_no_orf1ab$first = substring(last_thousand_no_orf1ab$codons, 1, 1)
last_thousand_no_orf1ab$second = substring(last_thousand_no_orf1ab$codons, 2, 2)
last_thousand_no_orf1ab$three = substring(last_thousand_no_orf1ab$codons, 3, 3)
last_thousand_no_orf1ab$first = factor(last_thousand_no_orf1ab$first, levels = c("C","T","G","A"))
last_thousand_no_orf1ab$second = factor(last_thousand_no_orf1ab$second, levels = c("C","T","G","A"))
last_thousand_no_orf1ab$three = factor(last_thousand_no_orf1ab$three, levels = c("C","T","G","A"))

df_last_thousand_no_orf1ab <-last_thousand_no_orf1ab[order(last_thousand_no_orf1ab$second, last_thousand_no_orf1ab$first, last_thousand_no_orf1ab$three, last_thousand_no_orf1ab$Aa),]
rownames(df_last_thousand_no_orf1ab) = 1: nrow(df_last_thousand_no_orf1ab)

df_last_thousand_no_orf1ab$normed_codons_in_aa = 0
df_last_thousand_no_orf1ab$normed_aa = 0
for (i in unique(df_last_thousand_no_orf1ab$Aa)){
  df_last_thousand_no_orf1ab[df_last_thousand_no_orf1ab$Aa == i,]$normed_codons_in_aa = df_last_thousand_no_orf1ab[df_last_thousand_no_orf1ab$Aa == i,]$codon_one_seq / sum(df_last_thousand_no_orf1ab[df_last_thousand_no_orf1ab$Aa == i,]$codon_one_seq)
  df_last_thousand_no_orf1ab[df_last_thousand_no_orf1ab$Aa == i,]$normed_aa = sum(df_last_thousand_no_orf1ab[df_last_thousand_no_orf1ab$Aa == i,]$codon_one_seq) / sum(df_last_thousand_no_orf1ab$codon_one_seq)
}

#orf1ab codons 1000 start
A_Nuc = df_first_thousand_orf1ab[df_first_thousand_orf1ab$second=="A",]$normed_codon
C_Nuc = df_first_thousand_orf1ab[df_first_thousand_orf1ab$second=="C",]$normed_codon
G_Nuc = df_first_thousand_orf1ab[df_first_thousand_orf1ab$second=="G",]$normed_codon
T_Nuc = df_first_thousand_orf1ab[df_first_thousand_orf1ab$second=="T",]$normed_codon

heat_df = data.frame("C" = C_Nuc,"T" = T_Nuc,"G" = G_Nuc,"A" = A_Nuc)
row.names(heat_df) = c("C1","C2","C3","C4","T1","T2","T3","T4","G1","G2","G3","G4","A1","A2","A3","A4")
df_first_thousand_orf1ab$codon_aa_count = paste(df_first_thousand_orf1ab$codons, "(",df_first_thousand_orf1ab$Aa,")", round(df_first_thousand_orf1ab$normed_codon,4), sep = " ")

A_nuc_text = df_first_thousand_orf1ab[df_first_thousand_orf1ab$second=="A",]$codon_aa_count
C_nuc_text = df_first_thousand_orf1ab[df_first_thousand_orf1ab$second=="C",]$codon_aa_count
G_nuc_text = df_first_thousand_orf1ab[df_first_thousand_orf1ab$second=="G",]$codon_aa_count
T_nuc_text = df_first_thousand_orf1ab[df_first_thousand_orf1ab$second=="T",]$codon_aa_count
text_df = data.frame("C" = C_nuc_text,"T" = T_nuc_text,"G" = G_nuc_text,"A" = A_nuc_text)
heatmaply(
  as.matrix(heat_df),
  colors = viridis(n = 256,  option = "magma"),
  fontsize_row = 15,
  fontsize_col = 15,
  Rowv = NULL,
  Colv = NULL,
  xlab = "Second Nucleotide",
  ylab = 'First Nucleotide',
  main = 'Codon Usage for first 1000 sequences in ORF1ab',
  show_dendrogram = c(FALSE, FALSE),
  plot_method = "plotly",
  draw_cellnote = T,
  cellnote = text_df,
  cellnote_size = 13 ,
  cellnote_textposition = "middle center",
  file = "D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/R_figure/HeatMap/orf_no_orf/20.Codon_Usage_orf1ab_start1000.html"
)

#no orf1ab codons 1000 start
A_Nuc = df_first_thousand_no_orf1ab[df_first_thousand_no_orf1ab$second=="A",]$normed_codon
C_Nuc = df_first_thousand_no_orf1ab[df_first_thousand_no_orf1ab$second=="C",]$normed_codon
G_Nuc = df_first_thousand_no_orf1ab[df_first_thousand_no_orf1ab$second=="G",]$normed_codon
T_Nuc = df_first_thousand_no_orf1ab[df_first_thousand_no_orf1ab$second=="T",]$normed_codon

heat_df = data.frame("C" = C_Nuc,"T" = T_Nuc,"G" = G_Nuc,"A" = A_Nuc)
row.names(heat_df) = c("C1","C2","C3","C4","T1","T2","T3","T4","G1","G2","G3","G4","A1","A2","A3","A4")
df_first_thousand_no_orf1ab$codon_aa_count = paste(df_first_thousand_no_orf1ab$codons, "(",df_first_thousand_no_orf1ab$Aa,")", round(df_first_thousand_no_orf1ab$normed_codon,4), sep = " ")

A_nuc_text = df_first_thousand_no_orf1ab[df_first_thousand_no_orf1ab$second=="A",]$codon_aa_count
C_nuc_text = df_first_thousand_no_orf1ab[df_first_thousand_no_orf1ab$second=="C",]$codon_aa_count
G_nuc_text = df_first_thousand_no_orf1ab[df_first_thousand_no_orf1ab$second=="G",]$codon_aa_count
T_nuc_text = df_first_thousand_no_orf1ab[df_first_thousand_no_orf1ab$second=="T",]$codon_aa_count
text_df = data.frame("C" = C_nuc_text,"T" = T_nuc_text,"G" = G_nuc_text,"A" = A_nuc_text)
heatmaply(
  as.matrix(heat_df),
  colors = viridis(n = 256,  option = "magma"),
  fontsize_row = 15,
  fontsize_col = 15,
  Rowv = NULL,
  Colv = NULL,
  xlab = "Second Nucleotide",
  ylab = 'First Nucleotide',
  main = 'Codon Usage for first 1000 sequences beside orf1ab',
  show_dendrogram = c(FALSE, FALSE),
  plot_method = "plotly",
  draw_cellnote = T,
  cellnote = text_df,
  cellnote_size = 13 ,
  cellnote_textposition = "middle center",
  file = "D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/R_figure/HeatMap/orf_no_orf/20.Codon_Usage_no_orf1ab_start1000.html"
)
#orf1ab codons in aa 1000 start
A_Nuc = df_first_thousand_orf1ab[df_first_thousand_orf1ab$second=="A",]$normed_codons_in_aa
C_Nuc = df_first_thousand_orf1ab[df_first_thousand_orf1ab$second=="C",]$normed_codons_in_aa
G_Nuc = df_first_thousand_orf1ab[df_first_thousand_orf1ab$second=="G",]$normed_codons_in_aa
T_Nuc = df_first_thousand_orf1ab[df_first_thousand_orf1ab$second=="T",]$normed_codons_in_aa

heat_df = data.frame("C" = C_Nuc,"T" = T_Nuc,"G" = G_Nuc,"A" = A_Nuc)
row.names(heat_df) = c("C1","C2","C3","C4","T1","T2","T3","T4","G1","G2","G3","G4","A1","A2","A3","A4")
df_first_thousand_orf1ab$codon_aa_count = paste(df_first_thousand_orf1ab$codons, "(",df_first_thousand_orf1ab$Aa,")", round(df_first_thousand_orf1ab$normed_codons_in_aa,4), sep = " ")

A_nuc_text = df_first_thousand_orf1ab[df_first_thousand_orf1ab$second=="A",]$codon_aa_count
C_nuc_text = df_first_thousand_orf1ab[df_first_thousand_orf1ab$second=="C",]$codon_aa_count
G_nuc_text = df_first_thousand_orf1ab[df_first_thousand_orf1ab$second=="G",]$codon_aa_count
T_nuc_text = df_first_thousand_orf1ab[df_first_thousand_orf1ab$second=="T",]$codon_aa_count
text_df = data.frame("C" = C_nuc_text,"T" = T_nuc_text,"G" = G_nuc_text,"A" = A_nuc_text)
heatmaply(
  as.matrix(heat_df),
  colors = viridis(n = 256,  option = "magma"),
  fontsize_row = 15,
  fontsize_col = 15,
  Rowv = NULL,
  Colv = NULL,
  xlab = "Second Nucleotide",
  ylab = 'First Nucleotide',
  main = 'Codon Usage in Aa for first 1000 sequences in ORF1ab',
  show_dendrogram = c(FALSE, FALSE),
  plot_method = "plotly",
  draw_cellnote = T,
  cellnote = text_df,
  cellnote_size = 13 ,
  cellnote_textposition = "middle center",
  file = "D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/R_figure/HeatMap/orf_no_orf/20.Codon_Usage_in_Aa_orf1ab_start1000.html"
)
#no orf1ab codons in aa 1000 start
A_Nuc = df_first_thousand_no_orf1ab[df_first_thousand_no_orf1ab$second=="A",]$normed_codons_in_aa
C_Nuc = df_first_thousand_no_orf1ab[df_first_thousand_no_orf1ab$second=="C",]$normed_codons_in_aa
G_Nuc = df_first_thousand_no_orf1ab[df_first_thousand_no_orf1ab$second=="G",]$normed_codons_in_aa
T_Nuc = df_first_thousand_no_orf1ab[df_first_thousand_no_orf1ab$second=="T",]$normed_codons_in_aa

heat_df = data.frame("C" = C_Nuc,"T" = T_Nuc,"G" = G_Nuc,"A" = A_Nuc)
row.names(heat_df) = c("C1","C2","C3","C4","T1","T2","T3","T4","G1","G2","G3","G4","A1","A2","A3","A4")
df_first_thousand_no_orf1ab$codon_aa_count = paste(df_first_thousand_no_orf1ab$codons, "(",df_first_thousand_no_orf1ab$Aa,")", round(df_first_thousand_no_orf1ab$normed_codons_in_aa,4), sep = " ")

A_nuc_text = df_first_thousand_no_orf1ab[df_first_thousand_no_orf1ab$second=="A",]$codon_aa_count
C_nuc_text = df_first_thousand_no_orf1ab[df_first_thousand_no_orf1ab$second=="C",]$codon_aa_count
G_nuc_text = df_first_thousand_no_orf1ab[df_first_thousand_no_orf1ab$second=="G",]$codon_aa_count
T_nuc_text = df_first_thousand_no_orf1ab[df_first_thousand_no_orf1ab$second=="T",]$codon_aa_count
text_df = data.frame("C" = C_nuc_text,"T" = T_nuc_text,"G" = G_nuc_text,"A" = A_nuc_text)
heatmaply(
  as.matrix(heat_df),
  colors = viridis(n = 256,  option = "magma"),
  fontsize_row = 15,
  fontsize_col = 15,
  Rowv = NULL,
  Colv = NULL,
  xlab = "Second Nucleotide",
  ylab = 'First Nucleotide',
  main = 'Codon Usage in Aa for first 1000 sequences beside orf1ab',
  show_dendrogram = c(FALSE, FALSE),
  plot_method = "plotly",
  draw_cellnote = T,
  cellnote = text_df,
  cellnote_size = 13 ,
  cellnote_textposition = "middle center",
  file = "D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/R_figure/HeatMap/orf_no_orf/20.Codon_Usage_in_Aa_no_orf1ab_start1000.html"
)
#orf1ab aa 1000 start
A_Nuc = df_first_thousand_orf1ab[df_first_thousand_orf1ab$second=="A",]$normed_aa
C_Nuc = df_first_thousand_orf1ab[df_first_thousand_orf1ab$second=="C",]$normed_aa
G_Nuc = df_first_thousand_orf1ab[df_first_thousand_orf1ab$second=="G",]$normed_aa
T_Nuc = df_first_thousand_orf1ab[df_first_thousand_orf1ab$second=="T",]$normed_aa

heat_df = data.frame("C" = C_Nuc,"T" = T_Nuc,"G" = G_Nuc,"A" = A_Nuc)
row.names(heat_df) = c("C1","C2","C3","C4","T1","T2","T3","T4","G1","G2","G3","G4","A1","A2","A3","A4")
df_first_thousand_orf1ab$codon_aa_count = paste(df_first_thousand_orf1ab$codons, "(",df_first_thousand_orf1ab$Aa,")", round(df_first_thousand_orf1ab$normed_aa,4), sep = " ")

A_nuc_text = df_first_thousand_orf1ab[df_first_thousand_orf1ab$second=="A",]$codon_aa_count
C_nuc_text = df_first_thousand_orf1ab[df_first_thousand_orf1ab$second=="C",]$codon_aa_count
G_nuc_text = df_first_thousand_orf1ab[df_first_thousand_orf1ab$second=="G",]$codon_aa_count
T_nuc_text = df_first_thousand_orf1ab[df_first_thousand_orf1ab$second=="T",]$codon_aa_count
text_df = data.frame("C" = C_nuc_text,"T" = T_nuc_text,"G" = G_nuc_text,"A" = A_nuc_text)
heatmaply(
  as.matrix(heat_df),
  colors = viridis(n = 256,  option = "magma"),
  fontsize_row = 15,
  fontsize_col = 15,
  Rowv = NULL,
  Colv = NULL,
  xlab = "Second Nucleotide",
  ylab = 'First Nucleotide',
  main = 'Aa Usage for first 1000 sequences in ORF1ab',
  show_dendrogram = c(FALSE, FALSE),
  plot_method = "plotly",
  draw_cellnote = T,
  cellnote = text_df,
  cellnote_size = 13 ,
  cellnote_textposition = "middle center",
  file = "D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/R_figure/HeatMap/orf_no_orf/20.Aa_Usage_orf1ab_start1000.html"
)
#no orf1ab aa 1000 start
A_Nuc = df_first_thousand_no_orf1ab[df_first_thousand_no_orf1ab$second=="A",]$normed_aa
C_Nuc = df_first_thousand_no_orf1ab[df_first_thousand_no_orf1ab$second=="C",]$normed_aa
G_Nuc = df_first_thousand_no_orf1ab[df_first_thousand_no_orf1ab$second=="G",]$normed_aa
T_Nuc = df_first_thousand_no_orf1ab[df_first_thousand_no_orf1ab$second=="T",]$normed_aa

heat_df = data.frame("C" = C_Nuc,"T" = T_Nuc,"G" = G_Nuc,"A" = A_Nuc)
row.names(heat_df) = c("C1","C2","C3","C4","T1","T2","T3","T4","G1","G2","G3","G4","A1","A2","A3","A4")
df_first_thousand_no_orf1ab$codon_aa_count = paste(df_first_thousand_no_orf1ab$codons, "(",df_first_thousand_no_orf1ab$Aa,")", round(df_first_thousand_no_orf1ab$normed_aa,4), sep = " ")

A_nuc_text = df_first_thousand_no_orf1ab[df_first_thousand_no_orf1ab$second=="A",]$codon_aa_count
C_nuc_text = df_first_thousand_no_orf1ab[df_first_thousand_no_orf1ab$second=="C",]$codon_aa_count
G_nuc_text = df_first_thousand_no_orf1ab[df_first_thousand_no_orf1ab$second=="G",]$codon_aa_count
T_nuc_text = df_first_thousand_no_orf1ab[df_first_thousand_no_orf1ab$second=="T",]$codon_aa_count
text_df = data.frame("C" = C_nuc_text,"T" = T_nuc_text,"G" = G_nuc_text,"A" = A_nuc_text)
heatmaply(
  as.matrix(heat_df),
  colors = viridis(n = 256,  option = "magma"),
  fontsize_row = 15,
  fontsize_col = 15,
  Rowv = NULL,
  Colv = NULL,
  xlab = "Second Nucleotide",
  ylab = 'First Nucleotide',
  main = 'Aa Usage for first 1000 sequences beside orf1ab',
  show_dendrogram = c(FALSE, FALSE),
  plot_method = "plotly",
  draw_cellnote = T,
  cellnote = text_df,
  cellnote_size = 13 ,
  cellnote_textposition = "middle center",
  file = "D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/R_figure/HeatMap/orf_no_orf/20.Aa_Usage_no_orf1ab_start1000.html"
)

#orf1ab codons 1000 end
A_Nuc = df_last_thousand_orf1ab[df_last_thousand_orf1ab$second=="A",]$normed_codon
C_Nuc = df_last_thousand_orf1ab[df_last_thousand_orf1ab$second=="C",]$normed_codon
G_Nuc = df_last_thousand_orf1ab[df_last_thousand_orf1ab$second=="G",]$normed_codon
T_Nuc = df_last_thousand_orf1ab[df_last_thousand_orf1ab$second=="T",]$normed_codon

heat_df = data.frame("C" = C_Nuc,"T" = T_Nuc,"G" = G_Nuc,"A" = A_Nuc)
row.names(heat_df) = c("C1","C2","C3","C4","T1","T2","T3","T4","G1","G2","G3","G4","A1","A2","A3","A4")
df_last_thousand_orf1ab$codon_aa_count = paste(df_last_thousand_orf1ab$codons, "(",df_last_thousand_orf1ab$Aa,")", round(df_last_thousand_orf1ab$normed_codon,4), sep = " ")

A_nuc_text = df_last_thousand_orf1ab[df_last_thousand_orf1ab$second=="A",]$codon_aa_count
C_nuc_text = df_last_thousand_orf1ab[df_last_thousand_orf1ab$second=="C",]$codon_aa_count
G_nuc_text = df_last_thousand_orf1ab[df_last_thousand_orf1ab$second=="G",]$codon_aa_count
T_nuc_text = df_last_thousand_orf1ab[df_last_thousand_orf1ab$second=="T",]$codon_aa_count
text_df = data.frame("C" = C_nuc_text,"T" = T_nuc_text,"G" = G_nuc_text,"A" = A_nuc_text)
heatmaply(
  as.matrix(heat_df),
  colors = viridis(n = 256,  option = "magma"),
  fontsize_row = 15,
  fontsize_col = 15,
  Rowv = NULL,
  Colv = NULL,
  xlab = "Second Nucleotide",
  ylab = 'First Nucleotide',
  main = 'Codon Usage for last 1000 sequences in ORF1ab',
  show_dendrogram = c(FALSE, FALSE),
  plot_method = "plotly",
  draw_cellnote = T,
  cellnote = text_df,
  cellnote_size = 13 ,
  cellnote_textposition = "middle center",
  file = "D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/R_figure/HeatMap/orf_no_orf/20.Codon_Usage_orf1ab_end1000.html"
)

#no orf1ab codons 1000 end
A_Nuc = df_last_thousand_no_orf1ab[df_last_thousand_no_orf1ab$second=="A",]$normed_codon
C_Nuc = df_last_thousand_no_orf1ab[df_last_thousand_no_orf1ab$second=="C",]$normed_codon
G_Nuc = df_last_thousand_no_orf1ab[df_last_thousand_no_orf1ab$second=="G",]$normed_codon
T_Nuc = df_last_thousand_no_orf1ab[df_last_thousand_no_orf1ab$second=="T",]$normed_codon

heat_df = data.frame("C" = C_Nuc,"T" = T_Nuc,"G" = G_Nuc,"A" = A_Nuc)
row.names(heat_df) = c("C1","C2","C3","C4","T1","T2","T3","T4","G1","G2","G3","G4","A1","A2","A3","A4")
df_last_thousand_no_orf1ab$codon_aa_count = paste(df_last_thousand_no_orf1ab$codons, "(",df_last_thousand_no_orf1ab$Aa,")", round(df_last_thousand_no_orf1ab$normed_codon,4), sep = " ")

A_nuc_text = df_last_thousand_no_orf1ab[df_last_thousand_no_orf1ab$second=="A",]$codon_aa_count
C_nuc_text = df_last_thousand_no_orf1ab[df_last_thousand_no_orf1ab$second=="C",]$codon_aa_count
G_nuc_text = df_last_thousand_no_orf1ab[df_last_thousand_no_orf1ab$second=="G",]$codon_aa_count
T_nuc_text = df_last_thousand_no_orf1ab[df_last_thousand_no_orf1ab$second=="T",]$codon_aa_count
text_df = data.frame("C" = C_nuc_text,"T" = T_nuc_text,"G" = G_nuc_text,"A" = A_nuc_text)
heatmaply(
  as.matrix(heat_df),
  colors = viridis(n = 256,  option = "magma"),
  fontsize_row = 15,
  fontsize_col = 15,
  Rowv = NULL,
  Colv = NULL,
  xlab = "Second Nucleotide",
  ylab = 'First Nucleotide',
  main = 'Codon Usage for last 1000 sequences beside orf1ab',
  show_dendrogram = c(FALSE, FALSE),
  plot_method = "plotly",
  draw_cellnote = T,
  cellnote = text_df,
  cellnote_size = 13 ,
  cellnote_textposition = "middle center",
  file = "D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/R_figure/HeatMap/orf_no_orf/20.Codon_Usage_no_orf1ab_end1000.html"
)
#orf1ab codons in aa 1000 end
A_Nuc = df_last_thousand_orf1ab[df_last_thousand_orf1ab$second=="A",]$normed_codons_in_aa
C_Nuc = df_last_thousand_orf1ab[df_last_thousand_orf1ab$second=="C",]$normed_codons_in_aa
G_Nuc = df_last_thousand_orf1ab[df_last_thousand_orf1ab$second=="G",]$normed_codons_in_aa
T_Nuc = df_last_thousand_orf1ab[df_last_thousand_orf1ab$second=="T",]$normed_codons_in_aa

heat_df = data.frame("C" = C_Nuc,"T" = T_Nuc,"G" = G_Nuc,"A" = A_Nuc)
row.names(heat_df) = c("C1","C2","C3","C4","T1","T2","T3","T4","G1","G2","G3","G4","A1","A2","A3","A4")
df_last_thousand_orf1ab$codon_aa_count = paste(df_last_thousand_orf1ab$codons, "(",df_last_thousand_orf1ab$Aa,")", round(df_last_thousand_orf1ab$normed_codons_in_aa,4), sep = " ")

A_nuc_text = df_last_thousand_orf1ab[df_last_thousand_orf1ab$second=="A",]$codon_aa_count
C_nuc_text = df_last_thousand_orf1ab[df_last_thousand_orf1ab$second=="C",]$codon_aa_count
G_nuc_text = df_last_thousand_orf1ab[df_last_thousand_orf1ab$second=="G",]$codon_aa_count
T_nuc_text = df_last_thousand_orf1ab[df_last_thousand_orf1ab$second=="T",]$codon_aa_count
text_df = data.frame("C" = C_nuc_text,"T" = T_nuc_text,"G" = G_nuc_text,"A" = A_nuc_text)
heatmaply(
  as.matrix(heat_df),
  colors = viridis(n = 256,  option = "magma"),
  fontsize_row = 15,
  fontsize_col = 15,
  Rowv = NULL,
  Colv = NULL,
  xlab = "Second Nucleotide",
  ylab = 'First Nucleotide',
  main = 'Codon Usage in Aa for last 1000 sequences in ORF1ab',
  show_dendrogram = c(FALSE, FALSE),
  plot_method = "plotly",
  draw_cellnote = T,
  cellnote = text_df,
  cellnote_size = 13 ,
  cellnote_textposition = "middle center",
  file = "D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/R_figure/HeatMap/orf_no_orf/20.Codon_Usage_in_Aa_orf1ab_end1000.html"
)
#no orf1ab codons in aa 1000 end
A_Nuc = df_last_thousand_no_orf1ab[df_last_thousand_no_orf1ab$second=="A",]$normed_codons_in_aa
C_Nuc = df_last_thousand_no_orf1ab[df_last_thousand_no_orf1ab$second=="C",]$normed_codons_in_aa
G_Nuc = df_last_thousand_no_orf1ab[df_last_thousand_no_orf1ab$second=="G",]$normed_codons_in_aa
T_Nuc = df_last_thousand_no_orf1ab[df_last_thousand_no_orf1ab$second=="T",]$normed_codons_in_aa

heat_df = data.frame("C" = C_Nuc,"T" = T_Nuc,"G" = G_Nuc,"A" = A_Nuc)
row.names(heat_df) = c("C1","C2","C3","C4","T1","T2","T3","T4","G1","G2","G3","G4","A1","A2","A3","A4")
df_last_thousand_no_orf1ab$codon_aa_count = paste(df_last_thousand_no_orf1ab$codons, "(",df_last_thousand_no_orf1ab$Aa,")", round(df_last_thousand_no_orf1ab$normed_codons_in_aa,4), sep = " ")

A_nuc_text = df_last_thousand_no_orf1ab[df_last_thousand_no_orf1ab$second=="A",]$codon_aa_count
C_nuc_text = df_last_thousand_no_orf1ab[df_last_thousand_no_orf1ab$second=="C",]$codon_aa_count
G_nuc_text = df_last_thousand_no_orf1ab[df_last_thousand_no_orf1ab$second=="G",]$codon_aa_count
T_nuc_text = df_last_thousand_no_orf1ab[df_last_thousand_no_orf1ab$second=="T",]$codon_aa_count
text_df = data.frame("C" = C_nuc_text,"T" = T_nuc_text,"G" = G_nuc_text,"A" = A_nuc_text)
heatmaply(
  as.matrix(heat_df),
  colors = viridis(n = 256,  option = "magma"),
  fontsize_row = 15,
  fontsize_col = 15,
  Rowv = NULL,
  Colv = NULL,
  xlab = "Second Nucleotide",
  ylab = 'First Nucleotide',
  main = 'Codon Usage in Aa for last 1000 sequences beside orf1ab',
  show_dendrogram = c(FALSE, FALSE),
  plot_method = "plotly",
  draw_cellnote = T,
  cellnote = text_df,
  cellnote_size = 13 ,
  cellnote_textposition = "middle center",
  file = "D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/R_figure/HeatMap/orf_no_orf/20.Codon_Usage_in_Aa_no_orf1ab_end1000.html"
)
#orf1ab aa 1000 end
A_Nuc = df_last_thousand_orf1ab[df_last_thousand_orf1ab$second=="A",]$normed_aa
C_Nuc = df_last_thousand_orf1ab[df_last_thousand_orf1ab$second=="C",]$normed_aa
G_Nuc = df_last_thousand_orf1ab[df_last_thousand_orf1ab$second=="G",]$normed_aa
T_Nuc = df_last_thousand_orf1ab[df_last_thousand_orf1ab$second=="T",]$normed_aa

heat_df = data.frame("C" = C_Nuc,"T" = T_Nuc,"G" = G_Nuc,"A" = A_Nuc)
row.names(heat_df) = c("C1","C2","C3","C4","T1","T2","T3","T4","G1","G2","G3","G4","A1","A2","A3","A4")
df_last_thousand_orf1ab$codon_aa_count = paste(df_last_thousand_orf1ab$codons, "(",df_last_thousand_orf1ab$Aa,")", round(df_last_thousand_orf1ab$normed_aa,4), sep = " ")

A_nuc_text = df_last_thousand_orf1ab[df_last_thousand_orf1ab$second=="A",]$codon_aa_count
C_nuc_text = df_last_thousand_orf1ab[df_last_thousand_orf1ab$second=="C",]$codon_aa_count
G_nuc_text = df_last_thousand_orf1ab[df_last_thousand_orf1ab$second=="G",]$codon_aa_count
T_nuc_text = df_last_thousand_orf1ab[df_last_thousand_orf1ab$second=="T",]$codon_aa_count
text_df = data.frame("C" = C_nuc_text,"T" = T_nuc_text,"G" = G_nuc_text,"A" = A_nuc_text)
heatmaply(
  as.matrix(heat_df),
  colors = viridis(n = 256,  option = "magma"),
  fontsize_row = 15,
  fontsize_col = 15,
  Rowv = NULL,
  Colv = NULL,
  xlab = "Second Nucleotide",
  ylab = 'First Nucleotide',
  main = 'Aa Usage for last 1000 sequences in ORF1ab',
  show_dendrogram = c(FALSE, FALSE),
  plot_method = "plotly",
  draw_cellnote = T,
  cellnote = text_df,
  cellnote_size = 13 ,
  cellnote_textposition = "middle center",
  file = "D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/R_figure/HeatMap/orf_no_orf/20.Aa_Usage_orf1ab_end1000.html"
)
#no orf1ab aa 1000 end
A_Nuc = df_last_thousand_no_orf1ab[df_last_thousand_no_orf1ab$second=="A",]$normed_aa
C_Nuc = df_last_thousand_no_orf1ab[df_last_thousand_no_orf1ab$second=="C",]$normed_aa
G_Nuc = df_last_thousand_no_orf1ab[df_last_thousand_no_orf1ab$second=="G",]$normed_aa
T_Nuc = df_last_thousand_no_orf1ab[df_last_thousand_no_orf1ab$second=="T",]$normed_aa

heat_df = data.frame("C" = C_Nuc,"T" = T_Nuc,"G" = G_Nuc,"A" = A_Nuc)
row.names(heat_df) = c("C1","C2","C3","C4","T1","T2","T3","T4","G1","G2","G3","G4","A1","A2","A3","A4")
df_last_thousand_no_orf1ab$codon_aa_count = paste(df_last_thousand_no_orf1ab$codons, "(",df_last_thousand_no_orf1ab$Aa,")", round(df_last_thousand_no_orf1ab$normed_aa,4), sep = " ")

A_nuc_text = df_last_thousand_no_orf1ab[df_last_thousand_no_orf1ab$second=="A",]$codon_aa_count
C_nuc_text = df_last_thousand_no_orf1ab[df_last_thousand_no_orf1ab$second=="C",]$codon_aa_count
G_nuc_text = df_last_thousand_no_orf1ab[df_last_thousand_no_orf1ab$second=="G",]$codon_aa_count
T_nuc_text = df_last_thousand_no_orf1ab[df_last_thousand_no_orf1ab$second=="T",]$codon_aa_count
text_df = data.frame("C" = C_nuc_text,"T" = T_nuc_text,"G" = G_nuc_text,"A" = A_nuc_text)
heatmaply(
  as.matrix(heat_df),
  colors = viridis(n = 256,  option = "magma"),
  fontsize_row = 15,
  fontsize_col = 15,
  Rowv = NULL,
  Colv = NULL,
  xlab = "Second Nucleotide",
  ylab = 'First Nucleotide',
  main = 'Aa Usage for last 1000 sequences beside orf1ab',
  show_dendrogram = c(FALSE, FALSE),
  plot_method = "plotly",
  draw_cellnote = T,
  cellnote = text_df,
  cellnote_size = 13 ,
  cellnote_textposition = "middle center",
  file = "D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/R_figure/HeatMap/orf_no_orf/20.Aa_Usage_no_orf1ab_end1000.html"
)
#last_orf - first_orf
df_first_thousand_orf1ab$l_f_normed_codons = df_last_thousand_orf1ab$normed_codon - df_first_thousand_orf1ab$normed_codon
df_first_thousand_orf1ab$l_f_normed_codons_aa = df_last_thousand_orf1ab$normed_codons_in_aa - df_first_thousand_orf1ab$normed_codons_in_aa
df_first_thousand_orf1ab$l_f_normed_aa = df_last_thousand_orf1ab$normed_aa - df_first_thousand_orf1ab$normed_aa

#last_no_orf - first_no_orf
df_first_thousand_no_orf1ab$l_f_normed_codons = df_last_thousand_no_orf1ab$normed_codon - df_first_thousand_no_orf1ab$normed_codon
df_first_thousand_no_orf1ab$l_f_normed_codons_aa = df_last_thousand_no_orf1ab$normed_codons_in_aa - df_first_thousand_no_orf1ab$normed_codons_in_aa
df_first_thousand_no_orf1ab$l_f_normed_aa = df_last_thousand_no_orf1ab$normed_aa - df_first_thousand_no_orf1ab$normed_aa

#no_orf - orf first
df_first_thousand_orf1ab$o_no_normed_codons_first = df_first_thousand_no_orf1ab$normed_codon - df_first_thousand_orf1ab$normed_codon
df_first_thousand_orf1ab$o_no_normed_codons_aa_first = df_first_thousand_no_orf1ab$normed_codons_in_aa - df_first_thousand_orf1ab$normed_codons_in_aa
df_first_thousand_orf1ab$o_no_normed_aa_first = df_first_thousand_no_orf1ab$normed_aa - df_first_thousand_orf1ab$normed_aa

#no_orf - orf last
df_last_thousand_orf1ab$o_no_normed_codons_last = df_last_thousand_no_orf1ab$normed_codon - df_last_thousand_orf1ab$normed_codon
df_last_thousand_orf1ab$o_no_normed_codons_aa_last = df_last_thousand_no_orf1ab$normed_codons_in_aa - df_last_thousand_orf1ab$normed_codons_in_aa
df_last_thousand_orf1ab$o_no_normed_aa_last = df_last_thousand_no_orf1ab$normed_aa - df_last_thousand_orf1ab$normed_aa

#codon last_orf - first_orf
A_Nuc = df_first_thousand_orf1ab[df_first_thousand_orf1ab$second=="A",]$l_f_normed_codons
C_Nuc = df_first_thousand_orf1ab[df_first_thousand_orf1ab$second=="C",]$l_f_normed_codons
G_Nuc = df_first_thousand_orf1ab[df_first_thousand_orf1ab$second=="G",]$l_f_normed_codons
T_Nuc = df_first_thousand_orf1ab[df_first_thousand_orf1ab$second=="T",]$l_f_normed_codons

heat_df = data.frame("C" = C_Nuc,"T" = T_Nuc,"G" = G_Nuc,"A" = A_Nuc)
row.names(heat_df) = c("C1","C2","C3","C4","T1","T2","T3","T4","G1","G2","G3","G4","A1","A2","A3","A4")
df_first_thousand_orf1ab$codon_aa_count = paste(df_first_thousand_orf1ab$codons, "(",df_first_thousand_orf1ab$Aa,")", round(df_first_thousand_orf1ab$l_f_normed_codons,4), sep = " ")

A_nuc_text = df_first_thousand_orf1ab[df_first_thousand_orf1ab$second=="A",]$codon_aa_count
C_nuc_text = df_first_thousand_orf1ab[df_first_thousand_orf1ab$second=="C",]$codon_aa_count
G_nuc_text = df_first_thousand_orf1ab[df_first_thousand_orf1ab$second=="G",]$codon_aa_count
T_nuc_text = df_first_thousand_orf1ab[df_first_thousand_orf1ab$second=="T",]$codon_aa_count
text_df = data.frame("C" = C_nuc_text,"T" = T_nuc_text,"G" = G_nuc_text,"A" = A_nuc_text)
heatmaply(
  as.matrix(heat_df),
  colors = viridis(n = 256,  option = "magma"),
  fontsize_row = 15,
  fontsize_col = 15,
  Rowv = NULL,
  Colv = NULL,
  xlab = "Second Nucleotide",
  ylab = 'First Nucleotide',
  main = 'Change of codon usage last-first 1000 in ORF1ab',
  show_dendrogram = c(FALSE, FALSE),
  plot_method = "plotly",
  draw_cellnote = T,
  cellnote = text_df,
  cellnote_size = 13 ,
  cellnote_textposition = "middle center",
  file = "D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/R_figure/HeatMap/orf_no_orf/Change/time/20.Change_of_codons_end-start_orf1ab.html"
)
#codon in aa last_orf - first_orf
A_Nuc = df_first_thousand_orf1ab[df_first_thousand_orf1ab$second=="A",]$l_f_normed_codons_aa
C_Nuc = df_first_thousand_orf1ab[df_first_thousand_orf1ab$second=="C",]$l_f_normed_codons_aa
G_Nuc = df_first_thousand_orf1ab[df_first_thousand_orf1ab$second=="G",]$l_f_normed_codons_aa
T_Nuc = df_first_thousand_orf1ab[df_first_thousand_orf1ab$second=="T",]$l_f_normed_codons_aa

heat_df = data.frame("C" = C_Nuc,"T" = T_Nuc,"G" = G_Nuc,"A" = A_Nuc)
row.names(heat_df) = c("C1","C2","C3","C4","T1","T2","T3","T4","G1","G2","G3","G4","A1","A2","A3","A4")
df_first_thousand_orf1ab$codon_aa_count = paste(df_first_thousand_orf1ab$codons, "(",df_first_thousand_orf1ab$Aa,")", round(df_first_thousand_orf1ab$l_f_normed_codons_aa,4), sep = " ")

A_nuc_text = df_first_thousand_orf1ab[df_first_thousand_orf1ab$second=="A",]$codon_aa_count
C_nuc_text = df_first_thousand_orf1ab[df_first_thousand_orf1ab$second=="C",]$codon_aa_count
G_nuc_text = df_first_thousand_orf1ab[df_first_thousand_orf1ab$second=="G",]$codon_aa_count
T_nuc_text = df_first_thousand_orf1ab[df_first_thousand_orf1ab$second=="T",]$codon_aa_count
text_df = data.frame("C" = C_nuc_text,"T" = T_nuc_text,"G" = G_nuc_text,"A" = A_nuc_text)
heatmaply(
  as.matrix(heat_df),
  colors = viridis(n = 256,  option = "magma"),
  fontsize_row = 15,
  fontsize_col = 15,
  Rowv = NULL,
  Colv = NULL,
  xlab = "Second Nucleotide",
  ylab = 'First Nucleotide',
  main = 'Change of codon usage among amino acids last-first 1000 in ORF1ab',
  show_dendrogram = c(FALSE, FALSE),
  plot_method = "plotly",
  draw_cellnote = T,
  cellnote = text_df,
  cellnote_size = 13 ,
  cellnote_textposition = "middle center",
  file = "D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/R_figure/HeatMap/orf_no_orf/Change/time/20.Change_of_codons_in_Aa_end-start_orf1ab.html"
)
#aa last_orf - first_orf
A_Nuc = df_first_thousand_orf1ab[df_first_thousand_orf1ab$second=="A",]$l_f_normed_aa
C_Nuc = df_first_thousand_orf1ab[df_first_thousand_orf1ab$second=="C",]$l_f_normed_aa
G_Nuc = df_first_thousand_orf1ab[df_first_thousand_orf1ab$second=="G",]$l_f_normed_aa
T_Nuc = df_first_thousand_orf1ab[df_first_thousand_orf1ab$second=="T",]$l_f_normed_aa

heat_df = data.frame("C" = C_Nuc,"T" = T_Nuc,"G" = G_Nuc,"A" = A_Nuc)
row.names(heat_df) = c("C1","C2","C3","C4","T1","T2","T3","T4","G1","G2","G3","G4","A1","A2","A3","A4")
df_first_thousand_orf1ab$codon_aa_count = paste(df_first_thousand_orf1ab$codons, "(",df_first_thousand_orf1ab$Aa,")", round(df_first_thousand_orf1ab$l_f_normed_aa,4), sep = " ")

A_nuc_text = df_first_thousand_orf1ab[df_first_thousand_orf1ab$second=="A",]$codon_aa_count
C_nuc_text = df_first_thousand_orf1ab[df_first_thousand_orf1ab$second=="C",]$codon_aa_count
G_nuc_text = df_first_thousand_orf1ab[df_first_thousand_orf1ab$second=="G",]$codon_aa_count
T_nuc_text = df_first_thousand_orf1ab[df_first_thousand_orf1ab$second=="T",]$codon_aa_count
text_df = data.frame("C" = C_nuc_text,"T" = T_nuc_text,"G" = G_nuc_text,"A" = A_nuc_text)
heatmaply(
  as.matrix(heat_df),
  colors = viridis(n = 256,  option = "magma"),
  fontsize_row = 15,
  fontsize_col = 15,
  Rowv = NULL,
  Colv = NULL,
  xlab = "Second Nucleotide",
  ylab = 'First Nucleotide',
  main = 'Change of amino acid usage last-first 1000 in ORF1ab',
  show_dendrogram = c(FALSE, FALSE),
  plot_method = "plotly",
  draw_cellnote = T,
  cellnote = text_df,
  cellnote_size = 13 ,
  cellnote_textposition = "middle center",
  file = "D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/R_figure/HeatMap/orf_no_orf/Change/time/20.Change_of_Aa_end-start_orf1ab.html"
)

#codon last_no-orf - first_no-orf
A_Nuc = df_first_thousand_no_orf1ab[df_first_thousand_no_orf1ab$second=="A",]$l_f_normed_codons
C_Nuc = df_first_thousand_no_orf1ab[df_first_thousand_no_orf1ab$second=="C",]$l_f_normed_codons
G_Nuc = df_first_thousand_no_orf1ab[df_first_thousand_no_orf1ab$second=="G",]$l_f_normed_codons
T_Nuc = df_first_thousand_no_orf1ab[df_first_thousand_no_orf1ab$second=="T",]$l_f_normed_codons

heat_df = data.frame("C" = C_Nuc,"T" = T_Nuc,"G" = G_Nuc,"A" = A_Nuc)
row.names(heat_df) = c("C1","C2","C3","C4","T1","T2","T3","T4","G1","G2","G3","G4","A1","A2","A3","A4")
df_first_thousand_no_orf1ab$codon_aa_count = paste(df_first_thousand_no_orf1ab$codons, "(",df_first_thousand_no_orf1ab$Aa,")", round(df_first_thousand_no_orf1ab$l_f_normed_codons,4), sep = " ")

A_nuc_text = df_first_thousand_no_orf1ab[df_first_thousand_no_orf1ab$second=="A",]$codon_aa_count
C_nuc_text = df_first_thousand_no_orf1ab[df_first_thousand_no_orf1ab$second=="C",]$codon_aa_count
G_nuc_text = df_first_thousand_no_orf1ab[df_first_thousand_no_orf1ab$second=="G",]$codon_aa_count
T_nuc_text = df_first_thousand_no_orf1ab[df_first_thousand_no_orf1ab$second=="T",]$codon_aa_count
text_df = data.frame("C" = C_nuc_text,"T" = T_nuc_text,"G" = G_nuc_text,"A" = A_nuc_text)
heatmaply(
  as.matrix(heat_df),
  colors = viridis(n = 256,  option = "magma"),
  fontsize_row = 15,
  fontsize_col = 15,
  Rowv = NULL,
  Colv = NULL,
  xlab = "Second Nucleotide",
  ylab = 'First Nucleotide',
  main = 'Change of codon usage last-first 1000 beside ORF1ab',
  show_dendrogram = c(FALSE, FALSE),
  plot_method = "plotly",
  draw_cellnote = T,
  cellnote = text_df,
  cellnote_size = 13 ,
  cellnote_textposition = "middle center",
  file = "D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/R_figure/HeatMap/orf_no_orf/Change/time/20.Change_of_codons_end-start_no_orf1ab.html"
)
#codon in aa last_no-orf - first_no-orf
A_Nuc = df_first_thousand_no_orf1ab[df_first_thousand_no_orf1ab$second=="A",]$l_f_normed_codons_aa
C_Nuc = df_first_thousand_no_orf1ab[df_first_thousand_no_orf1ab$second=="C",]$l_f_normed_codons_aa
G_Nuc = df_first_thousand_no_orf1ab[df_first_thousand_no_orf1ab$second=="G",]$l_f_normed_codons_aa
T_Nuc = df_first_thousand_no_orf1ab[df_first_thousand_no_orf1ab$second=="T",]$l_f_normed_codons_aa

heat_df = data.frame("C" = C_Nuc,"T" = T_Nuc,"G" = G_Nuc,"A" = A_Nuc)
row.names(heat_df) = c("C1","C2","C3","C4","T1","T2","T3","T4","G1","G2","G3","G4","A1","A2","A3","A4")
df_first_thousand_no_orf1ab$codon_aa_count = paste(df_first_thousand_no_orf1ab$codons, "(",df_first_thousand_no_orf1ab$Aa,")", round(df_first_thousand_no_orf1ab$l_f_normed_codons_aa,4), sep = " ")

A_nuc_text = df_first_thousand_no_orf1ab[df_first_thousand_no_orf1ab$second=="A",]$codon_aa_count
C_nuc_text = df_first_thousand_no_orf1ab[df_first_thousand_no_orf1ab$second=="C",]$codon_aa_count
G_nuc_text = df_first_thousand_no_orf1ab[df_first_thousand_no_orf1ab$second=="G",]$codon_aa_count
T_nuc_text = df_first_thousand_no_orf1ab[df_first_thousand_no_orf1ab$second=="T",]$codon_aa_count
text_df = data.frame("C" = C_nuc_text,"T" = T_nuc_text,"G" = G_nuc_text,"A" = A_nuc_text)
heatmaply(
  as.matrix(heat_df),
  colors = viridis(n = 256,  option = "magma"),
  fontsize_row = 15,
  fontsize_col = 15,
  Rowv = NULL,
  Colv = NULL,
  xlab = "Second Nucleotide",
  ylab = 'First Nucleotide',
  main = 'Change of codon usage among amino acids last-first 1000 beside ORF1ab',
  show_dendrogram = c(FALSE, FALSE),
  plot_method = "plotly",
  draw_cellnote = T,
  cellnote = text_df,
  cellnote_size = 13 ,
  cellnote_textposition = "middle center",
  file = "D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/R_figure/HeatMap/orf_no_orf/Change/time/20.Change_of_codons_in_Aa_end-start_no_orf1ab.html"
)
#aa last_no-orf - first_no-orf
A_Nuc = df_first_thousand_no_orf1ab[df_first_thousand_no_orf1ab$second=="A",]$l_f_normed_aa
C_Nuc = df_first_thousand_no_orf1ab[df_first_thousand_no_orf1ab$second=="C",]$l_f_normed_aa
G_Nuc = df_first_thousand_no_orf1ab[df_first_thousand_no_orf1ab$second=="G",]$l_f_normed_aa
T_Nuc = df_first_thousand_no_orf1ab[df_first_thousand_no_orf1ab$second=="T",]$l_f_normed_aa

heat_df = data.frame("C" = C_Nuc,"T" = T_Nuc,"G" = G_Nuc,"A" = A_Nuc)
row.names(heat_df) = c("C1","C2","C3","C4","T1","T2","T3","T4","G1","G2","G3","G4","A1","A2","A3","A4")
df_first_thousand_no_orf1ab$codon_aa_count = paste(df_first_thousand_no_orf1ab$codons, "(",df_first_thousand_no_orf1ab$Aa,")", round(df_first_thousand_no_orf1ab$l_f_normed_aa,4), sep = " ")

A_nuc_text = df_first_thousand_no_orf1ab[df_first_thousand_no_orf1ab$second=="A",]$codon_aa_count
C_nuc_text = df_first_thousand_no_orf1ab[df_first_thousand_no_orf1ab$second=="C",]$codon_aa_count
G_nuc_text = df_first_thousand_no_orf1ab[df_first_thousand_no_orf1ab$second=="G",]$codon_aa_count
T_nuc_text = df_first_thousand_no_orf1ab[df_first_thousand_no_orf1ab$second=="T",]$codon_aa_count
text_df = data.frame("C" = C_nuc_text,"T" = T_nuc_text,"G" = G_nuc_text,"A" = A_nuc_text)
heatmaply(
  as.matrix(heat_df),
  colors = viridis(n = 256,  option = "magma"),
  fontsize_row = 15,
  fontsize_col = 15,
  Rowv = NULL,
  Colv = NULL,
  xlab = "Second Nucleotide",
  ylab = 'First Nucleotide',
  main = 'Change of amino acid usage last-first 1000 beside ORF1ab',
  show_dendrogram = c(FALSE, FALSE),
  plot_method = "plotly",
  draw_cellnote = T,
  cellnote = text_df,
  cellnote_size = 13 ,
  cellnote_textposition = "middle center",
  file = "D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/R_figure/HeatMap/orf_no_orf/Change/time/20.Change_of_Aa_end-start_no_orf1ab.html"
)

#codon no_orf - orf first
A_Nuc = df_first_thousand_orf1ab[df_first_thousand_orf1ab$second=="A",]$o_no_normed_codons_first
C_Nuc = df_first_thousand_orf1ab[df_first_thousand_orf1ab$second=="C",]$o_no_normed_codons_first
G_Nuc = df_first_thousand_orf1ab[df_first_thousand_orf1ab$second=="G",]$o_no_normed_codons_first
T_Nuc = df_first_thousand_orf1ab[df_first_thousand_orf1ab$second=="T",]$o_no_normed_codons_first

heat_df = data.frame("C" = C_Nuc,"T" = T_Nuc,"G" = G_Nuc,"A" = A_Nuc)
row.names(heat_df) = c("C1","C2","C3","C4","T1","T2","T3","T4","G1","G2","G3","G4","A1","A2","A3","A4")
df_first_thousand_orf1ab$codon_aa_count = paste(df_first_thousand_orf1ab$codons, "(",df_first_thousand_orf1ab$Aa,")", round(df_first_thousand_orf1ab$o_no_normed_codons_first,4), sep = " ")

A_nuc_text = df_first_thousand_orf1ab[df_first_thousand_orf1ab$second=="A",]$codon_aa_count
C_nuc_text = df_first_thousand_orf1ab[df_first_thousand_orf1ab$second=="C",]$codon_aa_count
G_nuc_text = df_first_thousand_orf1ab[df_first_thousand_orf1ab$second=="G",]$codon_aa_count
T_nuc_text = df_first_thousand_orf1ab[df_first_thousand_orf1ab$second=="T",]$codon_aa_count
text_df = data.frame("C" = C_nuc_text,"T" = T_nuc_text,"G" = G_nuc_text,"A" = A_nuc_text)
heatmaply(
  as.matrix(heat_df),
  colors = viridis(n = 256,  option = "magma"),
  fontsize_row = 15,
  fontsize_col = 15,
  Rowv = NULL,
  Colv = NULL,
  xlab = "Second Nucleotide",
  ylab = 'First Nucleotide',
  main = 'Change of codon usage beside ORF1ab - in ORF1ab first 1000',
  show_dendrogram = c(FALSE, FALSE),
  plot_method = "plotly",
  draw_cellnote = T,
  cellnote = text_df,
  cellnote_size = 13 ,
  cellnote_textposition = "middle center",
  file = "D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/R_figure/HeatMap/orf_no_orf/Change/orf_no_orf/20.Change_of_codons_o_no_o_first.html"
)
#codon in aa no_orf - orf first
A_Nuc = df_first_thousand_orf1ab[df_first_thousand_orf1ab$second=="A",]$o_no_normed_codons_aa_first
C_Nuc = df_first_thousand_orf1ab[df_first_thousand_orf1ab$second=="C",]$o_no_normed_codons_aa_first
G_Nuc = df_first_thousand_orf1ab[df_first_thousand_orf1ab$second=="G",]$o_no_normed_codons_aa_first
T_Nuc = df_first_thousand_orf1ab[df_first_thousand_orf1ab$second=="T",]$o_no_normed_codons_aa_first

heat_df = data.frame("C" = C_Nuc,"T" = T_Nuc,"G" = G_Nuc,"A" = A_Nuc)
row.names(heat_df) = c("C1","C2","C3","C4","T1","T2","T3","T4","G1","G2","G3","G4","A1","A2","A3","A4")
df_first_thousand_orf1ab$codon_aa_count = paste(df_first_thousand_orf1ab$codons, "(",df_first_thousand_orf1ab$Aa,")", round(df_first_thousand_orf1ab$o_no_normed_codons_aa_first,4), sep = " ")

A_nuc_text = df_first_thousand_orf1ab[df_first_thousand_orf1ab$second=="A",]$codon_aa_count
C_nuc_text = df_first_thousand_orf1ab[df_first_thousand_orf1ab$second=="C",]$codon_aa_count
G_nuc_text = df_first_thousand_orf1ab[df_first_thousand_orf1ab$second=="G",]$codon_aa_count
T_nuc_text = df_first_thousand_orf1ab[df_first_thousand_orf1ab$second=="T",]$codon_aa_count
text_df = data.frame("C" = C_nuc_text,"T" = T_nuc_text,"G" = G_nuc_text,"A" = A_nuc_text)
heatmaply(
  as.matrix(heat_df),
  colors = viridis(n = 256,  option = "magma"),
  fontsize_row = 15,
  fontsize_col = 15,
  Rowv = NULL,
  Colv = NULL,
  xlab = "Second Nucleotide",
  ylab = 'First Nucleotide',
  main = 'Change of codon usage among amino acids beside ORF1ab - in ORF1ab first 1000',
  show_dendrogram = c(FALSE, FALSE),
  plot_method = "plotly",
  draw_cellnote = T,
  cellnote = text_df,
  cellnote_size = 13 ,
  cellnote_textposition = "middle center",
  file = "D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/R_figure/HeatMap/orf_no_orf/Change/orf_no_orf/20.Change_of_codons_in_Aa_o_no_o_first.html"
)
#aa no_orf - orf first
A_Nuc = df_first_thousand_orf1ab[df_first_thousand_orf1ab$second=="A",]$o_no_normed_aa_first
C_Nuc = df_first_thousand_orf1ab[df_first_thousand_orf1ab$second=="C",]$o_no_normed_aa_first
G_Nuc = df_first_thousand_orf1ab[df_first_thousand_orf1ab$second=="G",]$o_no_normed_aa_first
T_Nuc = df_first_thousand_orf1ab[df_first_thousand_orf1ab$second=="T",]$o_no_normed_aa_first

heat_df = data.frame("C" = C_Nuc,"T" = T_Nuc,"G" = G_Nuc,"A" = A_Nuc)
row.names(heat_df) = c("C1","C2","C3","C4","T1","T2","T3","T4","G1","G2","G3","G4","A1","A2","A3","A4")
df_first_thousand_orf1ab$codon_aa_count = paste(df_first_thousand_orf1ab$codons, "(",df_first_thousand_orf1ab$Aa,")", round(df_first_thousand_orf1ab$o_no_normed_aa_first,4), sep = " ")

A_nuc_text = df_first_thousand_orf1ab[df_first_thousand_orf1ab$second=="A",]$codon_aa_count
C_nuc_text = df_first_thousand_orf1ab[df_first_thousand_orf1ab$second=="C",]$codon_aa_count
G_nuc_text = df_first_thousand_orf1ab[df_first_thousand_orf1ab$second=="G",]$codon_aa_count
T_nuc_text = df_first_thousand_orf1ab[df_first_thousand_orf1ab$second=="T",]$codon_aa_count
text_df = data.frame("C" = C_nuc_text,"T" = T_nuc_text,"G" = G_nuc_text,"A" = A_nuc_text)
heatmaply(
  as.matrix(heat_df),
  colors = viridis(n = 256,  option = "magma"),
  fontsize_row = 15,
  fontsize_col = 15,
  Rowv = NULL,
  Colv = NULL,
  xlab = "Second Nucleotide",
  ylab = 'First Nucleotide',
  main = 'Change of amino acid usage beside ORF1ab - in ORF1ab first 1000',
  show_dendrogram = c(FALSE, FALSE),
  plot_method = "plotly",
  draw_cellnote = T,
  cellnote = text_df,
  cellnote_size = 13 ,
  cellnote_textposition = "middle center",
  file = "D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/R_figure/HeatMap/orf_no_orf/Change/orf_no_orf/20.Change_of_Aa_o_no_o_first.html"
)

#codon no_orf - orf last
A_Nuc = df_last_thousand_orf1ab[df_last_thousand_orf1ab$second=="A",]$o_no_normed_codons_last
C_Nuc = df_last_thousand_orf1ab[df_last_thousand_orf1ab$second=="C",]$o_no_normed_codons_last
G_Nuc = df_last_thousand_orf1ab[df_last_thousand_orf1ab$second=="G",]$o_no_normed_codons_last
T_Nuc = df_last_thousand_orf1ab[df_last_thousand_orf1ab$second=="T",]$o_no_normed_codons_last

heat_df = data.frame("C" = C_Nuc,"T" = T_Nuc,"G" = G_Nuc,"A" = A_Nuc)
row.names(heat_df) = c("C1","C2","C3","C4","T1","T2","T3","T4","G1","G2","G3","G4","A1","A2","A3","A4")
df_last_thousand_orf1ab$codon_aa_count = paste(df_last_thousand_orf1ab$codons, "(",df_last_thousand_orf1ab$Aa,")", round(df_last_thousand_orf1ab$o_no_normed_codons_last,4), sep = " ")

A_nuc_text = df_last_thousand_orf1ab[df_last_thousand_orf1ab$second=="A",]$codon_aa_count
C_nuc_text = df_last_thousand_orf1ab[df_last_thousand_orf1ab$second=="C",]$codon_aa_count
G_nuc_text = df_last_thousand_orf1ab[df_last_thousand_orf1ab$second=="G",]$codon_aa_count
T_nuc_text = df_last_thousand_orf1ab[df_last_thousand_orf1ab$second=="T",]$codon_aa_count
text_df = data.frame("C" = C_nuc_text,"T" = T_nuc_text,"G" = G_nuc_text,"A" = A_nuc_text)
heatmaply(
  as.matrix(heat_df),
  colors = viridis(n = 256,  option = "magma"),
  fontsize_row = 15,
  fontsize_col = 15,
  Rowv = NULL,
  Colv = NULL,
  xlab = "Second Nucleotide",
  ylab = 'First Nucleotide',
  main = 'Change of codon usage beside ORF1ab - in ORF1ab last 1000',
  show_dendrogram = c(FALSE, FALSE),
  plot_method = "plotly",
  draw_cellnote = T,
  cellnote = text_df,
  cellnote_size = 13 ,
  cellnote_textposition = "middle center",
  file = "D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/R_figure/HeatMap/orf_no_orf/Change/orf_no_orf/20.Change_of_codons_o_no_o_last.html"
)
#codon in aa no_orf - orf last
A_Nuc = df_last_thousand_orf1ab[df_last_thousand_orf1ab$second=="A",]$o_no_normed_codons_aa_last
C_Nuc = df_last_thousand_orf1ab[df_last_thousand_orf1ab$second=="C",]$o_no_normed_codons_aa_last
G_Nuc = df_last_thousand_orf1ab[df_last_thousand_orf1ab$second=="G",]$o_no_normed_codons_aa_last
T_Nuc = df_last_thousand_orf1ab[df_last_thousand_orf1ab$second=="T",]$o_no_normed_codons_aa_last

heat_df = data.frame("C" = C_Nuc,"T" = T_Nuc,"G" = G_Nuc,"A" = A_Nuc)
row.names(heat_df) = c("C1","C2","C3","C4","T1","T2","T3","T4","G1","G2","G3","G4","A1","A2","A3","A4")
df_last_thousand_orf1ab$codon_aa_count = paste(df_last_thousand_orf1ab$codons, "(",df_last_thousand_orf1ab$Aa,")", round(df_last_thousand_orf1ab$o_no_normed_codons_aa_last,4), sep = " ")

A_nuc_text = df_last_thousand_orf1ab[df_last_thousand_orf1ab$second=="A",]$codon_aa_count
C_nuc_text = df_last_thousand_orf1ab[df_last_thousand_orf1ab$second=="C",]$codon_aa_count
G_nuc_text = df_last_thousand_orf1ab[df_last_thousand_orf1ab$second=="G",]$codon_aa_count
T_nuc_text = df_last_thousand_orf1ab[df_last_thousand_orf1ab$second=="T",]$codon_aa_count
text_df = data.frame("C" = C_nuc_text,"T" = T_nuc_text,"G" = G_nuc_text,"A" = A_nuc_text)
heatmaply(
  as.matrix(heat_df),
  colors = viridis(n = 256,  option = "magma"),
  fontsize_row = 15,
  fontsize_col = 15,
  Rowv = NULL,
  Colv = NULL,
  xlab = "Second Nucleotide",
  ylab = 'First Nucleotide',
  main = 'Change of codon usage among amino acids beside ORF1ab - in ORF1ab last 1000',
  show_dendrogram = c(FALSE, FALSE),
  plot_method = "plotly",
  draw_cellnote = T,
  cellnote = text_df,
  cellnote_size = 13 ,
  cellnote_textposition = "middle center",
  file = "D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/R_figure/HeatMap/orf_no_orf/Change/orf_no_orf/20.Change_of_codons_in_Aa_o_no_o_last.html"
)
#aa no_orf - orf last
A_Nuc = df_last_thousand_orf1ab[df_last_thousand_orf1ab$second=="A",]$o_no_normed_aa_last
C_Nuc = df_last_thousand_orf1ab[df_last_thousand_orf1ab$second=="C",]$o_no_normed_aa_last
G_Nuc = df_last_thousand_orf1ab[df_last_thousand_orf1ab$second=="G",]$o_no_normed_aa_last
T_Nuc = df_last_thousand_orf1ab[df_last_thousand_orf1ab$second=="T",]$o_no_normed_aa_last

heat_df = data.frame("C" = C_Nuc,"T" = T_Nuc,"G" = G_Nuc,"A" = A_Nuc)
row.names(heat_df) = c("C1","C2","C3","C4","T1","T2","T3","T4","G1","G2","G3","G4","A1","A2","A3","A4")
df_last_thousand_orf1ab$codon_aa_count = paste(df_last_thousand_orf1ab$codons, "(",df_last_thousand_orf1ab$Aa,")", round(df_last_thousand_orf1ab$o_no_normed_aa_last,4), sep = " ")

A_nuc_text = df_last_thousand_orf1ab[df_last_thousand_orf1ab$second=="A",]$codon_aa_count
C_nuc_text = df_last_thousand_orf1ab[df_last_thousand_orf1ab$second=="C",]$codon_aa_count
G_nuc_text = df_last_thousand_orf1ab[df_last_thousand_orf1ab$second=="G",]$codon_aa_count
T_nuc_text = df_last_thousand_orf1ab[df_last_thousand_orf1ab$second=="T",]$codon_aa_count
text_df = data.frame("C" = C_nuc_text,"T" = T_nuc_text,"G" = G_nuc_text,"A" = A_nuc_text)
heatmaply(
  as.matrix(heat_df),
  colors = viridis(n = 256,  option = "magma"),
  fontsize_row = 15,
  fontsize_col = 15,
  Rowv = NULL,
  Colv = NULL,
  xlab = "Second Nucleotide",
  ylab = 'First Nucleotide',
  main = 'Change of amino acid usage beside ORF1ab - in ORF1ab last 1000',
  show_dendrogram = c(FALSE, FALSE),
  plot_method = "plotly",
  draw_cellnote = T,
  cellnote = text_df,
  cellnote_size = 13 ,
  cellnote_textposition = "middle center",
  file = "D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/R_figure/HeatMap/orf_no_orf/Change/orf_no_orf/20.Change_of_Aa_o_no_o_last.html"
)