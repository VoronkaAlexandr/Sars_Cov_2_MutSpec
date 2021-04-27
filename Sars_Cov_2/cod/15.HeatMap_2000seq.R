library(plotly)
library(data.table)
library("heatmaply")
library(seqinr)

first_thousand = read.csv('D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/data_obtained/12.Codon_Usage_begin.csv')
second_thousand = read.csv('D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/data_obtained/12.Codon_Usage_end.csv')
first_thousand[[1]] <- NULL
second_thousand[[1]] <- NULL
first_thousand[[1]] <- NULL
second_thousand[[1]] <- NULL

first_thousand$Aa = ""
second_thousand$Aa = ""

for (i in first_thousand$Codons){
  first_thousand[first_thousand$Codons == i,]$Aa = translate(s2c(i))
}
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

first_thousand$first = substring(first_thousand$Codons, 1, 1)
first_thousand$second = substring(first_thousand$Codons, 2, 2)
first_thousand$three = substring(first_thousand$Codons, 3, 3)

second_thousand$first = substring(second_thousand$Codons, 1, 1)
second_thousand$second = substring(second_thousand$Codons, 2, 2)
second_thousand$three = substring(second_thousand$Codons, 3, 3)

first_thousand$first = factor(first_thousand$first, levels = c("C","T","G","A"))
first_thousand$second = factor(first_thousand$second, levels = c("C","T","G","A"))
first_thousand$three = factor(first_thousand$three, levels = c("C","T","G","A"))

second_thousand$first = factor(second_thousand$first, levels = c("C","T","G","A"))
second_thousand$second = factor(second_thousand$second, levels = c("C","T","G","A"))
second_thousand$three = factor(second_thousand$three, levels = c("C","T","G","A"))

df_codon_first <-first_thousand[order(first_thousand$second, first_thousand$first, first_thousand$three, first_thousand$Aa),]
df_codon_last <-second_thousand[order(second_thousand$second, second_thousand$first, first_thousand$three, second_thousand$Aa),]
rownames(df_codon_first) = 1: nrow(df_codon_first)
rownames(df_codon_last) = 1: nrow(df_codon_last)

df_codon_first$normed_codons = as.numeric(df_codon_first$Mean_codon_number) / sum(as.numeric(df_codon_first$Mean_codon_number))
df_codon_first$normed_codons_in_aa = 0
df_codon_first$normed_aa = 0
for (i in unique(df_codon_first$Aa)){
  df_codon_first[df_codon_first$Aa == i,]$normed_codons_in_aa = df_codon_first[df_codon_first$Aa == i,]$Number_of_codons / sum(df_codon_first[df_codon_first$Aa == i,]$Number_of_codons)
  df_codon_first[df_codon_first$Aa == i,]$normed_aa = sum(df_codon_first[df_codon_first$Aa == i,]$Number_of_codons) / sum(df_codon_first$Number_of_codons)
}

df_codon_last$normed_codons = as.numeric(df_codon_last$Mean_codon_number) / sum(as.numeric(df_codon_last$Mean_codon_number))
df_codon_last$normed_codons_in_aa = 0
df_codon_last$normed_aa = 0
for (i in unique(df_codon_last$Aa)){
  df_codon_last[df_codon_last$Aa == i,]$normed_codons_in_aa = df_codon_last[df_codon_last$Aa == i,]$Number_of_codons / sum(df_codon_last[df_codon_last$Aa == i,]$Number_of_codons)
  df_codon_last[df_codon_last$Aa == i,]$normed_aa = sum(df_codon_last[df_codon_last$Aa == i,]$Number_of_codons) / sum(df_codon_last$Number_of_codons)
}
#codons 1000 start
A_Nuc = df_codon_first[df_codon_first$second=="A",]$normed_codons
C_Nuc = df_codon_first[df_codon_first$second=="C",]$normed_codons
G_Nuc = df_codon_first[df_codon_first$second=="G",]$normed_codons
T_Nuc = df_codon_first[df_codon_first$second=="T",]$normed_codons

heat_df = data.frame("C" = C_Nuc,"T" = T_Nuc,"G" = G_Nuc,"A" = A_Nuc)
row.names(heat_df) = c("C1","C2","C3","C4","T1","T2","T3","T4","G1","G2","G3","G4","A1","A2","A3","A4")
df_codon_first$codon_aa_count = paste(df_codon_first$Codons, "(",df_codon_first$Aa,")", round(df_codon_first$normed_codons,4), sep = " ")

A_nuc_text = df_codon_first[df_codon_first$second=="A",]$codon_aa_count
C_nuc_text = df_codon_first[df_codon_first$second=="C",]$codon_aa_count
G_nuc_text = df_codon_first[df_codon_first$second=="G",]$codon_aa_count
T_nuc_text = df_codon_first[df_codon_first$second=="T",]$codon_aa_count
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
  main = 'Codon Usage for first 1000 sequences',
  show_dendrogram = c(FALSE, FALSE),
  plot_method = "plotly",
  draw_cellnote = T,
  cellnote = text_df,
  cellnote_size = 13 ,
  cellnote_textposition = "middle center",
  file = "D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/R_figure/HeatMap/15.Codon_Usage_start1000.html"
)
#codons 1000 end
A_Nuc = df_codon_last[df_codon_last$second=="A",]$normed_codons
C_Nuc = df_codon_last[df_codon_last$second=="C",]$normed_codons
G_Nuc = df_codon_last[df_codon_last$second=="G",]$normed_codons
T_Nuc = df_codon_last[df_codon_last$second=="T",]$normed_codons

heat_df = data.frame("C" = C_Nuc,"T" = T_Nuc,"G" = G_Nuc,"A" = A_Nuc)
row.names(heat_df) = c("C1","C2","C3","C4","T1","T2","T3","T4","G1","G2","G3","G4","A1","A2","A3","A4")
df_codon_last$codon_aa_count = paste(df_codon_last$Codons, "(",df_codon_last$Aa,")", round(df_codon_last$normed_codons,4), sep = " ")

A_nuc_text = df_codon_last[df_codon_last$second=="A",]$codon_aa_count
C_nuc_text = df_codon_last[df_codon_last$second=="C",]$codon_aa_count
G_nuc_text = df_codon_last[df_codon_last$second=="G",]$codon_aa_count
T_nuc_text = df_codon_last[df_codon_last$second=="T",]$codon_aa_count
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
  main = 'Codon Usage for last 1000 sequences',
  show_dendrogram = c(FALSE, FALSE),
  plot_method = "plotly",
  draw_cellnote = T,
  cellnote = text_df,
  cellnote_size = 13 ,
  cellnote_textposition = "middle center",
  file = "D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/R_figure/HeatMap/15.Codon_Usage_end1000.html"
)
#codon in aa first 1000
A_Nuc = df_codon_first[df_codon_first$second=="A",]$normed_codons_in_aa
C_Nuc = df_codon_first[df_codon_first$second=="C",]$normed_codons_in_aa
G_Nuc = df_codon_first[df_codon_first$second=="G",]$normed_codons_in_aa
T_Nuc = df_codon_first[df_codon_first$second=="T",]$normed_codons_in_aa

heat_df = data.frame("C" = C_Nuc,"T" = T_Nuc,"G" = G_Nuc,"A" = A_Nuc)
row.names(heat_df) = c("C1","C2","C3","C4","T1","T2","T3","T4","G1","G2","G3","G4","A1","A2","A3","A4")
df_codon_first$codon_aa_count = paste(df_codon_first$Codons, "(",df_codon_first$Aa,")", round(df_codon_first$normed_codons_in_aa,4), sep = " ")

A_nuc_text = df_codon_first[df_codon_first$second=="A",]$codon_aa_count
C_nuc_text = df_codon_first[df_codon_first$second=="C",]$codon_aa_count
G_nuc_text = df_codon_first[df_codon_first$second=="G",]$codon_aa_count
T_nuc_text = df_codon_first[df_codon_first$second=="T",]$codon_aa_count
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
  main = 'Codon in amino acinds Usage for first 1000 sequences',
  show_dendrogram = c(FALSE, FALSE),
  plot_method = "plotly",
  draw_cellnote = T,
  cellnote = text_df,
  cellnote_size = 13 ,
  cellnote_textposition = "middle center",
  file = "D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/R_figure/HeatMap/15.Codon_in_aa_Usage_start1000.html"
)
#codon in aa last 1000
A_Nuc = df_codon_last[df_codon_last$second=="A",]$normed_codons_in_aa
C_Nuc = df_codon_last[df_codon_last$second=="C",]$normed_codons_in_aa
G_Nuc = df_codon_last[df_codon_last$second=="G",]$normed_codons_in_aa
T_Nuc = df_codon_last[df_codon_last$second=="T",]$normed_codons_in_aa

heat_df = data.frame("C" = C_Nuc,"T" = T_Nuc,"G" = G_Nuc,"A" = A_Nuc)
row.names(heat_df) = c("C1","C2","C3","C4","T1","T2","T3","T4","G1","G2","G3","G4","A1","A2","A3","A4")
df_codon_last$codon_aa_count = paste(df_codon_last$Codons, "(",df_codon_last$Aa,")", round(df_codon_last$normed_codons_in_aa,4), sep = " ")

A_nuc_text = df_codon_last[df_codon_last$second=="A",]$codon_aa_count
C_nuc_text = df_codon_last[df_codon_last$second=="C",]$codon_aa_count
G_nuc_text = df_codon_last[df_codon_last$second=="G",]$codon_aa_count
T_nuc_text = df_codon_last[df_codon_last$second=="T",]$codon_aa_count
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
  main = 'Codon in amino acinds Usage for last 1000 sequences',
  show_dendrogram = c(FALSE, FALSE),
  plot_method = "plotly",
  draw_cellnote = T,
  cellnote = text_df,
  cellnote_size = 13 ,
  cellnote_textposition = "middle center",
  file = "D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/R_figure/HeatMap/15.Codon_in_aa_Usage_end1000.html"
)
#aa usage first 1000
A_Nuc = df_codon_first[df_codon_first$second=="A",]$normed_aa
C_Nuc = df_codon_first[df_codon_first$second=="C",]$normed_aa
G_Nuc = df_codon_first[df_codon_first$second=="G",]$normed_aa
T_Nuc = df_codon_first[df_codon_first$second=="T",]$normed_aa

heat_df = data.frame("C" = C_Nuc,"T" = T_Nuc,"G" = G_Nuc,"A" = A_Nuc)
row.names(heat_df) = c("C1","C2","C3","C4","T1","T2","T3","T4","G1","G2","G3","G4","A1","A2","A3","A4")
df_codon_first$codon_aa_count = paste(df_codon_first$Codons, "(",df_codon_first$Aa,")", round(df_codon_first$normed_aa,4), sep = " ")

A_nuc_text = df_codon_first[df_codon_first$second=="A",]$codon_aa_count
C_nuc_text = df_codon_first[df_codon_first$second=="C",]$codon_aa_count
G_nuc_text = df_codon_first[df_codon_first$second=="G",]$codon_aa_count
T_nuc_text = df_codon_first[df_codon_first$second=="T",]$codon_aa_count
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
  main = 'Amino acinds Usage for first 1000 sequences',
  show_dendrogram = c(FALSE, FALSE),
  plot_method = "plotly",
  draw_cellnote = T,
  cellnote = text_df,
  cellnote_size = 13 ,
  cellnote_textposition = "middle center",
  file = "D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/R_figure/HeatMap/15.Aa_Usage_start1000.html"
)
#aa usage last 1000
A_Nuc = df_codon_last[df_codon_last$second=="A",]$normed_aa
C_Nuc = df_codon_last[df_codon_last$second=="C",]$normed_aa
G_Nuc = df_codon_last[df_codon_last$second=="G",]$normed_aa
T_Nuc = df_codon_last[df_codon_last$second=="T",]$normed_aa

heat_df = data.frame("C" = C_Nuc,"T" = T_Nuc,"G" = G_Nuc,"A" = A_Nuc)
row.names(heat_df) = c("C1","C2","C3","C4","T1","T2","T3","T4","G1","G2","G3","G4","A1","A2","A3","A4")
df_codon_last$codon_aa_count = paste(df_codon_last$Codons, "(",df_codon_last$Aa,")", round(df_codon_last$normed_aa,4), sep = " ")

A_nuc_text = df_codon_last[df_codon_last$second=="A",]$codon_aa_count
C_nuc_text = df_codon_last[df_codon_last$second=="C",]$codon_aa_count
G_nuc_text = df_codon_last[df_codon_last$second=="G",]$codon_aa_count
T_nuc_text = df_codon_last[df_codon_last$second=="T",]$codon_aa_count
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
  main = 'Amino acinds Usage for last 1000 sequences',
  show_dendrogram = c(FALSE, FALSE),
  plot_method = "plotly",
  draw_cellnote = T,
  cellnote = text_df,
  cellnote_size = 13 ,
  cellnote_textposition = "middle center",
  file = "D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/R_figure/HeatMap/15.Aa_Usage_end1000.html"
)
#last-first
df_codon_first$l_f_normed_codons = df_codon_last$normed_codons - df_codon_first$normed_codons
df_codon_first$l_f_normed_codons_aa = df_codon_last$normed_codons_in_aa - df_codon_first$normed_codons_in_aa
df_codon_first$l_f_normed_aa = df_codon_last$normed_aa - df_codon_first$normed_aa
#codon
A_Nuc = df_codon_first[df_codon_first$second=="A",]$l_f_normed_codons
C_Nuc = df_codon_first[df_codon_first$second=="C",]$l_f_normed_codons
G_Nuc = df_codon_first[df_codon_first$second=="G",]$l_f_normed_codons
T_Nuc = df_codon_first[df_codon_first$second=="T",]$l_f_normed_codons

heat_df = data.frame("C" = C_Nuc,"T" = T_Nuc,"G" = G_Nuc,"A" = A_Nuc)
row.names(heat_df) = c("C1","C2","C3","C4","T1","T2","T3","T4","G1","G2","G3","G4","A1","A2","A3","A4")
df_codon_first$codon_aa_count = paste(df_codon_first$Codons, "(",df_codon_first$Aa,")", round(df_codon_first$l_f_normed_codons,4), sep = " ")

A_nuc_text = df_codon_first[df_codon_first$second=="A",]$codon_aa_count
C_nuc_text = df_codon_first[df_codon_first$second=="C",]$codon_aa_count
G_nuc_text = df_codon_first[df_codon_first$second=="G",]$codon_aa_count
T_nuc_text = df_codon_first[df_codon_first$second=="T",]$codon_aa_count
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
  main = 'Change of codon usage',
  show_dendrogram = c(FALSE, FALSE),
  plot_method = "plotly",
  draw_cellnote = T,
  cellnote = text_df,
  cellnote_size = 13 ,
  cellnote_textposition = "middle center",
  file = "D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/R_figure/HeatMap/15.Change_of_codons_end-start.html"
)
#codon in aa
A_Nuc = df_codon_first[df_codon_first$second=="A",]$l_f_normed_codons_aa
C_Nuc = df_codon_first[df_codon_first$second=="C",]$l_f_normed_codons_aa
G_Nuc = df_codon_first[df_codon_first$second=="G",]$l_f_normed_codons_aa
T_Nuc = df_codon_first[df_codon_first$second=="T",]$l_f_normed_codons_aa

heat_df = data.frame("C" = C_Nuc,"T" = T_Nuc,"G" = G_Nuc,"A" = A_Nuc)
row.names(heat_df) = c("C1","C2","C3","C4","T1","T2","T3","T4","G1","G2","G3","G4","A1","A2","A3","A4")
df_codon_first$codon_aa_count = paste(df_codon_first$Codons, "(",df_codon_first$Aa,")", round(df_codon_first$l_f_normed_codons_aa,4), sep = " ")

A_nuc_text = df_codon_first[df_codon_first$second=="A",]$codon_aa_count
C_nuc_text = df_codon_first[df_codon_first$second=="C",]$codon_aa_count
G_nuc_text = df_codon_first[df_codon_first$second=="G",]$codon_aa_count
T_nuc_text = df_codon_first[df_codon_first$second=="T",]$codon_aa_count
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
  main = 'Change of codon usage among amino acids',
  show_dendrogram = c(FALSE, FALSE),
  plot_method = "plotly",
  draw_cellnote = T,
  cellnote = text_df,
  cellnote_size = 13 ,
  cellnote_textposition = "middle center",
  file = "D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/R_figure/HeatMap/15.Change_of_codons_in_aa_end-start.html"
)
#aa
A_Nuc = df_codon_first[df_codon_first$second=="A",]$l_f_normed_aa
C_Nuc = df_codon_first[df_codon_first$second=="C",]$l_f_normed_aa
G_Nuc = df_codon_first[df_codon_first$second=="G",]$l_f_normed_aa
T_Nuc = df_codon_first[df_codon_first$second=="T",]$l_f_normed_aa

heat_df = data.frame("C" = C_Nuc,"T" = T_Nuc,"G" = G_Nuc,"A" = A_Nuc)
row.names(heat_df) = c("C1","C2","C3","C4","T1","T2","T3","T4","G1","G2","G3","G4","A1","A2","A3","A4")
df_codon_first$codon_aa_count = paste(df_codon_first$Codons, "(",df_codon_first$Aa,")", round(df_codon_first$l_f_normed_aa,4), sep = " ")

A_nuc_text = df_codon_first[df_codon_first$second=="A",]$codon_aa_count
C_nuc_text = df_codon_first[df_codon_first$second=="C",]$codon_aa_count
G_nuc_text = df_codon_first[df_codon_first$second=="G",]$codon_aa_count
T_nuc_text = df_codon_first[df_codon_first$second=="T",]$codon_aa_count
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
  main = 'Change of amino acid usage',
  show_dendrogram = c(FALSE, FALSE),
  plot_method = "plotly",
  draw_cellnote = T,
  cellnote = text_df,
  cellnote_size = 13 ,
  cellnote_textposition = "middle center",
  file = "D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/R_figure/HeatMap/15.Change_of_aa_end-start.html"
)

four_vir_start = df_codon_first[df_codon_first$Aa != 'F' & df_codon_first$Aa != 'I' & df_codon_first$Aa != 'M' & df_codon_first$Aa != 'W' & df_codon_first$Aa != '*_TA' & df_codon_first$Aa != '*_TG' & df_codon_first$Aa != 'C' & df_codon_first$Aa != 'H' & df_codon_first$Aa != 'Q' & df_codon_first$Aa != 'Y' & df_codon_first$Aa != 'D' & df_codon_first$Aa != 'E' & df_codon_first$Aa != 'N' & df_codon_first$Aa != 'K' & df_codon_first$Codons != 'AGC' & df_codon_first$Codons != 'AGT' & df_codon_first$Codons != 'TTG' & df_codon_first$Codons != 'TTA' & df_codon_first$Codons != 'AGG' & df_codon_first$Codons != 'AGA',]
four_vir_1000_start = aggregate(four_vir_start$Mean_codon_number, by = list(four_vir_start$three), FUN = sum);
names(four_vir_1000_start) = c('Nucleotides','four_fold_1000_start')

four_vir_end = df_codon_last[df_codon_last$Aa != 'F' & df_codon_last$Aa != 'I' & df_codon_last$Aa != 'M' & df_codon_last$Aa != 'W' & df_codon_last$Aa != '*_TA' & df_codon_last$Aa != '*_TG' & df_codon_last$Aa != 'C' & df_codon_last$Aa != 'H' & df_codon_last$Aa != 'Q' & df_codon_last$Aa != 'Y' & df_codon_last$Aa != 'D' & df_codon_last$Aa != 'E' & df_codon_last$Aa != 'N' & df_codon_last$Aa != 'K' & df_codon_last$Codons != 'AGC' & df_codon_last$Codons != 'AGT' & df_codon_last$Codons != 'TTG' & df_codon_last$Codons != 'TTA' & df_codon_last$Codons != 'AGG' & df_codon_last$Codons != 'AGA',]
four_vir_1000_end = aggregate(four_vir_end$Mean_codon_number, by = list(four_vir_end$three), FUN = sum);
names(four_vir_1000_end) = c('Nucleotides','four_fold_1000_end')

four_vir_1000_start$Normed_mean_four_fold_1000_start = four_vir_1000_start$four_fold_1000_start / sum(four_vir_1000_start$four_fold_1000_start)
four_vir_1000_end$Normed_mean_four_fold_1000_end = four_vir_1000_end$four_fold_1000_end / sum(four_vir_1000_end$four_fold_1000_end)

nucl = c('C', 'T', 'G', 'A')
normed_four_vir_1000_start = data.frame(nucl)
names(normed_four_vir_1000_start) = 'Nucleotides'
normed_four_vir_1000_start = merge(normed_four_vir_1000_start, four_vir_1000_start, by = 'Nucleotides')
normed_four_vir_1000_start = merge(normed_four_vir_1000_start, four_vir_1000_end, by = 'Nucleotides')

syn_1000_start = aggregate(df_codon_first$Mean_codon_number, by = list(df_codon_first$three), FUN = sum);
syn_1000_end = aggregate(df_codon_last$Mean_codon_number, by = list(df_codon_last$three), FUN = sum);
names(syn_1000_start) = c('Nucleotides', 'syn_1000_start')
names(syn_1000_end) = c('Nucleotides', 'syn_1000_end')

syn_1000_start$Normed_mean_syn_1000_start = syn_1000_start$syn_1000_start / sum(syn_1000_start$syn_1000_start)
syn_1000_end$Normed_mean_syn_1000_end = syn_1000_end$syn_1000_end / sum(syn_1000_end$syn_1000_end)

normed_four_vir_1000_start = merge(normed_four_vir_1000_start, syn_1000_start, by = 'Nucleotides')
normed_four_vir_1000_start = merge(normed_four_vir_1000_start, syn_1000_end, by = 'Nucleotides')

nuc_change_view = normed_four_vir_1000_start[c(1, 3, 5, 7, 9)]
write.csv(nuc_change_view, 'D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/data_obtained/15.Nuc_Usage_end_start.csv')
