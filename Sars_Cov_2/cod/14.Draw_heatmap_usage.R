library(plotly)
library(data.table)
library("heatmaply")

uhani_codon_usage = read.csv('D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/data_obtained/Usage/13.Codon_usage_NC_045512.2.csv')
uhani_codon_usage$first = substring(uhani_codon_usage$Codons, 1, 1)
uhani_codon_usage$second = substring(uhani_codon_usage$Codons, 2, 2)
df <-uhani_codon_usage[order(uhani_codon_usage$second, uhani_codon_usage$first, uhani_codon_usage$Aacids),]
rownames(df) = 1: nrow(df)
A_nuc = df[df$second=="A",]$normed_number
C_nuc = df[df$second=="C",]$normed_number
G_Nuc = df[df$second=="G",]$normed_number
T_Nuc = df[df$second=="T",]$normed_number
heat_df = data.frame("A" = A_nuc,"C" = C_nuc,"G" = G_Nuc,"T" = T_Nuc)
row.names(heat_df) = c("A1","A2","A3","A4","C1","C2","C3","C4","G1","G2","G3","G4","T1","T2","T3","T4")

heatmaply(
  as.matrix(heat_df),
  Rowv = NULL,
  Colv = NULL,
  xlab = "First Nucleotide",
  ylab = 'Second Nucleotide',
  main = 'Codon Usage NC_045512.2',
  show_dendrogram = c(FALSE, FALSE),
  plot_method = "plotly",
  draw_cellnote = T,
  cellnote_textposition = "middle right",
  file = "D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/R_figure/HeatMap/14.Codon_Usage_Uhani.html"
)
df$Aa_normed = 0
df$forAa = 0
for (i in unique(df$Aacids)){
  df[df$Aacids == i,]$Aa_normed = df[df$Aacids == i,]$Number_of_codons / sum(df[df$Aacids == i,]$Number_of_codons)
  df[df$Aacids == i,]$forAa = sum(df[df$Aacids == i,]$Number_of_codons) / sum(df$Number_of_codons)
}
A_nuc = df[df$second=="A",]$Aa_normed
C_nuc = df[df$second=="C",]$Aa_normed
G_Nuc = df[df$second=="G",]$Aa_normed
T_Nuc = df[df$second=="T",]$Aa_normed

heat_df_aa = data.frame("A" = A_nuc,"C" = C_nuc,"G" = G_Nuc,"T" = T_Nuc)
row.names(heat_df_aa) = c("A1","A2","A3","A4","C1","C2","C3","C4","G1","G2","G3","G4","T1","T2","T3","T4")

heatmaply(
  as.matrix(heat_df_aa),
  Rowv = NULL,
  Colv = NULL,
  xlab = "First Nucleotide",
  ylab = 'Second Nucleotide',
  main = 'Codon Usage among Amino Acids NC_045512.2',
  show_dendrogram = c(FALSE, FALSE),
  plot_method = "plotly",
  draw_cellnote = T,
  cellnote_textposition = "middle right",
  file = "D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/R_figure/HeatMap/14.Codon_Usage_among_aa_Uhani.html"
)

A_nuc = df[df$second=="A",]$forAa
C_nuc = df[df$second=="C",]$forAa
G_Nuc = df[df$second=="G",]$forAa
T_Nuc = df[df$second=="T",]$forAa

heat_df_aa_only = data.frame("A" = A_nuc,"C" = C_nuc,"G" = G_Nuc,"T" = T_Nuc)
row.names(heat_df_aa_only) = c("A1","A2","A3","A4","C1","C2","C3","C4","G1","G2","G3","G4","T1","T2","T3","T4")

heatmaply(
  as.matrix(heat_df_aa_only),
  Rowv = NULL,
  Colv = NULL,
  xlab = "First Nucleotide",
  ylab = 'Second Nucleotide',
  main = 'Amino Acids usage NC_045512.2',
  show_dendrogram = c(FALSE, FALSE),
  plot_method = "plotly",
  draw_cellnote = T,
  cellnote_textposition = "middle right",
  file = "D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/R_figure/HeatMap/14.AA_Usage_Uhani.html"
)