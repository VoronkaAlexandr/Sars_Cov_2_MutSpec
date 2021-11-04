library(ggplot2)
library(plotly)
library(reshape2)
library(seqinr)

ref_seq = read.fasta('D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/data/Covid_ref.fasta',forceDNAtolower = FALSE)

df_it = read.csv('D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/new_data/data/U_ideal_table.csv')
ff_aa = c('S', 'P', 'S_UC', 'A', 'T', 'L_CU', 'V', 'R_CG')
ff_ref_pos = unique(df_it[df_it$RefAa %in% ff_aa & df_it$NucInCodon == 3,]$Pos)
ff_ref_seq = c()
for (i in ff_ref_pos){
  ff_ref_seq <- c(ff_ref_seq, ref_seq$NC_045512.2[i])
}
ff_ref_nuc = data.frame(matrix(vector(), 1, 4,
                               dimnames=list(c(), c('A', 'U', 'G', 'C'))),
                        stringsAsFactors=F)
ff_ref_nuc$A = table(ff_ref_seq)[1][[1]]
ff_ref_nuc$C = table(ff_ref_seq)[2][[1]]
ff_ref_nuc$G = table(ff_ref_seq)[3][[1]]
ff_ref_nuc$U = table(ff_ref_seq)[4][[1]]
norm_ff_ref_nuc = ff_ref_nuc
norm_ff_ref_nuc$A = ff_ref_nuc$A / rowSums(ff_ref_nuc)
norm_ff_ref_nuc$U = ff_ref_nuc$U / rowSums(ff_ref_nuc)
norm_ff_ref_nuc$G = ff_ref_nuc$G / rowSums(ff_ref_nuc)
norm_ff_ref_nuc$C = ff_ref_nuc$C / rowSums(ff_ref_nuc)
norm_ff_ref_nuc$date = 'ref'
norm_ff_ref_nuc.m <- melt(norm_ff_ref_nuc, id.var = "date")

ref_seq = ref_seq$NC_045512.2[265:29674]
ref_nuc = data.frame(matrix(vector(), 1, 4,
                            dimnames=list(c(), c('A', 'U', 'G', 'C'))),
                     stringsAsFactors=F)
ref_nuc$A = table(ref_seq)[1][[1]]
ref_nuc$C = table(ref_seq)[2][[1]]
ref_nuc$G = table(ref_seq)[3][[1]]
ref_nuc$U = table(ref_seq)[4][[1]]
norm_ref_nuc = ref_nuc
norm_ref_nuc$A = ref_nuc$A / rowSums(ref_nuc)
norm_ref_nuc$U = ref_nuc$U / rowSums(ref_nuc)
norm_ref_nuc$G = ref_nuc$G / rowSums(ref_nuc)
norm_ref_nuc$C = ref_nuc$C / rowSums(ref_nuc)
norm_ref_nuc$date = 'ref'
norm_ref_nuc.m <- melt(norm_ref_nuc, id.var = "date")


#Nuc

nuc_start = read.csv('D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/new_data/data_obtained/19.Nuc_Start.csv')
norm_nuc_start = nuc_start
norm_nuc_start$A = nuc_start$A / rowSums(nuc_start)
norm_nuc_start$U = nuc_start$U / rowSums(nuc_start)
norm_nuc_start$G = nuc_start$G / rowSums(nuc_start)
norm_nuc_start$C = nuc_start$C / rowSums(nuc_start)

nuc_inter = read.csv('D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/new_data/data_obtained/19.Nuc_Inter.csv')
nuc_inter <- nuc_inter[-5765, ]
norm_nuc_inter = nuc_inter
norm_nuc_inter$A = nuc_inter$A / rowSums(nuc_inter)
norm_nuc_inter$U = nuc_inter$U / rowSums(nuc_inter)
norm_nuc_inter$G = nuc_inter$G / rowSums(nuc_inter)
norm_nuc_inter$C = nuc_inter$C / rowSums(nuc_inter)

nuc_end = read.csv('D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/new_data/data_obtained/19.Nuc_End.csv')

norm_nuc_end = nuc_end
norm_nuc_end$A = nuc_end$A / rowSums(nuc_end)
norm_nuc_end$U = nuc_end$U / rowSums(nuc_end)
norm_nuc_end$G = nuc_end$G / rowSums(nuc_end)
norm_nuc_end$C = nuc_end$C / rowSums(nuc_end)

norm_nuc_start$date = 'start'
norm_nuc_inter$date = 'inter'
norm_nuc_end$date = 'end'

nuc = rbind(norm_nuc_start, norm_nuc_inter)
nuc = rbind(nuc, norm_nuc_end)

nuc.m <- melt(nuc, id.var = "date")


p <- ggplot(data = nuc.m, aes(x=variable, y=value)) + 
  geom_boxplot(aes(fill=factor(date, level = c('start', 'inter', 'end'))), )
p <- p + guides(fill=guide_legend(title="Sequencing Time")) + facet_wrap( ~ variable, scales="free") +
  geom_hline(data = norm_ref_nuc.m, aes(yintercept = value))
ggsave('D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/new_data/figures/19.Nuc_changes.svg', plot = p, device = 'svg',  limitsize = F)


#FF

ff_start = read.csv('D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/new_data/data_obtained/19.FF_Start.csv')
norm_ff_start = ff_start
norm_ff_start$A = ff_start$A / rowSums(ff_start)
norm_ff_start$U = ff_start$U / rowSums(ff_start)
norm_ff_start$G = ff_start$G / rowSums(ff_start)
norm_ff_start$C = ff_start$C / rowSums(ff_start)

ff_inter = read.csv('D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/new_data/data_obtained/19.FF_Inter.csv')
ff_inter <- ff_inter[-5765, ]
norm_ff_inter = ff_inter
norm_ff_inter$A = ff_inter$A / rowSums(ff_inter)
norm_ff_inter$U = ff_inter$U / rowSums(ff_inter)
norm_ff_inter$G = ff_inter$G / rowSums(ff_inter)
norm_ff_inter$C = ff_inter$C / rowSums(ff_inter)

ff_end = read.csv('D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/new_data/data_obtained/19.FF_End.csv')

norm_ff_end = ff_end
norm_ff_end$A = ff_end$A / rowSums(ff_end)
norm_ff_end$U = ff_end$U / rowSums(ff_end)
norm_ff_end$G = ff_end$G / rowSums(ff_end)
norm_ff_end$C = ff_end$C / rowSums(ff_end)

norm_ff_start$date = 'start'
norm_ff_inter$date = 'inter'
norm_ff_end$date = 'end'

ff = rbind(norm_ff_start, norm_ff_inter)
ff = rbind(ff, norm_ff_end)

ff.m <- melt(ff, id.var = "date")

p <- ggplot(data = ff.m, aes(x=variable, y=value)) + 
  geom_boxplot(aes(fill=factor(date, level = c('start', 'inter', 'end'))), )
p <- p + guides(fill=guide_legend(title="Sequencing Time")) + facet_wrap( ~ variable, scales="free") +
  geom_hline(data = norm_ff_ref_nuc.m, aes(yintercept = value))
ggsave('D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/new_data/figures/19.FF_changes.svg', plot = p, device = 'svg',  limitsize = F)
