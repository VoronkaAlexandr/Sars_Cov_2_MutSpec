ann = read.csv("D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/data_obtained/ideal_table.csv")

second_str = read.csv('D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/data/secondary_structure_on_genome.csv')
ann$IsStem = NaN
ann$SsPairs = NaN
for (i in 1:nrow(second_str)){
  ann[ann$Pos == second_str$Pos[i], ]$IsStem = second_str$IsStem[i]
  ann[ann$Pos == second_str$Pos[i], ]$SsPairs = second_str$SsPairs[i]
}
write.csv(ann, "D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/data_obtained/ideal_table.csv")