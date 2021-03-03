rm(list=ls(all=TRUE))

ann = read.csv("../../Sars_Cov_2/data_obtained/ideal_table.csv")
names(ann)
table(ann$GenName)
table(ann$RefAa) # it seems too much '*'
table(ann[ann$RefAa == '*',]$GenName) # everythere is 9 (3*3) - GOOD! BUT ORF1ab is NOT GOOD
# E      M      N  ORF10 ORF1ab  ORF3a   ORF6  ORF7b   ORF8      S 
# 9      9      9      9   1548      9      9      9      9      9 


