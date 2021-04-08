rm(list=ls(all=TRUE))

ColorA = rgb(1,0,0,0.3)
ColorT = rgb(0,1,0,0.3)
ColorG = rgb(0,0,1,0.3)
ColorC = rgb(0,1,1,0.3)

Final = read.table("../../Sars_Cov_2/R_function/data/11.A.FromMutSpecToNucContent.R.FinalTable.txt", header = TRUE)


###: PLOT
# real data for cold-water fishes: A ~ 0.24, C ~ 0.08, G ~ 0.34, T ~ 0.34 # total = 1 0.24+0.08+0.34+0.34
# real data for warm-water fishes: A ~ 0.20, C ~ 0.06, G ~ 0.36, T ~ 0.38 # total = 1
# summary(Final$InitGenome)
pdf("../../Sars_Cov_2/R_function/data/11.B.FromMutSpecToNucContent.R.FinalTable.pdf",  width = 20, height = 10) # dev.off()
# par(mfrow=c(2,2))
for (InitGenome in 1:10)
{
  plot(Final[Final$InitGenome == InitGenome,]$Gener,Final[Final$InitGenome == InitGenome,]$FrA, ylim = c(0,1), pch = 16, col = ColorA, main = '', xlab = '', ylab = ''); par(new=TRUE)
  plot(Final[Final$InitGenome == InitGenome,]$Gener,Final[Final$InitGenome == InitGenome,]$FrT, ylim = c(0,1), pch = 16, col = ColorT, main = '', xlab = '', ylab = ''); par(new=TRUE)
  plot(Final[Final$InitGenome == InitGenome,]$Gener,Final[Final$InitGenome == InitGenome,]$FrG, ylim = c(0,1), pch = 16, col = ColorG, main = '', xlab = '', ylab = ''); par(new=TRUE)
  plot(Final[Final$InitGenome == InitGenome,]$Gener,Final[Final$InitGenome == InitGenome,]$FrC, ylim = c(0,1), pch = 16, col = ColorC, main = 'covid, neutral equilibrium', xlab = '', ylab = ''); par(new=TRUE)
}
# observed whole genome nucleotides
A	= 8734
C	= 5360
G	= 5735
T	= 9426
Sum = A+C+G+T; A = A/Sum; C = C/Sum; G = G/Sum; T = T/Sum; 
abline(h = A, col = ColorA, lwd = 5, lt = 2)
abline(h = T, col = ColorT, lwd = 5, lt = 2)
abline(h = G, col = ColorG, lwd = 5, lt = 2)
abline(h = C, col = ColorC, lwd = 5, lt = 2)

# observed syn nucleotides
A	= 2744
C	= 1521
G	= 1236
T	= 4250
Sum = A+C+G+T; A = A/Sum; C = C/Sum; G = G/Sum; T = T/Sum; 
abline(h = A, col = ColorA, lwd = 5, lt = 1)
abline(h = T, col = ColorT, lwd = 5, lt = 1)
abline(h = G, col = ColorG, lwd = 5, lt = 1)
abline(h = C, col = ColorC, lwd = 5, lt = 1)

legend(1,0.9, legend = c('A','T','G','C'), col = c(ColorA, ColorT, ColorG, ColorC), pch = 16)

dev.off()

### make final boxplots and print out summary (mean, median)

### IF EQUILIBRIUM SENSITIVE TO STARTING CONDITIONS?
### IF EQUILIBRIUM IS SENSITIVE TO MUTSPEC (cold versus warm fishes)
### how the same derive with analogous functions (Lynch, Valerian, Stepan)


