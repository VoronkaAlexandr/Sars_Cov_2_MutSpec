rm(list=ls(all=TRUE))

ColorA = rgb(1,0,0,0.3)
ColorT = rgb(0,1,0,0.3)
ColorG = rgb(0,0,1,0.3)
ColorC = rgb(0,1,1,0.3)

Final = read.table("../../Sars_Cov_2/R_function/data/11.A.FromMutSpecToNucContent.R.FinalTable.txt", header = TRUE)

###: PLOT
# summary(Final$InitGenome)
pdf("../../Sars_Cov_2/R_function/data/11.B.FromMutSpecToNucContent.R.FinalTable.pdf",  width = 20, height = 10) # dev.off()
# par(mfrow=c(2,2))
for (InitGenome in 1:10)
{
  plot(Final[Final$InitGenome == InitGenome,]$Gener,Final[Final$InitGenome == InitGenome,]$FrA, ylim = c(0,1), pch = 16, col = ColorA, main = '', xlab = '', ylab = '', yaxt='n', xaxt = 'n'); par(new=TRUE)
  plot(Final[Final$InitGenome == InitGenome,]$Gener,Final[Final$InitGenome == InitGenome,]$FrT, ylim = c(0,1), pch = 16, col = ColorT, main = '', xlab = '', ylab = '', yaxt='n', xaxt = 'n'); par(new=TRUE)
  plot(Final[Final$InitGenome == InitGenome,]$Gener,Final[Final$InitGenome == InitGenome,]$FrG, ylim = c(0,1), pch = 16, col = ColorG, main = '', xlab = '', ylab = '', yaxt='n', xaxt = 'n'); par(new=TRUE)
  plot(Final[Final$InitGenome == InitGenome,]$Gener,Final[Final$InitGenome == InitGenome,]$FrC, ylim = c(0,1), pch = 16, col = ColorC, main = '', xlab = '', ylab = '', yaxt='n', xaxt = 'n'); par(new=TRUE)
  if (InitGenome == 10) {plot(Final[Final$InitGenome == InitGenome,]$Gener,Final[Final$InitGenome == InitGenome,]$FrC, ylim = c(0,1), pch = 16, col = ColorC, main = 'covid, neutral equilibrium', xlab = 'PseudoTime', ylab = 'NucleotideFrequencies')}; par(new=TRUE)
}

###### data from Valerian:
# input for Valerian: https://github.com/VoronkaAlexandr/Sars_Cov_2_MutSpec/blob/main/Sars_Cov_2/R_function/data/07.MutSpec12_ForFullGenome.csv
#Костя, для коронавируса получаются следующие предельные значения:
A= 0.1401067180
G= 0.02243383779
T= 0.7611186393
C= 0.0763408032
X = 505000

plot(X,A,col = "red", ylim = c(0,1), xlim = c(0,500000), pch = 15, cex = 3,  xlab = '', ylab = '', yaxt='n', xaxt = 'n'); par(new=TRUE)
plot(X,G,col = "blue", ylim = c(0,1), xlim = c(0,500000), pch = 15, cex = 3,  xlab = '', ylab = '', yaxt='n', xaxt = 'n'); par(new=TRUE)
plot(X,T,col = "green", ylim = c(0,1), xlim = c(0,500000), pch = 15, cex = 3,  xlab = '', ylab = '', yaxt='n', xaxt = 'n'); par(new=TRUE)
plot(X,C,col = "cyan", ylim = c(0,1), xlim = c(0,500000), pch = 15, cex = 3,  xlab = '', ylab = '', yaxt='n', xaxt = 'n'); par(new=TRUE)


###### observed whole genome nucleotides
A	= 8734
C	= 5360
G	= 5735
T	= 9426
Sum = A+C+G+T; A = A/Sum; C = C/Sum; G = G/Sum; T = T/Sum; 
#abline(h = A, col = ColorA, lwd = 5, lt = 2)
#abline(h = T, col = ColorT, lwd = 5, lt = 2)
#abline(h = G, col = ColorG, lwd = 5, lt = 2)
#abline(h = C, col = ColorC, lwd = 5, lt = 2)

# observed syn nucleotides
A	= 2744
C	= 1521
G	= 1236
T	= 4250
Sum = A+C+G+T; A = A/Sum; C = C/Sum; G = G/Sum; T = T/Sum; 
X = 20000

plot(X,A,col = "red", ylim = c(0,1), xlim = c(0,500000), pch = 19, cex = 3,  xlab = '', ylab = '', yaxt='n', xaxt = 'n'); par(new=TRUE)
plot(X,G,col = "blue", ylim = c(0,1), xlim = c(0,500000), pch = 19, cex = 3,  xlab = '', ylab = '', yaxt='n', xaxt = 'n'); par(new=TRUE)
plot(X,T,col = "green", ylim = c(0,1), xlim = c(0,500000), pch = 19, cex = 3,  xlab = '', ylab = '', yaxt='n', xaxt = 'n'); par(new=TRUE)
plot(X,C,col = "cyan", ylim = c(0,1), xlim = c(0,500000), pch = 19, cex = 3,  xlab = '', ylab = '', yaxt='n', xaxt = 'n'); par(new=TRUE)

legend(1,0.9, legend = c('A','T','G','C'), col = c(ColorA, ColorT, ColorG, ColorC), pch = 16)

dev.off()

### make final boxplots and print out summary (mean, median)

### IF EQUILIBRIUM SENSITIVE TO STARTING CONDITIONS?
### IF EQUILIBRIUM IS SENSITIVE TO MUTSPEC (cold versus warm fishes)
### how the same derive with analogous functions (Lynch, Valerian, Stepan)



