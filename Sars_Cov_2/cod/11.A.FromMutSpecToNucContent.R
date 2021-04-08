rm(list=ls(all=TRUE))

##### INITIALIZE PARAMETERS:

Final = data.frame()
GenomeLength = 10000  
SimulationLengthNumberOfGenerations = 500000

##### A: INITIALIZE GENOME

for (MutSpecProb in c('GlobalCovidSyn'))
{
### choose initial nucleotide frequencies: equal to 25% if InitGenome == 1 or random if InitGenome > 1
for (InitGenome in 1:10)
{  # InitGenome = 2
  if (InitGenome == 1) {frA = frG = frC = frT = 0.25} # frA = 0.145342343
  if (InitGenome >  1) 
    {
    frA = runif(1); frG = runif(1);
    frC = runif(1); frT = runif(1);  
    Summa = frA+frT+frG+frC
    frA = frA/Summa; frG = frG/Summa; frC = frC/Summa; frT = frT/Summa; 
    frA+frG+frC+frT # should be 1
    }
  
### make a genome
genome = sample(c(rep('A',round(frA*GenomeLength)),rep('T',round(frT*GenomeLength)),rep('G',round(frG*GenomeLength)),rep('C',round(frC*GenomeLength))))
length(genome) # should be equal GenomeLength

##### B: DEFINE MUTATIONAL SPECTRUM
# Take MutSpec from this file:
# https://github.com/VoronkaAlexandr/Sars_Cov_2_MutSpec/blob/main/Sars_Cov_2/R_function/data/07.MutSpec12_ForFullGenome.csv

if (MutSpecProb == 'GlobalCovidSyn')
{ 
VecMutSpec = read.table("../../Sars_Cov_2/R_function/data/07.MutSpec12_ForFullGenome.csv", header = TRUE, sep = ',')
MutSpec = data.frame(VecMutSpec[c(2,6)])
MutSpec$From = gsub(">(.*)","",MutSpec$NucSubst)
MutSpec$To = gsub("(.*)>","",MutSpec$NucSubst)
MutSpec$Prob = as.numeric(MutSpec$MutSpec)
sum(MutSpec$Prob)

ExpectedFrA = sum(MutSpec[MutSpec$To == 'A',]$Prob)/sum(MutSpec[MutSpec$From == 'A',]$Prob)
ExpectedFrT = sum(MutSpec[MutSpec$To == 'T',]$Prob)/sum(MutSpec[MutSpec$From == 'T',]$Prob)
ExpectedFrG = sum(MutSpec[MutSpec$To == 'G',]$Prob)/sum(MutSpec[MutSpec$From == 'G',]$Prob)
ExpectedFrC = sum(MutSpec[MutSpec$To == 'C',]$Prob)/sum(MutSpec[MutSpec$From == 'C',]$Prob)
Summa = ExpectedFrA + ExpectedFrT + ExpectedFrG + ExpectedFrC
ExpectedFrA = ExpectedFrA/Summa
ExpectedFrT = ExpectedFrT/Summa
ExpectedFrG = ExpectedFrG/Summa
ExpectedFrC = ExpectedFrC/Summa
MutSpecProb
ExpectedFrA # 0.09887561
ExpectedFrT # 0.6645893
ExpectedFrG # 0.21168
ExpectedFrC # 0.02485505

sum(MutSpec$Prob)
for (i in 1:nrow(MutSpec))
{ # i = 1
if (i == 1) {MutSpec$RulletFrom[i] = 0}
MutSpec$RulletTo[i] = sum(MutSpec[seq(1:i),]$Prob)
if (i  > 1) MutSpec$RulletFrom[i] = MutSpec$RulletTo[i-1]
}

##### C: MUTATE AND SAVE NUCLEOTIDE CONTENT EVERY 100 GENERATIONS
for (gener in 1:SimulationLengthNumberOfGenerations)
{
### 3: choose a random position in a genome
 RandomPos = sample(1:length(genome), 1)  
 NucInRandomPos = genome[RandomPos];
 
### 2: choose a random mutation
 Rullet = runif(1)
 PotentialSubstitution = MutSpec[MutSpec$RulletFrom <= Rullet & MutSpec$RulletTo > Rullet,]
 
### 3: mutation happens
 if (nrow(PotentialSubstitution) == 1)
 {
  if (NucInRandomPos == PotentialSubstitution$From)  {genome[RandomPos] = PotentialSubstitution$To}
 }

### 4: print out every 100 generations
if ((gener %% 1000) == 0)
  {
  Res = data.frame(table(genome))
  Res = data.frame(t(Res[order(Res$genome),]))
  Res = Res[2,]
  Res$Gener = gener
  Res$InitGenome = InitGenome
  Res$MutSpecProb = MutSpecProb
  names(Res)= c('A','C','G','T','Gener','InitGenome','MutSpecProb')
  Final = rbind(Final,Res)
  }
}}}}

### 5: DERIVE FRACTIONS AND SAVE

Final$FrA = as.numeric(Final$A)/length(genome)
Final$FrT = as.numeric(Final$T)/length(genome)
Final$FrG = as.numeric(Final$G)/length(genome)
Final$FrC = as.numeric(Final$C)/length(genome)

write.table(Final, "../../Sars_Cov_2/R_function/data/11.A.FromMutSpecToNucContent.R.FinalTable.txt")

