rm(list=ls(all=TRUE))

##### 0. settings:

MutationType = 'synonymous'
AnalyzedGenes = 'ORF1ab' # 'All' # 

##### 1. read dataset with mutations and prepare a data-frame ObservedMutations with four columns: ("Pos","RefNuc","AltNuc","ObservedMutations")
#### only A T G C should be in RefNuc and AltNuc
df = read.csv('../../Sars_Cov_2/data/norm_data_modernized_nextstrain.csv')
names(df)
DfShort = df[colnames(df) %in% c("Pos","RefNuc","AltNuc")]
DfShort$MutNumber = 1
ObservedMutations = aggregate(DfShort$MutNumber, by = list(DfShort$Pos,DfShort$RefNuc,DfShort$AltNuc),FUN = sum)
names(ObservedMutations) = c("Pos","RefNuc","AltNuc","ObservedMutations")
table(ObservedMutations$RefNuc)
table(ObservedMutations$AltNuc)
ObservedMutations = ObservedMutations[ObservedMutations$RefNuc != '-' & ObservedMutations$AltNuc != '-',]

##### 2. read genome annotation
ann = read.csv("../../Sars_Cov_2/data_obtained/ideal_table.csv")
names(ann)
table(ann$GenName)

##### 3. merge ObservedMutations with ann keeping all ann
AnnMut = merge(ann,ObservedMutations, by = c("Pos","RefNuc","AltNuc"), all.x = TRUE)
AnnMut[is.na(AnnMut$ObservedMutations),] <-0
summary(AnnMut$ObservedMutations)

##### 4. observed, but not expected - how many, who they are and what to do with them?
##### these variatns with RefNuc different from the Ref - they are not so numerous (1.6%) forget about them for now
ObservedMutationsId = paste(ObservedMutations$Pos,ObservedMutations$RefNuc,ObservedMutations$AltNuc, sep = '_')
ExpectedMutationsId = paste(ann$Pos,ann$RefNuc,ann$AltNuc, sep = '_')
ObservedButNotExpected = setdiff(ObservedMutationsId,ExpectedMutationsId)
length(ObservedButNotExpected)/length(ObservedMutationsId)
# an example
head(ObservedButNotExpected)
AnnMut[AnnMut$Pos == 6851,]

##### 5. Filter AnnMut according to the goal (syn/nons, genes)

if (MutationType == 'synonymous')
{ AnnMut = AnnMut[AnnMut$GenType == 'translated' & AnnMut$RefAa != '*' & AnnMut$RefAa == AnnMut$AltAa,] }

if (AnalyzedGenes != 'All') 
{ AnnMut = AnnMut[AnnMut$GenName == AnalyzedGenes,] }

##### 6. Estimate MutSpec: Expected, Observed and 
AnnMut$MutExp <- 1
AnnMut$NucSubst = paste(AnnMut$RefNuc,AnnMut$AltNuc,sep='>')
Exp12Comp = aggregate(AnnMut$MutExp, by = list(AnnMut$NucSubst), FUN = sum); 
names(Exp12Comp)=c('NucSubst','ExpFr')

Obs12Comp = aggregate(AnnMut$ObservedMutations, by = list(AnnMut$NucSubst), FUN = sum); 
names(Obs12Comp)=c('NucSubst','ObsFr')

MutSpec12Comp = merge(Exp12Comp,Obs12Comp, by = 'NucSubst')
MutSpec12Comp$ObsToExp = MutSpec12Comp$ObsFr/MutSpec12Comp$ExpFr
MutSpec12Comp$MutSpec = MutSpec12Comp$ObsToExp/sum(MutSpec12Comp$ObsToExp)
sum(MutSpec12Comp$MutSpec) # 1 == GOOD
MutSpec12Comp

barplot(MutSpec12Comp$MutSpec, names =  MutSpec12Comp$NucSubst, main = 'mut spec 12 comp syn mut', cex.names = 0.7)

# knit: Ctrl + Shift + K => html

