rm(list=ls(all=TRUE))

ann = read.csv("../../Sars_Cov_2/data_obtained/ideal_table.csv")
names(ann)

##### 1: modify MutExp and MutObsNC
ann$MutExp <-1
ann[is.na(ann$MutObsNC),]$MutObsNC <- 0
summary(ann$MutObsNC) # max is 36347.00 -> I think on NC we may have count of mutations (mutation frequency), not the count of events => check it

##### 2: filter out everything except synonymous mutations
table(ann$GenType)  # there are 42 empty!!! what is it??

ann = ann[ann$GenType == 'translated',]
ann = ann[ann$RefAa != '*',]
ann = ann[ann$RefAa == ann$AltAa,]
table(ann$AaSub) # only 'S' => GOOD

##### 3: derive 12 component MutSpec: expected and observed
ann$NucSubst = paste(ann$RefNuc,ann$AltNuc,sep='>')
Exp12Comp = aggregate(ann$MutExp, by = list(ann$NucSubst), FUN = sum); 
names(Exp12Comp)=c('NucSubst','ExpFr')
barplot(Exp12Comp$ExpFr, names =  Exp12Comp$NucSubst, main = 'mut spec expected all syn mut', cex.names = 0.7)

Obs12Comp = aggregate(ann$MutObsNC, by = list(ann$NucSubst), FUN = sum); 
names(Obs12Comp)=c('NucSubst','ObsFr')
barplot(Obs12Comp$ObsFr, names =  Obs12Comp$NucSubst, main = 'mut spec observed all syn mut', cex.names = 0.7)

MutSpec12Comp = merge(Exp12Comp,Obs12Comp, by = 'NucSubst')
MutSpec12Comp$ObsToExp = MutSpec12Comp$ObsFr/MutSpec12Comp$ExpFr
MutSpec12Comp$MutSpec = MutSpec12Comp$ObsToExp/sum(MutSpec12Comp$ObsToExp)
sum(MutSpec12Comp$MutSpec) # 1 == GOOD
MutSpec12Comp

barplot(MutSpec12Comp$MutSpec, names =  MutSpec12Comp$NucSubst, main = 'mut spec 12 comp syn mut', cex.names = 0.7)

# knit: Ctrl + Shift + K => html









