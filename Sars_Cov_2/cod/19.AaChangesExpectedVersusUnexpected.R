rm(list=ls(all=TRUE))

FromTo = read.table("../../Sars_Cov_2/R_function/gisaid_mutspec/data/17.Speed_aa_mut.csv", sep = ',',header = TRUE)
FromTo$UniqueIDForAaPair = 'NA'
for (i in 1:nrow(FromTo))
{ # i = 1
FromTo$UniqueIDForAaPair[i] = paste(sort(c(unlist(strsplit(FromTo$Aa_Mutation[i],'>'))[1],unlist(strsplit(FromTo$Aa_Mutation[i],'>'))[2])), collapse = '>')
}

Expected = FromTo[FromTo$Group == 1,]
for (i in 1:ncol(Expected))
{
  names(Expected)[i] = paste('Exp',names(Expected)[i],sep='.')  
}
unique(Expected$Exp.UniqueIDForAaPair)
VecOfExp.UniqueIDForAaPair.ShapedByCToT = c("A>V","H>Y","F>L_CT","L_CT>P","P>S_TC","F>S_TC","L_TT>S_TC","I>T","M>T")

Unexpected = FromTo[FromTo$Group == 0 & FromTo$UniqueIDForAaPair %in% Expected$Exp.UniqueIDForAaPair,]
for (i in 1:ncol(Unexpected))
{
  names(Unexpected)[i] = paste('Unexp',names(Unexpected)[i],sep='.')  
}

Compara = merge(Expected,Unexpected, by.x = 'Exp.UniqueIDForAaPair',  by.y = 'Unexp.UniqueIDForAaPair', all = TRUE)
Compara$RatioOfRatesExpToUnexp = Compara$Exp.Speed_Mut/Compara$Unexp.Speed_Mut
Compara = Compara[order(Compara$RatioOfRatesExpToUnexp),]
summary(Compara$RatioOfRatesExpToUnexp)
wilcox.test(Compara$RatioOfRatesExpToUnexp, mu = 1)

ComparaNoStopCodons = Compara[!grepl("\\*",Compara$Exp.UniqueIDForAaPair),]
summary(ComparaNoStopCodons$RatioOfRatesExpToUnexp)

### only due to C>T
summary(Compara[Compara$Exp.UniqueIDForAaPair %in% VecOfExp.UniqueIDForAaPair.ShapedByCToT,]$RatioOfRatesExpToUnexp)
wilcox.test(Compara[Compara$Exp.UniqueIDForAaPair %in% VecOfExp.UniqueIDForAaPair.ShapedByCToT,]$RatioOfRatesExpToUnexp,mu=1)

summary(ComparaNoStopCodons[ComparaNoStopCodons$Exp.UniqueIDForAaPair %in% VecOfExp.UniqueIDForAaPair.ShapedByCToT,]$RatioOfRatesExpToUnexp)


wilcox.test(ComparaNoStopCodons$RatioOfRatesExpToUnexp, mu = 1)

##### only due to G>T NoStopCodons
summary(ComparaNoStopCodons[!ComparaNoStopCodons$Exp.UniqueIDForAaPair %in% VecOfExp.UniqueIDForAaPair.ShapedByCToT,]$RatioOfRatesExpToUnexp) #  !!!!! stronger effect!!
wilcox.test(ComparaNoStopCodons[!ComparaNoStopCodons$Exp.UniqueIDForAaPair %in% VecOfExp.UniqueIDForAaPair.ShapedByCToT,]$RatioOfRatesExpToUnexp, mu = 1)

##### Try fisher test for some of them and try regions of proteins less constrained



# from one letter to three letter code AA:
A = c("G","A","L","M","F","W","K","Q","E","S","P","V","I","C","Y","H","R","N","D","T","*")
AAA = c("Gly","Ala","Leu","Met","Phe","Trp","Lys","Gln","Glu","Ser","Pro","Val","Ile","Cys","Tyr","His","Arg","Asn","Asp","Thr","Stop")
AA = data.frame(A,AAA)
