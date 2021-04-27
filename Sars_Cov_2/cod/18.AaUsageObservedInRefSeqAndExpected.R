rm(list=ls(all=TRUE))

### aminoacid usage in all protein coding genes on SarsCov2:
CodonUsageRefSeq = read.table("../../Sars_Cov_2/data_obtained/Usage/13.Codon_usage_NC_045512.2.csv", sep = ',',header = TRUE)

### expected amino-acid usage at the end of Vitya evolution
data  = read.table("../../Sars_Cov_2/data_obtained/AminoSimFromVictor/generations1.csv", sep = ';',header = TRUE)

data = data[data$generations == 1000000,]
data = data[,-c(1)]
sum(data[c(1:64)]) # 6400
data = data/6400 # 

#### from CodonUsage to AaUsage
DataT = data.frame(t(data))
DataT$Codons = row.names(DataT)
DataT$Aa = DataT$Codons

TranslateStandardCodonsIntoThreeLetterAa<-function(x)
{
  if (x %in% c('TTT','TTC')) {return ("Phe")}
  if (x %in% c('TTA','TTG')) {return ("LeuTT")}
  if (x %in% c('CTT','CTC','CTA','CTG')) {return ("LeuCT")}
  if (x %in% c('ATT','ATC','ATA')) {return ("Ile")} # !!!
  if (x %in% c('ATG')) {return ("Met")} # !!!
  if (x %in% c('GTC','GTA','GTG','GTT')) {return ("Val")}
  
  if (x %in% c('TCT','TCC','TCA','TCG')) {return ("SerTC")}
  if (x %in% c('CCT','CCC','CCA','CCG')) {return ("Pro")}
  if (x %in% c('ACT','ACC','ACA','ACG')) {return ("Thr")}
  if (x %in% c('GCT','GCC','GCA','GCG')) {return ("Ala")}
  
  if (x %in% c('TAT','TAC')) {return ("Tyr")}
  if (x %in% c('TAA','TAG','TGA')) {return ("Stop")} # !!!
  if (x %in% c('CAT','CAC')) {return ("His")}
  if (x %in% c('CAA','CAG')) {return ("Gln")}
  if (x %in% c('AAT','AAC')) {return ("Asn")}
  if (x %in% c('AAA','AAG')) {return ("Lys")}
  if (x %in% c('GAT','GAC')) {return ("Asp")}
  if (x %in% c('GAA','GAG')) {return ("Glu")}
  
  if (x %in% c('TGT','TGC')) {return ("Cys")}
  if (x %in% c('TGG')) {return ("Trp")} # !!!!
  if (x %in% c('CGT','CGC','CGA','CGG')) {return ("Arg")}
  if (x %in% c('AGT','AGC')) {return ("SerAG")}
  if (x %in% c('AGA','AGG')) {return ("Arg")}  # !!! 
  if (x %in% c('GGT','GGC','GGA','GGG')) {return ("Gly")}
}

DataT$Aa = apply(as.matrix(DataT$Codons),1,FUN = TranslateStandardCodonsIntoThreeLetterAa)
names(DataT)[1]=c('ExpCodonFreq')
agg = aggregate(DataT$ExpCodonFreq, by = list(DataT$Aa), FUN = sum)
names(agg)=c('Aa','Freq')

#### 
ExpAndObsCodonUsage = merge(CodonUsageRefSeq, DataT, by = 'Codons')
cor.test(ExpAndObsCodonUsage$ExpCodonFreq,ExpAndObsCodonUsage$normed_number, method= 'spearman')
plot(log2(ExpAndObsCodonUsage$ExpCodonFreq),log2(ExpAndObsCodonUsage$normed_number))

#### 
Agg = aggregate(list(ExpAndObsCodonUsage$ExpCodonFreq,ExpAndObsCodonUsage$normed_number),by=list(ExpAndObsCodonUsage$Aa), FUN = sum)
names(Agg)=c('Aa','ExpCodonFreq','ObsCodonFreq')
cor.test(Agg$ExpCodonFreq,Agg$ObsCodonFreq, method = 'spearman')
plot(Agg$ExpCodonFreq,Agg$ObsCodonFreq)


