library(ggplot2)
library(lubridate)
library(dplyr)
library(plotly)
library(stringr)
library(seqinr)
library(reshape2)
'%not_in%' <- Negate('%in%')
rm(list=ls(all=TRUE))

#reading data from GISAID and rename columns for function working


mut_df = read.csv("D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/new_data/data/U_mutations.csv")

#reading reference data (Wuhan genome NC_045512.2) for normalization
ideal_table = read.csv("D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/new_data/data/U_ideal_table.csv")

#code - working with nextstrain information

# ann_gens = ann[, c(3,5,6,10,15,16)]
# 
# df = merge(x = df, y = ann_gens, by = "Pos", all.x = TRUE)
# 
# nstrain = subset(df, RefNuc != '-' & AltNuc != '-')
# 
# nstrain$NeighL = ''
# nstrain$NeighR = ''
# 
# nstrain[nstrain$Pos != 1 & nstrain$Pos != 2 & nstrain$Pos != 29902 & nstrain$Pos != 29903,]$NeighL = toupper(str_split_fixed(nstrain[nstrain$Pos != 1 & nstrain$Pos != 2 & nstrain$Pos != 29902 & nstrain$Pos != 29903,]$parent_nucl_context, '', 5)[,2])
# nstrain[nstrain$Pos != 1 & nstrain$Pos != 2 & nstrain$Pos != 29902 & nstrain$Pos != 29903,]$NeighR = toupper(str_split_fixed(nstrain[nstrain$Pos != 1 & nstrain$Pos != 2 & nstrain$Pos != 29902 & nstrain$Pos != 29903,]$parent_nucl_context, '', 5)[,4])
# 
# nstrain = subset(nstrain, NeighL != '-' & NeighR != '-')
# nstrain$CodonRef = ''
# nstrain$CodonAlt = ''
# 
# nstrain[nstrain$NucInCodon == 1,]$CodonRef = toupper(paste(str_split_fixed(nstrain[nstrain$NucInCodon == 1,]$parent_nucl_context, '', 5)[,3],str_split_fixed(nstrain[nstrain$NucInCodon == 1,]$parent_nucl_context, '', 5)[,4],str_split_fixed(nstrain[nstrain$NucInCodon == 1,]$parent_nucl_context, '', 5)[,5],sep=''))
# nstrain[nstrain$NucInCodon == 1,]$CodonAlt = toupper(paste(str_split_fixed(nstrain[nstrain$NucInCodon == 1,]$child_nucl_context, '', 5)[,3],str_split_fixed(nstrain[nstrain$NucInCodon == 1,]$parent_nucl_context, '', 5)[,4],str_split_fixed(nstrain[nstrain$NucInCodon == 1,]$parent_nucl_context, '', 5)[,5],sep=''))
# 
# nstrain[nstrain$NucInCodon == 2,]$CodonRef = toupper(paste(str_split_fixed(nstrain[nstrain$NucInCodon == 2,]$parent_nucl_context, '', 5)[,2],str_split_fixed(nstrain[nstrain$NucInCodon == 2,]$parent_nucl_context, '', 5)[,3],str_split_fixed(nstrain[nstrain$NucInCodon == 2,]$parent_nucl_context, '', 5)[,4],sep=''))
# nstrain[nstrain$NucInCodon == 2,]$CodonAlt = toupper(paste(str_split_fixed(nstrain[nstrain$NucInCodon == 2,]$child_nucl_context, '', 5)[,2],str_split_fixed(nstrain[nstrain$NucInCodon == 2,]$parent_nucl_context, '', 5)[,3],str_split_fixed(nstrain[nstrain$NucInCodon == 2,]$parent_nucl_context, '', 5)[,4],sep=''))
# 
# nstrain[nstrain$NucInCodon == 3,]$CodonRef = toupper(paste(str_split_fixed(nstrain[nstrain$NucInCodon == 3,]$parent_nucl_context, '', 5)[,1],str_split_fixed(nstrain[nstrain$NucInCodon == 3,]$parent_nucl_context, '', 5)[,2],str_split_fixed(nstrain[nstrain$NucInCodon == 3,]$parent_nucl_context, '', 5)[,3],sep=''))
# nstrain[nstrain$NucInCodon == 3,]$CodonAlt = toupper(paste(str_split_fixed(nstrain[nstrain$NucInCodon == 3,]$child_nucl_context, '', 5)[,1],str_split_fixed(nstrain[nstrain$NucInCodon == 3,]$parent_nucl_context, '', 5)[,2],str_split_fixed(nstrain[nstrain$NucInCodon == 3,]$parent_nucl_context, '', 5)[,3],sep=''))
# 
# nstrain = subset(nstrain, '-' %not_in% CodonRef & '-' %not_in% CodonAlt)
# 
# nstrain$AaRef = ''
# nstrain$AaAlt = ''
# ref_string = toString(nstrain[nstrain$GenType == 'translated',]$CodonRef)
# ref_string = str_replace_all(ref_string, ', ', '')
# ref_string = s2c(ref_string)
# nstrain[nstrain$GenType == 'translated',]$AaRef = translate(seq = ref_string)
# 
# alt_string = toString(nstrain[nstrain$GenType == 'translated',]$CodonAlt)
# alt_string = str_replace_all(alt_string, ', ', '')
# alt_string = s2c(alt_string)
# nstrain[nstrain$GenType == 'translated',]$AaAlt = translate(seq = alt_string)
# 
# nstrain$AaSub = ''
# nstrain[nstrain$GenType == 'translated',]$AaSub = ifelse(nstrain[nstrain$GenType == 'translated',]$AaRef == nstrain[nstrain$GenType == 'translated',]$AaAlt, 'S', 'NS')
# 
# write.csv(x = nstrain, 'D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/new_data/data/uppdate_mutations.csv')



#Function for drawing Mutation spectrums
draw_mutspec = function(df, ann, n_mutspec = 12, mut_type = 'S', spec_type = 'general', region_type = 'translated'){
  ann = as.data.frame(ann)
  df = as.data.frame(df)
  if (region_type == 'translated'){
    df = df[df$GenType == 'translated',]
    ann = ann[ann$GenType == 'translated',]
  }else if(region_type == 'untranslated'){
    df = df[df$GenType == 'untranslated',]
    ann = ann[ann$GenType == 'untranslated',]
  }
  if (mut_type == 'S'){
    df = df[df$AaSub == 'S',]  #synonyms - without changing AA
    ann = ann[ann$AaSub == 'S',]
  }else if(mut_type == 'NS'){
    df = df[df$AaSub == 'NS',]
    ann = ann[ann$AaSub == 'NS',]
  }
  if (spec_type == 'general'){
    # Draw MutSpec 12 components for all genome 
    if(n_mutspec == 12){
      mut_list = c('A>C','A>G','A>U','C>A','C>G','C>U','G>A','G>C','G>U','U>A','U>C','U>G')
      mut_list = data.frame(mut_list)
      names(mut_list)=c('NucSubst')
      ann$MutExp <- 1
      ann$NucSubst = paste(ann$RefNuc,ann$AltNuc,sep='>')
      Exp12Comp = aggregate(ann$MutExp, by = list(ann$NucSubst), FUN = sum);
      names(Exp12Comp) = c('NucSubst','ExpFr')
      
      df$MutObs = 1
      df$NucSubst = paste(df$RefNuc,df$AltNuc,sep='>')
      Obs12Comp = aggregate(df$MutObs, by = list(df$NucSubst), FUN = sum);
      names(Obs12Comp) = c('NucSubst','ObsFr')
      
      MutSpec12Comp = merge(mut_list,Exp12Comp, by = 'NucSubst', all.x = TRUE)
      MutSpec12Comp = merge(MutSpec12Comp,Obs12Comp, by = 'NucSubst',all.x = TRUE)
      MutSpec12Comp[is.na(MutSpec12Comp)] <- 0
      MutSpec12Comp$MutSpec = MutSpec12Comp$ObsFr/MutSpec12Comp$ExpFr  #Normalization on all possible same mutations in RefSeq
      MutSpec12Comp$MutSpec[is.na(MutSpec12Comp$MutSpec)] <- 0
      #MutSpec12Comp$MutSpec = MutSpec12Comp$ObsToExp/sum(MutSpec12Comp$ObsToExp)
      write.csv(x = MutSpec12Comp, 'D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/new_data/data_obtained/07.S_MutSpec12_ForFullGenome.csv')
      
      png('D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/new_data/figures/07.MutSpec12_ForFullGenome_S.png')
      barplot(MutSpec12Comp$MutSpec, names =  MutSpec12Comp$NucSubst, main = 'MutSpec 12 components for full genome S', cex.names = 0.7)
      dev.off()
    }else if(n_mutspec == 192){
      #Draw MutSpec 192 components for all Genome
      ann$MutExp = 1
      ann$FromWithNeigh = paste(ann$NeighL,ann$RefNuc,ann$NeighR, sep='')
      ann$ToWithNeigh = paste(ann$NeighL,ann$AltNuc,ann$NeighR, sep='')
      ann$MutSpecWithNeigh = paste(ann$FromWithNeigh, ann$ToWithNeigh, sep = '>')
      Exp192Comp = aggregate(ann$MutExp, by = list(ann$MutSpecWithNeigh), FUN = sum);
      names(Exp192Comp) = c('NucSubst','ExpFr')
      ann$Colour = paste(ann$RefNuc, ann$AltNuc, sep = '>')
      Exp_with_colour = aggregate(ann$Colour, by = list(ann$MutSpecWithNeigh), FUN=unique)
      names(Exp_with_colour) = c('NucSubst','Colour')
      Exp192Comp = merge(Exp_with_colour, Exp192Comp, by = 'NucSubst', all.x = TRUE)
      Exp192Comp = Exp192Comp[order(Exp192Comp$Colour),]
      rownames(Exp192Comp) <- 1:nrow(Exp192Comp)
      Exp192Comp$Colour = factor(Exp192Comp$Colour)
      
      df$MutObs = 1
      df$FromWithNeigh = paste(df$NeighL,df$RefNuc,df$NeighR, sep='')
      df$ToWithNeigh = paste(df$NeighL,df$AltNuc,df$NeighR, sep='')
      df$MutSpecWithNeigh = paste(df$FromWithNeigh, df$ToWithNeigh, sep = '>')
      Obs192Comp = aggregate(df$MutObs, by = list(df$MutSpecWithNeigh), FUN = sum);
      names(Obs192Comp) = c('NucSubst','ObsFr')
      
      MutSpec192Comp = merge(Obs192Comp,Exp192Comp, by = 'NucSubst', all.x = TRUE)
      
      MutSpec192Comp$Colour = paste(str_split_fixed(MutSpec192Comp$NucSubst, '', 7)[,2],str_split_fixed(MutSpec192Comp$NucSubst, '', 7)[,6],sep='>')
      MutSpec192Comp[is.na(MutSpec192Comp)] <- 0
      
      MutSpec192Comp$ObsToExp = MutSpec192Comp$ObsFr/MutSpec192Comp$ExpFr
      MutSpec192Comp$ObsToExp[is.na(MutSpec192Comp$ObsToExp)] <- 0
      MutSpec192Comp$MutSpec = MutSpec192Comp$ObsToExp/sum(MutSpec192Comp$ObsToExp)
      MutSpec192Comp$Context = paste(str_split_fixed(MutSpec192Comp$NucSubst, '', 7)[,1], 'x', str_split_fixed(MutSpec192Comp$NucSubst, '', 7)[,3],sep='')
      write.csv(x = MutSpec192Comp, 'D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/new_data/data_obtained/07.MutSpec192_ForFullGenome.csv')
      
      fig <- ggplot(MutSpec192Comp,aes(x = reorder(NucSubst, Colour, FUN = unique, order = is.ordered(NucSubst)), y = MutSpec, ymin = 0, fill = Colour))+
        geom_bar(stat = "identity", width = 0.6)+ 
        #coord_flip() + 
        xlab("Mutations") +
        ylab("Share Of Mutations")+
        theme(legend.position="none")+
        theme(axis.text.x = element_text(size=3.2, angle = 90))+
        ggtitle('Covid Mut spec 192 components for full genome')+
        #theme(axis.text.y = element_text(hjust = 3.4))+
        geom_text(aes(label=NucSubst), size = 1.4, angle = 90)
      fig %>% ggplotly
      ggsave('D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/new_data/figures/07.MutSpec192_ForFullGenome.svg', plot = fig, device = 'svg',  limitsize = F)
    }
  }else if(spec_type == 'genes'){
    #Draw MutSpec 12 components per Gen
    if(n_mutspec == 12){
      for (gen in unique(df$GenName)){
        gen_df = df[df$GenName == gen,]
        gen_ann = ann[ann$GenName == gen,]
        mut_list = c('A>C','A>G','A>T','C>A','C>G','C>T','G>A','G>C','G>T','T>A','T>C','T>G')
        mut_list = data.frame(mut_list)
        names(mut_list)=c('NucSubst')
        gen_ann$MutExp = 1
        gen_ann$NucSubst = paste(gen_ann$RefNuc,gen_ann$AltNuc,sep='>')
        Exp12Comp = aggregate(gen_ann$MutExp, by = list(gen_ann$NucSubst), FUN = sum);
        names(Exp12Comp) = c('NucSubst','ExpFr')
        
        gen_df$MutObs = 1
        gen_df$NucSubst = paste(gen_df$RefNuc,gen_df$AltNuc,sep='>')
        Obs12Comp = aggregate(gen_df$MutObs, by = list(gen_df$NucSubst), FUN = sum);
        names(Obs12Comp) = c('NucSubst','ObsFr')
        
        MutSpec12Comp = merge(mut_list,Exp12Comp, by = 'NucSubst', all.x = TRUE)
        MutSpec12Comp = merge(MutSpec12Comp,Obs12Comp, by = 'NucSubst',all.x = TRUE)
        MutSpec12Comp[is.na(MutSpec12Comp)] <- 0
        MutSpec12Comp$ObsToExp = MutSpec12Comp$ObsFr/MutSpec12Comp$ExpFr
        MutSpec12Comp$ObsToExp[is.na(MutSpec12Comp$ObsToExp)] <- 0
        MutSpec12Comp$MutSpec = MutSpec12Comp$ObsToExp/sum(MutSpec12Comp$ObsToExp)
        write.csv(x = MutSpec12Comp, sprintf('D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/new_data/data_obtained/07.MutSpec12_%s.csv', gen))
        
        png(sprintf('D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/new_data/figures/07.MutSpec12_%s.png', gen))
        barplot(MutSpec12Comp$MutSpec, names =  MutSpec12Comp$NucSubst, main = sprintf('MutSpec 12 components per Gen %s', gen), cex.names = 0.7)
        dev.off()
      }
    }else if(n_mutspec == 192){
      #Draw MutSpec 192 components per Gen
      for (gen in unique(df$GenName)){
        gen_df = df[df$GenName == gen,]
        gen_ann = ann[ann$GenName == gen,]
        gen_ann$MutExp = 1
        gen_ann$FromWithNeigh = paste(gen_ann$NeighL,gen_ann$RefNuc,gen_ann$NeighR, sep='')
        gen_ann$ToWithNeigh = paste(gen_ann$NeighL,gen_ann$AltNuc,gen_ann$NeighR, sep='')
        gen_ann$MutSpecWithNeigh = paste(gen_ann$FromWithNeigh, gen_ann$ToWithNeigh, sep = '>')
        Exp192Comp = aggregate(gen_ann$MutExp, by = list(gen_ann$MutSpecWithNeigh), FUN = sum);
        names(Exp192Comp) = c('NucSubst','ExpFr')
        gen_ann$Colour = paste(gen_ann$RefNuc, gen_ann$AltNuc, sep = '>')
        Exp_with_colour = aggregate(gen_ann$Colour, by = list(gen_ann$MutSpecWithNeigh), FUN=unique)
        names(Exp_with_colour) = c('NucSubst','Colour')
        Exp192Comp = merge(Exp_with_colour, Exp192Comp, by = 'NucSubst', all.x = TRUE)
        Exp192Comp = Exp192Comp[order(Exp192Comp$Colour),]
        rownames(Exp192Comp) <- 1:nrow(Exp192Comp)
        Exp192Comp$Colour = factor(Exp192Comp$Colour)
        
        gen_df$MutObs = 1
        gen_df$FromWithNeigh = paste(gen_df$NeighL,gen_df$RefNuc,gen_df$NeighR, sep='')
        gen_df$ToWithNeigh = paste(gen_df$NeighL,gen_df$AltNuc,gen_df$NeighR, sep='')
        gen_df$MutSpecWithNeigh = paste(gen_df$FromWithNeigh, gen_df$ToWithNeigh, sep = '>')
        Obs192Comp = aggregate(gen_df$MutObs, by = list(gen_df$MutSpecWithNeigh), FUN = sum);
        names(Obs192Comp) = c('NucSubst','ObsFr')
        
        MutSpec192Comp = merge(Obs192Comp,Exp192Comp, by = 'NucSubst', all.x = TRUE)
        MutSpec192Comp[is.na(MutSpec192Comp)] <- 0
        MutSpec192Comp$ObsToExp = MutSpec192Comp$ObsFr/MutSpec192Comp$ExpFr
        MutSpec192Comp$ObsToExp[is.na(MutSpec192Comp$ObsToExp)] <- 0
        MutSpec192Comp$MutSpec = MutSpec192Comp$ObsToExp/sum(MutSpec192Comp$ObsToExp)
        write.csv(x = MutSpec192Comp, sprintf('D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/new_data/data_obtained/07.MutSpec192_%s.csv', gen))
        
        fig <- ggplot(MutSpec192Comp,aes(x = reorder(NucSubst, Colour, FUN = unique, order = is.ordered(NucSubst)), y = MutSpec, ymin = 0, fill = Colour))+
          geom_bar(stat = "identity", width = 0.6)+ 
          #coord_flip() + 
          xlab("Mutations") +
          ylab("Share Of Mutations")+
          theme(legend.position="none")+
          theme(axis.text.x = element_text(size=3.2, angle = 90))+
          ggtitle(sprintf('Covid Mut spec 192 components for %s', gen))+
          #theme(axis.text.y = element_text(hjust = 3.4))+
          geom_text(aes(label=NucSubst), size = 1, angle = 90)
        fig %>% ggplotly
        ggsave(sprintf('D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/new_data/figures/07.MutSpec192_%s.svg', gen), plot = fig, device = 'svg',  limitsize = F)
      }
    }
  }else if(spec_type == 'data'){
    if(n_mutspec == 12){
      date_df <- df %>%
        arrange(Time)
      ann$MutExp = 1
      ann$NucSubst = paste(ann$RefNuc,ann$AltNuc,sep='>')
      Exp12Comp = aggregate(ann$MutExp, by = list(ann$NucSubst), FUN = sum);
      names(Exp12Comp) = c('NucSubst','ExpFr')
      for (i in seq(1, length(as.Date(date_df$Time)), 500)){
        if ((i+500) <= length(date_df$Time)){
          time_set = date_df[i:(i+500),]
          srez = paste(date_df$Time[i], date_df$Time[i+500], sep = ' to ')
        }
        else{
          time_set = date_df[i:length(date_df$Time),]
          srez = paste(date_df$Time[i], date_df$Time[length(date_df$Time)], sep = ' to ')
        }
        mut_list = c('A>C','A>G','A>T','C>A','C>G','C>T','G>A','G>C','G>T','T>A','T>C','T>G')
        mut_list = data.frame(mut_list)
        names(mut_list)=c('NucSubst')
        
        time_set$MutObs = 1
        time_set$NucSubst = paste(time_set$RefNuc,time_set$AltNuc,sep='>')
        Obs12Comp = aggregate(time_set$MutObs, by = list(time_set$NucSubst), FUN = sum);
        names(Obs12Comp) = c('NucSubst','ObsFr')
        
        MutSpec12Comp = merge(mut_list,Exp12Comp, by = 'NucSubst', all.x = TRUE)
        MutSpec12Comp = merge(MutSpec12Comp,Obs12Comp, by = 'NucSubst',all.x = TRUE)
        MutSpec12Comp[is.na(MutSpec12Comp)] <- 0
        MutSpec12Comp$ObsToExp = MutSpec12Comp$ObsFr/MutSpec12Comp$ExpFr
        MutSpec12Comp$ObsToExp[is.na(MutSpec12Comp$ObsToExp)] <- 0
        MutSpec12Comp$MutSpec = MutSpec12Comp$ObsToExp/sum(MutSpec12Comp$ObsToExp)
        write.csv(x = MutSpec12Comp, sprintf('D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/new_data/data_obtained/07.Date_MutSpec12_%s.csv', srez))
        
        png(sprintf('D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/new_data/figures/07.Date_MutSpec12_%s.png', srez))
        barplot(MutSpec12Comp$MutSpec, names =  MutSpec12Comp$NucSubst, main = sprintf('MutSpec 12 components for Date %s', srez), cex.names = 0.7)
        dev.off()
      }
    }else if(n_mutspec == 192){
      date_df <- df %>%
        arrange(Time)
      ann$MutExp = 1
      ann$FromWithNeigh = paste(ann$NeighL,ann$RefNuc,ann$NeighR, sep='')
      ann$ToWithNeigh = paste(ann$NeighL,ann$AltNuc,ann$NeighR, sep='')
      ann$MutSpecWithNeigh = paste(ann$FromWithNeigh, ann$ToWithNeigh, sep = '>')
      Exp192Comp = aggregate(ann$MutExp, by = list(ann$MutSpecWithNeigh), FUN = sum);
      names(Exp192Comp) = c('NucSubst','ExpFr')
      ann$Colour = paste(ann$RefNuc, ann$AltNuc, sep = '>')
      Exp_with_colour = aggregate(ann$Colour, by = list(ann$MutSpecWithNeigh), FUN=unique)
      names(Exp_with_colour) = c('NucSubst','Colour')
      Exp192Comp = merge(Exp_with_colour, Exp192Comp, by = 'NucSubst', all.x = TRUE)
      Exp192Comp = Exp192Comp[order(Exp192Comp$Colour),]
      rownames(Exp192Comp) <- 1:nrow(Exp192Comp)
      Exp192Comp$Colour = factor(Exp192Comp$Colour)
      for (i in seq(1, length(as.Date(date_df$Time)), 500)){
        if ((i+500) <= length(date_df$Time)){
          time_set = date_df[i:(i+500),]
          srez = paste(date_df$Time[i], date_df$Time[i+500], sep = ' to ')
        }else{
          time_set = date_df[i:length(date_df$Time),]
          srez = paste(date_df$Time[i], date_df$Time[length(date_df$Time)], sep = ' to ')
        }
        time_set$MutObs = 1
        time_set$FromWithNeigh = paste(time_set$NeighL,time_set$RefNuc,time_set$NeighR, sep='')
        time_set$ToWithNeigh = paste(time_set$NeighL,time_set$AltNuc,time_set$NeighR, sep='')
        time_set$MutSpecWithNeigh = paste(time_set$FromWithNeigh, time_set$ToWithNeigh, sep = '>')
        Obs192Comp = aggregate(time_set$MutObs, by = list(time_set$MutSpecWithNeigh), FUN = sum);
        names(Obs192Comp) = c('NucSubst','ObsFr')
        
        MutSpec192Comp = merge(Obs192Comp,Exp192Comp, by = 'NucSubst', all.x = TRUE)
        MutSpec192Comp[is.na(MutSpec192Comp)] <- 0
        MutSpec192Comp$ObsToExp = MutSpec192Comp$ObsFr/MutSpec192Comp$ExpFr
        MutSpec192Comp$ObsToExp[is.na(MutSpec192Comp$ObsToExp)] <- 0
        MutSpec192Comp$MutSpec = MutSpec192Comp$ObsToExp/sum(MutSpec192Comp$ObsToExp)
        write.csv(x = MutSpec192Comp, sprintf('D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/new_data/data_obtained/07.Date_MutSpec192_%s.csv', srez))
        
        fig <- ggplot(MutSpec192Comp,aes(x = reorder(NucSubst, Colour, FUN = unique, order = is.ordered(NucSubst)), y = MutSpec, ymin = 0, fill = Colour))+
          geom_bar(stat = "identity", width = 0.6)+ 
          #coord_flip() + 
          xlab("Mutations") +
          ylab("Share Of Mutations")+
          theme(legend.position="none")+
          theme(axis.text.x = element_text(size=3.2, angle = 90))+
          ggtitle(sprintf('MutSpec 192 components for Date %s', srez))+
          #theme(axis.text.y = element_text(hjust = 3.4))+
          geom_text(aes(label=NucSubst), size = 1, angle = 90)
        fig %>% ggplotly
        ggsave(sprintf('D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/new_data/figures/07.Date_MutSpec192_%s.svg', srez), plot = fig, device = 'svg',  limitsize = F)
      }
    }
  }
}

#This function takes at the entres two data frames, df - data frame with information from nextstrain, ann - table with mutations, neighbours, aminoacids and etc.
# Possible variable values for drawing n_mutspec = 12/192, mut_type = 'S'/'NS', spec_type = 'general'/'genes'/'data', region_type = 'translated'/'untranslated'
draw_mutspec(df = mut_df, ann = ideal_table, n_mutspec = 192, spec_type = 'general', mut_type='S')




df = mut_df[mut_df$AaSub == 'NS',]
ann = ideal_table[ideal_table$AaSub == 'NS',]
###  Unique mutspec for EDGE
mut_list = c('A>C','A>G','A>U','C>A','C>G','C>U','G>A','G>C','G>U','U>A','U>C','U>G')
mut_list = data.frame(mut_list)
names(mut_list)=c('NucSubst')
ann$MutExp <- 1
ann$NucSubst = paste(ann$RefNuc,ann$AltNuc,sep='>')
Exp12Comp = aggregate(ann$MutExp, by = list(ann$NucSubst), FUN = sum);
names(Exp12Comp) = c('NucSubst','ExpFr')

Obs12Comp = mut_list
for (i in sort(unique(df$edge_level))){
  edge_set = df[df$edge_level == i,]
  edge_set$MutObs = 1
  edge_set$NucSubst = paste(edge_set$RefNuc,edge_set$AltNuc,sep='>')
  edge_Obs12Comp = aggregate(edge_set$MutObs, by = list(edge_set$NucSubst), FUN = sum);
  names(edge_Obs12Comp) = c('NucSubst',sprintf('ObsFr_edge%s',i))
  Obs12Comp = merge(Obs12Comp, edge_Obs12Comp, by = 'NucSubst', all.x = TRUE)
}

MutSpec12Comp['NucSubst']
MutSpec12Comp = merge(mut_list,Exp12Comp, by = 'NucSubst', all.x = TRUE)
MutSpec12Comp = merge(MutSpec12Comp,Obs12Comp, by = 'NucSubst',all.x = TRUE)
MutSpec12Comp[is.na(MutSpec12Comp)] <- 0
for (i in sort(unique(df$edge_level))){
  MutSpec12Comp[sprintf('ObstoExp_edge%s',i)] = MutSpec12Comp[sprintf('ObsFr_edge%s',i)]/MutSpec12Comp$ExpFr
  MutSpec12Comp[sprintf('ObstoExp_edge%s',i)][is.na(MutSpec12Comp[sprintf('ObstoExp_edge%s',i)])] <- 0
  MutSpec12Comp[sprintf('MutSpec_edge%s',i)] = MutSpec12Comp[sprintf('ObstoExp_edge%s',i)]/sum(MutSpec12Comp[sprintf('ObstoExp_edge%s',i)])
}
Names = names(MutSpec12Comp)

plot_MutSpec12Comp = melt(MutSpec12Comp[c('NucSubst',Names[grepl("Mut", Names)])], id.vars="NucSubst")
write.csv(x = MutSpec12Comp, 'D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/new_data/data_obtained/07.NS_Edge_MutSpec_12.csv')
fig = ggplot(data=plot_MutSpec12Comp, aes(x=variable, y=value, group=NucSubst, colour = NucSubst)) +
  geom_line(size = 2)+
  geom_point(color="red", size=3)+
  scale_color_brewer(palette="Paired")
ggsave('D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/new_data/figures/07.NS_Edge_MutSpec_12.svg', plot = fig, device = 'svg',  width=15, height=8)



#library(RColorBrewer)

#mycolors = c(brewer.pal(name="Dark2", n = 8), brewer.pal(name="Paired", n = 8))
#                             192 MutSpec with points
# fig <- ggplot(MutSpec192Comp,aes(x = reorder(NucSubst, Colour, FUN = unique, order = is.ordered(NucSubst)), y = MutSpec, ymin = 0, fill = Colour))+
#   geom_bar(stat = "identity", width = 0.6)+ 
#   #coord_flip() + 
#   xlab("Mutations") +
#   ylab("Share Of Mutations")+
#   theme(legend.position="none")+
#   theme(axis.text.x = element_text(size=3.2, angle = 90))+
#   ggtitle('Covid Mut spec 192 components for full genome')+
#   #theme(axis.text.y = element_text(hjust = 3.4))+
#   geom_point(aes(colour = Context), size = 1)+
#   scale_color_manual(values=mycolors)+
#   theme(legend.position="right")+
#   theme(legend.title = element_text(size=5,face="bold"))+
#   theme(legend.text = element_text(size=5,face="bold"))+
#   theme(legend.key.size = unit(0.5, 'cm'))
# fig %>% ggplotly
# ggsave('D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/new_data/figures/07.point_MutSpec192_ForFullGenome.svg', plot = fig, device = 'svg',  limitsize = F)
