library(ggplot2)
library(lubridate)
library(dplyr)
library(plotly)
rm(list=ls(all=TRUE))
#df = read.csv("../../Sars_Cov_2/data_obtained/ideal_table.csv")
df = read.csv('D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/data/norm_data_modernized_nextstrain.csv')
#nstrain = read.csv("../../Sars_Cov_2/data/norm_data_modernized_nextstrain.csv")
ann = read.csv("D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/data_obtained/ideal_table.csv")
nstrain = subset(df, RefNuc != '-' & AltNuc != '-')

for (i in 1:nrow(nstrain)){
  temporary_subset = subset(ann, Pos == nstrain$Pos[i] & AltNuc == nstrain$AltNuc[i])
  if (length(temporary_subset$RefNuc) != 0){
    nstrain$GenName[i] = temporary_subset$GenName
    nstrain$GenType[i] = temporary_subset$GenType
    nstrain$AaSub[i] = temporary_subset$AaSub
    nstrain$NeighL[i] = temporary_subset$NeighL
    nstrain$NeighR[i] = temporary_subset$NeighR
  }
}

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
    df = df[df$AaSub == 'S',]
    ann = ann[ann$AaSub == 'S',]
  }else if(mut_type == 'NS'){
    df = df[df$AaSub == 'NS',]
    ann = ann[ann$AaSub == 'NS',]
  }
  if (spec_type == 'general'){
    # Draw MutSpec 12 components for all genome 
    if(n_mutspec == 12){
      mut_list = c('A>C','A>G','A>T','C>A','C>G','C>T','G>A','G>C','G>T','T>A','T>C','T>G')
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
      MutSpec12Comp$ObsToExp = MutSpec12Comp$ObsFr/MutSpec12Comp$ExpFr
      MutSpec12Comp$ObsToExp[is.na(MutSpec12Comp$ObsToExp)] <- 0
      MutSpec12Comp$MutSpec = MutSpec12Comp$ObsToExp/sum(MutSpec12Comp$ObsToExp)
      write.csv(x = MutSpec12Comp, 'D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/R_function/data/07.MutSpec12_ForFullGenome.csv')
      
      png('D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/R_function/figures/07.MutSpec12_ForFullGenome.png')
      barplot(MutSpec12Comp$MutSpec, names =  MutSpec12Comp$NucSubst, main = 'MutSpec 12 components for full genome', cex.names = 0.7)
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
      MutSpec192Comp[is.na(MutSpec192Comp)] <- 0
      MutSpec192Comp$ObsToExp = MutSpec192Comp$ObsFr/MutSpec192Comp$ExpFr
      MutSpec192Comp$ObsToExp[is.na(MutSpec192Comp$ObsToExp)] <- 0
      MutSpec192Comp$MutSpec = MutSpec192Comp$ObsToExp/sum(MutSpec192Comp$ObsToExp)
      write.csv(x = MutSpec192Comp, 'D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/R_function/data/07.MutSpec192_ForFullGenome.csv')
      
      fig <- ggplot(MutSpec192Comp,aes(x = reorder(NucSubst, Colour, FUN = unique, order = is.ordered(NucSubst)), y = MutSpec, ymin = 0, fill = Colour))+
        geom_bar(stat = "identity", width = 0.6)+ 
        #coord_flip() + 
        xlab("Mutations") +
        ylab("Share Of Mutations")+
        theme(legend.position="none")+
        theme(axis.text.x = element_text(size=3.2, angle = 90))+
        ggtitle('Covid Mut spec 192 components for full genome')+
        #theme(axis.text.y = element_text(hjust = 3.4))+
        geom_text(aes(label=NucSubst), size = 1, angle = 90)
      fig %>% ggplotly
      ggsave('D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/R_function/figures/07.MutSpec192_ForFullGenome.svg', plot = fig, device = 'svg',  limitsize = F)
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
        write.csv(x = MutSpec12Comp, sprintf('D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/R_function/data/07.MutSpec12_%s.csv', gen))
        
        png(sprintf('D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/R_function/figures/07.MutSpec12_%s.png', gen))
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
        write.csv(x = MutSpec192Comp, sprintf('D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/R_function/data/07.MutSpec192_%s.csv', gen))
        
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
        ggsave(sprintf('D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/R_function/figures/07.MutSpec192_%s.svg', gen), plot = fig, device = 'svg',  limitsize = F)
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
        write.csv(x = MutSpec12Comp, sprintf('D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/R_function/data/07.Date_MutSpec12_%s.csv', srez))
        
        png(sprintf('D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/R_function/figures/07.Date_MutSpec12_%s.png', srez))
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
        write.csv(x = MutSpec192Comp, sprintf('D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/R_function/data/07.Date_MutSpec192_%s.csv', srez))
        
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
        ggsave(sprintf('D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/R_function/figures/07.Date_MutSpec192_%s.svg', srez), plot = fig, device = 'svg',  limitsize = F)
      }
    }
  }
}

#This function takes at the entres two data frames, df - data frame with information from nextstrain, ann - table with mutations, neighbours, aminoacids and etc.
# Possible variable values for drawing n_mutspec = 12/192, mut_type = 'S'/'NS', spec_type = 'general'/'genes'/'data', region_type = 'translated'/'untranslated'
draw_mutspec(df = nstrain, ann = ann, n_mutspec = 192, spec_type = 'data')
