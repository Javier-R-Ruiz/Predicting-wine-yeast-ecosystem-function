##PHYLOGENETIC_FUNCTION_ANALYSIS

rm(list=ls())
library(data.table)
library(msa)
library(phangorn)
library(ggtree)
library(phytools)
library(Biostrings)
library(ape)
library(ggplot2)
library(operators)
library(ggpubr)
library(Rphylopars)
library(scales)
library(DECIPHER)

#Build_Phylogenetic_Tree

fasta = readDNAStringSet('../Data/26S_sequences_113_strains.txt')
mult <- msa(fasta, method="ClustalW", type="dna", order="input")
writeXStringSet(as(unmasked(mult), "XStringSet"), file='../Data/Temporary/26s.fasta')
write.phyDat(phy_mult,file='../Data/Temporary/26s.phy')
system('trimal -in ../Data/Temporary/26s.phy -out ../Data/Temporary/Trimmed_26s.phy -automated1 ')
system('iqtree -s  ../Data/Temporary/Trimmed_26s.phy -bb 1000 -redo')

rm(list=ls())

pdata = fread('../Data/Phenotypes_60.csv')
Phylogenetic_Tree = read.newick('../Data/Temporary/Trimmed_26s.phy.contree')

#Root to basidiomycota
Phylogenetic_Tree <- root(Phylogenetic_Tree,outgroup=c('Rt_1','Rm_1','Ng_1','Ca_1','Ca_2','Ca_3','Ca_4','Ca_5','Ca_6'),resolve.root = TRUE)
Phylogenetic_Tree113 = chronos(Phylogenetic_Tree)
write.tree(Phylogenetic_Tree113,'../Data/Ultrametric_Tree_113.phy')
Phylogenetic_Tree60 =keep.tip(Phylogenetic_Tree113,pdata$Strains)
write.tree(Phylogenetic_Tree60,'../Data/Ultrametric_Tree_60.phy')


#Calculate phylogenetic distance and phylogenetic signal

data = fread('../Data/Phenotypes_60.csv')
data = data[ ,-c('Replicates')]
data =  data[,lapply(.SD, function(x) as.numeric(mean(x,na.rm=TRUE))),by=Strains ]
t1 =read.newick('../Data/Ultrametric_Tree_60.phy')


phylo_dist = cophenetic(t1)

  #Calculate distance in scaled trait value for environmental preference, fermentation performances and both
  data_a =scale(data[,c('ETHANOL',
                        'ACETIC ACID',
                        'SUCCINIC ACID',
                        'MALIC ACID',
                        'LACTIC ACID',
                        'TARTARIC ACID',
                        'CITRIC ACID',
                        'GLUCOSE',
                        'FRUCTOSE',
                        'GLICERINE',
                        'AMMONIA',
                        "pH",
                        'Total Acidity',
                        'PAN',
                        'SUGARS')])
  data_b= scale(data[,colnames(data) %!in% c(colnames(data_a),'Strains'),with=FALSE])
  data_c = cbind(data_a,data_b)
  rownames(data_a) = data$Strains
  rownames(data_b) = data$Strains
  rownames(data_c) = data$Strains
  
  F_Dist = as.matrix(dist(scale(data_a)))
  E_Dist = as.matrix(dist(scale(data_b)))
  Pheno_Dist = as.matrix(dist(scale(data_c)))
  xy <- t(combn(colnames(F_Dist), 2))
  xy <- data.frame(xy,
                   F_dist=F_Dist[xy],
                   E_dist = E_Dist[xy],
                   P_dist = Pheno_Dist[xy],
                   Phy_dist = phylo_dist[xy])
  cor.test(xy$Phy_dist,xy$P_dist)
  fwrite(xy,'../Data/Distances.csv')
  lambda = c()
  pval = c()
  pval2 = c()
  for(i in 1:ncol(data_c)){
    print(i)
    lambda=c(lambda,phylosig(t1,data_c[,i][!is.na(data_c[,i])],method='lambda')$lambda)
    blambda = c()
    for(j in 1:1000){
      t_d = sample(data_c[,i],replace=FALSE)
      names(t_d) = rownames(data_c)
      blambda = c(blambda,phylosig(t1,t_d,method='lambda')$lambda)
    }
    pval = c(pval,length(blambda[blambda>lambda])/length(blambda))
  }
  lambda_df = data.frame(trait =colnames(data_c),
                         lambda = lambda,
                         p_permutation = pval)
  
  fwrite(lambda_df,'../Data/Phylosig.csv')
  
#Phylogeny imputation_Strain
  
  t1 = read.newick('../Data/Ultrametric_Tree_60.phy')
  pheno = fread('../Data/Phenotypes_60.csv')[,-c(2)]
  
  colnames(pheno)[1] = 'species'
  pheno = pheno[pheno$species %in% t1$tip.label,]
  t1 = keep.tip(t1,(unique(pheno$species)))
  pheno = pheno[,-c(grep('^R\\-',colnames(pheno)),
                    grep('RATIO',colnames(pheno))),with=FALSE]
  observed_mean = pheno
  predicted_mean = data.frame()
  
  for(i in unique(observed_mean$species)){
    print(i)
    temp_observed_mean = observed_mean
    temp_observed_mean[temp_observed_mean$species==i,2:ncol(temp_observed_mean)]<-NA #Remove strain data
    imputed = c()
    for(k in 2:ncol(temp_observed_mean)){
      pobj =phylopars(trait_data = temp_observed_mean[,c(1,k),with=FALSE],tree = t1)
      imputed = c(imputed,pobj$anc_recon[i,])
    }  
    predicted_mean = rbind(predicted_mean,imputed)
  }
  colnames(predicted_mean) = colnames(observed_mean)[-1]
  predicted_mean$species = unique(observed_mean$species)
  final = merge(melt(predicted_mean,value.name='predicted',variable.name='trait'),
                melt(observed_mean,value.name='observed',variable.name='trait'))
  
  fwrite(final,'../Data/60_LOOCV_Pred_Strain.csv')
  
#Phylogeny imputation_Outgroup
  
  t1 = read.newick('../Data/Ultrametric_Tree_113.phy')
  pheno = fread('../Data/Phenotypes_60.csv')
  pheno2 = fread('../Data/Phenotypes_113.csv')
  pheno = pheno[,c('Strains','SGM')]
  pheno2 = pheno2[,c('Strains','E-SGM')]
  
  colnames(pheno)[1] = c('species','SGM')
  colnames(pheno2)[1] = c('species','SGM')
  
  pheno = pheno[pheno$species %in% t1$tip.label,]
  pheno2 = pheno2[pheno2$species %in% t1$tip.label,]
  
  observed_mean = pheno
  observed_mean = observed_mean[,c('species','SGM')]
  observed_mean2 = pheno2
  observed_mean = observed_mean[,c('species','SGM')]
  
  pobj =phylopars(trait_data = observed_mean,tree = t1)
  observed_mean2$predicted = pobj$anc_recon[observed_mean2$species,]
  observed_mean2$Type='Outgroup' #New strains
  observed_mean2[observed_mean2$species%in% observed_mean$species,]$Type = 'Ingroup' #Ingroup
  colnames(observed_mean2)[2] ='observed'
  fwrite(observed_mean2,'../Data/Outgroup_Pred.csv')
  
 #Phylogeny imputation_26S 

  t = readDNAStringSet('../Data/Temporary/26s.fasta')
  t1 = read.newick('../Data/Ultrametric_Tree_60.phy')
  t = t[names(t) %in% t1$tip.label]
  #Pairwise distance to determine identical 26s sequences
  dd = melt(DistanceMatrix(t))
  nn = data.table(dd[dd$value==0 & dd$Var1!=dd$Var2,]) #No mismatches
  pheno = fread('../Data/Phenotypes_60.csv')[,-c(2)]
  
  colnames(pheno)[1] = 'species'
  pheno = pheno[pheno$species %in% t1$tip.label,]
  t1 = keep.tip(t1,(unique(pheno$species)))
  pheno = pheno[,-c(grep('^R\\-',colnames(pheno)),
                    grep('RATIO',colnames(pheno))),with=FALSE]
  observed_mean = pheno
  predicted_mean = data.frame()
  
  for(i in unique(observed_mean$species)){
    print(i)
    temp_observed_mean = observed_mean
    t_i = c(i,as.character(nn[Var1==i,]$Var2),as.character(nn[Var2==i,]$Var1))
    temp_observed_mean[temp_observed_mean$species%in% t_i,2:ncol(temp_observed_mean)]<-NA   #Remove closest relatives (identical 26s) 
    imputed = c()
    for(k in 2:ncol(temp_observed_mean)){
      pobj =phylopars(trait_data = temp_observed_mean[,c(1,k),with=FALSE],tree = t1)
      imputed = c(imputed,pobj$anc_recon[i,])
    }  
    predicted_mean = rbind(predicted_mean,imputed)
  }
  colnames(predicted_mean) = colnames(observed_mean)[-1]
  predicted_mean$species = unique(observed_mean$species)
  final = merge(melt(predicted_mean,value.name='predicted',variable.name='trait'),
                melt(observed_mean,value.name='observed',variable.name='trait'))
  
  fwrite(final,'../Data/60_LOOCV_Pred_26s.csv')
  
#Phylogeny imputation_Species
  
  t1 = read.newick('../Data/Ultrametric_Tree_60.phy')
  taxonomy = fread('../Data/strain_info.csv')
  pheno = fread('../Data/Phenotypes_60.csv')[,-c(2)]
  
  colnames(pheno)[1] = 'species'
  pheno = pheno[pheno$species %in% t1$tip.label,]
  t1 = keep.tip(t1,(unique(pheno$species)))
  pheno = pheno[,-c(grep('^R\\-',colnames(pheno)),
                    grep('RATIO',colnames(pheno))),with=FALSE]
  observed_mean = pheno
  predicted_mean = data.frame()
  
  for(i in unique(observed_mean$species)){
    print(i)
    temp_observed_mean = observed_mean
    t_i = temp_observed_mean$species[grep(strsplit(i,'_')[[1]][1],temp_observed_mean$species)]
    temp_observed_mean[temp_observed_mean$species%in% t_i,2:ncol(temp_observed_mean)]<-NA
    imputed = c()
    for(k in 2:ncol(temp_observed_mean)){
      pobj =phylopars(trait_data = temp_observed_mean[,c(1,k),with=FALSE],tree = t1)
      imputed = c(imputed,pobj$anc_recon[i,])
    }  
    predicted_mean = rbind(predicted_mean,imputed)
  }
  colnames(predicted_mean) = colnames(observed_mean)[-1]
  predicted_mean$species = unique(observed_mean$species)
  final = merge(melt(predicted_mean,value.name='predicted',variable.name='trait'),
                melt(observed_mean,value.name='observed',variable.name='trait'))
  
  fwrite(final,'../Data/60_LOOCV_Pred_Species.csv')
  
#Phylogeny imputation_Genus
  
  t1 = read.newick('../Data/Ultrametric_Tree_60.phy')
  taxonomy = fread('../Data/strain_info.csv')
  taxonomy$Genus = sapply(strsplit(taxonomy$Species,split=' '),function(x) x[[1]])
  pheno = fread('../Data/Phenotypes_60.csv')[,-c(2)]
  
  colnames(pheno)[1] = 'species'
  pheno = pheno[pheno$species %in% t1$tip.label,]
  t1 = keep.tip(t1,(unique(pheno$species)))
  pheno = pheno[,-c(grep('^R\\-',colnames(pheno)),
                    grep('RATIO',colnames(pheno))),with=FALSE]
  observed_mean = pheno
  predicted_mean = data.frame()
  
  for(i in unique(observed_mean$species)){
    print(i)
    temp_observed_mean = observed_mean
    t_i = taxonomy[taxonomy$Genus ==  taxonomy[taxonomy$Code==i,]$Genus,]$Code
    temp_observed_mean[temp_observed_mean$species%in% t_i,2:ncol(temp_observed_mean)]<-NA
    imputed = c()
    for(k in 2:ncol(temp_observed_mean)){
      pobj =phylopars(trait_data = temp_observed_mean[,c(1,k),with=FALSE],tree = t1)
      imputed = c(imputed,pobj$anc_recon[i,])
    }  
    predicted_mean = rbind(predicted_mean,imputed)
  }
  colnames(predicted_mean) = colnames(observed_mean)[-1]
  predicted_mean$species = unique(observed_mean$species)
  final = merge(melt(predicted_mean,value.name='predicted',variable.name='trait'),
                melt(observed_mean,value.name='observed',variable.name='trait'))
  
  fwrite(final,'../Data/60_LOOCV_Pred_Genus.csv')
  
#Phylogeny imputation_FEE
  
  t1 = read.newick('../Data/Ultrametric_Tree_113.phy')
  fees = fread('../Data/fees.txt')
  fees = fees[fees$measurement=='sugars_consumption']
  fees$species = c('Ap_1','Hop_1','Ku_1','Lt_2','Mp_3','Pk_3','Sc_5','Sc_8','Sp_2','Td_8','Wa_5','Zb_1')
  fees = data.frame(fees[,c(3,4,5)])
  t1 = keep.tip(t1,fees$species)
  fees = fees[match(t1$tip.label,fees$species),]
  fees_plot = data.frame(
    Slope = as.numeric(fees$fee_slope),
    Intercept = as.numeric(fees$fee_intercept))
  rownames(fees_plot) = fees$species
  p1a <- ggtree(t1) + geom_tiplab() +geom_rootedge(0.1) 
  p1 <- gheatmap(p1a,scale(fees_plot),low='yellow',high='blue',
                 legend_title='Standardized Value',colnames_position = 'top',width=0.5,offset=0.3,colnames_angle=45,
                 font.size = 3.5,
                 colnames_offset_y =0.35,colnames_offset_x = 0.05) + theme(legend.position = 'bottom',
                                                                           legend.box="horizontal") +
    guides(fill = guide_colourbar(title.position="top", title.hjust = 0.5))
  predicted_mean = data.frame()
  for(i in unique(fees$species)){
    print(i)
    temp_observed_mean = data.frame(fees)
    temp_observed_mean[temp_observed_mean$species==i,2:ncol(temp_observed_mean)]<-NA
    imputed = c()
    pobj_intercept =phylopars(trait_data = temp_observed_mean[,c(1,2)],
                              tree = t1)
    
    pobj_slope =phylopars(trait_data = temp_observed_mean[,c(1,3)],
                          tree = t1)
    imputed = c(pobj_intercept$anc_recon[i,],pobj_slope$anc_recon[i,])
    
    predicted_mean = rbind(predicted_mean,imputed)
  }
  
  
  colnames(predicted_mean) = colnames(fees)[-1]
  predicted_mean$species = unique(fees$species)
  final = merge(melt(fees,value.name='observed',variable.name='trait'),
                melt(predicted_mean,value.name='predicted',variable.name='trait'))
  p2  <- ggplot(final[final$trait=='fee_slope',],aes(x=predicted, y =observed)) + 
    geom_point(size=2,stroke=1,shape=1) +
    theme_classic() +
    geom_smooth(method='lm',se=FALSE,col='Red') +
    stat_cor() + 
    labs(x='Predicted Fee Slope',y='Observed Fee Slope')
  p3  <- ggplot(final[final$trait=='fee_intercept',],aes(x=predicted, y =observed)) + 
    geom_point(size=2,stroke=1,shape=1) +
    theme_classic() +
    geom_smooth(method='lm',se=FALSE,col='Red') +
    stat_cor() + 
    labs(x='Predicted Fee intercept',y='Observed Fee intercept')
  ggsave('../Plots/FigS15.png',ggarrange(p1,ggarrange(p2,p3,ncol=1,nrow=2,labels=c('B','C')),labels=c('A','')),height=7,width=7)
  