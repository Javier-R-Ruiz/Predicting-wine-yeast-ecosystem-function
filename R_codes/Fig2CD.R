rm(list=ls())
library(data.table)
library(ggplot2)
library(operators)
library(ggpubr)
xy = fread('../Data/Distances (C).csv')

p1<-ggplot(xy,aes(x=Phy_dist,y=P_dist)) + geom_point(shape=1,size=3,stroke=1.0) + 
  theme_classic() + 
  labs(x='Phylogenetic Distance',y = 'Phenotypic Distance') + 
  geom_smooth(method='lm',col='Red',se=FALSE,size=1.5) + 
  stat_cor(method = "pearson", label.x = 0.6, label.y = max(xy$P_dist),size=5) +
  scale_y_continuous(limits=c(0,21),breaks=c(0,20)) +
  scale_x_continuous(limits=c(0,2.1),breaks=c(0,2)) +
  
  theme(axis.line=element_line(size=1),
        axis.text=element_text(size=12),
        axis.title=element_text(size=20))



p2<-ggplot(xy,aes(x=Phy_dist,y=F_dist)) + 
  geom_point(shape=1,size=3,stroke=1.0) +
  theme_classic() +  
  labs(x='Phylogenetic Distance',y = 'Phenotypic Distance \n (Functions Only)') + 
  geom_smooth(method='lm',col='Red',se=FALSE,size=1.5) + 
  stat_cor(method = "pearson", label.x = 0.6, label.y = max(xy$P_dist),size=5) +
  theme(axis.line=element_line(size=1),
        axis.text=element_text(size=12),
        axis.title=element_text(size=15))



p3<-ggplot(xy,aes(x=Phy_dist,y=E_dist))  + geom_point(shape=1,size=3,stroke=1.0) +
  theme_classic() +  
  labs(x='Phylogenetic Distance',y = 'Phenotypic Distance \n (Efficiency only)') + 
  geom_smooth(method='lm',col='Red',se=FALSE,size=1.5) + 
  stat_cor(method = "pearson", label.x = 0.6, label.y = max(xy$P_dist),size=5) +
  theme(axis.line=element_line(size=1),
        axis.text=element_text(size=12),
        axis.title=element_text(size=15))

ggsave('../Plots/FigS3.png',ggarrange(p3,p2,ncol=1,nrow=2,labels=c('A','B')),height=8,width=5)

lambda_df = fread('../Data/Phylosig.csv')
lambda_df = lambda_df[rev(order(lambda_df$lambda)),]

lambda_df$trait = factor(lambda_df$trait,levels=lambda_df$trait)
lambda_df$Signif = lambda_df$p_permutation<0.05
lambda_df$Signif = factor(lambda_df$Signif,levels=c(TRUE,FALSE))
p4 <-ggplot(lambda_df,aes(x=trait,y=lambda,fill=Signif)) +
  geom_bar(stat='identity',color='black') +theme_classic() + 
  geom_hline(yintercept=1,linetype=2,size=2) +
  labs(x = '', y = expression("Phylogenetic Signal ("*lambda*")"),fill=expression(P<0.05)) +
  scale_fill_manual(values=c('Grey','Red'),drop=FALSE) +
  scale_y_continuous(breaks=c(0,1),expand=c(0,0),limits=c(0,1.02)) +
  theme(axis.text.x=element_text(angle=90,size=8),
        legend.position = c(0.9, 0.88),
        legend.background = element_rect(fill = "white", color = "black"),
        legend.title = element_text( size=8), legend.text=element_text(size=8),
        legend.key.size = unit(0.5, "lines")) +
  theme(axis.line=element_line(size=1),
        axis.text.y=element_text(size=12),
        axis.title=element_text(size=20)) 

ggsave('../Plots/Fig2CD.png',ggarrange(p1,p4,ncol=1,nrow=2,labels=c('C','D')),height=10,width=6)


