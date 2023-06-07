rm(list=ls())
library(data.table)
library(ggplot2)
library(patchwork)
dt =  fread('../Data/60_LOOCV_Pred_Strain (EFGH).csv')
phylosig = fread('../Data/Phylosig.csv')
phylo_dist=fread('../Data/Distances.csv')
dt = dt[,lapply(.SD,mean,na.rm=TRUE),by=list(trait,species)]
#Scale traits to calculate predicted vs observed distance for each species pair.
dt[,scaled_predicted:= (predicted-mean(predicted,na.rm=TRUE))/sd(predicted,na.rm=TRUE), by='trait']
dt[,scaled_observed:= (observed-mean(observed,na.rm=TRUE))/sd(observed,na.rm=TRUE), by='trait']
newdt = dt[,  cor.test(predicted, observed,use='complete.obs')[c('estimate','p.value')], by='trait']
newdt2 = dt[, sqrt(sum(abs(scaled_predicted- scaled_observed)^2,na.rm=TRUE)),by='species']
newdt = merge(newdt,phylosig)
newdt = newdt[rev(order(newdt$estimate)),]

newdt$trait = factor(newdt$trait,levels=newdt$trait)

for(i in 1:nrow(newdt2)){
  S = newdt2$species[i]
  newdt2$NN[i] =  min(phylo_dist[X1==S | X2 == S]$Phy_dist,na.rm = TRUE)
  newdt2$MN[i] =  mean(phylo_dist[X1==S | X2 == S]$Phy_dist,na.rm=TRUE)
}
newdt$trait = factor(newdt$trait,levels=newdt$trait)
newdt$signif = factor(newdt$p.value<0.05,levels=c(FALSE,TRUE))
p1 <-ggplot(newdt,aes(x=trait,y=estimate,fill=signif)) +geom_bar(stat='identity',color='black') +theme_classic() + 
  geom_hline(yintercept=1,linetype=2,size=2) +
  labs(x = '', y = expression("Pearson "*r*" (LOOCV)"),fill=expression(P<0.05)) +
  scale_fill_manual(values=c('Red','Grey'),drop=FALSE) +
  scale_y_continuous(breaks=c(0,1),expand=c(0,0),limits=c(0,1)) +
  theme(axis.text.x=element_text(angle=90,size=8),
        legend.position = c(0.9, 0.88),
        legend.background = element_rect(fill = "white", color = "black"),
        legend.title = element_text( size=8), legend.text=element_text(size=8),
        legend.key.size = unit(0.5, "lines")) +
  theme(axis.line=element_line(size=1),
        axis.text.y=element_text(size=12),
        axis.title=element_text(size=20))



p2<-ggplot(newdt,aes(x=lambda,y=estimate))  +
  geom_point(shape=1,size=3,stroke=2.0) + 
  theme_classic() + 
  labs(x=expression('Phylogenetic Signal '*lambda),y = expression("Pearson "*r*" (LOOCV)")) + 
  geom_smooth(method='glm',col='Red',se=FALSE,size=1.5) + 
  scale_y_continuous(limits=c(0.1,1.1),breaks=c(0.3,1)) +
  scale_x_continuous(limits=c(0.2,1.1),breaks=c(0.3,1))+
  stat_cor(method = "pearson", label.x = 0.5, label.y =1.1,size=5) +
  
  theme(axis.line=element_line(size=1),
        axis.text=element_text(size=12),
        axis.title=element_text(size=20))

p3<-ggplot(newdt2,aes(x=NN,y=V1)) + 
  geom_point(shape=1,size=1.5,stroke=1.0) + 
  theme_classic() + 
  labs(x=expression('PNND'),y = expression("    Phenotypic Distance \n (Predicted vs Observed)")) + 
  geom_smooth(method='glm',col='Red',se=FALSE,size=1.5) + 
  scale_x_continuous(limits=c(0.0,1.75),breaks=c(0.0,1.5)) +
  scale_y_continuous(limits=c(0,16),breaks=c(0,15)) +
  stat_cor(method = "pearson", label.x = 0., label.y =15,size=2.5) +
  theme(axis.line=element_line(size=1),
        axis.text=element_text(size=6),
        axis.title=element_text(size=10))
p2 <- p2 +inset_element(p3,left=0.05,right=0.4,bottom=0.5,top=0.9)
ggsave('../Plots/Fig2GH.png',ggarrange(p1,p2,ncol=1,nrow=2,labels=c('G','H')),height=10,width=6)
