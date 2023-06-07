rm(list=ls())
library(scales)
library(ggplot2)
library(data.table)
data = fread('../Data/60_LOOCV_Pred_Strain (EFGH).csv')
data = data[, sapply(.SD, function(x) list(mean = mean(x,na.rm=TRUE))), .SDcols = c('predicted','observed'), by = c('species','trait')]
colnames(data) = c('species','trait','predicted','observed')
newdt = data[,  cor.test(predicted, observed,use='complete.obs')[c('estimate','p.value')], by='trait']
newdt = newdt[order(newdt$estimate)]
data$trait = factor(data$trait,levels=rev(newdt$trait))
data$Signif = TRUE
data[data$trait %in% newdt[newdt$p.value>0.05,]$trait]$Signif=FALSE


p1 <- ggplot(data)+
  geom_point(mapping = aes(x=predicted,y=observed,col=Signif),
             size=3,stroke=1.5,shape=1) +
  scale_color_manual(values=c('Red','Black')) + guides(col=FALSE) +
  geom_point(mapping = aes(x = observed*1.2, y = predicted*1.2), alpha = 0) + 
  geom_point(mapping = aes(x = predicted*1.2, y = observed*1.2), alpha = 0) + 
  scale_x_continuous(breaks = trans_breaks(identity, identity, n = 2)) +
  scale_y_continuous(breaks = trans_breaks(identity, identity, n = 2)) +
  stat_cor(mapping = aes(x=predicted,y=observed),
           label.y.npc=1, label.x.npc = 0, method = "pearson",size=3) +
  labs(x='Predicted Trait Value (LOOCV)', y= 'Observed Trait Value')+
  theme_classic()  + 
  geom_abline(intercept=0,slope=1,col='Red',linetype=2,size=1) +
  theme(aspect.ratio = 1,axis.line=element_line(size=1.0),axis.text=element_text(size=12),axis.title = element_text(size=20))+
  facet_wrap(~trait,scales='free',nrow=7,ncol=7) +guides(col=FALSE)
ggsave('../Plots/FigS5.png',p1,height=15,width=15)

p2 <-ggplot( data[data$trait=='SGM'],aes(x=predicted,y=observed))+
  geom_point(shape=1,size=3,stroke=2.0) + 
  theme_classic() + 
  labs(x='Predicted  \n Efficiency in SGM (LOOCV)', y= 'Observed  \n Efficiency in SGM')+
  # geom_smooth(method='glm',col='Red',se=FALSE,size=1.5) + 
  geom_abline(yintercept=0,slope=1,col='Red',linetype=2,size=1.5)+
  stat_cor(method = "pearson", label.x = 0.4, label.y =1.5,size=5) +
  scale_y_continuous(limits=c(0,1.6),breaks=c(0,1.5)) +
  scale_x_continuous(limits=c(0,1.6),breaks=c(0,1.5))+
  theme(axis.line=element_line(size=1),
        axis.text=element_text(size=12),
        axis.title=element_text(size=20))

p3 <-ggplot( data[data$trait=='SUGARS'],aes(x=predicted,y=observed)) +
  geom_point(shape=1,size=3,stroke=2.0) + 
  theme_classic() + 
  labs(x='Predicted   \n Sugar Concentration (LOOCV)',y= 'Observed   \n Sugar Concentration')+
  geom_abline(yintercept=0,slope=1,col='Red',linetype=2,size=1.5)+
  stat_cor(method = "pearson", label.x = 50, label.y =205,size=5) +
  scale_y_continuous(limits=c(0,205),breaks=c(0,200)) +
  scale_x_continuous(limits=c(0,205),breaks=c(0,200)) +
  theme(axis.line=element_line(size=1),
        axis.text=element_text(size=12),
        axis.title=element_text(size=20))


ggsave('../Plots/Fig2EF.png',ggarrange(p2,p3,ncol=1,nrow=2,labels=c('E','F')),height=10,width=6)
