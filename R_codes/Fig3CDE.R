library(reshape2)
library(scales)
library(ggplot2)

# load data (panel C)
data <- read.csv('../Source data_figure3/communities_sugars (panel C).csv')

data<- subset(data, data$Sc!="woSc",)

rich <- as.factor(data$`order_interaction (with Sc)`)
obs <- as.numeric(data$fraction_sugars_consumed)

#plot (panel C)
ggplot(data) +
  geom_jitter(aes(x=rich, y=obs, color=rich), width=0.3, size =2, alpha =0.3, show.legend=F) +
  geom_boxplot(aes(x=rich, y=obs), width=0.3, size =0.9, alpha=0.2, outlier.shape = NA) +
  scale_y_continuous(limits=c(0,1), breaks = seq(0, 1, by=0.2)) + 
  ylab("Fraction of sugars consumed") +
  xlab(expression(plain('# of strains co-inoculated with')~italic('Sc'))) + 
  scale_color_manual(values = c("#8c3092","#8c3092","#8c3092","#8c3092","#8c3092","#8c3092")) +
  theme_bw() +
  theme(axis.text.x=element_text(size=16, angle = 0, face="plain"),axis.title.x=element_text(size=16),
        axis.text.y=element_text(size=16, angle = 0, face="plain"),axis.title.y=element_text(size=16),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position = "none", legend.title =element_text(size=12, face='bold'), 
        legend.text = element_text(size=8, face='bold', color="black"))  

ggsave("Figure3C.png", dpi = 600, width = 100, height = 120, units = 'mm',
       limitsize = F)


# load data (panel D)
df <- read.csv('../Source data_figure3/df_sugars (panel D).csv')

pw$Sc_5 <- NULL
pw$Sc_8 <- NULL

head(pw)

strains <- as.factor(pw$strains)
Sc5 <- as.numeric(pw$Sc5)
Sc8 <- as.numeric(pw$Sc8)


std_pw <- pw
std_pw$Sc5 <- (pw$Sc5 - mean(pw$Sc5))/sd(pw$Sc5)
std_pw$Sc8 <- (pw$Sc8 - mean(pw$Sc8))/sd(pw$Sc8)


df_pw <- melt(std_pw, id.vars = 'strains')


df_pw$strains <- factor(df_pw$strains,
                        levels = c('Ap','Hop','Ku','Mp','Lt','Pk','Sp','Td','Wa','Zb'))


ggplot(df_pw, aes(x = variable, y = strains, fill = value)) +
  geom_tile(color = 'black') +
  scale_x_discrete(name = '',
                   position = 'top') +
  scale_fill_gradientn(colours = brewer.pal(4, "BuPu"))+
  scale_y_discrete(limits = rev) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(size = 16,
                                   angle = 45,
                                   hjust = 0),
        axis.ticks.y = element_blank(),
        axis.text.y = element_text(size = 16,
                                   face = 'italic'),
        axis.title.y = element_blank(),
        panel.border = element_blank(),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 16),
        aspect.ratio = 5.5)

ggsave(filename = 'Figure3D.png',
       device = 'png',
       dpi = 1000,
       width = 200,
       height = 150,
       units = 'mm',
       limitsize = F)


# load data (panel E)

df <- read.csv('../Source data_figure3/df_sugars (panel D).csv')

head(df)        

df<- subset(df, df$Sc!="-",)
df<- subset(df, df$prediction_add_model!="-",)
df<- subset(df,df$`Order_interaction (Sc)`!="0",)
df<- subset(df,df$`Order_interaction (Sc)`!="1",)

se <- as.numeric(df$`(pre-obs)^2`)
rich <- as.factor(df$`Order_interaction (Sc)`)


#plot (panel E)

ggplot(DF) +
  geom_jitter(aes(x=rich, y=se, color=rich), width=0.3, size =2, alpha =0.3, show.legend=F) +
  geom_boxplot(aes(x=rich, y=se), width=0.3, size =0.9, alpha=0.2, outlier.shape = NA) +
  scale_y_continuous(limits=c(0,0.6), breaks=c(0,0.2,0.4,0.6)) +
  ylab(bquote((F[pred] - F[obs])^2)) +
  xlab(expression(plain('# of strains co-inoculated with')~italic('Sc'))) + 
  scale_color_manual(values = c("#8c3092","#8c3092","#8c3092","#8c3092","#8c3092","#8c3092")) +
  theme_bw() +
  theme(axis.text.x=element_text(size=16, angle = 0, face="plain"),axis.title.x=element_text(size=16),
        axis.text.y=element_text(size=16, angle = 0, face="plain"),axis.title.y=element_text(size=16),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position = "none", legend.title =element_text(size=12, face='bold'), 
        legend.text = element_text(size=8, face='bold', color="black"))  

ggsave("Figure3E.png", dpi = 600, width = 105, height = 120, units = 'mm',
       limitsize = F)

