rm(list=ls())
library(data.table)
library(ggplot2)
library(operators)
library(ggtree)
library(phytools)
library(pheatmap)
library(ggplotify)
library(RColorBrewer)
library(gtable)
library(grid)
library(ggpubr)

  
#Read in Data
data = fread('../Data/Phenotypes_60(AB).csv')
data = data[ ,-c('Replicates')]
data =  data[,lapply(.SD, function(x) as.numeric(mean(x,na.rm=TRUE))),by=Strains ]
taxonomy = fread('../Data/strain_info.csv')
t1 =read.newick('../Data/Ultrametric_Tree_60.phy')
p1 <- ggtree(t1) +geom_rootedge(0.1) 

# tdf = merge(data.frame(ID = t1$tip.label,'Initials' = substr(t1$tip.label,1,2)),taxonomy)
mydf <- data.frame(row.names = taxonomy$Code , category =taxonomy$Family)
strains = data$Strains
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
data_b =scale(data[,which(colnames(data) %!in% colnames(data_a))[-1],with=FALSE])
  
  
annotation_colors = list(category= c(Dothioraceae    ='LightGreen',
                                       Pichiaceae  ='Orange',
                                       Cryptococcaceae  ='DarkGreen',
                                       Phaffomycetaceae  ='Orange',
                                       Debaryomycetaceae ='Orange',
                                       Saccharomycetaceae = 'Purple',
                                     Metschnikowiaceae = 'Orange',
                                     Metchsnikowiaceae = 'Orange',
                                     
                                       Filobasidiaceae = 'DarkGreen',
                                       Sporidiobolaceae = 'DarkGreen',
                                        Tremellaceae ='DarkGreen',
                                       "Incertae sedis" ='Orange',
                                       Schizosaccharomycetaceae ='LightGreen'))
  
  # 
rownames(data_a) = strains
rownames(data_b)=strains
plotdata = p1$data[1:60,]
data_a = data_a[plotdata[rev(order(plotdata$y)),]$label,]
data_b = data_b[plotdata[rev(order(plotdata$y)),]$label,]
  
  
p_a <- pheatmap(data_a,
                cluster_rows=FALSE,
                cluster_cols=FALSE,
                 show_rownames = FALSE,na_col='Grey', 
                 annotation_row = mydf,
                 annotation_colors = annotation_colors,
                 annotation_names_row=FALSE,
                 annotation_legend=FALSE,
                 color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
                 breaks=seq(-2,2,length.out=100), angle_col = 90)
  
  
p_b <- pheatmap(data_b,
                 cluster_rows=FALSE,
                 cluster_cols=FALSE,
                 
                 # clustering_distance_cols = 'euclidean',
                 show_rownames = FALSE,na_col='Grey', 
                 annotation_row = mydf,
                 # clustering_method = 'average',
                 annotation_colors = annotation_colors,
                 annotation_names_row=FALSE,
                 annotation_legend=FALSE,
                 color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
                 breaks=seq(-2,2,length.out=100), angle_col = 90)
  
  
#Get grobs we want - will use these to create own plot later
plot_a.grob <- p_a$gtable$grob[[1]]
xlab_a.grob <- p_a$gtable$grob[[2]]  
ylab_a.grob <- p_a$gtable$grob[[3]]  
legend_a.grob <- p_a$gtable$grob[[4]]  
  
  
  
#Get grobs we want - will use these to create own plot later
plot_b.grob <- p_b$gtable$grob[[1]]
xlab_b.grob <- p_b$gtable$grob[[2]]  
ylab_b.grob <- p_b$gtable$grob[[3]]  
legend_b.grob <- p_b$gtable$grob[[4]]  

  
#Shift both down by 1 inch
legend_a.grob$children[[1]]$y <- legend_a.grob$children[[1]]$y - unit(2.85,"inches") 
legend_a.grob$children[[2]]$y <- legend_a.grob$children[[2]]$y - unit(2.85,"inches") 
legend_a.grob$children[[1]]$x <- legend_a.grob$children[[1]]$x + unit(0.4,"inches") 
legend_a.grob$children[[2]]$x <- legend_a.grob$children[[2]]$x + unit(0.4,"inches") 


#Shift both down by 1 inch
legend_b.grob$children[[1]]$y <- legend_b.grob$children[[1]]$y - unit(2.85,"inches") 
legend_b.grob$children[[2]]$y <- legend_b.grob$children[[2]]$y - unit(2.85,"inches") 
legend_b.grob$children[[1]]$x <- legend_b.grob$children[[1]]$x + unit(0.4,"inches") 
legend_b.grob$children[[2]]$x <- legend_b.grob$children[[2]]$x + unit(0.4,"inches") 


#New legend label grob
leg_label_a <- textGrob(expression(frac(F[i*alpha] - bar(F[i*alpha]),sigma['  '*F[i*alpha]])),
                      x=0,y=0.7,hjust=-0.6,vjust=0,gp=gpar(fontsize=10,fontface="bold"))

#New legend label grob
leg_label_b <- textGrob(expression(frac(E[i*alpha] - bar(E[i*alpha]),sigma['  '*E[i*alpha]])),
                    x=0,y=0.7,hjust=-0.6,vjust=0,gp=gpar(fontsize=10,fontface="bold"))


#Add label to legend grob
legend_a.grob2 <- addGrob(legend_a.grob,leg_label_a)
legend_b.grob2 <- addGrob(legend_b.grob,leg_label_b)

  

twidths = p_a$gtable$widths
twidths[5] = unit.c(max(unit(1,"grobwidth",legend_a.grob),unit(12,"bigpts")+1.2*unit(1.1,"grobwidth",legend_a.grob)) + unit(1,"inches"))
theights = p_a$gtable$heights
my_new_gt <- gtable(widths = twidths,
                    heights=theights)

gtable_a <- gtable_add_grob(my_new_gt,plot_a.grob,4,3)
# gtable <- gtable_add_grob(gtable,tree.grob,2,3)
gtable_a <- gtable_add_grob(gtable_a,ylab_a.grob,4,2)
gtable_a <- gtable_add_grob(gtable_a,legend_a.grob2,4,5)
gtable_a <- gtable_add_grob(gtable_a,xlab_a.grob,5,3)

gtable_b <- gtable_add_grob(my_new_gt,plot_b.grob,4,3)
# gtable <- gtable_add_grob(gtable,tree.grob,2,3)
gtable_b <- gtable_add_grob(gtable_b,ylab_b.grob,4,2)
gtable_b <- gtable_add_grob(gtable_b,legend_b.grob2,4,5)
gtable_b <- gtable_add_grob(gtable_b,xlab_b.grob,5,3)

p2 <- as.ggplot(gtable_a) +theme(plot.margin = margin(t = 1,b = 1,l = 10,r = 1))
p3 <- as.ggplot(gtable_b) +theme(plot.margin = margin(t = 1,b = 1,l = 10,r = 1))
p1 <- p1 +theme(plot.margin = margin(t = 5,b = 82,l = 0,r = 2))


ggsave('../Plots/Fig2AB.png',ggarrange(p1,p3,p1,p2,ncol=4,labels=c('A','','B',''),
                                                       widths=c(0.7,1.5,0.7,1.5),
                                                       font.label=list(size=20)),
       height=10,
     width=15,
     bg='White')
