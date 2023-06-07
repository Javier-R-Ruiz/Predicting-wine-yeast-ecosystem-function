library(reshape2)
library(scales)

# load data
data <- read.csv('../data/functional_effects_final.csv')
#data <- data[!(data$background == 'Ap,Pk,Sc5,Sp,Zb' & data$knock_in == 'Ku'), ] # is this an outlier?
data$knock_in <- factor(data$knock_in,
                        levels = c('Sc5', 'Sc8',
                                   'Sp', 'Lt', 'Zb',
                                   'Td', 'Ku', 'Hop',
                                   'Wa', 'Mp', 'Pk', 'Ap'))

# get list of community functions
data_list <- data.frame(community = sapply(1:nrow(data),
                                           FUN = function(i) {
                                             x <- c(as.character(data$knock_in[i]), strsplit(data$background[i], split = ',')[[1]])
                                             x <- sort(x)
                                             x <- paste(x, collapse = ',')
                                             return(x)
                                           }),
                        fun = data$background_fun.mean + data$delta_fun.mean)
data_list <- aggregate(formula = fun ~ community,
                       data = data_list,
                       FUN = mean)

# plot histogram
data_list$sc <- 'No Sc'
data_list$sc[grepl('Sc5', data_list$community)] <- 'Sc5'
data_list$sc[grepl('Sc8', data_list$community)] <- 'Sc8'
data_list$sc <- c('Without Saccharomyces', 'With Saccharomyces')[1 + grepl('Sc', data_list$community)]

ggplot(data_list, aes(x = fun, group = sc, color = sc, fill = sc)) +
  stat_density(geom = 'line',
               position = 'identity',
               size = 1) +
  geom_histogram(aes(y = ..density..),
                 binwidth = 0.025,
                 color = NA,
                 alpha = 0.25,
                 position = 'identity') +
  scale_x_continuous(name = 'Fraction of sugars consumed',
                     breaks = seq(0, 1, by=0.2)) +
  scale_y_continuous(name = '# of communities') +
  scale_color_manual(values = c('#f26522', '#8c3092')) +
  scale_fill_manual(values = c('#f26522', '#8c3092')) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        panel.border = element_blank(),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        legend.title = element_blank(),
        legend.text = element_text(size = 16),
        legend.position = c(0.42, 0.9),
        aspect.ratio = 0.6) +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, size=0.5) +
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf,size=0.5)

if (save_plots) {
  ggsave(filename = '../plots/sugars_consumption_histogram.pdf',
         device = 'pdf',
         dpi = 600,
         width = 100,
         height = 80,
         units = 'mm',
         limitsize = F)
}
