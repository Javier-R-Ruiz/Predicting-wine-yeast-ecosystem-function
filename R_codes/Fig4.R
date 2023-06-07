m(list = ls())
library(reshape2)
library(scales)
save_plots <- T

# load data
data <- read.csv('../Source data_figure4/functional_effects_final (panels A, B and C).csv')

#data <- data[!(data$background == 'Ap,Pk,Sc5,Sp,Zb' & data$knock_in == 'Ku'), ] # is this an outlier?
data$knock_in <- factor(data$knock_in,
                        levels = c('Sc5', 'Sc8',
                                   'Sp', 'Lt', 'Zb',
                                   'Td', 'Ku', 'Hop',
                                   'Wa', 'Mp', 'Pk', 'Ap'))

# plot FEEs of Sc5 and Sc8 (PANEL A)
ggplot(data[data$knock_in %in% c('Sc5', 'Sc8'), ],
       aes(x = background_fun.mean, y = delta_fun.mean,
           xmin = background_fun.mean - background_fun.sd, xmax = background_fun.mean + background_fun.sd,
           ymin = delta_fun.mean - delta_fun.sd, ymax = delta_fun.mean + delta_fun.sd)) +
  geom_abline(slope = 0, intercept = 0,
              color = '#d1d3d4') +
  geom_abline(slope = -1, intercept = 1,
              color = '#d1d3d4',
              linetype = 'dashed') +
  geom_point(color = 'black',
             cex = 1.5,
             shape = 16) +
  geom_errorbar(alpha = 0.25) +
  geom_errorbarh(alpha = 0.25) +
  geom_smooth(method = 'lm',
              formula = y~x,
              color = 'firebrick1',
              se = F,
              fullrange = T) +
  facet_wrap(~ knock_in) +
  scale_x_continuous(name = 'F (background)',
                     limits = c(0, 1),
                     breaks = pretty_breaks(n = 3)) +
  scale_y_continuous(name = 'dF',
                     breaks = pretty_breaks(n = 3)) +
  theme_bw() +
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 16, face = 'italic'),
        panel.grid = element_blank(),
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
  ggsave(filename = '../plots/sugars_consumption_scFEEs.pdf',
         device = 'pdf',
         dpi = 600,
         width = 200,
         height = 80,
         units = 'mm',
         limitsize = F)
}


# plot every other fee, split backgrounds by presence/absence of Sc (Panel B)
plot_this <- data[!(data$knock_in %in% c('Sc5', 'Sc8')), ]
plot_this$Sc <- c('No', 'Yes')[1 + grepl('Sc', plot_this$background)]
plot_this$Sc <- factor(plot_this$Sc, levels = c('Yes', 'No'))
ggplot(plot_this,
       aes(x = background_fun.mean, y = delta_fun.mean, color = Sc,
           xmin = background_fun.mean - background_fun.sd, xmax = background_fun.mean + background_fun.sd,
           ymin = delta_fun.mean - delta_fun.sd, ymax = delta_fun.mean + delta_fun.sd)) +
  geom_abline(slope = 0, intercept = 0,
              color = '#d1d3d4') +
  geom_abline(slope = -1, intercept = 1,
              color = '#d1d3d4',
              linetype = 'dashed') +
  geom_point(cex = 1.5,
             shape = 16) +
  geom_errorbar(alpha = 0.25) +
  geom_errorbarh(alpha = 0.25) +
  geom_smooth(method = 'lm',
              formula = y~x,
              color = 'firebrick1',
              se = F,
              fullrange = T) +
  facet_wrap(~ knock_in,
             nrow = 2) +
  scale_x_continuous(name = 'F (background)',
                     limits = c(0, 1),
                     breaks = pretty_breaks(n = 3)) +
  scale_y_continuous(name = 'dF',
                     breaks = pretty_breaks(n = 2)) +
  scale_color_manual(name = expression(paste(italic(Sc), ' in background?')),
                     values = c('deepskyblue', 'gray')) +
  theme_bw() +
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 16, face = 'italic'),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        aspect.ratio = 0.6) +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, size=0.5) +
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf,size=0.5)

if (save_plots) {
  ggsave(filename = '../plots/sugars_consumption_FEEs_split-by-Sc.pdf',
         device = 'pdf',
         dpi = 600,
         width = 230,
         height = 80,
         units = 'mm',
         limitsize = F)
}


# plot of model accuracy vs community richness (Panel C)
# keep only communities that contain Sc
plot_this <- loo[grepl('Sc', loo$community) & loo$n_species > 1, ]
plot_this$n_species <- plot_this$n_species - 1
plot_this$n_species <- factor(plot_this$n_species,
                              levels = as.character(c(1, 3:6)))
ggplot(plot_this, aes(x = n_species, y = sq_err, group = n_species)) +
  geom_jitter(width = 0.15,
              alpha = 0.25,
              shape = 16) +
  geom_boxplot(outlier.shape = NA,
               fill = NA) +
  scale_x_discrete(name = '# of species\nco-inoculated with S. cerevisiae') +
  scale_y_continuous(name = expression((italic(F)[pred] - italic(F)[obs])^2),
                     limits = c(0, 0.6),
                     breaks = c(0, 0.2, 0.4, 0.6)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        legend.title = element_blank(),
        legend.text = element_text(size = 16),
        aspect.ratio = 1.4)

if (save_plots) {
  ggsave(filename = '../plots/squared_error_vs_richness.pdf',
         device = 'pdf',
         dpi = 600,
         width = 100,
         height = 120,
         units = 'mm',
         limitsize = F)
}


# load data
data <- read.csv('../Source data_figure4/observed_vs_predicted (panel D).csv')


#Plot correlation observed vs predicted_community function (Panel D)

predicted_value <- as.numeric(data$Fun)
value <- as.numeric(data$PROMEDIO)


ggplot(data,
       aes(x = predicted_value, y = value)) +
  geom_smooth(method = 'lm',
              formula = y~x,
              se = F,
              color = 'red',
              fullrange = T) +
  geom_point(cex = 3,
             shape = 1) +
  scale_x_continuous(name = 'Predicted fraction of sugars consumed',
                     limits = c(0, 1), breaks=c(0,0.2,0.4,0.6,0.8,1)) +
  scale_y_continuous(name = 'Observed fraction of sugars consumed',
                     limits = c(0, 1), breaks=c(0,0.2,0.4,0.6,0.8,1)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 16,
                                  hjust = 0),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 16),
        aspect.ratio = 1)


ggsave("OBS vs Pred_out of samples.png", 
       dpi = 600,
       width = 120, 
       height = 120, 
       units = 'mm',
       limitsize = F)

