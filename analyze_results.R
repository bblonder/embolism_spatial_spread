library(ggplot2)
library(dplyr)
library(tidyr)
library(sjPlot)
library(RColorBrewer)

dir_results = 'results'

files_vulnerability = dir(dir_results, pattern='*ts\\.csv$')

process_file_vulnerability <- function(file_this) {
  
  file_parts = strsplit(gsub('ts\\.csv','',file_this), split='_')[[1]]
  
  df = read.csv(file.path(dir_results, file_this))
  
  df_long = df %>% 
    apply(2, cumsum) %>%
    as.data.frame %>%
    mutate(time=row_number()) %>%
    pivot_longer(!time, values_to = 'n_pixels_cumulative') %>%
    mutate(species=paste(file_parts[1:2],collapse='_'),
           treatment=file_parts[3],
           focal_distance_max=file_parts[5],
           leaf_replicate=as.numeric(file_parts[4])) %>%
    select(-name)
  
  return(df_long)
}

df_all_vulnerability = do.call('rbind',lapply(files_vulnerability, process_file_vulnerability)) %>%
  group_by(species, treatment, focal_distance_max) %>%
  # renumber to 1-3
  mutate(leaf_replicate=as.numeric(factor(leaf_replicate))) %>%
  ungroup %>%
  as.data.frame
# do normalization
df_all_vulnerability_for_plotting = df_all_vulnerability %>%
  ungroup %>%
  group_by(species, treatment, leaf_replicate, focal_distance_max) %>%
  mutate(frac_pixels_cumulative = n_pixels_cumulative / max(n_pixels_cumulative)) %>%
  mutate(frac_time = time/max(time)) %>%
  mutate(time_hours = time*15/60) %>%
  ungroup %>%
  mutate(species=factor(species, labels=c('A. rubrum (palmate)', 'C. erythrophyllum (pinnate)'))) %>%
  mutate(treatment=factor(treatment, labels=c('control','lamina','midrib')))

# 
# ggplot(df_all_vulnerability_for_plotting, aes(x=time,y=n_pixels_cumulative, color=focal_distance_max)) +
#   geom_line() +
#   theme_bw() + 
#   facet_grid(species~treatment + leaf_replicate,labeller=label_both) +
#   theme(legend.position='bottom')

g_frac_pixels_cumulative = ggplot(df_all_vulnerability_for_plotting, aes(x=time_hours,y=frac_pixels_cumulative,
                                                            linetype=factor(leaf_replicate),
                                                            color=treatment)) +
  geom_line() +
  theme_bw() + 
  facet_grid(~species,labeller=label_both) +
  theme(legend.position='bottom') +
  scale_linetype(name='Leaf replicate') +
  scale_color_discrete(name='Treatment') +
  xlab('Time (hours)') +
  ylab('Cumulative sum # pixels embolized / Total # pixels embolized')
ggsave(g_frac_pixels_cumulative, file='results/g_frac_pixels_cumulative_whole_leaf.png',width=6,height=4)



# get total lamina pixels
df_px_counts = do.call("rbind",lapply(dir(path='results',pattern='*px_count\\.csv',full.names = TRUE), function(file_this) {
  
  file_parts = strsplit(gsub('px_count\\.csv','',basename(file_this)), split='_')[[1]]
  
  
  z = read.csv(file_this,header = FALSE)
  names(z) = c('px_mask','px_veins','px_total')
  
  z = z %>%
    mutate(species=paste(file_parts[1:2],collapse='_'),
           treatment=file_parts[3],
           focal_distance_max=file_parts[5],
           leaf_replicate=as.numeric(file_parts[4])) %>%
  
  return(z)
  })) %>%
  group_by(species, treatment, focal_distance_max) %>%
  # renumber to 1-3
  mutate(leaf_replicate=as.numeric(factor(leaf_replicate))) %>%
  ungroup %>%
  as.data.frame


# get cumulative embolized pixels
df_px_counts_joined = df_all_vulnerability %>% 
  group_by(species, treatment, leaf_replicate, focal_distance_max) %>% 
  summarize(n_pixels_embolized=max(n_pixels_cumulative)) %>%
  left_join(df_px_counts,by=c('species','treatment','leaf_replicate','focal_distance_max')) %>%
  # now drop the focal distance part - not needed
  group_by(species, treatment, leaf_replicate) %>%
  slice_head(n=1) %>%
  ungroup %>%
  mutate(species=factor(species, labels=c('A. rubrum (palmate)', 'C. erythrophyllum (pinnate)'))) %>%
  mutate(treatment=factor(treatment, labels=c('control','lamina','midrib')))

g_frac_px_total = ggplot(df_px_counts_joined, aes(x=species,y=n_pixels_embolized/px_mask,color=treatment)) +
  geom_boxplot() +
  theme_bw() + 
  xlab('Species') +
  ylab('Total pixels embolized / total pixels') +
  scale_color_discrete(name='Treatment')
ggsave(g_frac_px_total, file='results/g_frac_px_total_whole_leaf.png',width=6,height=4)


# do linear model
m_total_pixels_embolized = lm(n_pixels_embolized/px_mask ~ species+treatment, data=df_px_counts_joined %>% filter(focal_distance_max==200))
summary(m_total_pixels_embolized)
anova(m_total_pixels_embolized)
# write it out
tab_model(m_total_pixels_embolized, file='results/tab_m_total_pixels_embolized.html')



# get t50 for each leaf
df_t50_whole_leaf = df_all_vulnerability %>%
  group_by(species, treatment, focal_distance_max, leaf_replicate) %>%
  mutate(frac_cumulative = n_pixels_cumulative / max(n_pixels_cumulative)) %>%
  mutate(at_50 = frac_cumulative >= 0.5) %>%
  filter(at_50) %>%
  arrange(time) %>%
  slice_head(n=1) %>%
  ungroup %>%
  mutate(species=factor(species, labels=c('A. rubrum (palmate)', 'C. erythrophyllum (pinnate)'))) %>%
  mutate(treatment=factor(treatment, labels=c('control','lamina','midrib'))) %>%
  mutate(time_hours = time*15/60)

g_t50_whole_leaf = ggplot(df_t50_whole_leaf %>% filter(focal_distance_max==200), aes(x=species,y=time_hours,color=treatment)) +
  geom_boxplot() +
  theme_bw() +
  scale_color_discrete(name='Treatment') +
  xlab('Species') +
  ylab('T50 (hours)')
ggsave(g_t50_whole_leaf, file='results/g_t50_whole_leaf.png',width=6,height=4)

m_t50_whole_leaf = lm(time_hours~species+treatment,data=df_t50_whole_leaf %>% filter(focal_distance_max==200))
tab_model(m_t50_whole_leaf, file='results/tab_m_t50_whole_leaf.html')



# start spatial analyses
files_spatial_focal = dir(dir_results, pattern='*focal\\.csv$')
files_spatial_random_big = dir(dir_results, pattern='*random_big\\.csv$')
files_spatial_random_small = dir(dir_results, pattern='*random_small\\.csv$')

process_file_spatial <- function(file_this) {
  
  file_parts = strsplit(gsub('\\.csv','',file_this), split='_')[[1]]

  df = read.csv(file.path(dir_results, file_this),header=FALSE)
  # pull out the normalization (by available lamina area)
  df_normalization = as.numeric(df[1,])
  # subset to the non-normalization part
  df = df[-1,,drop=FALSE]
  
  print(str(df[,1]))
  print(str(df_normalization))
  for (i in 1:ncol(df))
  {
    df[,i] = df[,i] / df_normalization[i]
  }
  
  names(df) = c(1:ncol(df))
  
  df_long = df %>% 
    apply(2, cumsum) %>%
    as.data.frame %>%
    mutate(time=row_number()) %>%
    pivot_longer(!time, names_to='simulation_replicate') %>%
    #arrange(name) %>%
    mutate(species=paste(file_parts[1:2],collapse='_'),
           treatment=file_parts[3],
           leaf_replicate=as.numeric(file_parts[4]),
           focal_distance_max=file_parts[5],
           simulation_treatment=paste(file_parts[6:length(file_parts)],collapse='_')) %>%
    mutate(simulation_replicate=as.numeric(simulation_replicate))
  
  return(df_long)
}

df_all = do.call('rbind',lapply(c(files_spatial_focal, files_spatial_random_big, files_spatial_random_small), 
                                process_file_spatial)) %>%
  group_by(species, treatment) %>%
  # renumber to 1-3
  mutate(leaf_replicate=as.numeric(factor(leaf_replicate))) %>%
  ungroup %>%
  as.data.frame

df_all = df_all %>% 
  group_by(species, treatment, leaf_replicate, simulation_replicate, simulation_treatment, focal_distance_max) %>% 
  mutate(value_normalized_by_px_embolized_max = value / max(value)) %>%
  mutate(time_hours = time * 15/60)



# calculate t50 metric
df_t50_normalized = df_all %>%
  ungroup %>%
  group_by(species, treatment, leaf_replicate, simulation_replicate, simulation_treatment, focal_distance_max) %>%
  mutate(at_50 = value_normalized_by_px_embolized_max>=0.5) %>%
  filter(at_50) %>%
  arrange(time) %>%
  slice_head(n=1) %>%
  ungroup %>%
  mutate(species=factor(species, labels=c('A. rubrum (palmate)', 'C. erythrophyllum (pinnate)'))) %>%
  mutate(treatment=factor(treatment, labels=c('control','lamina','midrib')))

g_t50_normalized = ggplot(df_t50_normalized %>% filter(focal_distance_max==200), aes(x=species,y=time_hours,color=simulation_treatment)) +
  geom_boxplot() +
  facet_grid(treatment~leaf_replicate,labeller=label_both) +
  theme_bw() +
  #ggtitle('t50 = value_normalized_by_px_embolized_max') +
  xlab('Species') +
  ylab('T50 normalized (hours)') +
  theme(legend.position='bottom') +
  scale_color_discrete(name='Simulation treatment')
ggsave(g_t50_normalized, file='results/g_t50_normalized_200.png',width=10,height = 7)

df_t01 = df_all %>%
  ungroup %>%
  group_by(species, treatment, leaf_replicate, simulation_replicate, simulation_treatment, focal_distance_max) %>%
  mutate(at_01 = value>0.01) %>%
  filter(at_01) %>%
  arrange(time) %>%
  slice_head(n=1)

g_t01_normalized = ggplot(df_t01 %>% filter(focal_distance_max==200), aes(x=species,y=time_hours,color=simulation_treatment)) +
  geom_boxplot() +
  facet_grid(treatment~leaf_replicate,labeller=label_both) +
  theme_bw() +
  xlab('Species') +
  ylab('T01 (hours)') +
  theme(legend.position='bottom') +
  scale_color_discrete(name='Simulation treatment')
ggsave(g_t01_normalized, file='results/g_t01_200.png',width=10,height = 7)




#library(lme4)
#summary(lmer(time ~ species*treatment*simulation_treatment - species:treatment:simulation_treatment + (1|species/leaf_replicate), data=df_t05 %>% filter(focal_distance_max==200)))
m_t01_200 = lm(time_hours ~ species+treatment+simulation_treatment + species*treatment + species*simulation_treatment, data=df_t01 %>% filter(focal_distance_max==200))
tab_model(m_t01_200, file='results/tab_m_t01_200.html')

m_t50_200 = lm(time_hours ~ species+treatment+simulation_treatment + species*treatment + species*simulation_treatment, data=df_t50_normalized %>% filter(focal_distance_max==200))
tab_model(m_t50_200, file='results/tab_m_t50_200.html')





focal_distances = as.numeric(unique(df_all$focal_distance_max))
lapply(focal_distances, function(fd) 
{
  g_ts = ggplot(df_all %>% filter(focal_distance_max==fd) %>%
                  ungroup %>%
                  mutate(species=factor(species, labels=c('A. rubrum (palmate)', 'C. erythrophyllum (pinnate)'))) %>%
                  mutate(treatment=factor(treatment, labels=c('control','lamina','midrib')))
                , aes(x=time_hours,y=value,
                     color=simulation_treatment, 
                     alpha=simulation_treatment %in% c("focal"),
                     linewidth=simulation_treatment %in% c("focal"),
                     group=paste(simulation_replicate,simulation_treatment))) + 
    geom_line() +
    theme_bw() +
    scale_color_brewer(palette='Accent', name='Simulation treatment') +
    scale_alpha_discrete(range=c(0.3,1),guide='none') +
    scale_linewidth_discrete(range=c(0.1,1),guide='none') +
    facet_grid(species~treatment + leaf_replicate, labeller=label_both) +
    theme(legend.position='bottom') +
    ylim(0,0.5) +
    xlab('Time (hours)') +
    ylab('Cumulative # pixels embolized / # available leaf pixels')
  ggsave(g_ts, file=sprintf('results/g_ts_focal_distance_max_%d.png', fd),width=13,height=7)
})

lapply(focal_distances, function(fd) 
{
  g_ts = ggplot(df_all %>% filter(focal_distance_max==fd) %>%
                  ungroup %>%
                  mutate(species=factor(species, labels=c('A. rubrum (palmate)', 'C. erythrophyllum (pinnate)'))) %>%
                  mutate(treatment=factor(treatment, labels=c('control','lamina','midrib')))
                , aes(x=time_hours,y=value_normalized_by_px_embolized_max,
                      color=simulation_treatment, 
                      alpha=simulation_treatment %in% c("focal"),
                      linewidth=simulation_treatment %in% c("focal"),
                      group=paste(simulation_replicate,simulation_treatment))) + 
    geom_line() +
    theme_bw() +
    scale_color_brewer(palette='Accent', name='Simulation treatment') +
    scale_alpha_discrete(range=c(0.3,1),guide='none') +
    scale_linewidth_discrete(range=c(0.1,1),guide='none') +
    facet_grid(species~treatment + leaf_replicate, labeller=label_both) +
    theme(legend.position='bottom') +
    ylim(0,1) +
    xlab('Time (hours)') +
    ylab('Cumulative # pixels embolized / Total # embolized pixels')
  ggsave(g_ts, file=sprintf('results/g_ts_focal_distance_max_normalized_%d.png', fd),width=13,height=7)
})






 

calculate_quantiles <- function(df, var)
{
  # cat('.')
  vals_focal = df %>%
    filter(simulation_treatment=='focal') %>%
    pull(!!var)
  #stopifnot(length(vals_focal) <= 1)
  vals_random_big = df %>%
    filter(simulation_treatment=='random_big') %>%
    pull(!!var)
  vals_random_small = df %>%
    filter(simulation_treatment=='random_small') %>%
    pull(!!var)
  #print(str(vals_focal))
  # print(length(vals_random_big))
  # print(length(vals_random_small))
  # p.small=NA
  # p.big=NA
  # try(p.small <- t.test(vals_random_small, paired=FALSE, mu=vals_focal)$p.value)
  # try(p.big <- t.test(vals_random_big, paired=FALSE, mu=vals_focal)$p.value)
  # 
  return(data.frame(
    # p.small=p.small,
    # p.big=p.big,
    z.small=(vals_focal - mean(vals_random_small))/sd(vals_random_small),
    z.big=(vals_focal - mean(vals_random_big))/sd(vals_random_big)))
  # return(data.frame(q.small=length(which(vals_focal < vals_random_small)) / length(vals_random_small),
  #                    q.big=length(which(vals_focal < vals_random_big)) / length(vals_random_big)))
}





df_quantiles = df_all %>% 
  #filter(species=='Ace_rub' & treatment=='L' & leaf_replicate==1) %>%
  group_by(species, treatment, leaf_replicate, time, focal_distance_max) %>%
  do(calculate_quantiles(.,'value')) %>%
  as.data.frame %>% 
  ungroup %>%
  mutate(species=factor(species, labels=c('A. rubrum (palmate)', 'C. erythrophyllum (pinnate)'))) %>%
  mutate(treatment=factor(treatment, labels=c('control','lamina','midrib')))

df_quantiles_normalized = df_all %>% 
  #filter(species=='Ace_rub' & treatment=='L' & leaf_replicate==1) %>%
  group_by(species, treatment, leaf_replicate, time, focal_distance_max) %>%
  do(calculate_quantiles(.,'value_normalized_by_px_embolized_max')) %>%
  as.data.frame %>% 
  ungroup %>%
  mutate(species=factor(species, labels=c('A. rubrum (palmate)', 'C. erythrophyllum (pinnate)'))) %>%
  mutate(treatment=factor(treatment, labels=c('control','lamina','midrib')))

focal_distances = as.numeric(unique(df_all$focal_distance_max))
lapply(focal_distances, function(fd) 
{
  g_quantiles = ggplot(df_quantiles %>% 
                         filter(focal_distance_max==fd) %>%
                         filter(treatment!='control') %>%
                         pivot_longer(!species:focal_distance_max) %>%
                         ungroup %>%
                         mutate(time_hours=time*15/60) %>%
                         mutate(name=factor(name,labels=c('relative to random_big','relative to random_small'))), 
                                aes(x=time,y=value,color=name,linetype=factor(leaf_replicate))) +
    geom_hline(yintercept = 0) +
    geom_line() + 
    facet_grid(species~treatment,labeller=label_both) +
    theme_bw() +
    geom_hline(yintercept=c(-1.96,1.96),color='gray') +
    scale_color_manual(values=brewer.pal(n=3,'Accent')[2:3],name='Type') +
    scale_linetype_discrete(name='Leaf replicate') +
    theme(legend.position='bottom') +
    ylab('z-score') +
    xlab('Time (hours)')
  ggsave(g_quantiles, file=sprintf('results/g_quantiles_focal_distance_max_%d.png',fd),width=7,height=6)
})

lapply(focal_distances, function(fd) 
{
  g_quantiles = ggplot(df_quantiles_normalized %>% 
                         filter(focal_distance_max==fd) %>%
                         filter(treatment!='control') %>%
                         pivot_longer(!species:focal_distance_max) %>%
                         ungroup %>%
                         mutate(time_hours=time*15/60) %>%
                         mutate(name=factor(name,labels=c('relative to random_big','relative to random_small'))), 
                       aes(x=time_hours,y=value,color=name,linetype=factor(leaf_replicate))) +
    geom_hline(yintercept = 0) +
    geom_line() + 
    facet_grid(species~treatment,labeller=label_both) +
    theme_bw() +
    geom_hline(yintercept=c(-1.96,1.96),color='gray') +
    scale_color_manual(values=brewer.pal(n=3,'Accent')[2:3],name='Type') +
    scale_linetype_discrete(name='Leaf replicate') +
    theme(legend.position='bottom') +
    ylab('z-score') +
    xlab('Time (hours)')
  ggsave(g_quantiles, file=sprintf('results/g_quantiles_focal_distance_max_normalized_%d.png',fd),width=7,height=6)
})
