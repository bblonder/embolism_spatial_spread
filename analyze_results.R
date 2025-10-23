library(ggplot2)
library(dplyr)
library(tidyr)

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
           leaf_replicate=as.numeric(file_parts[4])) %>%
    select(-name)
  
  return(df_long)
}

df_all_vulnerability = do.call('rbind',lapply(files_vulnerability, process_file_vulnerability)) %>%
  group_by(species, treatment) %>%
  # renumber to 1-3
  mutate(leaf_replicate=as.numeric(factor(leaf_replicate))) %>%
  ungroup %>%
  as.data.frame
# do normalization
df_all_vulnerability = df_all_vulnerability %>%
  ungroup %>%
  group_by(species, treatment, leaf_replicate) %>%
  mutate(frac_pixels_cumulative = n_pixels_cumulative / max(n_pixels_cumulative)) %>%
  mutate(frac_time = time/max(time))


ggplot(df_all_vulnerability, aes(x=time,y=n_pixels_cumulative)) +
  geom_line() +
  theme_bw() + 
  facet_grid(species~treatment + leaf_replicate) +
  theme(legend.position='bottom')

g_frac_pixels_cumulative = ggplot(df_all_vulnerability, aes(x=time,y=frac_pixels_cumulative,
                                                            linetype=factor(leaf_replicate),
                                                            color=treatment)) +
  geom_line() +
  theme_bw() + 
  facet_grid(~species) +
  theme(legend.position='bottom')
ggsave(g_frac_pixels_cumulative, file='results/g_frac_pixels_cumulative.pdf',width=8,height=6)




files_spatial_focal = dir(dir_results, pattern='*focal\\.csv$')
files_spatial_random_big = dir(dir_results, pattern='*random_big\\.csv$')
files_spatial_random_small = dir(dir_results, pattern='*random_small\\.csv$')

process_file_spatial <- function(file_this) {
  
  file_parts = strsplit(gsub('\\.csv','',file_this), split='_')[[1]]

  df = read.csv(file.path(dir_results, file_this),header=FALSE)
  # pull out the normalization
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
           simulation_treatment=paste(file_parts[5:length(file_parts)],collapse='_')) %>%
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

g_ts = ggplot(df_all, aes(x=time,y=value,
                   color=simulation_treatment, 
                   alpha=simulation_treatment %in% c("focal"),
                   linewidth=simulation_treatment %in% c("focal"),
                   group=paste(simulation_replicate,simulation_treatment))) + 
  geom_line() +
  theme_bw() +
  scale_alpha_discrete(range=c(0.25,1),guide='none') +
  scale_linewidth_discrete(range=c(0.1,0.5),guide='none') +
  facet_grid(species~treatment + leaf_replicate) +
  theme(legend.position='bottom')
ggsave(g_ts, file='results/g_ts.pdf',width=13,height=7)



 

calculate_quantiles <- function(df)
{
  # cat('.')
  vals_focal = df %>%
    filter(simulation_treatment=='focal') %>%
    pull(value)
  #stopifnot(length(vals_focal) <= 1)
  vals_random_big = df %>%
    filter(simulation_treatment=='random_big') %>%
    pull(value)
  vals_random_small = df %>%
    filter(simulation_treatment=='random_small') %>%
    pull(value)
  #print(str(vals_focal))
  # print(length(vals_random_big))
  # print(length(vals_random_small))
  p.small=NA
  p.big=NA
  try(p.small <- t.test(vals_random_small, paired=FALSE, mu=vals_focal)$p.value)
  try(p.big <- t.test(vals_random_big, paired=FALSE, mu=vals_focal)$p.value)
  
  return(data.frame(
    p.small=p.small,
    p.big=p.big,
    z.small=(vals_focal - mean(vals_random_small))/sd(vals_random_small),
    z.big=(vals_focal - mean(vals_random_big))/sd(vals_random_big)))
  # return(data.frame(q.small=length(which(vals_focal < vals_random_small)) / length(vals_random_small),
  #                    q.big=length(which(vals_focal < vals_random_big)) / length(vals_random_big)))
}

df_quantiles = df_all %>% 
  #filter(species=='Ace_rub' & treatment=='L' & leaf_replicate==1) %>%
  group_by(species, treatment, leaf_replicate, time) %>%
  do(calculate_quantiles(.)) %>%
  as.data.frame

g_quantiles = ggplot(df_quantiles %>% 
         pivot_longer(!species:time) %>%
         filter(treatment!='C' & name %in% c('z.big','z.small')), aes(x=time,y=value,color=name,linetype=factor(leaf_replicate))) +
  geom_hline(yintercept = 0) +
  geom_line() + 
  facet_grid(species~treatment) +
  theme_bw() +
  geom_hline(yintercept=c(-2,2),color='green')
ggsave(g_quantiles, file='results/g_quantiles.png',width=8,height=8)


# not very useful below
df_quantiles_categorized = df_quantiles %>%
  mutate(outcome.small = factor(
           ifelse(p.small < 0.05, ifelse(z.small < 0, 'negative','positive'),'insignificant'),levels=c('negative','insignificant','positive')))

ggplot(df_quantiles_categorized %>% 
         filter(treatment!='C') %>%
         filter(!is.na(outcome.small)), 
       aes(x=time,y=outcome.small, color=factor(leaf_replicate),
           group=paste(species,treatment,leaf_replicate))) + 
  geom_line() + 
  facet_grid(species~treatment) +
  theme_bw()
