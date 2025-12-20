library(ggplot2)
library(dplyr)
library(tidyr)
library(terra)
library(RStoolbox)
library(ggspatial)

dir_results = 'results'

files_images = dir(dir_results, pattern='*first_time\\.tif$')

process_file <- function(file_this) {
  
  file_parts = strsplit(gsub('first_time\\.tif','',file_this), split='_')[[1]]
  
  df = data.frame(
           species=paste(file_parts[1:2],collapse='_'),
           treatment=file_parts[3],
           focal_distance_max=file_parts[5],
           leaf_replicate=as.numeric(file_parts[4]))
  
  return(df)
}

df_all_images = do.call('rbind',lapply(files_images, process_file)) %>%
  group_by(species, treatment, focal_distance_max) %>%
  # renumber to 1-3
  mutate(leaf_replicate=as.numeric(factor(leaf_replicate))) %>%
  ungroup %>%
  as.data.frame %>%
  mutate(index=row_number()) %>%
  filter(focal_distance_max==200)

damage_loci = read.csv('data/damage_loci.csv')

df_all_images = df_all_images %>% left_join(damage_loci, by=c('species','treatment','leaf_replicate'))

list_images = lapply(file.path(dir_results, files_images), rast)

max_vals = lapply(list_images, function(x) { max(x[], na.rm=TRUE)  })
max_max_vals = max(unlist(max_vals)) * 15/60 # scale to hours

indices = df_all_images$index

for (i in 1:length(indices))
{
  cat(i)
  
  im_this = list_images[[ indices[i] ]]
  # mask out the NA values
  im_this[im_this==1] = NA
  # convert from time slice to miniutes
  im_this = im_this * 15/60
  title_this = paste(df_all_images$species[i], df_all_images$treatment[i], df_all_images$leaf_replicate[i])
  
  g_this = ggplot() + 
    ggR(im_this, geom_raster=TRUE, ggLayer = TRUE, stretch='none', maxpixels = 1e6) +
    theme_void() +
    theme(axis.title.x=element_blank()) +
    theme(axis.title.y=element_blank()) +
    coord_equal() +
    geom_point(data=data.frame(x=df_all_images$y[i]/2,y=df_all_images$x[i]/2), aes(x=x, y=y), color = "white", size=2,shape=3, inherit.aes = FALSE) +
    theme(legend.position='bottom') +
    scale_fill_viridis_c(limits=c(0,max_max_vals),name='Time (hours)',option=3, na.value = 'lightgray') +
    #ggtitle(title_this) +
    scale_x_continuous(breaks=c(1,dim(im_this)[2])) +
    scale_y_continuous(breaks=c(1,dim(im_this)[1]))
  rm(im_this)
  ggsave(g_this, file=sprintf('results/image_%s.png', title_this),width=7,height=5)
}
