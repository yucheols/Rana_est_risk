#####  use ENMwrap and other packages to generate models  ----------

# load packages
library(ENMwrap)
library(megaSDM)
library(SDMtune)
library(raster)
library(plyr)
library(dplyr)
library(readr)

# clean working environment
rm(list = ls(all.names = T))
gc()


#####  part 1 ::: env data ----------

### prep calibration range
# import China polygon
ch_poly <- rgdal::readOGR('E:/Asia_shp/China/CHN_adm0.shp')

# climate 
clim <- raster::stack(list.files(path = 'E:/env layers/worldclim', pattern = '.tif$', full.names = T))
clim <- raster::crop(clim, extent(ch_poly))
clim <- raster::mask(clim, ch_poly)

# topo
topo <- raster('E:/env layers/elev_worldclim/wc2.1_30s_elev.tif')
topo <- raster::crop(topo, extent(ch_poly))
topo <- raster::mask(topo, ch_poly)

# land cover
land <- raster::stack(list.files(path = 'E:/env layers/land cover', pattern = '.tif$', full.names = T))
land <- raster::stack(subset(land , c('cultivated', 'water')))
land <- raster::crop(land, extent(ch_poly))
land <- raster::mask(land, ch_poly)

# stack them
envs <- raster::stack(clim, topo, land)
print(envs)
plot(envs[[1]])

# rename layers
names(envs) = c('bio1','bio10','bio11','bio12','bio13','bio14','bio15','bio16','bio17','bio18','bio19',
                'bio2','bio3','bio4','bio5','bio6','bio7','bio8','bio9','elev','cultivated','water')

# export processed calibration layers
for (i in 1:nlayers(envs)) {
  r <- envs[[i]]
  file_name <- paste0('data/envs/calibration/', names(envs)[i], '.bil')
  writeRaster(r, filename = file_name, overwrite = T)
}

# import shortcut 
envs <- raster::stack(list.files(path = 'data/envs/calibration', pattern = '.bil', full.names = T))

#####  part 2 ::: occurrence data ----------
# acquire occs data
megaSDM::OccurrenceCollection(spplist = c('Rana amurensis', 'Rana chensinensis', 'Rana coreana', 
                                          'Rana dybowskii', 'Rana huanrenensis', 'Rana kukunoris'),
                              output = 'data/occs',
                              trainingarea = extent(ch_poly))

# compile occurrence points
occs <- list.files(path = 'data/occs', pattern = '.csv', full.names = T) %>%
  lapply(read_csv) %>%
  rbind.fill() %>%
  dplyr::select(4,6,5)

colnames(occs) = c('species', 'long', 'lat')
head(occs)

# check the list of species
unique(occs$species)

# create a list of species
occs.list <- list(occs %>% dplyr::filter(species == unique(occs$species)[1]) %>% dplyr::select(2,3),
                  occs %>% dplyr::filter(species == unique(occs$species)[2]) %>% dplyr::select(2,3),
                  occs %>% dplyr::filter(species == unique(occs$species)[3]) %>% dplyr::select(2,3),
                  occs %>% dplyr::filter(species == unique(occs$species)[4]) %>% dplyr::select(2,3),
                  occs %>% dplyr::filter(species == unique(occs$species)[5]) %>% dplyr::select(2,3),
                  occs %>% dplyr::filter(species == unique(occs$species)[6]) %>% dplyr::select(2,3))

# thin occurrence points
occs.thin <- occs_thinner(occs_list = occs.list, envs = envs[[1]], long = 'long', lat = 'lat', 
                          spp_list = list('Rana amurensis', 'Rana chensinensis', 'Rana coreana', 'Rana dybowskii', 'Rana huanrenensis', 'Rana kukunoris'))


#####  part 3 ::: background data ----------

# make buffers per species == 700 km circular buffers
buff <- buffMaker(occs_list = occs.thin, envs = envs, buff_dist = 700000)

# sample bg
bg.list <- bgSampler(envs = envs[[1]], n = 10000, occs_list = occs.thin, buffer_list = buff, excludep = T, method = 'buffer')

# export bg
for (i in 1:length(bg.list)) {
  file <- bg.list[[i]]
  write.csv(file, paste0('data/bg/', unique(occs$species)[i], '.csv'))
}

#####  part 4 ::: select environmental data ----------

# sample 10000 for layer selection
bg.env.sel <- dismo::randomPoints(mask = envs[[1]], n = 10000) %>% as.data.frame()
write.csv(bg.env.sel, 'data/bg/bg_env_sel.csv')

# extract env values from bg
env_val <- raster::extract(envs, bg.env.sel)
head(env_val)

# generate correlation matrix
env_cor <- cor(env_val)
print(env_cor)

# run ntbox
ntbox::correlation_finder(cor_mat = env_cor, threshold = 0.7, verbose = T)

# subset environmental data 
envs.sub <- raster::stack(subset(envs, c('bio1','bio3','bio4','bio12','bio15','cultivated','water')))


#####  part 4 ::: prep data for candidate model testing ----------

# use SDMtune for speed

# function for automated SWD prep
swd.maker <- function(occs.list, bg.list, spp.list, env, categ = NULL) {
  out.swd <- list()
  
  for (i in 1:length(occs.list)) {
    swd <- SDMtune::prepareSWD(species = spp.list[[i]], env = terra::rast(env), 
                               p = occs.list[[i]], a = bg.list[[i]], categorical = categ, verbose = T)
    
    out.swd[[i]] <- swd
    print('==============   Done!   ==============')
  }
  return(out.swd)
}

# generate SWD
swd <- swd.maker(occs.list = occs.thin, bg.list = bg.list, 
                 spp.list = list('Rana amurensis', 'Rana chensinensis', 'Rana coreana', 'Rana dybowskii', 'Rana huanrenensis', 'Rana kukunoris'),
                 env = envs.sub)

print(swd)


# function to automate block CV folds generation
block.folds <- function(occs.list, bg.list, orientation) {
  folds.out <- list()
  
  for (i in 1:length(occs.list)) {
    folds <- ENMeval::get.block(occs = occs.list[[i]], bg = bg.list[[i]], orientation = orientation)
    folds.out[[i]] <- folds
  }
  return(folds.out)
}

# generate folds
cvfolds <- block.folds(occs.list = occs.thin, bg.list = bg.list, orientation = 'lat_lon')
print(cvfolds)

# automate the generation of default models
auto.model <- function(swd.list, method, folds, progress, iter, type) {
  models <- list()
  
  for (i in 1:length(swd.list)) {
    model <- SDMtune::train(method = method, data = swd.list[[i]], folds = folds[[i]], progress = progress, iter = iter, type = type)
    models[[i]] <- model
  }
  return(models)
}

# generate default models
def.mods <- auto.model(swd.list = swd, method = 'Maxent', folds = cvfolds, progress = T, iter = 5000, type = 'cloglog')


#####  part 5 ::: model tuning ----------

# function for automation
model.tune <- function(list.models, hypers, metric, test = NULL, env = NULL, save_models, interactive, progress) {
  tuned.models <- list()
  results <- list()
  
  for (i in 1:length(list.models)) {
    tune <- gridSearch(model = list.models[[i]], hypers = hypers, metric = metric, test = test, env = env, save_models = save_models, 
                       interactive = interactive, progress = progress)
    
    results[[i]] <- tune@results
    tuned.models[[i]] <- tune@models
  }
  return(list(results = results, tuned.models = tuned.models))
}

# run model tuning
test.models <- model.tune(list.models = def.mods, 
                          hypers = list(fc = c('l', 'lq', 'h', 'lqh', 'lqhp', 'lqhpt'),
                                        reg = seq(1,5, by = 0.5)),
                          metric = 'auc',
                          save_models = T,
                          interactive = F,
                          progress = T)
