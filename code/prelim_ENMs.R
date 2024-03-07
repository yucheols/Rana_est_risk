#####  load packages
library(megaSDM)
library(raster)
library(plyr)
library(dplyr)
library(readr)
library(SDMtune)
library(ggplot2)
library(rasterVis)
library(pals)

# prevent encoding error
Sys.getlocale()
Sys.setlocale("LC_CTYPE", ".1251")
Sys.getlocale()

##### Part 1 ::: envs data -----------------------------------------------------------------------------------------------------
### calibration range
# import China polygon
ch_poly <- rgdal::readOGR('E:/Asia_shp/China/CHN_adm0.shp')

# process envs
clim <- raster::stack(list.files(path = 'E:/env layers/worldclim', pattern = '.tif$', full.names = T))
clim <- raster::crop(clim, extent(ch_poly))
clim <- raster::mask(clim, ch_poly)

topo <- raster('E:/env layers/elev_worldclim/wc2.1_30s_elev.tif')
topo <- raster::crop(topo, extent(ch_poly))
topo <- raster::mask(topo, ch_poly)

land <- raster::stack(list.files(path = 'E:/env layers/land cover', pattern = '.tif$', full.names = T))
land <- raster::stack(subset(land, c('cultivated', 'water')))
land <- raster::crop(land, extent(ch_poly))
land <- raster::mask(land, ch_poly)

# stack env
envs <- raster::stack(clim, topo, land)
names(envs) = c('bio1', 'bio10', 'bio11', 'bio12', 'bio13', 'bio14', 'bio15', 'bio16', 'bio17', 'bio18', 'bio19', 
                'bio2', 'bio3', 'bio4', 'bio5', 'bio6', 'bio7', 'bio8', 'bio9', 'elev', 'cultivated', 'water')

print(envs)
plot(envs[[1]])


##### Part 2 ::: occs data -----------------------------------------------------------------------------------------------------
# acquire occs data
#megaSDM::OccurrenceCollection(spplist = c('Rana amurensis', 'Rana chensinensis', 'Rana coreana', 
#                                          'Rana dybowskii', 'Rana huanrenensis', 'Rana kukunoris'),
#                              output = 'data/occs',
#                              trainingarea = extent(ch_poly))

# compile occs
occs <- list.files(path = 'data/occs', pattern = '.csv', full.names = T) %>%
  lapply(read_csv) %>%
  rbind.fill %>%
  dplyr::select(4,6,5)

colnames(occs) = c('species', 'long', 'lat')
head(occs)

# sort per species
unique(sort(occs[, 1]))

amurensis <- occs %>% filter(species == 'Rana amurensis')
chensinensis <- occs %>% filter(species == 'Rana chensinensis')
coreana <- occs %>% filter(species == 'Rana coreana')
dybowskii <- occs %>% filter(species == 'Rana dybowskii')
huanrenensis <- occs %>% filter(species == 'Rana huanrenensis')
kukunoris <- occs %>% filter(species == 'Rana kukunoris')

### thin occurrence points
# automate
thinner <- function(occs_list, coords, env, x, y){
  output <- list()
  
  for (i in 1:length(occs_list)) {
    occs <- occs_list[[i]] %>% dplyr::select(2,3)
    thin <- thinData(coords = occs, env = terra::rast(env), x = x, y = y, verbose = T, progress = T)
    output[[i]] <- thin
  }
  
  return(output)
}

# run function
thin_occ <- thinner(occs_list = list(amurensis, chensinensis, coreana, dybowskii, huanrenensis, kukunoris),
                    coords = list(amurensis, chensinensis, coreana, dybowskii, huanrenensis, kukunoris),
                    env = envs, x = 'long', y = 'lat')

print(thin_occ)

# sort per species
amurensis <- thin_occ[[1]]
chensinensis <- thin_occ[[2]]
coreana <- thin_occ[[3]]    # only 2 occs.....drop this sp from modeling
dybowskii <- thin_occ[[4]]
huanrenensis <- thin_occ[[5]]
kukunoris <- thin_occ[[6]]

# plot out occs
plot(envs[[1]])

points(amurensis, col = 'red')
points(chensinensis, col = 'blue')
points(dybowskii, col = 'green')
points(huanrenensis, col = 'black')
points(kukunoris, col = 'orange')


##### Part 4 ::: background data -----------------------------------------------------------------------------------------------------
### automate buffer making
buffMaker <- function(occs_list, envs, buff_dist) {
  output <- list()
  
  for (i in 1:length(occs_list)) {
    occs.sf <- sf::st_as_sf(occs_list[[i]], coords = c('long', 'lat'), crs = raster::crs(envs))
    eq_area = '+proj=eck4 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs'
    occs.sf <- sf::st_transform(occs.sf, crs = eq_area)
    
    buff <- sf::st_buffer(occs.sf, dist = buff_dist) %>%
      sf::st_union() %>%
      sf::st_sf() %>%
      sf::st_transform(crs = raster::crs(envs))
    
    output[[i]] <- buff
  }
  
  return(output)
}

# run == 700 km buffers
buffers <- buffMaker(occs_list = list(amurensis, chensinensis, dybowskii, huanrenensis, kukunoris), 
                     envs = envs, buff_dist = 700000)

plot(envs[[1]])
plot(buffers[[4]], border = 'blue', lwd = 2, add = T)

### automate bg sampling 
bgSampler <- function(ex_bio, buffers, n, occs_list, excludep) {
  bg.out <- list()
  prog_bar = progress::progress_bar$new(total = length(buffers))
  
  for (i in 1:length(buffers)) {
    prog_bar$tick()
    
    mask <- raster::mask(ex_bio, buffers[[i]])
    bg <- dismo::randomPoints(mask = mask, n = n, p = occs_list[[i]], excludep = excludep)
    colnames(bg) = c('long', 'lat')
    
    bg.out[[i]] <- bg
  }
  return(bg.out)
}


# sample background points
bg <- bgSampler(ex_bio = envs[[1]], buffers = buffers, n = 10000, 
                occs_list = list(amurensis, chensinensis, dybowskii, huanrenensis, kukunoris), excludep = T)

# lets see if this is done correctly
plot(envs[[1]])

plot(buffers[[1]], add = T, lwd = 2, border = 'blue')
points(bg[[1]])


##### Part 5 ::: select envs data -----------------------------------------------------------------------------------------------------
# use ENMTools to sort out highly correlated rasters
cor.mat <- ENMTools::raster.cor.matrix(env = envs, method = 'pearson')
print(cor.mat)

find.cor <- caret::findCorrelation(data.matrix(cor.mat), cutoff = abs(0.7))
envs <- raster::dropLayer(envs, sort(find.cor))
print(envs)

##### Part 6 ::: model fitting -----------------------------------------------------------------------------------------------------
### make SWD
swd.maker <- function(occs.list, bg.list, spp.list, env, categ) {
  out.swd <- list()
  
  for (i in 1:length(occs.list)) {
    swd <- SDMtune::prepareSWD(species = spp.list[[i]], env = terra::rast(env), 
                               p = occs.list[[i]], a = bg.list[[i]], categorical = categ, verbose = T)
    
    out.swd[[i]] <- swd
    print('==============   Done!   ==============')
  }
  return(out.swd)
}

# make SWD per species
swd <- swd.maker(occs.list = list(amurensis, chensinensis, dybowskii, huanrenensis, kukunoris),
                 bg.list = bg, spp.list = c('amurensis', 'chensinensis', 'dybowskii', 'huanrenensis', 'kukunoris'),
                 env = envs, categ = NULL)

### generate CV folds
folds <- function(occs.list, bg.list, kfolds) {
  folds.out <- list()
  
  for (i in 1:length(occs.list)) {
    folds <- ENMeval::get.randomkfold(occs = occs.list[[i]], bg = bg.list[[i]], kfolds = kfolds)
    folds.out[[i]] <- folds
  }
  return(folds.out)
}

# get folds
cv.folds <- folds(occs.list = list(amurensis, chensinensis, dybowskii, huanrenensis, kukunoris), bg.list = bg, kfolds = 10)
print(cv.folds)


### generate default models
# modeling automation
auto.model <- function(swd.list, method, folds, progress, iter, type) {
  models <- list()
  
  for (i in 1:length(swd.list)) {
    model <- SDMtune::train(method = method, data = swd.list[[i]], folds = folds[[i]], progress = progress, iter = iter, type = type)
    models[[i]] <- model
  }
  return(models)
}

# make models
rana.enms <- auto.model(swd.list = swd, method = 'Maxent', folds = cv.folds, progress = T, iter = 5000, type = 'cloglog')
print(rana.enms)

# evaluate models
for (i in 1:length(rana.enms)) {
  print(auc(model = rana.enms[[i]], test = T))
}

##### Part 6 ::: model prediction -----------------------------------------------------------------------------------------------------
# automate predictions
model.pred <- function(models, data, type, clamp, progress) {
  pred.out <- list()
  
  for (i in 1:length(models)) {
    pred <- SDMtune::predict(object = models[[i]], data = terra::rast(data), 
                             type = type, clamp = clamp, progress = progress) %>% raster::raster()
    pred.out[[i]] <- pred 
  }
  return(pred.out)
}

# make pred
preds <- model.pred(models = rana.enms, data = envs, type = 'cloglog', clamp = T, progress = T)
preds <- raster::stack(preds)

plot(preds[[1]])
plot(preds[[2]])

# plot models

##### Part 7 ::: response curves -----------------------------------------------------------------------------------------------------





##### Part 8 ::: model projections -----------------------------------------------------------------------------------------------------
# ROK polygon
rok <- rgdal::readOGR('E:/Asia_shp/South Korea/KOR_adm0.shp')

# prep projections layers
clim2 <- raster::stack(list.files(path = 'E:/env layers/worldclim', pattern = '.tif$', full.names = T))
clim2 <- raster::crop(clim2, extent(rok))
clim2 <- raster::mask(clim2, rok)

topo2 <- raster('E:/env layers/elev_worldclim/wc2.1_30s_elev.tif')
topo2 <- raster::crop(topo2, extent(rok))
topo2 <- raster::mask(topo2, rok)

land2 <- raster::stack(list.files(path = 'E:/env layers/land cover', pattern = '.tif$', full.names = T))
land2 <- raster::stack(subset(land2, c('cultivated', 'water')))
land2 <- raster::crop(land2, extent(rok))
land2 <- raster::mask(land2, rok)

# stack
envs2 <- raster::stack(clim2, topo2, land2)
names(envs2) = c('bio1', 'bio10', 'bio11', 'bio12', 'bio13', 'bio14', 'bio15', 'bio16', 'bio17', 'bio18', 'bio19', 
                'bio2', 'bio3', 'bio4', 'bio5', 'bio6', 'bio7', 'bio8', 'bio9', 'elev', 'cultivated', 'water')

plot(envs2[[1]])

# make predictions to ROK
proj <- model.pred(models = rana.enms, data = envs2, type = 'cloglog', clamp = T, progress = T)
proj.out <- raster::stack(proj)
plot(proj.out)                       # model order == amurensis, chensinensis, dybowskii, huanrenensis, kukunoris

names(proj.out) = c('R.amurensis', 'R.chensinensis', 'R.dybowskii', 'R.huanrenensis', 'R.kukunoris')


##### Part 9 ::: plot outputs -----------------------------------------------------------------------------------------------------
# plot
gplot(proj.out) +
  geom_tile(aes(fill = value)) +
  facet_wrap(~ variable) +
  coord_equal() +
  scale_fill_gradientn(colors = rev(as.vector(ocean.thermal(1000))),
                       na.value = NA,
                       name = 'Suitability') +
  xlab('Long') + ylab('Lat') +
  theme_bw() +
  theme(strip.text = element_text(face = 'italic', size = 14),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12)) 


# save
ggsave('output/plot/prelim_output.png', width = 30, height = 22, dpi = 600, units = 'cm')

