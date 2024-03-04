#####  load packages
library(megaSDM)
library(raster)
library(plyr)
library(dplyr)
library(readr)
library(SDMtune)

# prevent encoding error
Sys.getlocale()
Sys.setlocale("LC_CTYPE", ".1251")
Sys.getlocale()

##### Part 1 ::: envs data -----------------------------------------------------------------------------------------------------
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
# import China polygon
ch_poly <- rgdal::readOGR('E:/Asia_shp/Chinas/CHN_adm0.shp')

# acquire occs data
megaSDM::OccurrenceCollection(spplist = c('Rana amurensis', 'Rana chensinensis', 'Rana coreana', 
                                          'Rana dybowskii', 'Rana huanrenensis', 'Rana kukunoris'),
                              output = 'data/occs',
                              trainingarea = extent(ch_poly))

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

# run
buffers <- buffMaker(occs_list = list(amurensis, chensinensis, coreana, dybowskii, huanrenensis, kukunoris), 
                     envs = envs, buff_dist = 300000)

plot(envs[[1]])
plot(buffers[[1]], border = 'blue', lwd = 2, add = T)

### sample background points
bgSampler <- function() {
  
}



