##### Korean native Rana species ENMs
# load packages
library(raster)
library(megaSDM)
library(plyr)
library(dplyr)
library(readr)
library(SDMtune)

# prevent encoding error
Sys.getlocale()
Sys.setlocale("LC_CTYPE", ".1251")
Sys.getlocale()

##### Part 1 ::: env data prep -----------------------------------------------------------------------------------------------------
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

# subset rasters
print(find.cor)

envs2 <- dropLayer(envs2, sort(find.cor))
print(envs2)
plot(envs2[[1]])


##### Part 2 ::: get occurrence data -----------------------------------------------------------------------------------------------------
OccurrenceCollection(spplist = c('Rana uenoi', 'Rana coreana', 'Rana huanrenensis'),
                     output = 'data/occs/Korea',
                     trainingarea = extent(envs2))

# compile occurrence points
rok.occs <- list.files(path = 'data/occs/Korea', pattern = '.csv', full.names = T) %>%
  lapply(read_csv) %>%
  rbind.fill %>%
  dplyr::select(4,6,5)

colnames(rok.occs) = c('species', 'long', 'lat')
head(rok.occs)

# sort occurrence points
coreana.rok <- rok.occs %>% filter(species == 'Rana coreana')
huanrenensis.rok <- rok.occs %>% filter(species == 'Rana huanrenensis')
uenoi.rok <- rok.occs %>% filter(species == 'Rana uenoi')

# thin occurrence data
coreana.rok <- thinData(coords = coreana.rok, env = terra:: rast(envs2[[1]]), x = 'long', y = 'lat', verbose = T, progress = T)
huanrenensis.rok <- thinData(coords = huanrenensis.rok, env = terra::rast(envs2[[1]]), x = 'long', y = 'lat', verbose = T, progress = T)
uenoi.rok <- thinData(coords = uenoi.rok, env = terra::rast(envs2[[1]]), x = 'long', y = 'lat', verbose = T, progress = T)

##### Part 3 ::: sample background data -----------------------------------------------------------------------------------------------------
rok.bg <- dismo::randomPoints(mask = envs2[[1]], n = 10000, p = rok.occs[, c(2,3)], excludep = T) %>% as.data.frame()
colnames(rok.bg) = c('long', 'lat')
plot(envs2[[1]])
points(rok.bg)

