#####  use ENMwrap to generate models  ----------

# load packages
library(ENMwrap)
library(megaSDM)
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

###