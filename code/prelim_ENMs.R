#####  load packages
library(megaSDM)
library(raster)
library()

# prevent encoding error
Sys.getlocale()
Sys.setlocale("LC_CTYPE", ".1251")
Sys.getlocale()

# Part 1 ::: data processing --------------------------------------------------------------------------------------------------------------
# import China polygon
ch_poly <- rgdal::readOGR('E:/Asia_shp/Chinas/CHN_adm0.shp')

# acquire occs data
megaSDM::OccurrenceCollection(spplist = c('Rana amurensis', 'Rana chensinensis', 'Rana coreana', 
                                          'Rana dybowskii', 'Rana huanrenensis', 'Rana kukunoris'),
                              output = 'data/occs',
                              trainingarea = extent(ch_poly))

# compile occs
