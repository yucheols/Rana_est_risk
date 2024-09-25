#####  use ENMwrap and other packages to generate models  ----------

# set random seed
set.seed(333)

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

# import shortcut == calibration range
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
tail(occs)

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
buff <- buff_maker(occs_list = occs.thin, envs = envs, buff_dist = 700000)

# sample bg
bg.list <- bg_sampler(envs = envs[[1]], n = 10000, occs_list = occs.thin, buffer_list = buff, excludep = T, method = 'buffer')

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

# check results
glimpse(test.models$results)
glimpse(test.models$tuned.models)

# export model object
saveRDS(test.models, 'output/models/model_test_output.rds')

# function to retrieve optimal parameter combinations
get.opt <- function(list.results) {
  output <- list()
  
  for (i in 1:length(list.results)) {
    opt <- list.results[[i]] %>% dplyr::filter(diff_AUC == min(diff_AUC)) %>%
      dplyr::filter(test_AUC == max(test_AUC))
    
    output[[i]] <- opt
  }
  return(output)
}

# get optimal parameter combinations
opt.params <- get.opt(list.results = test.models$results)
print(opt.params)


#' Rana amurensis      == lqhpt 1.0 == 6th model
#' Rana chensinensis   == lq 1.5    == 8th model
#' Rana coreana        == discard   
#' Rana dybowskii      == h 1.0     == 3rd model
#' Rana huanrenensis   == lqh 1.0   == 4th model
#' Rana kukunoris      == h 2.5     == 21st model

# assign models to new objects
amurensis_model <- test.models$tuned.models[[1]][[6]]
chensinensis_model <- test.models$tuned.models[[2]][[8]]
dybowskii_model <- test.models$tuned.models[[4]][[3]]
huanrenensis_model <- test.models$tuned.models[[5]][[4]]
kukunoris_model <- test.models$tuned.models[[6]][[21]]


#####  part 6 ::: variable importance ----------

# get variable importance
amurensis_varimp <- maxentVarImp(amurensis_model)
chensinensis_varimp <- maxentVarImp(chensinensis_model)
dybowskii_varimp <- maxentVarImp(dybowskii_model)
huanrenensis_varimp <- maxentVarImp(huanrenensis_model)
kukunoris_varimp <- maxentVarImp(kukunoris_model)

# group into list
varimp_list <- list(amurensis_varimp, chensinensis_varimp, dybowskii_varimp, huanrenensis_varimp, kukunoris_varimp)
sp.names = list('R.amurensis', 'R.chensinensis', 'R.dybowskii', 'R.huanrenensis', 'R.kukunoris')

# export
for (i in 1:length(varimp_list)) {
  file <- varimp_list[[i]]
  write.csv(file, paste0('output/varimp/', sp.names[[i]], '_varimp.csv'))
}

#####  part 7 ::: response curve ----------

###' function to pull out response data from SDMtune model outputs
respDataPull <- function(model, var, type, only_presence, marginal, species_name) {
  
  plotdata.list <- list()
  
  for (i in 1:length(var)) {
    plotdata <- plotResponse(model = model, var = var[[i]], type = type, only_presence = only_presence, marginal = marginal)
    plotdata <- ggplot2::ggplot_build(plotdata)$data
    plotdata <- plotdata[[1]]
    
    plotdata <- plotdata[, c(1:4)]
    plotdata$species <- species_name
    plotdata$var <- var[[i]]
    
    plotdata.list[[i]] <- plotdata
  }
  plotdata.df <- dplyr::bind_rows(plotdata.list) 
  return(plotdata.df)
}

###' amurensis
# get data
amurensis_resp <- respDataPull(species_name = 'R.amurensis', model = amurensis_model, var = names(envs.sub), 
                               type = 'cloglog', only_presence = F, marginal = F)

# plot
plot_response(amurensis_resp) +
  theme(axis.title = element_text(size = 16, face = 'bold'),
        axis.title.x = element_text(margin = margin(t = 20)),
        axis.title.y = element_text(margin = margin(r = 20)),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 16))

# export
ggsave('output/plot/amurensis_resp.jpg', width = 30, height = 22, dpi = 600, units = 'cm')


###' chensinensis
# get data
chensinensis_resp <- respDataPull(species_name = 'R.chensinensis', model = chensinensis_model, var = names(envs.sub),
                                  type = 'cloglog', only_presence = F, marginal = F)

# plot
plot_response(chensinensis_resp) +
  theme(axis.title = element_text(size = 16, face = 'bold'),
        axis.title.x = element_text(margin = margin(t = 20)),
        axis.title.y = element_text(margin = margin(r = 20)),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 16))

# export
ggsave('output/plot/chensinensis_resp.jpg', width = 30, height = 22, dpi = 600, units = 'cm')


###' dybowskii
# get data
dybowskii_resp <- respDataPull(species_name = 'R.dybowskii', model = dybowskii_model, var = names(envs.sub),
                               type = 'cloglog', only_presence = F, marginal = F)

# plot
plot_response(dybowskii_resp) +
  theme(axis.title = element_text(size = 16, face = 'bold'),
        axis.title.x = element_text(margin = margin(t = 20)),
        axis.title.y = element_text(margin = margin(r = 20)),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 16))

# export
ggsave('output/plot/dybowskii_resp.jpg', width = 30, height = 22, dpi = 600, units = 'cm')


###' huanrenensis 
# get data
huanrenensis_resp <- respDataPull(species_name = 'R.huanrenensis', model = huanrenensis_model, var = names(envs.sub),
                                  type = 'cloglog', only_presence = F, marginal = F)

# plot
plot_response(huanrenensis_resp) +
  theme(axis.title = element_text(size = 16, face = 'bold'),
        axis.title.x = element_text(margin = margin(t = 20)),
        axis.title.y = element_text(margin = margin(r = 20)),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 16))

# export
ggsave('output/plot/huanrenensis_resp.jpg', width = 30, height = 22, dpi = 600, units = 'cm')


###' kukunoris
# get data
kukunoris_resp <- respDataPull(species_name = 'R.kukunoris', model = kukunoris_model, var = names(envs.sub),
                               type = 'cloglog', only_presence = F, marginal = F)

# plot
plot_response(kukunoris_resp) +
  theme(axis.title = element_text(size = 16, face = 'bold'),
        axis.title.x = element_text(margin = margin(t = 20)),
        axis.title.y = element_text(margin = margin(r = 20)),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 16))

# export
ggsave('output/plot/kukunoris_resp.jpg', width = 30, height = 22, dpi = 600, units = 'cm')


#####  part 8 ::: model thresholds ----------

# automate
getTh <- function(model.list) {
  thresh.out <- list()
  
  for (i in 1:length(model.list)) {
    th <- SDMtune::maxentTh(combineCV(model = model.list[[i]]))
    thresh.out[[i]] <- th
  }
  return(thresh.out)
}

# get thresholds
model.th <- getTh(model.list = list(amurensis_model, chensinensis_model, dybowskii_model, huanrenensis_model, kukunoris_model))
print(model.th)

# export
for (i in 1:length(model.th)) {
  file <- model.th[[i]]
  write.csv(model.th, paste0('output/threshold/', sp.names[[i]], '_thresh.csv'))
}

#####  part 9 ::: model prediction ----------

# put models into a new list
opt.models <- list(models = c(amurensis_model, chensinensis_model, dybowskii_model, huanrenensis_model, kukunoris_model),
                   sp.names = c('R.amurensis', 'R.chensinensis', 'R.dybowskii', 'R.huanrenensis', 'R.kukunoris'))

# automate predictions
model.pred2 <- function(models, sp.names, data, type, clamp, progress) {
  pred.out <- list()
  
  for (i in 1:length(models)) {
    pred <- SDMtune::predict(object = models[[i]], data = terra::rast(data), 
                             type = type, clamp = clamp, progress = progress) %>% raster::raster()
    pred.out[[i]] <- pred
  }
  pred.out <- raster::stack(pred.out)
  names(pred.out) = sp.names
  return(pred.out)
}

# make predictions
opt.pred <- model.pred2(models = opt.models$models, sp.names = opt.models$sp.names, data = envs.sub, type = 'cloglog', clamp = T, progress = T)
print(opt.pred)

# plot predictions
plot_preds(preds = opt.pred, poly = ch_poly, pred.names = names(opt.pred), colors = rev(as.vector(pals::ocean.thermal(1000))))


#####  part 10 ::: model projection to R Korea ----------

### prep projection range
# import R Korea polygon
kr_poly <- rgdal::readOGR('E:/Asia_shp/South Korea/KOR_adm0.shp')

# climate 
clim <- raster::stack(list.files(path = 'E:/env layers/worldclim', pattern = '.tif$', full.names = T))
clim <- raster::crop(clim, extent(kr_poly))
clim <- raster::mask(clim, kr_poly)

# topo
topo <- raster('E:/env layers/elev_worldclim/wc2.1_30s_elev.tif')
topo <- raster::crop(topo, extent(kr_poly))
topo <- raster::mask(topo, kr_poly)

# land cover
land <- raster::stack(list.files(path = 'E:/env layers/land cover', pattern = '.tif$', full.names = T))
land <- raster::stack(subset(land , c('cultivated', 'water')))
land <- raster::crop(land, extent(kr_poly))
land <- raster::mask(land, kr_poly)

# stack them
envs.prj <- raster::stack(clim, topo, land)
print(envs.prj)
plot(envs.prj[[1]])

# rename layers
names(envs.prj) = c('bio1','bio10','bio11','bio12','bio13','bio14','bio15','bio16','bio17','bio18','bio19',
                    'bio2','bio3','bio4','bio5','bio6','bio7','bio8','bio9','elev','cultivated','water')

# export processed projection layers
for (i in 1:nlayers(envs.prj)) {
  r <- envs.prj[[i]]
  file_name <- paste0('data/envs/projection/', names(envs.prj)[i], '.bil')
  writeRaster(r, filename = file_name, overwrite = T)
}

# subset the projection layers to needed layers
envs.prj <- raster::stack(list.files(path = 'data/envs/projection', pattern = '.bil', full.names = T))
envs.prj <- raster::stack(subset(envs.prj, c('bio1','bio3','bio4','bio12','bio15','cultivated','water')))

print(envs.prj)
plot(envs.prj[[1]])

# model projection
kor.pred <- model.pred2(models = opt.models$models, sp.names = opt.models$sp.names, data = envs.prj, type = 'cloglog', clamp = T, progress = T)
print(kor.pred)                                                    

# plot projected models
plot_preds(preds = kor.pred, poly = kr_poly, colors = rev(as.vector(pals::ocean.thermal(1000))), pred.names = names(kor.pred), ncol = 5, nrow = 1) 
ggsave('output/plot/result_output.png', width = 35, height = 8, dpi = 600, units = 'cm')


#####  part 11 ::: binary map ----------
print(model.th)

# reformat threshold data
th.reform <- function(th.list) {
  th.reform.out <- list()
  
  for (i in 1:length(th.list)) {
    th.df <- data.frame(th.list[[i]]$threshold)
    colnames(th.df) = 'threshold'
    th.reform.out[[i]] <- th.df
  }
  return(th.reform.out)
}

# reformat
model.th <- th.reform(th.list = model.th)
print(model.th)

# get MTSS threshold == 7th value
mtss.list <- list(model.th[[1]][7, ], model.th[[2]][7, ], model.th[[3]][7, ], model.th[[4]][7, ], model.th[[5]][7, ])

# now generate binary maps
bin <- bin_maker(preds = kor.pred, th = mtss.list)
print(bin)

# plot binary
plot_preds(preds = bin, poly = kr_poly, ncol = 5, nrow = 1, colors = rev(terrain.colors(1000)), pred.names = names(bin))
ggsave('output/plot/binary_output.png', width = 35, height = 8, dpi = 600, units = 'cm')


#####  part 12 ::: extrapolation risk ----------

# conduct MESS
mess <- ntbox::ntb_mess(M_stack = envs.sub, G_stack = envs.prj)

print(mess)
plot(mess)

# plot 
                                                          