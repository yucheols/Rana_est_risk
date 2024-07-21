######  follow the biomod2 tutorial in the biomod2 website (https://biomodhub.github.io/biomod2/index.html) 

# clear working environment
rm(list = ls(all.names = T))

# load packages
library(biomod2)
library(terra)

# attach data
data('DataSpecies')
head(DataSpecies)

# select the name of the study species
myRespName <- 'GuloGulo'

# get corresponding p/a info
myResp <- as.numeric(DataSpecies[, myRespName])
head(myResp)

# get corresponding coordinate data
myRespXY <- DataSpecies[, c('X_WGS84', 'Y_WGS84')]
head(myRespXY)

# load environmental variables extracted from BIOCLIM (bio_3, bio_4, bio_7, bio_11 & bio_12)
data('bioclim_current')
myExpl <- rast(bioclim_current)
print(myExpl)


#####  part 1 ::: data prep ------------------------------ 
## format data with true absences
myBiomodData <- BIOMOD_FormatingData(resp.var = myResp,
                                     expl.var = myExpl,
                                     resp.xy = myRespXY,
                                     resp.name = myRespName)

## pseudoabsence extraction == turn 0 into NA
myResp.PA <- ifelse(myResp == 1, 1, NA)

## format data with pseudoabsences: random method
myBiomodData.r <- BIOMOD_FormatingData(resp.var = myResp.PA,
                                       expl.var = myExpl,
                                       resp.xy = myRespXY,
                                       resp.name = myRespName,
                                       PA.nb.rep = 4,
                                       PA.nb.absences = 1000,
                                       PA.strategy = 'random')

print(myBiomodData.r)
plot(myBiomodData.r)

## select multiple sets of pseudo-absences
myBiomodData.multi <- BIOMOD_FormatingData(resp.var = myResp.PA,
                                           expl.var = myExpl,
                                           resp.xy = myRespXY,
                                           resp.name = myRespName,
                                           PA.nb.rep = 4,
                                           PA.nb.absences = c(1000, 500, 500, 200),
                                           PA.strategy = 'random')

summary(myBiomodData.multi)
plot(myBiomodData.multi)


#####  part 2 ::: data partitioning == cross-validation ------------------------------  
## k-fold method
cv.k <- bm_CrossValidation(bm.format = myBiomodData,
                           strategy = 'kfold',
                           nb.rep = 2,
                           k = 10)

head(cv.k)

## geographically stratified method
cv.s <- bm_CrossValidation(bm.format = myBiomodData,
                           strategy = 'strat',
                           k = 2,
                           balance = 'presences',
                           strat = 'x')

head(cv.s)


#####  part 3 ::: check modeling options ------------------------------  
## bigboss parameters
opt.b <- bm_ModelingOptions(data.type = 'binary',
                            models = c('SRE','XGBOOST'),
                            strategy = 'bigboss')

print(opt.b)

## tuned parameters with formated data
opt.t <- bm_ModelingOptions(data.type = 'binary',
                            models = c('SRE', 'XGBOOST'),
                            strategy = 'tuned',
                            bm.format = myBiomodData)

print(opt.t)


#####  part 4 ::: run single models ------------------------------ 
# run models
myBiomodModelOut <- BIOMOD_Modeling(bm.format = myBiomodData,
                                    modeling.id = 'AllModels',
                                    CV.strategy = 'random',
                                    CV.nb.rep = 2,
                                    CV.perc = 0.8,
                                    OPT.strategy = 'bigboss',
                                    var.import = 3,
                                    metric.eval = c('ROC', 'TSS', 'BOYCE'),
                                    do.progress = T)

# Get evaluation scores & variables importance
get_evaluations(myBiomodModelOut)
get_variables_importance(myBiomodModelOut)

# plot evaluation scores
bm_PlotEvalBoxplot(bm.out = myBiomodModelOut)
bm_PlotEvalBoxplot(bm.out = myBiomodModelOut, group.by = c('algo', 'algo'))
bm_PlotEvalBoxplot(bm.out = myBiomodModelOut, group.by = c('algo', 'run'))

# plot variables importance
bm_PlotVarImpBoxplot(bm.out = myBiomodModelOut, group.by = c('expl.var', 'algo', 'algo'))
bm_PlotVarImpBoxplot(bm.out = myBiomodModelOut, group.by = c('expl.var', 'algo', 'run'))
bm_PlotVarImpBoxplot(bm.out = myBiomodModelOut, group.by = c('algo', 'expl.var', 'run'))

# plot response curves
bm_PlotResponseCurves(bm.out = myBiomodModelOut, 
                      models.chosen = get_built_models(myBiomodModelOut)[c(1:3, 12:14)],
                      fixed.var = 'median')

bm_PlotResponseCurves(bm.out = myBiomodModelOut, 
                      models.chosen = get_built_models(myBiomodModelOut)[c(1:3, 12:14)],
                      fixed.var = 'min')

bm_PlotResponseCurves(bm.out = myBiomodModelOut, 
                      models.chosen = get_built_models(myBiomodModelOut)[3],
                      fixed.var = 'median',
                      do.bivariate = TRUE)


#####  part 5 ::: project single models ------------------------------ 
myBiomodProj <- BIOMOD_Projection(bm.mod = myBiomodModelOut,
                                  proj.name = 'Current',
                                  new.env = myExpl,
                                  models.chosen = 'all',
                                  metric.binary = 'all',
                                  metric.filter = 'all',
                                  build.clamping.mask = T)

plot(myBiomodProj)


#####  part 6 ::: run ensemble models ------------------------------ 
myBiomodEM <- BIOMOD_EnsembleModeling(bm.mod = myBiomodModelOut,
                                      models.chosen = 'all',
                                      em.by = 'all',
                                      em.algo = c('EMmean', 'EMcv', 'EMci', 'EMmedian', 'EMca', 'EMwmean'),
                                      metric.select = 'TSS',
                                      metric.select.thresh = c(0.7),
                                      metric.eval = c('ROC', 'TSS', 'BOYCE'),
                                      var.import = 3,
                                      EMci.alpha = 0.05,
                                      EMwmean.decay = 'proportional')

print(myBiomodEM)

# Get evaluation scores & variables importance
get_evaluations(myBiomodEM)
get_variables_importance(myBiomodEM)

# plot evaluation scores
bm_PlotEvalMean(bm.out = myBiomodEM, group.by = 'full.name')
bm_PlotEvalBoxplot(bm.out = myBiomodEM, group.by = c('full.name', 'full.name'))

# plot variables importance
bm_PlotVarImpBoxplot(bm.out = myBiomodEM, group.by = c('expl.var', 'full.name', 'full.name'))
bm_PlotVarImpBoxplot(bm.out = myBiomodEM, group.by = c('expl.var', 'algo', 'merged.by.run'))
bm_PlotVarImpBoxplot(bm.out = myBiomodEM, group.by = c('algo', 'expl.var', 'merged.by.run'))

# plot response curves
bm_PlotResponseCurves(bm.out = myBiomodEM, 
                      models.chosen = get_built_models(myBiomodEM)[c(1, 6, 7)],
                      fixed.var = 'median')

bm_PlotResponseCurves(bm.out = myBiomodEM, 
                      models.chosen = get_built_models(myBiomodEM)[c(1, 6, 7)],
                      fixed.var = 'min')

bm_PlotResponseCurves(bm.out = myBiomodEM, 
                      models.chosen = get_built_models(myBiomodEM)[7],
                      fixed.var = 'median',
                      do.bivariate = TRUE)


#####  part 7 ::: project ensemble models ------------------------------ 

# Project ensemble models (from single projections)
myBiomodEMProj <- BIOMOD_EnsembleForecasting(bm.em = myBiomodEM, 
                                             bm.proj = myBiomodProj,
                                             models.chosen = 'all',
                                             metric.binary = 'all',
                                             metric.filter = 'all')

plot(myBiomodEMProj)

# Project ensemble models (building single projections)
myBiomodEMProj2 <- BIOMOD_EnsembleForecasting(bm.em = myBiomodEM, 
                                              bm.proj = myBiomodProj,
                                              models.chosen = 'all',
                                              metric.binary = 'all',
                                              metric.filter = 'all')

plot(myBiomodEMProj2)
