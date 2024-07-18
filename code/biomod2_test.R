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
cv.k <- bm_CrossValidation(bm.format = myBiomodData.multi,
                           strategy = 'kfold',
                           nb.rep = 2,
                           k = 10)

head(cv.k)

## geographically stratified method
cv.s <- bm_CrossValidation(bm.format = myBiomodData.multi,
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
                            bm.format = myBiomodData.multi)

print(opt.t)


#####  part 4 ::: run single models ------------------------------ '
# run models
myBiomodModelOut <- BIOMOD_Modeling(bm.format = myBiomodData.multi,
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

# Represent evaluation scores & variables importance