
# Load libraries ----------------------------------------------------------
require(pacman)
pacman::p_load(terra, tidyverse,  randomForestExplainer, caret, randomForest, geodata, sf, fs, glue, usdm, rnaturalearthdata, rnaturalearth)
pacman::p_load(raster, rgdal, pROC, cptcity, randomForest, mlbench, caret, rgeos, stringr, sf, dismo, tidyverse, geodata, glue, fs)

g <- gc(reset = T)
rm(list = ls())
options(scipen = 999, warn = -1)

# Functions to use --------------------------------------------------------
raster01 <- function(r){
  minmax_r <- range(values(r), na.rm=TRUE) # get the min max values
  return( (r-minmax_r[1]) / (diff(minmax_r))) # rescale 
}

# Load data ---------------------------------------------------------------
crna <- terra::rast('./data/tif/fires/Corona_Gg.tif')
vles <- read_csv('./data/tbl/fires/crds_fires-values_sample_vars.csv')
vars <- readRDS(file = './rds/vars_vif.rds')
fles <- dir_ls('./data/tif/vars_crnt') %>% as.character
nmes <- fles %>% basename %>% gsub('.tif', '', .)
stck <- fles %>% rast()
names(stck) <- nmes

# Select the variables from the stack -------------------------------------
stck <- stck[[grep(paste0(vars, collapse = '|'), names(stck))]]
nlyr(stck)

# Tunning random forest ---------------------------------------------------

set.seed(123)
vles <- drop_na(vles)
vles <- dplyr::select(vles, x, y, fire_corona, ID, vars)
mtrx <- vles[,c(3, 5:ncol(vles))]
mtrx <- drop_na(mtrx)

# MTRY
control <- trainControl(method = 'repeatedcv', number = 5, repeats = 3, search = 'random')
mtry <- round(sqrt(length(vars)), 0)
rf_random <- train(fire_corona ~., data = mtrx, method = "rf", metric = 'RMSE', tuneLength = 5, trControl = control)# MTRY = 3
saveRDS(object = rf_random, file = './rds/rf_random_run1.rds')

# Number of trees
tunegrid <- expand.grid(.mtry=c(sqrt(length(vars))))
ntrees <- c(100, 250, 500, 750, 1000, 1500, 2000, 2500)
# fit <- train(fire_corona ~ ., data = as.data.frame(mtrx), method = "rf", metric = 'RMSE', tuneGrid = tunegrid, trControl = control, ntree = 100)
mdls <- map(.x = 1:length(ntrees), .f = function(i){
  
  cat("Ntree: ", ntrees[i], '\n')
  set.seed(123)
  fit <- train(fire_corona ~ ., data = as.data.frame(mtrx), method = "rf", metric = 'RMSE', tuneGrid = tunegrid, trControl = control, ntree = ntrees[i])
  return(fit)
  
})
mdls
saveRDS(object = mdls, file = './rds/models.rds')
rslt <- resamples(mdls)
rslt <- summary(rslt)
saveRDS(object = rslt, file = './rds/models_summary.rds')

# MAE: mean absolute error (Model 8 - Mean, value lowest: 6.07)
# RMSE: root mean squared error (Model 8 - Mean, value lowest 11.54)
# RSquared: r squared (Model 4 - Mean, value lowest 11.1)

# Final model 8 (2500 trees)

# -------------------------------------------------------------------------
# To fit the model --------------------------------------------------------
# -------------------------------------------------------------------------
vles <- dplyr::select(vles, -ID)
fold <- kfold(vles, k = 25)
rm(rslt)

frml <- fire_corona ~ aei_rf_x10_025 + DEM_025 + DistRivers_km + DistRoads_km + DistRural_km + DistUrban_km + DistWaterBodies_km + FuelLoad_Mgha + GDP_2011 + PopDen_paxkm2_2011 + Population_2011 + slope_025 + Srad_mean_025 + wc2.1_30s_bio_15 + wc2.1_30s_bio_18 + wc2.1_30s_bio_19 + wc2.1_30s_bio_2 + wc2.1_30s_bio_3 + wc2.1_30s_bio_8 + wind_mean_025_x10 

rslt <- map(.x = 1:25, .f = function(i){
  
  cat('>>> To process: ', i, '<<<\n')
  tst <- vles[fold == i,]
  trn <- vles[fold != i,]
  
  tst <- drop_na(tst)
  trn <- drop_na(trn)
  
  # Presences and pseudo-absences
  env <- trn$fire_corona
  
  # Way 1
  rf1 <- randomForest::randomForest(frml, data = trn[,c(3:ncol(trn))], ntree = 2500, mtry = 3, na.action = na.omit, nodesize = 2,
                                    confusionMatrix = TRUE,
                                    do.trace = TRUE, keep.forest = TRUE, importance = TRUE, verbose = TRUE)
  
  rst <- terra::predict(stck, rf1, prob = TRUE, progress = 'text')
  
  # Importance variables 
  imp <- importance(rf1)
  imp <- as.data.frame(imp)
  imp <- rownames_to_column(imp)
  
  # AUC 
  prd <- as.numeric(predict(rf1, tst[,3:ncol(tst)]))
  obs <- as.vector(tst[,'fire_corona'])
  obs <- as.numeric(unlist(obs))
  auc <- auc(obs, prd) 
  
  # List the results
  lst <- list(rf1, rst, imp, auc)
  cat('Done!\n')
  return(lst)
  
})

rflist <- map(rslt, 1)
rffr <- do.call(randomForest::combine, rflist)
saveRDS(object = rffr, file = './rds/randomforest_model.rds')

rffr <- readRDS(file = './rds/randomforest_model.rds')

fnal <- predict(stck, rffr)
dir_create('./rf/run_1')
terra::writeRaster(x = fnal, filename = './rf/run_1/rf_combine-avg.tif', overwrite = T)

# To normalice
rstr <- raster01(r = fnal)

# AUCs
aucs <- map(rslt, 4)
aucs <- unlist(aucs)
aucs <- tibble(model = 1:25, auc = aucs)

write.csv(aucs, './rf/run_1/aucs.csv', row.names = FALSE)

impr <- map(rslt, 3)
impr <- map(.x = 1:25, .f = function(i){
  impr[[i]] %>% mutate(run = glue('run_{i}'))
})
impr <- bind_rows(impr)
impr <- as_tibble(impr)
impr <- impr %>% group_by(rowname) %>% summarise(`%IncMSE` = mean(`%IncMSE`), IncNodePurity = mean(IncNodePurity)) %>% ungroup()

gimpr <- ggplot(impr, aes(x=rowname, y=`%IncMSE`)) +
  geom_segment( aes(x=rowname, xend=rowname, y=0, yend=`%IncMSE`), color="skyblue") +
  geom_point(aes(size = IncNodePurity), color="blue", alpha=0.6) +
  labs(x = '%IncMSE', y = 'Variable') +
  theme_light() +
  coord_flip() +
  theme(
    legend.position="bottom",
    panel.grid.major.y = element_blank(),
    # panel.border = element_blank(),
    axis.ticks.y = element_blank()
  )

gimpr
dir_create('./png/graphs')
ggsave(plot = gimpr, filename = './png/graphs/run1_rfImp.png', units = 'in', width = 9, height = 7, dpi = 300)


# Check future climate ----------------------------------------------------
ftre.dirs <- dir_ls('./data/tif/ftre/avg') %>% as.character()

ssp126.prds <- ftre.dirs[1] %>% dir_ls() %>% map(., rast) %>% reduce(., c)
ssp245.prds <- ftre.dirs[2] %>% dir_ls() %>% map(., rast) %>% reduce(., c)
ssp370.prds <- ftre.dirs[3] %>% dir_ls() %>% map(., rast) %>% reduce(., c)
ssp575.prds <- ftre.dirs[4] %>% dir_ls() %>% map(., rast) %>% reduce(., c)

# Check bio 8 (same SSP)
ssp126.prds.b08 <- ssp126.prds[[grep('bio08', names(ssp126.prds), value = F)]]
ssp245.prds.b08 <- ssp245.prds[[grep('bio08', names(ssp245.prds), value = F)]]
ssp370.prds.b08 <- ssp370.prds[[grep('bio08', names(ssp370.prds), value = F)]]
ssp575.prds.b08 <- ssp575.prds[[grep('bio08', names(ssp575.prds), value = F)]]
plot(ssp575.prds.b08)

# Predict for the future --------------------------------------------------
rffr <- readRDS('./rds/randomforest_model.rds')
ftre.dirs <- dir_ls('./data/tif/ftre/avg') %>% as.character()
pos.bios <- names(stck) %>% grep('bio_', ., value = T) %>% str_sub(., nchar(.) - 1, nchar(.)) %>% parse_number()

# Stack with variables that are not climate1
stck.1 <- stck[[grep('_bio_', names(stck), value = F)]]
stck.2 <- stck[[-grep('_bio_', names(stck), value = F)]]

# Future parameters
ssps <- c('126', '245', '370', '585')
prds <- c('2021-2040', '2041-2060', '2061-2080', '2081-2100')

predict.ftre <- function(ssp, prd){
  
  cat('To process: ', ssp, prd, '\n')
  fle <- as.character(grep(prd, dir_ls(grep(ssp, ftre.dirs, value = T)), value = T))
  rst <- rast(fle)  
  sub <- rst[[pos.bios]]
  names(sub) <- names(stck.1)
  stk <- c(stck.2, sub)
  rsl <- terra::predict(stk, rffr)
  names(rsl) <- glue('rfrs_{ssp}_{prd}')
  return(rsl)
  
}

# To make the predict
prdc.ftre <- map(.x = 1:4, .f = function(s){
  rsl <- map(.x = 1:4, .f = function(p){
    rst <- predict.ftre(ssp = ssps[s], prd = prds[p])
    return(rst)
  })
  return(rsl)
}) %>% 
  reduce(., c)

names(prdc.ftre)

# Join future with the baseline 
fnal <- terra::rast('./rf/run_1/rf_combine-avg.tif')
names(fnal) <- glue('rfrs_current')
alls <- c(fnal, prdc.ftre)

dir_create('./rf/run_1')
terra::writeRaster(x = alls, filename = './rf/run_1/predict_alls.tif', overwrite = TRUE)

# To scale the rasters ----------------------------------------------------
raster01 <- function(r){
  minmax_r <- range(values(r), na.rm=TRUE) # get the min max values
  return( (r-minmax_r[1]) / (diff(minmax_r))) # rescale 
}
alls <- raster01(alls)

# To make the map  --------------------------------------------------------

# Administrative data 
mex0 <- gadm(country = 'MEX', level = 0, path = 'tmpr')
mex1 <- gadm(country = 'MEX', level = 1, path = 'tmpr')
wrld <- ne_countries(returnclass = 'sf', scale = 50)

# Tibble
tble <- terra::as.data.frame(alls, xy = T) %>% 
  as_tibble %>% 
  gather(var, value, -x, -y) 

lbls <- tble %>% 
  distinct(var) %>% 
  mutate(var2 = gsub('rfrs_', 'SSP', var)) %>% 
  mutate(var2 = gsub('SSPcurrent', 'Current', var2), 
         var2 = gsub('_', ' ', var2))

tble <- inner_join(tble, lbls, by = 'var') %>% 
  mutate(var2 = factor(var2, levels = lbls$var2))

gmap <- ggplot() + 
  geom_tile(data = tble, aes(x = x, y = y, fill = value)) + 
  facet_wrap(~var2) +
  scale_fill_gradientn(colors = cpt(pal = 'imagej_gyr_centre', n = 10, rev = FALSE)) +
  geom_sf(data = st_as_sf(mex0), fill = NA, col = 'grey40') +
  # geom_sf(data = st_as_sf(mex1), fill = NA, col = 'grey40') + 
  geom_sf(data = wrld, fill = NA, col = 'grey40') + 
  coord_sf(xlim = ext(mex0)[1:2], ylim = ext(mex0)[3:4]) + 
  labs(x = 'Lon', y = 'Lat', fill = 'Fire score') +
  theme_void() +
  theme(strip.text = element_text(face = 'bold'), 
        axis.text.y = element_text(size = 6, angle = 90, hjust = 0.5),
        axis.text.x = element_text(size = 6), 
        legend.position = 'bottom',
        legend.key.width = unit(3, 'line')) +
  guides(fill = guide_legend( 
    direction = 'horizontal',
    keyheight = unit(1.15, units = "mm"),
    keywidth = unit(15, units = "mm"),
    title.position = 'top',
    title.hjust = 0.5,
    label.hjust = .5,
    nrow = 1,
    byrow = T,
    reverse = F,
    label.position = "bottom"
  ))

ggsave(plot = gmap, filename = './png/maps/fire_allperiods.jpg', units = 'in', width = 15, height = 9.5, dpi = 300)

# To make a histogram for checking the results ----------------------------

ghist <- ggplot(data = tble, aes(x = value)) + 
  geom_histogram() +
  facet_wrap(~var2) +
  scale_x_continuous(trans = 'log10') +
  labs(x = 'Prob value', y = 'Pixel frequency') +
  theme_light() +
  theme(strip.text = element_text(face = 'bold'))

ggsave(plot = ghist, filename = './png/graphs/histogram_prob.jpg', units = 'in', width = 9, height = 8, dpi = 300)

# Finish ------------------------------------------------------------------
