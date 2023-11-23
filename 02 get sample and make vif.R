
# Load libraries ----------------------------------------------------------
require(pacman)
pacman::p_load(terra, tidyverse, randomForest, geodata, sf, fs, glue, usdm, rnaturalearthdata, rnaturalearth)

g <- gc(reset = T)
rm(list = ls())
options(scipen = 999, warn = -1)

# Load data ---------------------------------------------------------------

# Administrative data 
mex0 <- gadm(country = 'MEX', level = 0, path = 'tmpr')
mex1 <- gadm(country = 'MEX', level = 1, path = 'tmpr')
wrld <- ne_countries(returnclass = 'sf', scale = 50)

# Incendios
fire <- as.character(dir_ls('./data/tif/fires', regexp = '.tif$'))
nmes.fire <- gsub('.tif$', '', basename(fire))
fire <- rast(fire)
names(fire) <- nmes.fire
nmes.fire

# Variables 
vars <- as.character(dir_ls('./data/tif/vars_crnt', regexp = '.tif$'))
nmes.vars <- gsub('.tif', '', basename(vars))
vars <- terra::rast(vars)
names(vars) <- nmes.vars

names(vars)

# To select  --------------------------------------------------------------
crna <- fire[[1]]
crds <- terra::as.data.frame(crna, xy = T)
crds <- as_tibble(crds)
names(crds) <- c('x', 'y', 'fire_corona')

dir_create('./data/tbl/fires')
write.csv(crds, './data/tbl/fires/crds_fires-values.csv', row.names = FALSE)

# To select a sample for working with RF ----------------------------------
crds.smpl <- crds %>% sample_n(size = nrow(crds) * 0.33, replace = FALSE)
write.csv(crds.smpl, './data/tbl/fires/crds_fires-values_sample.csv', row.names = FALSE)

# To extract the values  --------------------------------------------------
vles <- terra::extract(vars, crds.smpl[,1:2]) %>% 
  cbind(crds.smpl, .)

# To make the VIF  --------------------------------------------------------
vif.step <- vifstep(x = vles[,5:ncol(vles)], th = 5)
vrs <- vif.step@results$Variables

nlyr(vars) 
length(vrs)

dir.create('./rds')
saveRDS(object = vrs, file = './rds/vars_vif.rds')
write.csv(vles, './data/tbl/fires/crds_fires-values_sample_vars.csv', row.names = FALSE)










