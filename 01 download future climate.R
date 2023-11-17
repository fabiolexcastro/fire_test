
# Load libraries ----------------------------------------------------------
require(pacman)
pacman::p_load(terra, geodata, tidyverse,  randomForestExplainer, caret, randomForest, geodata, sf, fs, glue, usdm, rnaturalearthdata, rnaturalearth)

g <- gc(reset = T)
rm(list = ls())
options(scipen = 999, warn = -1)

# Load data ---------------------------------------------------------------

# Variables
fles <- dir_ls('./data/tif/vars_crnt') %>% as.character
nmes <- fles %>% basename %>% gsub('.tif', '', .)
stck <- fles %>% rast()
names(stck) <- nmes
mask <- stck[[1]] * 0 + 1

# Administrative data 
mex0 <- gadm(country = 'MEX', level = 0, path = 'tmpr')
mex1 <- gadm(country = 'MEX', level = 1, path = 'tmpr')

# Parameters
gcms <- c("ACCESS-CM2", "ACCESS-ESM1-5", "AWI-CM-1-1-MR", "BCC-CSM2-MR", "CanESM5", "CanESM5-CanOE", "CMCC-ESM2", "CNRM-CM6-1", "CNRM-CM6-1-HR", "CNRM-ESM2-1", "EC-Earth3-Veg", "EC-Earth3-Veg-LR", "FIO-ESM-2-0", "GFDL-ESM4", "GISS-E2-1-G", "GISS-E2-1-H", "HadGEM3-GC31-LL", "INM-CM4-8", "INM-CM5-0", "IPSL-CM6A-LR", "MIROC-ES2L", "MIROC6", "MPI-ESM1-2-HR", "MPI-ESM1-2-LR", "MRI-ESM2-0", "UKESM1-0-LL")
ssps <- c('126', '245', '370', '585')
prds <- c('2021-2040', '2041-2060', '2061-2080', '2081-2100')

# Function to use ---------------------------------------------------------
down.clma <- function(ssp, gcm, prd){
  
  # gcm <- gcms[1]
  # ssp <- ssps[1]
  # prd <- prds[1]
  
  cat('>>> To process: ', ssp, ' ', gcm, ' ', prd, '<<<\n')
  
  bio <- cmip6_world(model = gcm, ssp = ssp, time = prd, var = 'bioc', res = 10, path = './tmpr')
  bio <- terra::crop(bio, mex0)
  bio <- terra::mask(bio, mex0)
  bio <- terra::resample(bio, mask, method = 'bilinear')
  
  dir <- glue('./data/tif/ftre/raw/{ssp}/{prd}')
  dir_create(dir)
  terra::writeRaster(x = bio, filename = glue('{dir}/bio_{gcm}.tif'), overwrite = TRUE)
  cat('Done!\n')
  
}

# To apply the function ---------------------------------------------------


map(.x = 1:length(ssps), .f = function(s){
  
  map(.x = 1:length(prds), .f = function(p){
    
    map(.x = 1:length(gcms), .f = function(g){
      
      cat('To process: ', ssps[s], ' ', prds[p], ' ', gcms[g], '\n')    
      
      try(expr = {
        
        down.clma(ssp = ssps[s], gcm = gcms[g], prd = prds[p])
        
      })
      
    })
    
  })
  
})

