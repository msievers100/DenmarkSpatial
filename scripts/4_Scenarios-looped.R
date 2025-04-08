## Denmark Scenarios
## Sievers and Brown

rm(list =ls())
library(ggplot2)
library(mgcv)
library(tidyr)
library(patchwork)
library(sf)
library(tmap)
library(dplyr)
library(viridis)
library(visreg)
library(gratia)

source("scripts/Scenario-Functions.R")

#
# Parameters 
#
#date to be used in file name for saving output
# suggest you update in case you want to go back and 
# check older results
run_date <- "2024-08-19"

#number of steps for transplant scenarios, from 0 to max bare area
nsteps <- 20

#number of iterations for random transplantation scenario
#use more to get cleaner result, but will take longer
niter <- 10

#MICHAEL - update this to be the ratio of grid_area (ie area of our modelled grids)
# over the area of the survey grids (about 3.6 m^2?)
grid_area_ratio <- 20.32833/3.6


theme_set(theme_classic())

#save all the scenario maps (will slow it down)? 
#not recommended to do if you are running all species
# if set to TRUE, recommending just running loop for 1 species
doplots <- FALSE 

ndraws <- 500

#
# Data 
#

buffer_input <- readRDS("outputs/buffers.rds")
grid_area_input <- readRDS("outputs/grid_area.rds") #!

#Read in spatial grid object
grid2_cropped <- readRDS("outputs/grid2_cropped.rds")
boulder <- readRDS("outputs/boulder.rds")
mussel_dissolved <- readRDS("outputs/mussel_dissolved.rds")

#Read in GAMs for all species
gam_list <- readRDS("outputs/gam-selected-species.rds")
species_names <- names(gam_list)

######
# Model predictions for seagrass restoration 
# Note: Each scenario now calculates abundances for all species simultaneously.
# Results dataframes will have columns for each species abundance (species_n)
# plus scenario metadata (area, transplanted_grids, etc)
#####

# Calculate initial predictions for all species
for(species in species_names) {
  grid2_cropped[[species]] <- predict(gam_list[[species]], newdata = grid2_cropped, type = "response")
  
  #Map initial predictions
  tmap_obj <- tm_shape(grid2_cropped) +
    tm_dots(col = species, style = "cont", palette = "viridis", size = 0.1)
  tmap_obj
  tmap_save(tmap_obj, paste0("outputs/plots/",species,"-InitialPred.png"), dpi = 300)
}

#Scenarios of random transplantation
xrand <- lapply(1:niter, 
                function(j){
                  transplant_random(nsteps, gam_list, grid2_cropped,
                                    buffers = buffer_input, grid_area = grid_area_input,
                                    grid_area_ratio = grid_area_ratio,
                                    save_plots = doplots,
                                    scenario_name = j, seed = j*14)}
                )

rand_scnrs <- bind_rows(xrand) 
    
# Example plot for random scenarios (uncomment and modify for specific species)
 ggplot(rand_scnrs, aes(x = area, y = Lesser.pipefish_n, color = scenario_name, group = scenario_name)) +
   geom_line() +
   labs(title = "Random Transplantation Scenarios for Lesser Pipefish",
        x = "Area", y = "Predicted Abundance") +
   theme_minimal()

#From bottom right ("starting away from seagrass")
bbox <- st_bbox(grid2_cropped)
bottom_right_point <- st_point(c(bbox["xmax"], bbox["ymin"]))
bottom_right_sfc <- st_sfc(bottom_right_point, crs = st_crs(grid2_cropped))
seagrass_scnrs_awaySG <- transplant_scenario(nsteps, bottom_right_sfc, 
                                             gam_list, grid2_cropped,
                                             buffers = buffer_input, grid_area = grid_area_input,
                                             grid_area_ratio = grid_area_ratio,
                                             doplots, "start_away_from_seagrass")

# Example with Scenario with CIs
seagrass_scnrs_awaySG_cis <- transplant_scenario_cis(nsteps, bottom_right_sfc, 
                                                     gam_list, grid2_cropped,
                                                     buffers = buffer_input, grid_area = grid_area_input,
                                                     grid_area_ratio = grid_area_ratio,
                                                     doplots, "start_away_from_seagrass", 
                                                     ndraws = ndraws, seed = 4)
#From top ("starting near seagrass")
coords <- st_coordinates(grid2_cropped)
max_y <- max(coords[, "Y"])
top_points <- grid2_cropped %>% 
  filter(st_coordinates(.)[, "Y"] == max_y)
top_point <- top_points[1, ]
seagrass_scnrs_nearSG <- transplant_scenario(nsteps, top_point, gam_list, grid2_cropped,
                                             buffers = buffer_input, grid_area = grid_area_input,
                                             grid_area_ratio = grid_area_ratio,
                                             doplots, "start_near_seagrass")

# Example with Scenario with CIs
seagrass_scnrs_nearSG_cis <- transplant_scenario_cis(nsteps, top_point, 
                                                     gam_list, grid2_cropped,
                                                     buffers = buffer_input, grid_area = grid_area_input,
                                                     grid_area_ratio = grid_area_ratio,
                                                     doplots, "start_near_seagrass", 
                                                     ndraws = ndraws, seed = 4)   
    
#create a scenario where we start at three points close to boulder reef
#find points closest to boulder reef
nearest_boulder <- st_distance(grid2_cropped, boulder)
iclosest <- order(nearest_boulder)[1:3]
#find three nearest points and get their coordinates
grid_close_boulder <- grid2_cropped[iclosest,]
seagrass_scnrs_near_reef <- transplant_scenario(nsteps, grid_close_boulder, 
                                                gam_list, grid2_cropped,
                                                buffers = buffer_input, grid_area = grid_area_input,
                                                grid_area_ratio = grid_area_ratio,
                                                doplots, "start_at_boulder")

# Example with Scenario with CIs
seagrass_scnrs_boulder_cis <- transplant_scenario_cis(nsteps, grid_close_boulder, 
                                                     gam_list, grid2_cropped,
                                                     buffers = buffer_input, grid_area = grid_area_input,
                                                     grid_area_ratio = grid_area_ratio,
                                                     doplots, "start_at_boulder", 
                                                     ndraws = ndraws, seed = 4)   
    
#scenario starting near mussel reef
nearest_mussel <- st_distance(grid2_cropped, mussel_dissolved)
iclosest <- order(nearest_mussel)[1:30] ## I had to make this larger than 3, so some sites were near both of the 2 reefs
grid_close_mussel <- grid2_cropped[iclosest,]
seagrass_scnrs_near_mussel <- transplant_scenario(nsteps, grid_close_mussel, 
                                                  gam_list, grid2_cropped,
                                                  buffers = buffer_input, grid_area = grid_area_input,
                                                  grid_area_ratio = grid_area_ratio,
                                                  doplots, "start_at_mussel")

# Example with Scenario with CIs
seagrass_scnrs_mussel_cis <- transplant_scenario_cis(nsteps, grid_close_mussel, 
                                                     gam_list, grid2_cropped,
                                                     buffers = buffer_input, grid_area = grid_area_input,
                                                     grid_area_ratio = grid_area_ratio,
                                                     doplots, "start_at_mussel", 
                                                     ndraws = ndraws, seed = 4)   

# Bind all scenarios together
all_scnrs <- bind_rows(seagrass_scnrs_awaySG, seagrass_scnrs_nearSG, seagrass_scnrs_near_reef, seagrass_scnrs_near_mussel)

# Save output of runs
saveRDS(all_scnrs, file = paste0("outputs/",run_date,"scenarios-distance.rds"))
saveRDS(rand_scnrs, file = paste0("outputs/",run_date,"scenarios-random.rds"))

# Save confidence interval scenarios
all_scnrs_cis <- bind_rows(seagrass_scnrs_awaySG_cis, seagrass_scnrs_nearSG_cis, 
                          seagrass_scnrs_boulder_cis, seagrass_scnrs_mussel_cis)
saveRDS(all_scnrs_cis, file = paste0("outputs/",run_date,"scenarios-distance-cis.rds"))
