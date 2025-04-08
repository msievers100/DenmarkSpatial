## Denmark Scenarios
## Sievers and Brown
## Greedy optimisation for seagrass restoration

# DONT RUN THIS SCRIPT ON ITS OWN
# It should be sourced from 3_source-greedy-new.R

library(ggplot2)
library(mgcv)
library(tidyr)
library(patchwork)
library(sf)
library(tmap)
library(dplyr)
library(viridis)
library(visreg)
library(patchwork)
library(gratia)


############

source("scripts/optimisation-Functions-new2-parallel.R")
grid_area_ratio <- readRDS("outputs/grid_area_ratio.rds")
grid_area_input <- readRDS("outputs/grid_area.rds") #!

#
# Parameters 
#
#date to be used in file name for saving output
# suggest you update in case you want to go back and 
# check older results

#number of iterations for the greedy algorithm
# higher values will take longer to run but will give more accurate results
niter <- 100
nstarts <- 100 #number of random starts for the greedy algorithm #20
ntol <- 1 #tolerance for improvement in outcomes (in units of total species n)
# . higher values will stop sooner b
# but will be further from the optimal solution
start_seeds <-  1:nstarts #random seeds for the greedy algorithm. So it can be 
# replicated
n_transplants <- 3307 #number of transplants to be made. We will optimise 
# arrangement of these transplants 
## 3973 total
## 3307 bare (10% of this is 331)

theme_set(theme_classic())


#
# Data 
#

buffer_input <- readRDS("outputs/buffers.rds")
grid_area <- readRDS("outputs/grid_area.rds") #!

#Read in spatial grid object
grid2_cropped <- readRDS("outputs/grid2_cropped.rds")
boulder <- readRDS("outputs/boulder.rds")
mussel_dissolved <- readRDS("outputs/mussel_dissolved.rds")

species_names[ispp_obj]
# calculate scaling values for each species in species_names
scalevals <- df %>% 
  select(all_of(species_names)) %>% 
  summarise_all(sd) %>% 
  as.numeric()
scalevalsdf <- data.frame(species = species_names, scalevals = scalevals)
write.csv(scalevalsdf, "outputs/scalevals.csv")

###################################
# STEP 2: Calculate marginal value of every grid cell 
###################################

# Calculate current summed abundance across target species
grid_temp <- compute_summed_abundances(grid2_cropped, gam_list, scalevals, ispp_obj, grid_area_ratio)
grid2_cropped$current_sum <- grid_temp$summed_abundance

# Create copy of data and change all 'Bare' treatments to '4-years-old' restored
grid2_treated <- grid2_cropped
grid2_treated$Treatment[grid2_treated$Treatment == "Bare"] <- "4-years-old"

# Calculate new summed abundance
grid2_treated <- compute_summed_abundances(grid2_treated, gam_list, scalevals, ispp_obj, grid_area_ratio)
grid2_cropped$treated_sum <- grid2_treated$summed_abundance

# Calculate marginal value
grid2_cropped$marginal_value <- with(grid2_cropped, treated_sum - current_sum)
grid2_cropped$vals <- 0
grid2_cropped$vals[order(grid2_cropped$marginal_value, decreasing = TRUE)[1:n_transplants]] <- 1
# Visualize marginal values
tmarg <- tm_shape(grid2_cropped) +
  tm_dots(col = "marginal_value", 
          style = "cont", 
          palette = "viridis", 
          size = 0.1,
          title = "Marginal Value\nof Restoration") +
  tm_layout(legend.outside = TRUE)
tmarg

#Plot to change the colour of Natural cells (which have a marginal value of 0, but perhaps are better shown differently...)
#tmarg2 <- tm_shape(grid2_cropped[grid2_cropped$Treatment == "Bare",]) +
#  tm_dots(col = "marginal_value", 
#          style = "cont", 
#          palette = "viridis", 
#          size = 0.1,
#          title = "Marginal Value\nof Restoration") +
#  tm_shape(grid2_cropped[grid2_cropped$Treatment == "Natural",]) +
#  tm_dots(col = "Treatment",
#          palette = c("darkolivegreen"),
#          size = 0.1) +
#  tm_layout(legend.outside = TRUE)

#Setting scale for palette so all objectives are on the same colour scale
tmarg2 <- tm_shape(grid2_cropped[grid2_cropped$Treatment == "Bare",]) +
  tm_dots(col = "marginal_value", 
          style = "cont", 
          palette = "viridis",
          breaks = seq(-2, 6, by = 0.5),
          size = 0.1,
          title = "Marginal Value\nof Restoration") +
  tm_shape(grid2_cropped[grid2_cropped$Treatment == "Natural",]) +
  tm_dots(col = "Treatment",
          palette = c("darkolivegreen"),
          size = 0.1) +
  tm_layout(legend.outside = TRUE)
tmarg2
# Save the plot
tmap_save(tmarg2, filename = paste0("outputs/plots/marginal-values_",obj_name,"_", run_date,".png")#,
        #  width = 3, height = 3, dpi = 300
        )

###################################
# STEP 3: Implement optimization
###################################

# Run the best guess optimization
best_guess_results <- best_guess_optimize(grid2_cropped, n_transplants, update_area10m = update_area10m)
best_guess_grid <- best_guess_results$final_grid

# Plot progress
progress_df <- do.call(rbind, lapply(best_guess_results$progress, function(x) {
  data.frame(n_transplants = x$n_transplants, abundance = x$abundance)
}))

###################################
# STEP 4: SAVE RESULTS
###################################

saveRDS(best_guess_results, file = paste0("outputs/greedy-results/best_guess_results_",obj_name,"_", run_date,".rds"))

