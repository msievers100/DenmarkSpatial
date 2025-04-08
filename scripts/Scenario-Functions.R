transplant_scenario <- function(nareas, grid_center, gam_list, grid2_cropped, buffers, grid_area,
                                grid_area_ratio = 1,
                                save_plots = TRUE, 
                                scenario_name = "scenario") {

#Scenario of transplants starting from selected points and expanding by nearest distance to those
# starting points. 
# Inputs:
# nareas: number of steps for area of transplants
# grid_center: center of the grid for where transplants start. Can be multiple points
# gam_list: named list of GAM models, one for each species
# grid2_cropped: grid data with the grid points
# save_plots: logical, whether to save plots
# scenario_name: name of the scenario used for saving plots
  nbare <- sum(grid2_cropped$Treatment == "Bare")
  species_names <- names(gam_list)
  
  # Initialize results dataframe with columns for each species
  seagrass_scnrs <- data.frame(scenario_name = scenario_name,
                              area = 1:nareas,
                              transplanted_grids = round(seq(0, nbare, length.out = nareas)))
  seagrass_scnrs$transplanted_area <- seagrass_scnrs$transplanted_grids * grid_area
  
  # Add columns for each species abundance
  for(sp in species_names) {
    seagrass_scnrs[[paste0(sp,"_n")]] <- NA
    # Calculate initial abundance for each species
    seagrass_scnrs[[paste0(sp,"_n")]][1] <- sum(predict(gam_list[[sp]], 
                                                        newdata = grid2_cropped, 
                                                        type = "response"))*grid_area_ratio
  }
  
  igrids_no_seagrass <- which(grid2_cropped$Treatment == "Bare")
  
  # Find the grid point that is closest to the center
  middle <- which.min(st_distance(grid2_cropped, grid_center))
  dist_to_middle <- st_distance(grid2_cropped, grid_center) %>% 
    apply(1, min) 
  #find nearest distance
  
  
  dist_to_middle <- data.frame(dist = dist_to_middle, 
                               id = 1:length(dist_to_middle), 
                               InSeagrass = grid2_cropped$InSeagrass)
  dist_to_middle_InSeagrass <- filter(dist_to_middle, !InSeagrass)
  
  for (i in 2:nareas) {
    ngrids_transplanted <- seagrass_scnrs$transplanted_grids[i]
    # Select grid points starting with the closest
    transplanted <- order(dist_to_middle_InSeagrass$dist)[1:ngrids_transplanted]
    grid2_temp <- grid2_cropped
    grid2_temp$Treatment[dist_to_middle_InSeagrass$id[transplanted]] <- "4-years-old"
    
    # Update Area10m
    grid2_temp$InSeagrass[grid2_temp$Treatment == "4-years-old"] <- TRUE
    nseagrass <- st_intersects(buffers, grid2_temp[grid2_temp$InSeagrass,]) %>%
      lapply(length) %>% unlist()
    grid2_temp$Area10m <- nseagrass * grid_area
    
    # Calculate abundances for all species
    for(sp in species_names) {
      grid2_temp[[sp]] <- predict(gam_list[[sp]], newdata = grid2_temp, type = "response")*grid_area_ratio
      seagrass_scnrs[[paste0(sp,"_n")]][i] <- sum(grid2_temp[[sp]])
    }
    
    if (save_plots) {
      tm1 <- tm_shape(grid2_temp) +
        tm_dots(col = "Treatment", size = 0.5)
      tmap_save(tm1, filename = paste0("outputs/plots/scenarioplots/",scenario_name,"_", i, ".png"), dpi = 300)
    }
  }
  
  return(seagrass_scnrs)
}


transplant_random <- function(nareas, gam_list, grid2_cropped, 
                              buffers, grid_area,  grid_area_ratio = 1,
                              save_plots = TRUE,
                              scenario_name = "scenario", seed=42) {
#Scenario of transplants that selects random locations for transplants. 
# Runs multiple areas transplanted from zero to the total area of the grid.
# Inputs:
# nareas: number of steps for area of transplants
# gam_list: named list of GAM models, one for each species
# grid2_cropped: data frame of grid data with the grid points
# save_plots: logical, whether to save plots
# scenario_name: name of the scenario used for saving plots. Suggest you use
# numbers here for each random iteration
  nbare <- sum(grid2_cropped$Treatment == "Bare")
  species_names <- names(gam_list)
  
  # Initialize results dataframe with columns for each species
  seagrass_scnrs <- data.frame(scenario_name = scenario_name,
                              area = 1:nareas,
                              transplanted_grids = round(seq(0, nbare, length.out = nareas)))
  seagrass_scnrs$transplanted_area <- seagrass_scnrs$transplanted_grids * grid_area
  
  # Add columns for each species abundance
  for(sp in species_names) {
    seagrass_scnrs[[paste0(sp,"_n")]] <- NA
    # Calculate initial abundance for each species
    seagrass_scnrs[[paste0(sp,"_n")]][1] <- sum(predict(gam_list[[sp]], 
                                                        newdata = grid2_cropped, 
                                                        type = "response"))*grid_area_ratio
  }
    
  grid2_temp <- grid2_cropped
  set.seed(seed)
  for (i in 2:nareas) {
    igrids_no_seagrass <- which(grid2_temp$Treatment == "Bare")
    # Randomly select i grid points to be transplanted
    ngrids_transplanted <- seagrass_scnrs$transplanted_grids[i]
    transplanted <- sample(igrids_no_seagrass, 
                           ngrids_transplanted-seagrass_scnrs$transplanted_grids[i-1], 
                           replace = FALSE)
    
    grid2_temp$Treatment[transplanted] <- "4-years-old"

    # Update Area10m
    
    #New faster method
   
    #plot to check 
    # tm_shape(working_grid) + 
    #   tm_dots() + 
    #   tm_shape(working_grid[affected_cells,]) + 
    #   tm_dots(col = "red") + 
    #   tm_shape(working_grid[best_cell,]) + 
    #   tm_dots(col = "green")
    
    grid2_temp$InSeagrass[grid2_temp$Treatment == "4-years-old"] <- TRUE
    nseagrass <- st_intersects(buffers, grid2_temp[grid2_temp$InSeagrass,]) %>%
      lapply(length) %>% unlist()
    grid2_temp$Area10m <- nseagrass * grid_area
    
    # Calculate abundances for all species
    for(sp in species_names) {
      grid2_temp[[sp]] <- predict(gam_list[[sp]], newdata = grid2_temp, type = "response")*grid_area_ratio
      seagrass_scnrs[[paste0(sp,"_n")]][i] <- sum(grid2_temp[[sp]])
    }
    
    if (save_plots) {
      tm1 <- tm_shape(grid2_temp) +
        tm_dots(col = "Treatment", size = 0.5)
      tmap_save(tm1, filename = paste0("outputs/plots/scenarioplots/random",scenario_name,"_", i, ".png"), dpi = 300)
    }
  }
  
  return(seagrass_scnrs)
}


transplant_scenario_cis <- function(nareas, grid_center, gam_list, grid2_cropped, 
                                    buffers,  grid_area,  grid_area_ratio = 1,
                                    save_plots = TRUE, 
                                scenario_name = "scenario", ndraws, seed = seed) {
  #Scenario of transplants starting from selected points and expanding by nearest distance to those
  # starting points. 
  # Generates predictions for multiple draws from the posterior
  # Inputs:
  # nareas: number of steps for area of transplants
  # grid_center: center of the grid for where transplants start. Can be multiple points
  # gam_list: named list of GAM models, one for each species
  # grid2_cropped: grid data with the grid points
  # save_plots: logical, whether to save plots
  # scenario_name: name of the scenario used for saving plots
  # ndraws: number of draws from posterior
  # seed: random seed for posterior draws
  nbare <- sum(grid2_cropped$Treatment == "Bare")
  species_names <- names(gam_list)
  
  # Initialize base dataframe
  base_df <- data.frame(scenario_name = scenario_name,
                       area = 1:nareas,
                       transplanted_grids = round(seq(0, nbare, length.out = nareas)))
  
  # Create draws dataframe
  draws_df <- data.frame(.draw = 1:ndraws)
  
  # Cross join to get all combinations
  seagrass_scnrs <- base_df %>%
    cross_join(draws_df)
  
  seagrass_scnrs$transplanted_area <- seagrass_scnrs$transplanted_grids * grid_area
  
  # Add columns for each species abundance
  for(sp in species_names) {
    seagrass_scnrs[[paste0(sp,"_n")]] <- NA
    # Calculate initial abundance for each species with uncertainty
    dfits <- fitted_samples(gam_list[[sp]], data = grid2_cropped, 
                           n = ndraws, seed = seed, scale = "response") %>%
      group_by(.draw) %>%
      summarize(n = sum(.fitted)*grid_area_ratio) %>%
      arrange(.draw)
    
    #match on .draw and area
    seagrass_scnrs[[paste0(sp,"_n")]][seagrass_scnrs$area == 1] <- dfits$n
  }
  
  igrids_no_seagrass <- which(grid2_cropped$Treatment == "Bare")
  
  # Find the grid point that is closest to the center
  middle <- which.min(st_distance(grid2_cropped, grid_center))
  dist_to_middle <- st_distance(grid2_cropped, grid_center) %>% 
    apply(1, min) 
  #find nearest distance
  
  dist_to_middle <- data.frame(dist = dist_to_middle, 
                               id = 1:length(dist_to_middle), 
                               InSeagrass = grid2_cropped$InSeagrass)
  dist_to_middle_InSeagrass <- filter(dist_to_middle, !InSeagrass)
  
  for (i in 2:nareas) {
    ngrids_transplanted <- seagrass_scnrs$transplanted_grids[seagrass_scnrs$area ==i][1]
    # Select grid points starting with the closest
    transplanted <- order(dist_to_middle_InSeagrass$dist)[1:ngrids_transplanted]
    grid2_temp <- grid2_cropped
    grid2_temp$Treatment[dist_to_middle_InSeagrass$id[transplanted]] <- "4-years-old"
    
    # Update Area10m
    grid2_temp$InSeagrass[grid2_temp$Treatment == "4-years-old"] <- TRUE
    nseagrass <- st_intersects(buffers, grid2_temp[grid2_temp$InSeagrass,]) %>%
      lapply(length) %>% unlist()
    grid2_temp$Area10m <- nseagrass * grid_area
    
    # Calculate abundances for all species with uncertainty
    for(sp in species_names) {
      dfits <- fitted_samples(gam_list[[sp]], data = grid2_temp, 
                             n = ndraws, seed = seed, scale = "response") %>%
        group_by(.draw) %>%
        summarize(n = sum(.fitted)*grid_area_ratio) %>%
        arrange(.draw)
      
      #match on .draw and area
      seagrass_scnrs[[paste0(sp,"_n")]][seagrass_scnrs$area == i] <- dfits$n
    }
    
    
    if (save_plots) {
      tm1 <- tm_shape(grid2_temp) +
        tm_dots(col = "Treatment", size = 0.5)
      tmap_save(tm1, filename = paste0("outputs/plots/scenarioplots/",scenario_name,"_", i, ".png"), dpi = 300)
    }
  }
  
  return(seagrass_scnrs)
}


#Greedy algorithm to optimise transplants 

transplant_greedy <- function(n_transplants, gam_list, grid2_cropped, 
                              buffers, grid_area, grid_area_ratio = 1,
                              save_plots = FALSE,
                              scenario_name = "scenario", seed = 42, maxiter = 10,
                              tol = 10, scalevals = NULL) {
  # Tries to optimize locations of a fixed number of transplants to maximize
  # combined abundance of multiple species. 
  # Uses a simple greedy algorithm (only make better moves)
  # Inputs:
  # n_transplants: number of transplants to make. 
  # gam_list: named list of GAM models, one for each species
  # grid2_cropped: data frame of grid data with the grid points
  # save_plots: logical, whether to save plots
  # scenario_name: name of the scenario used for saving plots
  # seed: random seed to initialize first transplant arrangement
  # maxiter: maximum number of iterations to try to improve the arrangement
  # tol: tolerance for improvement in abundance required to keep iterating
  # scalevals: named vector of scaling factors for each species (optional)
  
  nbare <- sum(grid2_cropped$Treatment == "Bare")
  if (n_transplants > nbare) {
    stop("Number of transplants is greater than the number of bare grid points")
  }
  
  species_names <- names(gam_list)
  
  # If no scaling values provided, use 1 for all species
  if (is.null(scalevals)) {
    scalevals <- setNames(rep(1, length(species_names)), species_names)
  }
  
  # Initialize results dataframe with columns for each species
  seagrass_scnrs <- data.frame(scenario_name = scenario_name,
                              iter = 1:maxiter)
  
  # Add columns for each species abundance
  for(sp in species_names) {
    seagrass_scnrs[[paste0(sp,"_n")]] <- NA
  }
  
  # Add column for total scaled abundance
  seagrass_scnrs$total_scaled_n <- NA
  
  #Initial arrangement and initial abundances
  set.seed(seed)
  grid2_temp <- grid2_cropped
  igrids_no_seagrass <- which(grid2_temp$Treatment == "Bare")
  # Randomly select points to be transplanted
  transplanted <- sample(igrids_no_seagrass, n_transplants, replace = FALSE)
  
  grid2_temp$Treatment[transplanted] <- "4-years-old"
  
  # Update Area10m
  grid2_temp$InSeagrass[grid2_temp$Treatment == "4-years-old"] <- TRUE
  nseagrass <- st_intersects(buffers, grid2_temp[grid2_temp$InSeagrass,]) %>%
    lapply(length) %>% unlist()
  grid2_temp$Area10m <- nseagrass * grid_area
  
  # Calculate initial abundances for all species
  total_scaled_n <- 0
  for(sp in species_names) {
    grid2_temp[[sp]] <- predict(gam_list[[sp]], newdata = grid2_temp, type = "response")*grid_area_ratio
    seagrass_scnrs[[paste0(sp,"_n")]][1] <- sum(grid2_temp[[sp]])
    total_scaled_n <- total_scaled_n + seagrass_scnrs[[paste0(sp,"_n")]][1] / scalevals[sp]
  }
  seagrass_scnrs$total_scaled_n[1] <- total_scaled_n
  
  tol_count <- 0
  for (i in 2:maxiter) {
    igrids_no_seagrass <- which(grid2_temp$Treatment == "Bare")
    igrids_transplanted <- which(grid2_temp$Treatment == "4-years-old")
    # Randomly select 1 bare grid point to be transplanted and one
    # transplanted grid to be bare. This ensures area remains the same. 
    
    transplanted <- sample(igrids_no_seagrass, 1, 
                           replace = FALSE)
    untransplanted <- sample(igrids_transplanted, 1, 
                             replace = FALSE)
    
    grid2_temp$Treatment[transplanted] <- "4-years-old"
    grid2_temp$Treatment[untransplanted] <- "Bare"
    
    # Update Area10m
    grid2_temp$InSeagrass[grid2_temp$Treatment == "4-years-old"] <- TRUE
    
    nseagrass <- st_intersects(buffers, grid2_temp[grid2_temp$InSeagrass,]) %>%
      lapply(length) %>% unlist()
    grid2_temp$Area10m <- nseagrass * grid_area
    
    # Calculate abundances for all species
    total_scaled_n <- 0
    for(sp in species_names) {
      grid2_temp[[sp]] <- predict(gam_list[[sp]], newdata = grid2_temp, type = "response")*grid_area_ratio
      ntemp <- sum(grid2_temp[[sp]])
      total_scaled_n <- total_scaled_n + ntemp / scalevals[sp]
    }
    
    #See if total_scaled_n is greater than previous iteration
    if (total_scaled_n > seagrass_scnrs$total_scaled_n[i-1]) {
      for(sp in species_names) {
        seagrass_scnrs[[paste0(sp,"_n")]][i] <- sum(grid2_temp[[sp]])
      }
      seagrass_scnrs$total_scaled_n[i] <- total_scaled_n
    } else {
      grid2_temp$Treatment[transplanted] <- "Bare"
      grid2_temp$Treatment[untransplanted] <- "4-years-old"
      for(sp in species_names) {
        seagrass_scnrs[[paste0(sp,"_n")]][i] <- seagrass_scnrs[[paste0(sp,"_n")]][i-1]
      }
      seagrass_scnrs$total_scaled_n[i] <- seagrass_scnrs$total_scaled_n[i-1]
    }
    
    if (save_plots) {
      tm1 <- tm_shape(grid2_temp) +
        tm_dots(col = "Treatment", size = 0.5)
      tmap_save(tm1, filename = paste0("outputs/plots/scenarioplots/greedy/optimized",scenario_name,"_", i, ".png"), dpi = 50)
    }
    
    #This stops iterating if the random changes don't improve the abundance more than three times in a row
    if (ntemp - seagrass_scnrs$n[i-1] < tol) {
      tol_count <- tol_count + 1
    } else {
      tol_count <- 0
    }
    if (tol_count>(niter/5)){
      break
    }
  }
  
  return(list(seagrass_scnrs = seagrass_scnrs, grid2_temp = grid2_temp))
}



optimize_transplants <- function(n_transplants, grid2_cropped, gam_list, 
                                 scalevals, ispp_obj, n_starts = 100) {
  
  # Initialize best solution trackers
  best_abundance <- -Inf
  best_grid <- grid2_cropped
  run_optimization(1)
  # Run multiple starts
  results <- lapply(1:n_starts, function(i) run_optimization(i))
  
  # Find best result
  abundances <- sapply(results, function(x) x$abundance)
  best_idx <- which.max(abundances)
  
  list(
    best_abundance = results[[best_idx]]$abundance,
    best_grid = results[[best_idx]]$grid,
    all_abundances = abundances
  )
}

# Function to run one optimization from random start
run_optimization <- function(seed) {
  set.seed(seed)
  
  # Create working copy of grid
  grid_temp <- grid2_cropped
  
  # Initialize with random transplants
  bare_cells <- which(grid_temp$Treatment == "Bare")
  initial_cells <- sample(bare_cells, n_transplants)
  grid_temp$Treatment[initial_cells] <- "4-years-old"
  
  improved <- TRUE
  while(improved) {
    improved <- FALSE
    
    # Calculate current abundance
    grid_temp <- compute_summed_abundances(grid_temp, gam_list, scalevals, ispp_obj)
    current_abundance <- sum(grid_temp$summed_abundance)
    
    # Try moving each transplant to each empty cell
    transplant_cells <- which(grid_temp$Treatment == "4-years-old")
    empty_cells <- which(grid_temp$Treatment == "Bare")
    
    for(from in transplant_cells) {
      for(to in empty_cells) {
        test_grid <- grid_temp
        test_grid$Treatment[from] <- "Bare"
        test_grid$Treatment[to] <- "4-years-old"
        
        # Recalculate Area10m for affected cells
        radius <- 10  # 10m radius for Area10m calculation
        affected_cells <- st_is_within_distance(test_grid[c(from, to),], test_grid, radius)
        for(i in unique(unlist(affected_cells))) {
          neighbors <- st_is_within_distance(test_grid[i,], test_grid, radius)[[1]]
          test_grid$Area10m[i] <- sum(test_grid$Treatment[neighbors] == "4-years-old") * grid_area
        }
        
        # Calculate new abundance
        test_grid <- compute_summed_abundances(test_grid, gam_list, scalevals, ispp_obj)
        new_abundance <- sum(test_grid$summed_abundance)
        
        if(new_abundance > current_abundance) {
          grid_temp <- test_grid
          current_abundance <- new_abundance
          improved <- TRUE
          break
        }
      }
      if(improved) break
    }
  }
  
  list(abundance = current_abundance, grid = grid_temp)
}
