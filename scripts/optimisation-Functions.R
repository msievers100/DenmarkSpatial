

compute_summed_abundances <- function(df, gam_list, scalevals, ispp_obj, grid_area_ratio) {
  predictions <- lapply(ispp_obj, function(ispp) {
    (predict(gam_list[[ispp]], newdata = df, type = "response")* grid_area_ratio) / scalevals[ispp]
  })
  df$summed_abundance <- Reduce(`+`, predictions) / length(ispp_obj)
  return(df)
}


optimize_transplants <- function(n_transplants, grid2_cropped, gam_list, 
                                 scalevals, ispp_obj, n_starts = 100) {
  
  # Initialize best solution trackers
  best_abundance <- -Inf
  best_grid <- grid2_cropped
  
  # First run with marginal value-based initialization
  results <- list()
  results[[1]] <- run_optimization(1, grid2_cropped, initial = TRUE)
  
  # Run remaining starts with random initialization
  for(i in 1:n_starts) {
    results[[i]] <- run_optimization(i, grid2_cropped, initial = FALSE)
  }
  
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
run_optimization <- function(seed, grid2_cropped, initial = FALSE) {
  set.seed(seed)
  
  # Create working copy of grid
  grid_temp <- grid2_cropped
  
  if (!initial) {
    # Initialize with random transplants
    bare_cells <- which(grid_temp$Treatment == "Bare")
    initial_cells <- sample(bare_cells, n_transplants)
    grid_temp$Treatment[initial_cells] <- "4-years-old"
  } else {
    # Use marginal value-based initialization
    bare_cells <- which(grid_temp$Treatment == "Bare")
    marginal_values <- grid_temp$marginal_value[bare_cells]
    initial_cells <- bare_cells[order(marginal_values, decreasing = TRUE)][1:n_transplants]
    grid_temp$Treatment[initial_cells] <- "4-years-old"
  }
  
  improved <- TRUE
  k=0
  while(improved) {
    improved <- FALSE
    k <- k+1
    cat("Iteration:", k, "\n")
    # Calculate current abundance
    grid_temp <- compute_summed_abundances(grid_temp, gam_list, scalevals, ispp_obj, grid_area_ratio)
    current_abundance <- sum(grid_temp$summed_abundance)
    print(current_abundance)
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
        test_grid <- compute_summed_abundances(test_grid, gam_list, scalevals, ispp_obj, grid_area_ratio)
        new_abundance <- sum(test_grid$summed_abundance)
        
        if(round(new_abundance,1) > round(current_abundance,1)) {
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




#Optimize transplant locations using best guess
# update_area10m = TRUE will update the Area10m for each cell after each transplant then recalculate marginal values
# update_area10m = FALSE will not update the Area10m for each cell after each transplant and therefore the 
# order of priorities is fixed based on initial marginal values
# both optoins will update area10m for calculation of predicted abundances. 
#
best_guess_optimize <- function(grid_data, n_transplants, update_area10m = TRUE) {
  # Create a working copy of the grid data
  working_grid <- grid_data
  selected_cells <- c()
  best_abundance <- -Inf
  # Initialize progress tracking
  progress <- list()
  progress_counter <- 0
  
  while(length(selected_cells) < n_transplants) {
    # Find cell with highest marginal value among bare cells
    bare_indices <- which(working_grid$Treatment == "Bare")
    if(length(bare_indices) == 0) break
    
    marginal_values <- working_grid$marginal_value[bare_indices]
    best_cell <- bare_indices[which.max(marginal_values)]
    
    # Add to selected cells
    selected_cells <- c(selected_cells, best_cell)
    
    # Update Treatment
    working_grid$Treatment[best_cell] <- "4-years-old"
    # Update Area10m for all cells within 10m radius
    
    # ------------------------
    # Slow method, but same as scenario-functions.R
    working_grid$InSeagrass[working_grid$Treatment == "4-years-old"] <- TRUE
    nseagrass <- st_intersects(buffer_input, working_grid[working_grid$InSeagrass,]) %>%
      lapply(length) %>% unlist()
    working_grid$Area10m <- nseagrass * grid_area_input
    
    # ------------------------
    # Fast method, but different to scenario-functions.R
    # affected_cells <- st_is_within_distance(working_grid[best_cell,], working_grid, 10)[[1]]
    
    # Add an extra cell to the area10m for each affected cell
    # for(i in affected_cells) {
    #   working_grid$Area10m[i] <- working_grid$Area10m[i] + grid_area
    # }
    # ------------------------
#plot to check 
    # tm_shape(working_grid) + 
    #   tm_dots() + 
    #   tm_shape(working_grid[affected_cells,]) + 
    #   tm_dots(col = "red") + 
    #   tm_shape(working_grid[best_cell,]) + 
    #   tm_dots(col = "green")
    

    #Calculate new value
    working_grid <- compute_summed_abundances(working_grid, gam_list, scalevals, ispp_obj, grid_area_ratio)
    
    #Calculate new marginal values - what if each bare cell now had a transplant? 
    working_grid_treated <- working_grid
    working_grid_treated$Treatment[working_grid_treated$Treatment == "Bare"] <- "4-years-old"
    #Add to area10m just for the cell that is treated (so not accounting for neighours as in maringal values
    # we are assuming they are staying the same)
    working_grid_treated$Area10m[working_grid_treated$Treatment == "Bare"] <- 
      working_grid_treated$Area10m[working_grid_treated$Treatment == "Bare"] + grid_area
    working_grid_treated <- compute_summed_abundances(working_grid_treated, gam_list, scalevals, ispp_obj, grid_area_ratio)

    # Marginal value is the difference in summed abundance and is used to rank order transplants grids for priority
    if (update_area10m){
      working_grid$marginal_value <- working_grid_treated$summed_abundance - working_grid$summed_abundance
    }
    
    
    # Save progress every 10 transplants
    if(length(selected_cells) %% 10 == 0) {
      progress_counter <- progress_counter + 1
      working_grid <- compute_summed_abundances(working_grid, gam_list, scalevals, ispp_obj, grid_area_ratio)
      
      if(sum(working_grid$summed_abundance) > best_abundance) {
        best_abundance <- sum(working_grid$summed_abundance)
        best_grid <- working_grid
      }
      
      progress[[progress_counter]] <- list(
        n_transplants = length(selected_cells),
        abundance = sum(working_grid$summed_abundance),
        grid = working_grid, 
        best_grid = best_grid,
        best_abundance = best_abundance
      )
    }
  }
  
  return(list(
    final_grid = working_grid,
    progress = progress
  ))
}
