
library(future)
library(future.apply)

#### This part is unchanged - see comment for where changes start (at )
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



#######
## Changes start here ###
#######

#Optimize transplant locations using best guess with improved spatial handling
#The first version (my recommendation) includes diversification strategies that help escape local optima. 
#It sacrifices some short-term gains for potentially better global solutions by exploring more of the solution space.
#best_guess_optimize <- function(grid_data, n_transplants, update_area10m = TRUE, 
 #                               random_seed_pct = 0.1, random_jumps = TRUE, 
  #                              jump_frequency = 50, batch_size = 100) {
  #The second version is a pure greedy algorithm that always selects what looks best at each step. 
  #It's more deterministic but more likely to get trapped in local optima, especially in spatial problems where cluster effects are important.
  #best_guess_optimize <- function(grid_data, n_transplants, update_area10m = TRUE, 
  #                               random_seed_pct = 0, random_jumps = FALSE, 
  #                              jump_frequency = 100, batch_size = 100) {
  
  #Optimize transplant locations with better performance

#With parallel processing
best_guess_optimize <- function(grid_data, n_transplants, update_area10m = TRUE, 
                                random_seed_pct = 0.0, random_jumps = TRUE, 
                                jump_frequency = 100, batch_size = 30) {
  # Create a working copy of the grid data
  working_grid <- grid_data
  selected_cells <- c()
  best_abundance <- -Inf
  
  # Initialize progress tracking
  progress <- list()
  progress_counter <- 0
  
  # Optional: Random seed a percentage of cells first to improve starting configuration
  if (random_seed_pct > 0) {
    bare_indices <- which(working_grid$Treatment == "Bare")
    n_random_seed <- floor(n_transplants * random_seed_pct)
    
    # Randomly select cells to seed
    random_seed_cells <- sample(bare_indices, n_random_seed)
    
    # Add to selected cells
    selected_cells <- random_seed_cells
    
    # Update Treatment
    working_grid$Treatment[random_seed_cells] <- "4-years-old"
    working_grid$InSeagrass[random_seed_cells] <- TRUE
    
    # Update Area10m for all cells based on new seagrass configuration
    nseagrass <- st_intersects(buffer_input, working_grid[working_grid$InSeagrass,]) %>%
      lapply(length) %>% unlist()
    working_grid$Area10m <- nseagrass * grid_area_input
    
    # Calculate abundance after random seeding
    working_grid <- compute_summed_abundances(working_grid, gam_list, scalevals, ispp_obj, grid_area_ratio)
    
    cat("Randomly seeded", n_random_seed, "cells. Current abundance:", sum(working_grid$summed_abundance), "\n")
    
    # Add to progress
    progress_counter <- progress_counter + 1
    progress[[progress_counter]] <- list(
      n_transplants = length(selected_cells),
      abundance = sum(working_grid$summed_abundance),
      grid = working_grid
    )
  }
  
  # Counter for random jumps
  steps_since_last_jump <- 0
  
  # Create a spatial index for quick distance queries
  # This pre-calculation significantly speeds up spatial operations
  distance_cache <- list()
  
  # Pre-calculate cell-to-buffer relationships
  cat("Pre-calculating spatial relationships...\n")
  buffer_index <- st_intersects(buffer_input, grid_data)
  
  # OPTIMIZATION: Pre-calculate all distance relationships at once
  cat("Pre-calculating distance relationships (this might take a while but will speed up later steps)...\n")
  bare_indices <- which(working_grid$Treatment == "Bare")
  
  # Only calculate for a reasonable subset if there are too many
  if (length(bare_indices) > 1000) {
    cat("Large grid detected! Pre-calculating distances for the first 1000 bare cells...\n")
    precalc_indices <- sample(bare_indices, 1000)
  } else {
    precalc_indices <- bare_indices
  }
  
  # Create a list to store the pre-calculated distances
  cell_distances <- list()
  
  # Use mclapply from the parallel package if available
  if (requireNamespace("parallel", quietly = TRUE)) {
    cat("Using parallel processing for pre-calculation...\n")
    # Get number of cores (leave one for system)
    n_cores <- max(1, parallel::detectCores() - 1)
    cat("Using", n_cores, "cores\n")
    
    # Split the work
    chunks <- split(precalc_indices, cut(seq_along(precalc_indices), n_cores))
    
    # Process each chunk in parallel
    cell_distances_chunks <- parallel::mclapply(chunks, function(chunk_indices) {
      chunk_result <- list()
      for (idx in chunk_indices) {
        chunk_result[[as.character(idx)]] <- st_is_within_distance(working_grid[idx,], working_grid, 10)[[1]]
      }
      return(chunk_result)
    }, mc.cores = n_cores)
    
    # Combine results
    for (chunk in cell_distances_chunks) {
      cell_distances <- c(cell_distances, chunk)
    }
  } else {
    # Sequential fallback
    cat("Parallel package not available, using sequential pre-calculation...\n")
    pb <- txtProgressBar(min = 0, max = length(precalc_indices), style = 3)
    for (i in seq_along(precalc_indices)) {
      idx <- precalc_indices[i]
      cell_distances[[as.character(idx)]] <- st_is_within_distance(working_grid[idx,], working_grid, 10)[[1]]
      setTxtProgressBar(pb, i)
    }
    close(pb)
  }
  
  cat("Pre-calculation complete!\n")
  
  # Main optimization loop
  while(length(selected_cells) < n_transplants) {
    bare_indices <- which(working_grid$Treatment == "Bare")
    if(length(bare_indices) == 0) break
    
    # Progress reporting
    if(length(selected_cells) %% 10 == 0) {
      cat("Processing transplant", length(selected_cells), "of", n_transplants, "\n")
    }
    
    # Decide whether to make a random jump
    make_random_jump <- random_jumps && steps_since_last_jump >= jump_frequency
    
    if (make_random_jump) {
      # Random jump: Choose a random bare cell instead of the best one
      best_cell <- sample(bare_indices, 1)
      steps_since_last_jump <- 0
      cat("Making random jump at", length(selected_cells), "transplants\n")
    } else {
      # PERFORMANCE OPTIMIZATION: Use a small batch size
      # For very large grids, even smaller batches (10-20) can work well
      if (length(bare_indices) > batch_size) {
        candidate_indices <- sample(bare_indices, batch_size)
      } else {
        candidate_indices <- bare_indices
      }
      
      # Initialize vector to store marginal values for candidate cells
      marginal_values <- numeric(length(candidate_indices))
      
      # Get current abundance as baseline
      current_abundance <- sum(working_grid$summed_abundance)
      
      # PARALLEL OPTIMIZATION: Use parallel processing for this batch if available
      if (requireNamespace("parallel", quietly = TRUE) && length(candidate_indices) >= 4) {
        # Number of cores to use (leave one for system)
        n_cores <- min(parallel::detectCores() - 1, length(candidate_indices))
        
        # Prepare the data needed for each worker
        # We'll create a copy with only what's needed to reduce memory usage
        shared_data <- list(
          working_grid = working_grid,
          current_abundance = current_abundance,
          cell_distances = cell_distances,
          grid_area_input = grid_area_input,
          gam_list = gam_list,
          scalevals = scalevals,
          ispp_obj = ispp_obj,
          grid_area_ratio = grid_area_ratio
        )
        
        # Process each candidate cell in parallel
        results <- parallel::mclapply(seq_along(candidate_indices), function(i) {
          cell_index <- candidate_indices[i]
          w_grid <- shared_data$working_grid
          
          # Create a temporary copy of the treatment and InSeagrass status
          temp_treatment <- w_grid$Treatment
          temp_inseagrass <- w_grid$InSeagrass
          
          # Update for this candidate
          temp_treatment[cell_index] <- "4-years-old"
          temp_inseagrass[cell_index] <- TRUE
          
          # Create a temporary grid with these changes
          test_grid <- w_grid
          test_grid$Treatment <- temp_treatment
          test_grid$InSeagrass <- temp_inseagrass
          
          # Update Area10m for all cells
          test_grid$Area10m <- w_grid$Area10m  # Start with current values
          
          # Find cells that would be affected by this transplant (within 10m)
          cell_key <- as.character(cell_index)
          if (cell_key %in% names(shared_data$cell_distances)) {
            affected_cells <- shared_data$cell_distances[[cell_key]]
          } else {
            affected_cells <- st_is_within_distance(test_grid[cell_index,], test_grid, 10)[[1]]
          }
          
          # For affected cells, recalculate Area10m
          for (affected_idx in affected_cells) {
            # Get cells within 10m of the affected cell
            affected_key <- as.character(affected_idx)
            if (affected_key %in% names(shared_data$cell_distances)) {
              seagrass_cells <- shared_data$cell_distances[[affected_key]]
            } else {
              seagrass_cells <- st_is_within_distance(test_grid[affected_idx,], test_grid, 10)[[1]]
            }
            
            # Count seagrass points
            test_grid$Area10m[affected_idx] <- sum(temp_inseagrass[seagrass_cells]) * shared_data$grid_area_input
          }
          
          # Calculate new abundance with this cell transplanted
          test_grid <- compute_summed_abundances(test_grid, shared_data$gam_list, 
                                                 shared_data$scalevals, shared_data$ispp_obj, 
                                                 shared_data$grid_area_ratio)
          
          # Return marginal value
          sum(test_grid$summed_abundance) - shared_data$current_abundance
        }, mc.cores = n_cores)
        
        # Extract results
        marginal_values <- unlist(results)
      } else {
        # Sequential fallback
        for (i in seq_along(candidate_indices)) {
          cell_index <- candidate_indices[i]
          
          # Create a temporary copy of the treatment and InSeagrass status
          temp_treatment <- working_grid$Treatment
          temp_inseagrass <- working_grid$InSeagrass
          
          # Update for this candidate
          temp_treatment[cell_index] <- "4-years-old"
          temp_inseagrass[cell_index] <- TRUE
          
          # Create a temporary grid with these changes
          test_grid <- working_grid
          test_grid$Treatment <- temp_treatment
          test_grid$InSeagrass <- temp_inseagrass
          
          # Update Area10m for all cells
          test_grid$Area10m <- working_grid$Area10m  # Start with current values
          
          # Find cells that would be affected by this transplant (within 10m)
          cell_key <- as.character(cell_index)
          if (cell_key %in% names(cell_distances)) {
            affected_cells <- cell_distances[[cell_key]]
          } else {
            affected_cells <- st_is_within_distance(test_grid[cell_index,], test_grid, 10)[[1]]
            cell_distances[[cell_key]] <- affected_cells  # Cache for future use
          }
          
          # For affected cells, recalculate Area10m
          for (affected_idx in affected_cells) {
            # Get cells within 10m of the affected cell
            affected_key <- as.character(affected_idx)
            if (affected_key %in% names(cell_distances)) {
              seagrass_cells <- cell_distances[[affected_key]]
            } else {
              seagrass_cells <- st_is_within_distance(test_grid[affected_idx,], test_grid, 10)[[1]]
              cell_distances[[affected_key]] <- seagrass_cells  # Cache for future use
            }
            
            # Count seagrass points
            test_grid$Area10m[affected_idx] <- sum(temp_inseagrass[seagrass_cells]) * grid_area_input
          }
          
          # Calculate new abundance with this cell transplanted
          test_grid <- compute_summed_abundances(test_grid, gam_list, scalevals, ispp_obj, grid_area_ratio)
          
          # Marginal value is the improvement in total abundance
          marginal_values[i] <- sum(test_grid$summed_abundance) - current_abundance
        }
      }
      
      # Find cell with highest marginal value
      best_cell <- candidate_indices[which.max(marginal_values)]
      steps_since_last_jump <- steps_since_last_jump + 1
    }
    
    # Add selected cell to the transplant list
    selected_cells <- c(selected_cells, best_cell)
    
    # Update Treatment
    working_grid$Treatment[best_cell] <- "4-years-old"
    working_grid$InSeagrass[best_cell] <- TRUE
    
    # OPTIMIZATION: Only update Area10m for cells within 10m of the transplant
    cell_key <- as.character(best_cell)
    if (cell_key %in% names(cell_distances)) {
      affected_cells <- cell_distances[[cell_key]]
    } else {
      affected_cells <- st_is_within_distance(working_grid[best_cell,], working_grid, 10)[[1]]
      cell_distances[[cell_key]] <- affected_cells
    }
    
    # Update Area10m only for affected cells
    for (affected_idx in affected_cells) {
      affected_key <- as.character(affected_idx)
      if (affected_key %in% names(cell_distances)) {
        seagrass_cells <- cell_distances[[affected_key]]
      } else {
        seagrass_cells <- st_is_within_distance(working_grid[affected_idx,], working_grid, 10)[[1]]
        cell_distances[[affected_key]] <- seagrass_cells
      }
      
      working_grid$Area10m[affected_idx] <- sum(working_grid$InSeagrass[seagrass_cells]) * grid_area_input
    }
    
    # Calculate new abundance after this transplant
    working_grid <- compute_summed_abundances(working_grid, gam_list, scalevals, ispp_obj, grid_area_ratio)
    
    # Save progress every 10 transplants
    if(length(selected_cells) %% 10 == 0) {
      progress_counter <- progress_counter + 1
      current_abundance <- sum(working_grid$summed_abundance)
      
      cat("Transplants:", length(selected_cells), "Abundance:", current_abundance, "\n")
      
      if(current_abundance > best_abundance) {
        best_abundance <- current_abundance
        best_grid <- working_grid
      }
      
      progress[[progress_counter]] <- list(
        n_transplants = length(selected_cells),
        abundance = current_abundance,
        grid = working_grid
      )
    }
  }
  
  # Final abundance calculation to ensure consistency
  final_abundance <- sum(working_grid$summed_abundance)
  
  return(list(
    final_grid = working_grid,
    progress = progress,
    final_abundance = final_abundance
  ))
}