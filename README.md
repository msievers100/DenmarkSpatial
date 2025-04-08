# Readme for the project "Optimizing multi-habitat restoration for seascape-scale animal benefits in coastal ecosystems "

Sievers, Brown 

## Summary

Multi-habitat restoration has been advocated to best enhance biodiversity and ecosystem services in degraded environments (Vozzo et al. 2023), and the functional role of animals is increasingly recognized as a critical component of successful restoration (Sievers et al. 2022). Empirical studies examining how animals respond to different spatial configurations of restored habitats, however, remain scarce. We collect and analyse animal data from one of the world’s largest multi-habitat coastal restoration sites comprising restored seagrass (Zostera marina), boulder reefs, and mussel reefs (Mytilus edulis). We use distinct small-scale spatial patterns in population abundances to model outcomes for five hypothetical configurations of restored seagrass and produce a series of optimizations, demonstrating it is practical to configure restoration to optimize biodiversity objectives. Surveys of mobile fauna – primarily fishes, crustaceans and gastropods – across 132 sites in bare areas and natural and restored seagrass showed that diversity was highest near mussel reefs and in more mature seagrass, while individual species responded in complex ways to restored seagrass age and seascape characteristics. Species-specific responses translated to variable outcomes across restoration scenarios and optimizations. For example, pipefish (Sygnathus rostellatus) abundance was highest under random site selection, while flatfish abundance was greatest for a planting scenario starting near boulder reefs. The optimal number and arrangement of transplants varied significantly depending on the target species or species group, but a common arrangement was patchy planting, highlighting the importance of not homogenizing seascapes. Our approach provides a practical framework for incorporating animal monitoring data into restoration planning, helping practitioners design and optimize spatial planting configurations.

## For community metrics (richness, Shannon diversity, Evenness)
## Run SpeciesCommGAMAnalysis.R, using UTMdat.csv
## For dbRDA plot, run CommunityMV.R, using UTMdat1MV.csv
#Used columns are defined below:
UTM_E - UTM east
UTM_N - UTM north
DistBouldEdge - Distance to boulder reef edge in meters
DistMussEdge - Distance to mussel reef edge in meters
Area10m - Area of seagrass within a 10m radius around site
Perim10m - Perimeter of seagrass within a 10m radius around site
PerAreaRatio10m - Perimeter to area ratio within a 10m radius around site
RatioLog - Log of PerAreaRatio10m
Treatment - Site type: Natural seagrass, Bare sand, Restored seagrass (1, 2, 3 or 4-years-old)
Large.shrimp - Abundance of Baltic prawns
Lesser.pipefish - Abundance of pipefish
Other.shrimp - Abundance of brown shrimp (Crangon prawn)
Smooth.flatfish - Abundance of flatfish
Herb.snail.ab - Abundance of periwinkle
Whelk.ab - Abundance of whelk






## Spatial scenarios and optimization below ##


## Directory structure

### Data (`/data`)
- Survey data: `datCABMS.csv`
- Spatial data files:
  - Boulder reef data: `BoulderReefPoints.*`, `BoulderReefPoly.*`, `NewBould.*`
  - Mussel reef data: `MusselReef.*`
  - Seagrass data: `SeagrassRestNat.*`
  - Site locations: `Sites.*`
  - Additional data: `Seagrass polygon area.xlsx`

### Outputs (`/outputs`)
- Model outputs (`.rds` files):
  - `boulder.rds`, `buffers.rds`, `gam-selected-species.rds`
  - `grid_area.rds`, `grid2_cropped.rds`, `mussel_dissolved.rds`
  - `newbould.rds`, `ngrids_no_seagrass.rds`
- Analysis results: `deviance_Table.csv`
- Visualization outputs:
  - Species plots (`/plots`): Visreg plots for various species
  - RShiny data (`/RShinyrds`): Progress and results files for different scenarios
  - CSV outputs (`/csv`)

### Scripts (`/scripts`)

#### Main Workflow Scripts
1. `1_SpatialDataPrep.R` - Prepare spatial data for analysis
2. `2_GAMS-looped.R` - Run GAMS for each species in a loop
3. `3_Source-greedy-new.R` - Source the greedy optimisation script with parameters of choice
4. `4_Scenarios-looped.R` - Run restoration scenarios for each species, 
5. `5_plot-scenario-objectives-new.R` - Plot the objective functions for each species including random and distance based scenarios
6. `6_plot-greedy.R` - Plot greedy functions/maps

#### Analysis and Comparison Scripts
- `Greedy-optimisation-new.R` - Run the greedy optimisation, sourced from`3_Source-greedy.R`
- `GAMs.R` - Test generalized additive models for species abundance
- `Compare-Objectives.R`, `Compare-Solutions.R`
- `Plot-Comparison.R`
- `calculate-objective-functions.R`
- `Scenarios.R` - Run restoration scenarios for one species

#### Utility Scripts
- `Scenario-Functions.R` - Functions for running restoration scenarios
- `optimisation-functions-new2-parallel.R` - Functions for the greedy optimisation


## Optimisation methods

We optimised the location of a fixed number of transplants to maximize the combined abundance of three key species: crangon prawn, baltic prawn and flatfish. The objective function was:

$$
\text{Objective Function} = \sum{j=1}^{K} \sum_{i=1}^{N} (A_{j,i}/s_{j})
\text{Subject to:} \sum_{i=1}^{n} T_i = T
$$

where \( A_{j,i} \) was the abundances of the K species (crangon prawn, baltic prawn, and flatfish) at grid cell \( i \), respectively, and \( N \) was the total number of grid cells. the \s_{j} where scaling factors for each species, which rescaled their abundances by the standard deviations from the data.  \( T_i \) was an indicator (0/1) of transplanting at grid cell \( i \), and \( T \) was the target number of transplants.

We used a combination of heuristics to identify the optimal transplant locations. First, created a greedy heuristic. The heuristic iteratively applied these steps, starting with a seascape with no transplanted seagrass:
1. Calculate the benefit of replacing each individual bare cell with a seagrass transplant, calculated as difference between the objective function with that transplant versus without. 
2. Rank order all cells from high to low benefit. 
3. Implement the transplant in the highest benefit cell. 
4. Repeat steps 1-3 until the target number of transplants is reached. 

We also implemented a greedy heuristic that started at a random allocation of the target number of transplants. It then iteratively swapped a transplant cell with a bare sediment cell, accepting swaps that improved the objective function. This algorithm never found a better solution than the greedy heuristic. So results are presented for the greedy heuristic. 


# Seagrass Transplant Optimization Improvements

This document explains the optimizations made to the seagrass transplant location optimization algorithm and how to use the improved functions.

## Background

The original optimization algorithm (`best_guess_optimize`) was designed to find optimal locations for seagrass transplants to maximize fish abundance. However, it had performance issues and sometimes performed worse than random placement.

## Key Issues Identified

1. **Incorrect Marginal Value Calculation**: The original algorithm calculated marginal values incorrectly by:
   - Creating a grid where ALL bare cells were treated as transplants
   - Using the difference between this all-treated scenario and current abundance
   - This didn't capture the true individual contribution of each cell

2. **Performance Bottlenecks**:
   - Expensive spatial operations (`st_is_within_distance`) called repeatedly
   - Recalculation of Area10m for all cells after each transplant
   - Inefficient grid copying and state management

## Optimizations Implemented

The new `best_guess_optimize_fast` function includes the following improvements:

1. **Pre-calculated Spatial Relationships**:
   - One-time calculation of all cells within 10m of each cell
   - Stored in a lookup table for fast access
   - Eliminates repeated expensive spatial operations

2. **Direct Area10m Updates**:
   - Instead of recounting all seagrass points within 10m
   - Simply adds `grid_area_input` to Area10m for affected cells
   - Much faster and mathematically equivalent

3. **Parallel Processing**:
   - Evaluates candidate cells in parallel
   - Utilizes multiple CPU cores
   - Significantly speeds up batch evaluation

## Performance Comparison

The benchmark script (`benchmark-optimize.R`) compares the performance of the original and optimized functions. Typical results show:

- **Speed Improvement**: 5-10x faster execution
- **Result Quality**: Equivalent or better fish abundance
- **Memory Usage**: More efficient, especially for large grids

## How to Use

### Basic Usage

```r
# Source the optimized functions
source("scripts/optimisation-Functions-optimized.R")

# Run optimization
result <- best_guess_optimize_fast(
  grid_data = grid_data,
  n_transplants = 100,
  update_area10m = TRUE,
  random_jumps = TRUE,
  batch_size = 20,
  parallel = TRUE
)

# Access results
final_grid <- result$final_grid
final_abundance <- result$final_abundance
progress <- result$progress
```

### Parameters

- `grid_data`: Spatial grid data with Treatment and other attributes
- `n_transplants`: Number of transplants to place
- `update_area10m`: Whether to update Area10m after each transplant
- `random_seed_pct`: Percentage of transplants to randomly seed at start
- `random_jumps`: Whether to make random jumps to escape local optima
- `jump_frequency`: How often to make random jumps
- `batch_size`: Number of candidate cells to evaluate at each step
- `parallel`: Whether to use parallel processing
- `n_cores`: Number of cores to use (NULL = auto-detect)

### Example Script

See `example-usage-optimized.R` for a complete example of how to use the optimized function, including:
- Loading data
- Setting parameters
- Running the optimization
- Plotting progress
- Saving results
- Creating maps of optimized transplant locations

## Benchmarking

To compare the performance of the original and optimized functions:

```r
# Run the benchmark script
source("scripts/benchmark-optimize.R")
```

This will:
1. Run both functions with the same parameters
2. Measure execution time
3. Compare results for consistency
4. Calculate speedup factor
5. Save benchmark results

## Technical Details

### Spatial Relationship Pre-calculation

The most significant optimization is pre-calculating all spatial relationships:

```r
# Create a complete distance matrix (sparse representation)
n_cells <- nrow(grid_data)
distance_matrix <- list()

for (i in 1:n_cells) {
  # Find all cells within 10m of cell i
  affected_cells <- st_is_within_distance(grid_data[i,], grid_data, 10)[[1]]
  distance_matrix[[i]] <- affected_cells
}
```

### Direct Area10m Updates

Instead of recounting all seagrass points:

```r
# Original (slow) approach
seagrass_count <- sum(temp_inseagrass[st_is_within_distance(
  test_grid[affected_idx,], test_grid, 10)[[1]]])
test_grid$Area10m[affected_idx] <- seagrass_count * grid_area_input

# Optimized (fast) approach
test_grid$Area10m[affected_cells] <- test_grid$Area10m[affected_cells] + grid_area_input
```

### Parallel Processing

Using the parallel package for candidate evaluation:

```r
# Create a cluster
cl <- parallel::makeCluster(n_cores)

# Export necessary objects to the cluster
parallel::clusterExport(cl, varlist = c("working_grid", "current_abundance", 
                                       "distance_matrix", "grid_area_input",
                                       "gam_list", "scalevals", "ispp_obj", 
                                       "grid_area_ratio", "compute_summed_abundances"), 
                       envir = environment())

# Parallel evaluation of candidates
marginal_values <- parallel::parSapply(cl, candidate_indices, function(cell_index) {
  # Evaluation code...
})

# Stop the cluster
parallel::stopCluster(cl)
```

## Conclusion

The optimized algorithm correctly calculates the true marginal value of each potential transplant location and runs significantly faster than the original version. This allows for more extensive exploration of the solution space and better final results.
