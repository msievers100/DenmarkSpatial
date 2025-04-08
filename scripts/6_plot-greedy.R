####

library(tmap)
library(sf)
library(viridis)
library(dplyr)

# Load the optimization results
greedy_results <- readRDS("outputs/greedy-results/best_guess_results_ispp_pipe_2025-03-13.rds")
#greedy_results <- readRDS("outputs/greedy-results/best_guess_results_ispp_fisheries_2025-03-13.rds")
#greedy_results <- readRDS("outputs/greedy-results/best_guess_results_ispp_SGfunction_2025-03-13.rds")
#greedy_results <- readRDS("outputs/greedy-results/best_guess_results_ispp_all_2025-03-13.rds")

#greedy_results <- readRDS("outputs/greedy-results/best_guess_results_baltic_2025-03-13.rds")
#greedy_results <- readRDS("outputs/greedy-results/best_guess_results_browns_2025-03-13.rds")
#greedy_results <- readRDS("outputs/greedy-results/best_guess_results_flatfish_2025-03-13.rds")
#greedy_results <- readRDS("outputs/greedy-results/best_guess_results_peri_2025-03-13.rds")
#greedy_results <- readRDS("outputs/greedy-results/best_guess_results_whelk_2025-03-13.rds")
#greedy_results <- readRDS("outputs/greedy-results/best_guess_results_Figure2_2025-03-13.rds")

# Extract the progress data to find the point of maximum abundance
progress_data <- do.call(rbind, lapply(greedy_results$progress, function(x) {
  data.frame(n_transplants = x$n_transplants, abundance = x$abundance)
}))

# Find the transplant count with maximum abundance
optimal_transplants <- progress_data$n_transplants[which.max(progress_data$abundance)]
max_abundance <- max(progress_data$abundance)

cat("Maximum abundance of", max_abundance, "at", optimal_transplants, "transplants\n")

# Find the progress entry closest to this optimal transplant count
optimal_index <- which.min(abs(sapply(greedy_results$progress, function(x) x$n_transplants - optimal_transplants)))
optimal_grid <- greedy_results$progress[[optimal_index]]$grid

# Verify we have the right grid
cat("Selected grid has", sum(optimal_grid$Treatment == "4-years-old"), "transplanted sites\n")

# Create new Treatment factor with only the levels we want
optimal_grid$Treatment <- factor(
  as.character(optimal_grid$Treatment),
  levels = c("Bare", "4-years-old", "Natural")
)

# Create a map showing the optimized transplant locations
#tbest <- tm_shape(optimal_grid) +
#  tm_dots(col = "Treatment", 
#          palette = c("Bare" = "khaki1", "4-years-old" = "darkorchid1", "Natural" = "darkolivegreen"),
#          size = 0.22,
#          title = "Optimal Abundance\nTransplant Locations") +
#  tm_layout(legend.outside = TRUE,
#            main.title = paste0("Maximum Abundance at ", optimal_transplants, " Sites"),
#            main.title.position = "center",
#            frame = FALSE) #+
  #tm_scale_bar(position = c("left", "bottom"),
   #            text.size = 0.7) 

# Display the map
#tbest

# Save the map
#tmap_save(tbest, filename = paste0("outputs/plots/optimal-abundance-map-Figure2", Sys.Date(), ".png"),
#          width = 5, height = 3.5, dpi = 300)

#Adding mussel and boulder reefs

mussel_dissolved <- readRDS("outputs/mussel_dissolved.rds")
bould <- readRDS("outputs/boulder.rds")

# First get the bounding box of the optimal_grid layer
bbox_expanded <- bbox_optimal
bbox_expanded[1] <- bbox_expanded[1] - 20  # xmin - 50m
bbox_expanded[2] <- bbox_expanded[2] - 10  # ymin - 50m
bbox_expanded[3] <- bbox_expanded[3] + 50  # xmax + 50m
bbox_expanded[4] <- bbox_expanded[4] + 10  # ymax + 50m

tbest2 <- tm_shape(mussel_dissolved, bbox = bbox_expanded) +
  tm_fill(col = "#990000", alpha = 0.6) +
  
  tm_shape(bould) +
  tm_fill(col = "#FF7F00", alpha = 0.6) +
  
  tm_shape(optimal_grid) +
  tm_dots(col = "Treatment", 
          palette = c("Bare" = "khaki1", "4-years-old" = "darkorchid1", "Natural" = "darkolivegreen"),
          size = 0.22,
          title = "Optimal Abundance\nTransplant Locations") +
  
  tm_layout(legend.outside = TRUE,
            main.title = paste0("Maximum Abundance at ", optimal_transplants, " Sites"),
            main.title.position = "center") + #,
            #frame = FALSE) +
  tm_scale_bar(position = c("LEFT", "bottom"),  # Using uppercase LEFT to position it farther left
               text.size = 0.7)
tbest2
tmap_save(tbest2, filename = paste0("outputs/plots/optimal-abundance-map2-SGfunction", Sys.Date(), ".png"),
          width = 5, height = 3.5, dpi = 300)


### Top 10% of sites
# This is approximated by the maximum number of transplants in the progress data
total_available_sites <- max(progress_data$n_transplants)

# Calculate 10% of the available sites
ten_percent_sites <- round(total_available_sites * 0.1)

cat("Total available sites:", total_available_sites, "\n")
cat("10% of available sites:", ten_percent_sites, "\n")

# Find the progress entry closest to the 10% target
ten_percent_index <- which.min(abs(sapply(greedy_results$progress, function(x) x$n_transplants - ten_percent_sites)))
ten_percent_grid <- greedy_results$progress[[ten_percent_index]]$grid

# Verify we have the right grid
cat("Selected grid has", sum(ten_percent_grid$Treatment == "4-years-old"), "transplanted sites\n")

# Create new Treatment factor with only the levels we want
ten_percent_grid$Treatment <- factor(
  as.character(ten_percent_grid$Treatment),
  levels = c("Bare", "4-years-old", "Natural")
)

# Create a map showing the 10% restoration locations
tmap_10percent <- tm_shape(ten_percent_grid) +
  tm_dots(col = "Treatment", 
          palette = c("Bare" = "khaki1", "4-years-old" = "darkorchid1", "Natural" = "darkolivegreen"),
          size = 0.1,
          title = "Treatment Type") +
  tm_layout(legend.outside = TRUE,
            main.title = paste0("10% Restoration (", ten_percent_sites, " Sites)"),
            main.title.position = "center")

# Display the map
tmap_10percent

# Save the map
tmap_save(tmap_10percent, filename = paste0("outputs/plots/10-percent-restoration-map-fisheries", Sys.Date(), ".png"))#,
          #width = 8, height = 6, dpi = 300)











#### OLD code below. 


#saveRDS(best_guess_results, file = paste0("outputs/greedy-results/best_guess_results_",obj_name,"_", run_date,".rds"))

greedy_results <- readRDS("outputs/greedy-results/best_guess_results_baltic_2025-01-28.rds")
summary(greedy_results)
head(greedy_results)
#I don't know how to make maps of the optimisation for a given N transplants (which will be max abundance)

p_progress <- ggplot(progress_df) +
  aes(x = n_transplants, y = abundance) +
  geom_line() +
  #Compare to the greedy without updating
  # geom_point(aes(x = n_transplants[1], y = simple_greedy_best))
  labs(x = "Number of Transplants", y = "Total Abundance",
       title = "Optimization Progress") +
  theme_classic() #+ geom_vline(xintercept = 1810, linetype = "dashed")
p_progress
#ggsave("outputs/plots/optimization_progress-WithArea3973-Pipefish.png", p_progress, width = 5, height = 4)

## 3973 total
## 3307 bare (10% of this is 331)

progress_df$n_transplants[which.max(progress_df$abundance)]
#For fisheries spp. [1] 1680 - use for n_transplants to map optimal arrangement under unlimited resources
#For All 6 spp. [1] 930 - use for n_transplants to map optimal arrangement under unlimited resources
#For SG spp. [1] 2770 - use for n_transplants to map optimal arrangement under unlimited resources
#For Pipefish [1] 1810 - use for n_transplants to map optimal arrangement under unlimited resources
#Can use these values for i, below

# Visualize best guess placement for max transplants
tbest <- tm_shape(best_guess_grid) +
  tm_dots(col = "Treatment", 
          palette = c("Bare" = "khaki1", "4-years-old" = "darkorchid1", "Natural" = "darkolivegreen"),
          size = 0.1,
          title = "Best Guess\nTransplant Locations") +
  tm_layout(legend.outside = TRUE)
tbest






#
# Calculate summed abundance
#
grid_temp <- grid2_cropped

bare_cells <- which(grid_temp$Treatment == "Bare")
marginal_values <- grid_temp$marginal_value[bare_cells]
initial_cells <- bare_cells[order(marginal_values, decreasing = TRUE)][1:n_transplants]
grid_temp$Treatment[initial_cells] <- "4-years-old"

# update area10m for all cells within 10m radius
affected_cells <- st_is_within_distance(grid_temp[initial_cells,], grid_temp, 10)[[1]]
for(i in affected_cells) {
  grid_temp$Area10m[i] <- grid_temp$Area10m[i] + grid_area
}

grid_temp <- compute_summed_abundances(grid_temp, gam_list, scalevals, ispp_obj)

simple_greedy_best <- sum(grid_temp$summed_abundance)
##'simple_greedy_best' is just picking the top N sites
##'the best guess is picking the top site, recalculating Area10m, re-ranking, then picking the next top site

tm_shape(grid_temp) +
  tm_dots(col = "Treatment", 
          palette = c("Bare" = "khaki1", "4-years-old" = "darkorchid1", "Natural" = "darkolivegreen"),
          size = 0.1,
          title = "Best Guess\nTransplant Locations") +
  tm_layout(legend.outside = TRUE)


###################################
# STEP 4: Create a best guess optimisation
###################################

# Run the best guess optimization
best_guess_results <- best_guess_optimize(grid2_cropped, n_transplants)
best_guess_grid <- best_guess_results$final_grid

# Plot progress
progress_df <- do.call(rbind, lapply(best_guess_results$progress, function(x) {
  data.frame(n_transplants = x$n_transplants, abundance = x$abundance)
}))

p_progress <- ggplot(progress_df) +
  aes(x = n_transplants, y = abundance) +
  geom_line() +
  #Compare to the greedy without updating
  # geom_point(aes(x = n_transplants[1], y = simple_greedy_best))
  labs(x = "Number of Transplants", y = "Total Abundance",
       title = "Optimization Progress") +
  theme_classic() #+ geom_vline(xintercept = 1810, linetype = "dashed")
p_progress
#ggsave("outputs/plots/optimization_progress-WithArea3973-Pipefish.png", p_progress, width = 5, height = 4)

## 3973 total
## 3307 bare (10% of this is 331)

progress_df$n_transplants[which.max(progress_df$abundance)]
#For fisheries spp. [1] 1680 - use for n_transplants to map optimal arrangement under unlimited resources
#For All 6 spp. [1] 930 - use for n_transplants to map optimal arrangement under unlimited resources
#For SG spp. [1] 2770 - use for n_transplants to map optimal arrangement under unlimited resources
#For Pipefish [1] 1810 - use for n_transplants to map optimal arrangement under unlimited resources
#Can use these values for i, below

# Visualize best guess placement for max transplants
tbest <- tm_shape(best_guess_grid) +
  tm_dots(col = "Treatment", 
          palette = c("Bare" = "khaki1", "4-years-old" = "darkorchid1", "Natural" = "darkolivegreen"),
          size = 0.1,
          title = "Best Guess\nTransplant Locations") +
  tm_layout(legend.outside = TRUE)
tbest

#tmap_save(tbest, filename = paste0("outputs/plots/marginal-values-WithArea10percent-Pipefish", run_date, ".png"))



# Calculate effectiveness of best guess
best_guess_grid <- compute_summed_abundances(best_guess_grid, gam_list, scalevals, ispp_obj)
best_guess_abundance <- sum(best_guess_grid$summed_abundance)
print(paste("Best guess total abundance:", best_guess_abundance))
simple_greedy_best


#simple_greedy_best for Pipefish = 1459; 'best guess' for Pipefish = 1269.

#
#Visualize best guess placement for any number of transplants
# 
i <- 10
print(paste("This is for n_transplants = ", best_guess_results[[2]][[i]]$n_transplants))
my_grid <- best_guess_results[[2]][[i]]$grid
t1 <- tm_shape(my_grid) +
  tm_dots(col = "Treatment", 
          palette = c("Bare" = "khaki1", "4-years-old" = "darkorchid1", "Natural" = "darkolivegreen"),
          size = 0.1,
          title = "Best Guess\nTransplant Locations") +
  tm_layout(legend.outside = TRUE)
t1

#### Plotting all curves on the one plot ####
# Note, csv's for the four objectives saved in code below, in the RShiny section
pipefish <- read.csv("outputs/csv/progress_df_pipefish.csv")
all <- read.csv("outputs/csv/progress_df_all.csv")
fisheries <- read.csv("outputs/csv/progress_df_fisheries.csv")
sgfunction <- read.csv("outputs/csv/progress_df_sgfunction.csv")

combined_data <- bind_rows(
  mutate(pipefish, dataset = "Lesser pipefish"),
  mutate(all, dataset = "All species"),
  mutate(fisheries, dataset = "Fisheries spp."),
  mutate(sgfunction, dataset = "Seagrass specialists")
)

ggplot(all, aes(x = n_transplants, y = abundance)) +
  geom_line(size = 1.2) +  # Increased line thickness
  labs(x = "Number of Transplants",
       y = "Scaled abundance") +
  theme(legend.position = "right")

# First transform the data to create percentage values
combined_data_pct <- combined_data %>%
  group_by(dataset) %>%
  mutate(abundance_pct = (abundance / max(abundance) * 100)) %>%
  ungroup()


# Create the plot with thicker lines
progressplot <- ggplot(combined_data_pct, aes(x = n_transplants, y = abundance_pct, color = dataset)) +
  geom_line(size = 1.2) +  # Increased line thickness
  scale_color_manual(values = c("Lesser pipefish" = "#8884d8", 
                                "All species" = "#82ca9d",
                                "Fisheries spp." = "#ffc658",
                                "Seagrass specialists" = "#ff7300")) +
  geom_vline(xintercept = 2770, linetype = "dashed", colour = "#ff7300", size = 1) +
  geom_vline(xintercept = 930, linetype = "dashed", colour = "#82ca9d", size = 1) +
  geom_vline(xintercept = 1680, linetype = "dashed", colour = "#ffc658", size = 1) +
  geom_vline(xintercept = 1810, linetype = "dashed", colour = "#8884d8", size = 1) +
  labs(x = "Number of Transplants",
       y = "Scaled abundance (% of maximum)") +
  scale_y_continuous(labels = function(x) paste0(x, "%")) +  # Add % to y-axis labels
  theme(legend.position = "right")

#ggsave("outputs/plots/progress_comparison.png", progressplot, width = 4, height = 4, dpi = 300)


#### Mapping multi objectives at once ####
# Note, rds for the four objectives saved in code below, in the RShiny section
pipe_results_10percent <- readRDS("outputs/RShinyrds/pipe_results_10percent.rds")
fisheries_results_10percent <- readRDS("outputs/RShinyrds/fisheries_results_10percent.rds")
sgfunction_results_10percent <- readRDS("outputs/RShinyrds/sgfunction_results_10percent.rds")
all_results_10percent <- readRDS("outputs/RShinyrds/all_results_10percent.rds")


# Keep all sites from pipefish dataset
pipefish_sites <- pipe_results_10percent$final_grid
fisheries_sites <- fisheries_results_10percent$final_grid

# For 4-years-old sites, update categories
pipefish_4yr_indices <- which(pipefish_sites$Treatment == "4-years-old")
fisheries_4yr_indices <- which(fisheries_sites$Treatment == "4-years-old")

# Initialize category column for both datasets
pipefish_sites$category <- pipefish_sites$Treatment
fisheries_sites$category <- fisheries_sites$Treatment

# Update 4-years-old categories
pipefish_sites$category[pipefish_4yr_indices] <- "4-years-old (Pipefish Only)"
fisheries_sites$category[fisheries_4yr_indices] <- "4-years-old (Fisheries Only)"

# Find overlapping 4-year-old sites
matches <- st_equals(
  pipefish_sites[pipefish_4yr_indices,], 
  fisheries_sites[fisheries_4yr_indices,]
)

# Update categories for overlapping points
overlapping_indices <- which(sapply(matches, length) > 0)
if(length(overlapping_indices) > 0) {
  pipefish_sites$category[pipefish_4yr_indices[overlapping_indices]] <- "4-years-old (Both)"
}

# Get fisheries-only points
fisheries_4yr <- fisheries_sites[fisheries_4yr_indices,]
matches_all <- st_equals(fisheries_4yr, pipefish_sites[pipefish_4yr_indices,])
non_matching_indices <- which(sapply(matches_all, length) == 0)
fisheries_only <- fisheries_4yr[non_matching_indices,]

# Combine the datasets
all_sites <- rbind(pipefish_sites, fisheries_only)

# Create the map
tmPipeFisheries <- tm_shape(all_sites) +
  tm_dots(col = "category", 
          palette = c("Bare" = "khaki1",
                      "Natural" = "darkolivegreen",
                      "4-years-old (Pipefish Only)" = "darkorchid1",
                      "4-years-old (Both)" = "black",
                      "4-years-old (Fisheries Only)" = "skyblue3"),
          size = 0.1,
          title = "Transplant Locations") +
  tm_layout(legend.outside = TRUE)

tmap_save(tmPipeFisheries, filename = paste0("outputs/plots/tmPipeFisheries", run_date, ".png"))



