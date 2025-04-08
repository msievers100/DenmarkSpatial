## Denmark Scenarios - Plot Scenario Objectives
## Sievers and Brown
## 2025-02-21

library(ggplot2)
library(tidyr)
library(dplyr)
library(viridis)
library(patchwork)

# Load scenario outputs
scenarios_distance <- readRDS("outputs/2024-08-19scenarios-distance.rds")
scenarios_random <- readRDS("outputs/2024-08-19scenarios-random.rds")
#greedy_results <- readRDS("outputs/greedy-results/best_guess_results_ispp_fisheries_2025-03-13.rds")
#greedy_results <- readRDS("outputs/greedy-results/best_guess_results_ispp_pipe_2025-03-13.rds")
#greedy_results <- readRDS("outputs/greedy-results/best_guess_results_ispp_SGfunction_2025-03-13.rds")
#greedy_results <- readRDS("outputs/greedy-results/best_guess_results_ispp_all_2025-03-13.rds")

#greedy_results <- readRDS("outputs/greedy-results/best_guess_results_baltic_2025-03-13.rds")
#greedy_results <- readRDS("outputs/greedy-results/best_guess_results_browns_2025-01-28.rds")
#greedy_results <- readRDS("outputs/greedy-results/best_guess_results_flatfish_2025-01-28.rds")
#greedy_results <- readRDS("outputs/greedy-results/best_guess_results_peri_2025-01-28.rds")
greedy_results <- readRDS("outputs/greedy-results/best_guess_results_whelk_2025-01-28.rds")

scalevalsdf <- read.csv("outputs/scalevals.csv")

# Convert greedy results to scenario format
greedy_progress <- do.call(rbind, lapply(greedy_results$progress, function(x) {
  data.frame(area = x$n_transplants, abundance = x$abundance)
}))

scenarios_distance$area <- as.numeric(scenarios_distance$transplanted_grids)
scenarios_random$area <- as.numeric(scenarios_random$transplanted_grids)


# Get species names from columns ending in _n
species_names <- names(scenarios_distance)[grep("_n$", names(scenarios_distance))]

# Function to calculate scaled abundance
calc_scaled_abundance <- function(df, species_cols, scalevals) {
  scaled_vals <- mapply(function(x, scale) x/scale, 
                       df[species_cols], 
                       scalevals)
                       
  rowSums(scaled_vals)/length(species_cols)
}

# Define species groups
species_indices <- list(
  fisheries = c(2, 4),
  sgfunction = c(1, 2, 5),
  all = seq_along(species_names),
  pipe = 1, 
  baltic = 2,
  browns = 3,
  flatfish = 4,
  peri = 5,
  whelk = 6
)

#arrange order so it matches order of species in scenarios
scalevals <- scalevalsdf$scaleval[match(names(scenarios_distance)[grep("_n$", names(scenarios_distance))], 
  paste0(scalevalsdf$species, "_n"))]


# Calculate objective functions for each species group and scenario type
for (group_name in names(species_indices)) {
  group_name <- "whelk"
  group_species <- species_names[species_indices[[group_name]]]
  group_scalevals <- scalevals[species_indices[[group_name]]]
  
  scenarios_distance[[paste0("scaled_abundance_", group_name)]] <- 
    calc_scaled_abundance(scenarios_distance, group_species, group_scalevals)
  
  scenarios_random[[paste0("scaled_abundance_", group_name)]] <- 
    calc_scaled_abundance(scenarios_random, group_species, group_scalevals)
}

# Create plots for each species group
objective_plots <- lapply(names(species_indices), function(group_name) {
  ggplot() +
    geom_line(data = scenarios_random,
              aes(x = area, y = .data[[paste0("scaled_abundance_", group_name)]], 
                  group = scenario_name,
                  color = "Random"),
              alpha = 0.3) +
    geom_line(data = scenarios_distance,
              aes(x = area, y = .data[[paste0("scaled_abundance_", group_name)]],
                  group = scenario_name,
                  color = scenario_name)) +
    scale_color_viridis_d() +
    labs(x = "Sites Restored",
         y = "Scaled Abundance",
         title = paste("Species Group:", group_name),
         color = "Scenario") +
    theme_classic()
})



# Combine objective function plots and add greedy optimization results
objective_plots <- lapply(names(species_indices), function(group_name) {
  group_name <- "whelk"
  p <- ggplot() +
    geom_line(data = scenarios_random,
              aes(x = area, y = .data[[paste0("scaled_abundance_", group_name)]], 
                  group = scenario_name,
                  color = "Random"),
              alpha = 0.3) +
    geom_line(data = scenarios_distance,
              aes(x = area, y = .data[[paste0("scaled_abundance_", group_name)]],
                  group = scenario_name,
                  color = scenario_name), linewidth = 1)
  
  if(group_name == "whelk") {
    p <- p + geom_line(data = greedy_progress,
                       aes(x = area, y = abundance,
                           color = "Optimized"),
                       linewidth = 1)
  }
  
  p + scale_color_viridis_d() +
    labs(x = "Sites Restored",
         y = "Scaled Abundance",
         #title = paste("Species Group:", group_name),
         color = "Scenario") +
    theme_classic()
})

objective_plots[[9]]



#### Start Sievers new code - Using % of maximum instead ####
# First, find the maximum value across all datasets for this species group
group_name <- names(species_indices)[9]

max_value <- max(
  max(scenarios_random[[paste0("scaled_abundance_", group_name)]], na.rm = TRUE),
  max(scenarios_distance[[paste0("scaled_abundance_", group_name)]], na.rm = TRUE),
  if(group_name == "whelk") max(greedy_progress$abundance, na.rm = TRUE) else -Inf
)

# Modify the plot
p <- ggplot() +
  geom_line(data = scenarios_random,
            aes(x = area, 
                y = .data[[paste0("scaled_abundance_", group_name)]] / max_value * 100, 
                group = scenario_name,
                color = "Random"),
            alpha = 0.3) +
  geom_line(data = scenarios_distance,
            aes(x = area, 
                y = .data[[paste0("scaled_abundance_", group_name)]] / max_value * 100,
                group = scenario_name,
                color = scenario_name),
            linewidth = 1)

if(group_name == "whelk") {
  p <- p + geom_line(data = greedy_progress,
                     aes(x = area, 
                         y = abundance / max_value * 100,
                         color = "Optimized"),
                     linewidth = 1)
}

p2 <- p + scale_color_viridis_d() +
  labs(x = "Sites Restored",
       y = "Scaled abundance (% of max)",
       #title = paste("Species Group:", group_name),
       color = "Scenario") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_y_continuous(limits = c(0, 100)) +
  theme_classic()
p2
saveRDS(p2, "outputs/whelk0plot.rds")
ggsave("outputs/plots/whelk0plot.png", 
       plot = p2, 
       width = 4.8, 
       height = 2.5, 
       dpi = 300)


#data <- readRDS("outputs/sgfunction2650plot.rds")
#data2 <- data + scale_y_continuous(limits = c(33, 100))
#data2
#ggsave("outputs/plots/sgfunction2650plot.png", 
#       plot = data2, 
#       width = 4.8, 
#       height = 2.5, 
#       dpi = 300)








 #### End Sievers new code ####





p1 <- wrap_plots(objective_plots, ncol = 2)

# Plot individual species responses
species_plots <- lapply(species_names, function(sp) {
  sp_name <- gsub("_n$", "", sp)
  
  ggplot() +
    geom_line(data = scenarios_random,
              aes(x = area, y = .data[[sp]], 
                  group = scenario_name,
                  color = "Random"),
              alpha = 0.3) +
    geom_line(data = scenarios_distance,
              aes(x = area, y = .data[[sp]],
                  group = scenario_name,
                  color = scenario_name)) +
    scale_color_viridis_d() +
    labs(x = "Area Restored (m²)",
         y = "Abundance",
         title = sp_name,
         color = "Scenario") +
    theme_classic()
})

# Combine plots
p2 <- wrap_plots(species_plots, ncol = 2)

# Save plots
ggsave("outputs/plots/scenario-objective-functions.png", p1, width = 12, height = 10)
ggsave("outputs/plots/scenario-species-responses.png", p2, width = 12, height = 8)
