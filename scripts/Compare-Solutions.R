# Load required libraries
library(ggplot2)
library(tidyverse)

# Load results
greedy_results <- readRDS("outputs/greedy_objectives.rds")
heuristic_results <- readRDS("outputs/heuristic_objectives.rds")

# Combine results
all_results <- bind_rows(
  greedy_results %>% mutate(method = "Greedy"),
  heuristic_results %>% mutate(method = "Heuristic")
)

# Create comparison plot
ggplot(all_results, aes(x = area, y = objective, color = method)) +
  geom_line() +
  geom_point() +
  theme_minimal() +
  labs(x = "Area of Transplants",
       y = "Objective Function Value",
       title = "Comparison of Optimization Methods",
       color = "Method") +
  theme(legend.position = "bottom")

# Save plot
ggsave("outputs/optimization_comparison.png", width = 8, height = 6)

# Calculate and print summary statistics
summary_stats <- all_results %>%
  group_by(method) %>%
  summarise(
    mean_objective = mean(objective),
    max_objective = max(objective),
    min_objective = min(objective)
  )

print(summary_stats)
