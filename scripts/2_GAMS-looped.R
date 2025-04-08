## Denmark GAMs
## Loops over species 
## Sievers and Brown


library(mgcv)
library(visreg)
library(gratia)
library(ggplot2)

theme_set(theme_classic())

#Proportion of deviance explained:
# Define a function to calculate deviance explained by each term
calculate_deviance_explained <- function(gam2, reduced_formula) {
  reduced_model <- update(gam2, formula = reduced_formula)
  (deviance(reduced_model) - deviance(gam2)) / deviance(reduced_model) * 100
}


df <- read.csv("data/datCABMS.csv")
str(df)
df$Treatment <- as.factor(df$Treatment)
df$Treatment <- relevel(df$Treatment, ref = "Bare")

colnames(df)
selected_species <- c("Lesser.pipefish", "Baltic.prawn", "Crangon.prawn", "Flatfish", "Periwinkle", "Whelk")             
nspp_fit <- length(selected_species)
# Loop over species to create a list of GAM fits
gam_list <- lapply(selected_species, function(species) {
  print(species)
  gam_model <- gam(as.formula(paste(species, "~ s(DistBouldEdge, k = 3) + s(DistMussEdge, k = 3) + s(Area10m, k = 3) + Treatment")), data = df, family = nb())
  return(gam_model)
})


names(gam_list) <- selected_species

#Save GAM
saveRDS(gam_list, file = "outputs/gam-selected-species.rds")


#
# Interpretation of GAMs
#

deviance_df <- NULL

for (i in 1:nspp_fit) {
  species <- selected_species[i]
  
gam2 <- gam_list[[i]]

#Interpreting covariate relationship with response. Not using 'draw', to include Treatment
plot_DistBouldEdge <- visreg(gam2, "DistBouldEdge", gg = TRUE) + ggtitle("DistBouldEdge")
plot_DistMussEdge <- visreg(gam2, "DistMussEdge", gg = TRUE) + ggtitle("DistMussEdge")
plot_Area10m <- visreg(gam2, "Area10m", gg = TRUE) + ggtitle("Area10m")
plot_Treatment <- visreg(gam2, "Treatment", gg = TRUE) + ggtitle("Treatment")

combined_plot <- (plot_DistBouldEdge | plot_DistMussEdge) / (plot_Area10m | plot_Treatment)
ggsave(paste0("outputs/plots/",species,"_visreg_plots.png"), plot = combined_plot, width = 8, height = 7, dpi = 150)


# Calculate the deviance explained by each term
df_temp <- data.frame(Species = species,
                      Variable = c("DistBouldEdge", "DistMussEdge", "Area10m", "Treatment"),
                      Deviance_Explained = c(calculate_deviance_explained(gam2, . ~ . - s(DistBouldEdge, k = 3)),
                                             calculate_deviance_explained(gam2, . ~ . - s(DistMussEdge, k = 3)),
                                             calculate_deviance_explained(gam2, . ~ . - s(Area10m, k = 3)),
                                             calculate_deviance_explained(gam2, . ~ . - Treatment))
)

deviance_df <- rbind(deviance_df, df_temp)


}

write.csv(deviance_df, "outputs/deviance_Table.csv")



