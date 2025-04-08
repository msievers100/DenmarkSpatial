#03 July 2024
#Denmark data Cesar ideas 
#Michael Sievers

library(ggplot2)
library(Hmsc)
library(MASS)
library(tidyr)
library(patchwork)
library(sf)
library(tmap)
library(vegan)
library(dplyr)
library(ggcorrplot)
library(tidyr)
library(GGally)
library("viridis")

# distance-based redundancy analysis (dbRDA), which is functionally similar to DISTLM
##Where Ratio = 1 for no area
dat1MV <- read.csv("data/UTMdat1MV.csv")
str(dat1MV)
Species <- dat1MV[,-c(1:9)]

#vegdist Function: Calculates a distance matrix using the specified method (e.g., Bray-Curtis).
Speciesdist_matrix <- vegdist(Species, method = "bray") 
#capscale Function: Performs constrained ordination (dbRDA) on the distance matrix.
#dbRDA_model <- capscale(Speciesdist_matrix ~ Area10m + RatioLog + DistBouldEdge + DistMussEdge + UTM_E + UTM_N + Treatment, data = dat1MV)

#without lat and long (correlated)
dbRDA_model <- capscale(Speciesdist_matrix ~ Treatment + Area10m + RatioLog + DistBouldEdge + DistMussEdge, data = dat1MV)
summary(dbRDA_model)
#Partitioning of squared Bray distance:
#Inertia Proportion
#Total          32.676     1.0000
#Constrained     9.665     0.2958
#Unconstrained  23.012     0.7042

# ANOVA to assess significance
anova(dbRDA_model, permutations = 999)
# Check the significance of each predictor
anova(dbRDA_model, by = "term", permutations = 999)

# Extract residuals and fitted values
residuals_dbRDA <- residuals(dbRDA_model)
fitted_dbRDA <- fitted(dbRDA_model)

# Plot Residuals vs. Fitted Values
plot(fitted_dbRDA, residuals_dbRDA,
     xlab = "Fitted values",
     ylab = "Residuals",
     main = "Residuals vs Fitted Values")

# QQ Plot for Residuals
qqnorm(residuals_dbRDA)
qqline(residuals_dbRDA, col = "red")

# Screeplot for eigenvalues
screeplot(dbRDA_model, main = "Screeplot of Eigenvalues")

# Ordination plot
plot(dbRDA_model, display = "sites", main = "Ordination Plot")



##Removed Treatment for plotting; didn't change much but makes the plot clearer
dbRDA_model <- capscale(Speciesdist_matrix ~ Area10m + RatioLog + DistBouldEdge + DistMussEdge, data = dat1MV)

plot(dbRDA_model)
# Extract the site and biplot scores
site_scores <- scores(dbRDA_model, display = "sites")
biplot_scores <- scores(dbRDA_model, display = "bp")

# Convert to data frames for ggplot2
site_scores_df <- as.data.frame(site_scores)
biplot_scores_df <- as.data.frame(biplot_scores)

site_scores_df$Treatment <- dat1MV$Treatment
site_scores_df$Treatment <- relevel(site_scores_df$Treatment, ref = "Bare")

# Create a ggplot of the site scores and biplot arrows
# Create new columns for shifted text position
biplot_scores_df$CAP1_shifted <- biplot_scores_df$CAP1 * 1.8
biplot_scores_df$CAP2_shifted <- biplot_scores_df$CAP2 * 1.5

# Plot
dbRDA_plot <- ggplot() +
  geom_point(data = site_scores_df, aes(x = CAP1, y = CAP2, color = Treatment), size = 2) +
  geom_segment(data = biplot_scores_df, aes(x = 0, y = 0, xend = CAP1, yend = CAP2), 
               arrow = arrow(length = unit(0.3, "cm")), color = "red") +
  stat_ellipse(data = site_scores_df, aes(x = CAP1, y = CAP2, color = Treatment), 
               level = 0.95, geom = "polygon", alpha = 0.1) + 
  geom_text(data = biplot_scores_df, aes(x = CAP1_shifted, y = CAP2_shifted, label = rownames(biplot_scores_df)), 
            color = "black", size = 5) + 
  annotate("text", x = -Inf, y = Inf, label = "Constrained inertia propn = 0.3", hjust = -0.05, vjust = 2, size = 4) +
  scale_colour_viridis(discrete = TRUE, option = "D") +
  theme_bw() +
  labs(title = "", x = "dbRDA Axis 1 (74%)", y = "dbRDA Axis 2 (11%)")

print(dbRDA_plot)

