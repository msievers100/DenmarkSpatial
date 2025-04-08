#31 June 2024
#Denmark Spatial 
#Michael Sievers

library(mgcv)
library(MuMIn)
library(visreg)
library(ggplot2)
library("viridis") 
library(cowplot)
library(vegan)
library(ggrepel)
library(GGally)
library(gstat)
library(sf)
library(gratia)
library(patchwork)
library(spdep)
library(gstat)
library(sp)
library(emmeans)
library(ape)
library(dplyr)

dat <- read.csv("data/UTMdat.csv")
dat$Treatment <- as.factor(dat$Treatment)
dat$Treatment <- relevel(dat$Treatment, ref = "Bare")

##Row 58 is large outlier for the Per-Area ratio. 
dat2 <- dat[-58, ]
#Log transformed 
dat2$RatioLog <- log(dat2$PerAreaRatio10m + 1)

# Function to format the concurvity matrix
format_concurvity <- function(concur_mat, digits = 4) {
  # Convert scientific notation to decimal and round
  formatted_mat <- lapply(concur_mat, function(x) {
    # Convert matrix to data frame
    df <- as.data.frame(x)
    # Format all numeric columns
    df[] <- lapply(df, function(y) {
      format(round(y, digits), scientific = FALSE, nsmall = digits)
    })
    # Convert back to matrix
    as.matrix(df)
  })
  return(formatted_mat)
}

cor(dat2[, c("Area10m", "RatioLog", "DistBouldEdge", "DistMussEdge")])

#### Richness ####

Richm2a <- gam(Richness ~ Treatment + s(Area10m, k = 3) + s(RatioLog, k = 3) + 
                 s(DistBouldEdge, k = 3) + s(DistMussEdge, k =3) + s(UTM_E, UTM_N, bs="gp", k = 4), data = dat2, na.action = na.fail)
summary(Richm2a)
RichApp <- appraise(Richm2a)

draw(Richm2a, residuals = TRUE) & theme_cowplot(12)
options(scipen = 10)
concurvity(Richm2a, full = TRUE)
concurvity(Richm2a, full = FALSE)
gam.check(Richm2a)

Richm2a <- gam(Richness ~ Treatment + s(Area10m, k = 3) + s(RatioLog, k = 3) + 
                 s(DistBouldEdge, k = 3) + s(DistMussEdge, k =3) + s(UTM_E, UTM_N, bs="gp", k = 4), data = dat2, na.action = na.fail)
# Concurvity full model
conc_result <- concurvity(Richm2a, full = FALSE)
formatted_result <- format_concurvity(conc_result)
print(formatted_result$estimate, quote = FALSE)

anova(Richm2a)

#Proportion of deviance explained:
calculate_deviance_explained <- function(Richm2a, reduced_formula) {
  reduced_model <- update(Richm2a, formula = reduced_formula)
  (deviance(reduced_model) - deviance(Richm2a)) / deviance(reduced_model) * 100
}

calculate_deviance_explained(Richm2a, . ~ . - s(Area10m, k = 3))
calculate_deviance_explained(Richm2a, . ~ . - s(RatioLog, k = 3))
calculate_deviance_explained(Richm2a, . ~ . - s(DistBouldEdge, k = 3))
calculate_deviance_explained(Richm2a, . ~ . - s(DistMussEdge, k = 3))
calculate_deviance_explained(Richm2a, . ~ . - s(UTM_E,UTM_N, bs="gp", k = 4))
calculate_deviance_explained(Richm2a, . ~ . - Treatment)

#Seperate out plots to add more info to each
p1 <- draw(Richm2a, select = "s(Area10m)", residuals = TRUE) + ggtitle(NULL) + theme_cowplot(12) +
  annotate("text", x = -Inf, y = Inf, label = "p = 0.33; Dev = 2.3%", hjust = -0.05, vjust = 2, size = 4, color = "black") + labs(x = expression("Area of seagrass within 10m radius (m"^2*")"))
p2 <- draw(Richm2a, select = "s(RatioLog)", residuals = TRUE) + ggtitle(NULL) + theme_cowplot(12)+
  annotate("text", x = -Inf, y = Inf, label = "p = 0.16; Dev = 1.3%", hjust = -0.05, vjust = 2, size = 4, color = "black") + labs(x = "Seagrass perimeter:area (10m radius)")
p3 <- draw(Richm2a, select = "s(DistBouldEdge)", residuals = TRUE) + ggtitle(NULL) + theme_cowplot(12) +
  annotate("text", x = -Inf, y = Inf, label = "p = 0.42; Dev = 1.4%", hjust = -0.05, vjust = 2, size = 4, color = "black") + labs(x = "Distance to boulder reef (m)")
p4 <- draw(Richm2a, select = "s(DistMussEdge)", residuals = TRUE) + ggtitle(NULL) + theme_cowplot(12)+
  annotate("text", x = -Inf, y = Inf, label = "p = 0.01; Dev = 5.3%", hjust = -0.05, vjust = 2, size = 4, color = "black", fontface = "bold") + labs(x = "Distance to mussel reef (m)")
p5 <- draw(Richm2a, select = "s(UTM_E,UTM_N)", residuals = TRUE) + ggtitle(NULL) + theme_cowplot(12) +
  annotate("text", x = 542900, y = 6171620, label = "p = 0.11; \n Dev = 4.1%", size = 4, color = "black") + theme(axis.text.x = element_blank(),
                                                                                                                  axis.text.y = element_blank()) 

##Plotting Treatment effect

#Suppressing intercept for Treatment
Richm2bS <- gam(Richness ~ 0+Treatment + s(Area10m, k = 3) + s(RatioLog, k = 3) + 
                  s(DistBouldEdge, k = 3) + s(DistMussEdge, k = 3) + s(UTM_E, UTM_N, bs="gp", k = 4), data = dat2, na.action = na.fail)

# Extract the partial effects and standard errors for the Treatment factor
preds <- predict(Richm2bS, type = "terms", se.fit = TRUE, terms = "Treatment")
dat2$TreatmentEffect <- preds$fit[, "Treatment"]
dat2$SE <- preds$se.fit[, "Treatment"]

# Calculate the confidence intervals
dat2$CI_upper <- dat2$TreatmentEffect + 1.96 * dat2$SE
dat2$CI_lower <- dat2$TreatmentEffect - 1.96 * dat2$SE

# Plot the Treatment effect with confidence intervals
p6 <- ggplot(dat2, aes(x = Treatment, y = TreatmentEffect)) +
  geom_point() +
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper), width = 0.2) +
  theme_cowplot(12) +
  labs(x = "Treatment", y = "Partial effect") +
  annotate("text", x = -Inf, y = Inf, label = "p < 0.001; Dev = 16.9%", hjust = -0.05, vjust = 1.5, size = 4, color = "black", fontface = "bold") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

RichnessAllPlots <- plot_grid(p1, p2, p3, p4, p5, p6, labels = c('A', 'B', 'C', 'D','E', 'F'), label_size = 12, ncol = 2)
RichnessAllPlots

vis_smooth1 <- visreg(Richm2a, "DistMussEdge", plot = FALSE)
vis_cat <- visreg(Richm2a, "Treatment", by = "Treatment", plot = FALSE)

RichnessPlot1 <- ggplot(vis_smooth1$fit, aes(x = DistMussEdge, y = visregFit)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = visregLwr, ymax = visregUpr), alpha = 0.2) +
  labs(x = "Distance to mussel reef (m)", y = "Richness") +
  theme_minimal(base_size = 15)+
  annotate("text", x = -Inf, y = Inf, label = "p = 0.01", hjust = -0.05, vjust = 2, size = 4, color = "black", fontface = "bold") +
  scale_y_continuous(limits = c(4, 10)) +
  theme_cowplot(12)

RichnessPlot2 <- ggplot(vis_cat$fit, aes(x = as.factor(Treatment), y = visregFit)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = visregLwr, ymax = visregUpr), width = 0.2) +
  scale_fill_viridis(discrete = TRUE, option = "D") +
  labs(x = "", y = "Richness") +
  theme_cowplot(12) +
  annotate("text", x = -Inf, y = Inf, label = "p < 0.001", hjust = -0.05, vjust = 1.5, size = 4, color = "black", fontface = "bold") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.title=element_blank(), 
        legend.position = "none")

RichnessPlots <- plot_grid(RichnessPlot1, RichnessPlot2, labels = c('A', 'B'), label_size = 12)
RichnessPlots



dredged_models <- dredge(Richm2a)
print(dredged_models)
best_model <- get.models(dredged_models, 1)[[1]]
summary(best_model)
RichAppDredge <- appraise(best_model)

anova(best_model)

#Proportion of deviance explained:
calculate_deviance_explained <- function(best_model, reduced_formula) {
  reduced_model <- update(best_model, formula = reduced_formula)
  (deviance(reduced_model) - deviance(best_model)) / deviance(reduced_model) * 100
}

calculate_deviance_explained(best_model, . ~ . - s(DistMussEdge, k = 3))
calculate_deviance_explained(best_model, . ~ . - Treatment)

vis_smooth1 <- visreg(best_model, "DistMussEdge", scale="response",plot = FALSE)
#vis_smooth2 <- visreg(Richm2a, "Latitude", plot = FALSE)
#vis_smooth3 <- visreg(Richm2a, "Longitude", plot = FALSE)
vis_cat <- visreg(best_model, "Treatment", by = "Treatment", scale="response",plot = FALSE)

RichnessPlot1 <- ggplot(vis_smooth1$fit, aes(x = DistMussEdge, y = visregFit)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = visregLwr, ymax = visregUpr), alpha = 0.2) +
  labs(x = "Distance to mussel reef (m)", y = "Richness") +
  theme_minimal(base_size = 15)+
  annotate("text", x = -Inf, y = Inf, label = "p < 0.001", hjust = -0.05, vjust = 2, size = 4, color = "black", fontface = "bold") +
  #annotate("text", x = -Inf, y = Inf, label = "p < 0.001; Dev = 9.0%", hjust = -0.05, vjust = 2, size = 4, color = "black", fontface = "bold") +
  #scale_y_continuous(breaks = seq(0, 10, by = 2)) +
  #scale_y_continuous(limits = c(4, 10)) +
  theme_cowplot(12)

RichnessPlot2 <- ggplot(vis_cat$fit, aes(x = as.factor(Treatment), y = visregFit)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = visregLwr, ymax = visregUpr), width = 0.2) +
  scale_fill_viridis(discrete = TRUE, option = "D") +
  labs(x = "", y = "Richness") +
  theme_cowplot(12) +
  annotate("text", x = -Inf, y = Inf, label = "p < 0.001", hjust = -0.05, vjust = 1.5, size = 4, color = "black", fontface = "bold") +
  #annotate("text", x = -Inf, y = Inf, label = "p < 0.001; Dev = 22.5%", hjust = -0.05, vjust = 1.5, size = 4, color = "black", fontface = "bold") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.title=element_blank(), 
        legend.position = "none") +
  scale_y_continuous(limits = c(0, 7.5)) 

RichnessPlotsDredge <- plot_grid(RichnessPlot1, RichnessPlot2, label_size = 12)
RichnessPlotsDredge
ggsave(RichnessPlotsDredge, filename = "2024-06-31-RichnessPlotsDredge.png", width = 6, height = 2.5)

# Create a data frame with coordinates and residuals
resid_df <- data.frame(
  x = dat2$UTM_E,
  y = dat2$UTM_N,
  resid = residuals(best_model))

# Convert to SpatialPointsDataFrame
coordinates(resid_df) <- ~x+y

# Create variogram
var_model <- variogram(resid~1, data=resid_df)
plot(var_model)

# Create distance matrix from UTM coordinates
dists <- dist(cbind(dat2$UTM_E, dat2$UTM_N))
dists_matrix <- as.matrix(dists)

# Perform Moran's I test
Moran.I(residuals(best_model), dists_matrix)



#### Shannon Diversity ####

Shannonm2a <- gam(Shannon ~ Treatment + s(Area10m, k = 3) + s(RatioLog, k = 3) + 
                    s(DistBouldEdge, k = 3) + s(DistMussEdge, k =3) + s(UTM_E, UTM_N, bs="gp", k = 4), data = dat2, na.action = na.fail)
summary(Shannonm2a)
ShannonApp <- appraise(Shannonm2a)

draw(Shannonm2a, residuals = TRUE) & theme_cowplot(12)
options(scipen = 10)
concurvity(Shannonm2a, full = TRUE)
concurvity(Shannonm2a, full = FALSE)
gam.check(Shannonm2a)

# Concurvity full model
conc_result <- concurvity(Shannonm2a, full = FALSE)
formatted_result <- format_concurvity(conc_result)
print(formatted_result$estimate, quote = FALSE)

anova(Shannonm2a)

#Proportion of deviance explained:
calculate_deviance_explained <- function(Shannonm2a, reduced_formula) {
  reduced_model <- update(Shannonm2a, formula = reduced_formula)
  (deviance(reduced_model) - deviance(Shannonm2a)) / deviance(reduced_model) * 100
}

calculate_deviance_explained(Shannonm2a, . ~ . - s(Area10m, k = 3))
calculate_deviance_explained(Shannonm2a, . ~ . - s(RatioLog, k = 3))
calculate_deviance_explained(Shannonm2a, . ~ . - s(DistBouldEdge, k = 3))
calculate_deviance_explained(Shannonm2a, . ~ . - s(DistMussEdge, k = 3))
calculate_deviance_explained(Shannonm2a, . ~ . - s(UTM_E,UTM_N, bs="gp", k = 4))
calculate_deviance_explained(Shannonm2a, . ~ . - Treatment)

#Seperate out plots to add more info to each
p1 <- draw(Shannonm2a, select = "s(Area10m)", residuals = TRUE) + ggtitle(NULL) + theme_cowplot(12) +
  annotate("text", x = -Inf, y = Inf, label = "p = 0.50; Dev = 1.8%", hjust = -0.05, vjust = 2, size = 4, color = "black") + labs(x = expression("Area of seagrass within 10m radius (m"^2*")"))
p2 <- draw(Shannonm2a, select = "s(RatioLog)", residuals = TRUE) + ggtitle(NULL) + theme_cowplot(12)+
  annotate("text", x = -Inf, y = Inf, label = "p = 0.13; Dev = 1.6%", hjust = -0.05, vjust = 2, size = 4, color = "black") + labs(x = "Seagrass perimeter:area (10m radius)")
p3 <- draw(Shannonm2a, select = "s(DistBouldEdge)", residuals = TRUE) + ggtitle(NULL) + theme_cowplot(12) +
  annotate("text", x = -Inf, y = Inf, label = "p = 0.34; Dev = 1.5%", hjust = -0.05, vjust = 2, size = 4, color = "black") + labs(x = "Distance to boulder reef (m)")
p4 <- draw(Shannonm2a, select = "s(DistMussEdge)", residuals = TRUE) + ggtitle(NULL) + theme_cowplot(12)+
  annotate("text", x = -Inf, y = Inf, label = "p = 0.007; Dev = 5.1%", hjust = -0.05, vjust = 2, size = 4, color = "black", fontface = "bold") + labs(x = "Distance to mussel reef (m)")
p5 <- draw(Shannonm2a, select = "s(UTM_E,UTM_N)", residuals = TRUE) + ggtitle(NULL) + theme_cowplot(12) +
  annotate("text", x = 542900, y = 6171620, label = "p = 0.018; \n Dev = 6.0%", size = 4, color = "black", fontface = "bold") + 
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank()) 

##Plotting Treatment effect

#Suppressing intercept for Treatment
Shannonm2bS <- gam(Shannon ~ 0+Treatment + s(Area10m, k = 3) + s(RatioLog, k = 3) + 
                     s(DistBouldEdge, k = 3) + s(DistMussEdge, k = 3) + s(UTM_E, UTM_N, bs="gp", k = 4), data = dat2, na.action = na.fail)

# Extract the partial effects and standard errors for the Treatment factor
preds <- predict(Shannonm2bS, type = "terms", se.fit = TRUE, terms = "Treatment")
dat2$TreatmentEffect <- preds$fit[, "Treatment"]
dat2$SE <- preds$se.fit[, "Treatment"]

# Calculate the confidence intervals
dat2$CI_upper <- dat2$TreatmentEffect + 1.96 * dat2$SE
dat2$CI_lower <- dat2$TreatmentEffect - 1.96 * dat2$SE

# Plot the Treatment effect with confidence intervals
p6 <- ggplot(dat2, aes(x = Treatment, y = TreatmentEffect)) +
  geom_point() +
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper), width = 0.2) +
  theme_cowplot(12) +
  labs(x = "Treatment", y = "Partial effect") +
  annotate("text", x = -Inf, y = Inf, label = "p = 0.23; Dev = 5.5%", hjust = -0.05, vjust = 1.5, size = 4, color = "black") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ShannonAllPlots <- plot_grid(p1, p2, p3, p4, p5, p6, labels = c('A', 'B', 'C', 'D','E', 'F'), label_size = 12, ncol = 2)
ShannonAllPlots

dredged_models <- dredge(Shannonm2a)
print(dredged_models)
best_model <- get.models(dredged_models, 1)[[1]]
summary(best_model)
ShannonAppDredge <- appraise(best_model)

anova(best_model)

#Proportion of deviance explained:
calculate_deviance_explained <- function(best_model, reduced_formula) {
  reduced_model <- update(best_model, formula = reduced_formula)
  (deviance(reduced_model) - deviance(best_model)) / deviance(reduced_model) * 100
}

calculate_deviance_explained(best_model, . ~ . - s(RatioLog, k = 3))
calculate_deviance_explained(best_model, . ~ . - s(DistMussEdge, k = 3))

vis_smooth1 <- visreg(best_model, "RatioLog", scale="response",plot = FALSE)
vis_smooth3 <- visreg(best_model, "DistMussEdge", scale="response",plot = FALSE)

ShannonPlot1 <- ggplot(vis_smooth1$fit, aes(x = RatioLog, y = visregFit)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = visregLwr, ymax = visregUpr), alpha = 0.2) +
  labs(x = expression(atop("Seagrass perimeter:area", "(10m radius)")), y = "Shannon") +
  labs(x = "Seagrass perimeter:area (10m radius)", y = "Shannon") +
    theme_minimal(base_size = 15)+
  annotate("text", x = -Inf, y = Inf, label = "p = 0.05", hjust = -0.05, vjust = 2, size = 4, color = "black") +
  theme_cowplot(12)

ShannonPlot3 <- ggplot(vis_smooth3$fit, aes(x = DistMussEdge, y = visregFit)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = visregLwr, ymax = visregUpr), alpha = 0.2) +
  labs(x = "Distance to mussel reef (m)", y = "Shannon") +
  theme_minimal(base_size = 15)+
  annotate("text", x = -Inf, y = Inf, label = "p = 0.04", hjust = -0.05, vjust = 2, size = 4, color = "black", fontface = "bold") +
  theme_cowplot(12)

ShannonPlotsDredge <- plot_grid(ShannonPlot1, ShannonPlot3, label_size = 12)
ShannonPlotsDredge
ggsave(ShannonPlotsDredge, filename = "2024-06-31-ShannonPlotsDredge.png", width = 6, height = 2.5)

resid_df <- data.frame(
  x = dat2$UTM_E,
  y = dat2$UTM_N,
  resid = residuals(best_model))
coordinates(resid_df) <- ~x+y
var_model <- variogram(resid~1, data=resid_df)
plot(var_model)
dists <- dist(cbind(dat2$UTM_E, dat2$UTM_N))
dists_matrix <- as.matrix(dists)
Moran.I(residuals(best_model), dists_matrix)


#### Evenness ####
## One site missing evenness data because of 1 species
dat2$Evenness
#Row 54
datEven <- dat2[-54, ]
datEven$Evenness

Evennessm2a <- gam(Evenness ~ Treatment + s(Area10m, k = 3) + s(RatioLog, k = 3) + 
                     s(DistBouldEdge, k = 3) + s(DistMussEdge, k =3) + s(UTM_E, UTM_N, bs="gp", k = 4), data = datEven, na.action = na.fail)
summary(Evennessm2a)
EvennessApp <- appraise(Evennessm2a)

draw(Evennessm2a, residuals = TRUE) & theme_cowplot(12)
options(scipen = 10)
concurvity(Evennessm2a, full = TRUE)
concurvity(Evennessm2a, full = FALSE)
gam.check(Evennessm2a)

# Concurvity full model
conc_result <- concurvity(Evennessm2a, full = FALSE)
formatted_result <- format_concurvity(conc_result)
print(formatted_result$estimate, quote = FALSE)

anova(Evennessm2a)

#Proportion of deviance explained:
calculate_deviance_explained <- function(Evennessm2a, reduced_formula) {
  reduced_model <- update(Evennessm2a, formula = reduced_formula)
  (deviance(reduced_model) - deviance(Evennessm2a)) / deviance(reduced_model) * 100
}

calculate_deviance_explained(Evennessm2a, . ~ . - s(Area10m, k = 3))
calculate_deviance_explained(Evennessm2a, . ~ . - s(RatioLog, k = 3))
calculate_deviance_explained(Evennessm2a, . ~ . - s(DistBouldEdge, k = 3))
calculate_deviance_explained(Evennessm2a, . ~ . - s(DistMussEdge, k = 3))
calculate_deviance_explained(Evennessm2a, . ~ . - s(UTM_E,UTM_N, bs="gp", k = 4))
calculate_deviance_explained(Evennessm2a, . ~ . - Treatment)

#Seperate out plots to add more info to each
p1 <- draw(Evennessm2a, select = "s(Area10m)", residuals = TRUE) + ggtitle(NULL) + theme_cowplot(12) +
  annotate("text", x = -Inf, y = Inf, label = "p = 0.27; Dev = 1.0%", hjust = -0.05, vjust = 2, size = 4, color = "black") + labs(x = expression("Area of seagrass within 10m radius (m"^2*")"))
p2 <- draw(Evennessm2a, select = "s(RatioLog)", residuals = TRUE) + ggtitle(NULL) + theme_cowplot(12)+
  annotate("text", x = -Inf, y = Inf, label = "p = 0.85; Dev = 0.02%", hjust = -0.05, vjust = 2, size = 4, color = "black") + labs(x = "Seagrass perimeter:area (10m radius)")
p3 <- draw(Evennessm2a, select = "s(DistBouldEdge)", residuals = TRUE) + ggtitle(NULL) + theme_cowplot(12) +
  annotate("text", x = -Inf, y = Inf, label = "p = 0.60; Dev = 0.2%", hjust = -0.05, vjust = 2, size = 4, color = "black") + labs(x = "Distance to boulder reef (m)")
p4 <- draw(Evennessm2a, select = "s(DistMussEdge)", residuals = TRUE) + ggtitle(NULL) + theme_cowplot(12)+
  annotate("text", x = -Inf, y = Inf, label = "p = 0.19; Dev = 0.4%", hjust = -0.05, vjust = 2, size = 4, color = "black") + labs(x = "Distance to mussel reef (m)")
p5 <- draw(Evennessm2a, select = "s(UTM_E,UTM_N)", residuals = TRUE) + ggtitle(NULL) + theme_cowplot(12) +
  annotate("text", x = 542900, y = 6171620, label = "p = 0.91; \n Dev = 0.2%", size = 4, color = "black") + 
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank()) 

##Plotting Treatment effect

#Suppressing intercept for Treatment
Evennessm2bS <- gam(Evenness ~ 0+ Treatment + s(Area10m, k = 3) + s(RatioLog, k = 3) + 
                      s(DistBouldEdge, k = 3) + s(DistMussEdge, k =3) + s(UTM_E, UTM_N, bs="gp", k = 4), data = datEven, action = na.fail)

# Extract the partial effects and standard errors for the Treatment factor
preds <- predict(Evennessm2bS, type = "terms", se.fit = TRUE, terms = "Treatment")
datEven$TreatmentEffect <- preds$fit[, "Treatment"]
datEven$SE <- preds$se.fit[, "Treatment"]

# Calculate the confidence intervals
datEven$CI_upper <- datEven$TreatmentEffect + 1.96 * datEven$SE
datEven$CI_lower <- datEven$TreatmentEffect - 1.96 * datEven$SE

# Plot the Treatment effect with confidence intervals
p6 <- ggplot(datEven, aes(x = Treatment, y = TreatmentEffect)) +
  geom_point() +
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper), width = 0.2) +
  theme_cowplot(12) +
  labs(x = "Treatment", y = "Partial effect") +
  annotate("text", x = -Inf, y = Inf, label = "p = 0.003; Dev = 10.8%", hjust = -0.05, vjust = 1.5, size = 4, color = "black", fontface = "bold") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

EvennessAllPlots <- plot_grid(p1, p2, p3, p4, p5, p6, labels = c('A', 'B', 'C', 'D','E', 'F'), label_size = 12, ncol = 2)
EvennessAllPlots

vis_cat <- visreg(Evennessm2a, "Treatment", by = "Treatment", plot = FALSE)

EvennessPlot2 <- ggplot(vis_cat$fit, aes(x = as.factor(Treatment), y = visregFit)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = visregLwr, ymax = visregUpr), width = 0.2) +
  scale_fill_viridis(discrete = TRUE, option = "D") +
  labs(x = "", y = "Evenness") +
  theme_cowplot(12) +
  annotate("text", x = -Inf, y = Inf, label = "p = 0.003; Dev = 10.8%", hjust = -0.05, vjust = 1.5, size = 4, color = "black", fontface = "bold") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous(limits = c(0, 0.55)) +
  theme(legend.title=element_blank(), 
        legend.position = "none")

EvennessPlot2


dredged_models <- dredge(Evennessm2a)
print(dredged_models)
best_model <- get.models(dredged_models, 1)[[1]]
summary(best_model)

EvennessAppDredge <- appraise(best_model)

anova(best_model)

#Proportion of deviance explained:
calculate_deviance_explained <- function(best_model, reduced_formula) {
  reduced_model <- update(best_model, formula = reduced_formula)
  (deviance(reduced_model) - deviance(best_model)) / deviance(reduced_model) * 100
}

calculate_deviance_explained(best_model, . ~ . - s(DistMussEdge, k = 3))
calculate_deviance_explained(best_model, . ~ . - Treatment)


vis_smooth3 <- visreg(best_model, "DistMussEdge",scale="response",plot = FALSE)

EvennessPlot1 <- ggplot(vis_smooth3$fit, aes(x = DistMussEdge, y = visregFit)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = visregLwr, ymax = visregUpr), alpha = 0.2) +
  labs(x = expression("Distance to mussel reef (m)"), y = "Evenness") +
  theme_minimal(base_size = 15)+
  annotate("text", x = -Inf, y = Inf, label = "p = 0.10", hjust = -0.05, vjust = 2, size = 4, color = "black") +
  theme_cowplot(12)

vis_cat <- visreg(best_model, "Treatment", by = "Treatment", scale="response",plot = FALSE)

EvennessPlot2 <- ggplot(vis_cat$fit, aes(x = as.factor(Treatment), y = visregFit)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = visregLwr, ymax = visregUpr), width = 0.2) +
  scale_fill_viridis(discrete = TRUE, option = "D") +
  labs(x = "", y = "Evenness") +
  theme_cowplot(12) +
  annotate("text", x = -Inf, y = Inf, label = "p = 0.002", hjust = -2.5, vjust = 1.5, size = 4, color = "black", fontface = "bold") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous(limits = c(0, 0.55)) +
  theme(legend.title=element_blank(), 
        legend.position = "none")

EvennessPlotsDredge <- plot_grid(EvennessPlot1, EvennessPlot2, label_size = 12)
EvennessPlotsDredge
ggsave(EvennessPlotsDredge, filename = "2024-06-31-EvennessPlotsDredge.png", width = 6, height = 2.5)


resid_df <- data.frame(
  x = datEven$UTM_E,
  y = datEven$UTM_N,
  resid = residuals(best_model))
coordinates(resid_df) <- ~x+y
var_model <- variogram(resid~1, data=resid_df)
plot(var_model)
dists <- dist(cbind(datEven$UTM_E, datEven$UTM_N))
dists_matrix <- as.matrix(dists)
Moran.I(residuals(best_model), dists_matrix)




#### Lesser.pipefish ####
### Only best model ###
## negative binomial ##

Lesser.pipefishm2a <- gam(Lesser.pipefish ~ Treatment + s(Area10m, k = 3) + s(RatioLog, k = 3) + 
                            s(DistBouldEdge, k = 3) + s(DistMussEdge, k =3) + s(UTM_E, UTM_N, bs="gp", k = 4), 
                          data = dat2, na.action = na.fail, family = nb())

# Concurvity full model
conc_result <- concurvity(Lesser.pipefishm2a, full = FALSE)
formatted_result <- format_concurvity(conc_result)
print(formatted_result$estimate, quote = FALSE)


dredged_models <- dredge(Lesser.pipefishm2a)
print(dredged_models)
best_model <- get.models(dredged_models, 1)[[1]]
summary(best_model)
anova(best_model)
appraise(best_model)

Lesser.pipefishApp <- appraise(best_model)

## Best model
## Lesser.pipefish ~ Treatment + s(Area10m, k = 3) + s(DistBouldEdge, k = 3)
# Get estimated marginal means
emm_result <- emmeans(best_model, specs = "Treatment")
pairs_result <- pairs(emm_result, adjust = "tukey")
pairs_df <- as.data.frame(pairs_result)
cld_result <- multcomp::cld(emm_result)
cld_result

draw(best_model, residuals = TRUE) & theme_cowplot(12)
options(scipen = 10)
concurvity(best_model, full = TRUE)
concurvity(best_model, full = FALSE)
gam.check(best_model)

anova(best_model)

#Proportion of deviance explained:
### This has to be done using the full model

calculate_deviance_explained <- function(best_model, reduced_formula) {
  reduced_model <- update(best_model, formula = reduced_formula)
  (deviance(reduced_model) - deviance(best_model)) / deviance(reduced_model) * 100
}

calculate_deviance_explained(best_model, . ~ . - s(Area10m, k = 3))
#calculate_deviance_explained(best_model, . ~ . - s(RatioLog, k = 3))
calculate_deviance_explained(best_model, . ~ . - s(DistBouldEdge, k = 3))
#calculate_deviance_explained(best_model, . ~ . - s(DistMussEdge, k = 3))
#calculate_deviance_explained(best_model, . ~ . - s(UTM_E,UTM_N, bs="gp", k = 4))
calculate_deviance_explained(best_model, . ~ . - Treatment)

#Seperate out plots to add more info to each
p1 <- draw(best_model, select = "s(Area10m)", residuals = TRUE) + ggtitle(NULL) + theme_cowplot(12) +
  annotate("text", x = -Inf, y = Inf, label = "p = 0.02; Dev = 4.1%", hjust = -0.05, vjust = 2, size = 4, color = "black") + labs(x = expression("Area of seagrass within 10m radius (m"^2*")"))
#p2 <- draw(best_model, select = "s(RatioLog)", residuals = TRUE) + ggtitle(NULL) + theme_cowplot(12)+
#  annotate("text", x = -Inf, y = Inf, label = "p = 0.74; Dev = 0.09%", hjust = -0.05, vjust = 2, size = 4, color = "black") + labs(x = "Seagrass perimeter:area (10m radius)")
p3 <- draw(best_model, select = "s(DistBouldEdge)", residuals = TRUE) + ggtitle(NULL) + theme_cowplot(12) +
  annotate("text", x = -Inf, y = Inf, label = "p = 0.006; Dev = 3.9%", hjust = -0.05, vjust = 2, size = 4, color = "black", fontface = "bold") + labs(x = "Distance to boulder reef (m)")
#p4 <- draw(best_model, select = "s(DistMussEdge)", residuals = TRUE) + ggtitle(NULL) + theme_cowplot(12)+
#  annotate("text", x = -Inf, y = Inf, label = "p = 0.53; Dev = 1.4%", hjust = -0.05, vjust = 2, size = 4, color = "black") + labs(x = "Distance to mussel reef (m)")
#p5 <- draw(best_model, select = "s(UTM_E,UTM_N)", residuals = TRUE) + ggtitle(NULL) + theme_cowplot(12) +
#  annotate("text", x = 542900, y = 6171620, label = "p = 0.40; \n Dev = 2.2%", size = 4, color = "black") + 
#  theme(axis.text.x = element_blank(),
#        axis.text.y = element_blank()) 

##Plotting Treatment effect

#Suppressing intercept for Treatment
Lesser.pipefishm2bS <- gam(Lesser.pipefish ~ 0+Treatment + s(Area10m, k = 3) + 
                             s(DistBouldEdge, k = 3), data = dat2, na.action = na.fail, family = nb())

# Extract the partial effects and standard errors for the Treatment factor
preds <- predict(Lesser.pipefishm2bS, type = "terms", se.fit = TRUE, terms = "Treatment")
dat2$TreatmentEffect <- preds$fit[, "Treatment"]
dat2$SE <- preds$se.fit[, "Treatment"]

# Calculate the confidence intervals
dat2$CI_upper <- dat2$TreatmentEffect + 1.96 * dat2$SE
dat2$CI_lower <- dat2$TreatmentEffect - 1.96 * dat2$SE

# Plot the Treatment effect with confidence intervals
p6 <- ggplot(dat2, aes(x = Treatment, y = TreatmentEffect)) +
  geom_point() +
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper), width = 0.2) +
  theme_cowplot(12) +
  labs(x = "Treatment", y = "Partial effect") +
  annotate("text", x = -Inf, y = Inf, label = "p < 0.001; Dev = 12.2%", hjust = -0.05, vjust = 1.5, size = 4, color = "black", fontface = "bold") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

Lesser.pipefishAllPlots <- plot_grid(p1, p3, p6, labels = c('A', 'B', 'C', 'D','E', 'F'), label_size = 12, ncol = 2)
Lesser.pipefishAllPlots


vis_smooth1 <- visreg(best_model, "Area10m", scale="response",plot = FALSE)
vis_smooth2 <- visreg(best_model, "DistBouldEdge", scale="response",plot = FALSE)

Lesser.pipefishPlot1 <- ggplot(vis_smooth1$fit, aes(x = Area10m, y = visregFit)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = visregLwr, ymax = visregUpr), alpha = 0.2) +
  labs(x = expression(atop("Area of seagrass", "within 10m radius (m"^2*")")), 
       y = "Abundance - Lesser pipefish") +
  #theme_minimal(base_size = 15)+
  annotate("text", x = -Inf, y = Inf, label = "p = 0.02", hjust = -0.05, vjust = 1.5, size = 4, color = "black", fontface = "bold") +
  #annotate("text", x = -Inf, y = Inf, label = "p = 0.02; Dev = 4.1%", hjust = -0.05, vjust = 2, size = 4, color = "black", fontface = "bold") +
  #scale_y_continuous(breaks = seq(0, 10, by = 2)) +
  #scale_y_continuous(limits = c(4, 10)) +
  theme_cowplot(12)

Lesser.pipefishPlot2 <- ggplot(vis_smooth2$fit, aes(x = DistBouldEdge, y = visregFit)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = visregLwr, ymax = visregUpr), alpha = 0.2) +
  labs(x = expression(atop("Distance to", "boulder reef (m)")), y = "Abundance - Lesser pipefish") +
  #theme_minimal(base_size = 15)+
  annotate("text", x = -Inf, y = Inf, label = "p = 0.006", hjust = -0.05, vjust = 1.5, size = 4, color = "black", fontface = "bold") +
  #annotate("text", x = -Inf, y = Inf, label = "p = 0.006; Dev = 3.9%", hjust = -0.05, vjust = 2, size = 4, color = "black", fontface = "bold") +
  #scale_y_continuous(breaks = seq(0, 10, by = 2)) +
  #scale_y_continuous(limits = c(4, 10)) +
  theme_cowplot(12)

vis_cat <- visreg(best_model, "Treatment", by = "Treatment", scale="response",plot = FALSE)

Lesser.pipefishPlot3 <- ggplot(vis_cat$fit, aes(x = as.factor(Treatment), y = visregFit)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = visregLwr, ymax = visregUpr), width = 0.2) +
  scale_fill_viridis(discrete = TRUE, option = "D") +
  labs(x = "Treatment", y = "Abundance - Lesser pipefish") +
  theme_cowplot(12) +
  annotate("text", x = -Inf, y = Inf, label = "p < 0.001", hjust = -0.05, vjust = 1.5, size = 4, color = "black", fontface = "bold") +
  #annotate("text", x = -Inf, y = Inf, label = "p < 0.001; Dev = 12.2%", hjust = -0.05, vjust = 1.5, size = 4, color = "black", fontface = "bold") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.title=element_blank(), 
        legend.position = "none") #+

Lesser.pipefishPlots <- plot_grid(Lesser.pipefishPlot1, Lesser.pipefishPlot2, Lesser.pipefishPlot3, nrow=1)
Lesser.pipefishPlots

resid_df <- data.frame(
  x = dat2$UTM_E,
  y = dat2$UTM_N,
  resid = residuals(best_model))
coordinates(resid_df) <- ~x+y
var_model <- variogram(resid~1, data=resid_df)
plot(var_model)
dists <- dist(cbind(dat2$UTM_E, dat2$UTM_N))
dists_matrix <- as.matrix(dists)
Moran.I(residuals(best_model), dists_matrix)

conc_result_best <- concurvity(best_model, full=FALSE) # Pairwise concurvity
formatted_result_best <- format_concurvity(conc_result_best)
print(formatted_result_best$estimate, quote = FALSE)







#### Flatfish ####

dat2 <- dat2 %>% 
  rename(Flatfish = Smooth.flatfish)

Flatfishm2a <- gam(Flatfish ~ Treatment + s(Area10m, k = 3) + s(RatioLog, k = 3) + 
                     s(DistBouldEdge, k = 3) + s(DistMussEdge, k =3) + s(UTM_E, UTM_N, bs="gp", k = 4), data = dat2, na.action = na.fail, family = nb())
summary(Flatfishm2a)
FlatfishApp <- appraise(Flatfishm2a)


# Concurvity full model
conc_result <- concurvity(Flatfishm2a, full = FALSE)
formatted_result <- format_concurvity(conc_result)
print(formatted_result$estimate, quote = FALSE)


summary(dat2$Flatfish)
plot(dat2$Flatfish ~ dat2$Area10m)
#None in Bare or 1-year-old
dat2_filtered <- dat2 %>%
  filter(!(Treatment %in% c("Bare", "1-year-old")))
dat2_filtered$Treatment <- droplevels(dat2_filtered$Treatment)
summary(dat2_filtered$Treatment)

dat2$Flatfishlog <- log((dat2$Flatfish) + 1)


## Try dredging
dredged_models <- dredge(Flatfishm2a)
print(dredged_models)
best_model <- get.models(dredged_models, 1)[[1]]
summary(best_model)
appraise(best_model)

FlatfishAppDredge <- appraise(best_model)

anova(best_model)
#Flatfish ~ s(DistMussEdge, k = 3) + s(RatioLog, k = 3) + Treatment + 1

#Parametric Terms:
#  df Chi.sq p-value
#Treatment  5  18.32 0.00257

#Approximate significance of smooth terms:
#  edf Ref.df Chi.sq p-value
#s(DistMussEdge) 1.884  1.986  6.666  0.0303
#s(RatioLog)     1.000  1.000  3.582  0.0584
emm_result <- emmeans(best_model, specs = "Treatment")
pairs_result <- pairs(emm_result, adjust = "tukey")
pairs_df <- as.data.frame(pairs_result)
cld_result <- multcomp::cld(emm_result)
cld_result

#Proportion of deviance explained:

calculate_deviance_explained <- function(best_model, reduced_formula) {
  reduced_model <- update(best_model, formula = reduced_formula)
  (deviance(reduced_model) - deviance(best_model)) / deviance(reduced_model) * 100
}

#calculate_deviance_explained(Flatfishm2a, . ~ . - s(Area10m, k = 3))
calculate_deviance_explained(best_model, . ~ . - s(RatioLog, k = 3))
#calculate_deviance_explained(Flatfishm2a, . ~ . - s(DistBouldEdge, k = 3))
calculate_deviance_explained(best_model, . ~ . - s(DistMussEdge, k = 3))
#calculate_deviance_explained(Flatfishm2a, . ~ . - s(UTM_E,UTM_N, bs="gp", k = 4))
calculate_deviance_explained(best_model, . ~ . - Treatment)

#Seperate out plots to add more info to each
#p1 <- draw(Flatfishm2a, select = "s(Area10m)", residuals = TRUE) + ggtitle(NULL) + theme_cowplot(12) +
#  annotate("text", x = -Inf, y = Inf, label = "p = 0.21; Dev = 2.2%", hjust = -0.05, vjust = 2, size = 4, color = "black") + labs(x = expression("Area of seagrass within 10m radius (m"^2*")"))
#p2 <- draw(Flatfishm2a, select = "s(RatioLog)", residuals = TRUE) + ggtitle(NULL) + theme_cowplot(12)+
#  annotate("text", x = -Inf, y = Inf, label = "p = 0.3; Dev = 1.1%", hjust = -0.05, vjust = 2, size = 4, color = "black") + labs(x = "Seagrass perimeter:area (10m radius)")
#p3 <- draw(Flatfishm2a, select = "s(DistBouldEdge)", residuals = TRUE) + ggtitle(NULL) + theme_cowplot(12) +
#  annotate("text", x = -Inf, y = Inf, label = "p = 0.67; Dev = 0.71%", hjust = -0.05, vjust = 2, size = 4, color = "black") + labs(x = "Distance to boulder reef (m)")
#p4 <- draw(Flatfishm2a, select = "s(DistMussEdge)", residuals = TRUE) + ggtitle(NULL) + theme_cowplot(12)+
#  annotate("text", x = -Inf, y = Inf, label = "p = 0.02; Dev = 7.9%", hjust = -0.05, vjust = 2, size = 4, color = "black", fontface = "bold") + labs(x = "Distance to mussel reef (m)")
#p5 <- draw(Flatfishm2a, select = "s(UTM_E,UTM_N)", residuals = TRUE) + ggtitle(NULL) + theme_cowplot(12) +
#  annotate("text", x = 542900, y = 6171620, label = "p = 0.81; \n Dev = 2.1%", size = 4, color = "black") + 
#  theme(axis.text.x = element_blank(),
#        axis.text.y = element_blank()) 

##Plotting Treatment effect

#Suppressing intercept for Treatment
#Flatfishm2bS <- gam(Flatfish ~ 0+Treatment + s(Area10m, k = 3) + s(RatioLog, k = 3) + 
#                             s(DistBouldEdge, k = 3) + s(DistMussEdge, k = 3) + s(UTM_E, UTM_N, bs="gp", k = 4), data = dat2, na.action = na.fail, family = poisson)

# Extract the partial effects and standard errors for the Treatment factor
#preds <- predict(Flatfishm2bS, type = "terms", se.fit = TRUE, terms = "Treatment")
#dat2$TreatmentEffect <- preds$fit[, "Treatment"]
#dat2$SE <- preds$se.fit[, "Treatment"]

# Calculate the confidence intervals
#dat2$CI_upper <- dat2$TreatmentEffect + 1.96 * dat2$SE
#dat2$CI_lower <- dat2$TreatmentEffect - 1.96 * dat2$SE

# Plot the Treatment effect with confidence intervals
#p6 <- ggplot(dat2, aes(x = Treatment, y = TreatmentEffect)) +
#  geom_point() +
#  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper), width = 0.2) +
#  theme_cowplot(12) +
#  labs(x = "Treatment", y = "Partial effect") +
#  annotate("text", x = -Inf, y = Inf, label = "p < 0.04; Dev = 12.0%", hjust = -0.05, vjust = 1.5, size = 4, color = "black", fontface = "bold") +
#  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#FlatfishAllPlots <- plot_grid(p1, p2, p3, p4, p5, p6, labels = c('A', 'B', 'C', 'D','E', 'F'), label_size = 12, ncol = 2)
#FlatfishAllPlots
#ggsave(FlatfishAllPlots, filename = "2024-06-31-FlatfishAllPlots.png", width = 6, height = 8)


##Plotting signficant terms
#  df    F p-value
#Treatment  5  18.32 0.00257

#Approximate significance of smooth terms:
#  edf Ref.df Chi.sq p-value
#s(DistMussEdge) 1.884  1.986  6.666  0.0303
#s(RatioLog)     1.000  1.000  3.582  0.0584

vis_smooth1 <- visreg(best_model, "DistMussEdge", scale="response",plot = FALSE)
vis_smooth2 <- visreg(best_model, "RatioLog", scale="response",plot = FALSE)

FlatfishPlot1 <- ggplot(vis_smooth1$fit, aes(x = DistMussEdge, y = visregFit)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = visregLwr, ymax = visregUpr), alpha = 0.2) +
  labs(x = expression(atop("Distance to", "mussel reef (m)")), 
       y = "Abundance - Flatfish") +
  #theme_minimal(base_size = 15)+
  annotate("text", x = -Inf, y = Inf, label = "p = 0.03", hjust = -0.05, vjust = 2, size = 4, color = "black", fontface = "bold") +
  #annotate("text", x = -Inf, y = Inf, label = "p = 0.03; Dev = -3.1%", hjust = -0.05, vjust = 2, size = 4, color = "black", fontface = "bold") +
  #scale_y_continuous(breaks = seq(0, 10, by = 2)) +
  #scale_y_continuous(limits = c(4, 10)) +
  theme_cowplot(12)

FlatfishPlot2 <- ggplot(vis_smooth2$fit, aes(x = RatioLog, y = visregFit)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = visregLwr, ymax = visregUpr), alpha = 0.2) +
  labs(x = expression(atop("Seagrass perimeter:area", "(10m radius)")), 
       , y = "Abundance - Flatfish") +  #theme_minimal(base_size = 15)+
  annotate("text", x = -Inf, y = Inf, label = "p = 0.06", hjust = -0.05, vjust = 2, size = 4, color = "black") +
  #annotate("text", x = -Inf, y = Inf, label = "p = 0.058; Dev = -4.0%", hjust = -0.05, vjust = 2, size = 4, color = "black") +
  #scale_y_continuous(breaks = seq(0, 10, by = 2)) +
  #scale_y_continuous(limits = c(4, 10)) +
  theme_cowplot(12)

vis_cat <- visreg(best_model, "Treatment", by = "Treatment", scale="response",plot = FALSE)

FlatfishPlot3 <- ggplot(vis_cat$fit, aes(x = as.factor(Treatment), y = visregFit)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = visregLwr, ymax = visregUpr), width = 0.2) +
  scale_fill_viridis(discrete = TRUE, option = "D") +
  labs(x = "Treatment", y = "Abundance - Flatfish") +
  theme_cowplot(12) +
  annotate("text", x = -Inf, y = Inf, label = "p = 0.002", hjust = -0.05, vjust = 1.5, size = 4, color = "black", fontface = "bold") +
  #annotate("text", x = -Inf, y = Inf, label = "p = 0.002; Dev = -0.5%", hjust = -0.05, vjust = 1.5, size = 4, color = "black", fontface = "bold") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.title=element_blank(), 
        legend.position = "none") +
  scale_y_continuous(limits = c(0, 2.5))

#FlatfishPlots <- plot_grid(FlatfishPlot1, FlatfishPlot2, FlatfishPlot3, labels = c('A', 'B', 'C'), label_size = 12)
FlatfishPlots <- plot_grid(FlatfishPlot1, FlatfishPlot2, FlatfishPlot3, nrow=1)
FlatfishPlots

resid_df <- data.frame(
  x = dat2$UTM_E,
  y = dat2$UTM_N,
  resid = residuals(best_model))
coordinates(resid_df) <- ~x+y
var_model <- variogram(resid~1, data=resid_df)
plot(var_model)
dists <- dist(cbind(dat2$UTM_E, dat2$UTM_N))
dists_matrix <- as.matrix(dists)
Moran.I(residuals(best_model), dists_matrix)







#### Baltic.prawn ####
dat2 <- dat2 %>% 
  rename(Baltic.prawn = Large.shrimp)

#Baltic.prawnm2a <- gam(Baltic.prawn ~ Treatment + s(Area10m, k = 3) + s(RatioLog, k = 3) + 
#                           s(DistBouldEdge, k = 3) + s(DistMussEdge, k =3) + s(UTM_E, UTM_N, bs="gp", k = 4), data = dat2, action = na.fail, family = poisson)
#summary(Baltic.prawnm2a)
#Baltic.prawnApp <- appraise(Baltic.prawnm2a)
#ggsave("Baltic.prawnApp.png", Baltic.prawnApp, width = 7, height = 5)

##Not great plots
summary(dat2$Baltic.prawn)
plot(dat2$Baltic.prawn ~ dat2$Treatment)
#None in Bare or 1-year-old
dat2_filtered <- subset(dat2, !Treatment %in% c("Bare", "1-year-old"))
summary(dat2_filtered$Treatment)

#Baltic.prawnm2a <- gam(Baltic.prawn ~ Treatment + s(Area10m, k = 3) + s(RatioLog, k = 3) + 
#                    s(DistBouldEdge, k = 3) + s(DistMussEdge, k =3) + s(UTM_E, UTM_N, bs="gp", k = 4), data = dat2_filtered, action = na.fail, family = poisson)
#summary(Baltic.prawnm2a)
#Baltic.prawnApp <- 
#appraise(Baltic.prawnm2a)
#ggsave("Baltic.prawnApp.png", Baltic.prawnApp, width = 7, height = 5)

# Use a negative binomial

Baltic.prawnm2a <- gam(Baltic.prawn ~ Treatment + s(Area10m, k = 3) + s(RatioLog, k = 3) + 
                         s(DistBouldEdge, k = 3) + s(DistMussEdge, k = 3) + 
                         s(UTM_E, UTM_N, bs="gp", k = 4), 
                       data = dat2_filtered, 
                       na.action = na.fail, 
                       family = nb())
summary(Baltic.prawnm2a)
appraise(Baltic.prawnm2a)
Baltic.prawnApp <- appraise(Baltic.prawnm2a)

draw(Baltic.prawnm2a, residuals = TRUE) & theme_cowplot(12)
options(scipen = 10)
concurvity(Baltic.prawnm2a, full = TRUE)
concurvity(Baltic.prawnm2a, full = FALSE)
gam.check(Baltic.prawnm2a)

anova(Baltic.prawnm2a)

#Proportion of deviance explained:

calculate_deviance_explained <- function(Baltic.prawnm2a, reduced_formula) {
  reduced_model <- update(Baltic.prawnm2a, formula = reduced_formula)
  (deviance(reduced_model) - deviance(Baltic.prawnm2a)) / deviance(reduced_model) * 100
}

calculate_deviance_explained(Baltic.prawnm2a, . ~ . - s(Area10m, k = 3))
calculate_deviance_explained(Baltic.prawnm2a, . ~ . - s(RatioLog, k = 3))
calculate_deviance_explained(Baltic.prawnm2a, . ~ . - s(DistBouldEdge, k = 3))
calculate_deviance_explained(Baltic.prawnm2a, . ~ . - s(DistMussEdge, k = 3))
calculate_deviance_explained(Baltic.prawnm2a, . ~ . - s(UTM_E,UTM_N, bs="gp", k = 4))
calculate_deviance_explained(Baltic.prawnm2a, . ~ . - Treatment)

#Seperate out plots to add more info to each
p1 <- draw(Baltic.prawnm2a, select = "s(Area10m)", residuals = TRUE) + ggtitle(NULL) + theme_cowplot(12) +
  annotate("text", x = -Inf, y = Inf, label = "p = 0.44; Dev = 0.9%", hjust = -0.05, vjust = 2, size = 4, color = "black") + labs(x = expression("Area of seagrass within 10m radius (m"^2*")"))
p2 <- draw(Baltic.prawnm2a, select = "s(RatioLog)", residuals = TRUE) + ggtitle(NULL) + theme_cowplot(12)+
  annotate("text", x = -Inf, y = Inf, label = "p = 0.025; Dev = 1.3%", hjust = -0.05, vjust = 2, size = 4, color = "black", fontface = "bold") + labs(x = "Seagrass perimeter:area (10m radius)")
p3 <- draw(Baltic.prawnm2a, select = "s(DistBouldEdge)", residuals = TRUE) + ggtitle(NULL) + theme_cowplot(12) +
  annotate("text", x = -Inf, y = Inf, label = "p = 0.34; Dev = 0.76%", hjust = -0.05, vjust = 2, size = 4, color = "black") + labs(x = "Distance to boulder reef (m)")
p4 <- draw(Baltic.prawnm2a, select = "s(DistMussEdge)", residuals = TRUE) + ggtitle(NULL) + theme_cowplot(12)+
  annotate("text", x = -Inf, y = Inf, label = "p = 0.77; Dev = 1.1%", hjust = -0.05, vjust = 2, size = 4, color = "black") + labs(x = "Distance to mussel reef (m)")
p5 <- draw(Baltic.prawnm2a, select = "s(UTM_E,UTM_N)", residuals = TRUE) + ggtitle(NULL) + theme_cowplot(12) +
  annotate("text", x = 542900, y = 6171620, label = "p = 0.76; \n Dev = 2.1%", size = 4, color = "black") + 
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank()) 

##Plotting Treatment effect

#Suppressing intercept for Treatment
Baltic.prawnm2bS <- gam(Baltic.prawn ~ 0+Treatment + s(Area10m, k = 3) + s(RatioLog, k = 3) + 
                          s(DistBouldEdge, k = 3) + s(DistMussEdge, k = 3) + s(UTM_E, UTM_N, bs="gp", k = 4), data = dat2_filtered, na.action = na.fail, family = poisson)

# Extract the partial effects and standard errors for the Treatment factor
preds <- predict(Baltic.prawnm2bS, type = "terms", se.fit = TRUE, terms = "Treatment")
dat2_filtered$TreatmentEffect <- preds$fit[, "Treatment"]
dat2_filtered$SE <- preds$se.fit[, "Treatment"]

# Calculate the confidence intervals
dat2_filtered$CI_upper <- dat2_filtered$TreatmentEffect + 1.96 * dat2_filtered$SE
dat2_filtered$CI_lower <- dat2_filtered$TreatmentEffect - 1.96 * dat2_filtered$SE

# Plot the Treatment effect with confidence intervals
p6 <- ggplot(dat2_filtered, aes(x = Treatment, y = TreatmentEffect)) +
  geom_point() +
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper), width = 0.2) +
  theme_cowplot(12) +
  labs(x = "Treatment", y = "Partial effect") +
  annotate("text", x = -Inf, y = Inf, label = "p = 0.44; Dev = 3.1%", hjust = -0.05, vjust = 1.5, size = 4, color = "black") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_y_continuous(limits = c(-1, 2))

Baltic.prawnAllPlots <- plot_grid(p1, p2, p3, p4, p5, p6, labels = c('A', 'B', 'C', 'D','E', 'F'), label_size = 12, ncol = 2)
Baltic.prawnAllPlots

#vis_smooth1 <- visreg(Baltic.prawnm2a, "DistBouldEdge", plot = FALSE)
vis_smooth2 <- visreg(Baltic.prawnm2a, "RatioLog", plot = FALSE)

Baltic.prawnPlot1 <- ggplot(vis_smooth2$fit, aes(x = RatioLog, y = visregFit)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = visregLwr, ymax = visregUpr), alpha = 0.2) +
  labs(x = "Seagrass perimeter:area (10m radius)", y = "Baltic prawn") +
  theme_minimal(base_size = 15)+
  annotate("text", x = -Inf, y = Inf, label = "p = 0.02; Dev = 7.9%", hjust = -0.05, vjust = 2, size = 4, color = "black", fontface = "bold") +
  #scale_y_continuous(breaks = seq(0, 10, by = 2)) +
  #scale_y_continuous(limits = c(4, 10)) +
  theme_cowplot(12)

Baltic.prawnPlot1

## Try dredging
dredged_models <- dredge(Baltic.prawnm2a)
print(dredged_models)
best_model <- get.models(dredged_models, 1)[[1]]
summary(best_model)
appraise(best_model)

Baltic.prawnAppDredge <- appraise(best_model)

anova(best_model)

#Proportion of deviance explained:

calculate_deviance_explained <- function(Baltic.prawnm2a, reduced_formula) {
  reduced_model <- update(Baltic.prawnm2a, formula = reduced_formula)
  (deviance(reduced_model) - deviance(Baltic.prawnm2a)) / deviance(reduced_model) * 100
}

#calculate_deviance_explained(Baltic.prawnm2a, . ~ . - s(Area10m, k = 3))
calculate_deviance_explained(Baltic.prawnm2a, . ~ . - s(RatioLog, k = 3))
#calculate_deviance_explained(Baltic.prawnm2a, . ~ . - s(DistBouldEdge, k = 3))
calculate_deviance_explained(Baltic.prawnm2a, . ~ . - s(DistMussEdge, k = 3))
#calculate_deviance_explained(Baltic.prawnm2a, . ~ . - s(UTM_E,UTM_N, bs="gp", k = 4))
calculate_deviance_explained(Baltic.prawnm2a, . ~ . - Treatment)


vis_smooth1 <- visreg(best_model, "DistBouldEdge",scale="response", plot = FALSE)
vis_smooth2 <- visreg(best_model, "RatioLog",scale="response", plot = FALSE)

Baltic.prawnPlot1 <- ggplot(vis_smooth1$fit, aes(x = DistBouldEdge, y = visregFit)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = visregLwr, ymax = visregUpr), alpha = 0.2) +
  labs(x = expression(atop("Distance to", "boulder reef (m)")), 
       y = "Abundance - Baltic prawn") +  #theme_minimal(base_size = 15)+
  annotate("text", x = -Inf, y = Inf, label = "p = 0.09", hjust = -0.05, vjust = 2, size = 4, color = "black") +
  #annotate("text", x = -Inf, y = Inf, label = "p = 0.09; Dev = 1.3%", hjust = -0.05, vjust = 2, size = 4, color = "black") +
  #scale_y_continuous(breaks = seq(0, 10, by = 2)) +
  #scale_y_continuous(limits = c(4, 10)) +
  theme_cowplot(12)

Baltic.prawnPlot2 <- ggplot(vis_smooth2$fit, aes(x = RatioLog, y = visregFit)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = visregLwr, ymax = visregUpr), alpha = 0.2) +
  labs(x = expression(atop("Seagrass perimeter:area", "(10m radius)")), 
       , y = "Abundance - Baltic prawn") +  #theme_minimal(base_size = 15)+
  annotate("text", x = -Inf, y = Inf, label = "p = 0.01", hjust = -0.05, vjust = 2, size = 4, color = "black", fontface = "bold") +
  #annotate("text", x = -Inf, y = Inf, label = "p = 0.01; Dev = 1.1%", hjust = -0.05, vjust = 2, size = 4, color = "black", fontface = "bold") +
  #scale_y_continuous(breaks = seq(0, 10, by = 2)) +
  #scale_y_continuous(limits = c(-2, 2)) +
  theme_cowplot(12)


Baltic.prawnPlots <- plot_grid(Baltic.prawnPlot1, Baltic.prawnPlot2)
Baltic.prawnPlots

resid_df <- data.frame(
  x = dat2_filtered$UTM_E,
  y = dat2_filtered$UTM_N,
  resid = residuals(best_model))
coordinates(resid_df) <- ~x+y
var_model <- variogram(resid~1, data=resid_df)
plot(var_model)
dists <- dist(cbind(dat2_filtered$UTM_E, dat2_filtered$UTM_N))
dists_matrix <- as.matrix(dists)
Moran.I(residuals(best_model), dists_matrix)








#### Periwinkle ####
dat2 <- dat2 %>% 
  rename(Periwinkle = Herb.snail.ab)

## negative binomial ##
Periwinklem2a <- gam(Periwinkle ~ Treatment + s(Area10m, k = 3) + s(RatioLog, k = 3) + 
                       s(DistBouldEdge, k = 3) + s(DistMussEdge, k =3) + s(UTM_E, UTM_N, bs="gp", k = 4), data = dat2, na.action = na.fail, family = nb())
summary(Periwinklem2a)
PeriwinkleApp <- appraise(Periwinklem2a)

draw(Periwinklem2a, residuals = TRUE) & theme_cowplot(12)
options(scipen = 10)
concurvity(Periwinklem2a, full = TRUE)
concurvity(Periwinklem2a, full = FALSE)
gam.check(Periwinklem2a)

anova(Periwinklem2a)
#Proportion of deviance explained:

calculate_deviance_explained <- function(Periwinklem2a, reduced_formula) {
  reduced_model <- update(Periwinklem2a, formula = reduced_formula)
  (deviance(reduced_model) - deviance(Periwinklem2a)) / deviance(reduced_model) * 100
}

calculate_deviance_explained(Periwinklem2a, . ~ . - s(Area10m, k = 3))
calculate_deviance_explained(Periwinklem2a, . ~ . - s(RatioLog, k = 3))
calculate_deviance_explained(Periwinklem2a, . ~ . - s(DistBouldEdge, k = 3))
calculate_deviance_explained(Periwinklem2a, . ~ . - s(DistMussEdge, k = 3))
calculate_deviance_explained(Periwinklem2a, . ~ . - s(UTM_E,UTM_N, bs="gp", k = 4))
calculate_deviance_explained(Periwinklem2a, . ~ . - Treatment)

#Seperate out plots to add more info to each
p1 <- draw(Periwinklem2a, select = "s(Area10m)", residuals = TRUE) + ggtitle(NULL) + theme_cowplot(12) +
  annotate("text", x = -Inf, y = Inf, label = "p = 0.53; Dev = 0.78%", hjust = -0.05, vjust = 2, size = 4, color = "black") + labs(x = expression("Area of seagrass within 10m radius (m"^2*")"))
p2 <- draw(Periwinklem2a, select = "s(RatioLog)", residuals = TRUE) + ggtitle(NULL) + theme_cowplot(12)+
  annotate("text", x = -Inf, y = Inf, label = "p = 0.75; Dev = 0.62%", hjust = -0.05, vjust = 2, size = 4, color = "black") + labs(x = "Seagrass perimeter:area (10m radius)")
p3 <- draw(Periwinklem2a, select = "s(DistBouldEdge)", residuals = TRUE) + ggtitle(NULL) + theme_cowplot(12) +
  annotate("text", x = -Inf, y = Inf, label = "p = 0.15; Dev = 1.1%", hjust = -0.05, vjust = 2, size = 4, color = "black") + labs(x = "Distance to boulder reef (m)")
p4 <- draw(Periwinklem2a, select = "s(DistMussEdge)", residuals = TRUE) + ggtitle(NULL) + theme_cowplot(12)+
  annotate("text", x = -Inf, y = Inf, label = "p = 0.10; Dev = 1.0%", hjust = -0.05, vjust = 2, size = 4, color = "black") + labs(x = "Distance to mussel reef (m)")
p5 <- draw(Periwinklem2a, select = "s(UTM_E,UTM_N)", residuals = TRUE) + ggtitle(NULL) + theme_cowplot(12) +
  annotate("text", x = 542900, y = 6171620, label = "p = 0.91; \n Dev = 1.5%", size = 4, color = "black") + 
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank()) 

##Plotting Treatment effect

#Suppressing intercept for Treatment
Periwinklem2bS <- gam(Periwinkle ~ 0+Treatment + s(Area10m, k = 3) + s(RatioLog, k = 3) + 
                        s(DistBouldEdge, k = 3) + s(DistMussEdge, k = 3) + s(UTM_E, UTM_N, bs="gp", k = 4), data = dat2, na.action = na.fail, family = poisson)

# Extract the partial effects and standard errors for the Treatment factor
preds <- predict(Periwinklem2bS, type = "terms", se.fit = TRUE, terms = "Treatment")
dat2$TreatmentEffect <- preds$fit[, "Treatment"]
dat2$SE <- preds$se.fit[, "Treatment"]

# Calculate the confidence intervals
dat2$CI_upper <- dat2$TreatmentEffect + 1.96 * dat2$SE
dat2$CI_lower <- dat2$TreatmentEffect - 1.96 * dat2$SE

# Plot the Treatment effect with confidence intervals
p6 <- ggplot(dat2, aes(x = Treatment, y = TreatmentEffect)) +
  geom_point() +
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper), width = 0.2) +
  theme_cowplot(12) +
  labs(x = "Treatment", y = "Partial effect") +
  annotate("text", x = -Inf, y = Inf, label = "p < 0.001; Dev = 3.8%", hjust = -0.05, vjust = 1.5, size = 4, color = "black", fontface = "bold") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_y_continuous(limits = c(0, 4))

PeriwinkleAllPlots <- plot_grid(p1, p2, p3, p4, p5, p6, labels = c('A', 'B', 'C', 'D','E', 'F'), label_size = 12, ncol = 2)
PeriwinkleAllPlots


#vis_smooth1 <- visreg(Periwinklem2a, "DistBouldEdge", plot = FALSE)
#vis_smooth2 <- visreg(Periwinklem2a, "DistMussEdge", plot = FALSE)
vis_cat <- visreg(Periwinklem2a, "Treatment", by = "Treatment", plot = FALSE)

PeriwinklePlot2 <- ggplot(vis_cat$fit, aes(x = as.factor(Treatment), y = visregFit)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = visregLwr, ymax = visregUpr), width = 0.2) +
  scale_fill_viridis(discrete = TRUE, option = "D") +
  labs(x = "", y = "Periwinkle") +
  theme_cowplot(12) +
  annotate("text", x = -Inf, y = Inf, label = "p < 0.001; Dev = 3.8%", hjust = -0.05, vjust = 1.5, size = 4, color = "black", fontface = "bold") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.title=element_blank(), 
        legend.position = "none") +
  scale_y_continuous(limits = c(0, 4))

## Try dredging
dredged_models <- dredge(Periwinklem2a)
print(dredged_models)
best_model <- get.models(dredged_models, 1)[[1]]
summary(best_model)
appraise(best_model)

PeriwinkleAppDredge <- appraise(best_model)

anova(best_model)

# Get estimated marginal means
emm_result <- emmeans(best_model, specs = "Treatment")
pairs_result <- pairs(emm_result, adjust = "tukey")
pairs_df <- as.data.frame(pairs_result)
cld_result <- multcomp::cld(emm_result)
cld_result


#Proportion of deviance explained:
calculate_deviance_explained <- function(best_model, reduced_formula) {
  reduced_model <- update(best_model, formula = reduced_formula)
  (deviance(reduced_model) - deviance(best_model)) / deviance(reduced_model) * 100
}

#calculate_deviance_explained(best_model, . ~ . - s(Area10m, k = 3))
#calculate_deviance_explained(best_model, . ~ . - s(RatioLog, k = 3))
calculate_deviance_explained(best_model, . ~ . - s(DistBouldEdge, k = 3))
calculate_deviance_explained(best_model, . ~ . - s(DistMussEdge, k = 3))
#calculate_deviance_explained(best_model, . ~ . - s(UTM_E,UTM_N, bs="gp", k = 4))
calculate_deviance_explained(best_model, . ~ . - Treatment)

vis_smooth1 <- visreg(best_model, "DistBouldEdge", scale="response",plot = FALSE)
vis_smooth2 <- visreg(best_model, "DistMussEdge", scale="response",plot = FALSE)

PeriwinklePlot1 <- ggplot(vis_smooth1$fit, aes(x = DistBouldEdge, y = visregFit)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = visregLwr, ymax = visregUpr), alpha = 0.2) +
  labs(x = expression(atop("Distance to", "boulder reef (m)")), 
       y = "Abundance - Periwinkle") +  #theme_minimal(base_size = 15)+
  annotate("text", x = -Inf, y = Inf, label = "p = 0.07", hjust = -0.05, vjust = 2, size = 4, color = "black") +
  #annotate("text", x = -Inf, y = Inf, label = "p = 0.07; Dev = 1.3%", hjust = -0.05, vjust = 2, size = 4, color = "black") +
  #scale_y_continuous(breaks = seq(0, 10, by = 2)) +
  #scale_y_continuous(limits = c(4, 10)) +
  theme_cowplot(12)

PeriwinklePlot2 <- ggplot(vis_smooth2$fit, aes(x = DistMussEdge, y = visregFit)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = visregLwr, ymax = visregUpr), alpha = 0.2) +
  labs(x = expression(atop("Distance to", "mussel reef (m)")), 
       y = "Abundance - Periwinkle") +  #theme_minimal(base_size = 15)+
  #annotate("text", x = -Inf, y = Inf, label = "p = 0.038; Dev = 1.6%", hjust = -0.05, vjust = 2, size = 4, color = "black", fontface = "bold") +
  annotate("text", x = -Inf, y = Inf, label = "p = 0.038", hjust = -0.05, vjust = 2, size = 4, color = "black", fontface = "bold") +
  #scale_y_continuous(breaks = seq(0, 10, by = 2)) +
  #scale_y_continuous(limits = c(-2, 2)) +
  theme_cowplot(12)

vis_cat <- visreg(best_model, "Treatment", by = "Treatment",  scale="response",plot = FALSE)

PeriwinklePlot3 <- ggplot(vis_cat$fit, aes(x = as.factor(Treatment), y = visregFit)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = visregLwr, ymax = visregUpr), width = 0.2) +
  scale_fill_viridis(discrete = TRUE, option = "D") +
  labs(x = "Treatment", y = "Abundance – Periwinkle") +
  theme_cowplot(12) +
  annotate("text", x = -Inf, y = Inf, label = "p < 0.001", hjust = -0.05, vjust = 1.5, size = 4, color = "black", fontface = "bold") +
  #annotate("text", x = -Inf, y = Inf, label = "p < 0.001; Dev = 5.1%", hjust = -0.05, vjust = 1.5, size = 4, color = "black", fontface = "bold") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.title=element_blank(), 
        legend.position = "none") +
  scale_y_continuous(limits = c(0, 30))


PeriwinklePlotsDredge <- plot_grid(PeriwinklePlot1, PeriwinklePlot2,PeriwinklePlot3, nrow=1)
PeriwinklePlotsDredge

ggsave(PeriwinklePlotsDredge, filename = "2024-07-31-PeriwinklePlotsSigDredgeResponse.png", width = 9, height = 3)

resid_df <- data.frame(
  x = dat2$UTM_E,
  y = dat2$UTM_N,
  resid = residuals(best_model))
coordinates(resid_df) <- ~x+y
var_model <- variogram(resid~1, data=resid_df)
plot(var_model)
dists <- dist(cbind(dat2$UTM_E, dat2$UTM_N))
dists_matrix <- as.matrix(dists)
Moran.I(residuals(best_model), dists_matrix)













#### Whelk ####
dat2 <- dat2 %>% 
  rename(Whelk = Whelk.ab)
## negative binomial ##
Whelkm2a <- gam(Whelk ~ Treatment + s(Area10m, k = 3) + s(RatioLog, k = 3) + 
                  s(DistBouldEdge, k = 3) + s(DistMussEdge, k =3) + s(UTM_E, UTM_N, bs="gp", k = 4), data = dat2, na.action = na.fail, family = nb())
summary(Whelkm2a)
WhelkApp <- appraise(Whelkm2a)

draw(Whelkm2a, residuals = TRUE) & theme_cowplot(12)
options(scipen = 10)
concurvity(Whelkm2a, full = TRUE)
concurvity(Whelkm2a, full = FALSE)
gam.check(Whelkm2a)

anova(Whelkm2a)

#Proportion of deviance explained:

calculate_deviance_explained <- function(Whelkm2a, reduced_formula) {
  reduced_model <- update(Whelkm2a, formula = reduced_formula)
  (deviance(reduced_model) - deviance(Whelkm2a)) / deviance(reduced_model) * 100
}

calculate_deviance_explained(Whelkm2a, . ~ . - s(Area10m, k = 3))
calculate_deviance_explained(Whelkm2a, . ~ . - s(RatioLog, k = 3))
calculate_deviance_explained(Whelkm2a, . ~ . - s(DistBouldEdge, k = 3))
calculate_deviance_explained(Whelkm2a, . ~ . - s(DistMussEdge, k = 3))
calculate_deviance_explained(Whelkm2a, . ~ . - s(UTM_E,UTM_N, bs="gp", k = 4))
calculate_deviance_explained(Whelkm2a, . ~ . - Treatment)

#Seperate out plots to add more info to each
#p1 <- draw(Whelkm2a, select = "s(Area10m)", residuals = TRUE) + ggtitle(NULL) + theme_cowplot(12) +
#  annotate("text", x = -Inf, y = Inf, label = "p = 0.53; Dev = 0.78%", hjust = -0.05, vjust = 2, size = 4, color = "black") + labs(x = expression("Area of seagrass within 10m radius (m"^2*")"))
#p2 <- draw(Whelkm2a, select = "s(RatioLog)", residuals = TRUE) + ggtitle(NULL) + theme_cowplot(12)+
#  annotate("text", x = -Inf, y = Inf, label = "p = 0.75; Dev = 0.62%", hjust = -0.05, vjust = 2, size = 4, color = "black") + labs(x = "Seagrass perimeter:area (10m radius)")
#p3 <- draw(Whelkm2a, select = "s(DistBouldEdge)", residuals = TRUE) + ggtitle(NULL) + theme_cowplot(12) +
#  annotate("text", x = -Inf, y = Inf, label = "p = 0.15; Dev = 1.1%", hjust = -0.05, vjust = 2, size = 4, color = "black") + labs(x = "Distance to boulder reef (m)")
#p4 <- draw(Whelkm2a, select = "s(DistMussEdge)", residuals = TRUE) + ggtitle(NULL) + theme_cowplot(12)+
#  annotate("text", x = -Inf, y = Inf, label = "p = 0.10; Dev = 1.0%", hjust = -0.05, vjust = 2, size = 4, color = "black") + labs(x = "Distance to mussel reef (m)")
#p5 <- draw(Whelkm2a, select = "s(UTM_E,UTM_N)", residuals = TRUE) + ggtitle(NULL) + theme_cowplot(12) +
#  annotate("text", x = 542900, y = 6171620, label = "p = 0.91; \n Dev = 1.5%", size = 4, color = "black") + 
#  theme(axis.text.x = element_blank(),
#        axis.text.y = element_blank()) 

##Plotting Treatment effect

#Suppressing intercept for Treatment
#Whelkm2bS <- gam(Whelk ~ 0+Treatment + s(Area10m, k = 3) + s(RatioLog, k = 3) + 
#                       s(DistBouldEdge, k = 3) + s(DistMussEdge, k = 3) + s(UTM_E, UTM_N, bs="gp", k = 4), data = dat2, na.action = na.fail, family = poisson)

# Extract the partial effects and standard errors for the Treatment factor
#preds <- predict(Whelkm2bS, type = "terms", se.fit = TRUE, terms = "Treatment")
#dat2$TreatmentEffect <- preds$fit[, "Treatment"]
#dat2$SE <- preds$se.fit[, "Treatment"]

# Calculate the confidence intervals
#dat2$CI_upper <- dat2$TreatmentEffect + 1.96 * dat2$SE
#dat2$CI_lower <- dat2$TreatmentEffect - 1.96 * dat2$SE

# Plot the Treatment effect with confidence intervals
#p6 <- ggplot(dat2, aes(x = Treatment, y = TreatmentEffect)) +
# geom_point() +
#  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper), width = 0.2) +
#  theme_cowplot(12) +
#  labs(x = "Treatment", y = "Partial effect") +
#  annotate("text", x = -Inf, y = Inf, label = "p < 0.001; Dev = 3.8%", hjust = -0.05, vjust = 1.5, size = 4, color = "black", fontface = "bold") +
#  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
#  scale_y_continuous(limits = c(0, 4))

#WhelkAllPlots <- plot_grid(p1, p2, p3, p4, p5, p6, labels = c('A', 'B', 'C', 'D','E', 'F'), label_size = 12, ncol = 2)
#WhelkAllPlots
#ggsave(WhelkAllPlots, filename = "2024-06-31-WhelkAllPlots.png", width = 6, height = 8)


#vis_smooth1 <- visreg(Whelkm2a, "DistBouldEdge", plot = FALSE)
#vis_smooth2 <- visreg(Whelkm2a, "DistMussEdge", plot = FALSE)
#vis_cat <- visreg(Whelkm2a, "Treatment", by = "Treatment", plot = FALSE)

#WhelkPlot2 <- ggplot(vis_cat$fit, aes(x = as.factor(Treatment), y = visregFit)) +
#  geom_bar(stat = "identity", position = "dodge") +
#  geom_errorbar(aes(ymin = visregLwr, ymax = visregUpr), width = 0.2) +
#  scale_fill_viridis(discrete = TRUE, option = "D") +
#  labs(x = "", y = "Whelk") +
#  theme_cowplot(12) +
#  annotate("text", x = -Inf, y = Inf, label = "p < 0.001; Dev = 3.8%", hjust = -0.05, vjust = 1.5, size = 4, color = "black", fontface = "bold") +
# theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
# theme(legend.title=element_blank(), 
#        legend.position = "none") +
#  scale_y_continuous(limits = c(0, 4))

#ggsave(WhelkPlot2, filename = "2024-07-31-WhelkPlotsSig.png", width = 3, height = 3)


## Try dredging
dredged_models <- dredge(Whelkm2a)
print(dredged_models)
best_model <- get.models(dredged_models, 1)[[1]]
summary(best_model)
appraise(best_model)

WhelkAppDredge <- appraise(best_model)

anova(best_model)

calculate_deviance_explained <- function(best_model, reduced_formula) {
  reduced_model <- update(best_model, formula = reduced_formula)
  (deviance(reduced_model) - deviance(best_model)) / deviance(reduced_model) * 100
}

calculate_deviance_explained(best_model, . ~ . - s(Area10m, k = 3))
#calculate_deviance_explained(best_model, . ~ . - s(RatioLog, k = 3))
#calculate_deviance_explained(best_model, . ~ . - s(DistBouldEdge, k = 3))
#calculate_deviance_explained(best_model, . ~ . - s(DistMussEdge, k = 3))
calculate_deviance_explained(best_model, . ~ . - s(UTM_E,UTM_N,  bs="gp", k = 4))
#calculate_deviance_explained(best_model, . ~ . - Treatment)

vis_smooth1 <- visreg(best_model, "Area10m", scale="response",plot = FALSE)
#vis_smooth2 <- visreg(best_model, "DistMussEdge", plot = FALSE)

WhelkPlot1 <- ggplot(vis_smooth1$fit, aes(x = Area10m, y = visregFit)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = visregLwr, ymax = visregUpr), alpha = 0.2) +
  labs(x = expression(atop("Area of seagrass", "within 10m radius (m"^2*")")), 
       y = "Abundance - Whelk") + #theme_minimal(base_size = 15)+
  annotate("text", x = -Inf, y = Inf, label = "p = 0.05", hjust = -0.05, vjust = 2, size = 4, color = "black") +
  #annotate("text", x = -Inf, y = Inf, label = "p = 0.05; Dev = 0.6%", hjust = -0.05, vjust = 2, size = 4, color = "black") +
  #scale_y_continuous(breaks = seq(0, 10, by = 2)) +
  #scale_y_continuous(limits = c(4, 10)) +
  theme_cowplot(12)

Whelkp5 <- draw(best_model, select = "s(UTM_E,UTM_N)", residuals = TRUE) + ggtitle(NULL) + theme_cowplot(12) +
  annotate("text", x = 542900, y = 6171620, label = "p = 0.002", size = 4, color = "black", fontface = "bold") + 
  #annotate("text", x = 542900, y = 6171620, label = "p = 0.002; \n Dev = 1.6%", size = 4, color = "black", fontface = "bold") + 
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank()) +
  labs(fill = "Partial effect - Whelk",
       color = "Partial effect - Whelk")


WhelkPlotsDredge <- plot_grid(WhelkPlot1, Whelkp5)
WhelkPlotsDredge

resid_df <- data.frame(
  x = dat2$UTM_E,
  y = dat2$UTM_N,
  resid = residuals(best_model))
coordinates(resid_df) <- ~x+y
var_model <- variogram(resid~1, data=resid_df)
plot(var_model)
dists <- dist(cbind(dat2$UTM_E, dat2$UTM_N))
dists_matrix <- as.matrix(dists)
Moran.I(residuals(best_model), dists_matrix)








#### Crangon.prawn ####
dat2 <- dat2 %>% 
  rename(Crangon.prawn = Other.shrimp)
## negative binomial ##
Crangon.prawnm2a <- gam(Crangon.prawn ~ Treatment + s(Area10m, k = 3) + s(RatioLog, k = 3) + 
                          s(DistBouldEdge, k = 3) + s(DistMussEdge, k =3) + s(UTM_E, UTM_N, bs="gp", k = 4), data = dat2, na.action = na.fail, family = nb())
summary(Crangon.prawnm2a)
Crangon.prawnApp <- appraise(Crangon.prawnm2a)

draw(Crangon.prawnm2a, residuals = TRUE) & theme_cowplot(12)
options(scipen = 10)
concurvity(Crangon.prawnm2a, full = TRUE)
concurvity(Crangon.prawnm2a, full = FALSE)
gam.check(Crangon.prawnm2a)

anova(Crangon.prawnm2a)

#Proportion of deviance explained:

calculate_deviance_explained <- function(Crangon.prawnm2a, reduced_formula) {
  reduced_model <- update(Crangon.prawnm2a, formula = reduced_formula)
  (deviance(reduced_model) - deviance(Crangon.prawnm2a)) / deviance(reduced_model) * 100
}

calculate_deviance_explained(Crangon.prawnm2a, . ~ . - s(Area10m, k = 3))
calculate_deviance_explained(Crangon.prawnm2a, . ~ . - s(RatioLog, k = 3))
calculate_deviance_explained(Crangon.prawnm2a, . ~ . - s(DistBouldEdge, k = 3))
calculate_deviance_explained(Crangon.prawnm2a, . ~ . - s(DistMussEdge, k = 3))
calculate_deviance_explained(Crangon.prawnm2a, . ~ . - s(UTM_E,UTM_N, bs="gp", k = 4))
calculate_deviance_explained(Crangon.prawnm2a, . ~ . - Treatment)

#Seperate out plots to add more info to each
#p1 <- draw(Crangon.prawnm2a, select = "s(Area10m)", residuals = TRUE) + ggtitle(NULL) + theme_cowplot(12) +
#  annotate("text", x = -Inf, y = Inf, label = "p = 0.53; Dev = 0.78%", hjust = -0.05, vjust = 2, size = 4, color = "black") + labs(x = expression("Area of seagrass within 10m radius (m"^2*")"))
#p2 <- draw(Crangon.prawnm2a, select = "s(RatioLog)", residuals = TRUE) + ggtitle(NULL) + theme_cowplot(12)+
#  annotate("text", x = -Inf, y = Inf, label = "p = 0.75; Dev = 0.62%", hjust = -0.05, vjust = 2, size = 4, color = "black") + labs(x = "Seagrass perimeter:area (10m radius)")
#p3 <- draw(Crangon.prawnm2a, select = "s(DistBouldEdge)", residuals = TRUE) + ggtitle(NULL) + theme_cowplot(12) +
#  annotate("text", x = -Inf, y = Inf, label = "p = 0.15; Dev = 1.1%", hjust = -0.05, vjust = 2, size = 4, color = "black") + labs(x = "Distance to boulder reef (m)")
#p4 <- draw(Crangon.prawnm2a, select = "s(DistMussEdge)", residuals = TRUE) + ggtitle(NULL) + theme_cowplot(12)+
#  annotate("text", x = -Inf, y = Inf, label = "p = 0.10; Dev = 1.0%", hjust = -0.05, vjust = 2, size = 4, color = "black") + labs(x = "Distance to mussel reef (m)")
#p5 <- draw(Crangon.prawnm2a, select = "s(UTM_E,UTM_N)", residuals = TRUE) + ggtitle(NULL) + theme_cowplot(12) +
#  annotate("text", x = 542900, y = 6171620, label = "p = 0.91; \n Dev = 1.5%", size = 4, color = "black") + 
#  theme(axis.text.x = element_blank(),
#        axis.text.y = element_blank()) 

##Plotting Treatment effect

#Suppressing intercept for Treatment
#Crangon.prawnm2bS <- gam(Crangon.prawn ~ 0+Treatment + s(Area10m, k = 3) + s(RatioLog, k = 3) + 
#                       s(DistBouldEdge, k = 3) + s(DistMussEdge, k = 3) + s(UTM_E, UTM_N, bs="gp", k = 4), data = dat2, na.action = na.fail, family = poisson)

# Extract the partial effects and standard errors for the Treatment factor
#preds <- predict(Crangon.prawnm2bS, type = "terms", se.fit = TRUE, terms = "Treatment")
#dat2$TreatmentEffect <- preds$fit[, "Treatment"]
#dat2$SE <- preds$se.fit[, "Treatment"]

# Calculate the confidence intervals
#dat2$CI_upper <- dat2$TreatmentEffect + 1.96 * dat2$SE
#dat2$CI_lower <- dat2$TreatmentEffect - 1.96 * dat2$SE

# Plot the Treatment effect with confidence intervals
#p6 <- ggplot(dat2, aes(x = Treatment, y = TreatmentEffect)) +
# geom_point() +
#  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper), width = 0.2) +
#  theme_cowplot(12) +
#  labs(x = "Treatment", y = "Partial effect") +
#  annotate("text", x = -Inf, y = Inf, label = "p < 0.001; Dev = 3.8%", hjust = -0.05, vjust = 1.5, size = 4, color = "black", fontface = "bold") +
#  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
#  scale_y_continuous(limits = c(0, 4))

#Crangon.prawnAllPlots <- plot_grid(p1, p2, p3, p4, p5, p6, labels = c('A', 'B', 'C', 'D','E', 'F'), label_size = 12, ncol = 2)
#Crangon.prawnAllPlots
#ggsave(Crangon.prawnAllPlots, filename = "2024-06-31-Crangon.prawnAllPlots.png", width = 6, height = 8)


#vis_smooth1 <- visreg(Crangon.prawnm2a, "DistBouldEdge", plot = FALSE)
#vis_smooth2 <- visreg(Crangon.prawnm2a, "DistMussEdge", plot = FALSE)
#vis_cat <- visreg(Crangon.prawnm2a, "Treatment", by = "Treatment", plot = FALSE)

#Crangon.prawnPlot2 <- ggplot(vis_cat$fit, aes(x = as.factor(Treatment), y = visregFit)) +
#  geom_bar(stat = "identity", position = "dodge") +
#  geom_errorbar(aes(ymin = visregLwr, ymax = visregUpr), width = 0.2) +
#  scale_fill_viridis(discrete = TRUE, option = "D") +
#  labs(x = "", y = "Crangon.prawn") +
#  theme_cowplot(12) +
#  annotate("text", x = -Inf, y = Inf, label = "p < 0.001; Dev = 3.8%", hjust = -0.05, vjust = 1.5, size = 4, color = "black", fontface = "bold") +
# theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
# theme(legend.title=element_blank(), 
#        legend.position = "none") +
#  scale_y_continuous(limits = c(0, 4))

#ggsave(Crangon.prawnPlot2, filename = "2024-07-31-Crangon.prawnPlotsSig.png", width = 3, height = 3)


## Try dredging
dredged_models <- dredge(Crangon.prawnm2a)
print(dredged_models)
best_model <- get.models(dredged_models, 1)[[1]]
summary(best_model)
appraise(best_model)

Crangon.prawnAppDredge <- appraise(best_model)

anova(best_model)

# Get estimated marginal means
emm_result <- emmeans(best_model, specs = "Treatment")
pairs_result <- pairs(emm_result, adjust = "tukey")
pairs_df <- as.data.frame(pairs_result)
cld_result <- multcomp::cld(emm_result)
cld_result


calculate_deviance_explained <- function(best_model, reduced_formula) {
  reduced_model <- update(best_model, formula = reduced_formula)
  (deviance(reduced_model) - deviance(best_model)) / deviance(reduced_model) * 100
}

#calculate_deviance_explained(best_model, . ~ . - s(Area10m, k = 3))
#calculate_deviance_explained(best_model, . ~ . - s(RatioLog, k = 3))
calculate_deviance_explained(best_model, . ~ . - s(DistBouldEdge, k = 3))
calculate_deviance_explained(best_model, . ~ . - s(DistMussEdge, k = 3))
calculate_deviance_explained(best_model, . ~ . - s(UTM_E,UTM_N, bs="gp", k = 4))
calculate_deviance_explained(best_model, . ~ . - Treatment)

vis_smooth1 <- visreg(best_model, "DistBouldEdge", scale="response",plot = FALSE)
vis_smooth2 <- visreg(best_model, "DistMussEdge", scale="response",plot = FALSE)


CrangonPlot1 <- ggplot(vis_smooth1$fit, aes(x = DistBouldEdge, y = visregFit)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = visregLwr, ymax = visregUpr), alpha = 0.2) +
  labs(x = expression(atop("Distance to", "boulder reef (m)")), 
       y = "Abundance - Brown shrimp") +  theme_minimal(base_size = 15)+
  annotate("text", x = -Inf, y = Inf, label = "p = 0.06", hjust = -0.05, vjust = 2, size = 4, color = "black") +
  #annotate("text", x = -Inf, y = Inf, label = "p = 0.06; Dev = %", hjust = -0.05, vjust = 2, size = 4, color = "black") +
  #scale_y_continuous(breaks = seq(0, 10, by = 2)) +
  #scale_y_continuous(limits = c(4, 10)) +
  theme_cowplot(12)

CrangonPlot2 <- ggplot(vis_smooth2$fit, aes(x = DistMussEdge, y = visregFit)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = visregLwr, ymax = visregUpr), alpha = 0.2) +
  labs(x = expression(atop("Distance to", "mussel reef (m)")), 
       y = "Abundance - Brown shrimp") +  #theme_minimal(base_size = 15)+
  annotate("text", x = -Inf, y = Inf, label = "p = 0.001", hjust = -0.05, vjust = 2, size = 4, color = "black", fontface = "bold") +
  #annotate("text", x = -Inf, y = Inf, label = "p = 0.001; Dev = %", hjust = -0.05, vjust = 2, size = 4, color = "black") +
  #scale_y_continuous(breaks = seq(0, 10, by = 2)) +
  #scale_y_continuous(limits = c(4, 10)) +
  theme_cowplot(12)

vis_cat <- visreg(best_model, "Treatment", by = "Treatment", scale="response", plot = FALSE)

CrangonPlot3 <- ggplot(vis_cat$fit, aes(x = as.factor(Treatment), y = visregFit)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = visregLwr, ymax = visregUpr), width = 0.2) +
  scale_fill_viridis(discrete = TRUE, option = "D") +
  labs(x = "Treatment", y = "Abundance - Brown shrimp") +
  theme_cowplot(12) +
  annotate("text", x = -Inf, y = Inf, label = "p < 0.001", hjust = -0.05, vjust = 1.5, size = 4, color = "black", fontface = "bold") +
  #annotate("text", x = -Inf, y = Inf, label = "p < 0.001; Dev = 5.1%", hjust = -0.05, vjust = 1.5, size = 4, color = "black", fontface = "bold") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.title=element_blank(), 
        legend.position = "none") +
  scale_y_continuous(limits = c(0, 20))

CrangonPlot4 <- draw(best_model, select = "s(UTM_E,UTM_N)", residuals = TRUE) + ggtitle(NULL) + theme_cowplot(12) +
  annotate("text", x = 542900, y = 6171620, label = "p = 0.002", size = 4, color = "black", fontface = "bold") + 
  #annotate("text", x = 542900, y = 6171620, label = "p = 0.002; \n Dev = 1.6%", size = 4, color = "black", fontface = "bold") + 
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank()) 


Crangon.prawnPlotsDredge <- plot_grid(CrangonPlot1, CrangonPlot2, CrangonPlot3, CrangonPlot4, nrow=1)
Crangon.prawnPlotsDredge

resid_df <- data.frame(
  x = dat2$UTM_E,
  y = dat2$UTM_N,
  resid = residuals(best_model))
coordinates(resid_df) <- ~x+y
var_model <- variogram(resid~1, data=resid_df)
plot(var_model)
dists <- dist(cbind(dat2$UTM_E, dat2$UTM_N))
dists_matrix <- as.matrix(dists)
Moran.I(residuals(best_model), dists_matrix)
