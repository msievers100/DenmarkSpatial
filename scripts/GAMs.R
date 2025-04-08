## Denmark GAMs
## Sievers and Brown
## 07 Aug 2024

library(mgcv)
library(visreg)
library(gratia)

#setwd("~/Library/CloudStorage/OneDrive-GriffithUniversity/DenmarkHmscScenarioModelling/DenmarkHmscScenarioModelling")

df <- read.csv("data/datCABMS.csv")
str(df)
df$Treatment <- as.factor(df$Treatment)
df$Treatment <- relevel(df$Treatment, ref = "Bare")

colnames(df)
###All potential reponse variables:
[13] "Black.goby"              
[14] "Brittlestar"             
[15] "Crab"                    
[16] "Green.crab"              
[17] "Hermit.crab"             
[18] "Isopod"                  
[19] "Baltic.prawn"            
[20] "Lesser.pipefish"         
[21] "Mysid"                   
[22] "Crangon.prawn"           
[23] "Rock.gunnel"             
[24] "Sand.goby"               
[25] "Seastar"                 
[26] "Shortfin.sculpin"        
[27] "Flatfish"                
[28] "Spider.crab"             
[29] "Spinach.spinach"         
[30] "Straightnosed.pipefish"  
[31] "Three.spined.stickleback"
[32] "Tiger"                   
[33] "Urchin"                  
[34] "Periwinkle"              
[35] "Whelk"                   
[36] "Richness"                
[37] "Shannon"                 
[38] "Simpson"                 
[39] "InvSimpson"              
[40] "Evenness"                
[41] "DeltaPlus"         

## Realistically, it's probably only worth looking at:
[19] "Baltic.prawn"            
[20] "Lesser.pipefish"         
[22] "Crangon.prawn"           
[27] "Flatfish"                
[30] "Straightnosed.pipefish"  
[34] "Periwinkle"              
[35] "Whelk"                   
[36] "Richness"                
[37] "Shannon"                 



#### Lesser Pipefish as example ####
# Sticking with calling it gam2, for now, as unsure best way to do this given potential loop for other species/response vars

gam2 <- gam(Lesser.pipefish ~ s(DistBouldEdge, k = 3) + s(DistMussEdge, k = 3) +
              s(Area10m, k = 3) + Treatment, data = df, family = nb())
summary(gam2)

#Save GAM
saveRDS(gam2, file = "outputs/gam2.rds")

#Interpreting covariate relationship with response. Not using 'draw', to include Treatment
plot_DistBouldEdge <- visreg(gam2, "DistBouldEdge", gg = TRUE) + ggtitle("DistBouldEdge")
plot_DistMussEdge <- visreg(gam2, "DistMussEdge", gg = TRUE) + ggtitle("DistMussEdge")
plot_Area10m <- visreg(gam2, "Area10m", gg = TRUE) + ggtitle("Area10m")
plot_Treatment <- visreg(gam2, "Treatment", gg = TRUE) + ggtitle("Treatment")

combined_plot <- (plot_DistBouldEdge | plot_DistMussEdge) / (plot_Area10m | plot_Treatment)
ggsave("outputs/plots/Pipefish_visreg_plots.png", plot = combined_plot, width = 8, height = 7, dpi = 150)

#Proportion of deviance explained:
# Define a function to calculate deviance explained by each term
calculate_deviance_explained <- function(gam2, reduced_formula) {
  reduced_model <- update(gam2, formula = reduced_formula)
  (deviance(reduced_model) - deviance(gam2)) / deviance(reduced_model) * 100
}

# Calculate the deviance explained by each term
calculate_deviance_explained(gam2, . ~ . - s(DistBouldEdge, k = 3))
calculate_deviance_explained(gam2, . ~ . - s(DistMussEdge, k = 3))
calculate_deviance_explained(gam2, . ~ . - s(Area10m, k = 3))
calculate_deviance_explained(gam2, . ~ . - Treatment)


