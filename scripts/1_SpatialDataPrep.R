## Denmark spatial data prep
## Sievers and Brown


library(ggplot2)
library(mgcv)
library(tidyr)
library(patchwork)
library(sf)
library(tmap)
library(dplyr)
library(viridis)
library(visreg)
library(patchwork)

# Specie survey data
df <- read.csv("data/datCABMS.csv")

survey_area <- 3.6 #m^2

#
# Maps of study areas
#
# Load the shapefile
boulder <- st_read("data/BoulderReefPoly.shp")
saveRDS(boulder, file = "outputs/boulder.rds")
mussel <- st_read("data/MusselReef.shp")
## I had to dissolve the mussel layer for it to work for the scenarios
mussel_dissolved <- st_union(mussel)
saveRDS(mussel_dissolved, file = "outputs/mussel_dissolved.rds")

newbould <- st_read("data/NewBould.shp")
saveRDS(newbould, file = "outputs/newbould.rds")

seagrass <- st_read("data/SeagrassRestNat.shp")

tm_shape(filter(seagrass, Assigned_c == "Transplanted_eelgrass")) + 
  tm_polygons(col = "green")

tm_shape(boulder) + 
  tm_polygons(col = "red")
tm_shape(mussel_dissolved) +
  tm_polygons(col = "blue") 

eelgrass_transplanted <- filter(seagrass, Assigned_c == "Transplanted_eelgrass")

sites <- st_read("data/Sites.shp")
sites_transformed <- st_transform(sites, crs = 25832)

# Calculate the minimum convex polygon (MCP) around the sites + 10m
sites_mcp <- st_convex_hull(st_union(sites_transformed))
sites_mcp_expanded <- st_buffer(sites_mcp, dist = 10)

#Crop to boundaries of polygons - just using all seagrass to capture larger area to start with
region_bbox <- st_bbox(sites_mcp_expanded)

npts <- 100
xpts <- seq(region_bbox$xmin, region_bbox$xmax, length.out = npts)
ypts <- seq(region_bbox$ymin, region_bbox$ymax, length.out = npts)
grid <- expand.grid(xpts, ypts)
colnames(grid) <- c("x", "y")
#convert to spatial object
grid <- st_as_sf(grid, coords = c("x", "y"), crs = st_crs(boulder))

tm_shape(grid) + 
  tm_dots() +
  tm_shape(sites_mcp_expanded) +
  tm_polygons(alpha = 0.4)

#calculate grid area as dimensions of study area divided by number of grid points
grid_area <- (region_bbox$xmax - region_bbox$xmin) * (region_bbox$ymax - region_bbox$ymin) / npts^2
grid_area

saveRDS(grid_area, file = "outputs/grid_area.rds")
grid_area_ratio <- grid_area / survey_area
saveRDS(grid_area_ratio, file = "outputs/grid_area_ratio.rds")
#20.33m^2 for npts = 100
#9m^2 for npts = 150
#3.25 for npts = 250
#use this to standardize model predictions (which will be animals per survey area...)

#check if grid points are on boulder or mussel reefs and remove them
grid$OnBoulder <- st_intersects(grid, boulder, sparse = FALSE) %>% apply(1, any)
grid$OnMussel <- st_intersects(grid, mussel, sparse = FALSE) %>% apply(1, any)
#remove points on boulder or mussel
grid2 <- grid[!grid$OnBoulder & !grid$OnMussel, ]

# Crop the grid to only include points within the sites_mcp convex hull + 10m 
grid2_cropped <- st_intersection(grid2, sites_mcp_expanded)
print(nrow(grid2_cropped))
#3973 for npts = 100
#9007 for npts = 150
plot(st_geometry(grid2_cropped), col = 'blue', main = 'Cropped Grid within MCP')
plot(st_geometry(sites_mcp), border = 'red', add = TRUE)

#plot grid2 over boulder layer
tm_shape(sites_mcp_expanded) +
  tm_polygons() +
  tm_shape(boulder) +
  tm_polygons(col = "red") +
  tm_shape(mussel) +
  tm_polygons(col = "blue") +
  tm_shape(grid2_cropped) +
  tm_dots(col = "green", size = 0.05)

#
# Distance layers 
#

#calculate distance of each grid point of boulder and mussel reefs
nearest_boulder <- st_nearest_feature(grid2_cropped, boulder)

# Calculate the distance to the nearest boulder feature
grid2_cropped$DistBouldEdge <- st_distance(grid2_cropped, boulder[nearest_boulder, ], by_element = TRUE) %>%
  as.numeric()

#check, plot points coloured by distance
tm_shape(boulder, bbox = grid2_cropped) +
  tm_polygons(col = "black") +
  tm_shape(grid2_cropped) +
  tm_dots(col = "DistBouldEdge", size = 0.25, palette = viridis(100),  style = "cont")

#repeat for mussel reef
nearest_mussel <- st_nearest_feature(grid2_cropped, mussel_dissolved)
grid2_cropped$DistMussEdge <- st_distance(grid2_cropped, mussel_dissolved[nearest_mussel, ], by_element = TRUE) %>%
  as.numeric()

#check, plot points coloured by distance
tm_shape(mussel_dissolved, bbox = grid2_cropped) +
  tm_polygons(col = "black") +
  tm_shape(grid2_cropped) +
  tm_dots(col = "DistMussEdge", size = 0.25, palette = viridis(100),  style = "cont")

#locations of natural seagrass
seagrass_natural <- filter(seagrass, Assigned_c == "Eelgrass")

#get grids in natural eelgrass
grid2_cropped$InSeagrass <- st_intersects(grid2_cropped, seagrass_natural, sparse = FALSE) %>% apply(1, any)
nrow(grid2_cropped)
#plot InSeagrass to check it, green for grids with seagrass
tm_shape(grid2_cropped) +
  tm_dots(col = "InSeagrass", palette = c("red", "green"), size = 0.5) + 
  tm_shape(seagrass_natural) +
  tm_polygons(alpha = 0.2) +
  tm_shape(eelgrass_transplanted) +
  tm_polygons(alpha = 0.3, col = "black")

#Area of seagrass in 10m buffer around each grid point
buffers <- st_buffer(grid2_cropped, dist = 10)
saveRDS(buffers, file = "outputs/buffers.rds")

# Add an identifier to each buffer
buffers$id <- seq_len(nrow(buffers))
#speeding this up by doing point intersection instead. Will just be less accurate -  intersect the buffers just with the 
#grid pionts and count the number of grid points with SG per buffer, then multiply by grid area. 
#This isn't quite as accurate (will be more accurate if you increase number grid points), but its also 100s X faster
nseagrass <- st_intersects(buffers, grid2_cropped[grid2_cropped$InSeagrass,]) %>%
  lapply(length) %>% unlist()
grid2_cropped$Area10m <- nseagrass * grid_area
summary(grid2_cropped$Area10m)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.00    0.00    0.00   42.35   60.98  264.27 

#note this needs to be dynamically updated with restoration

tm_shape(grid2_cropped) +
  tm_dots(col = "Area10m", style = "cont", palette = "viridis", size =0.25)
#check seagrass appears to be in right places 

#
#Set-up the treatment factor
#

grid2_cropped$Treatment <- "Bare"
#grids in natural seagrass get the treatment
grid2_cropped$Treatment[grid2_cropped$InSeagrass] <- "Natural"
grid2_cropped$Treatment <- factor(grid2_cropped$Treatment, levels = levels(factor(df$Treatment)))

tm_shape(grid2_cropped) +
  tm_dots(col = "Treatment", size =0.25)

igrids_no_seagrass <- which(grid2_cropped$Treatment == "Bare")
ngrids_no_seagrass <- length(igrids_no_seagrass)

saveRDS(ngrids_no_seagrass, file = "outputs/ngrids_no_seagrass.rds")
saveRDS(grid2_cropped, file = "outputs/grid2_cropped.rds")
