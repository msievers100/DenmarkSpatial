# Run 3_Greedy-optimisation.R from here
# Use to select species to run the greedy optimisation on
#
# This implements two options for the optimisation
# 1. 'simple_greedy' - calculates the top N sites for your objective function. Then iteratively adds 
# transplants in that order and recalculates the objective function.
# 2. 'iterative_greedy' - As for 1, except we calculate the top site, add transplants, recalculate Area10m, then recalculate
# the next best site and so on. 
# Our expectation is that 2 would be more optimal than 1. In practice they are almost identical, with method 2 slightly better

run_date <- "2025-03-13"

#Read in GAM - note Lesser Pipefish for example...
gam_list <- readRDS("outputs/gam-selected-species.rds")

df <- read.csv("data/datCABMS.csv")
species_names <- names(gam_list)
nspp <- length(species_names)

species_names
#Species to sum for the objective function. 
ispp_fisheries <- c(2, 4)
ispp_SGfunction <- c(1, 2, 5)
ispp_all <- 1:length(species_names)
ispp_pipe <- 1

baltic <- 2
browns <- 3
flatfish <- 4
peri <- 5
whelk <- 6

Figure2 <- c(1, 4, 6)

#############################################
############ PARAMETERS TO SET ############
#############################################

# ****** SELECT SPECIES GROUP HERE ******* #
ispp_obj <- Figure2
obj_name <- "Figure2" #set the name for the output files

update_area10m <- TRUE 
#update the Area10m for each cell after each transplant then recalculate marginal values
# This is the simple_greedy
#update_area10m <- FALSE #will not update the Area10m for each cell after each transplant and therefore the
# this is iterative_greedy
# order of priorities is fixed based on initial marginal values

# **************************************** #

#save all the scenario maps (will slow it down)? 
#not recommended to do if you are running all species
# if set to TRUE, recommending just running loop for 1 species
doplots <- FALSE

#
# Run the optimisation
#

source("scripts/Greedy-optimisation-new.R")
