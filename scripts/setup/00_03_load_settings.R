# Set to location of Dropbox and GitHub folders
path_dropbox = "INSERT PATH TO DROPBOX FOLDER HERE"
path_github = "INSERT PATH TO GITHUB REPO HERE"








if (!("conus" %in% ls())) {
  conus = T  # set to T for CONUS, F for regional, this and below governs loading all necessary globals
}
if (!("pm_type" %in% ls())) {
  pm_type = "mean" # set to "mean" for average PM and "extremeN" for extreme
}

# global vars you need for function calls
fixef_rm = "perfect" # removes observations where fixed effect results in perfect fit (e.g. only 0s for poisson)
tol = 1e-8 # convergence tolerance for algorithm
corner_push = 0.95 # how far away to push algorithm as it approaches corner solution
break_alpha = ifelse(conus, 0.2, 0.15) # defines alpha and 1-alpha quantiles over which the algorithm will search for a solution
break_fitstat = "dev" # using deviance for both gaussian and poisson
max_restart = 2 # how many times you want the algorithm to restart after failing to reach convergence in max_it iterations
max_it = 5 # how many iterations before a restart
glm_it = 15 # maximum number of iterations for feglm to try (note 15 is the same as total number that can occur in fixest_segmented)
smart_init = T # use smart search bounds to initialize segmented algorithm
n_boot = 1000 # number of bootstraps
fitstat_dir = ifelse(break_fitstat %in% c("r2"), -1, 1) # maximize for r2, minimize for deviance
breaks_by = 0.5 # half-year increments
use_fun = feglm # what fixest function to use
interval_width = 2 # how big do you want the smart initialization interval to be? (note: half-year units means 2 = 1 year on each side of the smart search breakpoing estimate)
setFixest_notes(F) # turn off notes if you want to keep console clean but may throw errors in full bootstrap loop 
setFixest_estimation(lean = T, warn = T) # leave warnings on for trycatch code


# build sample information
if(conus){
  sample_obs <- c(rep(50, 3), 100, rep(50, 3))
  sample_years <- c(5, 10, 15, 15, rep(15, 3))
  sample_drop <- c(rep(0, 4), 2020:2022)
  sample_names <- c("n50_y5", "n50_y10", "n50_y15", 
                    "n100_y15", "drop2020", "drop2021", "drop2022")
}
if(!conus){
  sample_names <- c("southeast", "south", "southwest", "west", "northeast",  
                    "east_north_central", "northwest", "central", "west_north_central")
  sample_obs <- rep(50, length(sample_names))
  sample_years <- rep(15, length(sample_names))
  sample_drop <- rep(0, length(sample_names))
  
}

file_names <- file.path(path_dropbox, "output", 
                        paste0(ifelse(conus, "conus/conus_", "regional/regional_"), 
                               pm_type, "_", n_boot, "b", "_FErm", fixef_rm, 
                               "_alpha", break_alpha*100, "_tol", tol, "_corner", corner_push*100,  
                               "_fitstat", break_fitstat, ifelse(smart_init, "_smartinit_", "_randominit_"), 
                               "_maxit", max_it, "_maxrestart", max_restart, "_", sample_names, ".feather"))

comb_name <- file.path(path_dropbox, "output", 
                       paste0(ifelse(conus, "conus/conus_", "regional/regional_"), 
                              pm_type, "_", n_boot, "b", "_FErm", fixef_rm, 
                              "_alpha", break_alpha*100, "_tol", tol, "_corner", corner_push*100,  
                              "_fitstat", break_fitstat, ifelse(smart_init, "_smartinit_", "_randominit_"), 
                              "_maxit", max_it, "_maxrestart", max_restart, "_COMBINED.feather"))
