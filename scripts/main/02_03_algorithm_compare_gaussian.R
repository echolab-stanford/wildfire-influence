# ------------------------------------------------------------------------------
# Written by: Brandon de la Cuesta
# Runs simulations using Gaussian distribution.
# ------------------------------------------------------------------------------
source("scripts/setup/00_01_load_packages.R")
source("scripts/setup/00_02_load_functions.R")
source("scripts/setup/00_03_load_settings.R")

parallel = T
n_compare <- 200
glm_it <- 15

# Condition check to set other globals/global fixest options
# If parallelizing, prevent fixest from using more than 1 thread at a time
if (parallel) {setFixest_nthreads(1)}
if (!parallel) {setFixest_nthreads(n_cores)}
# If not parallelizing, use cores_prop% of cores
cores_prop <- 0.7; n_cores <- round(cores_prop*availableCores())
closeAllConnections(); registerDoFuture()
plan(multisession, workers = n_cores)
registerDoRNG(1234)

# ------------------------------------------------------------------------------
df_psi <- foreach(i = 1:n_compare, .combine = "rbind", 
                  .errorhandling = "remove") %dorng% { # errorhandling here is for segmented function
  
  # Generate random slopes, and break points
  b1 = runif(1, -1, 1)
  b2 = runif(1, -1, 1)
  psi_true = runif(1, 15, 85)
  smart_offset = runif(1, 2, 5)
  smart_bound = c(psi_true + rnorm(1) - smart_offset, 
                  psi_true + rnorm(1) + smart_offset)
  psi_init = runif(1, 15, 85)
  
  # Use them to simulate data
  test_dt <- data.frame(z = 1:100, 
                        e = rnorm(100)) %>% 
    mutate(value = b1*z + b2*(z - psi_true)*(z > psi_true) + e)
  
  # Fit models with segmented and our function 
  our_fit <- fixest_segmented(dt = test_dt, 
                              tol = 1e-8, 
                              fit_family_fn = feglm,  
                              Z_var = "z", 
                              alpha = 0.05, 
                              psi0 = smart_bound, 
                              dep_var = "value", 
                              return_all = T, 
                              glm_it = glm_it)
  
  seg_fit <- segmented(lm(value ~ z, data = test_dt), 
                       seg.Z = ~z, 
                       npsi = 1, 
                       psi = runif(1, smart_bound[1], smart_bound[2]),
                       control = seg.control(conv.psi = T, 
                                             alpha = 0.05, 
                                             fc = 0.95, 
                                             maxit.glm = glm_it, 
                                             tol = 1e-8))
  
  # Return the psi and epsilon 
  res_psi <- data.frame(psi_true = psi_true, 
                        psi_init = psi_init, 
                        seg_psi =  seg_fit %>% 
                          extract2("psi") %>% 
                          magrittr::extract(2), 
                        seg_epsilon = seg_fit$epsilon, 
                        seg_it = seg_fit$it,
                        our_psi = our_fit %>% extract2("psi"), 
                        out_epsilon = our_fit %>% 
                          extract2("epsilon_list") %>% 
                          last, 
                        our_it = our_fit$iterations, 
                        run = i)
  return(res_psi)
}

df_plot <- df_psi %>% 
  group_by(run) %>% 
  dplyr::select(ends_with("psi"), ends_with("_it"), psi_true) %>%
  pivot_longer(!psi_true, names_pattern = "(.*)_(.*)", 
               names_to = c("test", ".value")) %>% 
  filter(!is.na(test))

saveRDS(df_plot, file.path(path_dropbox, "output", "other", 
                           "algorithm_compare_gaussian.rds"))

df_plot %>% 
  ggplot(aes(x = psi_true, y = psi, color = test, shape = test)) + 
  geom_point() + 
  geom_abline(slope = 1, intercept = 0) + 
  scale_color_manual(name = "", 
                     values = c("red", "blue"), 
                     labels = c("ours", "segmented")) + 
  scale_shape_manual(name = "", 
                     values = c(16, 5), 
                     labels = c("ours", "segmented")) + 
  xlab("true breakpoint") + ylab("estimated breakpoint") + 
  theme_classic() +
  theme(legend.position = c(0.8, 0.18))
