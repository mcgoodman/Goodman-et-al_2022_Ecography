
###############################################################################
##' @title VAST Delta-models for Eastern Bering Sea Groundfish
##' @author Maurice Goodman
###############################################################################

### Load and install packages -------------------------------------------------

## VAST package to run models
if(!require("VAST")) {
  if(!require("devtools")) {
    install.packages("devtools")
    require(devtools)
  }
  install_github("james-thorson/VAST", INSTALL_opts="--no-staged-install")
  require(FishStatsUtils)
  require(TMB)
  require(VAST)
}

## Other packages
packages <- c("here", "tidyverse")
for (i in seq_along(packages)) {
  if(!do.call(require, list(package = packages[i]))) {
    do.call(install.packages, list(pkgs = packages[i]))
    do.call(require, list(package = packages[i]))
  }
}

### Read-in biomass data ------------------------------------------------------

## List biomass files
biomass_files <- list.files(here("data", "biomass-catch-rate"), full.names = TRUE)

## Extract species names from file names
species_names <- vapply(
  X = biomass_files, 
  FUN = function(x) {
    gsub(".csv", "", strsplit(basename(x), "rates_")[[1]][2], fixed = TRUE)
  }, 
  FUN.VALUE = character(1), 
  USE.NAMES = FALSE
)

## Read in biomass files
biomass <- lapply(biomass_files, read.csv)
names(biomass) <- species_names

## Directories to save output to
save_dirs <- here("model output", "delta_VAST_output", species_names)
lapply(save_dirs, dir.create, recursive = TRUE)

## Specify data and model structure -------------------------------------------

## format data files for VAST
for (i in seq(length(biomass))) {
  biomass[[i]] <- data.frame(
    Catch_KG = biomass[[i]]$Biomass_catch_rate,
    Year = biomass[[i]]$Year, 
    Lat = biomass[[i]]$Latitude, 
    Lon = biomass[[i]]$Longitude, 
    AreaSwept_km2 = 1, # data is already standardized by area swept
    Vessel = "missing", 
    Pass = 0
  )
}

## Eastern Bering Sea extrapolation grid
EBS_grid <- FishStatsUtils::make_extrapolation_info(Region = "eastern_bering_sea")

## Size of extrapolation grid cells
extrap_size <- EBS_grid$Data_Extrap$Area_in_survey_km2
extrap_size <- unique(extrap_size[extrap_size > 0])

## List to hold knots for spatiotemporal model
EBS_knots <- vector("list", length(biomass))
names(EBS_knots) <- species_names
grid_size_km <- 25 # size of grid cells for determining spatial mesh

## Define knots for spatiotemporal model
for (i in seq_along(biomass)) {
  EBS_knots[[i]] <- FishStatsUtils::make_spatial_info(
    n_x = 300, # Number of knots from Gruss 2020 
    Lon_i = biomass[[i]]$Lon, 
    Lat_i = biomass[[i]]$Lat,
    Extrapolation_List = EBS_grid, 
    Method = "Mesh", # geometric anisotropy
    grid_size_km = grid_size_km, 
    knot_method = "grid", # determine knots based on grid cell locations
    Save_Results = TRUE, 
    DirPath = save_dirs[i],
    fine_scale = TRUE
  )
}

## List to hold model data objects
model_data <- vector("list", length = length(biomass))
names(model_data) <- species_names

## Specify model data and structure
for (i in seq(length(biomass))) {
  model_data[[i]] <- with(biomass[[i]], VAST::make_data(
    b_i = Catch_KG, 
    a_i = AreaSwept_km2, # biomass is already standardized, 
    t_i = Year,
    ObsModel_ez = c(
      "PosDist" = 2, # positive catch rate ~ Gamma
      "Link" = 0 # encounter probability has logit link
    ),
    FieldConfig = c(
      "Omega1" = 1, # spatial variation in encounter probability
      "Epsilon1" = 1, # spatiotemporal variation in encounter probability
      "Omega2" = 1, # spatial variation in positive catch rate
      "Epsilon2" =1 # spatiotemporal variation in positive catch rate
    ), 
    spatial_list = EBS_knots[[i]], 
    Aniso = 1, # Geometric anisotropy,
    RhoConfig = c(
      "Beta1" = 0, # Fixed temporal effects in encounter probability
      "Beta2" = 0, # Fixed temporal effects in positive catch rate
      "Epsilon1" = 0, # encounter probability, iid MVN
      "Epsilon2" = 0 # catch rate, iid MVN
    )
  ))
}

## Compile and optimize models ------------------------------------------------

biomass_models <- vector("list", length = length(model_data))
names(biomass_models) <- species_names

## Compile Models
for (i in seq_along(biomass_models)) {
  biomass_models[[i]] <- VAST::make_model(
    TmbData = model_data[[i]],
    Version = FishStatsUtils::get_latest_version(package = "VAST")
  )
}

## Empty list to store fitted model objects
biomass_fit <- vector("list", length = length(biomass_models))
names(biomass_fit) <- species_names

## Obtain maximum likelihood estimates via gradient descent
for (i in seq_along(biomass_models)) {
  biomass_fit[[i]] <- TMBhelper::fit_tmb(
    obj = biomass_models[[i]][["Obj"]], 
    lower = biomass_models[[i]][["Lower"]], 
    upper = biomass_models[[i]][["Upper"]], 
    getsd = TRUE, 
    bias.correct = FALSE,
    getReportCovariance = TRUE,
    getJointPrecision = TRUE, 
    getHessian = TRUE, 
    newtonsteps = 10, 
    savedir = save_dirs[i]
  )  
}

## Model reports
biomass_reports <- lapply(biomass_models, function(x) x$Obj$report())

## Extract Fitted values ------------------------------------------------------

## Read in convenience functions for VAST models
source(here("src", "VAST_functions.R"))

## Initialize empty list for fitted grid-scale values
biomass_preds <- vector("list", length(biomass_reports))
names(biomass_preds) <- species_names

## Extract fitted probability of occurrence, positive catch rate, density
for (i in seq_along(biomass_reports)) {
  biomass_preds[[i]] <- VAST_fitted(
    report = biomass_reports[[i]], 
    spatial_list = EBS_knots[[i]], 
    extrap_list = EBS_grid,
    years = unique(biomass[[i]]$Year)
  )
}

## Save fitted model objects --------------------------------------------------

## Bind objects into list
delta_models <- list(
  data = model_data, 
  model = biomass_models, 
  fitted = biomass_fit, 
  report = biomass_reports, 
  grid_cell_fits = biomass_preds
)

## Save object list
saveRDS(
  delta_models, 
  file = here("model output", "delta_VAST_output", "delta_models.rds")
)

## Plots of observed vs. fitted values ----------------------------------------

obs_fit <- vector("list", length(species_names))
names(obs_fit) <- species_names

for (i in seq_along(species_names)) {
  
  obs_fit[[i]] <- data.frame(
    observed = as.numeric(model_data[[i]]$b_i),
    fitted = as.numeric(biomass_reports[[i]]$D_i), 
    fitted_p = as.numeric(biomass_reports[[i]]$R1_i)
  )
  
  resid_dir <- paste(save_dirs[i], "residuals", sep = "/")
  if(!dir.exists(resid_dir)) dir.create(resid_dir)
  
  obs_fit_scatter <- obs_fit[[i]] %>%
    ggplot(aes(observed, fitted)) +
    geom_point() + 
    geom_abline(intercept = 0, slope = 1, color = "red", size = 1) + 
    theme_bw()
  
  ggsave(paste0(resid_dir, "/obs_fit_scatter.png"), obs_fit_scatter)
  
  obs_fit_scatter_log <- obs_fit[[i]] %>%
    ggplot(aes(log(observed + 1), log(fitted + 1))) +
    geom_point() + 
    geom_abline(intercept = 0, slope = 1, color = "red", size = 1) + 
    theme_bw()
  
  ggsave(paste0(resid_dir, "/obs_fit_scatter_log.png"), obs_fit_scatter_log)
  
  obs_fit_hist <- obs_fit[[i]] %>% 
    ggplot() + 
    geom_histogram(aes(observed, fill = "observed"), color = "black",
                   alpha = 0.6, bins = 50, position = ) + 
    geom_histogram(aes(fitted, fill = "fitted"), color = "black",
                   alpha = 0.6, bins = 50, position = ) + 
    theme_bw() + 
    theme(legend.position = c(1, 1), 
          legend.justification = c(1, 1), 
          legend.background = element_rect(color = "black"), 
          legend.title = element_blank()) + 
    coord_cartesian(clip = "off")
  
  ggsave(paste0(resid_dir, "/obs_fit_hist.png"), obs_fit_hist)
  
  obs_fit_hist_log <- obs_fit[[i]] %>%
    ggplot() + 
    geom_histogram(aes(log(observed + 1), fill = "observed"), color = "black",
                   alpha = 0.6, bins = 50, position = ) + 
    geom_histogram(aes(log(fitted + 1), fill = "fitted"), color = "black", 
                   alpha = 0.6, bins = 50, position = ) + 
    theme_bw() + 
    theme(legend.position = c(1, 1), 
          legend.justification = c(1, 1), 
          legend.background = element_rect(color = "black"), 
          legend.title = element_blank()) + 
    coord_cartesian(clip = "off")
  
  ggsave(paste0(resid_dir, "/obs_fit_hist_log.png"), obs_fit_hist_log)
  
}

## Visualize model output -----------------------------------------------------

## Plot spatial distribution of data
for (i in seq_along(biomass_fit)) {
  plot_data(
    Extrapolation_List = EBS_grid, 
    Spatial_List = EBS_knots[[i]], 
    Data_Geostat = biomass[[i]], 
    PlotDir = paste0(save_dirs[i], "/")
  )
}

## Values for fitted values plot labels
vars <- list(
  name = c(
    expression("ln biomass catch rate kg km"^-2), 
    "Probability of Occurrence", 
    expression("ln positive catch rate kg km"^-2)
  ), 
  filename = c("Ln_Biomass_Catch_Rate", "Probability_of_Occurrence", "Ln_Positive_Catch_Rate"),
  var = c("D_gct", "R1_gct", "R2_gct"), 
  log = c(TRUE, FALSE, TRUE)
)

## Plot fitted values at grid-scale
for (i in seq_along(biomass_preds)) {
  for (j in seq_along(vars$var)) {
    
    print(paste("plotting:", species_names[i], vars$name[j]))
    
    ## Log catch rate data (and convert to kg/km2)
    fill_var <- biomass_preds[[i]][[vars$var[j]]]
    if(vars$log[j]) fill_var <- log(fill_var * 100)
    
    fit_map <- map_EBS_grid(
      E_km = biomass_preds[[i]]$E_km, 
      N_km = biomass_preds[[i]]$N_km, 
      fill = fill_var, 
      facet = biomass_preds[[i]]$year, 
      grid_size = extrap_size
    ) + labs(fill = vars$name[j], 
             color = vars$name[j])
    
    ggsave(
      filename = paste0(save_dirs[i], "/", vars$filename[j], ".png"), 
      plot = fit_map, width = 10, height = 10, units = "in", dpi = 500
    )
    
  }
}

### Overlap metrics at MLE ----------------------------------------------------

## Read in predator-prey overlap functions
source(here("src", "range_overlap_functions.R"))

years <- unique(biomass$Juvenile_Walleye_pollock$Year)

## Names of predators to calculate overlap for
predators <- c("Flounder", "Pacific_cod", "Pacific_halibut", "Walleye_pollock")

## Coordinates of grid cells
grid_x <- as.data.frame(EBS_knots$Flounder$loc_g)$E_km
grid_y <- as.data.frame(EBS_knots$Flounder$loc_g)$N_km
n_cells <- length(grid_x)

### Prediction thresholds for area-overlap
fit_thresholds <- rep(NA, length(species_names))
names(fit_thresholds) <- species_names

### Plot ROC curves for each species
### and identify threshold that maximizes true positive + true negative rate
for (i in seq_along(fit_thresholds)) {
  p <- obs_fit[[i]]$fitted_p[obs_fit[[i]]$observed > 0]
  a <- obs_fit[[i]]$fitted_p[obs_fit[[i]]$observed == 0]
  fit_eval <- dismo::evaluate(p, a)
  png(paste(save_dirs[i], "ROC.png", sep = "/")); plot(fit_eval, "ROC"); dev.off()
  fit_thresholds[i] <- dismo::threshold(fit_eval)$spec_sens
}

## Initialize empty vector
overlap_mle <- vector("list", length = length(predators)*length(years))
index <- 1 ## to keep track of position in list

## Calculate overlap metrics at MLE
for (i in seq_along(predators)) {
  for (j in seq_along(years)) {
    overlap_mle[[index]] <- data.frame(
      species = predators[i], 
      year = years[j],
      loc_colloc = loc_collocfn( # Local index of collocation
        prey = biomass_reports$Juvenile_Walleye_pollock$D_gct[,1,j], 
        pred = biomass_reports[[predators[i]]]$D_gct[,1,j]
      ),
      area_overlap = area_overlapfn( # Area overlap
        prey = biomass_reports$Juvenile_Walleye_pollock$R1_gct[,1,j] > fit_thresholds["Juvenile_Walleye_pollock"],
        pred = biomass_reports[[predators[i]]]$R1_gct[,1,j] > fit_thresholds[predators[i]],
        area = rep(extrap_size, n_cells)
      ), 
      global_colloc = glob_collocfn( # Global Index of collocation
        prey_x = grid_x,
        prey_y = grid_y, 
        prey = biomass_reports$Juvenile_Walleye_pollock$D_gct[,1,j], 
        pred_x = grid_x, 
        pred_y = grid_y, 
        pred = biomass_reports[[predators[i]]]$D_gct[,1,j]
      )
    )
    index <- index + 1
  }
}

## bind MLE results into data frame
overlap_mle <- do.call("rbind", overlap_mle)

write.csv(
  rownames_to_column(data.frame(threshold = fit_thresholds), "species"),
  here("model output", "delta_VAST_output", "fit_thresholds.csv"), 
  row.names = FALSE
)

## Total biomass at MLE -------------------------------------------------------

biomass_mle <- vector("list", length(biomass_fit))
names(biomass_mle) <- species_names

## divide by 10:
## est (kg/ha) * area conv. (ha/km^2) * mass conv. (tonnes/kg) * area (km^2) 
for (i in seq_along(biomass_fit)) {
  biomass_mle[[species_names[i]]] <- data.frame(
    species = species_names[i], 
    year = years,
    biomass = colSums(
      biomass_reports[[species_names[i]]]$D_gct[,1,] * extrap_size / 10
    )
  )
}

biomass_mle <- do.call("rbind", biomass_mle)


## Predictive distribution of biomass and overlap -----------------------------

## grid-scale samples take up a lot of memory - process in batches
## from first run - about 2 GB RAM per 100 samples
## total samples = n_batches * batch_n
n_batches <- 10 ## number of batches of predictive distribution samples
batch_n <- 100 ## number of samples per batch
n_samples <- n_batches * batch_n ## total number of samples

## Initialize empty list to store biomass predictions
biomass_dist <- vector("list", length(biomass_fit)*n_samples)
b_index <- 1

## Initialize empty vector to store overlap predictions
overlap <- vector("list", length = length(predators)*length(years)*n_samples)
o_index <- 1

for (i in seq(n_batches)) {
  
  print(paste0("batch ", i, " (samples ", (i-1)*batch_n + 1, " - ", (i-1)*batch_n + batch_n, ")"))
  
  ## Initialize empty list to hold predictive distribution samples
  fit_samples <- vector("list", length(biomass_fit))
  names(fit_samples) <- species_names
  
  ## Sample from predictive distributions of occurrence and biomass
  for(j in seq_along(biomass_fit)) {
    
    print(paste("sampling:", species_names[j], "occurrence"))
    
    fit_samples[[j]]$occurrence <- sample_variable(
      Sdreport = biomass_fit[[j]]$SD, 
      Obj = biomass_models[[j]]$Obj, 
      variable_name = "R1_gct", 
      n_samples = batch_n,
      seed = sample(1:1000, 1) # otherwise, sample_variable will use same seed every time
    )
    
    print(paste("sampling:", species_names[j], "biomass"))
    
    fit_samples[[j]]$biomass <- sample_variable(
      Sdreport = biomass_fit[[j]]$SD, 
      Obj = biomass_models[[j]]$Obj, 
      variable_name = "D_gct", 
      n_samples = batch_n,
      seed = sample(1:1000, 1) # otherwise, sample_variable will use same seed every time
    )
    
    ## Compute biomass for predictive distribution samples
    for (k in 1:batch_n) {
      biomass_dist[[b_index]] <- data.frame(
        species = species_names[j], 
        sample = (i - 1)*batch_n + k, 
        year = years,
        biomass = colSums(
          fit_samples[[species_names[j]]]$biomass[,1,,k] * extrap_size / 10
        )
      )
      b_index <- b_index + 1
    }
  }
  
  ## Calculate overlap metrics for samples from predictive distribution
  for (j in seq_along(predators)) {
    for (k in 1:batch_n) {
      for (l in seq_along(years)) {
        overlap[[o_index]] <- data.frame(
          species = predators[j], 
          sample = (i - 1)*batch_n + k, 
          year = years[l],
          loc_colloc = loc_collocfn( # Local index of collocation
            prey = fit_samples$Juvenile_Walleye_pollock$biomass[,1,l,k], 
            pred = fit_samples[[predators[j]]]$biomass[,1,l,k]
          ), 
          area_overlap = area_overlapfn( # Area overlap
            prey = fit_samples$Juvenile_Walleye_pollock$occurrence[,1,l,k] > fit_thresholds["Juvenile_Walleye_pollock"],
            pred = fit_samples[[predators[j]]]$occurrence[,1,l,k] > fit_thresholds[predators[j]],
            area = rep(extrap_size, n_cells)
          ), 
          global_colloc = glob_collocfn( # Global Index of collocation
            prey_x = grid_x,
            prey_y = grid_y, 
            prey = fit_samples$Juvenile_Walleye_pollock$biomass[,1,l,k], 
            pred_x = grid_x, 
            pred_y = grid_y, 
            pred = fit_samples[[predators[j]]]$biomass[,1,l,k]
          )
        )
        o_index <- o_index + 1
      }
    }
  }
}

## bind results into data frames
biomass_dist <- do.call("rbind", biomass_dist)
overlap <- do.call("rbind", overlap)

## Save output
dir.create(here("model output", "overlap_metrics"), recursive = TRUE)
saveRDS(
  list(MLE = overlap_mle, samples = overlap),
  file = here("model output", "overlap_metrics", "overlap_metrics.rds")
)
saveRDS(
  list(MLE = biomass_mle, samples = biomass_dist), 
  here("model output", "delta_VAST_output", "biomass.rds")
)

### Plot time series of biomass -----------------------------------------------

biomass_ts <- biomass_dist %>% 
  group_by(species, year) %>% 
  summarize(sd = sd(biomass)/1e+06) %>% 
  left_join(biomass_mle) %>% 
  mutate(biomass = biomass/1e+06, 
         species = gsub("_", " ", species)) %>% 
  ggplot(aes(year, biomass)) +
  geom_errorbar(aes(ymin = biomass - sd, ymax = biomass + sd), 
                width = 0.4) +
  geom_line() + 
  geom_point() +
  facet_wrap(~species, ncol = 1, scales = "free") + 
  scale_x_continuous(breaks = seq(1985, 2015, 5), 
                     labels = function(x) str_sub(as.character(x), 3, 4)) + 
  theme_bw() + 
  theme(strip.background = element_blank()) + 
  labs(y = "biomass (Mt)")

ggsave(here("model output", "delta_VAST_output", "biomass.png"), biomass_ts,
       height = 8, width = 5, units = "in", dpi = 500)

### Plot time series of overlap metrics ---------------------------------------

## One-time use function to plot overlap metrics
overlap_plot <- function(metric, ylab, samples = overlap, mle = overlap_mle) {
  
  metric <- enquo(metric)
  
  samples %>% 
    mutate(species = gsub("_", " ", species)) %>% 
    ggplot(aes(year, !!metric, group = sample)) +
    geom_line(color = "grey60", alpha = 0.2) + 
    facet_wrap(~species, nrow = 1) + 
    geom_line(aes(year, !!metric), color = "black", 
              inherit.aes = FALSE, size = 1, 
              data = mutate(mle, species = gsub("_", " ", species))) + 
    theme_bw() + 
    theme(strip.background = element_blank()) + 
    scale_x_continuous(breaks = seq(1985, 2015, 5), 
                       labels = function(x) str_sub(as.character(x), 3, 4)) + 
    labs(y = ylab)
  
}

## Area Overlap 
ggsave(
  here("model output", "overlap_metrics", "area_overlap.png"),
  plot = overlap_plot(area_overlap, "Area Overlap"),
  height = 4, width = 10, units = "in", dpi = 500
)

## Local Index of Collocation
ggsave(
  here("model output", "overlap_metrics", "loc_colloc.png"),
  plot = overlap_plot(loc_colloc, "Local Index of Collocation"),
  height = 4, width = 10, units = "in", dpi = 500
)

## Global Index of Collocation
ggsave(
  here("model output", "overlap_metrics", "global_colloc.png"),
  plot = overlap_plot(global_colloc, "Global Index of Collocation"),
  height = 4, width = 10, units = "in", dpi = 500
)

### Range Centroids -----------------------------------------------------------

## Initialize empty vector
range_centroids <- vector("list", length(biomass_preds))
names(range_centroids) <- species_names

## Calculate range centroids for each year
for(i in seq_along(range_centroids)) {
  range_centroids[[i]] <- biomass_preds[[i]] %>% 
    select(E_km, N_km, year, D_gct, R1_gct) %>% 
    pivot_longer(cols = c(D_gct, R1_gct), names_to = "metric") %>% 
    group_by(year, metric) %>% 
    summarize(E_km = weighted.mean(E_km, value), 
              N_km = weighted.mean(N_km, value))
}

## bind results into data frame
range_centroids <- bind_rows(range_centroids, .id = "species")

## Plot northings over time
northing_plot <- range_centroids %>% 
  mutate(metric = ifelse(metric == "D_gct", "biomass", "occurrence"), 
         species = gsub("_", " ", species)) %>% 
  ggplot(aes(year, N_km, color = metric)) + 
  geom_line() + 
  geom_point() + 
  facet_wrap(~species, nrow = 1, scales = "free_y") + 
  theme_bw() + 
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        legend.position = "bottom", 
        axis.ticks = element_line(color = "black"), 
        axis.text = element_text(color = "black"), 
        panel.border = element_rect(color = "black")) + 
  scale_x_continuous(breaks = seq(1985, 2015, 5), 
                     labels = function(x) str_sub(as.character(x), 3, 4)) + 
  scale_color_manual(values = c("#50A3A4", "#F95335")) + 
  labs(y = "centroid northings (km)", color = "weighting") + 
  coord_cartesian(clip = "off")

## Save range centroid plot
ggsave(
  here("model output", "delta_VAST_output", "centroid_northings.png"),
  plot = northing_plot, height = 3, width = 10, units = "in", dpi = 500, 
  scale = 0.95
)

## Distance between predator centroids and juvenile pollock centroids
centroid_dists <- vector("list", length(predators))
prey_xy <- range_centroids[range_centroids$species == "Juvenile_Walleye_pollock", c("metric", "E_km", "N_km")]
prey_xy <- split(prey_xy[,c("E_km", "N_km")], prey_xy$metric)

dist2d <- function(p1, p2) {
  as.numeric(sqrt((p1[,1] - p2[,1])^2 + (p1[,2] - p2[,2])^2)[,1])
}

for(i in seq_along(centroid_dists)) {
  pred_xy <- range_centroids[range_centroids$species == predators[i], c("metric", "E_km", "N_km")]
  pred_xy <- split(pred_xy[,c("E_km", "N_km")], pred_xy$metric)
  centroid_dists[[i]] <- data.frame(
    predator = predators[i], 
    year = years, 
    dist_biomass = dist2d(prey_xy$D_gct, pred_xy$D_gct), 
    dist_occurrence = dist2d(prey_xy$R1_gct, pred_xy$D_gct)
  )
}

centroid_dists <- do.call("rbind", centroid_dists)

## Plot distance between juvenile pollock and predator centroids over time
dist_plot <- centroid_dists %>% 
  mutate(predator = gsub("_", " ", predator)) %>% 
  ggplot(aes(year, dist_biomass)) + 
  geom_line() + 
  geom_point() + 
  facet_wrap(~predator, nrow = 1) + 
  theme_bw() + 
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 12, face = "bold"),
        legend.position = "bottom", 
        axis.ticks = element_line(color = "black"), 
        axis.text = element_text(color = "black"), 
        panel.border = element_rect(color = "black")) + 
  scale_x_continuous(breaks = seq(1985, 2015, 5), 
                     labels = function(x) str_sub(as.character(x), 3, 4)) + 
  coord_cartesian(clip = "off") + 
  expand_limits(y = 0) + 
  labs(y = "distance (km)")
  
## Save centroid distance plot
ggsave(
  here("model output", "delta_VAST_output", "centroid_distances.png"),
  plot = dist_plot, height = 3, width = 10, units = "in", dpi = 500, 
  scale = 0.95
)

### Read in cold pool data
cold_pool <- read.csv(here("data", "cold_pool_extent_annual.csv")) %>% 
  select(year = YEAR, cold_pool = AREA_SUM_KM2_LTE2) %>% 
  mutate(cold_pool_sd = as.numeric(scale(cold_pool)))

### Regressions of biomass-weighted centroid northings vs cold pool extent
range_centroids %>% 
  filter(metric == "D_gct") %>% 
  left_join(cold_pool, by = "year") %>% 
  group_by(species) %>% 
  mutate(N_km = as.numeric(scale(N_km))) %>% 
  nest() %>% 
  mutate(model = purrr::map(data, ~lm(N_km ~ cold_pool_sd, data = .x)), 
         coef = purrr::map(model, ~as.data.frame(summary(.x)$coefficients)), 
         coef = purrr::map(coef, rownames_to_column, var = "term")) %>% 
  select(species, coef) %>% 
  unnest(cols = coef) %>% 
  mutate_if(is.numeric, round, digits = 3) %>% 
  write_csv(here("model output", "delta_VAST_output", "northings_coldpool.csv"))

## Plot range centroid northings against cold pool extent
cold_pool_northings_plot <- range_centroids %>% 
  filter(metric == "D_gct") %>% 
  mutate(species = gsub("_", " ", species)) %>% 
  left_join(cold_pool, by = "year") %>% 
  ggplot(aes(cold_pool/1e+03, N_km)) + 
  geom_point() + 
  geom_smooth(method = "lm", color = "black") + 
  facet_wrap(~species, scales = "free", nrow = 1) + 
  theme_bw() + 
  theme(strip.background = element_blank(), 
        strip.text = element_text(size = 10), 
        axis.text = element_text(color = "black"), 
        axis.ticks = element_line(color = "black"), 
        panel.border = element_rect(color = "black")) + 
  labs(y = "centroid northings (km)", 
       x = expression(paste("cold pool extent ", "(", "km \U00D7 10"^{3},")"))) + 
  coord_cartesian(clip = "off")

ggsave(
  here("model output", "delta_VAST_output", "northings_cold_pool.png"),
  plot = cold_pool_northings_plot, height = 3, width = 10, units = "in", dpi = 500
)

## Plot range centroid eastings against cold pool extent
cold_pool_eastings_plot <- range_centroids %>% 
  filter(metric == "D_gct") %>% 
  mutate(species = gsub("_", " ", species)) %>% 
  left_join(cold_pool, by = "year") %>% 
  ggplot(aes(cold_pool/1e+03, E_km)) + 
  geom_point() + 
  geom_smooth(method = "lm", color = "black") + 
  facet_wrap(~species, scales = "free", nrow = 1) + 
  theme_bw() + 
  theme(strip.background = element_blank(), 
        strip.text = element_text(size = 10), 
        axis.text = element_text(color = "black"), 
        axis.ticks = element_line(color = "black"), 
        panel.border = element_rect(color = "black")) + 
  labs(y = "centroid eastings (km)", 
       x = expression(paste("cold pool extent ", "(", "km \U00D7 10"^{3},")"))) + 
  coord_cartesian(clip = "off")

ggsave(
  here("model output", "delta_VAST_output", "eastings_cold_pool.png"),
  plot = cold_pool_eastings_plot, height = 3, width = 10, units = "in", dpi = 500
)

