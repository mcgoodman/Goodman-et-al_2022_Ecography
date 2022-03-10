#########################################################################################
##' @title VAST spatio-temporal analysis of predator biomass and stomach contents
##' @details Predator considered: Walleye Pollock bin 1
##' @author Arnaud Gruss, with additions by Maurice Goodman
#########################################################################################

### Required Packages ---------------------------------------------------------

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

packages <- c("here")
for (i in seq_along(packages)) {
  if(!do.call(require, list(package = packages[i]))) {
    do.call(install.packages, list(pkgs = packages[i]))
    do.call(require, list(package = packages[i]))
  }
}

######## Define working directories
proj_dir <- here()
output_dir <- paste0(here("model output", "PESC", "Walleye_Pollock_bin1"), "/")
data_dir <- paste0(here("data", "PESC"), "/")
dir.create(output_dir, recursive = TRUE)

### Model Settings ------------------------------------------------------------

######## Define the CPP version to use
Version = get_latest_version( package = "VAST" )

######## Define the spatial resolution for the model and whether to use a grid or mesh approximation
Method = c( "Grid", "Mesh", "Spherical_mesh" )[2]
grid_size_km = 25
#### Specify the number of "knots"
n_x = 300
#### Specify k-means settings
Kmeans_Config = list( "randomseed" = 1, "nstart" = 100, "iter.max" = 1e3 )    

######## Define whether to include spatial and spatio-temporal variation, the rank of the covariance among species, 
######## whether its autocorrelated, and whether there is overdispersion
FieldConfig = c( "Omega1" = "IID", "Epsilon1" = "IID", "Omega2" = "IID", "Epsilon2" = "IID" )
RhoConfig = c( "Beta1" = 0, "Beta2" = 0, "Epsilon1" = 0, "Epsilon2" = 0 ) 
OverdispersionConfig = c( "Eta1"=  0, "Eta2" = 0 )
ObsModel = c( 2, 1 )

######## Decide on which post-hoc calculations to include in the output
Options =  c( "SD_site_density" = 0, "SD_site_logdensity" = 0, "Calculate_Range" = 0, "Calculate_evenness" = 0, 
	"Calculate_effective_area" = 0, "Calculate_Cov_SE" = 0, 'Calculate_Synchrony' = 0, 'Calculate_Coherence' = 0, 
	'Calculate_proportion' = 0 ) 

######## Define any potential stratification of results
strata.limits <- data.frame( 'STRATA' = "All_areas" )

######## Set the "Region" parameter to "User"
Region = "Eastern_Bering_Sea"

######## Define other settings
Use_REML = FALSE
Aniso = TRUE
Npool = 20

######## Function to bundle objects into named list
######## Issue with THorsonUtilities::bundlelist in Rstudio jobs
make_list <- function(...) {
  x <- list(...)
  names(x) <- sapply(substitute(list(...))[-1], deparse)
  return(x)
}

######## Save all settings for later reference
Record = make_list(
  Version, Method, grid_size_km, n_x, FieldConfig,
  RhoConfig, OverdispersionConfig, ObsModel, Options, 
  Use_REML, Aniso, Npool
)

save( Record, file = file.path( output_dir, "Record.RData" ) )
capture.output( Record, file = paste0( output_dir, "Record.txt" ) )

######## Read in and process data
DF = read.csv( paste0( data_dir, "Data_WalleyePollock_bin1_VAST.csv" ) )
Data_Geostat = data.frame( "Year" = DF[,'Year'], "Lat" = DF[,'Latitude'], "Lon" = DF[,'Longitude'] ,
	"AreaSwept_km2" = DF[,"Area_swept"], "Catch_KG" = DF[,'Biomass'],  "spp" = DF[,'Category'] )
Data_Geostat = na.omit( Data_Geostat )

######## Generate an extrapolation grid
Extrapolation_List = make_extrapolation_info( Region = Region, strata.limits = strata.limits )

######## Generate the information used for conducting spatio-temporal parameter estimation, bundled in list `Spatial_List`
Spatial_List = make_spatial_info( grid_size_km = grid_size_km, n_x = n_x, Method = Method, 
	Lon = Data_Geostat[,'Lon'], Lat = Data_Geostat[,'Lat'], fine_scale = TRUE,
	Extrapolation_List = Extrapolation_List, DirPath = output_dir, Save_Results = TRUE, "knot_method" = "grid" )
Data_Geostat = cbind( Data_Geostat, "knot_i"=Spatial_List$knot_i )
save( Data_Geostat, file = file.path( output_dir,"Data_Geostat.RData" ) )
MyPlotDF = cbind( Extrapolation_List[["Data_Extrap"]][,c( 'Lat', 'Lon' )],
                  'x2i' = Spatial_List$PolygonList$NN_Extrap$nn.idx[,1] )
save( MyPlotDF, file = file.path( output_dir, "Extrapolation_grid_information.RData" ) )

######## Examine the data
tapply( Data_Geostat[,'Catch_KG'], INDEX = list( Data_Geostat[,'Year'], Data_Geostat[,'spp'] ), FUN = mean )

######## To estimate parameters, build a list of data inputs used for parameter estimation
Expansion_cz <- matrix( c( 0, 1, 0, 0 ), nrow = 2, ncol = 2 )
Data_orig = make_data( "CheckForErrors" = FALSE, "Version" = Version, "FieldConfig" = FieldConfig, 
	"OverdispersionConfig" = OverdispersionConfig, "RhoConfig" = RhoConfig, "ObsModel" = ObsModel, 
	"c_i" = as.numeric( Data_Geostat[,'spp'] ) - 1, "b_i" = Data_Geostat[,'Catch_KG'],
	"a_i" = Data_Geostat[,'AreaSwept_km2'], "v_i" = rep( 0, nrow( Data_Geostat ) ), 
	"s_i" = Data_Geostat[,'knot_i'] - 1, "t_i" = Data_Geostat[,'Year'], "a_xl" = Spatial_List$a_xl, 
	"MeshList" = Spatial_List$MeshList, "GridList" = Spatial_List$GridList,"Method" = Spatial_List$Method, 
	"Options" = Options, "Aniso" = Aniso, "Expansion_cz" = Expansion_cz, spatial_list = Spatial_List )

### Compile and Optimize Model ------------------------------------------------

######## Build the TMB object 
TmbList_orig = make_model( "build_model" = TRUE, "TmbData" = Data_orig, "RunDir" = output_dir, "Version" = Version, 
	"RhoConfig" = RhoConfig, "loc_x" = Spatial_List$loc_x, "Method" = Method, "Use_REML" = Use_REML, "Npool" = Npool )
Obj_orig = TmbList_orig[["Obj"]]

######## Conduct simple bug checks
#### Check that no data is missing a likelihood component
Report = Obj_orig$report()
Which = which( Report$LogProb1_i + Report$LogProb2_i == 0 )
if( length(Which)>0 ) stop( "Something is wrong" )

#### Check that no fixed effect has a zero gradient
Gr = Obj_orig$gr( Obj_orig$par )
if( any( Gr == 0 ) ) stop( "Something is wrong" )

######## Use a gradient-based nonlinear minimizer to identify maximum likelihood estimates for fixed effects
Opt_orig = TMBhelper::fit_tmb( obj = Obj_orig, lower = TmbList_orig[["Lower"]], upper = TmbList_orig[["Upper"]], 
	getsd = TRUE, savedir = output_dir, bias.correct = FALSE, newtonsteps = 1, 
	getReportCovariance = TRUE, getJointPrecision = TRUE, 
	bias.correct.control = list( sd = FALSE, split = NULL, nsplit = 1, vars_to_correct = "Index_cyl" ) )

### Save Results and Plot Output ----------------------------------------------

######## Size of extrapolation grid cells
extrap_size <- unique(
  Extrapolation_List[["Data_Extrap"]]$Area_in_survey_km2[Extrapolation_List[["Data_Extrap"]]$Area_in_survey_km2 > 0]
)
if(length(extrap_size) > 1) stop("problem with extrapolation grid")

######## Visualize the spatial distribution of data
plot_data(
  Extrapolation_List = Extrapolation_List, 
  Spatial_List = Spatial_List, 
  Data_Geostat = Data_Geostat, 
  PlotDir = output_dir
)

######## Years to return fitted data for 
Year_Set = seq( min( Data_Geostat[,'Year'] ), max( Data_Geostat[,'Year'] ) )
Years2Include = Year_Set[which( Year_Set %in% sort( unique( Data_Geostat[,'Year'] ) ) )]

######## Source some helper functions
source(here("src", "VAST_functions.R"))

######## Extract fitted values
Report_orig = Obj_orig$report()
grid_cell_fits <- VAST_fitted(
  report = Report_orig, 
  spatial_list = Spatial_List,
  extrap_list = Extrapolation_List,
  years = Year_Set
)

######## Remove years without diet data
grid_cell_fits <- grid_cell_fits[grid_cell_fits$year %in% Years2Include,]

######## Save results
Save_orig = list(
  Opt = Opt_orig, 
  Report = Report_orig, 
  ParHat = Obj_orig$env$parList(Opt_orig$par), 
  Data = Data_orig, 
  Obj = Obj_orig, 
  grid_cell_fits = grid_cell_fits
)
saveRDS(Save_orig, file = paste0(output_dir, "Save_orig.RData"))


######## Plot biomass catch rate
biomass_map <- map_EBS_grid(
  E_km = grid_cell_fits$E_km, 
  N_km = grid_cell_fits$N_km,
  fill = log(grid_cell_fits$D_gct_1), 
  facet = grid_cell_fits$year, 
  grid_size = extrap_size
) + labs(fill = "Ln biomass catch rate (kg/ha)", 
         color = "Ln biomass catch rate (kg/ha)")

ggsave(
  filename = paste0(output_dir, "Ln-biomass-catch-rate.png"), 
  plot = biomass_map, width = 10, height = 10, units = "in", dpi = 500
)


######## Plot prey-biomass-per-predator-biomass
predation_map <- map_EBS_grid(
  E_km = grid_cell_fits$E_km, 
  N_km = grid_cell_fits$N_km,
  fill = log(grid_cell_fits$D_gct_2), 
  facet = grid_cell_fits$year, 
  grid_size = extrap_size
) + labs(fill = "Ln prey-biomass-per-predator-biomass (kg/ha)", 
         color = "Ln prey-biomass-per-predator-biomass (kg/ha)")

ggsave(
  filename = paste0(output_dir, "Ln-prey-biomass-per-predator-biomass.png"), 
  plot = predation_map, width = 10, height = 10, units = "in", dpi = 500
)

######## Plot PESC (predator-expanded stomach contents)
predation_map <- map_EBS_grid(
  E_km = grid_cell_fits$E_km, 
  N_km = grid_cell_fits$N_km,
  fill = log(grid_cell_fits$D_gct_1 * grid_cell_fits$D_gct_2), 
  facet = grid_cell_fits$year, 
  grid_size = extrap_size
) + labs(fill = "Ln PESC (kg/ha)",
         color = "Ln PESC (kg/ha)")

ggsave(
  filename = paste0(output_dir, "Ln-PESC.png"), 
  plot = predation_map, width = 10, height = 10, units = "in", dpi = 500
)

######## Plot PWSC (predator-weighted stomach contents)

grid_cell_fits <- grid_cell_fits %>%  
  group_by(year) %>% 
  mutate(sum_biomass_kg = sum(D_gct_1)*grid_size_km) %>% 
  ungroup() %>% 
  mutate(PWSC = (D_gct_1 * D_gct_2) / sum_biomass_kg)

predation_map <- map_EBS_grid(
  E_km = grid_cell_fits$E_km, 
  N_km = grid_cell_fits$N_km,
  fill = log(grid_cell_fits$PWSC), 
  facet = grid_cell_fits$year, 
  grid_size = extrap_size
) + labs(fill = expression("Ln PWSC km"^{-2}), 
         color = expression("Ln PWSC km"^{-2}))

ggsave(
  filename = paste0(output_dir, "Ln PWSC.png"), 
  plot = predation_map, width = 10, height = 10, units = "in", dpi = 500
)

### Sample from PESC predictive distribution ----------------------------------

## grid-scale samples take up a lot of memory - process in batches
## from first run - about 2 GB RAM per 100 samples
## total samples = n_batches * batch_n
n_batches <- 10 ## number of batches of predictive distribution samples
batch_n <- 100 ## number of samples per batch

## Initialize empty list to store PESC predictions
PESC_samples <- vector("list", n_batches*batch_n)

for (i in seq(n_batches)) {
  
  print(paste0("batch ", i, " (samples ", (i-1)*batch_n + 1, " - ", (i-1)*batch_n + batch_n, ")"))
  
  ## Sample from predictive distribution at grid-scale resolution
  ## Dimensions of output: [grid cells, categories, years, sample number]  
  samples <- sample_variable(
    Sdreport = Opt_orig$SD,
    Obj = Obj_orig, 
    variable_name = "D_gct",
    n_samples = batch_n, 
    seed = sample(1:1000, 1) # otherwise, sample_variable will use same seed every time
  )
  
  ## Calculate PESC from samples
  for (j in seq(batch_n)) {
    sample_number <- (i -1)*batch_n + j
    
    ## convert from kg to metric tonnes: 1 metric tonne = 1000 kg
    biomass_grid_scale <- (samples[,1,,j] * extrap_size / 1000) # tonnes
    PESC_grid_scale <- (samples[,1,,j] * samples[,2,,j] * extrap_size / 1000) # tonnes
    
    PESC_samples[[sample_number]] <- data.frame(
      sample = sample_number, 
      year =  Year_Set, 
      biomass_t = colSums(biomass_grid_scale), # sum grid-scales for each year
      PESC_t = colSums(PESC_grid_scale) # sum grid-scales for each year
    )
    
    PESC_samples[[sample_number]]$PWSC <- with(PESC_samples[[sample_number]], PESC_t/biomass_t)
    
    PESC_samples[[sample_number]] <- 
      PESC_samples[[sample_number]][PESC_samples[[sample_number]]$year %in% Years2Include,]
  }
}

## Bind samples into data frame
PESC_samples <- do.call("rbind", PESC_samples)

## Save predictive samples
write.csv(PESC_samples, paste0(output_dir, "PESC_predictive_samples.csv"), row.names = FALSE)

## PESC, PWSC, and biomass maximum likelihood estimates
PESC_MLE <- grid_cell_fits %>% 
  mutate(PESC = D_gct_1 * D_gct_2 * extrap_size / 1000) %>% 
  group_by(year) %>% 
  summarize(PESC_t = sum(PESC), 
            biomass_t = sum(D_gct_1 * extrap_size / 1000), 
            PWSC = PESC_t / biomass_t) %>% 
  pivot_longer(cols = c(PESC_t, PWSC, biomass_t), names_to = "metric", values_to = "MLE")

## Standard deviation of predictive distribution
PESC_SD <- PESC_samples %>% 
  select(-sample) %>% group_by(year) %>% 
  summarize_all(sd) %>% 
  pivot_longer(cols = c(PESC_t, PWSC, biomass_t), names_to = "metric", values_to = "SD")

## Bind data frames together
PESC_summary <- PESC_MLE %>% left_join(PESC_SD, by = c("year", "metric"))
PESC_summary$metric <- gsub("_t", " (tonnes)", PESC_summary$metric)

## Save maximum likelihood estimates and SD of predictive distribution
write.csv(PESC_summary, paste0(output_dir, "PESC_MLE_SD.csv"), row.names = FALSE)

## Plot PESC predictions over time
PESC_plot <- PESC_summary %>% 
  ggplot(aes(year, MLE)) +
  geom_line() + 
  geom_point() + 
  geom_errorbar(aes(ymin = MLE - SD, ymax = MLE + SD), width = 0.4) + 
  facet_wrap(~metric, ncol = 1, scales = "free_y", strip.position = "left") + 
  theme_bw() + 
  scale_x_continuous(breaks = seq(1985, 2015, 5), 
                     labels = function(x) str_sub(as.character(x), 3, 4)) + 
  theme(strip.placement = "outside",
        strip.background = element_blank(), 
        axis.title.y = element_blank())

## Save PESC predictive distribution plot
ggsave(
  paste0(output_dir, "time_series.png"), PESC_plot,
  height = 6, width = 9, units = "in", dpi = 500
)
