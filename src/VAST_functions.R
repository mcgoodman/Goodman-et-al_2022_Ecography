
###############################################################################
##' @title Functions for handling output from VAST models
##' @author Maurice Goodman
###############################################################################

## Required Packages and data -------------------------------------------------

if(!require("VAST")) {
  if(!require("devtools")) {
    install.packages("devtools")
    require(devtools)
  }
  install_github("james-thorson/VAST", INSTALL_opts="--no-staged-install")
  require(FishStatsUtils)
  require(VAST)
}

packages <- c("tidyverse", "sp", "mapdata")
for (i in seq_along(packages)) {
  if(!do.call(require, list(package = packages[i]))) {
    do.call(install.packages, list(pkgs = packages[i]))
    do.call(require, list(package = packages[i]))
  }
}

## Simplified Alaskan coastline, including islands
alaska_utm <- read.csv(here("data", "GIS", "Alaska Shoreline", "alaska_crop_utm_simplified.csv"))

### Function to return fitted values ------------------------------------------

#' @title VAST_fitted 
#' @description return grid-scale fitted values from VAST in a data frame
#' @param report output of calling report() from compiled VAST model
#' @param spatial_list output from FishStatsUtils::make_spatial_info
#' @param vars name of variables to return (found in report)
#' @param years optional, unique years from model data
#' @return a data frame with each var, coordinates (eastings, northings, lat, lon), year, and grid size
VAST_fitted <- function(report, spatial_list, extrap_list,
                        vars = c("D_gct", "R1_gct", "R2_gct"), 
                        years = NULL) {
  
  if(!spatial_list$fine_scale) stop("results are not at grid-scale")
  
  grid_km <- as.data.frame(spatial_list$loc_g) # Eastings & Northings
  grid_latlon <- as.data.frame(spatial_list$latlon_g) # Latitude and longitude
  n_years <- dim(report[[vars[1]]])[3] # number of years in report output
  if(is.null(years)) years <- 1:n_years
  
  # check that vars are in report object
  for (i in seq_along(vars)) {
    if(is.null(report[[vars[i]]])) {
      stop("variable ", vars[i], " is not in report", sep = " ")
    }
  }
  
  # check dimensions of report and spatial list
  if(nrow(grid_km) != dim(report[[vars[1]]])[1]) {
    stop("dimensions of report and spatial_list do not match")
  }
  
  # check input for years. If years argument is NA, will not throw error
  if(!is.numeric(years) | length(years) != n_years) {
    stop("years must be numeric and match report dimensions")
  }
  
  # Empty list to store data
  df <- vector("list", n_years)
  
  # size of extrapolation grid cells (km)
  grid_size <- extrap_list$Data_Extrap$Area_in_survey_km2
  grid_size <- unique(grid_size[grid_size > 0])
  if(length(grid_size) > 1) {
    stop("something is wrong with extrapolation grid")
  }
  
  for (i in seq(n_years)) {
    # Store coordinates and grid cell size
    df[[i]] <- data.frame(
      E_km = grid_km$E_km, 
      N_km = grid_km$N_km, 
      lat = grid_latlon$Lat, 
      lon = grid_latlon$Lon,
      grid_size_km = grid_size, 
      year = years[i]
    )
    # Extract values for each variable
    for(j in seq_along(vars)) {
      n_cat <- dim(report[[vars[j]]])[2] # number of categories
      if(n_cat > 1) {
        for (k in seq(n_cat)) {
          df[[i]][[paste(vars[j], k, sep = "_")]] <- report[[vars[j]]][,k,i]
        }
      } else {
        df[[i]][[vars[j]]] <- report[[vars[j]]][,,i]
      }
    }
  }
  
  # Collapse list to single data frame
  df <- do.call("rbind", df)
  df
}


### Functions to plot fitted values -------------------------------------------

#' @title theme_map
#' @description theme for ggplot maps
#' @param plot ggplot object
theme_map <- function(plot) {
  
  require(tidyverse)
  
  plot +     
    theme_bw() + 
    guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5, 
                                 frame.colour = "black", frame.linewidth = 1.5, 
                                 ticks.colour = "black", ticks.linewidth = 1.5, 
                                 barwidth = unit(14, "lines"))) +
    guides(color = guide_colorbar(title.position = "top", title.hjust = 0.5, 
                                  frame.colour = "black", frame.linewidth = 1.5, 
                                  ticks.colour = "black", ticks.linewidth = 1.5, 
                                  barwidth = unit(14, "lines"))) + 
    theme(panel.border = element_rect(color = "grey70", size = 1),
          panel.background = element_rect(fill = "grey80"),
          panel.spacing.x = unit(0, "lines"),
          strip.text = element_text(size = 12, face = "bold"), 
          strip.background = element_blank(), 
          axis.text = element_text(size = 12, color = "grey70", face = "bold"),
          legend.title = element_text(size = 14),
          axis.ticks.length = unit(0.4, "lines"),
          axis.ticks = element_line(size = 1, color = "grey70"),
          axis.title = element_blank(),
          panel.grid = element_blank(), 
          plot.title = element_text(size = 16, face = "bold"),
          legend.position = "bottom")
}


#' @title map_EBS_grid
#' @description plot grid-scale fitted values for the Eastern Bering Sea
#' @param E_km vector with location of grid cells, eastings
#' @param N_km vector with location of grid cells, northings
#' @param fill vector with variable to color grid cells by
#' @param facet vector with variable to facet by, e.g. year
#' @param grid_size height & width of grid cells in km
#' @return a map of EBS survey area, with grid cells colored by fill variable
map_EBS_grid <- function(E_km, N_km, fill, facet, grid_size) {
  
  require(tidyverse)
  
  # Bind data into data frame
  plot_data <- data.frame(E_km = E_km, N_km = N_km, fill = fill, facet = facet)
  
  ## polygons for unique locations
  extrap_diag <- sqrt(2*grid_size)
  loc_g <- as.data.frame(unique(cbind(E_km, N_km)))
  poly_data <- purrr::map2_dfr(loc_g$E_km, loc_g$N_km, ~data.frame(
    x = c(.x - extrap_diag/2, .x, .x + extrap_diag/2, .x), 
    y = c(.y, .y - extrap_diag/2, .y, .y + extrap_diag/2)
  ), .id = "group")
  poly_data$E_km <- rep(loc_g$E_km, each = 4)
  poly_data$N_km <- rep(loc_g$N_km, each = 4)
  
  ## merge polygons with fill and facet information
  poly_data <- poly_data %>% left_join(plot_data, by = c("E_km", "N_km"))
  
  ### Plot EBS map
  EBS_map <- poly_data %>% 
    ggplot() + 
    scale_fill_viridis_c(option = "magma") +
    scale_color_viridis_c(option = "magma") +
    facet_wrap(~facet, ncol = 6) + 
    coord_cartesian(xlim = c(0, 1350), ylim = c(6000, 7000), expand = FALSE, clip = "off") +
    scale_x_continuous(breaks = c(0, 500, 1000)) + 
    geom_polygon(aes(x, y, group = group, fill = fill, color = fill), size = 0.0005) + 
    geom_polygon(aes(E_km, N_km, group = group), data = alaska_utm, 
                 fill = "white", color = NA, 
                 inherit.aes = FALSE)
  
  ### Clean up and return map
  theme_map(EBS_map)
  
}
