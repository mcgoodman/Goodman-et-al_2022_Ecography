
###############################################################################
##' @title Example plot - arrowtooth flounder predation in 2009 and 2015
##' @author Maurice Goodman
###############################################################################

packages <- c("here", "tidyverse", "cowplot", "mgcv", "rgdal", "sp")
for (i in seq_along(packages)) {
  if(!do.call(require, list(package = packages[i]))) {
    do.call(install.packages, list(pkgs = packages[i]))
    do.call(require, list(package = packages[i]))
  }
}

### Read in polygons ----------------------------------------------------------

## Simplified Alaskan coastline, including islands
alaska <- read.csv(here("data", "GIS", "Alaska Shoreline", "alaska_crop_utm_simplified.csv"))

## EBS Survey grid, including sub sampling around St. Matthew and Pribilof islands
EBS_survey <- fortify(readOGR(here("data", "GIS", "EBS survey grid")), region = "STATIONID")
EBS_survey <- EBS_survey %>% 
  mutate(
    E_km = long/1000, 
    N_km = lat/1000, 
    STATIONID = id, 
    group = as.numeric(as.factor(group))
  ) %>% 
  select(
    E_km, N_km, STATIONID, group
  )

## Boundary of EBS Survey grid
EBS_boundary <- fortify(readOGR(here("data", "GIS", "EBS Survey Boundary")))
EBS_boundary <- EBS_boundary %>% 
  mutate(
    E_km = long/1000, 
    N_km = lat/1000, 
    group = 1
  ) %>% 
  select(
    E_km, N_km, group
  )

### Functions and theme for plots ---------------------------------------------

## Theme for map panels
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
    theme(panel.border = element_rect(color = "grey80", size = 1),
          panel.background = element_rect(fill = "grey80"),
          panel.spacing.x = unit(0, "lines"),
          strip.text = element_text(size = 12, face = "bold"), 
          strip.background = element_blank(), 
          axis.text = element_text(size = 12, color = "grey70", face = "bold"),
          legend.title = element_text(size = 14),
          axis.ticks.length = unit(0.4, "lines"),
          axis.ticks = element_line(size = 1, color = "grey80"),
          axis.title = element_blank(),
          panel.grid = element_blank(), 
          plot.title = element_text(size = 16, face = "bold"),
          legend.position = "bottom")
}

## Function to plot data at extrapolation grid cells
## Grid cells are aligned in a diamond mesh - necessitates creating polygons for each
map_EBS_extrap <- function(E_km, N_km, fill, grid_size) {
  
  require(tidyverse)
  
  # Bind data into data frame
  plot_data <- data.frame(E_km = E_km, N_km = N_km, fill = fill)
  
  ## polygons for unique locations
  extrap_diag <- sqrt(2*grid_size)
  loc_g <- as.data.frame(unique(cbind(E_km, N_km)))
  poly_data <- purrr::map2_dfr(loc_g$E_km, loc_g$N_km, ~data.frame(
    x = c(.x - extrap_diag/2, .x, .x + extrap_diag/2, .x), 
    y = c(.y, .y - extrap_diag/2, .y, .y + extrap_diag/2)
  ), .id = "group")
  poly_data$E_km <- rep(loc_g$E_km, each = 4)
  poly_data$N_km <- rep(loc_g$N_km, each = 4)
  
  ## merge polygons with fill and information
  poly_data <- poly_data %>% left_join(plot_data, by = c("E_km", "N_km"))
  
  ### Plot EBS map
  EBS_map <- poly_data %>% 
    ggplot() +
    coord_cartesian(xlim = c(0, 1350), ylim = c(6000, 7000), expand = FALSE, clip = "off") +
    scale_x_continuous(breaks = c(0, 500, 1000)) + 
    geom_polygon(aes(x, y, group = group, fill = fill, color = fill), size = 0.0005) + 
    geom_polygon(aes(E_km, N_km, group = group), data = alaska, 
                 fill = "white", color = NA, 
                 inherit.aes = FALSE)
  
}

## Function to plot map for data taken at individual stations
map_EBS_survey <- function(STATIONID, color, limits, palette = "magma") {
  
  require(tidyverse)
  
  # Bind data into data frame
  plot_data <- data.frame(STATIONID = STATIONID, color = color)
  plot_data <- plot_data %>% left_join(EBS_survey, by = "STATIONID")
  plot_data <- na.omit(plot_data)
  
  ### Plot EBS map
  EBS_map <- plot_data %>% 
    ggplot(aes(E_km, N_km, group = STATIONID, subgroup = group)) +
    geom_polygon(aes(group = group), data = EBS_boundary, 
                 fill = "grey60", color = NA) + 
    geom_polygon(aes(fill = color, color = color), size = 0.01) +
    scale_fill_viridis_c(option = "magma", na.value = "grey70") +
    scale_color_viridis_c(option = "magma", na.value = "grey70") +
    geom_polygon(aes(E_km, N_km, group = group), data = alaska,
                 fill = "white", color = NA,
                 inherit.aes = FALSE) +
    coord_cartesian(xlim = c(0, 1350), ylim = c(6000, 7000), expand = FALSE, clip = "off") + 
    scale_x_continuous(breaks = c(0, 500, 1000)) + 
    scale_y_continuous(breaks = seq(6000, 7000, 250))
  
  ### Clean up and return map
  theme_map(EBS_map)
  
}

### Read in and process data --------------------------------------------------

## read in delta model fits
delta_models <- readRDS(here("model output", "delta_VAST_output", "delta_models.rds"))

## read in bottom temperature data
bottom_temp <- read.csv(here("data", "bottom_temp_data.csv"))
bottom_temp <- bottom_temp[bottom_temp$YEAR == 2009 | bottom_temp$YEAR == 2015,]
bottom_temp_coords <- FishStatsUtils::project_coordinates(
  X = bottom_temp$START_LONGITUDE, 
  Y = bottom_temp$START_LATITUDE, 
  projargs = "+proj=utm +zone=2 +datum=WGS84"
)
bottom_temp$E_km <- bottom_temp_coords[,1]/1000
bottom_temp$N_km <- bottom_temp_coords[,2]/1000


## GAMS for 2009 and 2015 bottom temperature, to use for contours
btemp_gam <- list(
  "2009" = gam(
    GEAR_TEMPERATURE ~ s(E_km, N_km, k = 150),
    data = bottom_temp[bottom_temp$YEAR == 2009,]
  ), 
  "2015" = gam(
    GEAR_TEMPERATURE ~ s(E_km, N_km, k = 150),
    data = bottom_temp[bottom_temp$YEAR == 2015,]
  )
)

## Predict on 5km EBS grid
EBS_grid <- read.csv(here("data", "GIS", "EBS_grid_5km_UTM.csv"))
EBS_grid$temp_2009 <- predict(btemp_gam[["2009"]], newdata = EBS_grid)
EBS_grid$temp_2015 <- predict(btemp_gam[["2015"]], newdata = EBS_grid)

## read in PESC grid cell fits for arrowtooth flounder
atf_bin1 <- readRDS(here("model output", "PESC", "Arrowtooth_Flounder_bin1", "Save_orig.RData"))$grid_cell_fits
atf_bin2 <- readRDS(here("model output", "PESC", "Arrowtooth_Flounder_bin2", "Save_orig.RData"))$grid_cell_fits
extrap_size <- unique(atf_bin1$grid_size_km)

## bind data and calculate PESC for each grid cell
atf_PESC <- atf_bin1 %>% 
  left_join(atf_bin2, by = c("E_km", "N_km", "lat", "lon", "grid_size_km", "year"), 
            suffix = c("_bin1", "_bin2")) %>% 
  filter(year == 2009 | year == 2015) %>% 
  mutate(PESC_bin1 = D_gct_1_bin1 * D_gct_2_bin1, 
         PESC_bin2 = D_gct_1_bin2 * D_gct_2_bin2,
         PESC = PESC_bin1 + PESC_bin2, 
         PWSC = PESC/(D_gct_1_bin1 + D_gct_1_bin2))


### Create individual panel plots -------------------------------------------------------

## Plot bottom temperature
bt_limits <- range(bottom_temp$GEAR_TEMPERATURE)

bt_2009 <- with(
  bottom_temp[bottom_temp$YEAR == 2009,], 
  map_EBS_survey(STATIONID, color = GEAR_TEMPERATURE, limits = bt_limits)  
) + ggtitle("a. bottom temperature 2009") + 
  theme(axis.text.x = element_blank())

bt_2015 <- with(
  bottom_temp[bottom_temp$YEAR == 2015,], 
  map_EBS_survey(STATIONID, color = GEAR_TEMPERATURE, limits = bt_limits)  
) + ggtitle("e. bottom temperature 2015")

## Plot juvenile walleye pollock

pollock_data <- delta_models$grid_cell_fits$Juvenile_Walleye_pollock %>% filter(year == 2009 | year == 2015)

pollock_2009 <- (with(
  pollock_data[pollock_data$year == 2009,],
  map_EBS_extrap(E_km, N_km, fill = R1_gct, grid_size = extrap_size)
) + geom_contour(aes(x = E_km, y = N_km, z = temp_2009), 
                 breaks = c(-4, 2, 10), color = "white", size = 0.8, 
                 data = EBS_grid) + 
  geom_contour(aes(x = E_km, y = N_km, z = temp_2009), 
               breaks = c(-4, 0, 10), color = "white", size = 0.8, 
               data = EBS_grid, linetype = "dashed") + 
  ggtitle("b. juvenile pollock 2009") + 
  scale_fill_viridis_c(limits = c(0, 1)) +
  scale_color_viridis_c(limits = c(0, 1))) %>% 
  theme_map() + 
  guides(color = "none") + 
  theme(axis.text = element_blank())

pollock_2015 <- (with(
  pollock_data[pollock_data$year == 2015,],
  map_EBS_extrap(E_km, N_km, fill = R1_gct, grid_size = extrap_size)
) + geom_contour(aes(x = E_km, y = N_km, z = temp_2015), 
                 breaks = c(-4, 2, 10), color = "white", size = 0.8, 
                 data = EBS_grid) + 
  geom_contour(aes(x = E_km, y = N_km, z = temp_2015), 
               breaks = c(-4, 0, 10), color = "white", size = 0.8, 
               data = EBS_grid, linetype = "dashed") +  
  ggtitle("f. juvenile pollock 2015") + 
  scale_fill_viridis_c(limits = c(0, 1)) +
  scale_color_viridis_c(limits = c(0, 1))) %>%  
  theme_map() + 
  guides(color = "none") + 
  theme(axis.text.y = element_blank())


## Plot arrowtooth flounder

atf_data <- delta_models$grid_cell_fits$Flounder %>% filter(year == 2009 | year == 2015)

atf_2009 <- (with(
  atf_data[atf_data$year == 2009,],
  map_EBS_extrap(E_km, N_km, fill = R1_gct, grid_size = extrap_size)
) + geom_contour(aes(x = E_km, y = N_km, z = temp_2009), 
                 breaks = c(-4, 2, 10), color = "white", size = 0.8, 
                 data = EBS_grid) + 
  geom_contour(aes(x = E_km, y = N_km, z = temp_2009), 
               breaks = c(-4, 0, 10), color = "white", size = 0.8, 
               data = EBS_grid, linetype = "dashed") + 
  ggtitle("c. flounder 2009") + 
  scale_fill_viridis_c(limits = c(0, 1)) +
  scale_color_viridis_c(limits = c(0, 1))) %>%
  theme_map() + 
  guides(color = "none") + 
  theme(axis.text = element_blank())


atf_2015 <- (with(
  atf_data[atf_data$year == 2015,],
  map_EBS_extrap(E_km, N_km, fill = R1_gct, grid_size = extrap_size)
) + geom_contour(aes(x = E_km, y = N_km, z = temp_2015), 
                 breaks = c(-4, 2, 10), color = "white", size = 0.8, 
                 data = EBS_grid) + 
  geom_contour(aes(x = E_km, y = N_km, z = temp_2015), 
               breaks = c(-4, 0, 10), color = "white", size = 0.8, 
               data = EBS_grid, linetype = "dashed") + 
  ggtitle("g. flounder 2015") + 
  scale_fill_viridis_c(limits = c(0, 1)) +
  scale_color_viridis_c(limits = c(0, 1))) %>% 
  theme_map() + 
  guides(color = "none") + 
  theme(axis.text.y = element_blank())


## Plot area overlap

source(here("src", "range_overlap_functions.R"))
thresholds <- read.csv(here("model output", "delta_VAST_output", "fit_thresholds.csv"))
thresholds <- setNames(thresholds$threshold, thresholds$species)

overlap_value <- list(
  "2009" = round(area_overlapfn(
    pollock_data$R1_gct[pollock_data$year == 2009] > thresholds["Juvenile_Walleye_pollock"], 
    atf_data$R1_gct[atf_data$year == 2009] > thresholds["Flounder"], 
    area = rep(extrap_size, sum(atf_data$year == 2009))
  ), 2), 
  "2015" = round(area_overlapfn(
    pollock_data$R1_gct[pollock_data$year == 2015] > thresholds["Juvenile_Walleye_pollock"], 
    atf_data$R1_gct[atf_data$year == 2015] > thresholds["Flounder"], 
    area = rep(extrap_size, sum(atf_data$year == 2015))
  ), 2)
)

pollock_data$area_overlap <- pollock_data$R1_gct > thresholds["Juvenile_Walleye_pollock"] & atf_data$R1_gct > thresholds["Flounder"]

area_colors <- c("#04386b", "#0a6bc9")

area_2009 <- with(
  pollock_data[pollock_data$year == 2009,],
  map_EBS_extrap(E_km, N_km, fill = area_overlap, grid_size = extrap_size)
) + geom_contour(aes(x = E_km, y = N_km, z = temp_2009), 
                 breaks = c(-4, 2, 10), color = "white", size = 0.8, 
                 data = EBS_grid) + 
  geom_contour(aes(x = E_km, y = N_km, z = temp_2009), 
               breaks = c(-4, 0, 10), color = "white", size = 0.8, 
               data = EBS_grid, linetype = "dashed") + 
  ggtitle("d. area overlap 2009") + 
  scale_fill_manual(values = area_colors) + 
  scale_color_manual(values = area_colors) + 
  geom_text(aes(x = 1150, y = 6800, label = overlap_value[["2009"]]), size = 6, 
            color = "grey60") + 
  theme_bw() + 
  theme(panel.border = element_rect(color = "grey80", size = 1),
        panel.background = element_rect(fill = "grey80"),
        panel.spacing.x = unit(0, "lines"),
        strip.text = element_text(size = 12, face = "bold"), 
        strip.background = element_blank(), 
        axis.text = element_blank(),
        legend.title = element_text(size = 14),
        axis.ticks.length = unit(0.4, "lines"),
        axis.ticks = element_line(size = 1, color = "grey80"),
        axis.title = element_blank(),
        panel.grid = element_blank(), 
        plot.title = element_text(size = 16, face = "bold"),
        legend.position = "bottom") + 
  guides(fill = guide_legend(title.position = "top", title.hjust = 0.5), 
         color = "none")

area_2015 <- with(
  pollock_data[pollock_data$year == 2015,],
  map_EBS_extrap(E_km, N_km, fill = area_overlap, grid_size = extrap_size)
) + geom_contour(aes(x = E_km, y = N_km, z = temp_2015), 
                 breaks = c(-4, 2, 10), color = "white", size = 0.8, 
                 data = EBS_grid) + 
  geom_contour(aes(x = E_km, y = N_km, z = temp_2015), 
               breaks = c(-4, 0, 10), color = "white", size = 0.8, 
               data = EBS_grid, linetype = "dashed") +  
  ggtitle("h. area overlap 2015") + 
  scale_fill_manual(values = area_colors) + 
  scale_color_manual(values = area_colors) + 
  geom_text(aes(x = 1150, y = 6800, label = overlap_value[["2015"]]), size = 6, 
            color = "grey60") + 
  theme_bw() + 
  theme(panel.border = element_rect(color = "grey80", size = 1),
        panel.background = element_rect(fill = "grey80"),
        panel.spacing.x = unit(0, "lines"),
        strip.text = element_text(size = 12, face = "bold"), 
        strip.background = element_blank(), 
        axis.text.x = element_text(size = 12, color = "grey70", face = "bold"),
        axis.text.y = element_blank(),
        legend.title = element_text(size = 14),
        axis.ticks.length = unit(0.4, "lines"),
        axis.ticks = element_line(size = 1, color = "grey80"),
        axis.title = element_blank(),
        panel.grid = element_blank(), 
        plot.title = element_text(size = 16, face = "bold"),
        legend.position = "bottom") + 
  guides(fill = guide_legend(title.position = "top", title.hjust = 0.5, 
                             label.position = "bottom", keyheight = unit(1.45, "lines")), 
         color = "none")



### Bind plots together -------------------------------------------------------

## Bind full plot
cold_pool_plot <- cowplot::plot_grid(
  bt_2009 + theme(legend.position = "none"), 
  bt_2015 + labs(color = expression("bottom temperature (°C)"), # phantom for spacing
                 fill = expression("bottom temperature (°C)")), 
  pollock_2009 + theme(legend.position = "none"), 
  pollock_2015 + labs(fill = expression("probability of encounter")),
  atf_2009 + theme(legend.position = "none"), 
  atf_2015 + labs(fill = expression("probability of encounter")),
  area_2009 + theme(legend.position = "none"), 
  area_2015,
  ncol = 4, byrow = FALSE, 
  rel_heights = rep(c(1, 1.35), 4), 
  rel_widths = c(1.15, rep(1, 3), 1.15, rep(1, 3))
)

## Save  plot
dir.create(here("model output", "PESC_overlap"), recursive = TRUE)
ggsave(here("model output", "PESC_overlap", "atf_example_plot2.png"), 
       cold_pool_plot, height = 9, width = 16, units = "in", 
       dpi = 500, scale = 0.9)
