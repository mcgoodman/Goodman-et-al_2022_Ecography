
###############################################################################
##' @title Regression and plots of predation vs overlap metrics
##' @author Maurice Goodman
###############################################################################

## Note: "PESC," or "predator-expanded-stomach-contents," corresponds
## to the metric of "total predation" used in Goodman et al.
## "PWSC," or "predator-weighted-stomach-contents" corresponds to
## the metric of "relative predation" used in Goodman et al.

### Load and install packages -------------------------------------------------

packages <- c("here", "tidyverse", "cowplot", "mapdata")
for (i in seq_along(packages)) {
  if(!do.call(require, list(package = packages[i]))) {
    do.call(install.packages, list(pkgs = packages[i]))
    do.call(require, list(package = packages[i]))
  }
}

### Load and join PESC/biomass/PWSC MLE datasets ------------------------------

## list PESC MLE and SD datasets
PESC_MLE <- list.files(here("model output", "PESC"), recursive = TRUE, 
                       pattern = "PESC_MLE_SD.csv", full.names = TRUE)


### read in PESC MLE and SD datasets
PESC_MLE <- lapply(PESC_MLE, read.csv, stringsAsFactors = FALSE)

## names and size bins of PESC species
PESC_names <- list.files(here("model output", "PESC"))

## extract species names
PESC_species <- vapply(
  PESC_names, 
  function(x) paste(strsplit(x, "_")[[1]][1:2], collapse = " "), 
  FUN.VALUE = character(1), 
  USE.NAMES = FALSE
)

## extract size bins
PESC_bins <- vapply(
  PESC_names, 
  function(x) as.numeric(str_sub(x, nchar(x), nchar(x))), 
  FUN.VALUE = numeric(1), 
  USE.NAMES = FALSE
)

names(PESC_MLE) <- paste(PESC_species, PESC_bins, sep = "-")

## Add MLE across size bins, and recalculate PWSC
## SD has to be re-calculated from summed predictive distributions (below)
PESC_MLE <- PESC_MLE %>% 
  bind_rows(.id = "id") %>% 
  separate(id, c("species", "bin"), sep = "-") %>% 
  select(-SD) %>% 
  pivot_wider(names_from = metric, values_from = MLE) %>% 
  rename(PESC_t = `PESC (tonnes)`, biomass_t = `biomass (tonnes)`) %>% 
  group_by(species, year) %>% 
  ## check that predation models exist for both size bins, remove if not 
  mutate(n_bins = length(unique(bin))) %>% 
  filter(n_bins == 2) %>% select(-n_bins) %>%
  summarize_at(vars(PESC_t, biomass_t), sum) %>% 
  mutate(PWSC = PESC_t/biomass_t) %>% 
  pivot_longer(c(biomass_t, PESC_t, PWSC), names_to = "metric", values_to = "MLE")

### Load and join PESC/biomass predictive distribution datasets ---------------

## list PESC predictive sample datasets
PESC_preds <- list.files(
  here("model output", "PESC"), 
  pattern = "PESC_predictive_samples.csv", 
  recursive = TRUE, 
  full.names = TRUE
)

## read in PESC predictive samples
PESC_preds <- lapply(PESC_preds, read.csv, stringsAsFactors = FALSE)
names(PESC_preds) <- paste(PESC_species, PESC_bins, sep = "-")

## Sum biomass, PESC, and PWSC across size bins for each species, year, and sample
PESC_pred_df <- PESC_preds %>% 
  bind_rows(.id = "id") %>% 
  separate(id, c("species", "bin"), sep = "-") %>% 
  group_by(species, sample, year) %>%
  summarize(biomass_t = sum(biomass_t), 
            PESC_t = sum(PESC_t)) %>% 
  mutate(PWSC = PESC_t/biomass_t)

## Calculate SD for PESC/ biomass / PWSC from predictive 
## distribution and join SD dataset to MLE dataset
PESC_MLE_SD <- PESC_pred_df %>% 
  pivot_longer(c(biomass_t, PESC_t, PWSC), names_to = "metric") %>% 
  group_by(species, year, metric) %>%
  summarize(SD = sd(value)) %>% 
  left_join(PESC_MLE) %>% 
  na.omit()

## directory to save data and plots to
save_dir <- here("model output", "PESC_overlap")
dir.create(save_dir)

## write out PESC, biomass, and PWSC MLE and SD
write.csv(PESC_MLE_SD, paste(save_dir, "PESC_MLE_SD.csv", sep = "/"), row.names = FALSE)

### Load overlap predictive distributions and join with PESC/biomass ----------

## Read in overlap metrics predictive samples
delta_preds <- readRDS(here("model output", "overlap_metrics", "overlap_metrics.rds"))

## Species names differ between PESC and distribution models
## List species names for both datasets
delta_species <- c(
  "Arrowtooth Flounder", "Pacific Cod", "Pacific Halibut", "Walleye Pollock"
)
delta_species <- setNames(
  delta_species, c("Flounder", "Pacific_cod", "Pacific_halibut", "Walleye_pollock")
)

## Change overlap metrics species names to match PESC
delta_preds$samples$species <- plyr::revalue(delta_preds$samples$species, delta_species)
delta_preds$MLE$species <- plyr::revalue(delta_preds$MLE$species, delta_species)

## predictive samples of juvenile pollock biomass (for multiple regressions)\
biomass_samples <- readRDS(here("model output", "delta_VAST_output", "biomass.rds"))
pollock_samples <- biomass_samples$samples %>% 
  filter(species == "Juvenile_Walleye_pollock") %>% 
  select(-species, pollock = biomass)
pollock_MLE <- biomass_samples$MLE %>% 
  filter(species == "Juvenile_Walleye_pollock") %>% 
  select(-species, pollock = biomass)

## Join datasets
preds_join <- PESC_pred_df %>% 
  left_join(delta_preds$samples, by = c("species", "sample", "year")) %>%
  left_join(pollock_samples, by = c("sample", "year"))

## write out joined datasets
write.csv(preds_join, paste(save_dir, "PESC_overlap.csv", sep = "/"), row.names = FALSE)

### Join overlap MLE and SD from predictive distribution ----------------------

overlap_MLE_SD <- delta_preds$samples %>% 
  pivot_longer(-c(species, year, sample), names_to = "metric") %>% 
  group_by(species, year, metric) %>% 
  summarize(SD = sd(value)) %>% 
  left_join(delta_preds$MLE %>% 
              pivot_longer(-c(species, year), names_to = "metric", values_to = "MLE"), 
            by = c("species", "year", "metric"))

### Time series PESC, biomass, and PWSC plot for each species -----------------

pred_cats <- c("PESC_t", "biomass_t", "PWSC")

ts_plots <- vector("list", length(pred_cats))

for (i in seq_along(pred_cats)) {
  ts_plots[[i]] <- PESC_MLE_SD %>% 
    filter(metric == pred_cats[i]) %>% 
    mutate(MLE = ifelse(metric == "PWSC", MLE, MLE/1e+06), 
           SD = ifelse(metric == "PWSC", SD, SD/1e+06)) %>% 
    ggplot(aes(year, MLE)) + 
    geom_line(color = "black") + 
    geom_point(color = "black") + 
    geom_errorbar(aes(ymin = MLE - SD, ymax = MLE + SD), color = "black", width = 0.4) + 
    facet_wrap(~species, nrow = 1, scales = "free_y") + 
    theme_bw() +
    theme(strip.background = element_blank(), 
          strip.text = element_text(size = 12, face = "bold"),
          axis.ticks = element_line(color = "black"), 
          axis.text = element_text(color = "black"), 
          panel.border = element_rect(color = "black"), 
          panel.grid.minor = element_blank()) + 
    labs(y = ifelse(
      pred_cats[i] == "PESC_t", "total predation (Mt)", ifelse(
        pred_cats[i] == "biomass_t", "biomass (Mt)", 
        expression(paste("relative predation ", "(", "kg kg"^{-1}, ")"))
      )
    )) +  
    scale_x_continuous(breaks = seq(1985, 2015, 5), 
                       labels = function(x) str_sub(as.character(x), 3, 4))
  
  ## modifications to individal plots
  if (i > 1) ts_plots[[i]] <- ts_plots[[i]] + theme(strip.text = element_blank())
  if (i < 3) ts_plots[[i]] <- ts_plots[[i]] + theme(axis.title.x = element_blank())

  
  ts_plots[[i]] <- ggplotGrob(ts_plots[[i]])
}

ts_plot <- cowplot::plot_grid(plotlist = ts_plots, ncol = 1, align = "v", 
                              rel_heights = c(1, 0.9, 1))

## Save overlap vs standardized PESC plots
ggsave(paste(save_dir, "PESC_PWSC_biomass_ts.png", sep = "/"), ts_plot, 
       height = 7, width = 10, units = "in", dpi = 500)

### Time series overlap plot for each species ---------------------------------

### Read in cold pool data
cold_pool <- read.csv(here("data", "cold_pool_extent_annual.csv"))

cold_pool <- cold_pool %>%  select(
  year = YEAR,
  cold_pool = AREA_SUM_KM2_LTE2, 
) %>% mutate(
  cold_pool_sd = c(scale(cold_pool))
)

### List overlap metrics to plot
overlap_metrics <- c(
  "area_overlap", "global_colloc", "loc_colloc"
)
names(overlap_metrics) <- c(
  "area overlap", "global index of collocation", "local index of collocation"
)
overlap_labels <- c(
  "area overlap", "global index of\ncollocation", "local index of\ncollocation"
)

### Plot of overlap with points colored by standardized cold pool extent
overlap_plot <- overlap_MLE_SD %>%
  filter(metric %in% overlap_metrics) %>% 
  mutate(metric = plyr::revalue(metric, setNames(overlap_labels, overlap_metrics))) %>% 
  left_join(cold_pool, by = "year") %>% 
  ggplot(aes(year, MLE)) + 
  geom_line() + 
  geom_errorbar(aes(ymin = MLE - SD, ymax = MLE + SD), width = 0.4) + 
  geom_point(aes(fill = cold_pool_sd), color = "black", shape = 21, size = 2) + 
  facet_grid(metric ~ species, scales = "free", switch = "y") + 
  theme_bw() + 
  theme(strip.background = element_blank(), 
        strip.text = element_text(size = 12, face = "bold"),
        strip.placement = "outside", 
        axis.title.y = element_blank(), 
        legend.position = "bottom") +
  scale_x_continuous(breaks = seq(1985, 2015, 5), 
                     labels = function(x) str_sub(as.character(x), 3, 4)) + 
  scale_fill_fermenter(breaks = c(-Inf, -2, -1, 0, 1, 2, Inf), 
                       palette = "RdYlBu", direction = 1, 
                       guide = guide_colorsteps(title.position = "top", title.hjust = 0.5, ticks = TRUE,
                                                frame.colour = "black", frame.linewidth = 1.5, 
                                                ticks.colour = "black", ticks.linewidth = 1.5, 
                                                barwidth = unit(10, "lines")), 
                       labels = function(x) paste(x, "SD")) + 
  labs(fill = "cold pool extent")

ggsave(paste(save_dir, "overlap_ts_cold_pool.png", sep = "/"), overlap_plot, height = 7, 
       width = 10, units = "in", dpi = 500)

## Plot of overlap without cold pool
overlap_plot <- overlap_MLE_SD %>%
  filter(metric %in% overlap_metrics) %>% 
  mutate(metric = plyr::revalue(metric, setNames(overlap_labels, overlap_metrics))) %>% 
  ggplot(aes(year, MLE)) + 
  geom_line() + 
  geom_errorbar(aes(ymin = MLE - SD, ymax = MLE + SD), width = 0.4) + 
  geom_point() + 
  facet_grid(metric ~ species, scales = "free", switch = "y") + 
  theme_bw() + 
  theme(strip.background = element_blank(), 
        strip.text = element_text(size = 12, face = "bold"),
        strip.placement = "outside", 
        axis.title.y = element_blank()) +
  scale_x_continuous(breaks = seq(1985, 2015, 5), 
                     labels = function(x) str_sub(as.character(x), 3, 4))

ggsave(paste(save_dir, "overlap_ts.png", sep = "/"), overlap_plot, height = 6, 
       width = 10, units = "in", dpi = 500)

## Regressions of overlap against cold pool extent
overlap_coldpool_lm <- overlap_MLE_SD %>%
  filter(metric %in% overlap_metrics) %>% 
  mutate(metric = plyr::revalue(metric, setNames(overlap_labels, overlap_metrics))) %>% 
  left_join(cold_pool, by = "year") %>% ungroup() %>% 
  group_by(species, metric) %>% nest() %>% 
  mutate(lm = purrr::map(data, function(df) lm(scale(MLE) ~ cold_pool_sd, data = df)), 
         slope = purrr::map_dbl(lm, function(x) coef(x)[2]), 
         lower = purrr::map_dbl(lm, function(x) confint(x)[2,1]), 
         upper = purrr::map_dbl(lm, function(x) confint(x)[2,2]), 
         t = purrr::map_dbl(lm, function(x) summary(x)$coefficients[2,"t value"]),
         p = purrr::map_dbl(lm, function(x) summary(x)$coefficients[2,"Pr(>|t|)"])) %>% 
  select(-data, -lm) %>% 
  ungroup()

write.csv(overlap_coldpool_lm, paste(save_dir, "overlap_coldpool_coefs.csv", sep = "/"), 
          row.names = FALSE)

### Plot of overlap against cold pool extent
overlap_coldpool_xyplot <- overlap_MLE_SD %>%
  filter(metric %in% overlap_metrics) %>% 
  mutate(metric = plyr::revalue(metric, setNames(overlap_labels, overlap_metrics))) %>% 
  left_join(cold_pool, by = "year") %>% 
  ggplot(aes(cold_pool/1e+03, MLE, color = metric, fill = metric)) + 
  geom_point() + 
  facet_wrap(~species, nrow = 1) + 
  theme_bw() +
  geom_smooth(method = "lm") + 
  theme(strip.background = element_blank(), 
        strip.text = element_text(size = 10, face = "bold"),
        axis.text = element_text(color = "black"), 
        axis.ticks = element_line(color = "black"), 
        panel.border = element_rect(color = "black"),
        legend.position = "bottom") + 
  labs(x = expression(paste("cold pool extent ", "(", "km \U00D7 10"^{3},")")), 
       y = "overlap") +
  scale_color_manual(values = c("#50A3A4", "#FCAF38", "#F95335")) + 
  scale_fill_manual(values = c("#50A3A4", "#FCAF38", "#F95335")) + 
  coord_cartesian(clip = "off")

ggsave(
  paste(save_dir, "overlap_coldpool_xyplot.png", sep = "/"),
  plot = overlap_coldpool_xyplot, height = 4, width = 10, units = "in", dpi = 500, scale = 0.9
)

## coefficient plot of overlap against cold pool extent
overlap_coldpool_lm_plot <- overlap_coldpool_lm %>% 
  ggplot(aes(species, slope)) + 
  facet_wrap(~metric) + 
  geom_hline(aes(yintercept = 0), linetype = "dashed") + 
  geom_errorbar(aes(ymin = lower, ymax = upper), size = 1, width = 0, 
                position = position_dodge(width = 0.4)) + 
  geom_point(size = 3, shape = 21, stroke = 1.5, fill = "white", 
             position = position_dodge(width = 0.4)) + 
  coord_flip(clip = "off") + 
  theme_bw() + 
  theme(strip.background = element_blank(), 
        strip.text = element_text(face = "bold"),
        axis.title.y = element_blank(),  
        axis.text = element_text(color = "black"),
        axis.text.y = element_text(face = "bold", color = "black"), 
        panel.grid.major.y = element_blank()) 

ggsave(paste(save_dir, "overlap_coldpool_coefs.png", sep = "/"), overlap_coldpool_lm_plot, 
       height = 4, width = 8, units = "in", dpi = 500)

### Scaled Per-sample regression of PESC/PWSC on overlap ----------------------

## Datasets with PESC & overlap MLE & SD - for plotting and summary
eiv_pesc_data <- vector("list", length(overlap_metrics))
names(eiv_pesc_data) <- names(overlap_metrics)

for (i in seq_along(overlap_metrics)) {
  
  eiv_pesc_data[[i]] <- PESC_MLE_SD %>% 
    filter(metric == "PESC_t") %>% 
    select(species, year, y_MLE = MLE, y_SD = SD)
  
  eiv_pesc_data[[i]] <- overlap_MLE_SD %>% 
    filter(metric == overlap_metrics[i]) %>% 
    select(species, year, x_MLE = MLE, x_SD = SD) %>% 
    left_join(eiv_pesc_data[[i]], by = c("species", "year")) %>% 
    left_join(pollock_MLE, by = "year")
  
}

## Datasets with PWSC & overlap MLE & SD - for plotting and summary
eiv_pwsc_data <- vector("list", length(overlap_metrics))
names(eiv_pwsc_data) <- names(overlap_metrics)

for (i in seq_along(overlap_metrics)) {
  
  eiv_pwsc_data[[i]] <- PESC_MLE_SD %>% 
    filter(metric == "PWSC") %>% 
    select(species, year, y_MLE = MLE, y_SD = SD)
  
  eiv_pwsc_data[[i]] <- overlap_MLE_SD %>% 
    filter(metric == overlap_metrics[i]) %>% 
    select(species, year, x_MLE = MLE, x_SD = SD) %>% 
    left_join(eiv_pwsc_data[[i]], by = c("species", "year")) %>% 
    left_join(pollock_MLE, by = "year")
  
}

overlap_lm <- vector("list", length(overlap_metrics)*length(delta_species))

n_samples <- max(preds_join$sample)

index <- 1

for (i in seq_along(overlap_metrics)) {
  for (j in seq_along(delta_species)) {
    
    coef_MLE_pesc <- coef(
      lm(y_MLE ~ x_MLE, 
         data = eiv_pesc_data[[names(overlap_metrics)[i]]] %>% ungroup() %>% 
           filter(species == delta_species[j] & !is.na(y_MLE)) %>% 
           mutate(y_MLE = as.numeric(scale(y_MLE)), 
                  x_MLE = as.numeric(scale(x_MLE))))
    )
    
    coef_MLE_pwsc <- coef(
      lm(y_MLE ~ x_MLE, 
         data = eiv_pwsc_data[[names(overlap_metrics)[i]]] %>% ungroup() %>% 
           filter(species == delta_species[j] & !is.na(y_MLE)) %>% 
           mutate(y_MLE = as.numeric(scale(y_MLE)), 
                  x_MLE = as.numeric(scale(x_MLE))))
    )
    
    ## list to store regression coefficients
    coef_samples_pesc <- matrix(nrow = n_samples, ncol = 2)
    coef_samples_pwsc <- matrix(nrow = n_samples, ncol = 2)
    
    ## per-sample regressions for PESC and PWSC
    for (k in seq(n_samples)) {
      
      ## data for per-sample regressions
      lm_data <- preds_join %>% filter(species == delta_species[j] & sample == k)
      lm_data$PESC <- as.numeric(scale(lm_data$PESC_t))
      lm_data$PWSC <- as.numeric(scale(lm_data$PWSC))
      lm_data$overlap <- scale(lm_data[[overlap_metrics[i]]])
      
      coef_samples_pesc[k,] <- coef(lm(
        PESC ~ overlap, 
        data = lm_data[lm_data$sample == k, ]
      ))
      
      coef_samples_pwsc[k,] <- coef(lm(
        PWSC ~ overlap, 
        data = lm_data[lm_data$sample == k, ]
      ))
      
    }
    
    ## variance and covariance of PESC coefficients
    coef_var_pesc <- apply(coef_samples_pesc, MARGIN = 2, FUN = var)
    coef_cov_pesc <- cov(coef_samples_pesc[,1], coef_samples_pesc[,2])
    
    ## variance and covariance of PWSC coefficients
    coef_var_pwsc <- apply(coef_samples_pwsc, MARGIN = 2, FUN = var)
    coef_cov_pwsc <- cov(coef_samples_pwsc[,1], coef_samples_pwsc[,2])
    
    ## output summary
    overlap_lm[[index]] <- data.frame(
      species = rep(delta_species[j], 2), 
      y = c("PESC", "PWSC"),
      x = rep(names(overlap_metrics)[i], 2),
      intercept = c(coef_MLE_pesc[1], coef_MLE_pwsc[1]), 
      intercept_var = c(coef_var_pesc[1], coef_var_pwsc[1]),
      slope = c(coef_MLE_pesc[2], coef_MLE_pwsc[2]), 
      slope_var = c(coef_var_pesc[2], coef_var_pwsc[2]), 
      cov = c(coef_cov_pesc, coef_cov_pwsc)
    )  
    
    index <- index + 1
  }
}

overlap_lm <- do.call("rbind", overlap_lm)

coef_summary <- overlap_lm %>% 
  mutate(`2.5%` = qnorm(0.025, slope, sqrt(slope_var)), 
         `10%` = qnorm(0.1, slope, sqrt(slope_var)), 
         `90%` = qnorm(0.9, slope, sqrt(slope_var)), 
         `97.5%` = qnorm(0.975, slope, sqrt(slope_var)), 
         Z = slope/sqrt(slope_var), 
         p = 2 * pnorm(-abs(Z), lower.tail = TRUE)) %>% 
  mutate(y = ifelse(y == "PESC", "total", "relative"))

coef_plot <- coef_summary %>% 
  ggplot(aes(species, slope, color = y)) + 
  facet_wrap(~x) + 
  geom_hline(aes(yintercept = 0), linetype = "dashed") + 
  geom_errorbar(aes(ymin = `2.5%`, ymax = `97.5%`), size = 1, width = 0, 
                 position = position_dodge(width = 0.4)) + 
  geom_errorbar(aes(ymin = `10%`, ymax = `90%`), size = 2, width = 0, 
                 position = position_dodge(width = 0.4)) + 
  geom_point(size = 3, shape = 21, stroke = 1.5, fill = "white", 
             position = position_dodge(width = 0.4)) + 
  geom_text(aes(x = species, y = `97.5%` + 0.05, label = y), 
            data = coef_summary %>% 
              filter(species == "Walleye Pollock" & x == "area overlap"), 
            position = position_dodge(width = 0.4), 
            hjust = 0) + 
  coord_flip(clip = "off") + 
  theme_bw() + 
  theme(strip.background = element_blank(), 
        strip.text = element_text(face = "bold"),
        axis.title.y = element_blank(),  
        axis.text = element_text(color = "black"),
        axis.text.y = element_text(face = "bold", color = "black"), 
        panel.grid.major.y = element_blank(), 
        legend.position = "none") + 
  scale_color_manual(values = c("grey50", "black")) + 
  scale_y_continuous(breaks = seq(-0.6, 0.6, 0.3))

## Save coefficients plot
ggsave(paste(save_dir, "overlap_PWSC_PESC_coefs.png", sep = "/"), coef_plot, 
       height = 4, width = 8, units = "in", dpi = 500)

### Scaled multiple regressions of PESC/PWSC on overlap -----------------------

## regressions include juvenile pollock biomass

multiple_lm <- vector("list", length(overlap_metrics)*length(delta_species))

index <- 1

for (i in seq_along(overlap_metrics)) {
  for (j in seq_along(delta_species)) {
    
    coef_MLE_pesc <- coef(
      lm(y_MLE ~ x_MLE + pollock, 
         data = eiv_pesc_data[[names(overlap_metrics)[i]]] %>% ungroup() %>% 
           filter(species == delta_species[j] & !is.na(y_MLE)) %>% 
           mutate(y_MLE = as.numeric(scale(y_MLE)), 
                  x_MLE = as.numeric(scale(x_MLE)), 
                  pollock = as.numeric(scale(pollock))))
    )
    
    coef_MLE_pwsc <- coef(
      lm(y_MLE ~ x_MLE + pollock, 
         data = eiv_pwsc_data[[names(overlap_metrics)[i]]] %>% ungroup() %>% 
           filter(species == delta_species[j] & !is.na(y_MLE)) %>% 
           mutate(y_MLE = as.numeric(scale(y_MLE)), 
                  x_MLE = as.numeric(scale(x_MLE)), 
                  pollock = as.numeric(scale(pollock))))
    )
    
    ## list to store regression coefficients
    coef_samples_pesc <- matrix(nrow = n_samples, ncol = 3)
    coef_samples_pwsc <- matrix(nrow = n_samples, ncol = 3)
    
    ## per-sample regressions for PESC and PWSC
    for (k in seq(n_samples)) {
      
      ## data for per-sample regressions
      lm_data <- preds_join %>% filter(species == delta_species[j] & sample == k)
      lm_data$PESC <- as.numeric(scale(lm_data$PESC_t))
      lm_data$PWSC <- as.numeric(scale(lm_data$PWSC))
      lm_data$overlap <- scale(lm_data[[overlap_metrics[i]]])
      lm_data$pollock <- as.numeric(scale(lm_data$pollock))
      
      coef_samples_pesc[k,] <- coef(lm(
        PESC ~ overlap + pollock, 
        data = lm_data
      ))
      
      coef_samples_pwsc[k,] <- coef(lm(
        PWSC ~ overlap + pollock, 
        data = lm_data
      ))
      
    }
    
    ## variance and covariance of PESC coefficients
    coef_var_pesc <- apply(coef_samples_pesc, MARGIN = 2, FUN = var)
    coef_cov_pesc <- cov(coef_samples_pesc[,1], coef_samples_pesc[,2])
    
    ## variance and covariance of PWSC coefficients
    coef_var_pwsc <- apply(coef_samples_pwsc, MARGIN = 2, FUN = var)
    coef_cov_pwsc <- cov(coef_samples_pwsc[,1], coef_samples_pwsc[,2])
    
    ## output summary
    multiple_lm[[index]] <- data.frame(
      species = rep(delta_species[j], 2), 
      y = c("PESC", "PWSC"),
      x = rep(names(overlap_metrics)[i], 2),
      intercept = c(coef_MLE_pesc[1], coef_MLE_pwsc[1]), 
      intercept_var = c(coef_var_pesc[1], coef_var_pwsc[1]),
      slope = c(coef_MLE_pesc[2], coef_MLE_pwsc[2]), 
      slope_var = c(coef_var_pesc[2], coef_var_pwsc[2]), 
      cov = c(coef_cov_pesc, coef_cov_pwsc)
    )  
    
    index <- index + 1
  }
}

multiple_lm <- do.call("rbind", multiple_lm)

coef_summary_multiple <- multiple_lm %>% 
  mutate(`2.5%` = qnorm(0.025, slope, sqrt(slope_var)), 
         `10%` = qnorm(0.1, slope, sqrt(slope_var)), 
         `90%` = qnorm(0.9, slope, sqrt(slope_var)), 
         `97.5%` = qnorm(0.975, slope, sqrt(slope_var))) %>% 
  mutate(y = ifelse(y == "PESC", "total", "relative"))

coef_plot_multiple <- coef_summary_multiple %>% 
  ggplot(aes(species, slope, color = y)) + 
  facet_wrap(~x) + 
  geom_hline(aes(yintercept = 0), linetype = "dashed") + 
  geom_errorbar(aes(ymin = `2.5%`, ymax = `97.5%`), size = 1, width = 0, 
                position = position_dodge(width = 0.4)) + 
  geom_errorbar(aes(ymin = `10%`, ymax = `90%`), size = 2, width = 0, 
                position = position_dodge(width = 0.4)) + 
  geom_point(size = 3, shape = 21, stroke = 1.5, fill = "white", 
             position = position_dodge(width = 0.4)) + 
  geom_text(aes(x = species, y = `97.5%` + 0.05, label = y), 
            data = coef_summary_multiple %>% 
              filter(species == "Walleye Pollock" & x == "area overlap"), 
            position = position_dodge(width = 0.4), 
            hjust = 0) + 
  coord_flip(clip = "off") + 
  theme_bw() + 
  theme(strip.background = element_blank(), 
        strip.text = element_text(face = "bold"),
        axis.title.y = element_blank(),  
        axis.text = element_text(color = "black"),
        axis.text.y = element_text(face = "bold", color = "black"), 
        panel.grid.major.y = element_blank(), 
        legend.position = "none") + 
  scale_color_manual(values = c("grey50", "black")) + 
  scale_y_continuous(breaks = seq(-0.6, 0.6, 0.3))


## Save coefficients plot
ggsave(paste(save_dir, "overlap_PWSC_PESC_coefs_pollock.png", sep = "/"), 
       coef_plot_multiple, height = 4, width = 8, units = "in", 
       dpi = 500)

### Unscaled per-sample regression of PESC/PWSC on overlap --------------------

overlap_lm2 <- vector("list", length(overlap_metrics)*length(delta_species))

index <- 1

for (i in seq_along(overlap_metrics)) {
  for (j in seq_along(delta_species)) {
    
    coef_MLE_pesc <- coef(
      lm(y_MLE/1e+06 ~ x_MLE, 
         data = eiv_pesc_data[[names(overlap_metrics)[i]]] %>% ungroup() %>% 
           filter(species == delta_species[j] & !is.na(y_MLE)))
    )
    
    coef_MLE_pwsc <- coef(
      lm(y_MLE ~ x_MLE, 
         data = eiv_pwsc_data[[names(overlap_metrics)[i]]] %>% ungroup() %>% 
           filter(species == delta_species[j] & !is.na(y_MLE)))
    )
    
    ## list to store regression coefficients
    coef_samples_pesc <- matrix(nrow = n_samples, ncol = 2)
    coef_samples_pwsc <- matrix(nrow = n_samples, ncol = 2)
    
    ## per-sample regressions for PESC and PWSC
    for (k in seq(n_samples)) {
      
      ## data for per-sample regressions
      lm_data <- preds_join %>% filter(species == delta_species[j] & sample == k)
      lm_data$overlap <- lm_data[[overlap_metrics[i]]]
      
      coef_samples_pesc[k,] <- coef(lm(
        PESC_t/1e+06 ~ overlap, 
        data = lm_data[lm_data$sample == k, ]
      ))
      
      coef_samples_pwsc[k,] <- coef(lm(
        PWSC ~ overlap, 
        data = lm_data[lm_data$sample == k, ]
      ))
      
    }
    
    ## variance and covariance of PESC coefficients
    coef_var_pesc <- apply(coef_samples_pesc, MARGIN = 2, FUN = var)
    coef_cov_pesc <- cov(coef_samples_pesc[,1], coef_samples_pesc[,2])
    
    ## variance and covariance of PWSC coefficients
    coef_var_pwsc <- apply(coef_samples_pwsc, MARGIN = 2, FUN = var)
    coef_cov_pwsc <- cov(coef_samples_pwsc[,1], coef_samples_pwsc[,2])
    
    ## output summary
    overlap_lm2[[index]] <- data.frame(
      species = rep(delta_species[j], 2), 
      y = c("PESC", "PWSC"),
      x = rep(names(overlap_metrics)[i], 2),
      intercept = c(coef_MLE_pesc[1], coef_MLE_pwsc[1]), 
      intercept_var = c(coef_var_pesc[1], coef_var_pwsc[1]),
      slope = c(coef_MLE_pesc[2], coef_MLE_pwsc[2]), 
      slope_var = c(coef_var_pesc[2], coef_var_pwsc[2]), 
      cov = c(coef_cov_pesc, coef_cov_pwsc)
    )  
    
    index <- index + 1
  }
}

overlap_lm2 <- do.call("rbind", overlap_lm2)

### Fitted means and confidence bands -----------------------------------------

## simple function to return a model matrix for a single predictor
model_matrix <- function(x_min, x_max, x_int) {
  x <- seq(x_min, x_max, by = x_int)
  intercept <- rep(1, length(x))
  as.matrix(cbind(intercept, x))
}

## function to calculate fitted mean, se of fit, and confidence intervals
## takes a model matrix `mm` (can be more than 1 predictor)
## takes fitted coefficients `mu` (e.g. from `coef()`)
## takes variance-covariance matrix of coefficients (e.g. from `vcov()`)
predict_lm <- function(mm, mu, sigma, level = 0.95) {
  
  y_hat <- mm %*% mu
  y_se <- sqrt(diag(mm %*% sigma %*% t(mm)))
  
  q <- c(lower = (1 - level)/2, upper = 1 - (1 - level)/2)
  
  y_lower <- qnorm(q[["lower"]], mean = y_hat, sd = y_se)
  y_upper <- qnorm(q[["upper"]], mean = y_hat, sd = y_se)
  
  data.frame(fit = y_hat, se.fit = y_se, lower = y_lower, upper = y_upper)
  
} 

## empty list to store fitted means and confidence bands
fit <- vector("list", nrow(overlap_lm2))

## store fitted means and confidence bands
for (i in seq(nrow(overlap_lm2))) {
  
  species_i <- overlap_lm$species[i]
  x_i <- overlap_lm2$x[i]
  y_i  <- overlap_lm2$y[i]
  
  if (y_i == "PESC") {
    eiv_i <- eiv_pesc_data[[x_i]] %>% filter(species == species_i)
  } else {
    eiv_i <- eiv_pwsc_data[[x_i]] %>% filter(species == species_i)
  }
  
  ### design matrix to predict on
  x_range <- range(eiv_i$x_MLE[!is.na(eiv_i$y_MLE)], na.rm = TRUE)
  X <- model_matrix(x_range[1], x_range[2], (x_range[2] - x_range[1])/500)
  
  ### mean vector and covariance matrices
  mu <- c(overlap_lm2$intercept[i], overlap_lm2$slope[i])
  sigma <- matrix(overlap_lm2$cov[i], nrow = 2, ncol = 2)
  diag(sigma) <- c(overlap_lm2$intercept_var[i], overlap_lm2$slope_var[i])
  
  ### obtain fitted means and confidence bands
  y_fit <- predict_lm(X, mu, sigma)
  
  fit[[i]] <- cbind(species = species_i, overlap_metric = x_i,  
                    pred_metric = y_i, x = X[,2], y_fit)
  
}

fit <- do.call("rbind", fit)


### Plot overlap & PESC -------------------------------------------------------

## function for long tick marks on colorbar
## from user "teunbrand" on Stack overflow
guide_longticks <- function(...) {
  guide <- guide_colorbar(...)
  class(guide) <- c("guide", "guide_longticks", "colorbar")
  guide
}

## function for long tick marks on colorbar
## from user "teunbrand" on Stack overflow
guide_gengrob.guide_longticks <- function(guide, theme) {
  dir <- guide$direction
  guide <- NextMethod()
  is_ticks <- grep("^ticks$", guide$layout$name)
  ticks <- guide$grobs[is_ticks][[1]]
  if (dir == "vertical") {
    ticks$x1 <- rep(tail(ticks$x1, 1), length(ticks$x1))
  } else {
    ticks$y1 <- rep(tail(ticks$y1, 1), length(ticks$y1))
  }
  
  guide$grobs[[is_ticks]] <- ticks
  guide
}

## Empty list to store overlap plots
eiv_pesc_plot <- vector("list", length(overlap_metrics))
names(eiv_pesc_plot) <- names(overlap_metrics)

## Render overlap plots
for (i in seq_along(eiv_pesc_plot)) {
  
  ## extra margins for first plot - panel titles
  if (i == 1) {
    plot_margins <- unit(c(0, 0.5, 0, -0.5), "lines")
  } else {
    plot_margins <- unit(c(-5, 0.5, 0, -0.5), "lines")
  }
  
  fit_i <- fit %>% 
    filter(overlap_metric == names(overlap_metrics)[i] & pred_metric == "PESC")
  
  p_table <- coef_summary %>% 
    filter(x == names(overlap_metrics)[i] & y == "total") %>% 
    mutate(label = ifelse(round(p, 3) == 0, "p < 0.001", paste("p =", round(p, 3))))
  
  label_data <- eiv_pesc_data[[i]] %>% 
    filter(!is.na(y_MLE)) %>% 
    group_by(species) %>% 
    summarize(x = max(x_MLE + x_SD), 
              y = 0.95*max(y_MLE + y_SD)/1e+06) %>% 
    left_join(p_table %>% select(species, p, label), by = "species")
  
  eiv_pesc_plot[[i]] <- eiv_pesc_data[[i]] %>% 
    filter(!is.na(y_MLE)) %>% 
    mutate(y_MLE = y_MLE/1e+06, 
           y_SD = y_SD/1e+06) %>% 
    ggplot(aes(x_MLE, y_MLE, color = year)) + 
    geom_point(size = 1.5) +
    geom_errorbar(aes(ymin = y_MLE - y_SD, ymax = y_MLE + y_SD), size = 0.8) + 
    geom_errorbarh(aes(xmin = x_MLE - x_SD, xmax = x_MLE + x_SD), size = 0.8) + 
    geom_line(aes(x, fit), data = fit_i, inherit.aes = FALSE) + 
    geom_ribbon(aes(x, ymin = lower, ymax = upper), data = fit_i, inherit.aes = FALSE, 
                alpha = 0.3) +
    geom_text(aes(x, y, label = label), data = label_data %>% filter(p >= 0.05), 
              inherit.aes = FALSE, size = 4, color = "grey40",
              hjust = 1) +
    geom_text(aes(x, y, label = label), data = label_data %>% filter(p < 0.05), 
              inherit.aes = FALSE, size = 4, color = "black", fontface = "bold",
              hjust = 1) +
    facet_wrap(~species, scales = "free", nrow = 1) + 
    labs(x = names(overlap_metrics[i])) + 
    theme_bw() +
    theme(strip.background = element_blank(), 
          strip.text = element_text(size = 12, face = "bold"),
          plot.margin = plot_margins,
          axis.title.y = element_blank(), 
          legend.position = "none") + 
    scale_color_viridis_c(limits = range(preds_join$year))
  
  ## remove titles from all panels except the first
  if(i > 1) {
    eiv_pesc_plot[[i]] <- eiv_pesc_plot[[i]] + 
      theme(strip.text = element_blank())
  }
  
  eiv_pesc_plot[[i]] <- ggplotGrob(eiv_pesc_plot[[i]])
}

## Dummy ggplot to extract legend from
preds_legend <- data.frame(
  x = 1:length(unique(preds_join$year)), 
  year = seq(min(preds_join$year), max(preds_join$year))) %>% 
  ggplot(aes(x, x, color = year)) + 
  geom_point() + 
  scale_color_viridis_c(breaks = seq(1985, 2010, 5)) + 
  theme(legend.direction = "horizontal", 
        legend.position = "bottom", 
        legend.title = element_text(face = "bold")) + 
  guides(color = guide_longticks(ticks = TRUE, ticks.colour = "white",
                                 ticks.linewidth = 3, title.position = "top", 
                                 title.hjust = 0.5, barheight = unit(0.5, "lines"),
                                 barwidth = unit(20, "lines")))

## Extract legend
preds_legend <- cowplot::get_legend(preds_legend)

## Object for shared y axis label
preds_ylab <- grid::textGrob("total predation (Mt)", gp = grid::gpar(fontsize = 12), rot = 90)

## Construct full plot
eiv_pesc_plot <-cowplot::plot_grid(
  preds_ylab, 
  cowplot::plot_grid(
    cowplot::plot_grid(plotlist = eiv_pesc_plot, align = "hv", ncol = 1, rel_heights = c(1.05, rep(1, 5))), 
    preds_legend, rel_heights = c(10, 1), ncol = 1
  ),
  ncol = 2, rel_widths = c(0.05, 0.95)
) + theme(panel.border = element_blank())

## Save overlap vs standardized PESC plots
ggsave(paste(save_dir, "PESC_overlap.png", sep = "/"), eiv_pesc_plot, 
       height = 8, width = 10, units = "in", dpi = 500)

### Plot overlap & PWSC -------------------------------------------------------

## Empty vector to store overlap plots
eiv_pwsc_plot <- vector("list", length(overlap_metrics))
names(eiv_pwsc_plot) <- names(overlap_metrics)

## Render overlap plots
for (i in seq_along(eiv_pwsc_plot)) {
  
  ## extra margins for first plot - panel titles
  if (i == 1) {
    plot_margins <- unit(c(0, 0.5, 0, -0.5), "lines")
  } else {
    plot_margins <- unit(c(-5, 0.5, 0, -0.5), "lines")
  }
  
  fit_i <- fit %>% 
    filter(overlap_metric == names(overlap_metrics)[i] & pred_metric == "PWSC")
  
  p_table <- coef_summary %>% 
    filter(x == names(overlap_metrics)[i] & y == "relative") %>% 
    mutate(label = ifelse(round(p, 3) == 0, "p < 0.001", paste("p =", round(p, 3))))
  
  label_data <- eiv_pwsc_data[[i]] %>% 
    filter(!is.na(y_MLE)) %>% 
    group_by(species) %>% 
    summarize(x = max(x_MLE + x_SD), 
              y = 0.95*max(y_MLE + y_SD)) %>% 
    left_join(p_table %>% select(species, p, label), by = "species")
  
  eiv_pwsc_plot[[i]] <- eiv_pwsc_data[[i]] %>% 
    filter(!is.na(y_MLE)) %>% 
    ggplot(aes(x_MLE, y_MLE, color = year)) + 
    geom_point(size = 1.5) +
    geom_errorbar(aes(ymin = y_MLE - y_SD, ymax = y_MLE + y_SD), size = 0.8) + 
    geom_errorbarh(aes(xmin = x_MLE - x_SD, xmax = x_MLE + x_SD), size = 0.8) + 
    geom_line(aes(x, fit), data = fit_i, inherit.aes = FALSE) + 
    geom_ribbon(aes(x, ymin = lower, ymax = upper), data = fit_i, inherit.aes = FALSE, 
                alpha = 0.3) +
    geom_text(aes(x, y, label = label), data = label_data %>% filter(p >= 0.05), 
              inherit.aes = FALSE, size = 4, color = "grey40",
              hjust = 1) +
    geom_text(aes(x, y, label = label), data = label_data %>% filter(p < 0.05), 
              inherit.aes = FALSE, size = 4, color = "black", fontface = "bold",
              hjust = 1) +
    facet_wrap(~species, scales = "free", nrow = 1) + 
    labs(x = names(overlap_metrics[i])) + 
    theme_bw() +
    theme(strip.background = element_blank(), 
          strip.text = element_text(size = 12, face = "bold"),
          plot.margin = plot_margins,
          axis.title.y = element_blank(), 
          legend.position = "none") + 
    scale_color_viridis_c(limits = range(preds_join$year))
  
  ## remove titles from all panels except the first
  if(i > 1) {
    eiv_pwsc_plot[[i]] <- eiv_pwsc_plot[[i]] + 
      theme(strip.text = element_blank())
  }
  
  eiv_pwsc_plot[[i]] <- ggplotGrob(eiv_pwsc_plot[[i]])
}

## Object for shared y axis label
pwsc_ylab <- grid::textGrob(expression(paste("relative predation ", "(", "kg kg"^{-1}, ")")), gp = grid::gpar(fontsize = 12), rot = 90)

## Construct full plot
eiv_pwsc_plot <-cowplot::plot_grid(
  pwsc_ylab, 
  cowplot::plot_grid(
    cowplot::plot_grid(plotlist = eiv_pwsc_plot, align = "hv", ncol = 1, rel_heights = c(1.05, rep(1, 5))), 
    preds_legend, rel_heights = c(10, 1), ncol = 1
  ),
  ncol = 2, rel_widths = c(0.05, 0.95)
) + theme(panel.border = element_blank())

## Save overlap vs standardized PESC plots
ggsave(paste(save_dir, "PWSC_overlap.png", sep = "/"), eiv_pwsc_plot, 
       height = 8, width = 10, units = "in", dpi = 500)
