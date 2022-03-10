
###############################################################################
##' @title Supplemental information: summaries of predator diets in the EBS
##' @author Maurice Goodman
###############################################################################

## Read in and format data ----------------------------------------------------

library(here)
library(tidyverse)

pred_info <- list(
  sciname = c(
    "Atheresthes stomias", 
    "Gadus macrocephalus", 
    "Hippoglossus stenolepis", 
    "Gadus chalcogrammus"
  ), 
  common = c(
    "Arrowtooth Flounder", 
    "Pacific Cod", 
    "Pacific Halibut", 
    "Walleye Pollock"
  ), 
  min_size = c(30, 30, 50, 40)
)

## read in from directory with csvs for all four predators
files <- list.files(here("data", "stomach contents"), full.names = TRUE)
diet <- do.call("rbind", lapply(files, read.csv))
diet$Pred_name <- plyr::revalue(diet$Pred_name, setNames(pred_info$common, pred_info$sciname))

## minimum size for predators
diet$threshold <- pred_info$min_size[match(diet$Pred_name, pred_info$common)]

## Predator diets by frequency of occurrence in stomachs ----------------------

## Table containing one entry for each prey item in each unique stomach sample
## Requires adding zeroes for prey items which do not already have entries, 
## and combining rows for prey items which have multiple entries for same stomach
prey_pred_df <- diet %>% 
  filter(Year >= 1982 & Year <= 2015 & Pred_len >= threshold & Prey_Name != "Empty" & !is.na(Prey_cnt)) %>% 
  mutate(sample = paste(Pred_name, Year, Hauljoin, Pred_specn, sep = "-")) %>% 
  select(sample, Prey_Name, Prey_cnt) %>% 
  complete(sample, Prey_Name, fill = list(Prey_cnt = 0)) %>% 
  group_by(sample, Prey_Name) %>% 
  summarize(n = sum(Prey_cnt)) %>%
  ungroup()

## Aggregate stomach samples to obtain annual proportions
prey_freq_annual <- prey_pred_df %>% 
  separate(sample, into = c("pred_name", "year", "haul_join", "pred_specn"), sep = "-", remove = FALSE) %>%
  mutate(present = as.numeric(n > 0)) %>% 
  group_by(pred_name, Prey_Name, year) %>% 
  summarize(n_present = sum(present), 
            n_sample = length(unique(sample))) %>%
  ungroup() %>% group_by(pred_name) %>% 
  mutate(p = n_present/n_sample)

## Obtain overall proportions, ranks, and SD
prey_freq <- prey_freq_annual %>% 
  group_by(pred_name, Prey_Name) %>% 
  summarize(p_mean = weighted.mean(p, n_sample), 
            SD = sd(p)) %>% 
  ungroup() %>% group_by(pred_name) %>% 
  mutate(rank = dense_rank(desc(p_mean))) %>% 
  arrange(pred_name, rank)

## Violin plot of annual frequencies
freq_plot <- prey_freq_annual %>% 
  filter(Prey_Name == "Walleye pollock") %>% 
  ggplot(aes(pred_name, p)) +
  geom_violin(scale = "width", fill = "grey60", color = NA, alpha = 0.5) + 
  geom_boxplot(width = 0.1, size = 0.8, color = "black") + 
  theme_bw() +
  theme(panel.grid.minor = element_blank(), 
        panel.grid.major.x = element_blank(), 
        axis.title.x = element_blank(), 
        axis.text = element_text(color = "black"), 
        axis.ticks = element_line(color = "black"), 
        panel.border = element_rect(color = "black")) + 
  ylab("frequency in stomachs") + 
  coord_cartesian(clip = "off")

## Predator diets by biomass --------------------------------------------------

## Aggregate samples to obtain annual proportions
## don't need to add zeroes first, but do have to add them after summarizing
prey_mass_annual <- diet %>% 
  filter(Year >= 1982 & Year <= 2015 & Pred_len >= threshold & Prey_Name != "Empty") %>%
  group_by(Pred_name, Prey_Name, Year) %>% 
  summarize(weight = sum(Prey_twt)) %>% 
  ungroup() %>% 
  complete(nesting(Pred_name, Year), Prey_Name, fill = list(weight = 0)) %>% 
  group_by(Pred_name, Year) %>% 
  mutate(twt = sum(weight),
         p = weight/sum(weight)) %>% 
  ungroup()

## Obtain overall proportions, rank and SD
prey_mass <- prey_mass_annual %>% 
  group_by(Pred_name, Prey_Name) %>% 
  summarize(p_mean = weighted.mean(p, twt), 
            SD = sd(p)) %>% 
  ungroup() %>% group_by(Pred_name) %>% 
  mutate(rank = dense_rank(desc(p_mean))) %>% 
  arrange(Pred_name, rank)

## violin plot of annual proportions by biomass
weight_plot <- prey_mass_annual %>% 
  filter(Prey_Name == "Walleye pollock") %>% 
  ggplot(aes(Pred_name, p)) +
  geom_violin(scale = "width", fill = "grey60", color = NA, alpha = 0.5) + 
  geom_boxplot(width = 0.1, size = 0.8, color = "black") + 
  theme_bw() +
  theme(panel.grid.minor = element_blank(), 
        panel.grid.major.x = element_blank(), 
        axis.title.x = element_blank(), 
        axis.text = element_text(color = "black"), 
        axis.ticks = element_line(color = "black"), 
        panel.border = element_rect(color = "black")) + 
  ylab("proportion of\nstomach content weight") + 
  coord_cartesian(clip = "off")


## Combine plots --------------------------------------------------------------

full_plot <- cowplot::plot_grid(freq_plot, weight_plot, ncol = 1, align = "hv")

ggsave(here("model output", "stomach contents.png"), full_plot, height = 4, 
       width = 8, units = "in",  dpi = 500)

