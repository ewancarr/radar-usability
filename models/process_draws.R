# Title:        Process posterior draws from a/b/c-path models
# Author:       Ewan Carr
# Started:      2022-01-13

renv::activate()
library(tidyverse)
library(ggh4x)
library(here)
library(tidybayes)
library(broom)
library(broom.mixed)
library(brms)
library(emmeans)
library(tidybayes)
library(patchwork)
load(here("data", "clean", "sel.Rdata"), verbose = TRUE)
load(here("outputs", "posteriors.Rdata"), verbose = TRUE)

extract_draws <- function(list_of_brmsfit) {
  post <- map_dfr(list_of_brmsfit, function(x) {
                  f <- spread_draws(x, 
                                    `b_gad_1sd|b_ids_1sd|b_wsas_1sd|b_psurev_1sd|b_tamrev_1sd|b_ease_1sd|b_useful_1sd`, regex = TRUE)
                  if (ncol(f) == 4) {
                      names(f) <- c(".chain", ".iteration", ".draw", ".value")
                  } else if (ncol(f) == 5) {
                      names(f) <- c(".chain", ".iteration", ".draw", ".value", ".value_zi")
                  }
                  return(f) 
                  }, .id = "model")
  return(post)
}

map(list(a = a_path), extract_draws)
map(list(b = b_path), extract_draws)
map(list(c = c_path), extract_draws)

p <- map(list(a = a_path, b = b_path, c = c_path), extract_draws) %>%
       map(~ mutate(.x, model = str_replace_all(model, "1sd|l3|of", ""))) %>%
       map_dfr(~ separate(.x, model, c("left", "right")), .id = "path") %>%
       group_by(path, left, right) %>%
       summarise(median_qi(.value)) %>%
       mutate(cell = str_glue("{sprintf('%.2f', y)} [{sprintf('%.2f', ymin)}, {sprintf('%.2f', ymax)}]")) %>%
       select(path, left, right, cell) 

reshape_table <- function(tab, p = "a") { 
  tab <- filter(tab, path == p) %>% 
    spread(right, cell) %>%
    ungroup() %>% 
    select(-path)
}

###############################################################################
####                                                                      #####
####                     Models for 'FitBit wear time'                    #####
####                                                                      #####
###############################################################################

# Gather together models where 'wear_time' is the outcome
wear_time <- list(fit = list(b_path[[5]], b_path[[6]], b_path[[7]], b_path[[8]],
                             c_path[[1]], c_path[[2]], c_path[[3]]),
                  x = c(# b-paths
                        "psurev_1sd", "tamrev_1sd", "ease_1sd", "useful_1sd",
                        # c-paths
                        "ids_1sd", "gad_1sd", "wsas_1sd"),
                  label = c(rep("b_path", 4), rep("c_path", 3)))

# Extract model predictions
extract_effect <- function(fit, param, diff = c(-0.5, 0.5)) {
  at_spec <- list(diff); names(at_spec) <- param
  emmeans(fit,
          param,
          at = at_spec,
          epred = TRUE) %>%
  gather_emmeans_draws() %>%
  set_names(c("exposure", ".chain", ".iteration", ".draw", ".value")) 
}

post_wt <- pmap_dfr(wear_time, function(fit, x, label) {
                      comparison <- seq(-1, 1, 0.5)
                      extract_effect(fit, x, diff = comparison) %>%
                        mutate(path = label, x = x)
          })

integer_breaks <- function(x) seq(ceiling(x[1]), floor(x[2]), by = 1)

# Figure
plot_wt <- post_wt %>%
  ggplot(aes(x = exposure,
             y = .value,
             color = ordered(x))) +
  stat_lineribbon(alpha = 0.5, size = 0.75) +
  facet_nested_wrap(~ path + x,
                    scale = "free_x",
                    ncol = 4) +
  scale_fill_brewer(palette = "Greys") +
  scale_color_brewer(type = "qual") +
  theme_minimal(base_family = "Times New Roman") + 
  theme(legend.position = "none") +
  scale_x_continuous(breaks = integer_breaks) +
  lims(y = c(0, 1)) +
  labs(y = "Proportion FitBit wear time")

ggsave(plot_wt,
       filename = here("writing", "figures", "wear_time.png"),
       dev = "png",
       width = 6,
       height = 5,
       dpi = 300)


# Table
post_wt %>%
  group_by(path, x) %>%
  mutate(mid = median(exposure),
         lo = mid - 0.5,
         hi = mid + 0.5) %>%
  filter(exposure == lo | exposure == hi) %>%
  group_by(path, x) %>%
  mutate(hilo = case_when(exposure == min(exposure) ~ FALSE,
                          exposure == max(exposure) ~ TRUE)) %>%
  ungroup() %>%
  select(.draw, hilo, .value, path, x) %>%
  spread(hilo, .value) %>%
  mutate(diff = `TRUE` - `FALSE`) %>%
  group_by(path, x) %>%
  summarise(diff = median_qi(diff))


###############################################################################
####                                                                      #####
####                  Models for PHQ-8 completion counts                  #####
####                                                                      #####
###############################################################################

total_phq8 <- list(fit = list(b_path[[1]], b_path[[2]], 
                              b_path[[3]], b_path[[4]],
                              c_path[[4]], c_path[[5]],
                              c_path[[6]]),
                  x = c(# b-paths
                        "psurev_1sd", "tamrev_1sd", "ease_1sd", "useful_1sd",
                        # c-paths
                        "ids_1sd", "gad_1sd", "wsas_1sd"),
                  label = c(rep("b_path", 4), rep("c_path", 3)),
                  xr = as.logical(c(1, 1, 0, 0, 1, 1, 1)))

get_draws <- function(model,
                      var,
                      data,
                      use_1sd = TRUE, 
                      n_steps = 10) {
  gm <- function(x) { median(x, na.rm = TRUE) }
  xr <- seq(-1, 1, 0.5)
  at_medians <- data.frame(age_1sd = gm(data$age_1sd),
                           prevwear = gm(data$prevwear),
                           psurev_1sd = gm(data$psurev_1sd),
                           tamrev_1sd = gm(data$tamrev_1sd),
                           ease_1sd = gm(data$ease_1sd),
                           useful_1sd = gm(data$useful_1sd),
                           wsas = gm(data$wsas),
                           inwork = gm(data$inwork),
                           comorf = "0",
                           gad_1sd = gm(data$gad_1sd),
                           wsas_1sd = gm(data$wsas_1sd),
                           ids_1sd = gm(data$ids_1sd),
                           edyrs_1sd = gm(data$edyrs_1sd)) %>%
    select(-any_of(var))
  nd <- crossing({{var}} := xr, at_medians)
  draws <- epred_draws(model,
                       newdata = nd,
                       re_formula = NA,
                       summary = FALSE) 
  keep <- c(names(draws)[1], tail(names(draws), 3))
  draws <- draws[, keep]
  names(draws) <- c("exposure", ".draw", ".category", ".epred")
  draws %>%
    ungroup() %>%
    mutate(pk_count = .epred * as.double(.category)) %>%
    group_by(.draw, exposure) %>% 
    summarise(mean_count = sum(pk_count))
}


post_phq8 <- pmap_dfr(total_phq8, 
     function(fit, x, label, xr) {
       d <- get_draws(model = fit, 
                      var = x, 
                      data = sel,
                      use_1sd = xr)
       d$model <- label
       d$label <- x
       return(d) 
     }) %>%
  ungroup()

plot_phq8 <- post_phq8 %>%
  ggplot(aes(x = exposure,
             y = mean_count,
             color = label)) +
  stat_lineribbon(size = 1,
                  point_interval = "median_qi",
                  .width = c(0.5, 0.8, 0.95)) +
  theme_minimal(base_family = "Times New Roman") +
  scale_fill_brewer(palette = "Greys") +
  facet_nested_wrap(model ~ label,
             scales = "free_x", 
             ncol = 4) +
  theme(legend.position = "none") +
  scale_y_continuous(limits = c(0, 7), 
                     breaks = 0:7,
                     labels = ~ if_else(.x < 7,
                                        as.character(.x),
                                        "7+")) +
  scale_x_continuous(breaks = integer_breaks) +
  labs(y = "Predicted number of PHQ-8 questionnaires completed")

ggsave(plot_phq8,
       filename = here("writing", "figures", "phq8.png"),
       dev = "png",
       width = 6,
       height = 5,
       dpi = 300)


###############################################################################
####                                                                      #####
####                         Models for mediators                         #####
####                                                                      #####
###############################################################################

med <- list(fit = a_path,
            x = rep(c("wsas_1sd", "ids_1sd", "gad_1sd"), each = 4),
            y = rep(c("psu", "tam", "useful", "ease"), 3))

post_med <- pmap_dfr(med, function(fit, x, y) {
                       extract_effect(fit, x, seq(-1, 1, 0.5)) %>%
                         mutate(y = y,
                                x = x)
            }) %>%
  mutate(y_label = case_when(
            y == "ease" ~ "Perceived ease of use (0-7)",
            y == "useful" ~ "Perceived usefulness (0-7)",
            y == "tam" ~ "TAM-FF total (0-112; reversed)",
            y == "psu" ~ "PSSUQ total, (0-114; reversed)")) 


make_mediator_plot <- function(data) {
  data %>%
    ggplot(aes(x = exposure,
               y = .value,
               color = x)) +
    stat_lineribbon(size = 1,
                    point_interval = "median_qi",
                    .width = c(0.5, 0.8, 0.95)) +
    theme_minimal(base_family = "Times New Roman") +
    scale_fill_brewer(palette = "Greys") +
    facet_nested_wrap(y_label ~ x,
                      scales = "free_x", 
                      ncol = 6) +
    theme(legend.position = "none") +
    scale_color_brewer(type = "qual")
}

top_plot <- post_med %>%
  filter(y %in% c("ease", "useful")) %>%
  make_mediator_plot() +
  scale_x_continuous(limits = c(-1, 1), breaks = integer_breaks) +
  scale_y_continuous(limits = c(0, 7), breaks = integer_breaks) 

bottom_plot <- post_med %>%
  filter(y %in% c("tam", "psu")) %>%
  make_mediator_plot() +
  scale_x_continuous(limits = c(-1, 1), breaks = integer_breaks) +
  scale_y_continuous(limits = c(0, 114)) 

plot_med <- top_plot / bottom_plot

ggsave(plot_med,
       filename = here("writing", "figures", "mediators.png"),
       dev = "png",
       width = 6,
       height = 5,
       dpi = 300)

###############################################################################
####                                                                      #####
####         Make table showing 1SD effect sizes for all outcomes         #####
####                                                                      #####
###############################################################################

make_cell <- function(est, lo, hi) {
  return(paste0(str_glue("{sprintf('%.2f', est)} "),
                str_glue("[{sprintf('%.2f', lo)}, "),
                str_glue("{sprintf('%.2f', hi)}]")))
}

# Assumble results for a-paths
tab_med <- post_med %>%
  filter(exposure %in% c(-0.5, 0.5)) %>%
  select(.draw, .value, y, y_label, exposure, x)  %>%
  spread(exposure, .value) %>%
  mutate(contrast = `0.5` - `-0.5`) %>%
  group_by(x, y, y_label) %>%
  median_qi(contrast) %>%
  mutate(cell = make_cell(contrast, .lower, .upper)) %>%
  select(y, x, cell) %>%
  mutate(id = "a_path")

a <- median(sel$ease_1sd, na.rm = TRUE)
b <- median(sel$useful_1sd, na.rm = TRUE)

# Assemble results for "Total PHQ-8 completions"
tab_phq8 <- post_phq8 %>%
  mutate(id = case_when((str_ends(label, "_1sd") & exposure == -0.5) ~ FALSE,
                        (str_ends(label, "_1sd") & exposure == 0.5) ~ TRUE)) %>%
  drop_na(id) %>%
  select(.draw, id, mean_count, model, label) %>%
  spread(id, mean_count) %>%
  mutate(contrast = `TRUE` - `FALSE`) %>%
  group_by(model, label) %>%
  median_qi(contrast) %>% 
  mutate(id = "phq8",
         cell = make_cell(contrast, .lower, .upper)) %>%
  select(id, path = model, x = label, cell) 


# Assemble results for "FitBit wear time"
tab_wt <- post_wt %>%
  mutate(id = case_when((str_ends(x, "_1sd") & exposure == -0.5) ~ FALSE,
                        (str_ends(x, "_1sd") & exposure == 0.5) ~ TRUE)) %>%
  drop_na(id) %>%
  ungroup() %>%
  select(.draw, id, .value, path, x) %>%
  spread(id, .value) %>%
  mutate(contrast = `TRUE` - `FALSE`) %>%
  group_by(path, x) %>%
  median_qi(contrast) %>%
  mutate(id = "wt",
         cell = make_cell(contrast, .lower, .upper)) %>%
  select(id, path, x, cell) 

tab_med %>% spread(y, cell) %>%
  write_csv("~/med.csv")

tab <- full_join(tab_phq8, tab_wt) %>%
  pivot_wider(names_from = "id",
              values_from = "cell")

write_csv(tab, "~/tab.csv")
