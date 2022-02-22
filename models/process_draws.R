# Title:        Process posterior draws from H1/H2/H3 models
# Author:       Ewan Carr
# Started:      2022-01-13

renv::activate()
library(tidyverse)
library(ggh4x)
library(here)
library(tidybayes)
library(broom)
library(broom.mixed)
library(rstanarm)
library(brms)
library(emmeans)
library(tidybayes)
library(patchwork)
library(RColorBrewer)
library(colorspace)
load(here("data", "clean", "sel.Rdata"), verbose = TRUE)
ff <- "Arial"

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

load(here("outputs", "fit_h2.Rdata"), verbose = TRUE)
load(here("outputs", "fit_h3.Rdata"), verbose = TRUE)

# Gather together models where 'wear_time' is the outcome
wear_time <- list(fit = list(fit_h2[[5]], fit_h2[[6]], fit_h2[[7]], fit_h2[[8]],
                             fit_h3[[1]], fit_h3[[2]], fit_h3[[3]]),
                  x = c(# H2
                        "psurev_1sd", "tamrev_1sd", "ease_1sd", "useful_1sd",
                        # H3
                        "ids_1sd", "gad_1sd", "wsas_1sd"),
                  label = c(rep("h2", 4), rep("h3", 3)))

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

labels <- tribble(~x,           ~label,   ~sub,
                  "psurev_1sd", "PSSUQ",  "Total score",
                  "tamrev_1sd", "TAM-FF", "Total score",
                  "useful_1sd", "TAM-FF", "Perceived usefulness",
                  "ease_1sd",   "TAM-FF", "Perceived ease of use",
                  "gad_1sd",    "GAD-7",  "Total score",
                  "ids_1sd",    "IDS-SR", "Total score",
                  "wsas_1sd",   "WSAS",   "Total score")
labels$label <- factor(labels$label, levels = c("PSSUQ", "TAM-FF", "IDS-SR", "GAD-7", "WSAS"))
labels$sub <- factor(labels$sub, levels = c("Total score", "Perceived ease of use", "Perceived usefulness"))
plot_data <- left_join(post_wt, labels)

pc <- c(tail(brewer.pal(6, "Blues"), 2), 
        tail(brewer.pal(6, "Purples"), 2), 
        brewer.pal(3, "Set2"))
col <- c("PSSUQ"  = pc[2],
         "TAM-FF" = pc[1],
         "TAM-FF" = pc[3],
         "TAM-FF" = pc[4],
         "GAD-7"  = pc[5],
         "IDS-SR" = pc[6],
         "WSAS"   = pc[7])

plot_wt <- plot_data %>%
  ggplot() +
  aes(x = exposure,
      fill = label,
      color = label,
      y = .value) +
  stat_lineribbon(.width = c(0.95),
                  point_interval = "median_qi") +
  facet_nested_wrap(~ label + sub, ncol = 4) +
  theme_minimal(base_family = ff) + 
  theme(legend.position = "none") +
  scale_x_continuous(breaks = integer_breaks) +
  lims(y = c(0, 1)) +
  labs(y = "Proportion of time wearing FitBit, next 3 months",
       x = "Standard deviations difference in exposure") +
  theme(strip.background = element_rect(fill = "gray85", line = 0)) +
  scale_color_manual(values = col) +
  scale_fill_manual(values = map(col, lighten, 0.8))

ggsave(plot_wt,
       filename = here("writing", "figures", "figs3_weartime.png"),
       dev = "png",
       width = 8,
       height = 6,
       dpi = 300)

###############################################################################
####                                                                      #####
####                  Models for PHQ-8 completion counts                  #####
####                                                                      #####
###############################################################################

n_phq8 <- list(fit = list(fit_h2[[1]], fit_h2[[2]], 
                          fit_h2[[3]], fit_h2[[4]],
                          fit_h3[[4]], fit_h3[[5]],
                          fit_h3[[6]]),
               x = c(# H2
                     "psurev_1sd", "tamrev_1sd", "ease_1sd", "useful_1sd",
                     # H3
                     "ids_1sd", "gad_1sd", "wsas_1sd"),
               label = c(rep("h2", 4), rep("h3", 3)),
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
                           n_lte_f = "0",
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


post_phq8 <- pmap_dfr(n_phq8, 
     function(fit, x, label, xr) {
       d <- get_draws(model = fit, 
                      var = x, 
                      data = sel,
                      use_1sd = xr)
       d$model <- label
       d$x <- x
       return(d) 
     }) %>%
  ungroup()

plot_data <- left_join(post_phq8, labels)

plot_phq8 <- post_phq8 %>%
  left_join(labels) %>%
  ggplot() +
  aes(x = exposure,
      y = mean_count,
      fill = label,
      color = label) +
  stat_lineribbon(.width = 0.95,
                  point_interval = "median_qi") +
  facet_nested_wrap(~ label + sub, ncol = 4) +
  theme_minimal(base_family = ff) + 
  theme(legend.position = "none") +
  scale_x_continuous(breaks = integer_breaks) +
  scale_y_continuous(limits = c(0, 6), 
                     breaks = 0:6) +
  labs(y = "Number of PHQ-8 questionnaires completed",
       x = "Standard deviations difference in exposure") +
  theme(strip.background = element_rect(fill = "gray85", line = 0)) +
  scale_color_manual(values = col) +
  scale_fill_manual(values = map(col, lighten, 0.8))

ggsave(plot_phq8,
       filename = here("writing", "figures", "figs2_phq8.png"),
       dev = "png",
       width = 8,
       height = 6,
       dpi = 300)

###############################################################################
####                                                                      #####
####                 Models for usability and acceptability               #####
####                                                                      #####
###############################################################################

load(here("outputs", "fit_h1.Rdata"), verbose = TRUE)

usab <- list(fit = fit_h1,
           x = rep(c("wsas_1sd", "ids_1sd", "gad_1sd"), each = 4),
           y = rep(c("psu", "tam", "useful", "ease"), 3))

post_usab <- pmap_dfr(usab, function(fit, x, y) {
                       extract_effect(fit, x, seq(-1, 1, 0.5)) %>%
                         mutate(y = y,
                                x = x)
            }) %>%
  mutate(y_label = case_when(
            y == "ease" ~ "Perceived ease of use (0-7)",
            y == "useful" ~ "Perceived usefulness (0-7)",
            y == "tam" ~ "TAM-FF total (0-112; reversed)",
            y == "psu" ~ "PSSUQ total, (0-114; reversed)")) 

plot_data <- left_join(post_usab, labels)

make_mediator_plot <- function(data, colors) {
  data %>%
    ggplot(aes(x = exposure,
               y = .value,
               color = label,
               fill = label)) +
    stat_lineribbon(point_interval = "median_qi",
                    .width = 0.95) +
    theme_minimal(base_family = ff) +
    facet_nested_wrap(y_label ~ label,
                      scales = "free_x", 
                      ncol = 6) +
    theme(legend.position = "none",
          strip.background = element_rect(fill = "gray85", line = 0)) +
    scale_color_manual(values = colors) +
    scale_fill_manual(values = map(colors, lighten, 0.8))
}

bottom_plot <- plot_data %>%
  filter(y %in% c("ease", "useful")) %>%
  make_mediator_plot(col) +
  scale_x_continuous(limits = c(-1, 1), breaks = integer_breaks) +
  scale_y_continuous(limits = c(0, 7), breaks = integer_breaks) +
  theme(axis.title.x = element_text(margin = margin(0.4, 0, 0, 0, "cm"))) +
  labs(y = "Mean score",
       x = "Standard deviations difference in exposure")

top_plot <- plot_data %>%
  filter(y %in% c("tam", "psu")) %>%
  make_mediator_plot(col) +
  scale_x_continuous(limits = c(-1, 1), breaks = integer_breaks) +
  scale_y_continuous(limits = c(0, 114)) +
  labs(y = "Mean score") +
  theme(axis.title.x = element_blank())

plot_usab <- top_plot / bottom_plot

ggsave(plot_usab,
       filename = here("writing", "figures", "figs1_usab.png"),
       dev = "png",
       width = 8,
       height = 6,
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

# Assumble results for "usability and acceptability" outcomes
tab_usab <- post_usab %>%
  filter(exposure %in% c(-0.5, 0.5)) %>%
  select(.draw, .value, y, y_label, exposure, x)  %>%
  spread(exposure, .value) %>%
  mutate(contrast = `0.5` - `-0.5`) %>%
  group_by(x, y, y_label) %>%
  median_qi(contrast) %>%
  mutate(cell = make_cell(contrast, .lower, .upper)) %>%
  select(y, x, cell) %>%
  mutate(id = "usab")

a <- median(sel$ease_1sd, na.rm = TRUE)
b <- median(sel$useful_1sd, na.rm = TRUE)

# Assemble results for "Total PHQ-8 completions"
tab_phq8 <- post_phq8 %>%
  mutate(id = case_when((str_ends(x, "_1sd") & exposure == -0.5) ~ FALSE,
                        (str_ends(x, "_1sd") & exposure == 0.5) ~ TRUE)) %>%
  drop_na(id) %>%
  select(.draw, id, mean_count, model, x) %>%
  spread(id, mean_count) %>%
  mutate(contrast = `TRUE` - `FALSE`) %>%
  group_by(model, x) %>%
  median_qi(contrast) %>% 
  mutate(id = "phq8",
         cell = make_cell(contrast, .lower, .upper)) %>%
  select(id, path = model, x, cell) 


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

tab_usab %>% spread(y, cell) %>%
  write_csv("~/usab.csv")

tab <- full_join(tab_phq8, tab_wt) %>%
  pivot_wider(names_from = "id",
              values_from = "cell")

write_csv(tab, "~/tab.csv")
