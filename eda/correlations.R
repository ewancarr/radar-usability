# Title:        Calculate correlations for paper
# Author:       Ewan Carr
# Started:      2022-02-21

renv::activate()
library(tidyverse)
library(here)
library(gtsummary)
library(confintr)
library(patchwork)
library(furrr)
library(huxtable)
library(RColorBrewer)
library(boot)
library(survey)
library(jtools)
set.seed(42)
plan(multicore)
load(here("data", "clean", "merged.Rdata"), verbose = TRUE)
load(here("data", "clean", "selected.Rdata"), verbose = TRUE)
ff <- "Arial"

sel <- merged %>%
  filter(subject_id %in% picks) %>%
  mutate(t_label = factor(t, levels = c(0, 3, 12),
                          labels = c("0 months\n(Enrolment)",
                                     "3 months",
                                     "12 months")),
         across(c(ids, gad, wsas, psurev, tamrev, 
                  ease, useful, age, edyrs),
                scale,
                .names = "{.col}_1sd"),
         pid = as.numeric(as.factor(subject_id)))

print(length(unique(sel$subject_id)))

# Table 2: Correlations between X, M, Y ---------------------------------------
boot_corr <- function(data, ix) {
    d <- data[ix, ]
    # NOTE: using svycorr from jtools package to account for
    #       repeated measures per participant.
    d_clus <- svydesign(id = ~pid, data = d)
    results <- svycor(x ~ y, design = d_clus, na.rm = TRUE)$cors[2]
    return(c(r = results))
}

calculate_correlation <- function(x, y, data, reps = 1000) {
  x <- data[[x]]
  y <- data[[y]]
  pid = data$pid
  r <- cor.test(x, y)$estimate
  bs <- boot(data.frame(x = x, 
                        y = y,
                        pid = pid),
             boot_corr,
             R = reps)
  ci <- boot.ci(bs, type = "perc")
  return(c(r, ci$percent[c(4, 5)]))
}

x <- c("gad_1sd", "ids_1sd", "wsas_1sd",
       "psurev_1sd", "tamrev_1sd", "useful_1sd", "ease_1sd")
y <- c("psurev_1sd", "tamrev_1sd", "useful_1sd", "ease_1sd",
       "wear_time_l3", "n_phq")

opts <- cross2(x, y)
r <- future_map(opts, ~ calculate_correlation(.x[[1]], .x[[2]], sel),
                seed = TRUE)

results <- pmap_dfr(list(key = opts, r = r),
                    ~ c(measure = .x, result = .y))

# Select required comparisons [TODO: find better way of doing this]

keepers <- rbind(# a-paths
                 c("gad_1sd", "psurev_1sd"),
                 c("gad_1sd", "tamrev_1sd"),
                 c("gad_1sd", "useful_1sd"),
                 c("gad_1sd", "ease_1sd"),
                 c("ids_1sd", "psurev_1sd"),
                 c("ids_1sd", "tamrev_1sd"),
                 c("ids_1sd", "useful_1sd"),
                 c("ids_1sd", "ease_1sd"),
                 c("wsas_1sd", "psurev_1sd"),
                 c("wsas_1sd", "tamrev_1sd"),
                 c("wsas_1sd", "useful_1sd"),
                 c("wsas_1sd", "ease_1sd"),
                 # b-paths
                 c("psurev_1sd", "wear_time_l3"),
                 c("psurev_1sd", "n_phq"),
                 c("tamrev_1sd", "wear_time_l3"),
                 c("tamrev_1sd", "n_phq"),
                 c("ease_1sd", "wear_time_l3"),
                 c("ease_1sd", "n_phq"),
                 c("useful_1sd", "wear_time_l3"),
                 c("useful_1sd", "n_phq"),
                 # c-paths
                 c("gad_1sd", "wear_time_l3"),
                 c("gad_1sd", "n_phq"),
                 c("ids_1sd", "wear_time_l3"),
                 c("ids_1sd", "n_phq"),
                 c("wsas_1sd", "wear_time_l3"),
                 c("wsas_1sd", "n_phq")) %>%
  as.data.frame() %>%
  set_names(c("measure1", "measure2"))

results <- keepers %>% left_join(results)

results$measure1 <- factor(results$measure1,
                           levels = c("gad_1sd", "ids_1sd", "wsas_1sd",
                                      "psurev_1sd", "tamrev_1sd",
                                      "ease_1sd", "useful_1sd"),
                           labels = c("GAD-7", "IDS-SR", "WSAS",
                                      "PSSUQ\n(reversed)",
                                      "TAM-FF\n(reversed)",
                                      "Perceived\nease of use",
                                      "Perceived\nusefulness"))

results$measure2 <- factor(results$measure2,
                           levels = c("psurev_1sd", "tamrev_1sd",
                                      "ease_1sd", "useful_1sd",
                                      "wear_time_l3", "n_phq"),
                           labels = c("PSSUQ\n(reversed)",
                                      "TAM-FF\n(reversed)",
                                      "Perceived\nease of use",
                                      "Perceived\nusefulness",
                                      "FitBit\nwear time",
                                      "No. PHQ-8\nreturned"))

results$r <- sprintf("%.2f", results$result.cor)
results$ci <- paste0("[", sprintf("%.2f", results$result2), ", ",
                     sprintf("%.2f", results$result3), "]")

p_corr <- ggplot(results,
       aes(x = measure2, y = forcats::fct_rev(measure1))) +
  geom_tile(aes(fill = result.cor), color = "white", size = 1.5) +
  geom_text(aes(label = r),
            nudge_y = 0.1,
            family = ff) +
  geom_text(aes(label = ci),
            nudge_y = -0.1,
            size = 3,
            family = ff) +
  scale_fill_gradient2(low = "#e41a1c",
                       high = "#4daf4a",
                       mid = "white",
                       midpoint = 0,
                       limit = c(-1, 1),
                       space = "Lab") +
  theme_minimal(base_family = ff) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text = element_text(size = 11)) +
  labs(fill = "Pearson\ncorrelation")
       # title = "Pearson correlations between predictors and outcomes (n = 551)",
       # subtitle = "[95% bootstrap percentile confidence interval]")

ggsave(p_corr,
       filename = here("writing", "figures", "p_corr.png"),
       dev = "png",
       width = 7,
       height = 5,
       dpi = 300)
