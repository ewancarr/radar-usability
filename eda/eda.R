# Title:        Prepare descriptives (Table 1 and 2) for RADAR-MDD paper
# Author:       Ewan Carr
# Started:      2022-01-11

renv::activate()
library(tidyverse)
library(here)
library(gtsummary)
library(confintr)
library(patchwork)
library(furrr)
library(huxtable)
library(RColorBrewer)
plan(multicore)
load(here("data", "clean", "merged.Rdata"), verbose = TRUE)
load(here("data", "clean", "sel.Rdata"), verbose = TRUE)
ff <- "Arial"

sel <- left_join(sel,
                 select(merged, subject_id, t, male), 
                 by = c("subject_id", "t")) %>%
  mutate(t_label = case_when(t == 0 ~ "0 months\n(Baseline)",
                             t == 3 ~ "3 months",
                             t == 12 ~ "12 months"),
         t_label = factor(t_label))

# Table 1: Descriptives of X, M, Y at 3 and 12 months -------------------------

sel_fu <- filter(sel, t > 0)

# Demographic characteristics
demo <- sel %>%
  filter(t == 0) %>%
  select(age, male, comorf, edyrs, inwork, t_label) %>%
  drop_na() %>%
  tbl_summary(by = t_label)


x <- sel %>%
  select(t_label, ids, gad, wsas) %>%
  tbl_summary(by = t_label,
              missing_text = "Missing",
              label = list(ids ~ "Inventory of Depressive Symptomatology (IDS)",
                           gad ~ "Generalized Anxiety Disorder-7",
                           wsas ~ "Work and Social Adjustment Scale (WSAS)"))

m <- sel %>%
  select(t_label, psurev, tamrev, useful, ease) %>%
  drop_na() %>%
  tbl_summary(by = t_label,
              missing_text = "Missing",
              label = list(psurev ~ "Post-Study System Usability Questionnaire (PSSUQ; reversed)",
                           tamrev ~ " Technology Acceptance Model (TAM; reversed)",
                           useful ~ "Perceived usefulness",
                           ease ~ "Perceived ease of use"))

y <- sel %>%
  select(t_label, wear_time_l3, total_phq8_l3) %>%
  drop_na() %>%
  tbl_summary(by = t_label,
              missing_text = "Missing",
              label = list(
                  wear_time_l3 ~ "Proportion of time wearing FitBit, next 3 months",
                  total_phq8_l3 ~ "Number of PHQ-8 questionnaire completions, next 3 months"),
              type = list(total_phq8_l3 ~ "continuous"))

table1 <- tbl_stack(list(demo, x, m, y),
                    group_header = c("Demographic characteristics",
                                     "Exposures",
                                     "Mediators",
                                     "Outcomes")) %>%
   modify_header(list(label ~"")) %>%
   as_gt() %>%
   gt::tab_style(style = gt::cell_text(weight = "bold"),
                 locations = gt::cells_row_groups(groups = everything()))

# table1 %>% as_huxtable() %>% quick_xlsx()

# Calculate range for each variable

sel %>%
  select(age, male, comorf, edyrs, inwork, t_label,
         ids, gad, wsas,
         psurev, tamrev, useful, ease,
         wear_time_l3, total_phq8_l3) %>%
  drop_na() %>%
  summarise(across(everything(), function(x) { str_glue("{min(as.numeric(x), na.rm = TRUE)}, {max(as.numeric(x), na.rm = TRUE)}")})) %>%
  pivot_longer(everything(),
               names_to = "var",
               values_to = "range")


# Figure 2: Distributions of mediators and outcomes ---------------------------

# Define colours
pc <- c(tail(brewer.pal(6, "Blues"), 2), 
        tail(brewer.pal(6, "Purples"), 2), 
        tail(brewer.pal(8, "Greens"), 2))

# Make plot for mediators
labels <- data.frame(measure = c("psurev",
                                          "tamrev",
                                          "useful",
                                          "ease",
                                          "total_phq8_l3_of7",
                                          "wear_time_l3"),
                              label = factor(1:6,
                                             labels = c("Post-Study System Usability\n(PSSUQ; reversed)",
                                                        "Technology Acceptance\nModel (TAM; reversed)",
                                                        "Perceived usefulness",
                                                        "Perceived ease of use",
                                                        "No. PHQ-8\nreturned",
                                                        "FitBit\nwear time")))
fu <- sel %>%
  filter(t %in% c(3, 12))

theming <- 
  theme_minimal(base_family = ff) +
  theme(axis.title = element_blank()) 

p_t1 <- ggplot(fu, aes(x = tamrev)) +
  stat_bin(fill = pc[1],
           binwidth = 10) +
  labs(y = "Count",
       title = "Technology Acceptance\nModel (TAM; reversed)") +
  theming +
  theme(axis.title.y = element_text()) +
  scale_x_continuous(breaks = seq(0, 120, 20)) +
  coord_cartesian(ylim = c(0, 250))

p_t2 <- fu %>% drop_na(psurev) %>%
  ggplot(aes(x = psurev)) +
  geom_histogram(fill = pc[2],
                 binwidth = 10) +
  labs(title = "Post-Study System Usability\n(PSSUQ; reversed)") +
  theming +
  scale_x_continuous(limits = c(0, 120),
                     breaks = seq(0, 120, 20)) +
  coord_cartesian(ylim = c(0, 250))

p_t3 <- ggplot(fu, aes(x = round(useful))) +
  stat_count(fill = pc[3],
             binwidth = 1) +
  theming +
  labs(title = "Perceived\nusefulness") +
  coord_cartesian(ylim = c(0, 250)) +
  scale_x_continuous(breaks = 0:7)

p_t4 <- ggplot(fu, aes(x = round(pc_ease_r))) +
  stat_count(fill = pc[4],
             binwidth = 1) +
  labs(title = "Perceived\nease of use") +
  theming +
  coord_cartesian(ylim = c(0, 250)) +
  scale_x_continuous(breaks = 0:7)

# Make plot for wear time
p_wt <- fu %>%
  select(wear_time_l3) %>%
  drop_na() %>%
  ggplot(aes(x = wear_time_l3)) +
  geom_histogram(fill = pc[5]) +
  scale_y_continuous(limits = c(0, 300)) +
  labs(title = "Proportion of time wearing FitBit, next three months",
       y = "Count") +
  theme_minimal(base_family = ff) +
  theme(axis.title.x = element_blank())

# Make plot for PHQ-8 completions
p_phq8 <- fu %>%
  select(total_phq8_l3_of7) %>%
  drop_na() %>%
  ggplot(aes(x = total_phq8_l3_of7)) +
  stat_count(fill = pc[6]) +
  scale_x_discrete(breaks = 0:7,
                   labels = ~ if_else(.x < 7, .x, "7+")) +
  scale_y_continuous(limits = c(0, 300)) +
  labs(title = "Total PHQ-8 completions, next three months") +
  theme_minimal(base_family = ff) +
  theme(axis.title = element_blank())

p_dist <- (p_t1 + p_t2 + p_t3 + p_t4 + plot_layout(nrow = 1)) / (p_wt + p_phq8)
ggsave(p_dist, file = here("writing", "figures", "distributions.png"),
       dev = "png",
       width = 10,
       height = 6,
       dpi = 300)



# Table 2: Correlations between X, M, Y ---------------------------------------

library(boot)
library(survey)
library(jtools)
set.seed(42)

boot_corr <- function(data, ix) {
    d <- data[ix, ]
    # NOTE: using svycorr from jtools package to account for repeated measures per participant
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
       "wear_time_l3", "total_phq8_l3")

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
                 c("psurev_1sd", "total_phq8_l3"),
                 c("tamrev_1sd", "wear_time_l3"),
                 c("tamrev_1sd", "total_phq8_l3"),
                 c("ease_1sd", "wear_time_l3"),
                 c("ease_1sd", "total_phq8_l3"),
                 c("useful_1sd", "wear_time_l3"),
                 c("useful_1sd", "total_phq8_l3"),
                 # c-paths
                 c("gad_1sd", "wear_time_l3"),
                 c("gad_1sd", "total_phq8_l3"),
                 c("ids_1sd", "wear_time_l3"),
                 c("ids_1sd", "total_phq8_l3"),
                 c("wsas_1sd", "wear_time_l3"),
                 c("wsas_1sd", "total_phq8_l3")) %>%
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
                                      "wear_time_l3", "total_phq8_l3"),
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

