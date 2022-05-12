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
library(janitor)
library(naniar)
library(RColorBrewer)
plan(multicore)
load(here("data", "clean", "merged.Rdata"), verbose = TRUE)
load(here("data", "clean", "selected.Rdata"), verbose = TRUE)
ff <- "Arial"

sel <- merged %>%
  filter(subject_id %in% picks) %>%
  mutate(t_label = factor(t, levels = c(0, 3, 12),
                          labels = c("0 months\n(Enrolment)",
                                     "3 months",
                                     "12 months"))) %>%
  ungroup()

# Table 1: Descriptives of X, M, Y at 3 and 12 months -------------------------

stat <- all_continuous() ~ c("{median} ({p25}, {p75}) [{min}, {max}]",
                             "{N_miss}")

tab_demo <- sel %>%
  select(age, male, comorf, edyrs, inwork, n_lte_f, t_label) %>%
  tbl_summary(by = t_label,
              statistic = all_continuous() ~ c("{median} ({p25}, {p75}) [{min}, {max}]"),
              type = list(all_continuous() ~ "continuous",
                          n_lte_f ~ "categorical"),
              missing = "no",
              label = list(age = "Age",
                           male = "Male gender",
                           comorf = "Number of comorbid conditions",
                           edyrs = "Years of education",
                           inwork = "Currently working",
                           n_lte_f = "Lifetime events, last 3 months"))

tab_clinical <- sel %>%
  select(t_label, ids, gad, wsas) %>%
  tbl_summary(by = t_label,
              statistic = stat,
              type = all_continuous() ~ "continuous2",
              missing = "no",
              label = list(ids ~ "Inventory of Depressive Symptomatology (IDS)",
                           gad ~ "Generalized Anxiety Disorder-7",
                           wsas ~ "Work and Social Adjustment Scale (WSAS)"))

tab_usability <- sel %>%
  select(t_label, psurev, tamrev, useful, ease) %>%
  tbl_summary(by = t_label,
              statistic = stat,
              type = all_continuous() ~ "continuous2",
              missing = "no",
              label = list(psurev ~ "Post-Study System Usability Questionnaire (PSSUQ; reversed)",
                           tamrev ~ " Technology Acceptance Model (TAM; reversed)",
                           useful ~ "Perceived usefulness",
                           ease ~ "Perceived ease of use"))

tab_usage <- sel %>%
  select(t_label, wear_time_l3, n_phq_f) %>%
  drop_na() %>%
  tbl_summary(by = t_label,
              statistic = stat,
              type = list(all_continuous() ~ "continuous2",
                          n_phq_f ~ "categorical"),
              missing = "no",
              label = list(wear_time_l3 ~ "Proportion of time wearing FitBit, next 3 months",
                           n_phq_f ~ "Number of PHQ-8 questionnaire completions, next 3 months"))

table1 <- tbl_stack(list(tab_demo,
                         tab_clinical,
                         tab_usability,
                         tab_usage),
                    group_header = c("Demographic characteristics",
                                     "Clinical variables",
                                     "Usability and acceptability",
                                     "Usage")) %>%
   modify_header(list(label ~"")) %>%
   as_gt() %>%
   gt::tab_style(style = gt::cell_text(weight = "bold"),
                 locations = gt::cells_row_groups(groups = everything()))

table1

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
                                 "n_phq8_f",
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

length(unique(fu$subject_id))

theming <- 
  theme_minimal(base_family = ff) +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_text(margin = margin(0.5, 0, 0, 0, "cm")),
        plot.title = element_text(size = 12,
                                  margin = margin(0.5, 0, 0.1, 0, "cm")),
        plot.subtitle = element_text(size = 12, color = "gray50")) 


p_t1 <- fu %>% 
  drop_na(psurev) %>%
  ggplot(aes(x = psurev)) +
  geom_histogram(fill = pc[2],
                 binwidth = 10) +
  labs(title = "PSSUQ",
       x = "Total score",
       y = "Count",
       subtitle = "Total score") +
  theming +
  theme(axis.title.y = element_text(angle = 90,
                                    margin = margin(0, 0.4, 0, 0, "cm"))) +
  scale_x_continuous(limits = c(0, 120),
                     breaks = seq(0, 120, 20)) +
  coord_cartesian(ylim = c(0, 250))

p_t2 <- ggplot(fu, aes(x = tamrev)) +
  stat_bin(fill = pc[1],
           binwidth = 10) +
  labs(x = "Total score",
       title = "TAM-FF",
       subtitle = "Total score") +
  theming +
  scale_x_continuous(breaks = seq(0, 120, 20)) +
  coord_cartesian(ylim = c(0, 250))

p_t3 <- ggplot(fu, aes(x = round(useful))) +
  stat_count(fill = pc[3]) +
  theming +
  labs(title = "TAM-FF",
       x = "Mean score",
       subtitle = "Perceived usefulness") +
  coord_cartesian(ylim = c(0, 250)) +
  scale_x_continuous(breaks = 0:7)

p_t4 <- ggplot(fu, aes(x = round(ease))) +
  stat_count(fill = pc[4]) +
  labs(title = "TAM-FF",
       x = "Mean score",
       subtitle = "Perceived ease of use") +
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
       x = "Proportion",
       y = "Count") +
  theme_minimal(base_family = ff) +
  theming +
  theme(axis.title.y = element_text(angle = 90,
                                    margin = margin(0, 0.4, 0, 0, "cm")))

# Make plot for PHQ-8 completions
p_phq8 <- fu %>%
  select(n_phq_f) %>%
  drop_na() %>%
  ggplot(aes(x = n_phq_f)) +
  stat_count(fill = pc[6]) +
  scale_x_discrete(breaks = 0:6,
                   labels = ~ .x) +
  scale_y_continuous(limits = c(0, 300)) +
  labs(title = "Total PHQ-8 completions, next three months",
       x = "Number of PHQ-8 questionnaires returned") +
  theme_minimal(base_family = ff) +
  theming + 
  theme(axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank())

p_dist <- (p_t1 + p_t2 + p_t3 + p_t4 + plot_layout(nrow = 1)) / (p_wt + p_phq8)
ggsave(p_dist, file = here("writing", "figures", "distributions.png"),
       dev = "png",
       width = 10,
       height = 7,
       dpi = 300)

# Counts needed for the text

tabyl(fu$wear_time_l3 == 0) %>% adorn_totals()
summary(fu$wear_time_l3[fu$wear_time_l3 > 0])

tabyl(na.omit(fu$n_phq_f)) %>% adorn_totals()

np <- as.numeric(na.omit(fu$n_phq_f)) - 1
summary(np[np > 0])
as.numeric(np) - 1
summary(as.numeric(fu$n_phq_f[!is.na(fu$n_phq_f) & as.numeric(fu$n_phq_f) > 0]) - 1)

# Supplementary Figure: Describe outcomes by sociodemographics ----------------

# Covariates: age, 
#             male gender
#             number of comorbid conditions
#             years of education
#             currently working
#             recent life events
# Outcomes:   Proportion of time wearing FitBit, next 3 months
#             Total PHQ-8 completions

sup <- sel %>%
  select(subject_id, t,
         age, gender, comorf, edyrs, inwork, n_lte_f,
         n_phq_f, wear_time_l3)  %>%
  group_by(subject_id) %>%
  mutate(across(c(age, gender, comorf, edyrs, n_lte_f, inwork),
                ~ first(na.omit(.x))),
         edyrs_cat = cut(edyrs, 
                         breaks=c(-Inf, 10, 20, Inf), 
                         labels=c("0-10",
                                  "11-20",
                                  "21+")),
         age_cat = cut(age,
                       breaks=c(-Inf, 30, 45, 65, Inf), 
                       labels=c("18-30",
                                "31-45",
                                "46-65",
                                "66+")))

get_summary <- function(var) {
  sup %>%
    group_by({{var}}) %>%
    summarise(n_phq = n_complete(as.numeric(n_phq_f)),
              med_phq = median(as.numeric(n_phq_f), na.rm = TRUE),
              n_wt = n_complete(wear_time_l3),
              med_wt = median(wear_time_l3, na.rm = TRUE))
}


bind_rows(get_summary(age_cat),
          get_summary(gender),
          get_summary(comorf),
          get_summary(edyrs_cat),
          get_summary(inwork)) %>% 
  mutate(cell1 = str_glue("{round(med_phq)} [{n_phq}]"),
         cell2 = str_glue("{sprintf('%.2f', med_wt)} [{n_wt}]")) %>%
  write_csv("~/tidy.csv")

# END.
