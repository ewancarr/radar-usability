# Title:        Fitting separate models for each outcome
# Author:       Ewan Carr
# Started:      2021-11-24

renv::activate()
library(here)
library(tidyverse)
library(naniar)
library(dagitty)
library(brms)
library(cmdstanr)
library(tidybayes)
load(here("data", "clean", "merged.Rdata"), verbose = TRUE)
load(here("data", "clean", "selected.Rdata"), verbose = TRUE)
load(here("outputs", "dags.Rdata"), verbose = TRUE)
n_iter <- 2e4
set_cmdstan_path("~/.cmdstan/cmdstan-2.28.2")
options(mc.cores = 24,
        brms.backend = "cmdstanr",
        brms.chains = 4)

make_labels <- function(opts) {
  return(map(opts, function(i) {
               paste(str_replace_all(i[[1]], "_", ""),
                     str_replace_all(i[[2]], "_", ""), sep = "_") }))
}

# Select required variables  --------------------------------------------------
sel <- merged %>%
  filter(subject_id %in% picks) %>%
  select(subject_id, t,
         # Exposures
         wsas,
         ids,
         gad,
         # Mediators
         psurev,
         tamrev,
         ease,
         useful,
         # Outcomes
         wear_time_l3,
         n_phq_f,
         # Covariates
         age,
         comorf,
         edyrs,
         inwork,
         prevwear,
         n_lte_f) %>%
  ungroup() 
nrow(sel)
length(unique(sel$subject_id))

# Standardise continuous variables -------------------------------------------

sel <- sel %>%
  mutate(across(c(ids, gad, wsas, psurev, tamrev, 
                  ease, useful, age, edyrs),
                scale,
                .names = "{.col}_1sd"))

# Create numeric ID variable --------------------------------------------------
sel$pid <- as.numeric(as.factor(sel$subject_id))

###############################################################################
####                                                                      #####
####                 H1: Clinical variables --> Usability                 #####
####                                                                      #####
###############################################################################

"
1.   ids_1sd --> psurev
2.   ids_1sd --> tamrev
3.   ids_1sd --> ease
4.   ids_1sd --> useful
5.   gad_1sd --> psurev
6.   gad_1sd --> tamrev
7.   gad_1sd --> ease
8.   gad_1sd --> useful
9.   wsas_1sd --> psurev
10.  wsas_1sd --> tamrev
11.  wsas_1sd --> ease
12.  wsas_1sd --> useful
"

adjustmentSets(d1)

# Fit linear mixed models -----------------------------------------------------

spec_h1 <- list(
  # wsas
  wsas_psu    = "psurev ~ wsas_1sd + age_1sd + edyrs_1sd + (1 | pid)",
  wsas_tam    = "tamrev ~ wsas_1sd + age_1sd + edyrs_1sd + (1 | pid)",
  wsas_useful = "useful ~ wsas_1sd + age_1sd + edyrs_1sd + (1 | pid)",
  wsas_ease   = "ease   ~ wsas_1sd + age_1sd + edyrs_1sd + (1 | pid)",
  # ids
  ids_psu     = "psurev ~ ids_1sd + age_1sd + edyrs_1sd + (1 | pid)",
  ids_tam     = "tamrev ~ ids_1sd + age_1sd + edyrs_1sd + (1 | pid)",
  ids_useful  = "useful ~ ids_1sd + age_1sd + edyrs_1sd + (1 | pid)",
  ids_ease    = "ease   ~ ids_1sd + age_1sd + edyrs_1sd + (1 | pid)",
  # gad
  gad_psu     = "psurev ~ gad_1sd + age_1sd + edyrs_1sd + (1 | pid)",
  gad_tam     = "tamrev ~ gad_1sd + age_1sd + edyrs_1sd + (1 | pid)",
  gad_useful  = "useful ~ gad_1sd + age_1sd + edyrs_1sd + (1 | pid)",
  gad_ease    = "ease   ~ gad_1sd + age_1sd + edyrs_1sd + (1 | pid)"
)


prior_h1 <- list(
  set_prior("normal(-5, 15)", class = "b"),      # X --> psurev [0-112]
  set_prior("normal(-5, 15)", class = "b"),      # X --> tamrev [0-114]
  set_prior("normal(-1, 3)", class = "b"),       # X --> useful [0-7]
  set_prior("normal(-1, 3)", class = "b")        # X --> ease   [0-7]
)

fit_h1 <- map2(spec_h1,
               rep(prior_h1, 3),
                ~ brm(as.formula(.x),
                      data = sel,
                      iter = n_iter,
                      prior = .y,
                      threads = threading(6)))

###############################################################################
####                                                                      #####
####                      H2: Usability --> Usage                         #####
####                                                                      #####
###############################################################################

"
Wear time
13. wear_time, psurev
14. wear_time, tamrev
15. wear_time, ease
16. wear_time, useful

PHQ-8 completions
17. n_phq_of, psurev
18. n_phq_of, tamrev
19. n_phq_of, use
20. n_phq_of, useful
"

adjustmentSets(d2)

# Create ordinal measure of "Number of PHQ-8 completions" ---------------------

# NOTE: we're collapsing the last two categories, so we have 0, 1, ... 7+. This
# is reccomended by Min and Agresti (2005, p.13) "When grouping the count
# values together to form the K categories, our simulation studies suggest that
# when the number of groups is too small, one will lose some efficiency."

sel$n_phq_of = factor(sel$n_phq_f, ordered = TRUE, levels = 0:6)

# Fit model for each combination of mediator and outcome ----------------------

opts <- list(m = c("psurev_1sd",
                   "tamrev_1sd",
                   "ease_1sd",
                   "useful_1sd"),
             y = c("n_phq_of",
                   "wear_time_l3")) %>%
  cross()

adjust <- "age_1sd + edyrs_1sd + prevwear + wsas"

fit_h2 <- map(opts, function(i) {
               f <- as.formula(str_glue("{i$y} ~ {i$m} + {adjust} + (1 | pid)"))
               cat(str_glue("\nFitting model for: {i$m} --> {i$y}\n{f}\n\n"))
               if (i$y == "n_phq_of") {
                 # Fit sequential ordinal model for 'Total PHQ-8 completions'
                 fit <- brm(f,
                            family = sratio("cloglog"),
                            threads = threading(6),
                            data = sel)
               } else {
                 # Fit zero inflated beta for 'FitBit wear time'                     
                 f_zi <- as.formula(str_glue("zi ~ {i$m} + {adjust} + (1 | pid)"))
                 f_phi <- as.formula(str_glue("phi ~ {i$m}"))
                 fit <- brm(formula = bf(f, f_phi, f_zi),
                            family = zero_inflated_beta(),
                            threads = threading(6),
                            data = sel,
                            iter = n_iter,
                            inits = 0,
                            control = list(adapt_delta = 0.99))
               }
               return(fit)
             })

###############################################################################
####                                                                      #####
####                H3: Clinical variables --> Usage                      #####
####                                                                      #####
###############################################################################

"
Wear time
21. ids_1sd --> wear_time_l3
22. gad_1sd --> wear_time_l3
23. wsas_1sd --> wear_time_l3

PHQ-8 completions
24. ids_1sd --> n_phq_of
25. gad_1sd --> n_phq_of
26. wsas_1sd --> n_phq_of
"

spec_h3 <- cross(list(x = c("ids_1sd", "gad_1sd", "wsas_1sd"),
                      y = c("wear_time_l3", "n_phq_of")))

adj <- "age_1sd + comorf + edyrs_1sd + inwork + prevwear + n_lte_f"

fit_single_model <- function(formula, data, fam, iter) {
  fit <- brm(formula,
             family = fam,
             data = data,
             iter = iter,
             cores = 24,
             chains = 4,
             backend = "cmdstanr",
             threads = threading(6),
             inits = 0,
             control = list(adapt_delta = 0.99),
             seed = 42)
  return(fit)
}

fit_h3 <- map(spec_h3, function(i) {
                if (i$y == "wear_time_l3") {
                  fam = zero_inflated_beta()
                  rhs = str_glue("{i$x} + {adj} + (1 | pid)")
                  f <- bf(as.formula(str_glue("wear_time_l3 ~ {rhs}")),
                          as.formula(str_glue("zi ~ {rhs}")),
                          as.formula(str_glue("phi ~ {i$x}")))
                } else {
                  fam = sratio("cloglog")
                  f <- str_glue("{i$y} ~ {i$x} + {adj}")
                }
                return(fit_single_model(formula = f, 
                                        data = sel, 
                                        fam = fam,
                                        iter = n_iter))
})
names(fit_h3) <- make_labels(spec_h3)

# Save ------------------------------------------------------------------------

# Posteriors
save(fit_h1, file = here("outputs", "fit_h1.Rdata"))
save(fit_h2, file = here("outputs", "fit_h2.Rdata"))
save(fit_h3, file = here("outputs", "fit_h3.Rdata"))

# Data
save(sel, file = here("data", "clean", "sel.Rdata"))
haven::write_dta(sel, here("data", "clean", "sel.dta"))

# END.
