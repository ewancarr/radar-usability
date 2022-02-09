# title: "Fitting separate models for each outcome"
# author: "Ewan Carr"
# date: "24th November 2021"

renv::activate()
library(here)
library(tidyverse)
library(naniar)
library(broom)
library(gtsummary)
library(dagitty)
library(lme4)
library(broom.mixed)
library(brms)
library(sjPlot)
library(cmdstanr)
library(modelr)
library(tidybayes)
load(here("data", "clean", "merged.Rdata"), verbose = TRUE)
load(here("data", "clean", "selected.Rdata"), verbose = TRUE)
refit <- TRUE
n_iter <- 1e4
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
         pc_ease_r,
         pc_useful_r,
         # Outcomes
         wear_time_l3,
         total_phq8_l3,
         # Covariates
         age,
         comorf,
         edyrs,
         inwork,
         prevwear) %>%
  ungroup() 
nrow(sel)
length(unique(sel$subject_id))

# Standardise continuous variables -------------------------------------------

sel <- sel %>%
  rename(ease = pc_ease_r,
         useful = pc_useful_r) %>%
  mutate(across(c(ids, gad, wsas, psurev, tamrev, 
                  ease, useful, age, edyrs),
                scale,
                .names = "{.col}_1sd"))

# Create numeric ID variable --------------------------------------------------
sel$pid <- as.numeric(as.factor(sel$subject_id))

# Define models ---------------------------------------------------------------

"
We have three exposures (X):

- ids
- gad
- wsas

Four mediators (M):

- psurev
- tamrev
- pc_useful_r
- pc_ease_r

Two outcomes (Y):

- wear_time_l3
- total_phq8_l3

We'll fit a separate model for each 'a', 'b', and 'c' path.
"


###############################################################################
####                                                                      #####
####                                a-paths                               #####
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

d1 <- dagitty('dag {
              "benefit receipt" [pos="0.235,0.358"]
              "in work" [pos="0.343,0.601"]
              "previous wearable use" [pos="0.693,0.267"]
              "switch from iphone" [pos="0.567,0.252"]
              age [pos="0.396,0.140"]
              comorbidities [pos="0.246,0.253"]
              education [pos="0.233,0.460"]
              gender [pos="0.448,0.140"]
              partner [pos="0.329,0.246"]
              usability [outcome,pos="0.576,0.414"]
              wsas [exposure,pos="0.449,0.523"]
              "benefit receipt" -> wsas
              "in work" -> "benefit receipt"
              "in work" -> comorbidities [pos="0.297,0.401"]
              "in work" -> wsas
              "previous wearable use" -> usability
              "switch from iphone" -> usability
              age -> "in work"
              age -> comorbidities
              age -> partner
              age -> usability
              age -> wsas
              comorbidities -> "benefit receipt"
              comorbidities -> wsas
              education -> "benefit receipt"
              education -> "in work"
              education -> comorbidities [pos="0.171,0.374"]
              education -> usability
              education -> wsas
              gender -> wsas
              partner -> usability
              wsas -> usability
              }')

adjustmentSets(d1)
plot(d1)


# Fit linear mixed models -----------------------------------------------------

if (refit) {
  models <- list(
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
  a_path <- map(models, 
                ~ brm(as.formula(.x),
                      data = sel,
                      iter = n_iter,
                      threads = threading(6)))
}


###############################################################################
####                                                                      #####
####                                b-paths                               #####
####                                                                      #####
###############################################################################

"
Wear time
13. wear_time, psurev
14. wear_time, tamrev
15. wear_time, ease
16. wear_time, useful

PHQ-8 completions
17. total_phq8, psurev
18. total_phq8, tamrev
19. total_phq8, use
20. total_phq8, useful
"

d2 <- dagitty('dag {
bb="0,0,1,1"
"benefit receipt" [pos="0.248,0.369"]
"in work" [pos="0.300,0.598"]
"previous wearable use" [pos="0.596,0.333"]
"switch from iphone" [pos="0.530,0.204"]
age [pos="0.396,0.140"]
comorbidities [pos="0.260,0.243"]
education [pos="0.226,0.493"]
gender [pos="0.461,0.162"]
partner [pos="0.304,0.156"]
usability [exposure,pos="0.488,0.338"]
wear_time [outcome,pos="0.617,0.484"]
wsas [pos="0.365,0.495"]
"benefit receipt" -> wsas
"in work" -> "benefit receipt"
"in work" -> comorbidities [pos="0.297,0.401"]
"in work" -> wear_time [pos="0.483,0.573"]
"in work" -> wsas
"previous wearable use" -> usability
"previous wearable use" -> wear_time
"switch from iphone" -> usability
age -> "in work"
age -> comorbidities
age -> partner
age -> usability
age -> wsas
comorbidities -> "benefit receipt"
comorbidities -> wear_time
comorbidities -> wsas
education -> "benefit receipt"
education -> "in work"
education -> comorbidities [pos="0.171,0.374"]
education -> usability
education -> wsas
gender -> wsas
partner -> usability
usability -> wear_time
wsas -> usability
wsas -> wear_time
}')

plot(d2)
adjustmentSets(d2)

# Create ordinal measure of "Number of PHQ-8 completions" ---------------------

# NOTE: we're collapsing the last two categories, so we have 0, 1, ... 7+. This
# is reccomended by Min and Agresti (2005, p.13) "When grouping the count
# values together to form the K categories, our simulation studies suggest that
# when the number of groups is too small, one will lose some efficiency."

sel <- sel %>%
  mutate(total_phq8_l3_of = factor(total_phq8_l3,
                                   ordered = TRUE,
                                   levels = 0:9),
         total_phq8_l3_of7 = factor(if_else(total_phq8_l3 >= 7,
                                            7, total_phq8_l3),
                                    ordered = TRUE,
                                    levels = 0:7))

# Fit model for each combination of mediator and outcome ----------------------

opts <- list(m = c("psurev_1sd",
                   "tamrev_1sd",
                   "ease_1sd",
                   "useful_1sd"),
             y = c("total_phq8_l3_of7",
                   "wear_time_l3")) %>%
  cross()

adjust <- "age_1sd + edyrs_1sd + prevwear + wsas"

b_path <- map(opts, function(i) {
               f <- as.formula(str_glue("{i$y} ~ {i$m} + {adjust} + (1 | pid)"))
               cat(str_glue("\nFitting model for: {i$m} --> {i$y}\n{f}\n\n"))
               if (i$y == "total_phq8_l3_of7") {
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

names(b_path) <- make_labels(opts)

###############################################################################
####                                                                      #####
####                                c-paths                               #####
####                                                                      #####
###############################################################################

"
Wear time
21. ids_1sd --> wear_time_l3
22. gad_1sd --> wear_time_l3
23. wsas_1sd --> wear_time_l3

PHQ-8 completions
24. ids_1sd --> total_phq8_l3_of
25. gad_1sd ---> total_phq8_l3_of
26. wsas_1sd --> total_phq8_l3_of
"

d3 <- dagitty('dag {
bb="0,0,1,1"
"benefit receipt" [pos="0.248,0.369"]
"in work" [adjusted,pos="0.300,0.598"]
"previous wearable use" [adjusted,pos="0.587,0.286"]
"switch from iphone" [pos="0.530,0.204"]
age [adjusted,pos="0.396,0.140"]
comorbidities [adjusted,pos="0.260,0.243"]
education [adjusted,pos="0.228,0.504"]
gender [pos="0.448,0.140"]
partner [pos="0.304,0.156"]
usability [pos="0.488,0.338"]
wear_time [outcome,pos="0.617,0.484"]
wsas [exposure,pos="0.365,0.495"]
"benefit receipt" -> wsas
"in work" -> "benefit receipt"
"in work" -> comorbidities [pos="0.297,0.401"]
"in work" -> wear_time [pos="0.483,0.573"]
"in work" -> wsas
"previous wearable use" -> usability
"previous wearable use" -> wear_time
"switch from iphone" -> usability
age -> "in work"
age -> comorbidities
age -> partner
age -> usability
age -> wsas
comorbidities -> "benefit receipt"
comorbidities -> wear_time
comorbidities -> wsas
education -> "benefit receipt"
education -> "in work"
education -> comorbidities [pos="0.171,0.374"]
education -> usability
education -> wsas
gender -> wsas
partner -> usability
usability -> wear_time
wsas -> usability
wsas -> wear_time
}')

plot(d3)
adjustmentSets(d3)

opts <- cross(list(x = c("ids_1sd", "gad_1sd", "wsas_1sd"),
                   y = c("wear_time_l3", "total_phq8_l3_of7")))
adj <- "age_1sd + comorf + edyrs_1sd + inwork + prevwear"

fit_c <- function(formula, data, fam, iter) {
  fit <- brm(
    formula,
    family = fam,
    data = data,
    iter = iter,
    cores = 24,
    chains = 4,
    backend = "cmdstanr",
    threads = threading(6),
    inits = 0,
    control = list(adapt_delta = 0.99),
    seed = 42
  )
  return(fit)
}

if (refit) {
  c_path <- map(opts, function(i) {
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
                  return(fit_c(formula = f, 
                               data = sel, 
                               fam = fam,
                               iter = n_iter))
  })
}
names(c_path) <- make_labels(opts)

# Save ------------------------------------------------------------------------

save(a_path, b_path, c_path,
     file = here("outputs", "posteriors.Rdata"))

save(sel,
     file = here("data", "clean", "sel.Rdata"))

haven::write_dta(sel, here("data", "clean", "sel.dta"))

# END.
