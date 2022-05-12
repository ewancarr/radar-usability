# Title:        Select analytical sample
# Author:       Ewan Carr
# Started:      2022-01-11

library(tidyverse)
library(here)
library(janitor)
library(gtsummary)
library(readxl)
load(here("data", "clean", "merged.Rdata"), verbose = TRUE)

count_n <- function(d) { length(unique(d$subject_id)) }

noisy_drop <- function(d, ...) {
  n_pre <- count_n(d)
  d <- drop_na(d, ...)
  n_post <- count_n(d)
  cat(str_glue("\n{n_pre} --> {n_post} ({n_pre - n_post} dropped)\n"))
  return(d)
}

# Select sample ---------------------------------------------------------------

# Before dropping anything
sel <- merged
cat(count_n(sel))

# [Clinical variables]
sel <- sel %>% noisy_drop(gad, ids, wsas)

# [Usability] At least one measurement of usability and acceptability measures
sel <- sel %>% noisy_drop(psurev, tamrev, useful, ease)

# [Usage] At least one measurement of each usage outcome
sel <- sel %>% noisy_drop(wear_time_l3, n_phq_f)

# Remove missing on covariates
sel <- sel %>% noisy_drop(age, comorf, edyrs, inwork, prevwear, n_lte_f)

count_n(sel)

# As proportion of total sample
round(100 * (count_n(sel) / 623), 1)

# Compare included vs. excluded -----------------------------------------------

included <- filter(merged, subject_id %in% picks)
included$samp <- "included"
excluded <- filter(merged, !(subject_id %in% picks))
excluded$samp <- "excluded"

compare <- bind_rows(included, excluded)

compare %>%
  select(samp,
         ids,
         age,
         male,
         comorf,
         edyrs,
         inwork,
         prevwear) %>%
  summarise(across(everything(), ~ first(na.omit(.x)))) %>%
  ungroup() %>%
  select(-subject_id) %>%
  tbl_summary(by = "samp")

# Export chosen participants
picks <- unique(as.character(sel$subject_id))
sel %>%
  select(subject_id) %>%
  distinct() %>%
  write_csv(here("data", "clean", "selected.csv"))
save(picks, file = here("data", "clean", "selected.Rdata"))

# Save Stata dataset
merged %>% filter(subject_id %in% picks) %>%
  haven::write_dta("~/selected_sample.dta")
