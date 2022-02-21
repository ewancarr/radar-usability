# Title:        Select analytical sample
# Author:       Ewan Carr
# Started:      2022-01-11

library(tidyverse)
library(here)
library(janitor)
library(gtsummary)
library(readxl)
load(here("data", "clean", "merged.Rdata"), verbose = TRUE)

# Select sample ---------------------------------------------------------------

# Before dropping anything
sel <- merged
length(unique(sel$subject_id))

# [Usage] At least one measurement of each outcome
sel <- sel %>% drop_na(wear_time_l3, total_phq8_l3)
length(unique(sel$subject_id))

# [Usability] At least one measurement of each mediator
sel <- merged %>% drop_na(psurev, tamrev, useful, ease)
length(unique(sel$subject_id))

# [Clinical variables]
sel <- sel %>% drop_na(gad, ids, wsas)
length(unique(sel$subject_id))

# Remove missing on covariates
sel <- sel %>% drop_na(age, comorf, edyrs, inwork, prevwear, n_lte_f)
length(unique(sel$subject_id))

# Compare included vs. excluded -----------------------------------------------

compare <- merged %>%
  filter(!(pid %in% sel$pid)) %>%
  mutate(samp = "excluded") %>%
  bind_rows(mutate(sel, samp = "analytical"))

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
