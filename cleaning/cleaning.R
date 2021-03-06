# Title:        Data cleaning for RADAR-MDD usability analysis
# Author:       Ewan Carr
# Started:      2021-08-03

renv::load()
library(tidyverse)
library(fs)
library(here)
library(haven)
library(janitor)
library(lubridate)
library(naniar)
library(readxl)

parse_event <- function(event) {
  if_else(event == "enrolment_arm_1",
          "0",
          str_replace(event, "_month_assessmen[t]*_arm_1", ""))
}

# Load survey data ------------------------------------------------------------

ed <- read_dta(here("data",
                    "survey",
                    "extended_data_2021_09_30.dta")) %>%
    rename(fut = Followuptime,
           psu = PSSUQ_TOTAL,
           eth = ETHCAT2,
           tam = TAM_TOTAL) %>%
    clean_names() %>%
    mutate(event = parse_event(redcap_event_name),
           t = as.numeric(event),
           pid = as.numeric(as.factor(subject_id)),
           male = as.numeric(gender) == 0,
           partner = as.numeric(marital_status) == 1,
           edyrs = as.numeric(education_years),
           smoker = as.numeric(mh_smoking) == 2,
           anyben = as.numeric(benefits_type_10 == 0),
           wsas = as.numeric(wsas_total),
           gad = as.numeric(gad7_total),
           phq = as.numeric(total_phq8),
           ids = as.numeric(ids_total),
           site = as.numeric(recruitmentsite),
           survey_time = ymd_hms(na_if(pssuq_timestamp, "[not completed]")),
           debt = case_when(utilities_10 == 1 ~ FALSE,
                            utilities_10 == 0 ~ TRUE,
                            TRUE ~ NA),
           prevwear = as.numeric(fitness_tracker_used),
           bame = if_else(is.na(eth), NA, as.numeric(eth) %in% 2:5)) %>%
    drop_na(pid, t, fut) %>%
    select(-enrolment_date,
           -redcap_event_name,
           -ends_with("_timestamp"),
           -event,
           -survey_time)

# Get extra information in larger dataset -------------------------------------

td <- read_dta(here("data",
                    "survey",
                    "totaldataset.dta")) %>%
  clean_names() %>%
  select(subject_id,
         redcap_event_name,
         smartphone_type,
         enrolment_date,
         employ = csri_7,
         ids_date = id_sdate,
         lte_sr_timestamp, ids_sr_timestamp,
         cidi_sf_timestamp, gad7_timestamp,
         wsas_timestamp, bipq_timestamp,
         csri_timestamp, pssuq_timestamp,
         starts_with("mh_longst_illness_"),
         starts_with("lte"),
         # Remove 'Depression' comorbidity
         -mh_longst_illness_type_4a) %>%
  mutate(event = parse_event(redcap_event_name),
         t = as.numeric(event))

# Get date of each follow-up assessment ---------------------------------------

td <- td %>%
  mutate(# Get first non-missing value from various timestamps
         across(c(lte_sr_timestamp, ids_sr_timestamp,
                  cidi_sf_timestamp, gad7_timestamp,
                  wsas_timestamp, bipq_timestamp,
                  csri_timestamp, pssuq_timestamp),
                ~ case_when(.x == "" ~ NA_character_,
                            .x == "[not completed]" ~ NA_character_,
                            nchar(.x) < 15 ~ paste0(.x, " 12:00:00"),
                            TRUE ~ .x)),
         survey_timestamp = coalesce(lte_sr_timestamp, ids_sr_timestamp,
                                     cidi_sf_timestamp, gad7_timestamp,
                                     wsas_timestamp, bipq_timestamp,
                                     csri_timestamp, pssuq_timestamp),
         survey_timestamp = ymd_hms(survey_timestamp),
         # Set survey date as IDS date or timestamp
         survey_date = coalesce(ids_date, as_date(survey_timestamp)))

# Manually replace "survey_time" for two participants with data from Faith ----

td[
    td$subject_id == "35287589-fcfd-423f-8a08-500dbcec2922" &
    td$t == 12, 
  ]$survey_date <- ymd("2020-11-05")

td[
    td$subject_id == "deeeadb9-ead2-42c3-8efe-afebf1f08f13" &
    td$t == 3, 
  ]$survey_date <- ymd("2019-04-18")

# Derive comorbid conditions, other variables ---------------------------------

td <- td %>%
  mutate(across(starts_with("mh_longst_illness_type_"), as.numeric),
         comor = rowSums(across(starts_with("mh_longst_illness_type_"))),
         comorf = factor(if_else(comor > 2, 2, comor),
                         levels = 0:2,
                         labels = c("0", "1", "2+")),
         # Other variables
         switchedphone = as.numeric(smartphone_type),
         enrol_date = ymd(enrolment_date),
         inwork = case_when(employ %in% c(0, 2, 5) ~ TRUE,
                            employ %in% c(1, 3, 4, 6, 7, 8, 
                                          9, 10, 11, 12) ~ FALSE,
                            TRUE ~ NA))

# Lifetime events -------------------------------------------------------------

recode_lte <- function(dat) {
  group_by(dat, subject_id, t) %>%
  summarise(n_lte = sum(in_the_last),
            n_lte_f = factor(case_when(n_lte %in% 0:2 ~ n_lte,
                                       n_lte > 2 ~ as.integer(3),
                                       TRUE ~ NA_integer_),
                             levels = 0:3,
                             labels = c("0", "1", "2", "3+")))
}

# At 3 and 12 months: use "ltefu"
lte_fu <- td %>%
  select(subject_id, t, ltefu_1:ltefu_12) %>%
  pivot_longer(ltefu_1:ltefu_12,
               names_prefix = "ltfu_") %>%
  filter(t %in% c(3, 12),
         name != "ltefu_13") %>%
  mutate(in_the_last = value == 1,
         item = parse_number(name)) %>%
  select(subject_id, t, item, in_the_last) %>%
  recode_lte()

# At 0 months/enrolment: use lte == "in the last year"
lte_0m <- td %>%
  filter(t == 0) %>%
  select(subject_id, t, lte_1:lte_12a) %>%
  pivot_longer(lte_1:lte_12a,
               names_prefix = "lte_") %>%
  mutate(in_the_last = value == 1,
         item = parse_number(name)) %>%
  select(subject_id, t, item, in_the_last) %>%
  recode_lte()

lte <- bind_rows(lte_0m, lte_fu)
td <- full_join(td, lte, by = c("subject_id", "t"))

# Merge both survey datasets --------------------------------------------------

table(unique(ed$subject_id) %in% unique(td$subject_id)) 
survey <- full_join(ed, td, by = c("subject_id", "t"))

# Check: how many people are missing survey_time?
survey %>%
  filter(t %in% c(3, 12),
         is.na(survey_date)) %>%
  pluck("subject_id")

# Copy enrolment values of covariates into subsequent rows --------------------

survey <- survey %>%
    group_by(subject_id) %>%
    mutate(across(c(enrol_date, pid, age, male, bame,
                  anyben, smoker, partner, phq,
                  comorf, switchedphone, inwork,
                  edyrs, site, debt, prevwear),
                  ~ first(na.omit(.x))))

# Recode longitudinal measures ------------------------------------------------

survey <- survey %>%
  rename(future_1 = tam_usage_1,
         future_2 = tam_usage_2,
         future_3 = tam_usage_3,
         future_4 = tam_usage_4) %>%
    mutate(psurev = (19 * 7) - psu,             # 19 questions, each scored 0-7
           tamrev = (16 * 7) - tam,             # 16 questions, each scored 0-7
           tf = as.factor(t),
           across(starts_with("future_"), as.numeric)) %>%
    rowwise() %>%
    mutate(# Future use
           fut = sum(c_across(c(future_1, future_2, future_3, future_4))),
           # Perceived usefulness (mean score, reversed)
           pc_useful = mean(c_across(starts_with("tam_usefulness"))),
           useful = 7 - pc_useful,
           # Perceived ease of use (mean score, reversed)
           pc_ease = mean(c_across(starts_with("tam_ease_of_use"))),
           ease = 7 - pc_ease) %>%
    ungroup()

# Derive dates ----------------------------------------------------------------

survey <- survey %>%
    arrange(subject_id, t) %>%
    complete(subject_id, t) %>%
    group_by(subject_id) %>%
    mutate(enrol_date = first(na.omit(enrol_date))) %>%
    mutate(fu_date = enrol_date %m+% months(t),
           since = interval(enrol_date, fu_date) / days(1),
           t_d29 = round(floor(since / 29)))

# Manual fix: replace 25 with 24.
survey$t_d29[survey$t_d29 == 25] <- 24

###############################################################################
####                                                                      #####
####                Calculate change in clinical variables                #####
####                                                                      #####
###############################################################################

survey <- survey %>%
  complete(subject_id, t)  %>%
  arrange(subject_id, t) %>%
  group_by(subject_id) %>%
  mutate(across(c(ids, wsas, gad),
                ~ .x - lag(.x),
                .names = "{.col}_lag"))

###############################################################################
####                                                                      #####
####                        Prepare mobile app data                       #####
####                                                                      #####
###############################################################################

app <- dir_ls(here("data", "firebase", "2021-09-16"), glob = "*.csv") %>%
    map_dfr(read_csv,
        col_types = "ccccfccccc",
        skip = 1,
        col_names = c("index_time",
                      "index",
                      "event_date",
                      "event_timestamp",
                      "event_name",
                      "user_pseudo_id",
                      "subjectId",
                      "engagement_time_msec",
                      "message_type",
                      "projectId")) %>%
    clean_names() %>%
    mutate(e_ts = as_datetime(as.numeric(event_timestamp) / 1e6),
           e_durmsec = as.numeric(engagement_time_msec),
           e_month = month(e_ts),
           e_yr = year(e_ts)) %>%
    filter(subject_id %in% unique(survey$subject_id))

# Derive timepoints -----------------------------------------------------------

# Get enrolment date from survey data
app <- select(survey, subject_id, enrol_date) %>%
    distinct() %>%
    right_join(app)

# Calculate 't_d29' (i.e. Dan's time variable) by dividing number of days
# between enrolment and now by 29.
app <- app %>%
    mutate(since = interval(enrol_date, e_ts) / days(1),
        t_d29 = round(floor(since / 29)))

# Total duration using the app ------------------------------------------------

duration <- app %>%
    select(subject_id, t_d29, e_ts, e_durmsec, e_month, e_yr) %>%
    group_by(subject_id, t_d29) %>%
    summarise(e_totdurms = sum(e_durmsec, na.rm = TRUE)) %>%
    mutate(e_totdurmin = e_totdurms / 60000)

# Number of times accessing the app -------------------------------------------

access <- app %>%
    filter(event_name == "user_engagement",
           engagement_time_msec > 5000) %>%
    group_by(subject_id, t_d29) %>%
    count(name = "n_access")

app_data <- duration %>%
  full_join(access) %>%
  mutate(n_access = replace_na(n_access, 0))

###############################################################################
####                                                                      #####
####                         Prepare data from Dan                        #####
####                                                                      #####
###############################################################################

from_dan <- read_csv(here("data", "from_dan", "data_requested.csv")) %>%
    select(-`...1`) %>%
    clean_names()

# Check: what is the first timepoint for each person?
from_dan <- from_dan %>%
  select(subject_id, days_in_study, matches("_[0-9]+$")) %>%
  gather(k, value, -subject_id) %>%
  mutate(t_d29 = parse_number(str_extract(k, "[0-9]+$")),
         measure = str_match(k, "^(.*)_[0-9]+$")[, 2]) %>%
  select(-k) %>%
  spread(measure, value) %>%
  arrange(subject_id, t_d29)

###############################################################################
####                                                                      #####
####                    Prepare RADAR PHQ-8 data                          #####
####                                                                      #####
###############################################################################

lookup <- survey %>%
  select(subject_id, t, survey_date) %>%
  filter(t %in% c(3, 12),
         !is.na(survey_date)) %>%
  mutate(window_start = survey_date,
         window_end = survey_date + weeks(12))

radphq <- read_csv(here("data", "radar", "phq8_data.csv")) %>%
  clean_names()

phq8_timings <- radphq %>%
  select(subject_id = user_id, 
         time_completed_utc) %>%
  filter(subject_id %in% lookup$subject_id) %>%
  mutate(time_completed_utc = as_datetime(time_completed_utc)) %>%
  select(subject_id, phq_time = time_completed_utc)

phq8 <- lookup %>%
  group_by(subject_id, t) %>%
  group_split() %>%
  map_dfr(function(i) {
    x <- phq8_timings[(phq8_timings$subject_id == i$subject_id &
                       phq8_timings$phq_time > i$window_start &
                       phq8_timings$phq_time < i$window_end), ]
    x$wk <- week(x$phq_time)
    # Use slice_sample to select a single PHQ-8 response per week
    x <- group_by(x, wk) %>% slice_sample(n = 1)
    i$n_phq <- nrow(x)
    return(i) }) %>%
  select(subject_id, t, n_phq)

phq8$n_phq_f <- factor(if_else(phq8$n_phq > 6, as.integer(6), phq8$n_phq),
                       levels = 0:6,
                       labels = paste0(0:6))

###############################################################################
####                                                                      #####
####                 Derive measures of app usage over NEXT 3 months      #####
####                                                                      #####
###############################################################################

cum_lead <- function(vec) {
    lead <- dplyr::lead
    return(lead(vec, 1) + lead(vec, 2) + lead(vec, 3))
}
 
next_3m <- from_dan %>%
    arrange(subject_id, t_d29) %>%
    group_by(subject_id) %>%
    mutate(across(c(wear_time, total_phq8),
                  cum_lead, .names = "{.col}_l3")) %>%
    select(subject_id, t_d29, ends_with("_l3"))

###############################################################################
####                                                                      #####
####                          Merge all datasets                          #####
####                                                                      #####
###############################################################################

# NOTE: not merging "app_data", since we're not using these measures in paper.
merged <- survey %>%
  filter(t %in% c(0, 3, 12)) %>%
  left_join(next_3m, by = c("subject_id", "t_d29")) %>%
  left_join(phq8, by = c("subject_id", "t"))

# Rescale 'wear time' from 0-3 to 0-1
merged$wear_time_l3 <- merged$wear_time_l3 / 3

###############################################################################
####                                                                      #####
####                     Remove withdrawn participants                    #####
####                                                                      #####
###############################################################################

withdrawals <- read_xlsx(here("data", "survey",
                              "Withdrawal information.xlsx")) %>%
  clean_names() %>%
  mutate(last_app = case_when(apps_disconnected == "Between 3-6 months" ~ 3,
                              apps_disconnected == "Between 12-15 months" ~ 12))

merged <- merged %>%
  left_join(withdrawals, by = c("subject_id" = "participant_id"))  %>%
  mutate(excluded = t >= last_app) 

# Check/inspect excluded observations
merged %>%
  group_by(subject_id) %>%
  filter(any(!is.na(excluded))) %>% 
  select(subject_id, t, psu, last_app, excluded)

# Remove them
merged <- filter(merged, !excluded | is.na(excluded))

###############################################################################
####                                                                      #####
####                                 Save                                 #####
####                                                                      #####
###############################################################################

save(merged, file = here("data", "clean", "merged.Rdata"))

# END.
