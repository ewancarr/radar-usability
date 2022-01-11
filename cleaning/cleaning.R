# Title:        Data cleaning for RADAR-MDD analysis
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

# Load survey data
dat <- read_dta(here("data", "survey", "extended_data_2021_09_30.dta")) %>%
    rename(event = redcap_event_name,
           fut = Followuptime,
           psu = PSSUQ_TOTAL,
           eth = ETHCAT2,
           tam = TAM_TOTAL) %>%
    clean_names() %>%
    mutate(
        event = if_else(event == "enrolment_arm_1",
            "0",
            str_replace(
                event,
                "_month_assessmen[t]*_arm_1",
                ""
            )
        ),
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
        debt = case_when(
            utilities_10 == 1 ~ FALSE,
            utilities_10 == 0 ~ TRUE,
            TRUE ~ NA
        ),
        prevwear = as.numeric(fitness_tracker_used),
        bame = if_else(is.na(eth), NA, as.numeric(eth) %in% 2:5)
    ) %>%
    drop_na(pid, t, fut)

# Get extra information in larger dataset -------------------------------------

extra <- read_dta(here("data", "survey",
                       "SwitchandPhysicalHealthData.dta")) %>%
    clean_names() %>%
    select(subject_id,
        starts_with("mh_longst_illness_"),
        smartphone_type,
        employ = csri_7,
        enrol_date = enrolmentdate
    ) %>%
    # Remove 'Depression' comorbidity
    select(-mh_longst_illness_type_4a) %>%
    mutate(across(starts_with("mh_longst_illness_type_"), as.numeric),
        comor = rowSums(across(starts_with("mh_longst_illness_type_"))),
        comorf = factor(if_else(comor > 2, 2, comor),
            levels = 0:2,
            labels = c("0", "1", "2+")
        ),
        switchedphone = as.numeric(smartphone_type),
        enrol_date = ymd(enrol_date),
        inwork = case_when(employ %in% c(0, 2, 5) ~ TRUE,
                           employ %in% c(1, 3, 4, 6, 7, 8, 
                                         9, 10, 11, 12) ~ FALSE,
                           TRUE ~ NA)
    )

# # Select analytical sample ----------------------------------------------------

# # ==> Must have at least one measure of PSU and TAM.

# incl <- dat %>%
#     select(pid, psu, tam) %>%
#     drop_na() %>%
#     count(pid) %>%
#     pluck("pid")

# Select baseline variables ---------------------------------------------------

baseline <- dat %>%
    group_by(pid) %>%
    select(
        subject_id, pid, age, male, bame, anyben, smoker, partner, phq,
        edyrs, site, debt, prevwear
    ) %>%
    summarise(across(everything(), ~ first(na.omit(.x))))

# Merge with extra variables (smartphone, physical comorbidities) -------------

baseline <- baseline %>%
    left_join(extra, by = "subject_id")

# Select longitudinal variables -----------------------------------------------

longitudinal <- dat %>%
    select(pid, t, survey_time, psu, tam, wsas, gad,
           starts_with("tam_usefulness"),
           starts_with("tam_ease_of_use"),
           future_1 = tam_usage_1,
           future_2 = tam_usage_2,
           future_3 = tam_usage_3,
           future_4 = tam_usage_4
           ) %>%
    drop_na(pid, t)

# Merge baseline and longitudinal measures ------------------------------------

clean <- longitudinal %>%
    left_join(baseline, by = "pid")

# Recode some variables -------------------------------------------------------

clean <- clean %>%
    # Reverse direction of PSU and TAM
    mutate(
        psurev = max(psu) - psu,
        tamrev = max(tam) - tam,
        tf = as.factor(t),
        across(starts_with("future_"), as.numeric)
    ) %>%
    rowwise() %>%
    mutate(# Future use
           fut = sum(c_across(c(future_1, future_2, future_3, future_4))),
           # Perceived usefulness (mean score, reversed)
           pc_useful = mean(c_across(starts_with("tam_usefulness"))),
           pc_useful_r = 7 - pc_useful,
           # Perceived ease of use (mean score, reversed)
           pc_ease = mean(c_across(starts_with("tam_ease_of_use"))),
           pc_ease_r = 7 - pc_ease) %>%
    ungroup()

# Select required variables ---------------------------------------------------

clean <- clean %>%
    select(
        subject_id, pid, t, survey_time, psu, tam, wsas, gad, fut, age, male,
        smoker, partner, phq, edyrs, site, debt, prevwear, comor, comorf,
        enrol_date, smartphone_type, tf, tamrev, psurev, employ, inwork,
        pc_useful_r, pc_ease_r
    ) %>%
    arrange(pid, t)

# Derive dates ----------------------------------------------------------------

survey <- clean %>%
    arrange(subject_id, t) %>%
    complete(subject_id, t) %>%
    group_by(subject_id) %>%
    mutate(enrol_date = first(na.omit(enrol_date))) %>%
    mutate(fu_date = enrol_date %m+% months(t),
           since = interval(enrol_date, fu_date) / days(1),
           t_d29 = round(floor(since / 29)))

# Manual fix: replace 25 with 24.
survey$t_d29[survey$t_d29 == 25] <- 24

# Save ------------------------------------------------------------------------

save(survey, file = here("data", "clean", "survey.Rdata"))
write_dta(survey, path = here("data", "clean", "survey.dta"))

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
    filter(subject_id %in% unique(clean$subject_id))

# Derive timepoints -----------------------------------------------------------

# Get enrolment date from survey data
app <- select(clean, subject_id, enrol_date) %>%
    distinct() %>%
    right_join(app)

# Calculate 't_d29' (i.e. Dan's time variable) by dividing number of days
# between enrolment and now by 29.

app <- app %>%
    mutate(since = interval(enrol_date, e_ts) / days(1),
        t_d29 = round(floor(since / 29)))

app %>% filter(subject_id == '007751c5-d7ad-4bec-a58f-abf32500e2ae') %>%
    select(t_d29) %>% count(t_d29)

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
  select(subject_id, matches("_[0-9]+$")) %>%
  gather(k, value, -subject_id) %>%
  mutate(t_d29 = parse_number(str_extract(k, "[0-9]+$")),
         measure = str_match(k, "^(.*)_[0-9]+$")[, 2]) %>%
  select(-k) %>%
  spread(measure, value) %>%
  arrange(subject_id, t_d29)

###############################################################################
####                                                                      #####
####                          Merge all datasets                          #####
####                                                                      #####
###############################################################################

merge_keys <- c("subject_id", "t_d29")

merged <- survey %>%
    full_join(from_dan, by = merge_keys) %>%
    full_join(app_data, by = merge_keys)

###############################################################################
####                                                                      #####
####                 Derive measures of app usage over NEXT 3 months      #####
####                                                                      #####
###############################################################################

cum_lead <- function(vec) {
    lead <- dplyr::lead
    return(lead(vec, 1) + lead(vec, 2) + lead(vec, 3))
}
 
leads <- merged %>%
    arrange(subject_id, t_d29) %>%
    group_by(subject_id) %>%
    mutate(across(c(wear_time, e_totdurmin, n_access, total_phq8,
                    total_esm, total_esm28q),
                  cum_lead, .names = "{.col}_l3")) %>%
    filter(t_d29 %in% c(3, 12)) %>%
    drop_na(psu)

###############################################################################
####                                                                      #####
####                                 Save                                 #####
####                                                                      #####
###############################################################################

save(merged, leads, file = here("data", "clean", "clean.Rdata"))
