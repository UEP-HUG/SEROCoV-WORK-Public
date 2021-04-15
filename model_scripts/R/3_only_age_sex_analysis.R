# Preamble ----------------------------------------------------------------
library(tidyverse)
library(rstan)
library(uuid) # for generating a new Universally Unique Identifier
library(parallel)
library(foreach)
library(ggridges)
library(rlang)

source("model_scripts/R/utils_new.R")

# use a session ID for similar filenames from same code run
session_ID = substr(uuid::UUIDgenerate(), 1, 8)

## Stan control settings that can be changed for testing
# (for production use 4 chains and at least 1500 iter, 250 warmup)
n_chains = 4
n_iter = 1500
n_warmup = 250
Stan_script = "model_scripts/Stan/work_serosurvey_univariate.stan"
random_seed = 1
redo = TRUE

## Stan settings, don't need to change
options(mc.cores = 4)
p_delta = 0.99
n_treedepth = 20

# Define model --------------------------------------------------------------
model = "genre + ageCat"

only_mobilised = TRUE

cat("The model is ****", model, "****\n")

# WORK study data -------------------------------------------
survey_dat = readRDS(here::here("data", "example_WORK_survey_data_4Stan.rds"))

# Set the reference levels
sex_ref = "Female"
age_ref = "[18,35)"

survey_dat = survey_dat %>%
  mutate(
    genre = fct_relevel(genre, ref = sex_ref),
    ageCat = fct_relevel(ageCat, ref = age_ref)
  )

survey_dat = survey_dat %>%
  mutate(
    UEP_result = fct_recode(
      UEP_result,
      "0" = "NEGATIF", "1" = "POSITIF"
    ) %>%
      as.character() %>%
      as.integer()
  ) %>%
  drop_na(str_split(model, " \\+ ")[[1]])

if (only_mobilised) {
  survey_dat = survey_dat %>% filter(is_totally_confined2 == "FALSE")
  session_ID = paste0(session_ID, "_onlyMobilised")
}

strata = survey_dat %>%
  count(ageCat, genre, name = "var_strata_size") %>% # NB Also could do across(all_of(univariate)), or even !!as.name(univariate), instead but not sure they're much clearer??
  arrange(genre, ageCat)

# Control data ------------------------------------------------------------

## bring in lab validation data (from Meyer et al, 2020 doi: 10.1016/j.cmi.2020.06.024
## "Validation of a commercially available SARS-CoV-2 serological immunoassay")

## number of positive controls from validation data
pos_control = 181 # from Table 1 Samples "COVID-19 Patient" All
## number of negative controls
neg_control = 326 # from Table 1 Samples "Negative Controls"

## number of true positives for cases
control_tp = 154 # Table 1 ELISA IgG Positive for "COVID-19 Patient" All

## number of false positives for controls (1-specificity)
control_fp = 4 # Table 1 ELISA IgG Positive for "Negative Controls"

# Run models --------------------------------------------------------------

output_result_list_filename = paste0(
  "results/",
  session_ID,
  "_results-list_",
  format(Sys.time(), "%Y-%m-%d-%H-%M-%S"),
  ".rds"
)

## model the overall seropositivity
## NB Make sure you realise the line in this function after the probability
## calculation variable = strata[i,3][[1]] should work but is meaningless!!
calc_seropos = run_analysis_stan_uni(
  model_script = Stan_script,
  dat = survey_dat,
  coef_eqn = model,
  pos_control = pos_control,
  neg_control = neg_control,
  control_tp = control_tp,
  control_fp = control_fp,
  n_cores = getOption("mc.cores"),
  chains = n_chains,
  iter = n_iter,
  warmup = n_warmup,
  Stan_control = list(
    adapt_delta = p_delta,
    max_treedepth = n_treedepth
  ),
  seed = random_seed,
  redo = redo,
  session_ID = session_ID,
  strata
)
cat(paste0("Saving RDS ", output_result_list_filename, " at ", format(Sys.time(), "%Y-%m-%d-%H:%M:%S"), "\n"))
saveRDS(calc_seropos, output_result_list_filename)

# Relative risk for each variable table ----------------------------------------------------------

int_probs = calc_seropos$int_probs

## overall estimate
overall_probs = int_probs %>%
  mutate(
    var = "Overall",
    val = ""
  ) %>%
  group_by(sim, var, val) %>%
  summarize(p = weighted.mean(seropos, pop)) %>%
  ungroup()

## find age specific probabilities in order to make relative risks
age_probs = int_probs %>%
  filter(Sex == sex_ref) %>%
  mutate(var = "Age") %>%
  rename(val = ageCat) %>%
  group_by(sim, var, val) %>%
  summarize(p = weighted.mean(seropos, pop)) %>%
  ungroup()

# sex-specific probabilities
sex_probs = int_probs %>%
  filter(ageCat == age_ref) %>%
  mutate(var = "Sex") %>%
  rename(val = Sex) %>%
  group_by(sim, var, val) %>%
  summarize(p = weighted.mean(seropos, pop)) %>%
  ungroup()

subset_est = bind_rows(
  overall_probs,
  sex_probs,
  age_probs
)

## Compute age estimates
age_seroprev = computeSeroPrev(int_probs, "ageCat", age_ref) %>%
  bind_rows(
    subset_est %>%
      filter(var == "Overall") %>%
      mutate(var = "ageCat", val = "all") %>%
      mutate(coef_val = as.numeric(NA))
  )

age_counts = survey_dat %>%
  rbind(survey_dat %>% mutate(ageCat = "all")) %>%
  group_by(ageCat) %>%
  summarize(
    n = n(),
    pos = sum(UEP_result == 1),
    neg = n() - pos
  ) %>%
  mutate(ageCat = as.character(ageCat))

age_res = left_join(age_seroprev, age_counts, by = c("val" = "ageCat")) %>%
  mutate(var = "Age")

## Compute sex estimates
sex_seroprev = computeSeroPrev(int_probs, "Sex", sex_ref)

sex_counts = survey_dat %>%
  mutate(val = genre) %>%
  group_by(val) %>%
  summarize(
    n = n(),
    pos = sum(UEP_result == 1),
    neg = n() - pos
  )

sex_res = left_join(sex_seroprev, sex_counts)

# Combine results for seropositive estimate
seropos = bind_rows(age_res, sex_res) %>%
  group_by(var, val, n, pos, neg) %>%
  summarize(
    `Seroprevalence (95% CI)` = paste0(
      mean(100 * p) %>%
        formatC(1, format = "f"), " (",
      quantile(100 * p, probs = .025) %>%
        formatC(1, format = "f"), "-",
      quantile(100 * p, probs = .975) %>%
        formatC(1, format = "f"), ")"
    ),
    p = ifelse(is.na(mean(coef_val)), "--",
      min(2 * c(mean(coef_val > 0), mean(coef_val < 0))) %>%
        formatC(4, format = "f")
    )
  ) %>%
  ungroup() %>%
  mutate(
    pos = paste0(pos, " (", formatC(100 * pos / n, 2, format = "f"), "%)"),
    neg = paste0(neg, " (", formatC(100 * neg / n, 2, format = "f"), "%)")
  ) %>%
  rename(
    Category = val, Obs = n, `Test positive` = pos, `Test negative` = neg
  ) %>%
  mutate(Category = factor(Category, levels = c(levels(survey_dat$ageCat), "Male", "Female", "all"))) %>%
  arrange(Category)

# Compute relative risk
rrisk = subset_est %>%
  filter(var == "Age") %>%
  group_by(sim) %>%
  mutate(rr = ifelse(val == age_ref, NA, p / p[val == age_ref])) %>%
  ungroup() %>%
  left_join(
    survey_dat %>%
      group_by(ageCat) %>%
      summarize(
        n = n(),
        pos = sum(UEP_result == 1),
        neg = n() - pos
      ) %>% mutate(ageCat = as.character(ageCat)),
    by = c("val" = "ageCat")
  ) %>%
  bind_rows(
    subset_est %>%
      filter(var == "Sex") %>%
      group_by(sim) %>%
      mutate(rr = ifelse(val == sex_ref, NA, p / p[val == sex_ref])) %>%
      ungroup() %>%
      left_join(survey_dat %>%
        mutate(val = genre) %>%
        group_by(val) %>%
        summarize(
          n = n(),
          pos = sum(UEP_result == 1),
          neg = n() - pos
        ))
  ) %>%
  group_by(var, val, n, pos, neg) %>%
  summarise(
    `Relative risk (95% CI)` = ifelse(is.na(mean(rr)), "--",
      paste0(
        mean(rr, na.rm = T) %>%
          formatC(2, format = "f"),
        " (", quantile(rr, probs = .025, na.rm = T) %>%
          formatC(2, format = "f"), "-",
        quantile(rr, probs = .975, na.rm = T) %>%
          formatC(2, format = "f"), ")"
      )
    ),
    p = ifelse(is.na(mean(rr)), "--",
      min(2 * c(
        mean(rr > 1, na.rm = T),
        mean(rr < 1, na.rm = T)
      )) %>%
        formatC(4, format = "f")
    )
  ) %>%
  ungroup() %>%
  mutate(
    pos = paste0(pos, " (", formatC(100 * pos / n, 2, format = "f"), "%)"),
    neg = paste0(neg, " (", formatC(100 * neg / n, 2, format = "f"), "%)")
  ) %>%
  rename(
    `Test positive` = pos, `Test negative` = neg,
    Obs = n, Category = val
  ) %>%
  select(-var) %>%
  mutate(Category = factor(Category, levels = c(levels(survey_dat$ageCat), "Male", "Female"))) %>%
  arrange(Category)

# Write table
final_table = seropos %>%
  select(-p) %>%
  left_join(rrisk %>%
    select(Category, `Relative risk (95% CI)`, p)) %>%
  arrange(Category) %>%
  replace_na(list(`Relative risk (95% CI)` = "--", p = "--"))
print(final_table)

res_table_filename = paste0(
  "results/",
  session_ID,
  "_results_table_",
  str_replace_all(model, " \\+ ", "-"),
  format(Sys.time(), "_%Y-%m-%d-%H-%M-%S"),
  ".csv"
)
final_table %>% write_csv(res_table_filename)
