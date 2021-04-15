# Preamble ----------------------------------------------------------------
library(tidyverse)
library(rstan)
library(uuid) # for generating a new Universally Unique Identifier
library(parallel)
library(foreach)
library(ggridges)

source("model_scripts/R/utils_new.R")

# use a session ID for similar filenames from same code run
session_ID = substr(uuid::UUIDgenerate(), 1, 8)

## Stan control settings that can be changed for testing
# (for production use 4 chains and at least 1500 iter, 250 warmup)
n_chains = 4
n_iter = 1500
n_warmup = 250
Stan_script = "model_scripts/Stan/work_serosurvey_no_sector_pooling.stan"
random_seed = 1
redo = TRUE

## Stan settings, don't need to change
options(mc.cores = 4)
p_delta = 0.99
n_treedepth = 20

# Define model --------------------------------------------------------------

model = "genre + ageCat"
only_mobilised = FALSE

# WORK study data -------------------------------------------
survey_dat = readRDS(here::here("data", "example_WORK_survey_data_4Stan.rds"))

# Set the reference levels
survey_dat = survey_dat %>%
  mutate(
    genre = fct_relevel(genre, ref = "Female")
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
  select(str_split(model, " \\+ ")[[1]]) %>%
  unique()
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

## model the overall seropositivity with a random effect for company
calc_seropos = run_analysis_stan_re_final(
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
print(calc_seropos$failures)

int_vars_summary = int_probs %>%
  group_by(Sex, ageCat) %>%
  summarise(
    `Rrisk (95% CI)` = ifelse(is.na(mean(rr)), "--",
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
        formatC(3, format = "f")
    )
  )
int_vars_summary %>% print(n = Inf)

# Relative risk for each sector table ----------------------------------------------------------

sec_comp_dat = int_probs %>%
  filter(var_index == 1) %>%
  ungroup() %>%
  group_by(sim) %>%
  mutate(rr_sec = seropos / seropos[sector == 1])

sec_rr_summary <-
  sec_comp_dat %>%
  ungroup() %>%
  group_by(name_sector) %>%
  summarise(
    `Rrisk (95% CI)` = ifelse(is.na(mean(rr_sec)), "--",
      paste0(
        mean(rr_sec, na.rm = T) %>%
          formatC(2, format = "f"),
        " (", quantile(rr_sec, probs = .025, na.rm = T) %>%
          formatC(2, format = "f"), "-",
        quantile(rr_sec, probs = .975, na.rm = T) %>%
          formatC(2, format = "f"), ")"
      )
    ),
    p = ifelse(is.na(mean(rr_sec)), "--",
      min(2 * c(
        mean(rr_sec > 1, na.rm = T),
        mean(rr_sec < 1, na.rm = T)
      )) %>%
        formatC(3, format = "f")
    )
  )
sec_rr_summary %>% print(n = Inf)

cat("\n---- Done ", output_result_list_filename, " at ", format(Sys.time(), "%Y-%m-%d-%H:%M:%S"), "\n\n")

sec_raw = survey_dat %>%
  filter(!is.na(Secteur_group)) %>%
  count(UEP_result, Secteur_group) %>%
  complete(UEP_result, Secteur_group, fill = list(n = 0L)) %>%
  pivot_wider(names_from = c(UEP_result), values_from = n) %>%
  mutate(
    N = `0` + `1`,
    sec_N_pc_of_all = 100 * N / sum(N),
    seropos_N_pc_per_sector = 100 * `1` / N
  ) %>%
  rename("NEGATIF" = "0", "POSITIF" = "1", "name_sector" = "Secteur_group") %>%
  select(name_sector, N, sec_N_pc_of_all, NEGATIF, POSITIF, seropos_N_pc_per_sector)

target_strata = survey_dat %>% count(genre, ageCat, Secteur_group, name = "sec_strata_size")

results = int_probs %>%
  select(-c(var_index, sector, message, rr)) %>%
  inner_join(target_strata, by = c("ageCat", "Sex" = "genre", "name_sector" = "Secteur_group")) %>%
  group_by(sim, name_sector) %>%
  summarise(seroprev = 100 * weighted.mean(seropos, sec_strata_size)) %>%
  ungroup()

sec_seroprev_summary = results %>%
  group_by(name_sector) %>%
  summarise(mean(seroprev), quantile(seroprev, .025), 
            quantile(seroprev, .25), median(seroprev), 
            quantile(seroprev, .75), quantile(seroprev, .975))

sec_raw %>%
  full_join(sec_seroprev_summary) %>%
  full_join(sec_rr_summary) %>%
  write_csv(paste0("results/", session_ID, "-", 
                   format(Sys.time(), "%Y-%m-%d-%H-%M-%S"), 
                   "-SexAgeSectors-seroprev_summary.csv"))
