#' @description for given variable and reference
#'
#' @param df serprevalance draws dataframe
#' @param var variable of interest
#' @param ref reference value
#'
#' @return return
computeSeroPrev = function(df, var, ref) {
  # browser()
  colnames(df)[colnames(df) == var] = "val"
  df$var = var
  df %>%
    group_by(sim, val, var) %>%
    summarize(p = weighted.mean(seropos, pop)) %>%
    group_by(sim) %>%
    mutate(coef_val = ifelse(val == ref, NA,
      ifelse(p > p[val == ref], 1, -1)
    ))
}

#' @title Run Stan analysis with random effects
#' @description Runs the Stan seroprevalence model with company random effects
#'
#' @param model_script Stan model script
#' @param dat list with data to run the model
#' @param coef_eqn formula in character format expressing the probability of seropositivity on the logit scale
#' @param pos_control number of predicted positives in the control dataset
#' @param neg_control number of predicted negatives in the control dataset
#' @param control_tp number of true positives in the control dataset
#' @param control_tn number of true negatives in the control dataset
#' @param n_cores number of cores to use for parallel computation
#' @param chains number of chains to use in Stan
#' @param iter number of total iterations per chain in Stan
#' @param warmup number of warmup iterations per chain in Stan
#' @param Stan_control List of control parameters to be passed on to Stan
#' @param seed Random seed
#' @param redo If FALSE redo fit or if TRUE load pre-computed posteriors if available
#' @param session_ID Unique identifier for linking output files
#'
#' @return a list with parameter posteriors and results
run_analysis_stan_re_final = function(model_script,
                                       dat,
                                       coef_eqn,
                                       pos_control,
                                       neg_control,
                                       control_tp,
                                       control_fp,
                                       n_cores = getOption("mc.cores"),
                                       chains,
                                       iter,
                                       warmup,
                                       Stan_control,
                                       seed,
                                       redo = FALSE,
                                       session_ID,
                                       strata,
                                       ...) {
  ## Prepare data for analysis
  # Set analysis data
  dat = as_tibble(dat)

  # Set model matrix
  X = model.matrix(as.formula(paste("~", coef_eqn)), data = dat)

  if (!dir.exists("results")) {
    dir.create("results")
  }

  # name of the stan output file
  stan_sampling_filename = paste0(
    "results/", session_ID, "-stan_fit_",
    basename(model_script), "_", chains,
    "_", iter, "_", warmup, "_",
    format(Sys.time(), "%Y-%m-%d-%H-%M-%S"), ".rds"
  )

  if (!file.exists(stan_sampling_filename) | redo) {
    # Unique company ids from 1 to N_survey
    unique_entreprises = unique(dat$nom_entreprise) %>% sort()
    dat$unique_entreprise_ids = map_dbl(dat$nom_entreprise, ~ which(unique_entreprises == .))
    # Unique sector ids from 1 to N_companies
    sector_comp_dict = dat %>%
      distinct(Secteur_group, nom_entreprise) %>%
      arrange(Secteur_group)
    unique_sectors = unique(sector_comp_dict$Secteur_group)
    # is below really the best way of doing this??
    comp_to_sectors = map_dbl(unique_entreprises, ~ which(unique_sectors == sector_comp_dict$Secteur_group[sector_comp_dict$nom_entreprise == .]))

    # Check if mapping is correct
    tibble(
      nom_entreprise = unique_entreprises,
      Secteur_group = unique_sectors[comp_to_sectors]
    ) %>%
      inner_join(sector_comp_dict,
        by = c("nom_entreprise"),
        suffix = c(".mapped", ".dict")
      ) %>%
      summarise(n_diff = sum(Secteur_group.mapped != Secteur_group.dict)) %>%
      {
        if (.$n_diff != 0) {
          stop("Mapping between companies and sectors failed")
        }
      }

    # Run stan model
    cat(paste0("Starting sampling session ", session_ID, " at ", format(Sys.time(), "%Y-%m-%d-%H:%M:%S"), "\n"))

    stan_posterior = rstan::sampling(stan_model(model_script),
      data = list(
        N_survey = nrow(dat),
        N_companies = length(unique_entreprises),
        N_sectors = length(unique_sectors),
        company = dat$unique_entreprise_ids,
        company_sector_mapping = comp_to_sectors,
        p_vars = ncol(X),
        X = X,
        survey_pos = dat$UEP_result,
        N_pos_control = pos_control,
        control_tp = control_tp,
        N_neg_control = neg_control,
        control_fp = control_fp
      ),
      chains = chains,
      iter = iter,
      warmup = warmup,
      seed = seed,
      control = Stan_control
    )
    cat(paste0("Finished sampling at ", format(Sys.time(), "%Y-%m-%d-%H:%M:%S"), "\nSaving as RDS file ", stan_sampling_filename, "\n"))
    saveRDS(stan_posterior, stan_sampling_filename)
    cat(paste0("Finished saving RDS ", stan_sampling_filename, " at ", format(Sys.time(), "%Y-%m-%d-%H:%M:%S"), "\n"))
  } else {
    cat("Loading pre-computed posteriors from ", stan_sampling_filename, "\n")
    stan_posterior = readRDS(stan_sampling_filename)
  }

  ## extract log-likelihoods
  stan_ll = loo::extract_log_lik(stan_posterior) %>% loo::loo()

  ## extract parameters
  beta = rstan::extract(stan_posterior, pars = "beta")[[1]]
  sigma_sector = rstan::extract(stan_posterior, pars = "sigma_sector")[[1]]
  mu_sector = rstan::extract(stan_posterior, pars = "mu_sector")[[1]]

  sec_pop_mat = strata %>% model.matrix(as.formula(paste("~", coef_eqn)), data = .)

  # name of the integrated seroprev output file
  integrated_probs_filename = paste0(
    "results/",
    session_ID,
    "_integrated-probs_",
    format(Sys.time(), "%Y-%m-%d-%H-%M-%S"),
    ".rds"
  )

  cat(paste0("Starting integrating seroprevalences at ", format(Sys.time(), "%Y-%m-%d-%H:%M:%S"), "\n"))

  if (!file.exists(integrated_probs_filename) | redo) {
    ## setup parallel code
    cl = parallel::makeCluster(n_cores)
    doParallel::registerDoParallel(cl)

    int_probs = foreach(
      i = 1:nrow(sec_pop_mat),
      .combine = rbind,
      .inorder = FALSE,
      .packages = c("tidyverse", "foreach")
    ) %dopar% {
      foreach(
        j = 1:nrow(beta),
        .combine = rbind,
        .inorder = T
      ) %do% {
        foreach(
          k = 1:ncol(sigma_sector),
          .combine = rbind,
          .inorder = T
        ) %do% {
          # Compute probability integrating across company random effects
          prob = integrate(function(x) {
            plogis(qnorm(
              x, beta[j, , drop = F] %*% t(sec_pop_mat[i, , drop = F]) + mu_sector[j, k],
              sigma_sector[j, k]
            ))
          }, 0, 1, rel.tol = 1e-8, stop.on.error = FALSE)

          tibble(
            var_index = i,
            Sex = strata$genre[i],
            ageCat = strata$ageCat[i],
            seropos = prob[[1]],
            sim = j,
            sector = k,
            message = prob[[4]]
          )
        }
      }
    }
    parallel::stopCluster(cl)
    cat(paste0("Finished integrating seroprevalences at ", format(Sys.time(), "%Y-%m-%d-%H:%M:%S"), "\n"))
    saveRDS(int_probs, integrated_probs_filename)
  } else {
    cat("Loading pre-computed seroprevalence estimates from ", integrated_probs_filename, "\n")
    int_probs = readRDS(integrated_probs_filename)
  }

  int_probs = int_probs %>% mutate(name_sector = factor(levels(unique_sectors)[sector],
    levels = levels(unique_sectors)
  ))
  int_probs = int_probs %>%
    group_by(sim, sector) %>%
    mutate(rr = seropos / seropos[var_index == 1]) %>%
    ungroup()

  failures = int_probs[which(int_probs$message != "OK"), ]

  # results list
  res = list(
    beta = beta,
    model_mtx = X,
    mu_sector = mu_sector,
    sigma_sector = sigma_sector,
    sens = extract(stan_posterior, pars = "sens")[[1]],
    spec = extract(stan_posterior, pars = "spec")[[1]],
    N_obs = nrow(dat),
    pos = sum(dat$UEP_result == 1),
    neg = sum(dat$UEP_result == 0),
    stan_ll = stan_ll,
    int_probs = int_probs,
    failures = failures
  )

  return(res)
}

#' @title Run Stan analysis without random effects
#' @description Runs the Stan seroprevalence model
#'
#' @param model_script Stan model script
#' @param dat list with data to run the model
#' @param coef_eqn formula in character format expressing the probability of seropositivity on the logit scale
#' @param pos_control number of predicted positives in the control dataset
#' @param neg_control number of predicted negatives in the control dataset
#' @param control_tp number of true positives in the control dataset
#' @param control_tn number of true negatives in the control dataset
#' @param n_cores number of cores to use for parallel computation
#' @param chains number of chains to use in Stan
#' @param iter number of total iterations per chain in Stan
#' @param warmup number of warmup iterations per chain in Stan
#' @param Stan_control List of control parameters to be passed on to Stan
#' @param seed Random seed
#' @param redo If FALSE redo fit or if TRUE load pre-computed posteriors if available
#' @param session_ID Unique identifier for linking output files
#'
#' @return a list with parameter posteriors and results
run_analysis_stan_uni = function(model_script,
                                  dat,
                                  coef_eqn,
                                  pos_control,
                                  neg_control,
                                  control_tp,
                                  control_fp,
                                  n_cores = getOption("mc.cores"),
                                  chains,
                                  iter,
                                  warmup,
                                  Stan_control,
                                  seed,
                                  redo = FALSE,
                                  session_ID,
                                  strata,
                                  ...) {
  ## Prepare data for analysis
  # Set analysis data
  dat = as_tibble(dat)

  # Set model matrix
  X = model.matrix(as.formula(paste("~", coef_eqn)), data = dat)

  if (!dir.exists("results")) {
    dir.create("results")
  }

  # name of the stan output file
  stan_sampling_filename = paste0(
    "results/", session_ID, "-stan_fit_",
    basename(model_script), "_", chains,
    "_", iter, "_", warmup, "_",
    format(Sys.time(), "%Y-%m-%d-%H-%M-%S"), ".rds"
  )

  if (!file.exists(stan_sampling_filename) | redo) {

    # Run stan model
    cat(paste0("Starting sampling session ", session_ID, " at ", format(Sys.time(), "%Y-%m-%d-%H:%M:%S"), "\n"))

    stan_posterior = rstan::sampling(stan_model(model_script),
      data = list(
        N_survey = nrow(dat),
        p_vars = ncol(X),
        X = X,
        survey_pos = dat$UEP_result,
        N_pos_control = pos_control,
        control_tp = control_tp,
        N_neg_control = neg_control,
        control_fp = control_fp
      ),
      chains = chains,
      iter = iter,
      warmup = warmup,
      seed = seed,
      control = Stan_control
    )
    cat(paste0("Finished sampling at ", format(Sys.time(), "%Y-%m-%d-%H:%M:%S"), "\nSaving as RDS file ", stan_sampling_filename, "\n"))
    saveRDS(stan_posterior, stan_sampling_filename)
    cat(paste0("Finished saving RDS ", stan_sampling_filename, " at ", format(Sys.time(), "%Y-%m-%d-%H:%M:%S"), "\n"))
  } else {
    cat("Loading pre-computed posteriors from ", stan_sampling_filename, "\n")
    stan_posterior = readRDS(stan_sampling_filename)
  }

  ## extract log-likelihoods
  stan_ll = loo::extract_log_lik(stan_posterior) %>% loo::loo()

  ## extract parameters
  beta = rstan::extract(stan_posterior, pars = "beta")[[1]]

  sec_pop_mat = strata %>%
    select(str_split(coef_eqn, " \\+ ")[[1]]) %>%
    model.matrix(as.formula(paste("~", coef_eqn)), data = .)

  # name of the integrated seroprev output file
  integrated_probs_filename = paste0(
    "results/",
    session_ID,
    "_integrated-probs_",
    format(Sys.time(), "%Y-%m-%d-%H-%M-%S"),
    ".rds"
  )

  cat(paste0("Starting integrating seroprevalences at ", format(Sys.time(), "%Y-%m-%d-%H:%M:%S"), "\n"))

  if (!file.exists(integrated_probs_filename) | redo) {
    ## setup parallel code
    cl = parallel::makeCluster(n_cores)
    doParallel::registerDoParallel(cl)

    int_probs = foreach(
      i = 1:nrow(sec_pop_mat),
      .combine = rbind,
      .inorder = FALSE,
      .packages = c("tidyverse", "foreach")
    ) %dopar% {
      foreach(
        j = 1:nrow(beta),
        .combine = rbind,
        .inorder = T
      ) %do% {
        # Compute probability
        prob = plogis(beta[j, , drop = F] %*% t(sec_pop_mat[i, , drop = F]))

        tibble_row(
          var_index = i,
          ageCat = strata$ageCat[i],
          Sex = strata$genre[i],
          variable = strata[i, 3][[1]],
          seropos = as.numeric(prob),
          pop = strata$var_strata_size[i],
          sim = j
        )
      }
    }
    parallel::stopCluster(cl)
    cat(paste0("Finished integrating seroprevalences at ", format(Sys.time(), "%Y-%m-%d-%H:%M:%S"), "\n"))
    saveRDS(int_probs, integrated_probs_filename)
  } else {
    cat("Loading pre-computed seroprevalence estimates from ", integrated_probs_filename, "\n")
    int_probs = readRDS(integrated_probs_filename)
  }

  # results list
  res = list(
    beta = beta,
    model_mtx = X,
    sens = extract(stan_posterior, pars = "sens")[[1]],
    spec = extract(stan_posterior, pars = "spec")[[1]],
    N_obs = nrow(dat),
    pos = sum(dat$UEP_result == 1),
    neg = sum(dat$UEP_result == 0),
    stan_ll = stan_ll,
    int_probs = int_probs
  )

  return(res)
}
