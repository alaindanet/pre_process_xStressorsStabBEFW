prep_sim_for_analysis <- function (
  file = NULL,
  extinct_bm_threshold = extinction_threshold,
  param_var = var_param_sim_d,
  param_table = sim_param_d2
  ) {
  open_dataset(file, format = "arrow") %>%
  collect() %>%
  # Transform back interaction strength into matrices
  mutate(across(where(is.list), as.list)) %>%
  filter(!is.na(richness)) %>%
  mutate(
    alive_species = map(bm_sp, ~which(.x > extinct_bm_threshold)),
    richness = map_int(alive_species, ~length(.x)),
    bm_sp_alive = map2(bm_sp, alive_species, ~.x[.y])
    ) %>%
  # Get rid of the extinct species to compute per capita int_strength
  mutate(across(c(int_strength, max_int),
      ~map2(.x, alive_species,
        function(y, sp) {
          mat <- matrix(y, nrow = sqrt(length(y)))
          mat[sp, sp, drop = FALSE]
        }
      )
    )
    ) %>%
  mutate(omega = map(omega, ~matrix(.x, nrow = sqrt(length(.x))))) %>%
  # Per capita interaction strength
  mutate(across(c(int_strength, max_int),
      ~map2(.x, bm_sp_alive,
        function(y, sp) {
          y %*% diag(1 / sp, nrow = length(sp))
        }
      )
    )
    ) %>%
  # Put good S for parameter values
  select(-all_of(c(param_var, "S", "ct"))) %>%
  left_join(param_table %>% select(-A)) %>%
  select(sim_id, fw_id, S, ct, h, env_stoch, rho, everything())
} 

make_ts_vector_to_matrices <- function(
  file = NULL,
  ntimestep = NULL 
  ) {
  open_dataset(file, format = "arrow") %>%
    collect() %>%
    # Transform back timeseries strength into matrices
    mutate(across(where(is.list), as.list)) %>%
    filter(!map_lgl(species, is.null)) %>%
    mutate(
      species = map(species,
        ~matrix(.x, nrow = ntimestep)
        ),
      stoch = map(stoch,
        ~matrix(.x, nrow = ntimestep)
      )
    )
}

prep_param_table <- function(file = NULL) {
  open_dataset(file,
      format = "arrow") %>%
    collect() %>%
    rename(env_stoch = sigma) %>%
    mutate(
      ct = map_dbl(A, ~sum(.x) / ((sqrt(length(.x)) - 1)^2))
    )
}

compute_network_metrics <- function(
  sim = NULL
  ) {
  sim %>%
    mutate(
      persistence = richness / S,
      async = 1 / sync,
      stab_pop = 1 / avg_cv_sp,
      ct_alive = map_dbl(max_int,
        ~ifelse(sum(.x) == 0, 0,
          sum(.x > 0) / ((ncol(.x))^2))),
      max_tlvl = map2_dbl(tlvl, alive_species,
        function(lvl, alive) {
          x <- lvl[alive]
          ifelse(length(x) == 0, 0, max(x))
        }),
      max_max_int = map_dbl(max_int,
        ~ifelse(length(.x[.x > 0]) == 0, 0, max(.x[.x > 0]))
        ),
      avg_int_per_cap =  map_dbl(int_strength,
        ~ifelse(length(.x[.x > 0]) == 0, 0, mean(.x[.x > 0]))
        ),
      avg_max_int = map_dbl(max_int,
        ~ifelse(length(.x[.x > 0]) == 0, 0, mean(.x[.x > 0]))
        ),
      sd_max_int = map_dbl(max_int,
        ~ifelse(length(.x[.x > 0]) <= 1, 0, sd(.x[.x > 0]))
        ),
      cv_max_int = sd_max_int / avg_max_int,
      inv_sd_max_int = 1 / sd_max_int,
      gini_max_int = map_dbl(max_int,
        ~ifelse(length(.x[.x > 0]) == 0, 0, gini(.x[.x > 0]))),
      avg_omnivory = map_dbl(omnivory, mean),
      disconnected_prod = map_lgl(max_int,
        ~any(colSums(.x) == 0 & rowSums(.x) == 0))
      ) %>%
  select(!where(is.list))
}

filter_alive_species_ts_matrices <- function(
  ts_mat = NULL,
  prep_sim = NULL
  ) {
  tmp <- ts_mat %>%
    left_join(prep_sim %>%
      select(sim_id, d, alive_species),
    by = "sim_id") %>%
  filter(
    !map_lgl(alive_species, ~length(.x) == 0),
  )
  filter_pb <- tmp %>% 
    filter(
      !map_int(alive_species, ~length(.x)) >
      map_int(stoch, ~ncol(.x))
  )

  pb_sim <- tmp$sim_id[!tmp$sim_id %in% filter_pb$sim_id]
  if (length(pb_sim) > 1) {
    warning("those sim have higher number of alive species than columns in the
      timeseries matrix:", cat(pb_sim, sep = ","))

  }
  filter_pb %>%
    mutate(
      species = map2(species, alive_species,
      function (x, alive) {
        matrix(x, nrow = 500)[, alive, drop = FALSE]
      }
      ),
    stoch = map2(stoch, alive_species,
      function (x, alive) {
        matrix(x, nrow = 500)[, alive, drop = FALSE]
      }
      ),
    d = map2(d, alive_species, ~.x[.y]),
    stoch_d = pmap(list(bm = species, stoch = stoch, d = d),
      function (bm, stoch, d){
        bm %*% diag(d, nrow = length(d)) * exp(stoch)
      }))
  
}

compute_stability_metrics_ts <- function(sim_ts_prep = NULL) {
  sim_ts_prep %>%
    mutate(
      cpe = map_dbl(species, compensatory_effect),
      cpe_env = map_dbl(stoch_d, compensatory_effect),
      cpe_int = cpe / cpe_env,
      stab_com = map_dbl(species, community_stability),
      pop_stab = map_dbl(species, population_stability),
      bm_sp = map(species, ~colMeans(.x)),
      bm_total = map_dbl(bm_sp, ~sum(.x)),
      sd_sp = map(species, ~apply(.x, 2, sd)),
      sum_sd_sp = map_dbl(sd_sp, ~sum(.x)),
      async = map_dbl(species, asynchrony),
      sae = map(species, statistical_averaging_effect),
      sae_total = map_dbl(sae, ~.x["total"]),
      sae_even = map_dbl(sae, ~.x["even"]),
      evenness_sae = map_dbl(sae, ~.x["eveness"])
      ) %>%
  select(-species, -stoch_d, -stoch)

}

get_sim_dataset <- function(
  sim_stab = NULL,
  sim_net = NULL 
  ) {
  left_join(
    sim_net %>%
      select(-async, -stab_com),
    sim_stab,
    by = "sim_id")
} 

filter_sim_fw <- function(x = sim) {
  x %>%
    filter(max_tlvl > 1, !disconnected_prod) %>%
    select(!where(is.list)) %>%
    mutate(resp_div = 1 - rho)
}

get_sem_dataset <- function(
  x = NULL
  ) {
   as.data.frame(x) %>%
     mutate(across(c(stab_com, pop_stab, async, cpe, sae_total, evenness_sae,
           sae_even, cpe_int, cpe_env), log))
}

get_ts_file <- function(dir = NULL, full_names = TRUE) {
  list.files(dir, full.names = full_names) %>% #
    .[str_detect(., "_ts")]
}
get_no_ts_file <- function(dir = NULL, full_names = TRUE) {
  list.files(dir, full.names = full_names) %>%
    .[!str_detect(., "_ts")]
}
