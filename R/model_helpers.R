get_sem_model_list <- function (
  sim_fw_sem = sim_fw_sem
  ) {
  list(
    ct_alive = lm(ct_alive ~ ct + S, sim_fw_sem),
    richness = lm(richness ~ ct + S, sim_fw_sem),
    #max_tlvl = lm(max_tlvl ~ ct + S + Z, sim_fw_sem),
    w_avg_tlvl = lm(w_avg_tlvl ~ ct + S + Z, sim_fw_sem),
    avg_omnivory = lm(avg_omnivory ~ ct + S, sim_fw_sem),
    avg_int_strength = lm(avg_int_strength ~ ct + S + h + Z, sim_fw_sem),
    pop_stab = lm(pop_stab ~ richness + avg_int_strength + w_avg_tlvl + ct_alive + env_stoch + h + Z, sim_fw_sem),
    evenness_sae = lm(evenness_sae ~ richness + avg_int_strength + w_avg_tlvl + ct_alive + h + Z, sim_fw_sem),
    sae_even = lm(sae_even ~ richness, sim_fw_sem),
    cpe_int = lm(cpe_int ~ richness + avg_int_strength + w_avg_tlvl + ct_alive + h + Z, sim_fw_sem),
    cpe_env = lm(cpe_env ~ richness + resp_div, sim_fw_sem),
    async = lm(async ~ sae_even + evenness_sae + cpe_env + cpe_int, sim_fw_sem),
    stab_com = lm(stab_com ~ pop_stab + async, sim_fw_sem)
  )
}

get_sem_simplified_model_list <- function (
  sim_fw_sem = NULL 
  ) {
  list(
    ct_alive = lm(ct_alive ~ ct + S, sim_fw_sem),
    richness = lm(richness ~ ct + S, sim_fw_sem),
    #max_tlvl = lm(max_tlvl ~ ct + S + Z, sim_fw_sem),
    w_avg_tlvl = lm(w_avg_tlvl ~ ct + S + Z, sim_fw_sem),
    avg_omnivory = lm(avg_omnivory ~ ct + S + Z, sim_fw_sem),
    avg_int_strength = lm(avg_int_strength ~ ct + S + Z, sim_fw_sem),
    pop_stab = lm(pop_stab ~ richness + avg_int_strength + w_avg_tlvl + ct_alive + env_stoch + h, sim_fw_sem),
    evenness_sae = lm(evenness_sae ~ richness + avg_int_strength + w_avg_tlvl + ct_alive + h, sim_fw_sem),
    #sae_even = lm(sae_even ~ richness, sim_fw_sem),
    sae_total = lm(sae_total ~ richness + evenness_sae, sim_fw_sem),
    cpe = lm(cpe ~ richness + avg_int_strength + w_avg_tlvl + ct_alive + h + resp_div, sim_fw_sem),
    async = lm(async ~ sae_total + cpe, sim_fw_sem),
    stab_com = lm(stab_com ~ pop_stab + async, sim_fw_sem)
  )
}

get_sem <- function(model_list = NULL) {
    ti <- as.psem(model_list)
    update(ti, avg_omnivory %~~% w_avg_tlvl)
}

get_sem_coeff <- function(sem = NULL) {
  ti <- summary(sem)$coefficients[, c(1, 2, 8, 9)]
  colnames(ti)[4] <- "Signif"
  ti

}

get_sem_tot_effect <- function(
  sem = NULL, 
  sem_data = NULL, 
  bootstrap_nb = 100,
  nb_cores = future::availableCores(),
  ci_type = "perc"
  ) {
  semEff(sem,
    ci.type = ci_type,
    R = bootstrap_nb,
    seed = 13,
    parallel = "snow",
    data = sem_data,
    ncpus = nb_cores 
    )
}
