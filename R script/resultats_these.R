# ------------------------------------------- #
# Rscript for the simulation of the scenarios #
# ------------------------------------------- #


library(tidyverse)
library(multibrasBOP2) # github.com/GuillaumeMulier/multibrasBOP2

source("R script/metadata.R")


# I/ Getting the different thresholds ----

# TOP and BOP for efficacy alone
load("data/cutoff_planifie.Rdata")
seuil_effseule <- seuil_effseule[[1]]

# TOP and BOP for efficacy and toxicity
load("data/scenar20240522.Rdata")
# 2 analyses
seuil_toxclassique <- tableau_scenars %>% 
  filter(scenar == "H0", nom_corr == "R<sub>pos,1</sub>") %>% 
  select(C__classique, gamma_classique, alpha_calc_classique, puissance_calc_classique) %>% 
  as.numeric() %>% 
  'names<-'(c("C_", "gamma", "alpha_calc", "puissance_calc"))
seuil_toxrapprochee <- tableau_scenars %>% 
  filter(scenar == "H0", nom_corr == "R<sub>pos,1</sub>") %>% 
  select(C__toxrapp, gamma_toxrapp, alpha_calc_toxrapp, puissance_calc_toxrapp) %>% 
  as.numeric() %>% 
  'names<-'(c("C_", "gamma", "alpha_calc", "puissance_calc"))

# Simon design r1 and r are computed inside the function
# You can see in the script thresholds.R how it is done
# r1 = 5 and r = 16


# II/ Simulate the scenarios for every correlation ----

# Create the dataframe that will have all the scenarios and store the results
tableau_results_Rvar <- tableau_scenars[, -c(2:9)]
tableau_results_Rvar[, c("nb_pts_topiva", "temps_etude_topiva", "accept_ttt_topiva", "arret_prema_topiva", "arret_fut_topiva", "arret_tox_topiva", "nb_eff_topiva", "nb_tox_topiva",
                         "nb_pts_bopiva", "temps_etude_bopiva", "accept_ttt_bopiva", "arret_prema_bopiva", "arret_fut_bopiva", "arret_tox_bopiva", "nb_eff_bopiva", "nb_tox_bopiva",
                         "nb_pts_simon", "temps_etude_simon", "accept_ttt_simon", "arret_prema_simon", "arret_fut_simon", "arret_tox_simon", "nb_eff_simon", "nb_tox_simon",
                         "moy_pts_classique", "moy_duree_classique", "rejet_h0_classique", "arret_prema_classique", "arret_fut_classique", "arret_tox_classique", "nb_eff_classique", "nb_tox_classique",
                         "moy_pts_toxrapp", "moy_duree_toxrapp", "rejet_h0_toxrapp", "arret_prema_toxrapp", "arret_fut_toxrapp", "arret_tox_toxrapp", "nb_eff_toxrapp", "nb_tox_toxrapp",
                         "moy_pts_toxprior", "moy_duree_toxprior", "rejet_h0_toxprior", "arret_prema_toxprior", "arret_fut_toxprior", "arret_tox_toxprior", "nb_eff_toxprior", "nb_tox_toxprior")] <- NA_real_
tableau_results_Rvar[, c("analyse_med_topiva", "analyse_med_bopiva", "analyse_med_simon")] <- NA_character_

# Simulate the trials and store the results
for (i in seq_len(nrow(tableau_results_Rvar))) {
  cat(paste0(tableau_results_Rvar$scenar[i], "/", tableau_results_Rvar$nom_corr[i], " : ", i, " / ", nrow(tableau_results_Rvar), "\n"))
  p_reel <- tableau_results_Rvar$probas[i] %>% pluck(1) %>% unname()
  # The trial and the cohort of MET+ patients  
  essais <- gen_patients_multinomTOP(n_sim = 10000,
                                     ana_inter = c(30, 51),
                                     ana_inter_tox = c(rep(5, 4), rep(10, 5), 11),
                                     interpatient = 6, max_teff = 180, max_ttox = 42, rand_ratio = 1,
                                     multinom_ttt = list(p_reel),
                                     seed = 1993)
  cohorte_2 <- lapply(
    X = essais,
    FUN = function(x) {
      res <- data.frame(temps_recrutement = sort(runif(9, 0, max(x$temps_recrutement))))
      res$tox <- rbinom(9, 1, p_reel[1] + p_reel[3])
      res$temps_tox <- gen_temps_patient(n = 9, tmax = 42)
      res$temps_tox <- ifelse(res$tox == 1, res$temps_recrutement + res$temps_tox, 1e+12)
      res$temps_obstox <- res$temps_recrutement + 42
      return(res)
    }
  )
  
  # TOP + PPtox
  tableau_results_Rvar[i, c("nb_pts_topiva", "temps_etude_topiva", 
                            "arret_fut_topiva", "arret_tox_topiva", "nb_eff_topiva", 
                            "nb_tox_topiva", "accept_ttt_topiva", "arret_prema_topiva", "analyse_med_topiva")] <- 
    simu_topiva(
      ana_inter_eff = c(30, 51),
      ana_inter_tox = c(rep(5, 4), rep(10, 6), 1),
      liste_essais = essais,
      liste_cohortes = cohorte_2,
      phieff = .15, phitox = .25,
      prior_eff = .15, prior_tox = c(1, 1), critere_tox = .95,
      C_ = seuil_effseule[["C_"]], gamm = seuil_effseule[["gamma"]],
      interpatient = 6, max_teff = 180, max_ttox = 42
    ) %>% 
    `colnames<-`(paste0(colnames(.), "_topiva")) %>% 
    tibble_row()
  
  # BOP + PPtox
  tableau_results_Rvar[i, c("nb_pts_bopiva", "temps_etude_bopiva", 
                            "arret_fut_bopiva", "arret_tox_bopiva", "nb_eff_bopiva", 
                            "nb_tox_bopiva", "accept_ttt_bopiva", "arret_prema_bopiva", "analyse_med_bopiva")] <- 
    simu_bopiva(
      ana_inter_eff = c(30, 51),
      ana_inter_tox = c(rep(5, 4), rep(10, 6), 1),
      liste_essais = essais,
      liste_cohortes = cohorte_2,
      phitox = .25,
      prior_tox = c(1, 1), critere_tox = .95,
      C_ = seuil_effseule[["C_"]], gamm = seuil_effseule[["gamma"]],
      p_n = tableau_results_Rvar$probas[tableau_results_Rvar$scenar == "H0" & tableau_results_Rvar$nom_corr == tableau_results_Rvar$nom_corr[i]] %>% pluck(1) %>% unname()
    ) %>% 
    `colnames<-`(paste0(colnames(.), "_bopiva")) %>% 
    tibble_row()
  
  # Simon + PPtox
  tableau_results_Rvar[i, c("nb_pts_simon", "temps_etude_simon", 
                            "arret_fut_simon", "arret_tox_simon", "nb_eff_simon", 
                            "nb_tox_simon", "accept_ttt_simon", "arret_prema_simon", "analyse_med_simon")] <- 
    simu_simon(
      ana_inter_eff = c(30, 51),
      ana_inter_tox = c(rep(5, 4), rep(10, 6), 1),
      liste_essais = essais,
      liste_cohortes = cohorte_2,
      pu = .15, pa = .3, alpha_simon = .1, 
      phitox = .25, critere_tox = .95, prior_tox = c(1, 1)
    ) %>% 
    `colnames<-`(paste0(colnames(.), "_simon")) %>% 
    tibble_row()
  
  # TOP eff/tox with 2 analyses
  tableau_results_Rvar[i, c("moy_pts_classique", "moy_duree_classique", "rejet_h0_classique", "arret_prema_classique",
                            "arret_fut_classique", "arret_tox_classique", "nb_eff_classique", 
                            "nb_tox_classique")] <- 
    simu_essais(alpha = .05,
                ana_inter = c(30, 51),
                ana_inter_tox = c(30, 51),
                p_n = tableau_results_Rvar$probas[tableau_results_Rvar$scenar == "H0" & tableau_results_Rvar$nom_corr == tableau_results_Rvar$nom_corr[i]] %>% pluck(1) %>% unname(),
                p_a = tableau_results_Rvar$probas[tableau_results_Rvar$scenar == "H1" & tableau_results_Rvar$nom_corr == tableau_results_Rvar$nom_corr[i]] %>% pluck(1) %>% unname(),
                p_reel = p_reel,
                liste_tableaux = essais,
                cut_seq = seuil_toxclassique[["C_"]],
                power_seq = seuil_toxclassique[["gamma"]],
                interpatient = 6, max_teff = 180, max_ttox = 42,
                stock_matrice = FALSE) %>%
    `names<-`(paste0(names(.), "_classique")) %>% 
    as_tibble_row()
  
  # TOP eff/tox with close monitoring of toxicity
  tableau_results_Rvar[i, c("moy_pts_toxrapp", "moy_duree_toxrapp", "rejet_h0_toxrapp", "arret_prema_toxrapp",
                            "arret_fut_toxrapp", "arret_tox_toxrapp", "nb_eff_toxrapp", 
                            "nb_tox_toxrapp")] <- 
    simu_essais(alpha = .05,
                ana_inter = c(30, 51),
                ana_inter_tox = c(rep(5, 4), rep(10, 5), 11),
                p_n = tableau_results_Rvar$probas[tableau_results_Rvar$scenar == "H0" & tableau_results_Rvar$nom_corr == tableau_results_Rvar$nom_corr[i]] %>% pluck(1) %>% unname(),
                p_a = tableau_results_Rvar$probas[tableau_results_Rvar$scenar == "H1" & tableau_results_Rvar$nom_corr == tableau_results_Rvar$nom_corr[i]] %>% pluck(1) %>% unname(),
                p_reel = p_reel,
                liste_tableaux = essais,
                cut_seq = seuil_toxrapprochee[["C_"]],
                power_seq = seuil_toxrapprochee[["gamma"]],
                interpatient = 6, max_teff = 180, max_ttox = 42,
                stock_matrice = FALSE) %>% 
    `names<-`(paste0(names(.), "_toxrapp")) %>% 
    as_tibble_row()
  
  # TOP eff/tox with update of prior for toxicity
  tableau_results_Rvar[i, c("moy_pts_toxprior", "moy_duree_toxprior", "rejet_h0_toxprior", "arret_prema_toxprior",
                            "arret_fut_toxprior", "arret_tox_toxprior", "nb_eff_toxprior", 
                            "nb_tox_toxprior")] <- 
    simu_essais_prior(alpha = .05,
                      ana_inter = c(30, 51),
                      ana_inter_tox = c(rep(5, 4), rep(10, 5), 11),
                      ana_par = 9,
                      p_n = tableau_results_Rvar$probas[tableau_results_Rvar$scenar == "H0" & tableau_results_Rvar$nom_corr == tableau_results_Rvar$nom_corr[i]] %>% pluck(1) %>% unname(),
                      p_a = tableau_results_Rvar$probas[tableau_results_Rvar$scenar == "H1" & tableau_results_Rvar$nom_corr == tableau_results_Rvar$nom_corr[i]] %>% pluck(1) %>% unname(),
                      p_reel = p_reel,
                      liste_tableaux = essais,
                      liste_tableaux_par = cohorte_2,
                      cut_seq = seuil_toxrapprochee[["C_"]],
                      power_seq = seuil_toxrapprochee[["gamma"]],
                      interpatient = 6, max_teff = 180, max_ttox = 42,
                      stock_matrice = FALSE) %>% 
    `names<-`(paste0(names(.), "_toxprior")) %>% 
    as_tibble_row()
}

save(tableau_results_Rvar, file = paste0("data/resultats_tox30_rposvar_", format(Sys.Date(), "%Y%m%d"), ".Rdata"))
