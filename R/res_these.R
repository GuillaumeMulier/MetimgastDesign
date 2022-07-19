# ------------------------------------------------------------- #
# Script des résultats pour la thèse sur l'essai METIMGAST      #
# Créé le 28/09/2021, modifié le 18/10/2021                     #
# ------------------------------------------------------------- #


# Packages et helpers ----

# devtools::install_github("GuillaumeMulier/multibrasBOP2", auth_token = "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")
# La clef sera à actualiser tous les 30 jours sur github

library(tidyverse)
library(ggtext)
library(multibrasBOP2)
library(patchwork)

theme_set(theme_light() +
            theme(strip.background = element_rect(fill = NA),
                  strip.text = element_textbox(
                    size = 12, 
                    color = "white", fill = "#7888C0", box.color = "#000066",
                    halign = 0.5, linetype = 1, r = unit(3, "pt"), width = unit(0.75, "npc"),
                    padding = margin(2, 0, 1, 0), margin = margin(3, 3, 3, 3))))


# J'ai repris les 30% de toxicité comme valeur seuil +++

# I/ Récupération des différents seuils ----

# Pour le TOP uniquement pour l'efficacité
load("data/cutoff_planifie.Rdata")
seuil_effseule <- seuil[[1]]
rm(seuil)

# Pour le TOP pour efficacité et toxicité
# Récupération dans le même temps des différents scénarios
# On ne se place que pour une corrélation à minima légèrement positive
load("data/resultats_tox30_20210907.Rdata")
seuil_toxclassique <- tableau_scenars %>% 
  filter(scenar == "H0", nom_corr == "Rpos1") %>% 
  select(C__classique, gamma_classique, alpha_calc_classique, puissance_calc_classique) %>% 
  as.numeric() %>% 
  'names<-'(c("C_", "gamma", "alpha_calc", "puissance_calc"))
seuil_toxrapprochee <- tableau_scenars %>% 
  filter(scenar == "H0", nom_corr == "Rpos1") %>% 
  select(C__toxrapp, gamma_toxrapp, alpha_calc_toxrapp, puissance_calc_toxrapp) %>% 
  as.numeric() %>% 
  'names<-'(c("C_", "gamma", "alpha_calc", "puissance_calc"))
liste_scenars <- tableau_scenars %>% 
  filter(nom_corr == "Rpos1") %>% 
  select(scenar, probas, EffTox, EffNotox, NoeffTox, NoeffNotox, R)
rm(tableau_scenars)


# II/ Simulations des 10 scénarios ----

source("R/realisations_essais_metimgast.R", encoding = "utf-8")

tableau_results <- liste_scenars
tableau_results[, c("analyse_med_topiva", "nb_pts_topiva", "temps_etude_topiva", 
                    "accept_ttt_topiva", "arret_prema_topiva", "arret_fut_topiva", "arret_tox_topiva", "nb_eff_topiva", 
                    "nb_tox_topiva")] <- NA_real_
tableau_results[, c("analyse_med_bopiva", "nb_pts_bopiva", "temps_etude_bopiva", 
                    "accept_ttt_bopiva", "arret_prema_bopiva", "arret_fut_bopiva", "arret_tox_bopiva", "nb_eff_bopiva", 
                    "nb_tox_bopiva")] <- NA_real_
tableau_results[, c("moy_pts_classique", "moy_duree_classique", "rejet_h0_classique", "arret_prema_classique",
                    "arret_fut_classique", "arret_tox_classique", "nb_eff_classique", 
                    "nb_tox_classique")] <- NA_real_
tableau_results[, c("moy_pts_toxrapp", "moy_duree_toxrapp", "rejet_h0_toxrapp", "arret_prema_toxrapp",
                    "arret_fut_toxrapp", "arret_tox_toxrapp", "nb_eff_toxrapp", 
                    "nb_tox_toxrapp")] <- NA_real_
tableau_results[, c("moy_pts_toxprior", "moy_duree_toxprior", "rejet_h0_toxprior", "arret_prema_toxprior",
                    "arret_fut_toxprior", "arret_tox_toxprior", "nb_eff_toxprior", 
                    "nb_tox_toxprior")] <- NA_real_

for (i in seq_len(nrow(liste_scenars))) {
    
    cat(paste0(liste_scenars$scenar[i], " : ", i, " / ", nrow(liste_scenars), "\n"))
    
    p_reel <- liste_scenars$probas[i] %>% pluck(1) %>% unname()
    
    essais <- gen_patients_multinomTOP(n_sim = 10000,
                                       ana_inter = c(30, 51),
                                       ana_inter_tox = c(rep(5, 4), rep(10, 5), 11),
                                       interpatient = 6, max_teff = 180, max_ttox = 42, rand_ratio = 1,
                                       multinom_ttt = list(p_reel),
                                       seed = 1993)
    # Generation de la 2ème cohorte qui est cMET négative
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
    
    # TOP avec approche Ivanova au temps du TOP
    resultats_topiva <- map2_dfr(.x = essais,
                                 .y = cohorte_2,
                                 .f = ~ realisation_essai_topiva(ana_inter_eff = c(30, 51),
                                                                 ana_inter_tox = c(rep(5, 4), rep(10, 6), 1),
                                                                 tableau = .x,
                                                                 cohorte = .y,
                                                                 phitox = c(.15, .25), prior_eff = .15, prior_tox = c(1, 1), 
                                                                 C_ = seuil_effseule[["C_"]], gamm = seuil_effseule[["gamma"]], critere_tox = .95,
                                                                 interpatient = 6, max_teff = 180, max_ttox = 42),
                                 .id = "essai")
    tableau_results[i, c("nb_pts_topiva", "temps_etude_topiva", 
                      "arret_fut_topiva", "arret_tox_topiva", "nb_eff_topiva", 
                      "nb_tox_topiva", "accept_ttt_topiva", "arret_prema_topiva", "analyse_med_topiva")] <- 
      resultats_topiva %>% 
      summarise(across(c(nb_pts, temps_etude, arret_fut:nb_tox), ~ mean(.x)),
                accept_ttt = mean(decision == "Acceptation traitement"),
                arret_prema = mean(str_detect(decision, "prématuré")),
                analyse_med_topiva = paste0("med = ", median(analyse), " / ", 
                                            paste(names(table(.$analyse)), "=", round(100 * prop.table(table(.$analyse)), 1), "%", collapse = " - "))) %>% 
      `colnames<-`(paste0(colnames(.), "_topiva")) %>% 
      tibble_row()
    
    # BOP avec approche Ivanova donc observation complète à chaque analyse intermédiaire
    resultats_bopiva <- map2_dfr(.x = essais,
                                 .y = cohorte_2,
                                 .f = ~ realisation_essai_bopiva(ana_inter_eff = c(30, 51),
                                                                 ana_inter_tox = c(rep(5, 4), rep(10, 6), 1),
                                                                 tableau = .x,
                                                                 cohorte = .y,
                                                                 phitox = c(.15, .25), prior_eff = .15, prior_tox = c(1, 1), 
                                                                 p_n = liste_scenars$probas[liste_scenars$scenar == "H0"] %>% pluck(1) %>% unname(), 
                                                                 C_ = seuil_effseule[["C_"]], gamm = seuil_effseule[["gamma"]], critere_tox = .95,
                                                                 interpatient = 6, max_teff = 180, max_ttox = 42),
                                 .id = "essai")
    tableau_results[i, c("nb_pts_bopiva", "temps_etude_bopiva", 
                         "arret_fut_bopiva", "arret_tox_bopiva", "nb_eff_bopiva", 
                         "nb_tox_bopiva", "accept_ttt_bopiva", "arret_prema_bopiva", "analyse_med_bopiva")] <- 
      resultats_bopiva %>% 
      summarise(across(c(nb_pts, temps_etude, arret_fut:nb_tox), ~ mean(.x)),
                accept_ttt = mean(decision == "Acceptation traitement"),
                arret_prema = mean(str_detect(decision, "prématuré")),
                analyse_med_bopiva = paste0("med = ", median(analyse), " / ", 
                                            paste(names(table(.$analyse)), "=", round(100 * prop.table(table(.$analyse)), 1), "%", collapse = " - "))) %>% 
      `colnames<-`(paste0(colnames(.), "_bopiva")) %>% 
      tibble_row()
    
    # TOP classique avec 2 analyses intermédiaires
    resultats_topclas <- simu_essais(alpha = .05,
                                     ana_inter = c(30, 51),
                                     ana_inter_tox = c(30, 51),
                                     p_n = liste_scenars$probas[liste_scenars$scenar == "H0"] %>% pluck(1) %>% unname(),
                                     p_a = liste_scenars$probas[liste_scenars$scenar == "H1"] %>% pluck(1) %>% unname(),
                                     p_reel = p_reel,
                                     liste_tableaux = essais,
                                     cut_seq = seuil_toxclassique[["C_"]],
                                     power_seq = seuil_toxclassique[["gamma"]],
                                     interpatient = 6, max_teff = 180, max_ttox = 42,
                                     stock_matrice = FALSE)
    tableau_results[i, c("moy_pts_classique", "moy_duree_classique", "rejet_h0_classique", "arret_prema_classique",
                         "arret_fut_classique", "arret_tox_classique", "nb_eff_classique", 
                         "nb_tox_classique")] <- 
      resultats_topclas %>%
      `names<-`(paste0(names(resultats_topclas), "_classique")) %>% 
      as_tibble_row()
    
    # TOP avec monitoring rapproché de la toxicité
    resultats_toprapp <- simu_essais(alpha = .05,
                                     ana_inter = c(30, 51),
                                     ana_inter_tox = c(rep(5, 4), rep(10, 5), 11),
                                     p_n = liste_scenars$probas[liste_scenars$scenar == "H0"] %>% pluck(1) %>% unname(),
                                     p_a = liste_scenars$probas[liste_scenars$scenar == "H1"] %>% pluck(1) %>% unname(),
                                     p_reel = p_reel,
                                     liste_tableaux = essais,
                                     cut_seq = seuil_toxrapprochee[["C_"]],
                                     power_seq = seuil_toxrapprochee[["gamma"]],
                                     interpatient = 6, max_teff = 180, max_ttox = 42,
                                     stock_matrice = FALSE)
    tableau_results[i, c("moy_pts_toxrapp", "moy_duree_toxrapp", "rejet_h0_toxrapp", "arret_prema_toxrapp",
                         "arret_fut_toxrapp", "arret_tox_toxrapp", "nb_eff_toxrapp", 
                         "nb_tox_toxrapp")] <- 
      resultats_toprapp %>% 
      `names<-`(paste0(names(resultats_toprapp), "_toxrapp")) %>% 
      as_tibble_row()
    
    # TOP en incorporant les données de la cohorte 2 dans la toxicité
    resultats_toxprior <- simu_essais_prior(alpha = .05,
                                            ana_inter = c(30, 51),
                                            ana_inter_tox = c(rep(5, 4), rep(10, 5), 11),
                                            ana_par = 9,
                                            p_n = liste_scenars$probas[liste_scenars$scenar == "H0"] %>% pluck(1) %>% unname(),
                                            p_a = liste_scenars$probas[liste_scenars$scenar == "H1"] %>% pluck(1) %>% unname(),
                                            p_reel = p_reel,
                                            liste_tableaux = essais,
                                            liste_tableaux_par = cohorte_2,
                                            cut_seq = seuil_toxrapprochee[["C_"]],
                                            power_seq = seuil_toxrapprochee[["gamma"]],
                                            interpatient = 6, max_teff = 180, max_ttox = 42,
                                            stock_matrice = FALSE)
    
    tableau_results[i, c("moy_pts_toxprior", "moy_duree_toxprior", "rejet_h0_toxprior", "arret_prema_toxprior",
                         "arret_fut_toxprior", "arret_tox_toxprior", "nb_eff_toxprior", 
                         "nb_tox_toxprior")] <- 
      resultats_toxprior %>% 
      `names<-`(paste0(names(resultats_toxprior), "_toxprior")) %>% 
      as_tibble_row()

}

save(tableau_results, file = paste0("data/resultats_tox30_rpos_", format(Sys.Date(), "%Y%m%d"), ".Rdata"))
rm(tableau_results, resultats_toxprior, resultats_toprapp, resultats_topclas, resultats_bopiva, resultats_topiva, cohorte_2, essais, p_reel)


# III/ Reprise des corrélations avec corrélation constant pour le TOP ----

## Estimation des risques alpha et puissances pour chaque schéma dans chaque situation

source("R/realisations_essais_metimgast.R", encoding = "utf-8")

generer_corr_multi <- function(vec_inter,
                               eff,
                               tox) {
  
  data.frame(Pet = vec_inter) %>% 
    mutate(probas = map(Pet, ~ c(Pefftox = .x, Peffnotox = eff - .x, Pnoefftox = tox - .x, Pnoeffnotox = 1 - eff - tox + .x)),
           correlation = map_dbl(probas, ~ (.x[1] - eff * tox) / sqrt((eff - eff * eff) * (tox - tox * tox))))
  
}

# Calculer tous les risques pour toutes les hypothèses possibles
tab_probas <- bind_rows(
  generer_corr_multi(seq(0, .15, .015), .15, .3) %>% mutate(hypothese = "H0"),
  generer_corr_multi(seq(0, .2, .02), .3, .2) %>% mutate(hypothese = "H1")
) 
tab_probas[, c("rejet_topiva", "rejet_bopiva", "rejet_classique", "rejet_toxrapp", "rejet_toxprior")] <- NA_real_
for (i in seq_len(nrow(tab_probas))) {
  
  cat(paste0("Ligne ", i, " / 22 :\n"))
  
  cat("Génération des patients.\n")
  # Génération de la cohorte principale
  essais <- gen_patients_multinomTOP(n_sim = 10000,
                                     ana_inter = c(30, 51),
                                     ana_inter_tox = c(rep(5, 4), rep(10, 5), 11),
                                     interpatient = 6, max_teff = 180, max_ttox = 42, rand_ratio = 1,
                                     multinom_ttt = list(as.numeric(tab_probas$probas[[i]])),
                                     seed = 1024)
  # Generation de la 2ème cohorte qui est cMET négative
  set.seed(1024)
  cohorte_2 <- lapply(
    X = essais,
    FUN = function(x) {
      res <- data.frame(temps_recrutement = sort(runif(9, 0, max(x$temps_recrutement))))
      res$tox <- rbinom(9, 1, as.numeric(tab_probas$probas[[i]])[1] + as.numeric(tab_probas$probas[[i]])[3])
      res$temps_tox <- gen_temps_patient(n = 9, tmax = 42)
      res$temps_tox <- ifelse(res$tox == 1, res$temps_recrutement + res$temps_tox, 1e+12)
      res$temps_obstox <- res$temps_recrutement + 42
      return(res)
    }
  )
  
  cat("Réalisation des essais.\n")
  # TOP avec approche Ivanova au temps du TOP
  resultats_topiva <- map2_dfr(.x = essais,
                               .y = cohorte_2,
                               .f = ~ realisation_essai_topiva(ana_inter_eff = c(30, 51),
                                                               ana_inter_tox = c(rep(5, 4), rep(10, 6), 1),
                                                               tableau = .x,
                                                               cohorte = .y,
                                                               phitox = c(.15, .25), prior_eff = .15, prior_tox = c(1, 1), 
                                                               C_ = seuil_effseule[["C_"]], gamm = seuil_effseule[["gamma"]], critere_tox = .95,
                                                               interpatient = 6, max_teff = 180, max_ttox = 42),
                               .id = "essai")
  tab_probas$rejet_topiva[i] <- mean(resultats_topiva$decision == "Acceptation traitement")
  # BOP avec approche Ivanova donc observation complète à chaque analyse intermédiaire
  resultats_bopiva <- map2_dfr(.x = essais,
                               .y = cohorte_2,
                               .f = ~ realisation_essai_bopiva(ana_inter_eff = c(30, 51),
                                                               ana_inter_tox = c(rep(5, 4), rep(10, 6), 1),
                                                               tableau = .x,
                                                               cohorte = .y,
                                                               phitox = c(.15, .25), prior_eff = .15, prior_tox = c(1, 1), 
                                                               p_n = c(0, .15, 0, .85), 
                                                               C_ = seuil_effseule[["C_"]], gamm = seuil_effseule[["gamma"]], critere_tox = .95,
                                                               interpatient = 6, max_teff = 180, max_ttox = 42),
                               .id = "essai")
  tab_probas$rejet_bopiva[i] <- mean(resultats_bopiva$decision == "Acceptation traitement")
  # TOP classique avec 2 analyses intermédiaires
  resultats_topclas <- simu_essais(alpha = .05,
                                   ana_inter = c(30, 51),
                                   ana_inter_tox = c(30, 51),
                                   p_n = c(.08, .07, .22, .63),
                                   p_a = c(.107, .193, .093, .607),
                                   p_reel = as.numeric(tab_probas$probas[[i]]),
                                   liste_tableaux = essais,
                                   cut_seq = seuil_toxclassique[["C_"]],
                                   power_seq = seuil_toxclassique[["gamma"]],
                                   interpatient = 6, max_teff = 180, max_ttox = 42,
                                   stock_matrice = FALSE)
  tab_probas$rejet_classique[i] <- resultats_topclas["rejet_h0"]
  # TOP avec monitoring rapproché de la toxicité
  resultats_toprapp <- simu_essais(alpha = .05,
                                   ana_inter = c(30, 51),
                                   ana_inter_tox = c(rep(5, 4), rep(10, 5), 11),
                                   p_n = c(.08, .07, .22, .63),
                                   p_a = c(.107, .193, .093, .607),
                                   p_reel = as.numeric(tab_probas$probas[[i]]),
                                   liste_tableaux = essais,
                                   cut_seq = seuil_toxrapprochee[["C_"]],
                                   power_seq = seuil_toxrapprochee[["gamma"]],
                                   interpatient = 6, max_teff = 180, max_ttox = 42,
                                   stock_matrice = FALSE)
  tab_probas$rejet_toxrapp[i] <- resultats_toprapp["rejet_h0"]
  # TOP en incorporant les données de la cohorte 2 dans la toxicité
  resultats_toxprior <- simu_essais_prior(alpha = .05,
                                          ana_inter = c(30, 51),
                                          ana_inter_tox = c(rep(5, 4), rep(10, 5), 11),
                                          ana_par = 9,
                                          p_n = c(.08, .07, .22, .63),
                                          p_a = c(.107, .193, .093, .607),
                                          p_reel = as.numeric(tab_probas$probas[[i]]),
                                          liste_tableaux = essais,
                                          liste_tableaux_par = cohorte_2,
                                          cut_seq = seuil_toxrapprochee[["C_"]],
                                          power_seq = seuil_toxrapprochee[["gamma"]],
                                          interpatient = 6, max_teff = 180, max_ttox = 42,
                                          stock_matrice = FALSE)
  tab_probas$rejet_toxprior[i] <- resultats_toxprior["rejet_h0"]
  set.seed(NULL) # Reset de graine
  
}

save(tab_probas, file = "data/probas_h0h1.Rdata")


# IV/ Reprise des 10 scénarios avec un seuil constant pris pour le calcul des seuils ----

source("R/realisations_essais_metimgast.R", encoding = "utf-8")
load("data/scenar20210907.Rdata")
tableau_scenars <- tableau_scenars[, -c(2:9)]

tableau_results_Rvar <- tableau_scenars
tableau_results_Rvar[, c("nb_pts_topiva", "temps_etude_topiva", 
                    "accept_ttt_topiva", "arret_prema_topiva", "arret_fut_topiva", "arret_tox_topiva", "nb_eff_topiva", 
                    "nb_tox_topiva")] <- NA_real_
tableau_results_Rvar[, c("nb_pts_bopiva", "temps_etude_bopiva", 
                    "accept_ttt_bopiva", "arret_prema_bopiva", "arret_fut_bopiva", "arret_tox_bopiva", "nb_eff_bopiva", 
                    "nb_tox_bopiva")] <- NA_real_
tableau_results_Rvar[, c("moy_pts_classique", "moy_duree_classique", "rejet_h0_classique", "arret_prema_classique",
                    "arret_fut_classique", "arret_tox_classique", "nb_eff_classique", 
                    "nb_tox_classique")] <- NA_real_
tableau_results_Rvar[, c("moy_pts_toxrapp", "moy_duree_toxrapp", "rejet_h0_toxrapp", "arret_prema_toxrapp",
                    "arret_fut_toxrapp", "arret_tox_toxrapp", "nb_eff_toxrapp", 
                    "nb_tox_toxrapp")] <- NA_real_
tableau_results_Rvar[, c("moy_pts_toxprior", "moy_duree_toxprior", "rejet_h0_toxprior", "arret_prema_toxprior",
                    "arret_fut_toxprior", "arret_tox_toxprior", "nb_eff_toxprior", 
                    "nb_tox_toxprior")] <- NA_real_
tableau_results_Rvar[, c("analyse_med_topiva", "analyse_med_bopiva")] <- NA_character_

for (i in seq_len(nrow(tableau_results_Rvar))) {
  
  cat(paste0(tableau_results_Rvar$scenar[i], "/", tableau_results_Rvar$nom_corr[i], " : ", i, " / ", nrow(tableau_results_Rvar), "\n"))
  
  p_reel <- tableau_results_Rvar$probas[i] %>% pluck(1) %>% unname()
  
  essais <- gen_patients_multinomTOP(n_sim = 10000,
                                     ana_inter = c(30, 51),
                                     ana_inter_tox = c(rep(5, 4), rep(10, 5), 11),
                                     interpatient = 6, max_teff = 180, max_ttox = 42, rand_ratio = 1,
                                     multinom_ttt = list(p_reel),
                                     seed = 1993)
  # Generation de la 2ème cohorte qui est cMET négative
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
  
  # TOP avec approche Ivanova au temps du TOP
  resultats_topiva <- map2_dfr(.x = essais,
                               .y = cohorte_2,
                               .f = ~ realisation_essai_topiva(ana_inter_eff = c(30, 51),
                                                               ana_inter_tox = c(rep(5, 4), rep(10, 6), 1),
                                                               tableau = .x,
                                                               cohorte = .y,
                                                               phitox = c(.15, .25), prior_eff = .15, prior_tox = c(1, 1), 
                                                               C_ = seuil_effseule[["C_"]], gamm = seuil_effseule[["gamma"]], critere_tox = .95,
                                                               interpatient = 6, max_teff = 180, max_ttox = 42),
                               .id = "essai")
  tableau_results_Rvar[i, c("nb_pts_topiva", "temps_etude_topiva", 
                       "arret_fut_topiva", "arret_tox_topiva", "nb_eff_topiva", 
                       "nb_tox_topiva", "accept_ttt_topiva", "arret_prema_topiva", "analyse_med_topiva")] <- 
    resultats_topiva %>% 
    summarise(across(c(nb_pts, temps_etude, arret_fut:nb_tox), ~ mean(.x)),
              accept_ttt = mean(decision == "Acceptation traitement"),
              arret_prema = mean(str_detect(decision, "prématuré")),
              analyse_med_topiva = paste0("med = ", median(analyse), " / ", 
                                          paste(names(table(.$analyse)), "=", round(100 * prop.table(table(.$analyse)), 1), "%", collapse = " - "))) %>% 
    `colnames<-`(paste0(colnames(.), "_topiva")) %>% 
    tibble_row()
  
  # BOP avec approche Ivanova donc observation complète à chaque analyse intermédiaire
  resultats_bopiva <- map2_dfr(.x = essais,
                               .y = cohorte_2,
                               .f = ~ realisation_essai_bopiva(ana_inter_eff = c(30, 51),
                                                               ana_inter_tox = c(rep(5, 4), rep(10, 6), 1),
                                                               tableau = .x,
                                                               cohorte = .y,
                                                               phitox = c(.15, .25), prior_eff = .15, prior_tox = c(1, 1), 
                                                               p_n = tableau_results_Rvar$probas[tableau_results_Rvar$scenar == "H0" & tableau_results_Rvar$nom_corr == tableau_results_Rvar$nom_corr[i]] %>% pluck(1) %>% unname(), 
                                                               C_ = seuil_effseule[["C_"]], gamm = seuil_effseule[["gamma"]], critere_tox = .95,
                                                               interpatient = 6, max_teff = 180, max_ttox = 42),
                               .id = "essai")
  tableau_results_Rvar[i, c("nb_pts_bopiva", "temps_etude_bopiva", 
                       "arret_fut_bopiva", "arret_tox_bopiva", "nb_eff_bopiva", 
                       "nb_tox_bopiva", "accept_ttt_bopiva", "arret_prema_bopiva", "analyse_med_bopiva")] <- 
    resultats_bopiva %>% 
    summarise(across(c(nb_pts, temps_etude, arret_fut:nb_tox), ~ mean(.x)),
              accept_ttt = mean(decision == "Acceptation traitement"),
              arret_prema = mean(str_detect(decision, "prématuré")),
              analyse_med_bopiva = paste0("med = ", median(analyse), " / ", 
                                          paste(names(table(.$analyse)), "=", round(100 * prop.table(table(.$analyse)), 1), "%", collapse = " - "))) %>% 
    `colnames<-`(paste0(colnames(.), "_bopiva")) %>% 
    tibble_row()
  
  # TOP classique avec 2 analyses intermédiaires
  resultats_topclas <- simu_essais(alpha = .05,
                                   ana_inter = c(30, 51),
                                   ana_inter_tox = c(30, 51),
                                   p_n = tableau_results_Rvar$probas[tableau_results_Rvar$scenar == "H0" & tableau_results_Rvar$nom_corr == tableau_results_Rvar$nom_corr[i]] %>% pluck(1) %>% unname(),
                                   p_a = tableau_results_Rvar$probas[tableau_results_Rvar$scenar == "H1" & tableau_results_Rvar$nom_corr == tableau_results_Rvar$nom_corr[i]] %>% pluck(1) %>% unname(),
                                   p_reel = p_reel,
                                   liste_tableaux = essais,
                                   cut_seq = seuil_toxclassique[["C_"]],
                                   power_seq = seuil_toxclassique[["gamma"]],
                                   interpatient = 6, max_teff = 180, max_ttox = 42,
                                   stock_matrice = FALSE)
  tableau_results_Rvar[i, c("moy_pts_classique", "moy_duree_classique", "rejet_h0_classique", "arret_prema_classique",
                       "arret_fut_classique", "arret_tox_classique", "nb_eff_classique", 
                       "nb_tox_classique")] <- 
    resultats_topclas %>%
    `names<-`(paste0(names(resultats_topclas), "_classique")) %>% 
    as_tibble_row()
  
  # TOP avec monitoring rapproché de la toxicité
  resultats_toprapp <- simu_essais(alpha = .05,
                                   ana_inter = c(30, 51),
                                   ana_inter_tox = c(rep(5, 4), rep(10, 5), 11),
                                   p_n = tableau_results_Rvar$probas[tableau_results_Rvar$scenar == "H0" & tableau_results_Rvar$nom_corr == tableau_results_Rvar$nom_corr[i]] %>% pluck(1) %>% unname(),
                                   p_a = tableau_results_Rvar$probas[tableau_results_Rvar$scenar == "H1" & tableau_results_Rvar$nom_corr == tableau_results_Rvar$nom_corr[i]] %>% pluck(1) %>% unname(),
                                   p_reel = p_reel,
                                   liste_tableaux = essais,
                                   cut_seq = seuil_toxrapprochee[["C_"]],
                                   power_seq = seuil_toxrapprochee[["gamma"]],
                                   interpatient = 6, max_teff = 180, max_ttox = 42,
                                   stock_matrice = FALSE)
  tableau_results_Rvar[i, c("moy_pts_toxrapp", "moy_duree_toxrapp", "rejet_h0_toxrapp", "arret_prema_toxrapp",
                       "arret_fut_toxrapp", "arret_tox_toxrapp", "nb_eff_toxrapp", 
                       "nb_tox_toxrapp")] <- 
    resultats_toprapp %>% 
    `names<-`(paste0(names(resultats_toprapp), "_toxrapp")) %>% 
    as_tibble_row()
  
  # TOP en incorporant les données de la cohorte 2 dans la toxicité
  resultats_toxprior <- simu_essais_prior(alpha = .05,
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
                                          stock_matrice = FALSE)
  
  tableau_results_Rvar[i, c("moy_pts_toxprior", "moy_duree_toxprior", "rejet_h0_toxprior", "arret_prema_toxprior",
                       "arret_fut_toxprior", "arret_tox_toxprior", "nb_eff_toxprior", 
                       "nb_tox_toxprior")] <- 
    resultats_toxprior %>% 
    `names<-`(paste0(names(resultats_toxprior), "_toxprior")) %>% 
    as_tibble_row()
  
}

save(tableau_results_Rvar, file = paste0("data/resultats_tox30_rposvar_", format(Sys.Date(), "%Y%m%d"), ".Rdata"))
rm(tableau_results_Rvar, resultats_toxprior, resultats_toprapp, resultats_topclas, resultats_bopiva, resultats_topiva, cohorte_2, essais, p_reel)


# V/ Reprise des scénarios sous H0 et H1 avec 3 longueurs de fenêtre de toxicité différentes ----

# Pour le TOP uniquement pour l'efficacité
load("data/cutoff_planifie.Rdata")
seuil_effseule <- seuil[[1]]
rm(seuil)

# Pour le TOP pour efficacité et toxicité
# Récupération dans le même temps des différents scénarios
# On ne se place que pour une corrélation à minima légèrement positive
load("data/resultats_tox30_20210907.Rdata")
seuil_toxclassique <- tableau_scenars %>% 
  filter(scenar == "H0", nom_corr == "Rpos1") %>% 
  select(C__classique, gamma_classique, alpha_calc_classique, puissance_calc_classique) %>% 
  as.numeric() %>% 
  'names<-'(c("C_", "gamma", "alpha_calc", "puissance_calc"))
seuil_toxrapprochee <- tableau_scenars %>% 
  filter(scenar == "H0", nom_corr == "Rpos1") %>% 
  select(C__toxrapp, gamma_toxrapp, alpha_calc_toxrapp, puissance_calc_toxrapp) %>% 
  as.numeric() %>% 
  'names<-'(c("C_", "gamma", "alpha_calc", "puissance_calc"))
liste_scenars <- tableau_scenars %>% 
  filter(nom_corr == "Rpos1") %>% 
  select(scenar, probas, EffTox, EffNotox, NoeffTox, NoeffNotox, R)
rm(tableau_scenars)

liste_scenars <- liste_scenars %>% 
  filter(scenar %in% c("H0", "H1")) %>% 
  mutate(temps_tox = list(c(42, 84, 126), c(42, 84, 126))) %>% 
  unnest(temps_tox)

source("R/realisations_essais_metimgast.R", encoding = "utf-8")

tableau_results <- liste_scenars
tableau_results[, c("nb_pts_topiva", "temps_etude_topiva", 
                    "accept_ttt_topiva", "arret_prema_topiva", "arret_fut_topiva", "arret_tox_topiva", "nb_eff_topiva", 
                    "nb_tox_topiva")] <- NA_real_
tableau_results[, c("nb_pts_bopiva", "temps_etude_bopiva", 
                    "accept_ttt_bopiva", "arret_prema_bopiva", "arret_fut_bopiva", "arret_tox_bopiva", "nb_eff_bopiva", 
                    "nb_tox_bopiva")] <- NA_real_
tableau_results[, c("moy_pts_classique", "moy_duree_classique", "rejet_h0_classique", "arret_prema_classique",
                    "arret_fut_classique", "arret_tox_classique", "nb_eff_classique", 
                    "nb_tox_classique")] <- NA_real_
tableau_results[, c("moy_pts_toxrapp", "moy_duree_toxrapp", "rejet_h0_toxrapp", "arret_prema_toxrapp",
                    "arret_fut_toxrapp", "arret_tox_toxrapp", "nb_eff_toxrapp", 
                    "nb_tox_toxrapp")] <- NA_real_
tableau_results[, c("moy_pts_toxprior", "moy_duree_toxprior", "rejet_h0_toxprior", "arret_prema_toxprior",
                    "arret_fut_toxprior", "arret_tox_toxprior", "nb_eff_toxprior", 
                    "nb_tox_toxprior")] <- NA_real_

tableau_results[, c("analyse_med_topiva", "analyse_med_bopiva")] <- NA_character_

for (i in seq_len(nrow(liste_scenars))) {
  
  cat(paste0(liste_scenars$scenar[i], " : ", i, " / ", nrow(liste_scenars), "\n"))
  
  p_reel <- liste_scenars$probas[i] %>% pluck(1) %>% unname()
  
  essais <- gen_patients_multinomTOP(n_sim = 10000,
                                     ana_inter = c(30, 51),
                                     ana_inter_tox = c(rep(5, 4), rep(10, 5), 11),
                                     interpatient = 6, max_teff = 180, max_ttox = liste_scenars$temps_tox[i], rand_ratio = 1,
                                     multinom_ttt = list(p_reel),
                                     seed = 1993)
  # Generation de la 2ème cohorte qui est cMET négative
  cohorte_2 <- lapply(
    X = essais,
    FUN = function(x) {
      res <- data.frame(temps_recrutement = sort(runif(9, 0, max(x$temps_recrutement))))
      res$tox <- rbinom(9, 1, p_reel[1] + p_reel[3])
      res$temps_tox <- gen_temps_patient(n = 9, tmax = liste_scenars$temps_tox[i])
      res$temps_tox <- ifelse(res$tox == 1, res$temps_recrutement + res$temps_tox, 1e+12)
      res$temps_obstox <- res$temps_recrutement + liste_scenars$temps_tox[i]
      return(res)
    }
  )
  
  # TOP avec approche Ivanova au temps du TOP
  resultats_topiva <- map2_dfr(.x = essais,
                               .y = cohorte_2,
                               .f = ~ realisation_essai_topiva(ana_inter_eff = c(30, 51),
                                                               ana_inter_tox = c(rep(5, 4), rep(10, 6), 1),
                                                               tableau = .x,
                                                               cohorte = .y,
                                                               phitox = c(.15, .25), prior_eff = .15, prior_tox = c(1, 1), 
                                                               C_ = seuil_effseule[["C_"]], gamm = seuil_effseule[["gamma"]], critere_tox = .95,
                                                               interpatient = 6, max_teff = 180, max_ttox = liste_scenars$temps_tox[i]),
                               .id = "essai")
  tableau_results[i, c("nb_pts_topiva", "temps_etude_topiva", 
                       "arret_fut_topiva", "arret_tox_topiva", "nb_eff_topiva", 
                       "nb_tox_topiva", "accept_ttt_topiva", "arret_prema_topiva", "analyse_med_topiva")] <- 
    resultats_topiva %>% 
    summarise(across(c(nb_pts, temps_etude, arret_fut:nb_tox), ~ mean(.x)),
              accept_ttt = mean(decision == "Acceptation traitement"),
              arret_prema = mean(str_detect(decision, "prématuré")),
              analyse_med_topiva = paste0("med = ", median(analyse), " / ", 
                                          paste(names(table(.$analyse)), "=", round(100 * prop.table(table(.$analyse)), 1), "%", collapse = " - "))) %>% 
    `colnames<-`(paste0(colnames(.), "_topiva")) %>% 
    tibble_row()
  
  # BOP avec approche Ivanova donc observation complète à chaque analyse intermédiaire
  resultats_bopiva <- map2_dfr(.x = essais,
                               .y = cohorte_2,
                               .f = ~ realisation_essai_bopiva(ana_inter_eff = c(30, 51),
                                                               ana_inter_tox = c(rep(5, 4), rep(10, 6), 1),
                                                               tableau = .x,
                                                               cohorte = .y,
                                                               phitox = c(.15, .25), prior_eff = .15, prior_tox = c(1, 1), 
                                                               p_n = liste_scenars$probas[liste_scenars$scenar == "H0"] %>% pluck(1) %>% unname(), 
                                                               C_ = seuil_effseule[["C_"]], gamm = seuil_effseule[["gamma"]], critere_tox = .95,
                                                               interpatient = 6, max_teff = 180, max_ttox = liste_scenars$temps_tox[i]),
                               .id = "essai")
  tableau_results[i, c("nb_pts_bopiva", "temps_etude_bopiva", 
                       "arret_fut_bopiva", "arret_tox_bopiva", "nb_eff_bopiva", 
                       "nb_tox_bopiva", "accept_ttt_bopiva", "arret_prema_bopiva", "analyse_med_bopiva")] <- 
    resultats_bopiva %>% 
    summarise(across(c(nb_pts, temps_etude, arret_fut:nb_tox), ~ mean(.x)),
              accept_ttt = mean(decision == "Acceptation traitement"),
              arret_prema = mean(str_detect(decision, "prématuré")),
              analyse_med_bopiva = paste0("med = ", median(analyse), " / ", 
                                          paste(names(table(.$analyse)), "=", round(100 * prop.table(table(.$analyse)), 1), "%", collapse = " - "))) %>% 
    `colnames<-`(paste0(colnames(.), "_bopiva")) %>% 
    tibble_row()
  
  # TOP classique avec 2 analyses intermédiaires
  resultats_topclas <- simu_essais(alpha = .05,
                                   ana_inter = c(30, 51),
                                   ana_inter_tox = c(30, 51),
                                   p_n = liste_scenars$probas[liste_scenars$scenar == "H0"] %>% pluck(1) %>% unname(),
                                   p_a = liste_scenars$probas[liste_scenars$scenar == "H1"] %>% pluck(1) %>% unname(),
                                   p_reel = p_reel,
                                   liste_tableaux = essais,
                                   cut_seq = seuil_toxclassique[["C_"]],
                                   power_seq = seuil_toxclassique[["gamma"]],
                                   interpatient = 6, max_teff = 180, max_ttox = liste_scenars$temps_tox[i],
                                   stock_matrice = FALSE)
  tableau_results[i, c("moy_pts_classique", "moy_duree_classique", "rejet_h0_classique", "arret_prema_classique",
                       "arret_fut_classique", "arret_tox_classique", "nb_eff_classique", 
                       "nb_tox_classique")] <- 
    resultats_topclas %>%
    `names<-`(paste0(names(resultats_topclas), "_classique")) %>% 
    as_tibble_row()
  
  # TOP avec monitoring rapproché de la toxicité
  resultats_toprapp <- simu_essais(alpha = .05,
                                   ana_inter = c(30, 51),
                                   ana_inter_tox = c(rep(5, 4), rep(10, 5), 11),
                                   p_n = liste_scenars$probas[liste_scenars$scenar == "H0"] %>% pluck(1) %>% unname(),
                                   p_a = liste_scenars$probas[liste_scenars$scenar == "H1"] %>% pluck(1) %>% unname(),
                                   p_reel = p_reel,
                                   liste_tableaux = essais,
                                   cut_seq = seuil_toxrapprochee[["C_"]],
                                   power_seq = seuil_toxrapprochee[["gamma"]],
                                   interpatient = 6, max_teff = 180, max_ttox = liste_scenars$temps_tox[i],
                                   stock_matrice = FALSE)
  tableau_results[i, c("moy_pts_toxrapp", "moy_duree_toxrapp", "rejet_h0_toxrapp", "arret_prema_toxrapp",
                       "arret_fut_toxrapp", "arret_tox_toxrapp", "nb_eff_toxrapp", 
                       "nb_tox_toxrapp")] <- 
    resultats_toprapp %>% 
    `names<-`(paste0(names(resultats_toprapp), "_toxrapp")) %>% 
    as_tibble_row()
  
  # TOP en incorporant les données de la cohorte 2 dans la toxicité
  resultats_toxprior <- simu_essais_prior(alpha = .05,
                                          ana_inter = c(30, 51),
                                          ana_inter_tox = c(rep(5, 4), rep(10, 5), 11),
                                          ana_par = 9,
                                          p_n = liste_scenars$probas[liste_scenars$scenar == "H0"] %>% pluck(1) %>% unname(),
                                          p_a = liste_scenars$probas[liste_scenars$scenar == "H1"] %>% pluck(1) %>% unname(),
                                          p_reel = p_reel,
                                          liste_tableaux = essais,
                                          liste_tableaux_par = cohorte_2,
                                          cut_seq = seuil_toxrapprochee[["C_"]],
                                          power_seq = seuil_toxrapprochee[["gamma"]],
                                          interpatient = 6, max_teff = 180, max_ttox = liste_scenars$temps_tox[i],
                                          stock_matrice = FALSE)
  
  tableau_results[i, c("moy_pts_toxprior", "moy_duree_toxprior", "rejet_h0_toxprior", "arret_prema_toxprior",
                       "arret_fut_toxprior", "arret_tox_toxprior", "nb_eff_toxprior", 
                       "nb_tox_toxprior")] <- 
    resultats_toxprior %>% 
    `names<-`(paste0(names(resultats_toxprior), "_toxprior")) %>% 
    as_tibble_row()
  
}

save(tableau_results, file = paste0("data/resultats_tox30_rpos_fentox_", format(Sys.Date(), "%Y%m%d"), ".Rdata"))
rm(tableau_results, resultats_toxprior, resultats_toprapp, resultats_topclas, resultats_bopiva, resultats_topiva, cohorte_2, essais, p_reel)


# VI/ Risque alpha et puissance lorsque les seuils sont calculés pour une efficacité et une toxicité indépendantes ----

## Retrouver les seuils calculés

load("data/cutoff_planifie.Rdata")
seuil_effseule <- seuil[[1]]
rm(seuil)

# Pour le TOP pour efficacité et toxicité
# Récupération dans le même temps des différents scénarios
# On ne se place que pour une corrélation à minima légèrement positive
load("data/resultats_tox30_20210907.Rdata")
seuil_toxclassique <- tableau_scenars %>% 
  filter(scenar == "H0", nom_corr == "Rind") %>% 
  select(C__classique, gamma_classique, alpha_calc_classique, puissance_calc_classique) %>% 
  as.numeric() %>% 
  'names<-'(c("C_", "gamma", "alpha_calc", "puissance_calc"))
seuil_toxrapprochee <- tableau_scenars %>% 
  filter(scenar == "H0", nom_corr == "Rind") %>% 
  select(C__toxrapp, gamma_toxrapp, alpha_calc_toxrapp, puissance_calc_toxrapp) %>% 
  as.numeric() %>% 
  'names<-'(c("C_", "gamma", "alpha_calc", "puissance_calc"))
liste_scenars <- tableau_scenars %>% 
  filter(nom_corr == "Rind") %>% 
  select(scenar, probas, EffTox, EffNotox, NoeffTox, NoeffNotox, R)
rm(tableau_scenars)

## Estimation des risques alpha et puissances pour chaque schéma dans chaque situation

source("R/realisations_essais_metimgast.R", encoding = "utf-8")

generer_corr_multi <- function(vec_inter,
                               eff,
                               tox) {
  
  data.frame(Pet = vec_inter) %>% 
    mutate(probas = map(Pet, ~ c(Pefftox = .x, Peffnotox = eff - .x, Pnoefftox = tox - .x, Pnoeffnotox = 1 - eff - tox + .x)),
           correlation = map_dbl(probas, ~ (.x[1] - eff * tox) / sqrt((eff - eff * eff) * (tox - tox * tox))))
  
}

# Calculer tous les risques pour toutes les hypothèses possibles
tab_probas <- bind_rows(
  generer_corr_multi(seq(0, .15, .015), .15, .3) %>% mutate(hypothese = "H0"),
  generer_corr_multi(seq(0, .2, .02), .3, .2) %>% mutate(hypothese = "H1")
) 
tab_probas[, c("rejet_topiva", "rejet_bopiva", "rejet_classique", "rejet_toxrapp", "rejet_toxprior")] <- NA_real_
for (i in seq_len(nrow(tab_probas))) {
  
  cat(paste0("Ligne ", i, " / 22 :\n"))
  
  cat("Génération des patients.\n")
  # Génération de la cohorte principale
  essais <- gen_patients_multinomTOP(n_sim = 10000,
                                     ana_inter = c(30, 51),
                                     ana_inter_tox = c(rep(5, 4), rep(10, 5), 11),
                                     interpatient = 6, max_teff = 180, max_ttox = 42, rand_ratio = 1,
                                     multinom_ttt = list(as.numeric(tab_probas$probas[[i]])),
                                     seed = 1024)
  # Generation de la 2ème cohorte qui est cMET négative
  set.seed(1024)
  cohorte_2 <- lapply(
    X = essais,
    FUN = function(x) {
      res <- data.frame(temps_recrutement = sort(runif(9, 0, max(x$temps_recrutement))))
      res$tox <- rbinom(9, 1, as.numeric(tab_probas$probas[[i]])[1] + as.numeric(tab_probas$probas[[i]])[3])
      res$temps_tox <- gen_temps_patient(n = 9, tmax = 42)
      res$temps_tox <- ifelse(res$tox == 1, res$temps_recrutement + res$temps_tox, 1e+12)
      res$temps_obstox <- res$temps_recrutement + 42
      return(res)
    }
  )
  
  cat("Réalisation des essais.\n")
  # TOP avec approche Ivanova au temps du TOP
  resultats_topiva <- map2_dfr(.x = essais,
                               .y = cohorte_2,
                               .f = ~ realisation_essai_topiva(ana_inter_eff = c(30, 51),
                                                               ana_inter_tox = c(rep(5, 4), rep(10, 6), 1),
                                                               tableau = .x,
                                                               cohorte = .y,
                                                               phitox = c(.15, .25), prior_eff = .15, prior_tox = c(1, 1), 
                                                               C_ = seuil_effseule[["C_"]], gamm = seuil_effseule[["gamma"]], critere_tox = .95,
                                                               interpatient = 6, max_teff = 180, max_ttox = 42),
                               .id = "essai")
  tab_probas$rejet_topiva[i] <- mean(resultats_topiva$decision == "Acceptation traitement")
  # BOP avec approche Ivanova donc observation complète à chaque analyse intermédiaire
  resultats_bopiva <- map2_dfr(.x = essais,
                               .y = cohorte_2,
                               .f = ~ realisation_essai_bopiva(ana_inter_eff = c(30, 51),
                                                               ana_inter_tox = c(rep(5, 4), rep(10, 6), 1),
                                                               tableau = .x,
                                                               cohorte = .y,
                                                               phitox = c(.15, .25), prior_eff = .15, prior_tox = c(1, 1), 
                                                               p_n = c(0, .15, 0, .85), 
                                                               C_ = seuil_effseule[["C_"]], gamm = seuil_effseule[["gamma"]], critere_tox = .95,
                                                               interpatient = 6, max_teff = 180, max_ttox = 42),
                               .id = "essai")
  tab_probas$rejet_bopiva[i] <- mean(resultats_bopiva$decision == "Acceptation traitement")
  # TOP classique avec 2 analyses intermédiaires
  resultats_topclas <- simu_essais(alpha = .05,
                                   ana_inter = c(30, 51),
                                   ana_inter_tox = c(30, 51),
                                   p_n = c(.08, .07, .22, .63),
                                   p_a = c(.107, .193, .093, .607),
                                   p_reel = as.numeric(tab_probas$probas[[i]]),
                                   liste_tableaux = essais,
                                   cut_seq = seuil_toxclassique[["C_"]],
                                   power_seq = seuil_toxclassique[["gamma"]],
                                   interpatient = 6, max_teff = 180, max_ttox = 42,
                                   stock_matrice = FALSE)
  tab_probas$rejet_classique[i] <- resultats_topclas["rejet_h0"]
  # TOP avec monitoring rapproché de la toxicité
  resultats_toprapp <- simu_essais(alpha = .05,
                                   ana_inter = c(30, 51),
                                   ana_inter_tox = c(rep(5, 4), rep(10, 5), 11),
                                   p_n = c(.08, .07, .22, .63),
                                   p_a = c(.107, .193, .093, .607),
                                   p_reel = as.numeric(tab_probas$probas[[i]]),
                                   liste_tableaux = essais,
                                   cut_seq = seuil_toxrapprochee[["C_"]],
                                   power_seq = seuil_toxrapprochee[["gamma"]],
                                   interpatient = 6, max_teff = 180, max_ttox = 42,
                                   stock_matrice = FALSE)
  tab_probas$rejet_toxrapp[i] <- resultats_toprapp["rejet_h0"]
  # TOP en incorporant les données de la cohorte 2 dans la toxicité
  resultats_toxprior <- simu_essais_prior(alpha = .05,
                                          ana_inter = c(30, 51),
                                          ana_inter_tox = c(rep(5, 4), rep(10, 5), 11),
                                          ana_par = 9,
                                          p_n = c(.08, .07, .22, .63),
                                          p_a = c(.107, .193, .093, .607),
                                          p_reel = as.numeric(tab_probas$probas[[i]]),
                                          liste_tableaux = essais,
                                          liste_tableaux_par = cohorte_2,
                                          cut_seq = seuil_toxrapprochee[["C_"]],
                                          power_seq = seuil_toxrapprochee[["gamma"]],
                                          interpatient = 6, max_teff = 180, max_ttox = 42,
                                          stock_matrice = FALSE)
  tab_probas$rejet_toxprior[i] <- resultats_toxprior["rejet_h0"]
  set.seed(NULL) # Reset de graine
  
}

save(tab_probas, file = "data/probas_h0h1_Rind.Rdata")
