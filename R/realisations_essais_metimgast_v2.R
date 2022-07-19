# ------------------------------------------------------------------------------- #
# Script de fonctions pour réaliser l'essai selon ce qui est prévu pour METIMGAST #
# Créé le 02/09/2021, modifé le 27/01/2022                                        #
# ------------------------------------------------------------------------------- #

# Règles d'arrêt selon Ivanova pour la toxicité ----

regle_arret <- function(prior, ana_inter, seuil, critere) {
  
  map_dfr(
    .x = ana_inter,
    .f = function(x) {
      vec_x <- seq_len(x + 1) - 1
      proba <- map_dbl(vec_x, ~ 1 - pbeta(seuil, prior[1] + .x, prior[2] + x - .x))
      proba <- proba <= critere
      return(c(n_pat = x, tox_max = sum(proba) - 1))
    }
  )
  
}


# Réalisation de l'essai au temps du TOP ----

realisation_essai_topiva <- function(ana_inter_eff,
                                     ana_inter_tox = NULL,
                                     tableau, tab_tox,
                                     cohorte,
                                     phieff, prior_eff,
                                     C_, gamm,
                                     interpatient, max_teff) {
  
  if (!is.null(ana_inter_tox)) {
    ana_eff         <- cumsum(ana_inter_eff)
    ana_tox         <- cumsum(ana_inter_tox)
    ana_inter_cum   <- sort(union(ana_eff, ana_tox))
    ana_inter       <- c(ana_inter_cum[1], diff(ana_inter_cum))
  } else {
    ana_eff         <- cumsum(ana_inter_eff)
    ana_tox         <- cumsum(ana_inter_eff)
    ana_inter_cum   <- cumsum(ana_inter_eff)
    ana_inter       <- ana_inter_eff
  }
  nmax <- max(ana_inter_cum)
  decision <- NULL
  
  # # Really bad code but didn't solve the problem of sample size at the end
  # tab_tox <- regle_arret(prior = prior_tox,
  #                        ana_inter = ana_tox %>% '[<-'(11, 90),
  #                        seuil = phitox[2],
  #                        critere = critere_tox)
  
  # Dataset pour la toxicité (tableau reste le dataset pour l'efficacité)
  donnees_tox <- rbind(tableau, 
        cbind(cohorte, ttt = NA, V1 = NA, V2 = NA, V3 = NA, V4 = NA, nb_pat = NA, analyse = NA, eff = NA, temps_eff = NA, temps_obseff = NA)) %>% 
    arrange(temps_recrutement)
  
  for (i in ana_inter_cum) { # 81 patients for efficacy is the end of the study and 90 patients for toxicity will be recruited at that time
    
    if (i != nmax) {
      temps_decision_eff <- tableau$temps_recrutement[i + 1] # Observation at the time when we would enroll the next patient
      temps_decision_tox <- donnees_tox$temps_recrutement[i + 1] 
    } else {
      temps_decision_eff <- temps_decision_tox <- max(tableau$temps_obseff[i], tableau$temps_obstox[i], donnees_tox$temps_obstox[90]) # Observation at the end of the trial if it goes that far
    }
    
    obs_eff        <- (tableau$temps_eff[seq_len(i)] <= temps_decision_eff) | (tableau$temps_obseff[seq_len(i)] <= temps_decision_eff) # Boolean of observation of response criterion
    evt_eff        <- tableau$temps_eff[seq_len(i)] <= temps_decision_eff # Boolean for response events
    yeff           <- sum(evt_eff) # Number of responses
    PP_eff         <- pbeta(phieff, prior_eff + yeff, 1 - prior_eff + i - yeff)
    seuil          <- 1 - C_ * (i / sum(ana_inter)) ^ gamm
    addtime        <- 0
    
    if (i %in% ana_eff && PP_eff > seuil) {
      
      temps_decision_eff_bis <- 0
      
      while (1 - sum(obs_eff) / i >= i / nmax) {
        # More than n / N pending patients so wait until enough non pending patients
        addtime <- addtime + interpatient
        temps_decision_eff_bis <- temps_decision_eff + addtime
        obs_eff <- (tableau$temps_eff[seq_len(i)] <= temps_decision_eff_bis) | (tableau$temps_obseff[seq_len(i)] <= temps_decision_eff_bis)
        
      }
      
      if (i != nmax) tableau[(i + 1):nmax, 10:14] <- tableau[(i + 1):nmax, 10:14] + addtime
      if (i != nmax) temps_decision_eff <- tableau$temps_recrutement[i + 1]
      
      obs_eff        <- (tableau$temps_eff[seq_len(i)] <= temps_decision_eff) | (tableau$temps_obseff[seq_len(i)] <= temps_decision_eff)
      evt_eff        <- tableau$temps_eff[seq_len(i)] <= temps_decision_eff
      yeff           <- sum(evt_eff)
      suivi_eff      <- (temps_decision_eff - tableau$temps_recrutement[seq_len(i)][!obs_eff]) / max_teff
      TESS_eff       <- sum(obs_eff) + sum(suivi_eff)
      PP_eff         <- pbeta(phieff, prior_eff + yeff, 1 - prior_eff + TESS_eff - yeff)
      
    }
    
    if (i < nmax) nb_tox <- sum(donnees_tox$temps_tox[seq_len(i)] <= temps_decision_tox) else nb_tox <- sum(donnees_tox$temps_tox[seq_len(90)] <= temps_decision_tox)
    if (i %in% ana_tox) {
      if (i < nmax) seuil_tox <- tab_tox$tox_max[tab_tox$n_pat == i] else seuil_tox <- tab_tox$tox_max[tab_tox$n_pat == 90]
    } else {
      seuil_tox <- i + 1
    }
    
    if (PP_eff * (i %in% ana_eff) > seuil || nb_tox > seuil_tox) {
      if (i == nmax) {decision <- "Non acceptation traitement"} else {decision <- "Non acceptation prématurée"}
      break
    }
    if (i == nmax) decision <- "Acceptation traitement"
    
  }
  
  sortie <- data.frame(
    analyse     = match(i, ana_inter_cum),
    nb_pts      = i,
    temps_etude = if (PP_eff * (i %in% ana_eff) > seuil & nb_tox <= seuil_tox) temps_decision_eff else temps_decision_tox,
    decision    = decision,
    arret_fut   = if (i == 30 & nb_tox > seuil_tox & temps_decision_eff > temps_decision_tox) 0 else as.numeric(PP_eff * (i %in% ana_eff) > seuil),
    arret_tox   = as.numeric(nb_tox > seuil_tox),
    nb_eff      = yeff,
    nb_tox      = nb_tox
  )
  
  return(sortie)
  
}

simu_topiva <- function(ana_inter_eff,
                        ana_inter_tox = NULL,
                        liste_essais, 
                        liste_cohortes,
                        phieff, phitox, prior_eff, 
                        prior_tox, critere_tox,
                        C_, gamm,
                        interpatient, max_teff, max_ttox) {
  
  # Really bad code but didn't solve the problem of sample size at the end
  tab_tox <- regle_arret(prior = prior_tox,
                         ana_inter = cumsum(ana_inter_tox) %>% '[<-'(11, 90),
                         seuil = phitox,
                         critere = critere_tox)
  
  map2_dfr(.x = liste_essais,
           .y = liste_cohortes,
           .f = ~ realisation_essai_topiva(ana_inter_eff = ana_inter_eff,
                                           ana_inter_tox = ana_inter_tox,
                                           tableau = .x, tab_tox = tab_tox,
                                           cohorte = .y,
                                           phieff = phieff, prior_eff = prior_eff, 
                                           C_ = C_, gamm = gamm,
                                           interpatient = interpatient, max_teff = max_teff),
           .id = "essai") %>% 
    summarise(across(c(nb_pts, temps_etude, arret_fut:nb_tox), ~ mean(.x)),
              accept_ttt = mean(decision == "Acceptation traitement"),
              arret_prema = mean(str_detect(decision, "prématuré")),
              analyse_med = paste0("med = ", median(analyse), " / ", 
                                   paste(names(table(.$analyse)), "=", round(100 * prop.table(table(.$analyse)), 1), "%", collapse = " - ")))
  
}


# Réalisation de l'essai au temps du BOP ----

realisation_essai_bopiva <- function(ana_inter_eff,
                                     ana_inter_tox = NULL,
                                     tableau, tab_tox, tab_eff,
                                     cohorte) {
  
  if (!is.null(ana_inter_tox)) {
    ana_eff         <- cumsum(ana_inter_eff)
    ana_tox         <- cumsum(ana_inter_tox)
    ana_inter_cum   <- sort(union(ana_eff, ana_tox))
    ana_inter       <- c(ana_inter_cum[1], diff(ana_inter_cum))
  } else {
    ana_eff         <- cumsum(ana_inter_eff)
    ana_tox         <- cumsum(ana_inter_eff)
    ana_inter_cum   <- cumsum(ana_inter_eff)
    ana_inter       <- ana_inter_eff
  }
  nmax <- sum(ana_inter)
  decision <- NULL
  
  donnees_tox <- rbind(tableau, 
                       cbind(cohorte, ttt = NA, V1 = NA, V2 = NA, V3 = NA, V4 = NA, nb_pat = NA, analyse = NA, eff = NA, temps_eff = NA, temps_obseff = NA)) %>% 
    arrange(temps_recrutement)
  
  for (i in ana_inter_cum) {
    
    temps_decision_eff <- max(tableau$temps_obseff[i], tableau$temps_obstox[i]) # Observation at the end of the trial if it goes that far
    temps_decision_tox <- donnees_tox$temps_obstox[i]
    if (i != nmax) { # Update the time of the recruited patients
      if (i %in% ana_eff & i %in% ana_tox) {
        tableau[(i + 1):nmax, 10:14] <- tableau[(i + 1):nmax, 10:14] + max(temps_decision_eff, temps_decision_tox) - tableau$temps_recrutement[i + 1]
        donnees_tox[(i + 1):nmax, 10:14] <- donnees_tox[(i + 1):nmax, 10:14] + max(temps_decision_eff, temps_decision_tox) - donnees_tox$temps_recrutement[i + 1]
        temps_decision <- max(temps_decision_eff, temps_decision_tox)
      } else if (i %in% ana_eff) {
        tableau[(i + 1):nmax, 10:14] <- tableau[(i + 1):nmax, 10:14] + temps_decision_eff - tableau$temps_recrutement[i + 1]
        donnees_tox[(i + 1):nmax, 10:14] <- donnees_tox[(i + 1):nmax, 10:14] + temps_decision_eff - donnees_tox$temps_recrutement[i + 1]
        temps_decision <- temps_decision_eff
      } else if (i %in% ana_tox) {
        tableau[(i + 1):nmax, 10:14] <- tableau[(i + 1):nmax, 10:14] + temps_decision_tox - tableau$temps_recrutement[i + 1]
        donnees_tox[(i + 1):nmax, 10:14] <- donnees_tox[(i + 1):nmax, 10:14] + temps_decision_tox - donnees_tox$temps_recrutement[i + 1]
        temps_decision <- temps_decision_tox
      }
    } else {
      temps_decision <- max(tableau$temps_obseff[nmax], tableau$temps_obstox[nmax], donnees_tox$temps_obstox[90])
    }

    nb_eff <- sum(tableau$eff[seq_len(i)])
    if (i %in% ana_eff) seuil_eff <- tab_eff$seuil_eff[tab_eff$tot_pat == i] else seuil_eff <- -1
    
    if (i != nmax) nb_tox <- sum(donnees_tox$tox[seq_len(i)]) else nb_tox <- sum(donnees_tox$tox)
    if (i %in% ana_tox) {
      if (i != nmax) {
        seuil_tox <- tab_tox$tox_max[tab_tox$n_pat == i] 
      } else {
        seuil_tox <- tab_tox$tox_max[tab_tox$n_pat == 90] 
      }
    } else {
      seuil_tox <- i + 1
    }
      
    
    if (nb_eff <= seuil_eff || nb_tox > seuil_tox) {
      if (i == nmax) {decision <- "Non acceptation traitement"} else {decision <- "Non acceptation prématurée"}
      break
    }
    if (i == nmax) decision <- "Acceptation traitement"
    
  }
  
  sortie <- data.frame(
    analyse     = match(i, ana_inter_cum),
    nb_pts      = i,
    temps_etude = temps_decision,
    decision    = decision,
    arret_fut   = as.numeric(nb_eff <= seuil_eff),
    arret_tox   = as.numeric(nb_tox > seuil_tox),
    nb_eff      = nb_eff,
    nb_tox      = nb_tox
  )
  
  return(sortie)
  
}

simu_bopiva <- function(ana_inter_eff,
                        ana_inter_tox = NULL,
                        liste_essais, 
                        liste_cohortes,
                        phitox, 
                        prior_tox, critere_tox, p_n,
                        C_, gamm) {
  
  # Really bad code but didn't solve the problem of sample size at the end
  tab_tox <- regle_arret(prior = prior_tox,
                         ana_inter = cumsum(ana_inter_tox) %>% '[<-'(11, 90),
                         seuil = phitox,
                         critere = critere_tox)
  
  tab_eff <- get_stopbound(ana_inter = ana_inter_eff,
                           C_ = C_, gamm = gamm,
                           p_n = p_n) %>% 
    suppressWarnings() %>% suppressMessages()
  
  map2_dfr(.x = liste_essais,
           .y = liste_cohortes,
           .f = ~ realisation_essai_bopiva(ana_inter_eff = ana_inter_eff,
                                           ana_inter_tox = ana_inter_tox,
                                           tableau = .x, tab_tox = tab_tox, tab_eff = tab_eff,
                                           cohorte = .y),
           .id = "essai") %>% 
    summarise(across(c(nb_pts, temps_etude, arret_fut:nb_tox), ~ mean(.x)),
              accept_ttt = mean(decision == "Acceptation traitement"),
              arret_prema = mean(str_detect(decision, "prématuré")),
              analyse_med = paste0("med = ", median(analyse), " / ", 
                                   paste(names(table(.$analyse)), "=", round(100 * prop.table(table(.$analyse)), 1), "%", collapse = " - ")))
  
}

# Réalisation de l'essai avec le design de Simon pour l'efficacité au temps du BOP2 ----

alpha_puiss_simon <- function(r1, n1, r, n, pu, pa) {
  
  alpha1 <- pbinom(r1, n1, pu)
  alpha2 <- sum(map_dbl(seq(r1 + 1, min(n1, r)), ~ dbinom(.x, n1, pu) * pbinom(r - .x, n - n1, pu)))
  alpha <- 1 - (alpha1 + alpha2)
  
  puissance1 <- pbinom(r1, n1, pa)
  puissance2 <- sum(map_dbl(seq(r1 + 1, min(n1, r)), ~ dbinom(.x, n1, pa) * pbinom(r - .x, n - n1, pa)))
  puissance <- 1 - (puissance1 + puissance2)
  
  EN_p0 <- n1 + (1 - alpha1) * (n - n1)
  
  return(c("r1" = r1, "n1" = n1, "r" = r, "n" = n, "pu" = pu, "pa" = pa, "alpha" = alpha, "puissance" = puissance, "EN_p0" = EN_p0, "PET_po" = alpha1))
  
}

simu_simon <- function(ana_inter_eff,
                       ana_inter_tox = NULL,
                       liste_essais, 
                       liste_cohortes,
                       alpha_simon = .05, puissance_simon = .9, pu, pa,
                       phitox, critere_tox, prior_tox) {
  
  # Obtenir les règles d'arrêt de Simon
  donnees_test <- expand.grid(
    n1 = ana_inter_eff[1],
    r1 = seq(0, ana_inter_eff[1], 1),
    n = sum(ana_inter_eff),
    r = seq(0, ana_inter_eff[1], 1),
    pu = pu,
    pa = pa
  ) %>% 
    filter(r >= r1)
  result_simu_simon <- map_dfr(seq_len(nrow(donnees_test)), ~ alpha_puiss_simon(donnees_test$r1[.x], donnees_test$n1[.x], donnees_test$r[.x], donnees_test$n[.x], donnees_test$pu[.x], donnees_test$pa[.x]))
  result_simu_simon <- result_simu_simon %>% 
    filter(alpha <= alpha_simon, puissance >= puissance_simon) %>% 
    arrange(EN_p0) %>% 
    slice(1)
  if (nrow(result_simu_simon) == 0) stop("Aucun design de Simon ne remplit les conditions de risque alpha et de puissance aux effectifs spécifiés.", call. = FALSE)
  result_simu_simon <- setNames(object = as.numeric(result_simu_simon), nm = names(result_simu_simon))
  tab_eff <- data.frame(
    nb_ana = 1:2,
    tot_pat = cumsum(ana_inter_eff),
    seuil_eff = result_simu_simon[c("r1", "r")]
  )
  
  # Really bad code but didn't solve the problem of sample size at the end
  tab_tox <- regle_arret(prior = prior_tox,
                         ana_inter = cumsum(ana_inter_tox) %>% '[<-'(11, 90),
                         seuil = phitox,
                         critere = critere_tox)
  
  map2_dfr(.x = liste_essais,
           .y = liste_cohortes,
           .f = ~ realisation_essai_bopiva(ana_inter_eff = ana_inter_eff,
                                           ana_inter_tox = ana_inter_tox,
                                           tableau = .x, tab_tox = tab_tox, tab_eff = tab_eff,
                                           cohorte = .y),
           .id = "essai") %>% 
    summarise(across(c(nb_pts, temps_etude, arret_fut:nb_tox), ~ mean(.x)),
              accept_ttt = mean(decision == "Acceptation traitement"),
              arret_prema = mean(str_detect(decision, "prématuré")),
              analyse_med = paste0("med = ", median(analyse), " / ", 
                                   paste(names(table(.$analyse)), "=", round(100 * prop.table(table(.$analyse)), 1), "%", collapse = " - ")))
  
}
