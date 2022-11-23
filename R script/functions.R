# --------------------------------------------- #
# Rscript for functions usefull in the analysis #
# --------------------------------------------- #


# The different used packages in the functions :
# library(tidyverse)
# library(patchwork)
# library(ggtext)


# I/ Correlation between efficacy and toxicity ----

# Generate the 4 probabilities and the correlation between efficacy 
# and toxicity given p(eff), p(tox) and p(eff & tox)
generer_corr_multi <- function(vec_inter,
                               eff,
                               tox) {
  
  data.frame(Pet = vec_inter) %>% 
    mutate(probas = map(Pet, ~ c(Pefftox = .x, Peffnotox = eff - .x, Pnoefftox = tox - .x, Pnoeffnotox = 1 - eff - tox + .x)),
           correlation = map_dbl(probas, ~ (.x[1] - eff * tox) / sqrt((eff - eff * eff) * (tox - tox * tox))))
  
}


# Compute correlation between efficacy and toxicity given p(eff), p(tox) and p(eff & tox)
calculate_corr <- function(Peff, Ptox, Pinter) {
  (Pinter - Peff * Ptox) / sqrt((Peff - Peff ^ 2) * (Ptox - Ptox ^ 2))
}


# Compute p(eff & tox) given p(eff), p(tox) and correlation between efficacy and toxicity
calculate_pinter <- function(Peff, Ptox, corr) {
  corr * sqrt((Peff - Peff ^ 2) * (Ptox - Ptox ^ 2)) + Peff * Ptox
}


# Generate the 6 multinomial laws with different correlation between efficacy 
# and toxicity as described in the article
generer_6corr_hyp <- function(Peff, Ptox, nom_scenar = NULL) {
  
  if (is.null(nom_scenar)) nom_scenar <- paste0("Eff", round(Peff, 2), "/Tox", round(Ptox, 2))
  
  Pefftox  <- Peff * Ptox
  vec_Rind <- c(Pefftox, Peff - Pefftox, Ptox - Pefftox, 1 - Peff - Ptox + Pefftox)
  
  Rmin     <- calculate_corr(Peff, Ptox, max(0, Peff + Ptox - 1))
  vec_Rmin <- c(max(0, Peff + Ptox - 1), Peff, Ptox, 1 - Peff - Ptox)
  
  Rmax     <- calculate_corr(Peff, Ptox, min(Peff, Ptox))
  vec_Rmax <- c(min(Peff, Ptox), Peff - min(Peff, Ptox), Ptox - min(Peff, Ptox), 1 - Peff - Ptox + min(Peff, Ptox))
  
  Rneg       <- Rmin / 2
  Pinter_neg <- calculate_pinter(Peff, Ptox, Rneg)
  vec_Rneg   <- c(Pinter_neg, Peff - Pinter_neg, Ptox - Pinter_neg, 1 - Peff - Ptox + Pinter_neg)
  
  Rpos1       <- Rmax / 3
  Pinter_pos1 <- calculate_pinter(Peff, Ptox, Rpos1)
  vec_Rpos1   <- c(Pinter_pos1, Peff - Pinter_pos1, Ptox - Pinter_pos1, 1 - Peff - Ptox + Pinter_pos1)
  
  Rpos2       <- 2 * Rmax / 3
  Pinter_pos2 <- calculate_pinter(Peff, Ptox, Rpos2)
  vec_Rpos2   <- c(Pinter_pos2, Peff - Pinter_pos2, Ptox - Pinter_pos2, 1 - Peff - Ptox + Pinter_pos2)
  
  resultats <- tibble(
    scenar = nom_scenar,
    nom_corr = c("Rmin", "Rneg", "Rind", "Rpos1", "Rpos2", "Rmax"),
    R = c(Rmin, Rneg, 0, Rpos1, Rpos2, Rmax),
    probas = map(list(vec_Rmin, vec_Rneg, vec_Rind, vec_Rpos1, vec_Rpos2, vec_Rmax), ~ setNames(.x, c("EffTox", "EffNotox", "NoeffTox", "NoeffNotox")))
  ) %>% 
    mutate(proba_temp = probas) %>% 
    unnest_wider(proba_temp)
  
  return(resultats)
  
}


# II/ PPtox designs ----

# Stopping rules for toxicity with PPtox design
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

# 
decision_essai <- function(liste, tab_arret) {
  liste %>% 
    bind_rows(.id = "data") %>% 
    left_join(tab_arret, by = c("n" = "n_pat")) %>% 
    group_by(data) %>% 
    mutate(decision = case_when(ntox > tox_max & row_number() < n() ~ "Early stopping",
                                ntox > tox_max & row_number() == n() ~ "Stopping",
                                ntox <= tox_max ~ "Accept treatment",
                                TRUE ~ NA_character_),
           bool = str_detect(decision, "[S|s]topping"),
           sum_bool = sum(bool),
           keep = ifelse(sum_bool == 0, vec_eval_tox[length(vec_eval_tox)], first(n[bool]))) %>% 
    ungroup() %>% 
    filter(n == keep) %>% 
    select(data:decision)
}


## A/ TOP + PPtox ----

# Results for 1 TOP+PPtox trial
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
  
  # Toxicity dataset 
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

# Get results of the simulated trials in liste_essais (generated with package multibrasBOP2)
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


## B/ BOP + PPtox ----

# Results for a single trial
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

# Results of all simulated trials stored in liste_essai (generated with package multibrasBOP2)
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


## C/ Simon + PPtox ----

# Compute the operating characteristics of a Simon design with fixed r1, n1, r and n
# Loop over possible values to determine optimal and minimax design
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

# Results of all simulated trials stored in liste_essai (generated with package multibrasBOP2)
simu_simon <- function(ana_inter_eff,
                       ana_inter_tox = NULL,
                       liste_essais, 
                       liste_cohortes,
                       alpha_simon = .05, puissance_simon = .9, pu, pa,
                       phitox, critere_tox, prior_tox) {
  
  # Stopping rules for efficacy
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


# III/ Other functions used in the work ----

# Compute the minimum number of patients satisfying type I error risk (alpha) and power (seuil_puissance)
# with threshold for exact binomial comparison.
# Plot the alpha and power for the optimal rule for each number of patients until nmax.
nsn_binom <- function(nmax, alpha, seuil_puissance, p0, p1) {
  
  resultats <- data.frame(
    nb_pat = seq_len(nmax),
    succes = NA,
    typeI = NA,
    puissance = NA
  )
  
  for (i in seq_len(nmax)) {
    
    for (j in seq_len(i + 1) - 1) {
      
      prop <- 1 - pbinom(j, i, p0)
      
      if (prop <= alpha) {
        
        resultats[i, 2] <- j
        resultats[i, 3] <- 1 - pbinom(j, i, p0)
        resultats[i, 4] <- 1 - pbinom(j, i, p1)
        break
        
      }
      
    }
    
  }
  
  plot <- ggplot(resultats %>% 
                   group_by(succes) %>% 
                   mutate(succes = if_else(row_number() == 1, paste0("> ", succes, " succès"), "")),
                 aes(x = nb_pat)) +
    geom_line(aes(color = "&alpha;", y = typeI)) +
    geom_line(aes(color = "1 - &beta;", y = puissance)) +
    geom_point(aes(color = "&alpha;", y = typeI)) +
    geom_point(aes(color = "1 - &beta;", y = puissance)) +
    geom_hline(yintercept = seuil_puissance, color = "#630929", linetype = "dashed") +
    geom_hline(yintercept = alpha, color = "#630929", linetype = "dashed") +
    geom_text(aes(y = 1.2, label = succes), angle = 90, hjust = 1) +
    scale_y_continuous(label = scales::percent_format(), breaks = seq(0, 1, .25)) +
    labs(x = "Nombre de patients",
         y = "&alpha; ou 1 - &beta;",
         colour = "Risque représenté") +
    theme_light() +
    theme(legend.text = element_markdown(),
          axis.title.y = element_markdown())
  print(plot)
  
  resultats_opti <- resultats %>% 
    filter(puissance >= seuil_puissance) %>% 
    slice_head(n = 1)
  
  return(list(resultats = resultats, resultats_opti = resultats_opti))
  
}


# Find parameters of beta law satisfying a specified mean (moyenne) and width of confidence interval (largeur_ic)
trouver_beta_ic <- function(moyenne, largeur_ic, conf.level = .9, dec = 2) {
  
  vec_ic <- c((1 - conf.level) / 2, (1 + conf.level) / 2)
  soluce <- uniroot(f = function(x) diff(qbeta(vec_ic, x, x * (1 - moyenne) / moyenne)) - largeur_ic,
                    interval = c(0, 100))
  a <- round(soluce$root, dec)
  b <- a * (1 - moyenne) / moyenne
  
  return(list(alpha = a, beta = b))
  
}

# Plot the operating characteristics for the 10 scenarios
plot_opchar <- function(tableau_scenars_plot, 
                        metadata, 
                        var,
                        niveaux = c("topiva", "bopiva", "simon", "classique", "toxrapp", "toxprior"),
                        etiquettes = c("TOP<sub>eff</sub> + PP<sub>tox</sub>", "BOP<sub>eff</sub> + PP<sub>tox</sub>",
                                     "Simon + PP<sub>tox</sub>", "TOP<sub>eff/tox</sub> with 2 analyses", 
                                     "TOP<sub>eff/tox</sub> with close monitoring of toxicity",
                                     "TOP<sub>eff/tox</sub> with informative prior")) {
  
  # Mise au format long
  tableau_res <- tableau_scenars_plot %>% 
    select(!matches("C_")) %>% 
    select(!matches("gamma")) %>% 
    select(-c(probas)) %>% 
    select(!starts_with("analyse_med")) %>% 
    select(scenar, scenar_p, scenarios, probas_plot, R, everything()) %>% 
    pivot_longer(cols = -c(scenar:R),
                 names_pattern = "(.+_.+)_(.+)",
                 names_to = c("stat", "schemas")) %>% 
    mutate(schemas = factor(schemas, levels = niveaux, labels = etiquettes))
  
  layout <- 
    "AB
     CC"
  
  # Barplot
  plot1 <- tableau_res %>% 
    filter(stat == var) %>% 
    ggplot(aes(x = scenar_p, y = value, fill = schemas)) +
    geom_col(position = position_dodge2(), color = "black") +
    geom_vline(xintercept = c(seq(1.5, 9.5, 1)), color = "grey50", 
               size = 1, linetype = "dashed") +
    coord_flip() +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.major.x = element_line(colour = "darkgrey", size = .7),
          panel.border = element_rect(color = "transparent"),
          legend.position = "bottom",
          axis.title.x = element_text(size = 20),
          axis.text.x = element_text(size = 16),
          legend.text = element_markdown(size = 20),
          legend.title = element_text(size = 24, face = "bold"),
          plot.margin = margin(t = 0, r = 0, b = 0, l = 0))  +
    guides(fill = guide_legend(nrow = 2)) +
    scale_y_continuous(name = metadata$label[match(var, metadata$var)],
                       limits = c(0, metadata$lim_sup[match(var, metadata$var)]),
                       labels = eval(metadata$format_legende[[match(var, metadata$var)]]),
                       breaks = eval(metadata$breaks[[match(var, metadata$var)]])) +
    scale_fill_manual(values = c("#648FFF", "#785EF0", "#DC267F", "#FE6100", "#FFB000", "#22482F")) +
    labs(x = NULL,
         fill = "Design")
  legende <- get_legend(plot1)
  
  # Annotations pour les scénarios
  plot2 <- tableau_scenars_plot %>% 
    ggplot() +
    geom_richtext(aes(x = scenar_p, y = 1.1, label = scenarios), hjust = 1, fill = "transparent", label.color = "transparent", size = 7.5) +
    geom_richtext(aes(x = scenar_p, y = 1.15, label = probas_plot), hjust = 0, fill = "transparent", label.color = "transparent", size = 7.5) +
    geom_vline(xintercept = c(seq(1.5, 9.5, 1)), color = "grey50", 
               size = 1, linetype = "dashed") +
    coord_flip(clip = "off") +
    ylim(c(.8, 1.7)) +
    theme_void() +
    theme(strip.background = element_rect(fill = NA, color = NA),
          strip.text = element_textbox(
            size = 12, 
            color = "transparent", fill = "transparent", box.color = "transparent",
            halign = 0.5, linetype = 1, r = unit(3, "pt"), width = unit(.001, "in"),
            padding = margin(2, 0, 1, 0), margin = margin(3, 3, 3, 3)),
          plot.margin = margin(r = 0))
  
  plot_tot <- plot2 + (plot1 + theme(legend.position = "none")) + legende +
    plot_layout(design = layout, widths = c(1.5, 7), heights = c(10, 1))
  
  return(plot_tot)
  
}

# Same as previous function but for all correlations
plot_opchar_facet <- function(tableau_scenars_plot, metadata, var) {
  
  # Mise au format long
  tableau_res <- tableau_scenars_plot %>% 
    select(!matches("C_")) %>% 
    select(!starts_with("analyse_med")) %>% 
    select(!matches("gamma")) %>% 
    select(-c(probas, EffTox, EffNotox, NoeffTox, NoeffNotox)) %>% 
    select(scenar, scenar_p, scenarios, probas_plot, R, nom_corr, everything()) %>% 
    pivot_longer(cols = -c(scenar:nom_corr),
                 names_pattern = "(.+_.+)_(.+)",
                 names_to = c("stat", "schemas")) %>% 
    mutate(schemas = factor(schemas, 
                            levels = c("topiva", "bopiva", "simon", "classique", "toxrapp", "toxprior"),
                            labels = c("TOP<sub>eff</sub> + PP<sub>tox</sub>", "BOP<sub>eff</sub> + PP<sub>tox</sub>", "Simon + PP<sub>tox</sub>", "TOP<sub>eff/tox</sub> with 2 analyses", "TOP<sub>eff/tox</sub> with close monitoring of toxicity", "TOP<sub>eff/tox</sub> with informative prior")))
  
  layout <- 
    "AB
     CC"
  
  # Barplot
  plot1 <- tableau_res %>% 
    filter(stat == var) %>% 
    ggplot(aes(x = scenar_p, y = value, fill = schemas)) +
    geom_col(position = position_dodge2(), color = "black") +
    geom_vline(xintercept = c(seq(1.5, 9.5, 1)), color = "grey50", 
               size = 1, linetype = "dashed") +
    coord_flip() +
    facet_wrap(vars(nom_corr), ncol = 2, dir = "v") +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.major.x = element_line(colour = "darkgrey", size = .7),
          panel.border = element_rect(color = "transparent"),
          legend.position = "bottom",
          axis.title.x = element_text(size = 20),
          axis.text.x = element_text(size = 16),
          legend.text = element_markdown(size = 20),
          legend.title = element_text(size = 24, face = "bold"),
          plot.margin = margin(t = 0, r = 0, b = 0, l = 0))  +
    guides(fill = guide_legend(nrow = 2)) +
    scale_y_continuous(name = metadata$label[match(var, metadata$var)],
                       limits = c(0, metadata$lim_sup[match(var, metadata$var)]),
                       labels = eval(metadata$format_legende[[match(var, metadata$var)]]),
                       breaks = eval(metadata$breaks[[match(var, metadata$var)]])) +
    scale_fill_manual(values = c("#648FFF", "#785EF0", "#DC267F", "#FE6100", "#FFB000", "#22482F")) +
    labs(x = NULL,
         fill = "Design")
  legende <- get_legend(plot1)
  
  # Annotations pour les scénarios
  plot2 <- tableau_scenars_plot %>% 
    filter(nom_corr %in% c("R<sub>min</sub>", "R<sub>neg</sub>", "R<sub>ind</sub>")) %>% 
    ggplot() +
    geom_richtext(aes(x = scenar_p, y = 1.1, label = scenarios), hjust = 1, fill = "transparent", label.color = "transparent", size = 7.5) +
    geom_richtext(aes(x = scenar_p, y = 1.15, label = probas_plot), hjust = 0, fill = "transparent", label.color = "transparent", size = 7.5) +
    geom_vline(xintercept = c(seq(1.5, 9.5, 1)), color = "grey50", 
               size = 1, linetype = "dashed") +
    coord_flip(clip = "off") +
    ylim(c(.8, 1.7)) +
    theme_void() +
    facet_wrap(vars(nom_corr), ncol = 2, nrow = 3, dir = "v") +
    theme(strip.background = element_rect(fill = NA, color = NA),
          strip.text = element_textbox(
            size = 16, 
            color = "transparent", fill = "transparent", box.color = "transparent",
            halign = 0.5, linetype = 1, r = unit(3, "pt"), width = unit(.001, "in"),
            padding = margin(2, 0, 1, 0), margin = margin(3, 3, 3, 3)),
          plot.margin = margin(r = 0))
  
  plot_tot <- plot2 + (plot1 + theme(legend.position = "none")) + legende +
    plot_layout(design = layout, widths = c(1.5, 7), heights = c(10, 1))
  
  return(plot_tot)
  
}
