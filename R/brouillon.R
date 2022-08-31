library(multibrasBOP2)
library(tidyverse)
seuil <- deter_cutoff(alpha = .05,
                      ana_inter = c(30, 51),
                      p_n = c(.1, .05, .15, .7),
                      p_a = c(.12, .18, .03, .67),
                      mat_beta_xi = matrix(c(1, 1, 0, 0, 0, 1, 0, 1), nrow = 2, byrow = TRUE),
                      methode = 3L,
                      affich_mat = "Non")
seuil[[1]]
# Pas exactement le même gamma, mais avec le même donne les mêmes chiffres. Et en plus ça donne les mêmes chiffres que l'application.

test_pat <- gen_patients_multinomTOP(n_sim = 10000,
                                     ana_inter = c(30, 51),
                                     interpatient = 6, max_teff = 180, max_ttox = 180,
                                     multinom_ttt = c(.1, .05, .15, .7),
                                     mat_beta_xi = matrix(c(1, 1, 0, 0, 0, 1, 0, 1), nrow = 2, byrow = TRUE))
test_pat[[1]]
ana_inter <- c(30, 51)

mat_beta_xi <- matrix(c(1, 1, 0, 0, 0, 1, 0, 1), nrow = 2, byrow = TRUE)
p_n <- c(.1, .05, .15, .7)
phitox <- c(sum(p_n * mat_beta_xi[1, ]), 1 - sum(p_n * mat_beta_xi[2, ]))

tableau <- test_pat[[1]]
ana_inter_cum <- cumsum(ana_inter)
nmax <- sum(ana_inter)

PhiEff <- phitox[1]
PhiTox <- phitox[2]
prior_eff <- sum(p_n * mat_beta_xi[1, ])
prior_tox <- 1 - sum(p_n * mat_beta_xi[2, ])

C_ <- .69
gamm <- .88
max_teff <- 180
max_ttox <- 180
interpatient <- 6

decision <- NULL

for (i in ana_inter_cum) {
  
  if (i != nmax) {
    temps_decision <- tableau$temps_recrutement[i + 1] # Observation at the time when we would enroll the next patient
  } else {
    temps_decision <- max(tableau$temps_obseff[i], tableau$temps_obstox[i]) # Observation at the end of the trial if it goes that far
  }
  obs_eff        <- (tableau$temps_eff[seq_len(i)] <= temps_decision) | (tableau$temps_obseff[seq_len(i)] <= temps_decision) # Boolean of observation of response criterion
  obs_tox        <- (tableau$temps_tox[seq_len(i)] <= temps_decision) | (tableau$temps_obstox[seq_len(i)] <= temps_decision) # Boolean of observation of toxicity criterion
  evt_eff        <- tableau$temps_eff[seq_len(i)] <= temps_decision # Boolean for response events
  evt_tox        <- tableau$temps_tox[seq_len(i)] <= temps_decision # Boolean for toxicity events
  yeff           <- sum(evt_eff) # Number of responses
  ytox           <- sum(evt_tox) # Number of toxicities
  PP_eff         <- pbeta(phitox[1], prior_eff + yeff, 1 - prior_eff + i - yeff)
  PP_tox         <- 1 - pbeta(phitox[2], prior_tox + ytox, 1 - prior_tox + i - ytox)
  seuil          <- 1 - C_ * (i / sum(ana_inter)) ^ gamm
  addtime        <- 0
    
  # Suspension rules: when not enough responses to be sure of efficacy or not enough toxicities to be sure about toxicity
  if (PP_tox <= seuil) { # Case where no stopping for toxicity 
    
    if (PP_eff > seuil) { # Case where stopping for futility
      
      temps_decision_bis <- 0
      
      while ((1 - sum(obs_eff) / i >= i / nmax) || (1 - sum(obs_tox) / i >= i / nmax)) {
        # More than n / N pending patients so wait until enough non pending patients
        addtime <- addtime + interpatient
        temps_decision_bis <- temps_decision + addtime
        obs_eff <- (tableau$temps_eff[seq_len(i)] <= temps_decision_bis) | (tableau$temps_obseff[seq_len(i)] <= temps_decision_bis)
        obs_tox <- (tableau$temps_tox[seq_len(i)] <= temps_decision_bis) | (tableau$temps_obstox[seq_len(i)] <= temps_decision_bis)
        
      }
      
    } else { # Case where no stopping for futility
      
      while ((1 - sum(obs_tox) / i >= i / nmax)) {
        # More than n / N pending patientsfor toxicity so wait until enough non pending patients
        # We are here sure that if more responses occured we still won't stop for futility
        addtime <- addtime + interpatient
        temps_decision_bis <- temps_decision + addtime
        obs_tox <- (tableau$temps_tox[seq_len(i)] <= temps_decision_bis) | (tableau$temps_obstox[seq_len(i)] <= temps_decision_bis)
        
      }
      
    }
    
  }
  
  if (i != nmax) tableau[(i + 1):nmax, 10:14] <- tableau[(i + 1):nmax, 10:14] + addtime
  if (i != nmax) temps_decision <- tableau$temps_recrutement[i + 1]
  
  obs_eff        <- (tableau$temps_eff[seq_len(i)] <= temps_decision) | (tableau$temps_obseff[seq_len(i)] <= temps_decision)
  obs_tox        <- (tableau$temps_tox[seq_len(i)] <= temps_decision) | (tableau$temps_obstox[seq_len(i)] <= temps_decision)
  evt_eff        <- tableau$temps_eff[seq_len(i)] <= temps_decision
  evt_tox        <- tableau$temps_tox[seq_len(i)] <= temps_decision
  yeff           <- sum(evt_eff)
  ytox           <- sum(evt_tox)
  suivi_eff      <- (temps_decision - tableau$temps_recrutement[seq_len(i)][!obs_eff]) / max_teff
  suivi_tox      <- (temps_decision - tableau$temps_recrutement[seq_len(i)][!obs_tox]) / max_ttox
  TESS_eff       <- sum(obs_eff) + sum(suivi_eff)
  TESS_tox       <- sum(obs_tox) + sum(suivi_tox)
  PP_eff         <- pbeta(phitox[1], prior_eff + yeff, 1 - prior_eff + TESS_eff - yeff)
  PP_tox         <- 1 - pbeta(phitox[2], prior_tox + ytox, 1 - prior_tox + TESS_tox - ytox)
  
  if (PP_eff > seuil || PP_tox > seuil) {
    if (i == nmax) {decision <- "Non acceptation traitement"} else {decision <- "Non acceptation prématurée"}
    break
  }
  if (i == nmax) decision <- "Acceptation traitement"
  
}

data.frame(
  id_essai = n,
  analyse = match(i, ana_inter_cum), nb_pts = i, temps_etude = temps_decision,
  decision = decision, arret_fut = as.numeric(PP_eff > seuil), arret_tox = as.numeric(PP_tox > seuil),
  nb_eff = yeff, nb_tox = ytox,
  nb_eff_true = sum(tableau$eff[seq_len(i)]), nb_tox_true = sum(tableau$tox[seq_len(i)]),
  C_ = C_, Gamma = gamm,
  alpha_eff = prior_eff + yeff, beta_eff = 1 - prior_eff + TESS_eff - yeff,
  alpha_tox = prior_tox + ytox, beta_tox = 1 - prior_tox + TESS_tox - ytox
)


ggplot() +
  map(seq(5, 25, 2), function(x) geom_function(aes(color = paste0("TESS: ", x)), fun = ~ pbeta(.x, 1 + 3, 1 + x - 3))) +
  xlim(c(0, 1)) +
  scale_color_discrete(breaks = paste0("TESS: ", seq(5, 25, 2)))

