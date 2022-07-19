# ---------------------------------------------------------------------- #
# Script de simulation des données pour évaluer les différentes méthodes #
# Créé le 27/08/2021, modifié le 27/08/2021                              #
# ---------------------------------------------------------------------- #

# Packages et helpers ----

# devtools::install_github("GuillaumeMulier/multibrasBOP2", auth_token = "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")
library(multibrasBOP2)

# I/ Génération des listes d'essais ----

## A/ Simulation settings ----

n_simu       <- 10000
ana_inter    <- c(30, 51)
interpatient <- 6
max_teff     <- 180
max_ttox     <- 180
seed         <- 1993
liste_proba  <- list(
  H0_Rmin    = c(0, .15, .2, .65),
  H0_Rind    = c(.03, .12, .17, .68),
  H0_Rpos    = c(.08, .07, .12, .73),
  H0_Rmax    = c(.15, 0, .05, .8),
  H1_30_Rmin = c(0, .3, .3, .4),
  H1_30_Rneg = c(.05, .25, .25, .45),
  H1_30_Rind = c(.09, .21, .21, .49),
  H1_30_Rpos = c(.15, .15, .15, .55),
  H1_30_Rmax = c(.3, 0, 0, .7),
  H1_40_Rmin = c(0, .3, .4, .3),
  H1_40_Rneg = c(.05, .25, .35, .35),
  H1_40_Rind = c(.12, .18, .28, .42),
  H1_40_Rpos = c(.15, .15, .25, .45),
  H1_40_Rmax = c(.3, 0, .1, .6)
)
noms_liste <- names(liste_proba)

## B/ Simulation ----

essais_H0H1 <- gen_patients_multinomTOP(n_sim = n_simu,
                                        ana_inter = ana_inter,
                                        interpatient = interpatient,
                                        max_teff = max_teff,
                                        max_ttox = max_ttox,
                                        seed = seed,
                                        multinom_ttt = liste_proba)

## C/ Stockage des résultats ----

save(essais_H0H1, file = paste0("data/essais_H0H1_", format(Sys.Date(), "%Y%m%d"), ".Rdata"))

## D/ Ménage ----

rm(essais_H0H1)








map_dfr(noms_liste, 
        function(x) simu_essais(ana_inter = c(30, 51), 
                                p_n = c(0, .15, 0, .85), 
                                p_a = c(0, .3, 0, .7), 
                                p_reel = liste_proba[[x]], 
                                cut_seq = .92, 
                                power_seq = .97, 
                                stock_matrice = TRUE,
                                mat_beta_xi = matrix(c(1, 1, 0, 0, 0, 0, 0, 0), nrow = 2, byrow = TRUE), 
                                liste_tableaux = lapply(essais_H0H1, function(tab) tab[tab$ttt == x, ]), 
                                interpatient = 6, 
                                max_teff = 180, 
                                max_ttox = 180)$op_char, 
        .id = "hyp")
