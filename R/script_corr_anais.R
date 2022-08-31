library(tidyverse)

# Calculer la corrélation quand on connait Peff, Ptox et Peff+tox
calculate_corr <- function(Peff, Ptox, Pinter) {
  (Pinter - Peff * Ptox) / sqrt((Peff - Peff ^ 2) * (Ptox - Ptox ^ 2))
}

# Calculer Peff+tox quand on connait Peff, Ptox et la corrélation
calculate_pinter <- function(Peff, Ptox, corr) {
  corr * sqrt((Peff - Peff ^ 2) * (Ptox - Ptox ^ 2)) + Peff * Ptox
}

# Générer mes 6 hypothèses (paragraphe 3 du HTML).
# Pour adapter tu peux te servir de ce qu'il y a 
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

# Exemple : J'ai une probabilité de réponse de 60% et une probablité de toxicité de 40%
Peff <- .6
Ptox <- .4
(Rmin     <- calculate_corr(Peff, Ptox, max(0, Peff + Ptox - 1)))
(Rmax     <- calculate_corr(Peff, Ptox, min(Peff, Ptox)))
# Corrélation qui peut aller entre -1 et 0.67
# On veut avoir la loi multinomiale avec une corrélation de 0.3 entre efficacité et toxicité
(Pinter <- calculate_pinter(Peff, Ptox, .3))
(vec_multinom   <- c(Pinter, Peff - Pinter, Ptox - Pinter, 1 - Peff - Ptox + Pinter))
# Voilà les 4 paramètres de la loi multinomiale avec dans l'ordre : 
# (Peff+tox+, Peff+tox-, Peff-tox+, Peff-tox-)

# On peut essayer de vérifier (et j'espère que ça va marcher lol)
set.seed(121221)
tab_test <- rmultinom(1000, 1, vec_multinom) %>% 
  t() %>%
  as.data.frame() %>% 
  mutate(eff = V1 + V2,
         tox = V1 + V3)
cor(tab_test$eff, tab_test$tox)
# On retrouve environ ce qu'on a demandé, ouf !