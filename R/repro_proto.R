# ------------------------------------------------------------- #
# Script pour retrouver les chiffres de Lucie dans le protocole #
# et pour comparer les operating characteristics des designs    #
# Créé le 09/08/2021, modifié le 18/10/2021                     #
# ------------------------------------------------------------- #

# Packages et helpers ----

library(tidyverse)
library(ggtext)
# devtools::install_github("GuillaumeMulier/multibrasBOP2", auth_token = "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")
library(multibrasBOP2)
library(patchwork)

theme_set(theme_light(base_size = 13) +
            theme(strip.background = element_rect(fill = NA),
                  strip.text = element_textbox(
                    size = 14, 
                    color = "white", fill = "#7888C0", box.color = "#000066",
                    halign = 0.5, linetype = 1, r = unit(3, "pt"), width = unit(0.75, "npc"),
                    padding = margin(2, 0, 1, 0), margin = margin(3, 3, 3, 3)),
                  legend.title = element_markdown()))

# I/ Pour le nombre de sujets avec la loi binomiale ----

## A/ Fonction pour trouver le nombre de sujet nécessaires ----

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

nsn_binom(nmax = 200, alpha = .05, seuil_puissance = .8, p0 = .08, p1 = .2)$resultats_opti


## B/ Fonction des seuils de la toxicité ----

### Fonction pour trouver la loi beta selon l'IC spécifié

trouver_beta_ic <- function(moyenne, largeur_ic, conf.level = .9, dec = 2) {
  
  vec_ic <- c((1 - conf.level) / 2, (1 + conf.level) / 2)
  soluce <- uniroot(f = function(x) diff(qbeta(vec_ic, x, x * (1 - moyenne) / moyenne)) - largeur_ic,
                    interval = c(0, 100))
  a <- round(soluce$root, dec)
  b <- a * (1 - moyenne) / moyenne
  
  return(list(alpha = a, beta = b))
  
}

trouver_beta_ic(moyenne = .5, largeur_ic = .8)

### Fonction pour avoir le nombre max de toxicités avant d'arriver à l'arrêt
# On prend un critère de la forme Pr(Ptox>seuil|Dn)>critere

regle_arret <- function(prior, ana_inter, seuil, critere) {
  
  map_dfr(
    .x = ana_inter,
    .f = function(x) {
      vec_x <- seq_len(x + 1) - 1
      proba <- map_dbl(vec_x, ~ pbeta(seuil, prior[1] + .x, prior[2] + x - .x))
      proba <- proba >= critere
      return(c(n_pat = x, tox_max = sum(proba) - 1))
    }
  )
  
}

regle_arret(prior = c(1, 1), ana_inter = c(5, 10, 15, 20, 30, 40, 50, 60, 70, 81), seuil = .25, critere = .05)


## C/ TOP-design pour l'efficacité ----

seuil <- deter_cutoff(alpha = .05,
                      ana_inter = c(30, 51),
                      ana_inter_tox = c(81),
                      p_n = c(0, .15, 0, .85),
                      p_a = c(0, .3, 0, .7),
                      mat_beta_xi = matrix(c(1, 1, 0, 0, 0, 0, 0, 0), nrow = 2, byrow = TRUE),
                      methode = 3L,
                      affich_mat = "No")
seuil[[1]]
# Pas exactement le même gamma, mais avec le même donne les mêmes chiffres. Et en plus ça donne les mêmes chiffres que l'application.
save(seuil, file = paste0("data/cutoff_planifie.Rdata"))


# II/ Planification du TOP avec les settings de l'essai METIMGAST ----

## A/ Influence de la corrélation entre efficacité et toxicité sur un schéma TOP classique ----

generer_corr_multi <- function(vec_inter,
                               eff,
                               tox) {
  
  data.frame(Pet = vec_inter) %>% 
    mutate(probas = map(Pet, ~ c(Pefftox = .x, Peffnotox = eff - .x, Pnoefftox = tox - .x, Pnoeffnotox = 1 - eff - tox + .x)),
           correlation = map_dbl(probas, ~ (.x[1] - eff * tox) / sqrt((eff - eff * eff) * (tox - tox * tox))))
  
}

corr_1 <- generer_corr_multi(seq(0, .15, .015), .15, .4)
corr_2 <- generer_corr_multi(seq(0, .15, .015), .15, .3)
corr_3 <- generer_corr_multi(seq(0, .2, .02), .3, .2)
tab_probas <- expand.grid(pas_0 = seq(0, .15, .015), pas_1 = seq(0, .2, .02)) %>% 
  left_join(corr_1 %>% `colnames<-`(paste0(colnames(.), "_H0")), by = c("pas_0" = "Pet_H0")) %>% 
  left_join(corr_3 %>% `colnames<-`(paste0(colnames(.), "_H1")), by = c("pas_1" = "Pet_H1"))

liste_resultats <- list()
for (i in seq_len(nrow(tab_probas))) {
  liste_resultats[[i]] <- deter_cutoff(alpha = .05,
                                       ana_inter = c(30, 51),
                                       ana_inter_tox = c(30, 51),
                                       p_n = unlist(tab_probas$probas_H0[i]),
                                       p_a = unlist(tab_probas$probas_H1[i]),
                                       mat_beta_xi = matrix(c(1, 1, 0, 0, 0, 1, 0, 1), nrow = 2, byrow = TRUE),
                                       methode = 3L,
                                       seed = 1993,
                                       affich_mat = "No")[[1]]
}
table_resume <- cbind(tab_probas, do.call(rbind, liste_resultats))
save(table_resume, file = paste0("data/essais_cutoff_", format(Sys.Date(), "%Y%m%d"), ".Rdata"))
# load("data/essais_cutoff_20210829.Rdata")

heatmap_top <- function(data, colonne) {
  labelleur <- c("C_" = "&lambda;", "gamma" = "&gamma;", "alpha_calc" = "Planned &alpha;", "puissance_calc" = "Planned power")
  data %>% 
    pivot_longer(c(C_, gamma, alpha_calc, puissance_calc)) %>% 
    filter(name == colonne) %>% 
    ggplot(aes(x = correlation_H0, y = correlation_H1, fill = value, label = value)) +
    geom_tile() +
    geom_text(color = "white", size = 4) +
    facet_wrap(vars(name), labeller = labeller(name = labelleur)) +
    scale_x_continuous(name = expression(paste("Correlation Eff/Tox under ", H[0])), expand = c(0, 0)) +
    scale_y_continuous(name = expression(paste("Correlation Eff/Tox under ", H[1])), expand = c(0, 0)) +
    labs(fill = labs(fill = labelleur[colonne]))
}
heatmap_seuils_40 <- map(c("C_", "gamma", "alpha_calc", "puissance_calc"), heatmap_top, data = table_resume) %>% 
  wrap_plots()
ggsave(plot = heatmap_seuils_40,
       filename = "figures/brouillon/heatmap_seuils_40.png",
       device = "png", height = 8, width = 16)

tab_probas <- expand.grid(pas_0 = seq(0, .15, .015), pas_1 = seq(0, .2, .02)) %>% 
  left_join(corr_2 %>% `colnames<-`(paste0(colnames(.), "_H0")), by = c("pas_0" = "Pet_H0")) %>% 
  left_join(corr_3 %>% `colnames<-`(paste0(colnames(.), "_H1")), by = c("pas_1" = "Pet_H1"))
liste_resultats <- list()
for (i in seq_len(nrow(tab_probas))) {
  liste_resultats[[i]] <- deter_cutoff(alpha = .05,
                                       ana_inter = c(30, 51),
                                       ana_inter_tox = c(30, 51),
                                       p_n = unlist(tab_probas$probas_H0[i]),
                                       p_a = unlist(tab_probas$probas_H1[i]),
                                       mat_beta_xi = matrix(c(1, 1, 0, 0, 0, 1, 0, 1), nrow = 2, byrow = TRUE),
                                       methode = 3L,
                                       seed = 1993,
                                       affich_mat = "No")[[1]]
}
table_resume_2 <- cbind(tab_probas, do.call(rbind, liste_resultats))
save(table_resume_2, file = paste0("data/essais_cutoff_bis_", format(Sys.Date(), "%Y%m%d"), ".Rdata"))
# load("data/essais_cutoff_bis_20210830.Rdata")

heatmap_seuils_30 <- map(c("C_", "gamma", "alpha_calc", "puissance_calc"), heatmap_top, data = table_resume_2) %>% 
  wrap_plots()
ggsave(plot = heatmap_seuils_30,
       filename = "figures/brouillon/heatmap_seuils_30eng.png",
       device = "png", height = 8, width = 16)

rm(liste_resultats, corr_1, corr_2, corr_3, tab_probas, table_resume, table_resume_2, heatmap_seuils_30, heatmap_seuils_40)


## B/ Influence de la corrélation entre efficacité et toxicité lors d'un TOP design avec analyses rapprochées pour toxicité ----

corr_1 <- generer_corr_multi(seq(0, .15, .015), .15, .4)
corr_2 <- generer_corr_multi(seq(0, .15, .015), .15, .3)
corr_3 <- generer_corr_multi(seq(0, .2, .02), .3, .2)

tab_probas <- expand.grid(pas_0 = seq(0, .15, .015), pas_1 = seq(0, .2, .02)) %>% 
  left_join(corr_1 %>% `colnames<-`(paste0(colnames(.), "_H0")), by = c("pas_0" = "Pet_H0")) %>% 
  left_join(corr_3 %>% `colnames<-`(paste0(colnames(.), "_H1")), by = c("pas_1" = "Pet_H1"))

liste_resultats <- list()
for (i in seq_len(nrow(tab_probas))) {
  liste_resultats[[i]] <- deter_cutoff(alpha = .05,
                                       ana_inter = c(30, 51),
                                       ana_inter_tox = c(rep(10, 7), 11),
                                       p_n = unlist(tab_probas$probas_H0[i]),
                                       p_a = unlist(tab_probas$probas_H1[i]),
                                       mat_beta_xi = matrix(c(1, 1, 0, 0, 0, 1, 0, 1), nrow = 2, byrow = TRUE),
                                       methode = 3L,
                                       seed = 1993,
                                       affich_mat = "No")[[1]]
  
}
table_resume_tox <- cbind(tab_probas, do.call(rbind, liste_resultats))
save(table_resume_tox, file = paste0("data/essais_cutoff_toxrapp_", format(Sys.Date(), "%Y%m%d"), ".Rdata"))
# load("data/essais_cutoff_toxrapp_20210830.Rdata")

heatmap_seuils_40 <- map(c("C_", "gamma", "alpha_calc", "puissance_calc"), heatmap_top, data = table_resume_tox) %>% 
  wrap_plots()
ggsave(plot = heatmap_seuils_40,
       filename = "figures/brouillon/heatmap_seuils_tox_40.png",
       device = "png", height = 8, width = 16)

tab_probas <- expand.grid(pas_0 = seq(0, .15, .015), pas_1 = seq(0, .2, .02)) %>% 
  left_join(corr_2 %>% `colnames<-`(paste0(colnames(.), "_H0")), by = c("pas_0" = "Pet_H0")) %>% 
  left_join(corr_3 %>% `colnames<-`(paste0(colnames(.), "_H1")), by = c("pas_1" = "Pet_H1"))
liste_resultats <- list()
for (i in seq_len(nrow(tab_probas))) {
  liste_resultats[[i]] <- deter_cutoff(alpha = .05,
                                       ana_inter = c(30, 51),
                                       ana_inter_tox = c(rep(10, 7), 11),
                                       p_n = unlist(tab_probas$probas_H0[i]),
                                       p_a = unlist(tab_probas$probas_H1[i]),
                                       mat_beta_xi = matrix(c(1, 1, 0, 0, 0, 1, 0, 1), nrow = 2, byrow = TRUE),
                                       methode = 3L,
                                       seed = 1993,
                                       affich_mat = "No")[[1]]
}
table_resume_tox_2 <- cbind(tab_probas, do.call(rbind, liste_resultats))
save(table_resume_tox_2, file = paste0("~/essais_cutoff_tox_bis_", format(Sys.Date(), "%Y%m%d"), ".Rdata"))
# load("data/essais_cutoff_tox_bis_20210830.Rdata")

heatmap_seuils_30 <- map(c("C_", "gamma", "alpha_calc", "puissance_calc"), heatmap_top, data = table_resume_tox_2) %>% 
  wrap_plots()
ggsave(plot = heatmap_seuils_30,
       filename = "figures/brouillon/heatmap_seuils_tox_30eng.png",
       device = "png", height = 8, width = 16)


# III/ Réviser les operating characteristics selon les cas de figures pour les hypothèses ----

source("R/realisations_essais_metimgast.R", encoding = "utf-8")

calculate_corr <- function(Peff, Ptox, Pinter) {
  (Pinter - Peff * Ptox) / sqrt((Peff - Peff ^ 2) * (Ptox - Ptox ^ 2))
}

calculate_pinter <- function(Peff, Ptox, corr) {
  corr * sqrt((Peff - Peff ^ 2) * (Ptox - Ptox ^ 2)) + Peff * Ptox
}

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

## A/ Ptox = 30% ----

liste_scenars <- list(
  H0 = c(Peff = .15, Ptox = .3),
  H1 = c(Peff = .3, Ptox = .2),
  Int1 = c(Peff = .2, Ptox = .25),
  Int2 = c(Peff = .25, Ptox = .25),
  Paseff1 = c(Peff = .15, Ptox = .2),
  Paseff2 = c(Peff = .1, Ptox = .15),
  Effint1 = c(Peff = .2, Ptox = .2),
  Effint2 = c(Peff = .2, Ptox = .3),
  Tox1 = c(Peff = .3, Ptox = .3),
  Tox2 = c(Peff = .4, Ptox = .35)
)
noms_scenars <- names(liste_scenars)
tableau_scenars <- map2_dfr(
  liste_scenars,
  noms_scenars,
  ~ generer_6corr_hyp(.x["Peff"], .x["Ptox"], .y)
)

load("data/cutoff_planifie.Rdata")
seuil <- seuil[[1]]

tableau_scenars <- tableau_scenars %>% 
  pull(nom_corr) %>% 
  unique() %>% 
  map_dfr(function(x) {
    seuil <- deter_cutoff(alpha = .05,
                          ana_inter = c(30, 51),
                          ana_inter_tox = c(30, 51),
                          p_n = tableau_scenars %>% filter(nom_corr == x, scenar == "H0") %>% pull(probas) %>% .[[1]],
                          p_a = tableau_scenars %>% filter(nom_corr == x, scenar == "H1") %>% pull(probas) %>% .[[1]],
                          affich_mat = "No",
                          methode = 3L)[[1]]
    c(nom_corr = x, setNames(seuil, paste0(names(seuil), "_classique")))
  }) %>% 
  right_join(tableau_scenars, by = "nom_corr")
tableau_scenars <- tableau_scenars %>% 
  pull(nom_corr) %>% 
  unique() %>% 
  map_dfr(function(x) {
    seuil <- deter_cutoff(alpha = .05,
                          ana_inter = c(30, 51),
                          ana_inter_tox = c(rep(5, 4), rep(10, 5), 11),
                          p_n = tableau_scenars %>% filter(nom_corr == x, scenar == "H0") %>% pull(probas) %>% .[[1]],
                          p_a = tableau_scenars %>% filter(nom_corr == x, scenar == "H1") %>% pull(probas) %>% .[[1]],
                          affich_mat = "No",
                          methode = 3L)[[1]]
    c(nom_corr = x, setNames(seuil, paste0(names(seuil), "_toxrapp")))
  }) %>% 
  right_join(tableau_scenars, by = "nom_corr")
save(tableau_scenars, file = paste0("data/scenar", format(Sys.Date(), "%Y%m%d"), ".Rdata"))

load("data/scenar20210907.Rdata")

noms_correlations <- unique(tableau_scenars$nom_corr)
tableau_scenars[, c("analyse_med_topiva", "nb_pts_topiva", "temps_etude_topiva", 
                    "accept_ttt_topiva", "arret_fut_topiva", "arret_tox_topiva", "nb_eff_topiva", 
                    "nb_tox_topiva")] <- NA_real_
tableau_scenars[, c("analyse_med_bopiva", "nb_pts_bopiva", "temps_etude_bopiva", 
                    "accept_ttt_bopiva", "arret_fut_bopiva", "arret_tox_bopiva", "nb_eff_bopiva", 
                    "nb_tox_bopiva")] <- NA_real_
tableau_scenars[, c("moy_pts_classique", "moy_duree_classique", "rejet_h0_classique", 
                    "arret_fut_classique", "arret_tox_classique", "nb_eff_classique", 
                    "nb_tox_classique")] <- NA_real_
tableau_scenars[, c("moy_pts_toxrapp", "moy_duree_toxrapp", "rejet_h0_toxrapp", 
                    "arret_fut_toxrapp", "arret_tox_toxrapp", "nb_eff_toxrapp", 
                    "nb_tox_toxrapp")] <- NA_real_

for (i in noms_scenars) {
  
  for (j in noms_correlations) {
    
    cat(paste0(i, ":", j, "\n"))
    
    p_reel <- tableau_scenars$probas[tableau_scenars$scenar == i & tableau_scenars$nom_corr == j] %>% pluck(1) %>% unname()
    
    essais <- gen_patients_multinomTOP(n_sim = 10000,
                                       ana_inter = c(30, 51),
                                       ana_inter_tox = c(rep(5, 4), rep(10, 5), 11),
                                       interpatient = 6, max_teff = 180, max_ttox = 180, rand_ratio = 1,
                                       multinom_ttt = list(p_reel),
                                       seed = 1993)
    # Generate the cohort 2 with cMET amplification aside from the main cohort (same code as in the package)
    cohorte_2 <- lapply(
      X = essais,
      FUN = function(x) {
        res <- data.frame(temps_recrutement = sort(runif(9, 0, max(x$temps_recrutement))))
        res$tox <- rbinom(9, 1, p_reel[1] + p_reel[3])
        res$temps_tox <- gen_temps_patient(n = 9, tmax = 180)
        res$temps_tox <- ifelse(res$tox == 1, res$temps_recrutement + res$temps_tox, 1e+12)
        res$temps_obstox <- res$temps_recrutement + max_ttox
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
                                                                 C_ = seuil[["C_"]], gamm = seuil[["gamma"]], critere_tox = .95,
                                                                 interpatient = 6, max_teff = 180, max_ttox = 180),
                                 .id = "essai")
    tableau_scenars[tableau_scenars$nom_corr == j & tableau_scenars$scenar == i, 
                    c("nb_pts_topiva", "temps_etude_topiva", 
                      "arret_fut_topiva", "arret_tox_topiva", "nb_eff_topiva", 
                      "nb_tox_topiva", "accept_ttt_topiva")] <- resultats_topiva %>% 
      summarise(across(c(nb_pts, temps_etude, arret_fut:nb_tox), ~ mean(.x)),
                accept_ttt = mean(decision == "Acceptation traitement")) %>% 
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
                                                                 p_n = tableau_scenars$probas[tableau_scenars$scenar == "H0" & tableau_scenars$nom_corr == j] %>% pluck(1) %>% unname(), 
                                                                 C_ = seuil[["C_"]], gamm = seuil[["gamma"]], critere_tox = .95,
                                                                 interpatient = 6, max_teff = 180, max_ttox = 180),
                                 .id = "essai")
    tableau_scenars[tableau_scenars$nom_corr == j & tableau_scenars$scenar == i, 
                    c("nb_pts_bopiva", "temps_etude_bopiva", 
                      "arret_fut_bopiva", "arret_tox_bopiva", "nb_eff_bopiva", 
                      "nb_tox_bopiva", "accept_ttt_bopiva")] <- resultats_bopiva %>% 
      summarise(across(c(nb_pts, temps_etude, arret_fut:nb_tox), ~ mean(.x)),
                accept_ttt = mean(decision == "Acceptation traitement")) %>% 
      `colnames<-`(paste0(colnames(.), "_bopiva")) %>% 
      tibble_row()
    
    # TOP classique avec 2 analyses intermédiaires
    resultats_topclas <- simu_essais(alpha = .05,
                                     ana_inter = c(30, 51),
                                     ana_inter_tox = c(30, 51),
                                     p_n = tableau_scenars$probas[tableau_scenars$scenar == "H0" & tableau_scenars$nom_corr == j] %>% pluck(1) %>% unname(),
                                     p_a = tableau_scenars$probas[tableau_scenars$scenar == "H1" & tableau_scenars$nom_corr == j] %>% pluck(1) %>% unname(),
                                     p_reel = tableau_scenars$probas[tableau_scenars$scenar == i & tableau_scenars$nom_corr == j] %>% pluck(1) %>% unname(),
                                     liste_tableaux = essais,
                                     cut_seq = tableau_scenars$C__classique[tableau_scenars$scenar == i & tableau_scenars$nom_corr == j] %>% as.numeric(),
                                     power_seq = tableau_scenars$gamma_classique[tableau_scenars$scenar == i & tableau_scenars$nom_corr == j] %>% as.numeric(),
                                     interpatient = 6, max_teff = 180, max_ttox = 180,
                                     stock_matrice = FALSE)
    tableau_scenars[tableau_scenars$nom_corr == j & tableau_scenars$scenar == i, 
                    c("moy_pts_classique", "moy_duree_classique", "rejet_h0_classique", 
                        "arret_fut_classique", "arret_tox_classique", "nb_eff_classique", 
                        "nb_tox_classique")] <- resultats_topclas %>% 
      `names<-`(paste0(names(resultats_topclas), "_classique")) %>% 
      as_tibble_row()
    
    # TOP avec monitoring rapproché de la toxicité
    resultats_toprapp <- simu_essais(alpha = .05,
                                     ana_inter = c(30, 51),
                                     ana_inter_tox = c(rep(5, 4), rep(10, 5), 11),
                                     p_n = tableau_scenars$probas[tableau_scenars$scenar == "H0" & tableau_scenars$nom_corr == j] %>% pluck(1) %>% unname(),
                                     p_a = tableau_scenars$probas[tableau_scenars$scenar == "H1" & tableau_scenars$nom_corr == j] %>% pluck(1) %>% unname(),
                                     p_reel = tableau_scenars$probas[tableau_scenars$scenar == i & tableau_scenars$nom_corr == j] %>% pluck(1) %>% unname(),
                                     liste_tableaux = essais,
                                     cut_seq = tableau_scenars$C__toxrapp[tableau_scenars$scenar == i & tableau_scenars$nom_corr == j] %>% as.numeric(),
                                     power_seq = tableau_scenars$gamma_toxrapp[tableau_scenars$scenar == i & tableau_scenars$nom_corr == j] %>% as.numeric(),
                                     interpatient = 6, max_teff = 180, max_ttox = 180,
                                     stock_matrice = FALSE)
    tableau_scenars[tableau_scenars$nom_corr == j & tableau_scenars$scenar == i, 
                    c("moy_pts_toxrapp", "moy_duree_toxrapp", "rejet_h0_toxrapp", 
                      "arret_fut_toxrapp", "arret_tox_toxrapp", "nb_eff_toxrapp", 
                      "nb_tox_toxrapp")] <- resultats_toprapp %>% 
      `names<-`(paste0(names(resultats_toprapp), "_toxrapp")) %>% 
      as_tibble_row()
    
  }
  
}

save(tableau_scenars, file = paste0("data/resultats_tox30_", format(Sys.Date(), "%Y%m%d"), ".Rdata"))

## B/ Ptox = 40% ----

liste_scenars_bis <- list(
  H0 = c(Peff = .15, Ptox = .4),
  H1 = c(Peff = .3, Ptox = .2),
  Int1 = c(Peff = .2, Ptox = .3),
  Int2 = c(Peff = .25, Ptox = .35),
  Paseff1 = c(Peff = .15, Ptox = .2),
  Paseff2 = c(Peff = .1, Ptox = .15),
  Effint1 = c(Peff = .2, Ptox = .2),
  Effint2 = c(Peff = .2, Ptox = .4),
  Tox1 = c(Peff = .3, Ptox = .4),
  Tox2 = c(Peff = .4, Ptox = .45)
)
noms_scenars_bis <- names(liste_scenars_bis)
tableau_scenars_bis <- map2_dfr(
  liste_scenars_bis,
  noms_scenars_bis,
  ~ generer_6corr_hyp(.x["Peff"], .x["Ptox"], .y)
)

load("data/cutoff_planifie.Rdata")
seuil <- seuil[[1]]

tableau_scenars_bis <- tableau_scenars_bis %>% 
  pull(nom_corr) %>% 
  unique() %>% 
  map_dfr(function(x) {
    seuil <- deter_cutoff(alpha = .05,
                          ana_inter = c(30, 51),
                          ana_inter_tox = c(30, 51),
                          p_n = tableau_scenars_bis %>% filter(nom_corr == x, scenar == "H0") %>% pull(probas) %>% .[[1]],
                          p_a = tableau_scenars_bis %>% filter(nom_corr == x, scenar == "H1") %>% pull(probas) %>% .[[1]],
                          affich_mat = "No",
                          methode = 3L)[[1]]
    c(nom_corr = x, setNames(seuil, paste0(names(seuil), "_classique")))
  }) %>% 
  right_join(tableau_scenars_bis, by = "nom_corr")
tableau_scenars_bis <- tableau_scenars_bis %>% 
  pull(nom_corr) %>% 
  unique() %>% 
  map_dfr(function(x) {
    seuil <- deter_cutoff(alpha = .05,
                          ana_inter = c(30, 51),
                          ana_inter_tox = c(rep(5, 4), rep(10, 5), 11),
                          p_n = tableau_scenars_bis %>% filter(nom_corr == x, scenar == "H0") %>% pull(probas) %>% .[[1]],
                          p_a = tableau_scenars_bis %>% filter(nom_corr == x, scenar == "H1") %>% pull(probas) %>% .[[1]],
                          affich_mat = "No",
                          methode = 3L)[[1]]
    c(nom_corr = x, setNames(seuil, paste0(names(seuil), "_toxrapp")))
  }) %>% 
  right_join(tableau_scenars_bis, by = "nom_corr")
save(tableau_scenars_bis, file = paste0("data/scenar_bis", format(Sys.Date(), "%Y%m%d"), ".Rdata"))

load("data/scenar_bis20210907.Rdata")

noms_correlations_bis <- unique(tableau_scenars_bis$nom_corr)
tableau_scenars_bis[, c("nb_pts_topiva", "temps_etude_topiva", 
                        "accept_ttt_topiva", "arret_fut_topiva", "arret_tox_topiva", "nb_eff_topiva", 
                        "nb_tox_topiva")] <- NA_real_
tableau_scenars_bis[, c("nb_pts_bopiva", "temps_etude_bopiva", 
                        "accept_ttt_bopiva", "arret_fut_bopiva", "arret_tox_bopiva", "nb_eff_bopiva", 
                        "nb_tox_bopiva")] <- NA_real_
tableau_scenars_bis[, c("moy_pts_classique", "moy_duree_classique", "rejet_h0_classique", 
                        "arret_fut_classique", "arret_tox_classique", "nb_eff_classique", 
                        "nb_tox_classique")] <- NA_real_
tableau_scenars_bis[, c("moy_pts_toxrapp", "moy_duree_toxrapp", "rejet_h0_toxrapp", 
                        "arret_fut_toxrapp", "arret_tox_toxrapp", "nb_eff_toxrapp", 
                        "nb_tox_toxrapp")] <- NA_real_

for (i in noms_scenars_bis) {
  
  for (j in noms_correlations_bis) {
    
    cat(paste0(i, ":", j, "\n"))
    
    essais <- gen_patients_multinomTOP(n_sim = 10000,
                                       ana_inter = c(30, 51),
                                       ana_inter_tox = c(rep(5, 4), rep(10, 5), 11),
                                       interpatient = 6, max_teff = 180, max_ttox = 180, rand_ratio = 1,
                                       multinom_ttt = list(tableau_scenars_bis$probas[tableau_scenars_bis$scenar == i & tableau_scenars_bis$nom_corr == j] %>% pluck(1) %>% unname()),
                                       seed = 1993)
    
    # TOP avec approche Ivanova au temps du TOP
    resultats_topiva <- map_dfr(.x = essais,
                                .f = ~ realisation_essai_topiva(ana_inter_eff = c(30, 51),
                                                                ana_inter_tox = c(rep(5, 4), rep(10, 5), 11),
                                                                tableau = .x,
                                                                phitox = c(.15, .25), prior_eff = .15, prior_tox = c(1, 1), 
                                                                C_ = seuil[["C_"]], gamm = seuil[["gamma"]], critere_tox = .95,
                                                                interpatient = 6, max_teff = 180, max_ttox = 180),
                                .id = "essai")
    tableau_scenars_bis[tableau_scenars_bis$nom_corr == j & tableau_scenars_bis$scenar == i, 
                    c("nb_pts_topiva", "temps_etude_topiva", 
                      "arret_fut_topiva", "arret_tox_topiva", "nb_eff_topiva", 
                      "nb_tox_topiva", "accept_ttt_topiva")] <- resultats_topiva %>% 
      summarise(across(c(nb_pts, temps_etude, arret_fut:nb_tox), ~ mean(.x)),
                accept_ttt = mean(decision == "Acceptation traitement")) %>% 
      `colnames<-`(paste0(colnames(.), "_topiva")) %>% 
      tibble_row()
    
    # BOP avec approche Ivanova donc observation complète à chaque analyse intermédiaire
    resultats_bopiva <- map_dfr(.x = essais,
                                .f = ~ realisation_essai_bopiva(ana_inter_eff = c(30, 51),
                                                                ana_inter_tox = c(rep(5, 4), rep(10, 5), 11),
                                                                tableau = .x,
                                                                phitox = c(.15, .25), prior_eff = .15, prior_tox = c(1, 1), 
                                                                p_n = tableau_scenars_bis$probas[tableau_scenars_bis$scenar == "H0" & tableau_scenars_bis$nom_corr == j] %>% pluck(1) %>% unname(), 
                                                                C_ = seuil[["C_"]], gamm = seuil[["gamma"]], critere_tox = .95,
                                                                interpatient = 6, max_teff = 180, max_ttox = 180),
                                .id = "essai")
    tableau_scenars_bis[tableau_scenars_bis$nom_corr == j & tableau_scenars_bis$scenar == i, 
                    c("nb_pts_bopiva", "temps_etude_bopiva", 
                      "arret_fut_bopiva", "arret_tox_bopiva", "nb_eff_bopiva", 
                      "nb_tox_bopiva", "accept_ttt_bopiva")] <- resultats_bopiva %>% 
      summarise(across(c(nb_pts, temps_etude, arret_fut:nb_tox), ~ mean(.x)),
                accept_ttt = mean(decision == "Acceptation traitement")) %>% 
      `colnames<-`(paste0(colnames(.), "_bopiva")) %>% 
      tibble_row()
    
    # TOP classique avec 2 analyses intermédiaires
    resultats_topclas <- simu_essais(alpha = .05,
                                     ana_inter = c(30, 51),
                                     ana_inter_tox = c(30, 51),
                                     p_n = tableau_scenars_bis$probas[tableau_scenars_bis$scenar == "H0" & tableau_scenars_bis$nom_corr == j] %>% pluck(1) %>% unname(),
                                     p_a = tableau_scenars_bis$probas[tableau_scenars_bis$scenar == "H1" & tableau_scenars_bis$nom_corr == j] %>% pluck(1) %>% unname(),
                                     p_reel = tableau_scenars_bis$probas[tableau_scenars_bis$scenar == i & tableau_scenars_bis$nom_corr == j] %>% pluck(1) %>% unname(),
                                     liste_tableaux = essais,
                                     cut_seq = tableau_scenars_bis$C__classique[tableau_scenars_bis$scenar == i & tableau_scenars_bis$nom_corr == j] %>% as.numeric(),
                                     power_seq = tableau_scenars_bis$gamma_classique[tableau_scenars_bis$scenar == i & tableau_scenars_bis$nom_corr == j] %>% as.numeric(),
                                     interpatient = 6, max_teff = 180, max_ttox = 180,
                                     stock_matrice = FALSE)
    tableau_scenars_bis[tableau_scenars_bis$nom_corr == j & tableau_scenars_bis$scenar == i, 
                    c("moy_pts_classique", "moy_duree_classique", "rejet_h0_classique", 
                      "arret_fut_classique", "arret_tox_classique", "nb_eff_classique", 
                      "nb_tox_classique")] <- resultats_topclas %>% 
      `names<-`(paste0(names(resultats_topclas), "_classique")) %>% 
      as_tibble_row()
    
    # TOP avec monitoring rapproché de la toxicité
    resultats_toprapp <- simu_essais(alpha = .05,
                                     ana_inter = c(30, 51),
                                     ana_inter_tox = c(rep(5, 4), rep(10, 5), 11),
                                     p_n = tableau_scenars_bis$probas[tableau_scenars_bis$scenar == "H0" & tableau_scenars_bis$nom_corr == j] %>% pluck(1) %>% unname(),
                                     p_a = tableau_scenars_bis$probas[tableau_scenars_bis$scenar == "H1" & tableau_scenars_bis$nom_corr == j] %>% pluck(1) %>% unname(),
                                     p_reel = tableau_scenars_bis$probas[tableau_scenars_bis$scenar == i & tableau_scenars_bis$nom_corr == j] %>% pluck(1) %>% unname(),
                                     liste_tableaux = essais,
                                     cut_seq = tableau_scenars_bis$C__toxrapp[tableau_scenars_bis$scenar == i & tableau_scenars_bis$nom_corr == j] %>% as.numeric(),
                                     power_seq = tableau_scenars_bis$gamma_toxrapp[tableau_scenars_bis$scenar == i & tableau_scenars_bis$nom_corr == j] %>% as.numeric(),
                                     interpatient = 6, max_teff = 180, max_ttox = 180,
                                     stock_matrice = FALSE)
    tableau_scenars_bis[tableau_scenars_bis$nom_corr == j & tableau_scenars_bis$scenar == i, 
                    c("moy_pts_toxrapp", "moy_duree_toxrapp", "rejet_h0_toxrapp", 
                      "arret_fut_toxrapp", "arret_tox_toxrapp", "nb_eff_toxrapp", 
                      "nb_tox_toxrapp")] <- resultats_toprapp %>% 
      `names<-`(paste0(names(resultats_toprapp), "_toxrapp")) %>% 
      as_tibble_row()
    
  }
  
}

save(tableau_scenars_bis, file = paste0("data/resultats_tox40_", format(Sys.Date(), "%Y%m%d"), ".Rdata"))

