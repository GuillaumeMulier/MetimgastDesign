# ------------------------------------------------------------- #
# Script des résultats pour la thèse sur l'essai METIMGAST      #
# Après la réunion avec Ruitao Lin                              #
# Créé le 05/01/2022, modifié le 05/01/2022                     #
# ------------------------------------------------------------- #


# Packages et helpers ----

# devtools::install_github("GuillaumeMulier/multibrasBOP2", auth_token = "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")
# La clef sera à actualiser tous les 30 jours sur github

library(tidyverse)
library(ggtext)
library(multibrasBOP2)
library(patchwork)
library(clinfun)
library(cowplot)

theme_set(theme_light() +
            theme(strip.background = element_rect(fill = NA),
                  strip.text = element_textbox(
                    size = 12, 
                    color = "white", fill = "#7888C0", box.color = "#000066",
                    halign = 0.5, linetype = 1, r = unit(3, "pt"), width = unit(0.75, "npc"),
                    padding = margin(2, 0, 1, 0), margin = margin(3, 3, 3, 3))))

vec_eval_tox <- c(5, 10, 15, 20, 30, 40, 50, 60, 70, 80, 90)


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


# II/ Déterminer les règles d'arrêt optimales pour le design d'Ivanova ----

source("R/realisations_essais_metimgast_v2.R", encoding = "UTF8")

## A/ En ne prenant en compte que la toxicité ----

set.seed(1993)

# On génère des données de toxicité sous H0 et H1
donnees_h0 <- map(1:10000, ~ rbinom(90, size = 1, prob = .3))
donnees_h0 <- map(donnees_h0,
                  function(data) {
                    map_dbl(vec_eval_tox, ~ sum(data[seq_len(.x)])) %>% 
                      bind_cols(n = vec_eval_tox, ntox = .)
                  })
donnees_h1 <- map(1:10000, ~ rbinom(90, size = 1, prob = .2))
donnees_h1 <- map(donnees_h1,
                  function(data) {
                    map_dbl(vec_eval_tox, ~ sum(data[seq_len(.x)])) %>% 
                      bind_cols(n = vec_eval_tox, ntox = .)
                  })

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

table_choix <- expand.grid(
  pi = seq(.2, .3, .01),
  tau = seq(.86, .99, .005)
) 
table_choix$alpha <- NA_real_
table_choix$puissance <- NA_real_

for (i in seq_len(nrow(table_choix))) {
  
  cat(i, "/", nrow(table_choix), ".\n")
  
  table_regles <- regle_arret(c(1, 1), vec_eval_tox, table_choix$pi[i], table_choix$tau[i])
  
  table_choix$alpha[i] <- decision_essai(donnees_h0, table_regles) %>% 
    pull(decision) %>% 
    str_detect("Accept treatment") %>% 
    sum() %>% 
    '/'(10000)
  
  table_choix$puissance[i] <- decision_essai(donnees_h1, table_regles) %>% 
    pull(decision) %>% 
    str_detect("Accept treatment") %>% 
    sum() %>% 
    '/'(10000)
  
}

plot1 <- ggplot(table_choix, aes(x = pi, y = tau, fill = alpha, label = paste0(round(100 * alpha, 1), "%"))) +
  geom_tile() +
  geom_text(color = "white") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = expression(pi[tox]),
       y = expression(tau),
       fill = "Type I error rate")
plot2 <- ggplot(table_choix, aes(x = pi, y = tau, fill = puissance, label = paste0(round(100 * puissance, 1), "%"))) +
  geom_tile() +
  geom_text(color = "white") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = expression(pi[tox]),
       y = expression(tau),
       fill = "Power")
plot_tot <- (plot1 | plot2) +
  plot_annotation(title = expression(paste("Power and type I error rate for ", PP[tox], " design on toxicity")),
                  subtitle = expression(paste("Stopping rule : P(", p[tox], ">", pi[tox], "|", D[n], ")>", tau))) &
  theme(text = element_text(size = 20))

ggsave(plot = plot_tot,
       filename = "figures/rech_seuil/recherche_seuil_iva_tox.png",
       device = "png", height = 9, width = 16)
save(table_choix, file = "data/alpha_puiss_iva_tox.Rdata")

rm(donnees_h0, donnees_h1, plot1, plot2, plot_tot, table_choix)


## B/ En prenant en compte la toxicité et l'efficacité ----

donnees_h0 <- gen_patients_multinomTOP(n_sim = 10000,
                                       ana_inter = c(30, 51),
                                       ana_inter_tox = c(rep(5, 4), rep(10, 5), 11),
                                       interpatient = 6, max_teff = 180, max_ttox = 42, rand_ratio = 1,
                                       multinom_ttt = list(c(.08, .07, .22, .63)),
                                       seed = 1993)
cohorte_h0 <- lapply(
  X = donnees_h0,
  FUN = function(x) {
    res <- data.frame(temps_recrutement = sort(runif(9, 0, max(x$temps_recrutement))))
    res$tox <- rbinom(9, 1, .3)
    res$temps_tox <- gen_temps_patient(n = 9, tmax = 42)
    res$temps_tox <- ifelse(res$tox == 1, res$temps_recrutement + res$temps_tox, 1e+12)
    res$temps_obstox <- res$temps_recrutement + 42
    return(res)
  }
)

donnees_h1 <- gen_patients_multinomTOP(n_sim = 10000,
                                       ana_inter = c(30, 51),
                                       ana_inter_tox = c(rep(5, 4), rep(10, 5), 11),
                                       interpatient = 6, max_teff = 180, max_ttox = 42, rand_ratio = 1,
                                       multinom_ttt = list(c(.107, .193, .093, .607)),
                                       seed = 1993)
cohorte_h1 <- lapply(
  X = donnees_h1,
  FUN = function(x) {
    res <- data.frame(temps_recrutement = sort(runif(9, 0, max(x$temps_recrutement))))
    res$tox <- rbinom(9, 1, .2)
    res$temps_tox <- gen_temps_patient(n = 9, tmax = 42)
    res$temps_tox <- ifelse(res$tox == 1, res$temps_recrutement + res$temps_tox, 1e+12)
    res$temps_obstox <- res$temps_recrutement + 42
    return(res)
  }
)

tableau <- donnees_h0[[1]]
cohorte <- cohorte_h0[[1]]
tab_tox <- regle_arret(c(1, 1), vec_eval_tox, .25, .95)
ana_inter_eff <- c(30, 51)
ana_inter_tox <- c(rep(5, 4), rep(10, 6), 1)
phitox = c(.15, .25)
prior_eff = .15
# prior_tox = c(1, 1) 
C_ = seuil_effseule[["C_"]]
gamm = seuil_effseule[["gamma"]]
# critere_tox = .95
interpatient = 6
max_teff = 180
max_ttox = 42

simu_topiva(ana_inter_eff = c(30, 51),
            ana_inter_tox = c(rep(5, 4), rep(10, 6), 1),
            liste_essais = donnees_h0, 
            liste_cohortes = cohorte_h0,
            phitox = c(.15, .25), prior_eff = .15, 
            prior_tox = c(1, 1), critere_tox = .95,
            C_ = seuil_effseule[["C_"]], gamm = seuil_effseule[["gamma"]],
            interpatient = 6, max_teff = 180, max_ttox = 42) 



# III/ Design de Simon ----

## A/ Simulations ----

design_simon <- ph2simon(pu = .15, pa = .30, ep1 = .05, ep2 = .1)
regles_simon <- design_simon$out[18, ]

tab_eff <- tribble(
  ~nb_ana, ~tot_pat, ~seuil_eff,
  1, 30, 4,
  2, 81, 17
)

n1 <- 30
n <- 81
r1 <- 5
r <- 17
pu <- .15
pa <- .3

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

donnees_test <- expand.grid(
  n1 = 30,
  r1 = 0:30,
  n = 81,
  r = 0:81,
  pu = .15,
  pa = .3
) %>% 
  filter(r >= r1)
result_simu_simon <- map_dfr(seq_len(nrow(donnees_test)), ~ alpha_puiss_simon(donnees_test$r1[.x], donnees_test$n1[.x], donnees_test$r[.x], donnees_test$n[.x], donnees_test$pu[.x], donnees_test$pa[.x]))
print(result_simu_simon, n = Inf)
result_simu_simon %>% filter(alpha <= .05, puissance >= .9) %>% arrange(EN_p0)
result_simu_simon %>% filter(alpha <= .05, puissance >= .9) %>% arrange(desc(puissance))
result_simu_simon %>% filter(alpha <= .05, puissance >= .8) %>% arrange(EN_p0)
result_simu_simon %>% filter(alpha <= .05, puissance >= .8) %>% arrange(desc(puissance))

tab_tox <- regle_arret(c(1, 1), vec_eval_tox, .20, .95)

result_simon <- liste_scenars %>% 
  select(scenar, R, probas)
result_simon[, c("nb_pts", "temps_etude", "arret_fut", "arret_tox", "nb_eff", "nb_tox", "accept_ttt")] <- NA_real_
result_simon$analyse_med <- NA_character_

for (i in seq_len(nrow(result_simon))) {
  
  cat(i, "/", nrow(result_simon), ".\n")
  p_reel <- result_simon$probas[i] %>% pluck(1) %>% unname()
  tab_donnee <- gen_patients_multinomTOP(n_sim = 10000,
                                         ana_inter = c(30, 51),
                                         ana_inter_tox = c(rep(5, 4), rep(10, 5), 11),
                                         interpatient = 6, max_teff = 180, max_ttox = 42, rand_ratio = 1,
                                         multinom_ttt = list(p_reel),
                                         seed = 1993)
  set.seed(1993)
  cohorte <- lapply(
    X = tab_donnee,
    FUN = function(x) {
      res <- data.frame(temps_recrutement = sort(runif(9, 0, max(x$temps_recrutement))))
      res$tox <- rbinom(9, 1, sum(p_reel * c(1, 0, 1, 0)))
      res$temps_tox <- gen_temps_patient(n = 9, tmax = 42)
      res$temps_tox <- ifelse(res$tox == 1, res$temps_recrutement + res$temps_tox, 1e+12)
      res$temps_obstox <- res$temps_recrutement + 42
      return(res)
    }
  )
  result_simon[i, c("nb_pts", "temps_etude", "arret_fut", "arret_tox", "nb_eff", "nb_tox", "accept_ttt", "rret_prema", "analyse_med")] <- 
    map2_dfr(.x = tab_donnee,
             .y = cohorte,
             .f = ~ realisation_essai_bopiva(ana_inter_eff = c(30, 51),
                                             ana_inter_tox = c(rep(5, 4), rep(10, 6), 1),
                                             tableau = .x, tab_tox = tab_tox, tab_eff = tab_eff,
                                             cohorte = .y),
             .id = "essai") %>% 
    summarise(across(c(nb_pts, temps_etude, arret_fut:nb_tox), ~ mean(.x)),
              accept_ttt = mean(decision == "Acceptation traitement"),
              arret_prema = mean(str_detect(decision, "prématuré")),
              analyse_med = paste0("med = ", median(analyse), " / ", 
                                   paste(names(table(.$analyse)), "=", round(100 * prop.table(table(.$analyse)), 1), "%", collapse = " - ")))
  
}

save(result_simon, file = "data/result_simon_v3.Rdata")

rm(result_simon, tab_donnee, cohorte)

## B / Figures ----

load("data/resultats_tox30_rpos_20211005.Rdata")
load("data/result_simon_v3.Rdata")
tableau_results <- select(tableau_results, -c(analyse_med_bopiva, analyse_med_topiva)) %>% 
  mutate(arret_both_topiva = arret_fut_topiva + arret_tox_topiva - (1 - accept_ttt_topiva),
         arret_both_bopiva = arret_fut_bopiva + arret_tox_bopiva - (1 - accept_ttt_bopiva),
         arret_both_classique = arret_fut_classique + arret_tox_classique - (1 - rejet_h0_classique),
         arret_both_toxrapp = arret_fut_toxrapp + arret_tox_toxrapp - (1 - rejet_h0_toxrapp),
         arret_both_toxprior = arret_fut_toxprior + arret_tox_toxprior - (1 - rejet_h0_toxprior))
tableau_results_simon <- select(result_simon, -analyse_med) %>% 
  `names<-`(c(names(.)[1:3], paste0(names(.)[-c(1:3)], "_simon"))) %>% 
  mutate(arret_both_simon = arret_fut_simon + arret_tox_simon - (1 - accept_ttt_simon))
names(tableau_results) <- str_replace_all(names(tableau_results), "nb_pts", "moy_pts")
names(tableau_results) <- str_replace_all(names(tableau_results), "temps_etude", "moy_duree")
names(tableau_results) <- str_replace_all(names(tableau_results), "rejet_h0", "accept_ttt")
names(tableau_results_simon) <- str_replace_all(names(tableau_results_simon), "nb_pts", "moy_pts")
names(tableau_results_simon) <- str_replace_all(names(tableau_results_simon), "temps_etude", "moy_duree")
names(tableau_results_simon) <- str_replace_all(names(tableau_results_simon), "rejet_h0", "accept_ttt")
tableau_results$scenar_p <- factor(tableau_results$scenar,
                                   levels = rev(c("H0", "H1", "Int1", "Int2", "Paseff1", "Paseff2", "Effint1", "Effint2", "Tox1", "Tox2")),
                                   labels = rev(c("Sc 1 : (0.15;0.30) [H<sub>0</sub>]", "Sc 2 : (0.30;0.20) [H<sub>1</sub>]", "Sc 3: (0.20;0.25)", "Sc 4 : (0.25;0.25)", "Sc 5 : (0.15;0.20)",
                                                  "Sc 6 : (0.10;0.15)", "Sc 7 : (0.20;0.20)", "Sc 8 : (0.20;0.30)", "Sc 9 : (0.30;0.30)", "Sc 10 : (0.40;0.35)")))
tableau_results$scenarios <- factor(tableau_results$scenar,
                                    levels = rev(c("H0", "H1", "Int1", "Int2", "Paseff1", "Paseff2", "Effint1", "Effint2", "Tox1", "Tox2")),
                                    labels = rev(c("Sc 1 :", "Sc 2 :", "Sc 3 :", "Sc 4 :", "Sc 5 :",
                                                   "Sc 6 :", "Sc 7 :", "Sc 8 :", "Sc 9 :", "Sc 10 :")))
tableau_results$probas_plot <- factor(tableau_results$scenar,
                                      levels = rev(c("H0", "H1", "Int1", "Int2", "Paseff1", "Paseff2", "Effint1", "Effint2", "Tox1", "Tox2")),
                                      labels = rev(c("(0.15;0.30) [H<sub>0</sub>]", "(0.30;0.20) [H<sub>1</sub>]", "(0.20;0.25)", "(0.25;0.25)", "(0.15;0.20)",
                                                     "(0.10;0.15)", "(0.20;0.20)", "(0.20;0.30)", "(0.30;0.30)", "(0.40;0.35)")))
tableau_results_simon$scenar_p <- factor(tableau_results_simon$scenar,
                                         levels = rev(c("H0", "H1", "Int1", "Int2", "Paseff1", "Paseff2", "Effint1", "Effint2", "Tox1", "Tox2")),
                                         labels = rev(c("Sc 1 : (0.15;0.30) [H<sub>0</sub>]", "Sc 2 : (0.30;0.20) [H<sub>1</sub>]", "Sc 3: (0.20;0.25)", "Sc 4 : (0.25;0.25)", "Sc 5 : (0.15;0.20)",
                                                        "Sc 6 : (0.10;0.15)", "Sc 7 : (0.20;0.20)", "Sc 8 : (0.20;0.30)", "Sc 9 : (0.30;0.30)", "Sc 10 : (0.40;0.35)")))
tableau_results_simon$scenarios <- factor(tableau_results_simon$scenar,
                                          levels = rev(c("H0", "H1", "Int1", "Int2", "Paseff1", "Paseff2", "Effint1", "Effint2", "Tox1", "Tox2")),
                                          labels = rev(c("Sc 1 :", "Sc 2 :", "Sc 3 :", "Sc 4 :", "Sc 5 :",
                                                         "Sc 6 :", "Sc 7 :", "Sc 8 :", "Sc 9 :", "Sc 10 :")))
tableau_results_simon$probas_plot <- factor(tableau_results_simon$scenar,
                                            levels = rev(c("H0", "H1", "Int1", "Int2", "Paseff1", "Paseff2", "Effint1", "Effint2", "Tox1", "Tox2")),
                                            labels = rev(c("(0.15;0.30) [H<sub>0</sub>]", "(0.30;0.20) [H<sub>1</sub>]", "(0.20;0.25)", "(0.25;0.25)", "(0.15;0.20)",
                                                           "(0.10;0.15)", "(0.20;0.20)", "(0.20;0.30)", "(0.30;0.30)", "(0.40;0.35)")))
tableau_results <- left_join(tableau_results, tableau_results_simon)

metadata <- tribble(
  ~var, ~label, ~lim_sup, ~breaks, ~format_legende,
  "accept_ttt", "Proportion d'acceptation du traitement", 1, expression(seq(0, 1, .2)), expression(scales::percent_format()),
  "arret_fut", "Proportion d'arrêt pour futilité", 1, expression(seq(0, 1, .2)), expression(scales::percent_format()),
  "arret_tox", "Proportion d'arrêt pour toxicité", 1, expression(seq(0, 1, .2)), expression(scales::percent_format()),
  "moy_pts", "Nombre moyen de patients par essai", 85, expression(seq(0, 80, 20)), expression(function(x) paste0(x, c(" pt", rep(" pts", length(x) - 1)))),
  "moy_duree", "Durée moyenne de l'étude", 1100, expression(seq(0, 1000, 200)), expression(function(x) paste0(round(x / 365.25, 2), c(" an", " an", rep(" ans", length(x) - 2)))),
  "nb_eff", "Nombre de réponses en moyenne par essai", 30, expression(seq(0, 30, 5)), expression(function(x) x),
  "nb_tox", "Nombre de toxicités en moyenne par essai", 30, expression(seq(0, 30, 5)), expression(function(x) x),
  "arret_prema", "Proportion d'arrêts précoces", 1, expression(seq(0, 1, .2)), expression(scales::percent_format())
)

metadata_eng <- tribble(
  ~var, ~label, ~lim_sup, ~breaks, ~format_legende,
  "accept_ttt", "% of acceptation of treatment", 1, expression(seq(0, 1, .2)), expression(scales::percent_format()),
  "arret_fut", "% futility stopping", 1, expression(seq(0, 1, .2)), expression(scales::percent_format()),
  "arret_tox", "% toxicity stopping", 1, expression(seq(0, 1, .2)), expression(scales::percent_format()),
  "moy_pts", "Mean number of patients by trial", 85, expression(seq(0, 80, 20)), expression(function(x) paste0(x, c(" pt", rep(" pts", length(x) - 1)))),
  "moy_duree", "Mean duration of a trial", 1100, expression(seq(0, 1000, 200)), expression(function(x) paste0(round(x / 365.25, 2), c(" year", " year", rep(" years", length(x) - 2)))),
  "nb_eff", "Mean number of responses per trial", 30, expression(seq(0, 30, 5)), expression(function(x) x),
  "nb_tox", "Mean number of toxicities per trial", 30, expression(seq(0, 30, 5)), expression(function(x) x),
  "arret_prema", "% early stopping", 1, expression(seq(0, 1, .2)), expression(scales::percent_format())
)

plot_opchar <- function(tableau_scenars_plot, metadata, var) {
  
  # Mise au format long
  tableau_res <- tableau_scenars_plot %>% 
    select(!matches("C_")) %>% 
    select(!matches("gamma")) %>% 
    select(-c(probas, EffTox, EffNotox, NoeffTox, NoeffNotox)) %>% 
    select(scenar, scenar_p, scenarios, probas_plot, R, everything()) %>% 
    pivot_longer(cols = -c(scenar:R),
                 names_pattern = "(.+_.+)_(.+)",
                 names_to = c("stat", "schemas")) %>% 
    mutate(schemas = factor(schemas, 
                            levels = c("topiva", "bopiva", "classique", "toxrapp", "toxprior", "simon"),
                            labels = c("TOP<sub>eff</sub> + Ivanova", "BOP<sub>eff</sub> + Ivanova", "TOP<sub>eff/tox</sub> with 2 analyses", "TOP<sub>eff/tox</sub> with close monitoring of toxicity", "TOP<sub>eff/tox</sub> with update of the prior of toxicity", "Simon + Ivanova")))
  
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


plot_accept <- plot_opchar(tableau_results, metadata_eng, "accept_ttt")
ggsave(filename = "figures/rech_seuil/calib_simon.png",
       plot = plot_accept, device = "png",
       height = 14, width = 20)

plot_pts <- plot_opchar(tableau_results, metadata_eng, "moy_pts")
ggsave(filename = "figures/rech_seuil/calib_simon_pts.png",
       plot = plot_pts, device = "png",
       height = 14, width = 20)

# IV/ Ivanova + TOP avec quelques valeurs de seuil ----

## A/ Simulations ----

result_iva <- liste_scenars %>% 
  select(scenar, R, probas)
result_iva[, c("nb_pts_2095", "temps_etude_2095", "arret_fut_2095", "arret_tox_2095", "nb_eff_2095", "nb_tox_2095", "accept_ttt_2095")] <- NA_real_
result_iva$analyse_med_2095 <- NA_character_
result_iva[, c("nb_pts_22595", "temps_etude_22595", "arret_fut_22595", "arret_tox_22595", "nb_eff_22595", "nb_tox_22595", "accept_ttt_22595")] <- NA_real_
result_iva$analyse_med_22595 <- NA_character_
result_iva[, c("nb_pts_2090", "temps_etude_2090", "arret_fut_2090", "arret_tox_2090", "nb_eff_2090", "nb_tox_2090", "accept_ttt_2090")] <- NA_real_
result_iva$analyse_med_2090 <- NA_character_
result_iva[, c("nb_pts_22590", "temps_etude_22590", "arret_fut_22590", "arret_tox_22590", "nb_eff_22590", "nb_tox_22590", "accept_ttt_22590")] <- NA_real_
result_iva$analyse_med_22590 <- NA_character_
result_iva[, c("nb_pts_2590", "temps_etude_2590", "arret_fut_2590", "arret_tox_2590", "nb_eff_2590", "nb_tox_2590", "accept_ttt_2590")] <- NA_real_
result_iva$analyse_med_2590 <- NA_character_

for (i in seq_len(nrow(result_iva))) {
  
  cat(i, "/", nrow(result_iva), ".\n")
  p_reel <- result_iva$probas[i] %>% pluck(1) %>% unname()
  tab_donnee <- gen_patients_multinomTOP(n_sim = 10000,
                                         ana_inter = c(31, 50),
                                         ana_inter_tox = c(rep(5, 4), rep(10, 5), 11),
                                         interpatient = 6, max_teff = 180, max_ttox = 42, rand_ratio = 1,
                                         multinom_ttt = list(p_reel),
                                         seed = 1993)
  set.seed(1993)
  cohorte <- lapply(
    X = tab_donnee,
    FUN = function(x) {
      res <- data.frame(temps_recrutement = sort(runif(9, 0, max(x$temps_recrutement))))
      res$tox <- rbinom(9, 1, sum(p_reel * c(1, 0, 1, 0)))
      res$temps_tox <- gen_temps_patient(n = 9, tmax = 42)
      res$temps_tox <- ifelse(res$tox == 1, res$temps_recrutement + res$temps_tox, 1e+12)
      res$temps_obstox <- res$temps_recrutement + 42
      return(res)
    }
  )
  result_iva[i, c("nb_pts_2095", "temps_etude_2095", "arret_fut_2095", "arret_tox_2095", "nb_eff_2095", "nb_tox_2095", "accept_ttt_2095", "arret_prema_2095", "analyse_med_2095")] <- 
    simu_topiva(ana_inter_eff = c(30, 51),
                ana_inter_tox = c(rep(5, 4), rep(10, 6), 1),
                liste_essais = tab_donnee, 
                liste_cohortes = cohorte,
                phitox = c(.15, .2), prior_eff = .15, 
                prior_tox = c(1, 1), critere_tox = .95,
                C_ = seuil_effseule[["C_"]], gamm = seuil_effseule[["gamma"]],
                interpatient = 6, max_teff = 180, max_ttox = 42) 
  result_iva[i, c("nb_pts_22595", "temps_etude_22595", "arret_fut_22595", "arret_tox_22595", "nb_eff_22595", "nb_tox_22595", "accept_ttt_22595", "arret_prema_22595", "analyse_med_22595")] <- 
    simu_topiva(ana_inter_eff = c(30, 51),
                ana_inter_tox = c(rep(5, 4), rep(10, 6), 1),
                liste_essais = tab_donnee, 
                liste_cohortes = cohorte,
                phitox = c(.15, .225), prior_eff = .15, 
                prior_tox = c(1, 1), critere_tox = .95,
                C_ = seuil_effseule[["C_"]], gamm = seuil_effseule[["gamma"]],
                interpatient = 6, max_teff = 180, max_ttox = 42) 
  result_iva[i, c("nb_pts_2090", "temps_etude_2090", "arret_fut_2090", "arret_tox_2090", "nb_eff_2090", "nb_tox_2090", "accept_ttt_2090", "arret_prema_2090", "analyse_med_2090")] <- 
    simu_topiva(ana_inter_eff = c(30, 51),
                ana_inter_tox = c(rep(5, 4), rep(10, 6), 1),
                liste_essais = tab_donnee, 
                liste_cohortes = cohorte,
                phitox = c(.15, .2), prior_eff = .15, 
                prior_tox = c(1, 1), critere_tox = .9,
                C_ = seuil_effseule[["C_"]], gamm = seuil_effseule[["gamma"]],
                interpatient = 6, max_teff = 180, max_ttox = 42) 
  result_iva[i, c("nb_pts_22590", "temps_etude_22590", "arret_fut_22590", "arret_tox_22590", "nb_eff_22590", "nb_tox_22590", "accept_ttt_22590", "arret_prema_22590", "analyse_med_22590")] <- 
    simu_topiva(ana_inter_eff = c(30, 51),
                ana_inter_tox = c(rep(5, 4), rep(10, 6), 1),
                liste_essais = tab_donnee, 
                liste_cohortes = cohorte,
                phitox = c(.15, .225), prior_eff = .15, 
                prior_tox = c(1, 1), critere_tox = .9,
                C_ = seuil_effseule[["C_"]], gamm = seuil_effseule[["gamma"]],
                interpatient = 6, max_teff = 180, max_ttox = 42) 
  result_iva[i, c("nb_pts_2590", "temps_etude_2590", "arret_fut_2590", "arret_tox_2590", "nb_eff_2590", "nb_tox_2590", "accept_ttt_2590", "arret_prema_2590", "analyse_med_2590")] <- 
    simu_topiva(ana_inter_eff = c(30, 51),
                ana_inter_tox = c(rep(5, 4), rep(10, 6), 1),
                liste_essais = tab_donnee, 
                liste_cohortes = cohorte,
                phitox = c(.15, .25), prior_eff = .15, 
                prior_tox = c(1, 1), critere_tox = .9,
                C_ = seuil_effseule[["C_"]], gamm = seuil_effseule[["gamma"]],
                interpatient = 6, max_teff = 180, max_ttox = 42) 
  
}

save(result_iva, file = "data/result_iva.Rdata")

## B/ Figures ----

load("data/result_iva.Rdata")
load("data/resultats_tox30_rpos_20211005.Rdata")
tableau_results_iva <- tableau_results %>% 
  select(scenar, probas, R, ends_with("topiva")) %>% 
  left_join(result_iva) %>% 
  mutate(arret_both_topiva = arret_fut_topiva + arret_tox_topiva - (1 - accept_ttt_topiva),
         arret_both_2095 = arret_fut_2095 + arret_tox_2095 - (1 - accept_ttt_2095),
         arret_both_22595 = arret_fut_22595 + arret_tox_22595 - (1 - accept_ttt_22595),
         arret_both_2090 = arret_fut_2090 + arret_tox_2090 - (1 - accept_ttt_2090),
         arret_both_22590 = arret_fut_22590 + arret_tox_22590 - (1 - accept_ttt_22590),
         arret_both_2590 = arret_fut_2590 + arret_tox_2590 - (1 - accept_ttt_2590))
names(tableau_results_iva) <- str_replace_all(names(tableau_results_iva), "nb_pts", "moy_pts")
names(tableau_results_iva) <- str_replace_all(names(tableau_results_iva), "temps_etude", "moy_duree")
names(tableau_results_iva) <- str_replace_all(names(tableau_results_iva), "rejet_h0", "accept_ttt")
tableau_results_iva$scenar_p <- factor(tableau_results_iva$scenar,
                                   levels = rev(c("H0", "H1", "Int1", "Int2", "Paseff1", "Paseff2", "Effint1", "Effint2", "Tox1", "Tox2")),
                                   labels = rev(c("Sc 1 : (0.15;0.30) [H<sub>0</sub>]", "Sc 2 : (0.30;0.20) [H<sub>1</sub>]", "Sc 3: (0.20;0.25)", "Sc 4 : (0.25;0.25)", "Sc 5 : (0.15;0.20)",
                                                  "Sc 6 : (0.10;0.15)", "Sc 7 : (0.20;0.20)", "Sc 8 : (0.20;0.30)", "Sc 9 : (0.30;0.30)", "Sc 10 : (0.40;0.35)")))
tableau_results_iva$scenarios <- factor(tableau_results_iva$scenar,
                                    levels = rev(c("H0", "H1", "Int1", "Int2", "Paseff1", "Paseff2", "Effint1", "Effint2", "Tox1", "Tox2")),
                                    labels = rev(c("Sc 1 :", "Sc 2 :", "Sc 3 :", "Sc 4 :", "Sc 5 :",
                                                   "Sc 6 :", "Sc 7 :", "Sc 8 :", "Sc 9 :", "Sc 10 :")))
tableau_results_iva$probas_plot <- factor(tableau_results_iva$scenar,
                                      levels = rev(c("H0", "H1", "Int1", "Int2", "Paseff1", "Paseff2", "Effint1", "Effint2", "Tox1", "Tox2")),
                                      labels = rev(c("(0.15;0.30) [H<sub>0</sub>]", "(0.30;0.20) [H<sub>1</sub>]", "(0.20;0.25)", "(0.25;0.25)", "(0.15;0.20)",
                                                     "(0.10;0.15)", "(0.20;0.20)", "(0.20;0.30)", "(0.30;0.30)", "(0.40;0.35)")))

plot_opchar_iva <- function(tableau_scenars_plot, metadata, var) {
  
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
    mutate(schemas = factor(schemas, 
                            levels = c("topiva", "22595", "2095", "2590", "22590", "2090"),
                            labels = c("&phi;<sub> tox</sub>=0.25/&tau;=0.95", "&phi;<sub> tox</sub>=0.225/&tau;=0.95", "&phi;<sub> tox</sub>=0.20/&tau;=0.95", "&phi;<sub> tox</sub>=0.25/&tau;=0.90", "&phi;<sub> tox</sub>=0.225/&tau;=0.90", "&phi;<sub> tox</sub>=0.20/&tau;=0.90")))
  
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

plot_accept <- plot_opchar_iva(tableau_results_iva, metadata_eng, "accept_ttt")
ggsave(filename = "figures/rech_seuil/calib_iva.png",
       plot = plot_accept, device = "png",
       height = 14, width = 20)

plot_pts <- plot_opchar_iva(tableau_results_iva, metadata_eng, "moy_pts")
ggsave(filename = "figures/rech_seuil/calib_iva_pts.png",
       plot = plot_pts, device = "png",
       height = 14, width = 20)
