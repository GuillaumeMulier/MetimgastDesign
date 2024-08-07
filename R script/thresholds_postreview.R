# ------------------------------------------- #
# Rscript for the determination of thresholds #
# ------------------------------------------- #


# I/ Packages and source scripts ----

source("R script/metadata.R")

library(tidyverse)
library(patchwork)
library(cowplot)
library(ggtext)
# devtools::install_github("GuillaumeMulier/multibrasBOP2")
library(multibrasBOP2)

theme_set(theme_light() +
            theme(strip.background = element_rect(fill = NA),
                  strip.text = element_textbox(
                    size = 12, 
                    color = "white", fill = "#7888C0", box.color = "#000066",
                    halign = 0.5, linetype = 1, r = unit(3, "pt"), width = unit(0.75, "npc"),
                    padding = margin(2, 0, 1, 0), margin = margin(3, 3, 3, 3))))

remove <- TRUE # Set to FALSE if you want to keep the different final objects of thresholds

# II/ BOP2 and TOP ----

## A/ Threshold for the design only for efficacy ----

seuil_effseule <- deter_cutoff(alpha = .1,
                      ana_inter = c(30, 51),
                      ana_inter_tox = c(81),
                      p_n = c(0, .15, 0, .85),
                      p_a = c(0, .3, 0, .7),
                      mat_beta_xi = matrix(c(1, 1, 0, 0, 0, 0, 0, 0), nrow = 2, byrow = TRUE),
                      methode = 4L,
                      affich_mat = "No")
save(seuil_effseule, file = paste0("data/cutoff_planifie.Rdata"))
if (remove) rm(seuil)

## B/ Threshold for TOP with 2 analyses in total ----

# Add the thresholds for 2 analyses (additionnal analysis not in the main manuscript)
tableau_scenars <- liste_scenars %>% 
  pull(nom_corr) %>%
  as.character() %>% 
  unique() %>% 
  map_dfr(function(x) {
    seuil <- deter_cutoff(alpha = .05,
                          ana_inter = c(30, 51),
                          ana_inter_tox = c(30, 51),
                          p_n = liste_scenars %>% filter(nom_corr == x, scenar == "H0") %>% pull(probas) %>% .[[1]],
                          p_a = liste_scenars %>% filter(nom_corr == x, scenar == "H1") %>% pull(probas) %>% .[[1]],
                          affich_mat = "No",
                          methode = 4L)[[1]]
    c(nom_corr = x, setNames(seuil, paste0(names(seuil), "_classique")))
  }) %>% 
  right_join(liste_scenars %>% mutate(nom_corr = as.character(nom_corr)), by = "nom_corr")

## C/ Thresholds for TOP with close monitoring of toxicity ----

# Add the thresholds for close monitoring
tableau_scenars <- tableau_scenars %>% 
  pull(nom_corr) %>% 
  as.character() %>% 
  unique() %>% 
  map_dfr(function(x) {
    seuil <- deter_cutoff(alpha = .05,
                          ana_inter = c(30, 51),
                          ana_inter_tox = c(rep(5, 4), rep(10, 5), 11),
                          p_n = tableau_scenars %>% filter(nom_corr == x, scenar == "H0") %>% pull(probas) %>% .[[1]],
                          p_a = tableau_scenars %>% filter(nom_corr == x, scenar == "H1") %>% pull(probas) %>% .[[1]],
                          affich_mat = "No",
                          methode = 4L)[[1]]
    c(nom_corr = x, setNames(seuil, paste0(names(seuil), "_toxrapp")))
  }) %>% 
  right_join(tableau_scenars %>% mutate(nom_corr = as.character(nom_corr)), by = "nom_corr")
save(tableau_scenars, file = paste0("data/scenar", format(Sys.Date(), "%Y%m%d"), ".Rdata"))
if (remove) rm(tableau_scenars)

## D/ Influence of the correlation on computed thresholds ----

# function to plot the heatmap
heatmap_top <- function(data, colonne) {
  labelleur <- c("C_" = "&lambda;", "gamma" = "&gamma;", "alpha_calc" = "&alpha;", "puissance_calc" = "Power")
  data %>% 
    pivot_longer(c(C_, gamma, alpha_calc, puissance_calc)) %>% 
    filter(name == colonne) %>% 
    ggplot(aes(x = correlation_H0, y = correlation_H1, fill = value, label = value)) +
    geom_tile() +
    geom_text(color = "white", size = 4) +
    theme(legend.title = element_markdown()) +
    facet_wrap(vars(name), labeller = labeller(name = labelleur)) +
    scale_x_continuous(name = expression(paste("Correlation Eff/Tox under ", H[0])), expand = c(0, 0)) +
    scale_y_continuous(name = expression(paste("Correlation Eff/Tox under ", H[1])), expand = c(0, 0)) +
    labs(fill = labs(fill = labelleur[colonne]))
}

# TOP with 2 analyses
corr_2 <- generer_corr_multi(seq(0, .15, .015), .15, .3)
corr_3 <- generer_corr_multi(seq(0, .2, .02), .3, .2)
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
                                       methode = 4L,
                                       seed = 1993,
                                       affich_mat = "No")[[1]]
}
table_resume_2 <- cbind(tab_probas, do.call(rbind, liste_resultats))
save(table_resume_2, file = paste0("data/essais_cutoff_bis_", format(Sys.Date(), "%Y%m%d"), ".Rdata"))
load("data/essais_cutoff_bis_20240522.Rdata")
heatmap_seuils_30 <- map(c("C_", "gamma", "alpha_calc", "puissance_calc"), heatmap_top, data = table_resume_2) %>% 
  wrap_plots()
ggsave(plot = heatmap_seuils_30,
       filename = "figures/threshold/heatmap_seuils_30eng.png",
       device = "png", height = 8, width = 16)
ggsave(plot = heatmap_seuils_30,
       filename = "figures/threshold/heatmap_seuils_30eng.eps",
       device = cairo_ps, height = 7000, width = 13000, units = "px", dpi = 800)
rm(liste_resultats, tab_probas, table_resume_2, heatmap_seuils_30, corr_2, corr_3)

# TOP with close monitoring of toxicity
corr_2 <- generer_corr_multi(seq(0, .15, .015), .15, .3)
corr_3 <- generer_corr_multi(seq(0, .2, .02), .3, .2)
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
                                       methode = 4L,
                                       seed = 1993,
                                       affich_mat = "No")[[1]]
}
table_resume_tox_2 <- cbind(tab_probas, do.call(rbind, liste_resultats))
save(table_resume_tox_2, file = paste0("~/essais_cutoff_tox_bis_", format(Sys.Date(), "%Y%m%d"), ".Rdata"))
load("data/essais_cutoff_tox_bis_20210830.Rdata")
heatmap_seuils_30 <- map(c("C_", "gamma", "alpha_calc", "puissance_calc"), heatmap_top, data = table_resume_tox_2) %>% 
  wrap_plots()
ggsave(plot = heatmap_seuils_30,
       filename = "figures/threshold/heatmap_seuils_tox_30eng.png",
       device = "png", height = 8, width = 16)
ggsave(plot = heatmap_seuils_30,
       filename = "figures/threshold/heatmap_seuils_tox_30eng.eps",
       device = cairo_ps, height = 7000, width = 13000, units = "px", dpi = 800)
rm(liste_resultats, tab_probas, table_resume_2, heatmap_seuils_30, corr_2, corr_3)


# III/ TOP + Ivanova ----

## A/ Only Ivanova ----

set.seed(1993)

# Generate the data
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

# Simulate the decisions
table_choix <- expand.grid(
  pi = seq(.2, .3, .01),
  tau = seq(.86, .99, .005),
  alpha = NA_real_,
  puissance = NA_real_
) 
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
save(table_choix, file = "data/alpha_puiss_iva_tox.Rdata")

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
       filename = "figures/threshold/recherche_seuil_iva_tox.png",
       device = "png", height = 9, width = 16)
ggsave(plot = plot_tot,
       filename = "figures/threshold/recherche_seuil_iva_tox.eps",
       device = "eps", height = 7000, width = 14000, units = "px", dpi = 800)
if (remove) rm(donnees_h0, donnees_h1, plot1, plot2, plot_tot, table_choix)

## B/ TOP + Ivanova ----
# = Search a good decision rule

# Simulate results overs the 10 scenarios with the given values of thresholds
result_iva <- liste_scenars %>%
  mutate(nom_corr = str_remove(str_remove(nom_corr, "<sub>"), "</sub>")) %>% 
  filter(nom_corr == "Rpos,1") %>% 
  select(scenar, R, probas)
result_iva[, c("nb_pts_2595", "temps_etude_2595", "arret_fut_2595", "arret_tox_2595", "nb_eff_2595", "nb_tox_2595", "accept_ttt_2595",
               "nb_pts_2095", "temps_etude_2095", "arret_fut_2095", "arret_tox_2095", "nb_eff_2095", "nb_tox_2095", "accept_ttt_2095",
               "nb_pts_22595", "temps_etude_22595", "arret_fut_22595", "arret_tox_22595", "nb_eff_22595", "nb_tox_22595", "accept_ttt_22595",
               "nb_pts_2090", "temps_etude_2090", "arret_fut_2090", "arret_tox_2090", "nb_eff_2090", "nb_tox_2090", "accept_ttt_2090",
               "nb_pts_22590", "temps_etude_22590", "arret_fut_22590", "arret_tox_22590", "nb_eff_22590", "nb_tox_22590", "accept_ttt_22590",
               "nb_pts_2590", "temps_etude_2590", "arret_fut_2590", "arret_tox_2590", "nb_eff_2590", "nb_tox_2590", "accept_ttt_2590")] <- NA_real_
result_iva[, c("analyse_med_2095", "analyse_med_22595", "analyse_med_2090", "analyse_med_22590", "analyse_med_2590")] <- NA_character_
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
  result_iva[i, c("nb_pts_2595", "temps_etude_2595", "arret_fut_2595", "arret_tox_2595", "nb_eff_2595", "nb_tox_2595", "accept_ttt_2595", "arret_prema_2595", "analyse_med_2095")] <- 
    simu_topiva(ana_inter_eff = c(30, 51),
                ana_inter_tox = c(rep(5, 4), rep(10, 6), 1),
                liste_essais = tab_donnee, 
                liste_cohortes = cohorte,
                phieff = .15, phitox = .25, prior_eff = .15, 
                prior_tox = c(1, 1), critere_tox = .95,
                C_ = seuil_effseule[[1]][["C_"]], gamm = seuil_effseule[[1]][["gamma"]],
                interpatient = 6, max_teff = 180, max_ttox = 42) 
  result_iva[i, c("nb_pts_2095", "temps_etude_2095", "arret_fut_2095", "arret_tox_2095", "nb_eff_2095", "nb_tox_2095", "accept_ttt_2095", "arret_prema_2095", "analyse_med_2095")] <- 
    simu_topiva(ana_inter_eff = c(30, 51),
                ana_inter_tox = c(rep(5, 4), rep(10, 6), 1),
                liste_essais = tab_donnee, 
                liste_cohortes = cohorte,
                phieff = .15, phitox = .2, prior_eff = .15, 
                prior_tox = c(1, 1), critere_tox = .95,
                C_ = seuil_effseule[[1]][["C_"]], gamm = seuil_effseule[[1]][["gamma"]],
                interpatient = 6, max_teff = 180, max_ttox = 42) 
  result_iva[i, c("nb_pts_22595", "temps_etude_22595", "arret_fut_22595", "arret_tox_22595", "nb_eff_22595", "nb_tox_22595", "accept_ttt_22595", "arret_prema_22595", "analyse_med_22595")] <- 
    simu_topiva(ana_inter_eff = c(30, 51),
                ana_inter_tox = c(rep(5, 4), rep(10, 6), 1),
                liste_essais = tab_donnee, 
                liste_cohortes = cohorte,
                phieff = .15, phitox = .225, prior_eff = .15, 
                prior_tox = c(1, 1), critere_tox = .95,
                C_ = seuil_effseule[[1]][["C_"]], gamm = seuil_effseule[[1]][["gamma"]],
                interpatient = 6, max_teff = 180, max_ttox = 42) 
  result_iva[i, c("nb_pts_2090", "temps_etude_2090", "arret_fut_2090", "arret_tox_2090", "nb_eff_2090", "nb_tox_2090", "accept_ttt_2090", "arret_prema_2090", "analyse_med_2090")] <- 
    simu_topiva(ana_inter_eff = c(30, 51),
                ana_inter_tox = c(rep(5, 4), rep(10, 6), 1),
                liste_essais = tab_donnee, 
                liste_cohortes = cohorte,
                phieff = .15, phitox = .2, prior_eff = .15, 
                prior_tox = c(1, 1), critere_tox = .9,
                C_ = seuil_effseule[[1]][["C_"]], gamm = seuil_effseule[[1]][["gamma"]],
                interpatient = 6, max_teff = 180, max_ttox = 42) 
  result_iva[i, c("nb_pts_22590", "temps_etude_22590", "arret_fut_22590", "arret_tox_22590", "nb_eff_22590", "nb_tox_22590", "accept_ttt_22590", "arret_prema_22590", "analyse_med_22590")] <- 
    simu_topiva(ana_inter_eff = c(30, 51),
                ana_inter_tox = c(rep(5, 4), rep(10, 6), 1),
                liste_essais = tab_donnee, 
                liste_cohortes = cohorte,
                phitox = .225, phieff = .15, prior_eff = .15, 
                prior_tox = c(1, 1), critere_tox = .9,
                C_ = seuil_effseule[[1]][["C_"]], gamm = seuil_effseule[[1]][["gamma"]],
                interpatient = 6, max_teff = 180, max_ttox = 42) 
  result_iva[i, c("nb_pts_2590", "temps_etude_2590", "arret_fut_2590", "arret_tox_2590", "nb_eff_2590", "nb_tox_2590", "accept_ttt_2590", "arret_prema_2590", "analyse_med_2590")] <- 
    simu_topiva(ana_inter_eff = c(30, 51),
                ana_inter_tox = c(rep(5, 4), rep(10, 6), 1),
                liste_essais = tab_donnee, 
                liste_cohortes = cohorte,
                phieff = .15, phitox = .25, prior_eff = .15, 
                prior_tox = c(1, 1), critere_tox = .9,
                C_ = seuil_effseule[[1]][["C_"]], gamm = seuil_effseule[[1]][["gamma"]],
                interpatient = 6, max_teff = 180, max_ttox = 42) 
}
names(result_iva) <- str_replace_all(names(result_iva), "nb_pts", "moy_pts")
names(result_iva) <- str_replace_all(names(result_iva), "temps_etude", "moy_duree")
names(result_iva) <- str_replace_all(names(result_iva), "rejet_h0", "accept_ttt")
save(result_iva, file = "data/result_iva.Rdata")
result_iva <- mutate(result_iva,
                     arret_both_2595 = arret_fut_2595 + arret_tox_2595 - (1 - accept_ttt_2595),
                     arret_both_2095 = arret_fut_2095 + arret_tox_2095 - (1 - accept_ttt_2095),
                     arret_both_22595 = arret_fut_22595 + arret_tox_22595 - (1 - accept_ttt_22595),
                     arret_both_2090 = arret_fut_2090 + arret_tox_2090 - (1 - accept_ttt_2090),
                     arret_both_22590 = arret_fut_22590 + arret_tox_22590 - (1 - accept_ttt_22590),
                     arret_both_2590 = arret_fut_2590 + arret_tox_2590 - (1 - accept_ttt_2590))
plot_accept <- plot_opchar(result_iva %>% left_join(liste_scenars %>% filter(nom_corr == "R<sub>pos,1</sub>") %>% select(scenar, scenar_p, scenarios, probas_plot)), metadata_eng, "accept_ttt", 
                           niveaux = c("2595", "22595", "2095", "2590", "22590", "2090"),
                           etiquettes = c("&phi;<sub> tox</sub>=0.25/&tau;=0.95", "&phi;<sub> tox</sub>=0.225/&tau;=0.95", "&phi;<sub> tox</sub>=0.20/&tau;=0.95", 
                                          "&phi;<sub> tox</sub>=0.25/&tau;=0.90", "&phi;<sub> tox</sub>=0.225/&tau;=0.90", "&phi;<sub> tox</sub>=0.20/&tau;=0.90"))
ggsave(filename = "figures/threshold/calib_iva.png",
       plot = plot_accept, device = "png",
       height = 14, width = 20)
ggsave(filename = "figures/threshold/calib_iva.eps",
       plot = plot_accept, device = cairo_ps,
       height = 9000, width = 17000, units = "px", dpi = 800)
plot_pts <- plot_opchar(result_iva %>% left_join(liste_scenars %>% filter(nom_corr == "R<sub>pos,1</sub>") %>% select(scenar, scenar_p, scenarios, probas_plot)), metadata_eng, "moy_pts",
                        niveaux = c("2595", "22595", "2095", "2590", "22590", "2090"),
                        etiquettes = c("&phi;<sub> tox</sub>=0.25/&tau;=0.95", "&phi;<sub> tox</sub>=0.225/&tau;=0.95", "&phi;<sub> tox</sub>=0.20/&tau;=0.95", 
                                       "&phi;<sub> tox</sub>=0.25/&tau;=0.90", "&phi;<sub> tox</sub>=0.225/&tau;=0.90", "&phi;<sub> tox</sub>=0.20/&tau;=0.90"))
ggsave(filename = "figures/threshold/calib_iva_pts.png",
       plot = plot_pts, device = "png",
       height = 14, width = 20)
ggsave(filename = "figures/threshold/calib_iva_pts.eps",
       plot = plot_pts, device = cairo_ps,
       height = 9000, width = 17000, units = "px", dpi = 800)
if (remove) rm(result_iva, plot_pts, plot_accept)


# IV/ Simon ----

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
result_simu_simon %>% filter(alpha <= .1, puissance >= .9) %>% arrange(EN_p0)

alpha_puiss_simon(3, 23, 11, 55, .15, .3)
