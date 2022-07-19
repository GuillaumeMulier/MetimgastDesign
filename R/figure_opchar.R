# -------------------------------------------------------------------- #
# Script des figures pour les résultats des simulations de l'influence #
# de la corrélation entre efficacité et toxicité sur les               #
# operating characteristics                                            #
# Créé le 08/09/2021, modifié le 18/10/2021                            #
# -------------------------------------------------------------------- #

# Packages et helpers ----

library(tidyverse)
library(patchwork)
library(rlang)
library(ggtext)
library(cowplot)

theme_set(theme_light() +
            theme(strip.background = element_rect(fill = NA),
                  strip.text = element_textbox(
                    size = 12, 
                    color = "white", fill = "#7888C0", box.color = "#000066",
                    halign = 0.5, linetype = 1, r = unit(3, "pt"), width = unit(0.75, "npc"),
                    padding = margin(2, 0, 1, 0), margin = margin(3, 3, 3, 3))))

# Imports des bases et modifications pour faire les graphiques ----

load("data/resultats_tox30_20210907.Rdata")
tableau_scenars <- select(tableau_scenars, -c(analyse_med_bopiva, analyse_med_topiva)) %>% 
  mutate(arret_both_topiva = arret_fut_topiva + arret_tox_topiva - (1 - accept_ttt_topiva),
         arret_both_bopiva = arret_fut_bopiva + arret_tox_bopiva - (1 - accept_ttt_bopiva),
         arret_both_classique = arret_fut_classique + arret_tox_classique - (1 - rejet_h0_classique),
         arret_both_toxrapp = arret_fut_toxrapp + arret_tox_toxrapp - (1 - rejet_h0_toxrapp))
names(tableau_scenars) <- str_replace_all(names(tableau_scenars), "nb_pts", "moy_pts")
names(tableau_scenars) <- str_replace_all(names(tableau_scenars), "temps_etude", "moy_duree")
names(tableau_scenars) <- str_replace_all(names(tableau_scenars), "rejet_h0", "accept_ttt")
tableau_scenars$nom_corr <- factor(tableau_scenars$nom_corr,
                                   levels = c("Rmin", "Rneg", "Rind", "Rpos1", "Rpos2", "Rmax"),
                                   labels = c("R<sub>min</sub>", "R<sub>neg</sub>", "R<sub>ind</sub>", "R<sub>pos,1</sub>", "R<sub>pos,2</sub>", "R<sub>max</sub>"))
tableau_scenars$scenar_p <- factor(tableau_scenars$scenar,
                                   levels = rev(c("H0", "H1", "Int1", "Int2", "Paseff1", "Paseff2", "Effint1", "Effint2", "Tox1", "Tox2")),
                                   labels = rev(c("Sc 1 : (0.15;0.30) [H<sub>0</sub>]", "Sc 2 : (0.30;0.20) [H<sub>1</sub>]", "Sc 3: (0.20;0.25)", "Sc 4 : (0.25;0.25)", "Sc 5 : (0.15;0.20)",
                                                  "Sc 6 : (0.10;0.15)", "Sc 7 : (0.20;0.20)", "Sc 8 : (0.20;0.30)", "Sc 9 : (0.30;0.30)", "Sc 10 : (0.40;0.35)")))
tableau_scenars$scenarios <- factor(tableau_scenars$scenar,
                                    levels = rev(c("H0", "H1", "Int1", "Int2", "Paseff1", "Paseff2", "Effint1", "Effint2", "Tox1", "Tox2")),
                                    labels = rev(c("Sc 1 :", "Sc 2 :", "Sc 3 :", "Sc 4 :", "Sc 5 :",
                                                   "Sc 6 :", "Sc 7 :", "Sc 8 :", "Sc 9 :", "Sc 10 :")))
tableau_scenars$probas_plot <- factor(tableau_scenars$scenar,
                                      levels = rev(c("H0", "H1", "Int1", "Int2", "Paseff1", "Paseff2", "Effint1", "Effint2", "Tox1", "Tox2")),
                                      labels = rev(c("(0.15;0.30) [H<sub>0</sub>]", "(0.30;0.20) [H<sub>1</sub>]", "(0.20;0.25)", "(0.25;0.25)", "(0.15;0.20)",
                                                     "(0.10;0.15)", "(0.20;0.20)", "(0.20;0.30)", "(0.30;0.30)", "(0.40;0.35)")))

load("data/resultats_tox40_20210908.Rdata")
tableau_scenars_bis <- select(tableau_scenars_bis, -c(analyse_med_bopiva, analyse_med_topiva)) %>% 
  mutate(arret_both_topiva = arret_fut_topiva + arret_tox_topiva - (1 - accept_ttt_topiva),
         arret_both_bopiva = arret_fut_bopiva + arret_tox_bopiva - (1 - accept_ttt_bopiva),
         arret_both_classique = arret_fut_classique + arret_tox_classique - (1 - rejet_h0_classique),
         arret_both_toxrapp = arret_fut_toxrapp + arret_tox_toxrapp - (1 - rejet_h0_toxrapp))
names(tableau_scenars_bis) <- str_replace_all(names(tableau_scenars_bis), "nb_pts", "moy_pts")
names(tableau_scenars_bis) <- str_replace_all(names(tableau_scenars_bis), "temps_etude", "moy_duree")
names(tableau_scenars_bis) <- str_replace_all(names(tableau_scenars_bis), "rejet_h0", "accept_ttt")
tableau_scenars_bis$nom_corr <- factor(tableau_scenars_bis$nom_corr,
                                   levels = c("Rmin", "Rneg", "Rind", "Rpos1", "Rpos2", "Rmax"),
                                   labels = c("R<sub>min</sub>", "R<sub>neg</sub>", "R<sub>ind</sub>", "R<sub>pos,1</sub>", "R<sub>pos,2</sub>", "R<sub>max</sub>"))
tableau_scenars_bis$scenar_p <- factor(tableau_scenars_bis$scenar,
                                   levels = rev(c("H0", "H1", "Int1", "Int2", "Paseff1", "Paseff2", "Effint1", "Effint2", "Tox1", "Tox2")),
                                   labels = rev(c("Sc 1 : (0.15;0.40) [H<sub>0</sub>]", "Sc 2 : (0.30;0.20) [H<sub>1</sub>]", "Sc 3: (0.20;0.30)", "Sc 4 : (0.25;0.35)", "Sc 5 : (0.15;0.20)",
                                                  "Sc 6 : (0.10;0.15)", "Sc 7 : (0.20;0.20)", "Sc 8 : (0.20;0.40)", "Sc 9 : (0.30;0.40)", "Sc 10 : (0.40;0.45)")))
tableau_scenars_bis$scenarios <- factor(tableau_scenars_bis$scenar,
                                    levels = rev(c("H0", "H1", "Int1", "Int2", "Paseff1", "Paseff2", "Effint1", "Effint2", "Tox1", "Tox2")),
                                    labels = rev(c("Sc 1 :", "Sc 2 :", "Sc 3 :", "Sc 4 :", "Sc 5 :",
                                                   "Sc 6 :", "Sc 7 :", "Sc 8 :", "Sc 9 :", "Sc 10 :")))
tableau_scenars_bis$probas_plot <- factor(tableau_scenars_bis$scenar,
                                      levels = rev(c("H0", "H1", "Int1", "Int2", "Paseff1", "Paseff2", "Effint1", "Effint2", "Tox1", "Tox2")),
                                      labels = rev(c("(0.15;0.40) [H<sub>0</sub>]", "(0.30;0.20) [H<sub>1</sub>]", "(0.20;0.30)", "(0.25;0.35)", "(0.15;0.20)",
                                                     "(0.10;0.15)", "(0.20;0.20)", "(0.20;0.40)", "(0.30;0.40)", "(0.40;0.45)")))

metadata <- tribble(
  ~var, ~label, ~lim_sup, ~breaks, ~format_legende,
  "accept_ttt", "Proportion d'acceptation du traitement", 1, expression(seq(0, 1, .2)), expression(scales::percent_format()),
  "arret_fut", "Proportion d'arrêt pour futilité", 1, expression(seq(0, 1, .2)), expression(scales::percent_format()),
  "arret_tox", "Proportion d'arrêt pour toxicité", 1, expression(seq(0, 1, .2)), expression(scales::percent_format()),
  "moy_pts", "Nombre moyen de patients par essai", 85, expression(seq(0, 80, 20)), expression(function(x) paste0(x, c(" pt", rep(" pts", length(x) - 1)))),
  "moy_duree", "Durée moyenne de l'étude", 2400, expression(seq(0, 2400, 400)), expression(function(x) paste0(round(x / 365.25, 2), c(" an", rep(" ans", length(x) - 1)))),
  "nb_eff", "Nombre de réponses en moyenne par essai", 30, expression(seq(0, 30, 5)), expression(function(x) x),
  "nb_tox", "Nombre de toxicités en moyenne par essai", 30, expression(seq(0, 30, 5)), expression(function(x) x)
)

metadata_eng <- tribble(
  ~var, ~label, ~lim_sup, ~breaks, ~format_legende,
  "accept_ttt", "% of acceptation of treatment", 1, expression(seq(0, 1, .2)), expression(scales::percent_format()),
  "arret_fut", "% futility stopping", 1, expression(seq(0, 1, .2)), expression(scales::percent_format()),
  "arret_tox", "% toxicity stopping", 1, expression(seq(0, 1, .2)), expression(scales::percent_format()),
  "moy_pts", "Mean number of patients by trial", 85, expression(seq(0, 80, 20)), expression(function(x) paste0(x, c(" pt", rep(" pts", length(x) - 1)))),
  "moy_duree", "Mean duration of a trial", 2400, expression(seq(0, 2400, 400)), expression(function(x) paste0(round(x / 365.25, 2), c(" an", rep(" ans", length(x) - 1)))),
  "nb_eff", "Mean number of responses per trial", 30, expression(seq(0, 30, 5)), expression(function(x) x),
  "nb_tox", "Mean number of toxicities per trial", 30, expression(seq(0, 30, 5)), expression(function(x) x)
)

plot_opchar <- function(tableau_scenars_plot, metadata, var) {
  
  # Mise au format long
  tableau_res <- tableau_scenars_plot %>% 
    select(!matches("C_")) %>% 
    select(!matches("gamma")) %>% 
    select(-c(probas, EffTox, EffNotox, NoeffTox, NoeffNotox)) %>% 
    select(nom_corr, scenar, scenar_p, scenarios, probas_plot, R, everything()) %>% 
    pivot_longer(cols = -c(nom_corr:R),
                 names_pattern = "(.+_.+)_(.+)",
                 names_to = c("stat", "schemas")) %>% 
    mutate(schemas = factor(schemas, 
                            levels = c("topiva", "bopiva", "classique", "toxrapp"),
                            labels = c("TOP + Ivanova", "BOP + Ivanova", "TOP avec 2 analyses intermédiaires", "TOP avec monitoring rapproché de la toxicité")))
  
  layout <- 
    "AB
     CC"
  
  # Barplot
  plot1 <- tableau_res %>% 
    filter(stat == var) %>% 
    ggplot(aes(x = scenar_p, y = value, fill = schemas)) +
    geom_col(position = position_dodge2()) +
    geom_vline(xintercept = c(seq(1.5, 9.5, 1)), color = "grey50", 
               size = 1, linetype = "dashed") +
    coord_flip() +
    facet_wrap(vars(nom_corr), ncol = 2, dir = "v") +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.border = element_rect(color = "transparent"),
          legend.position = "bottom",
          legend.text = element_text(size = 13),
          legend.title = element_text(size = 16, face = "bold"))  +
    scale_y_continuous(name = metadata$label[match(var, metadata$var)],
                       limits = c(0, metadata$lim_sup[match(var, metadata$var)]),
                       labels = eval(metadata$format_legende[[match(var, metadata$var)]]),
                       breaks = eval(metadata$breaks[[match(var, metadata$var)]])) +
    labs(x = NULL,
         fill = "Type de schéma d'essai")
  legende <- get_legend(plot1)
  
  # Annotations pour les scénarios
  plot2 <- tableau_scenars_plot %>% 
    filter(nom_corr %in% c("R<sub>min</sub>", "R<sub>neg</sub>", "R<sub>ind</sub>")) %>% 
    ggplot() +
    geom_richtext(aes(x = scenar_p, y = 1.1, label = scenarios), hjust = 1, fill = "transparent", label.color = "transparent") +
    geom_richtext(aes(x = scenar_p, y = 1.3, label = probas_plot), hjust = 0, fill = "transparent", label.color = "transparent") +
    geom_vline(xintercept = c(seq(1.5, 9.5, 1)), color = "grey50", 
               size = 1, linetype = "dashed") +
    coord_flip(clip = "off") +
    ylim(c(.8, 1.7)) +
    theme_void() +
    facet_wrap(vars(nom_corr), ncol = 2, nrow = 3, dir = "v") +
    theme(strip.background = element_rect(fill = NA, color = NA),
          strip.text = element_textbox(
            size = 12, 
            color = "transparent", fill = "transparent", box.color = "transparent",
            halign = 0.5, linetype = 1, r = unit(3, "pt"), width = unit(.001, "in"),
            padding = margin(2, 0, 1, 0), margin = margin(3, 3, 3, 3)))
  
  plot_tot <- plot2 + (plot1 + theme(legend.position = "none")) + legende +
    plot_layout(design = layout, widths = c(1, 7), heights = c(10, 1))
  
  return(plot_tot)
  
}

# I/ Figures pour le cas avec 30% de toxicité comme hypothèse nulle ----

plot_rejet <- plot_opchar(tableau_scenars, metadata_eng, "accept_ttt")
# ggsave(plot = plot_rejet,
#        filename = "figures/eval_carac/accept_ttt_30.png",
#        device = "png", width = 20, height = 12)
rm(plot_rejet)


plot_opchar(tableau_scenars, metadata, "moy_pts")
plot_opchar(tableau_scenars, metadata, "moy_duree")
plot_opchar(tableau_scenars, metadata, "arret_fut")
plot_opchar(tableau_scenars, metadata, "arret_tox")
plot_opchar(tableau_scenars, metadata, "nb_eff")
plot_opchar(tableau_scenars, metadata, "nb_tox")

# II/ Avec 40% comme hypothèse nulle de toxicité ----

plot_opchar(tableau_scenars_bis, metadata, "accept_ttt")
plot_opchar(tableau_scenars_bis, metadata, "moy_pts")
plot_opchar(tableau_scenars_bis, metadata, "moy_duree")
plot_opchar(tableau_scenars_bis, metadata, "arret_fut")
plot_opchar(tableau_scenars_bis, metadata, "arret_tox")
plot_opchar(tableau_scenars_bis, metadata, "nb_eff")
plot_opchar(tableau_scenars_bis, metadata, "nb_tox")