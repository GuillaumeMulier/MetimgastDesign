# -------------------------------------------------------------------- #
# Script des figures pour les résultats des scénarios pour la thèse    #
# Créé le 01/10/2021, modifié le 02/02/2022                            #
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
                    size = 16, 
                    color = "white", fill = "#7888C0", box.color = "#000066",
                    halign = 0.5, linetype = 1, r = unit(3, "pt"), width = unit(0.75, "npc"),
                    padding = margin(2, 0, 1, 0), margin = margin(3, 3, 3, 3)),
                  legend.text = element_markdown()))


# Imports des bases et modifications pour faire les graphiques ----

load("data/resultats_tox30_rpos_20220128.Rdata")
tableau_results <- select(tableau_results, -c(analyse_med_bopiva, analyse_med_topiva, analyse_med_simon,
                                              moy_pts_simon, moy_duree_simon, rejet_h0_simon)) %>% 
  mutate(arret_both_topiva = arret_fut_topiva + arret_tox_topiva - (1 - accept_ttt_topiva),
         arret_both_bopiva = arret_fut_bopiva + arret_tox_bopiva - (1 - accept_ttt_bopiva),
         arret_both_classique = arret_fut_classique + arret_tox_classique - (1 - rejet_h0_classique),
         arret_both_toxrapp = arret_fut_toxrapp + arret_tox_toxrapp - (1 - rejet_h0_toxrapp),
         arret_both_simon = arret_fut_simon + arret_tox_simon - (1 - accept_ttt_simon),
         arret_both_toxprior = arret_fut_toxprior + arret_tox_toxprior - (1 - rejet_h0_toxprior))
names(tableau_results) <- str_replace_all(names(tableau_results), "nb_pts", "moy_pts")
names(tableau_results) <- str_replace_all(names(tableau_results), "temps_etude", "moy_duree")
names(tableau_results) <- str_replace_all(names(tableau_results), "rejet_h0", "accept_ttt")
tableau_results$scenar_p <- factor(tableau_results$scenar,
                                   levels = rev(c("H0", "H1", "Int1", "Int2", "Paseff1", "Paseff2", "Effint1", "Effint2", "Tox1", "Tox2")),
                                   labels = rev(c("Sc 1 : (0.15;0.30)", "Sc 2 : (0.30;0.20)", "Sc 3 : (0.20;0.25)", "Sc 4 : (0.25;0.25)", "Sc 5 : (0.15;0.20)",
                                                  "Sc 6 : (0.10;0.15)", "Sc 7 : (0.20;0.20)", "Sc 8 : (0.20;0.30)", "Sc 9 : (0.30;0.30)", "Sc 10 : (0.40;0.35)")))
tableau_results$scenarios <- factor(tableau_results$scenar,
                                    levels = rev(c("H0", "H1", "Int1", "Int2", "Paseff1", "Paseff2", "Effint1", "Effint2", "Tox1", "Tox2")),
                                    labels = rev(c("Sc 1 :", "Sc 2 :", "Sc 3 :", "Sc 4 :", "Sc 5 :",
                                                   "Sc 6 :", "Sc 7 :", "Sc 8 :", "Sc 9 :", "Sc 10 :")))
tableau_results$probas_plot <- factor(tableau_results$scenar,
                                      levels = rev(c("H0", "H1", "Int1", "Int2", "Paseff1", "Paseff2", "Effint1", "Effint2", "Tox1", "Tox2")),
                                      labels = rev(c("(0.15;0.30)", "(0.30;0.20)", "(0.20;0.25)", "(0.25;0.25)", "(0.15;0.20)",
                                                     "(0.10;0.15)", "(0.20;0.20)", "(0.20;0.30)", "(0.30;0.30)", "(0.40;0.35)")))

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
    select(!starts_with("analyse_med")) %>% 
    select(-c(probas, EffTox, EffNotox, NoeffTox, NoeffNotox)) %>% 
    select(scenar, scenar_p, scenarios, probas_plot, R, everything()) %>% 
    pivot_longer(cols = -c(scenar:R),
                 names_pattern = "(.+_.+)_(.+)",
                 names_to = c("stat", "schemas")) %>% 
    mutate(schemas = factor(schemas, 
                            levels = c("topiva", "bopiva", "simon", "classique", "toxrapp", "toxprior"),
                            labels = c("TOP<sub>eff</sub> + PP<sub>tox</sub>", "BOP<sub>eff</sub> + PP<sub>tox</sub>", "Simon + PP<sub>tox</sub>", "TOP<sub>eff/tox</sub> with 2 analyses", "TOP<sub>eff/tox</sub> with close monitoring of toxicity", "TOP<sub>eff/tox</sub> with update of the prior of toxicity")))
  
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

plot_opchar_facet <- function(tableau_scenars_plot, metadata, var) {
  
  # Mise au format long
  tableau_res <- tableau_scenars_plot %>% 
    select(!matches("C_")) %>% 
    select(!starts_with("analyse_med")) %>% 
    select(!matches("gamma")) %>% 
    select(-c(probas, EffTox, EffNotox, NoeffTox, NoeffNotox)) %>% 
    select(scenar, scenar_p, scenarios, probas_plot, R, everything()) %>% 
    pivot_longer(cols = -c(scenar:nom_corr),
                 names_pattern = "(.+_.+)_(.+)",
                 names_to = c("stat", "schemas")) %>% 
    mutate(schemas = factor(schemas, 
                            levels = c("topiva", "bopiva", "simon", "classique", "toxrapp", "toxprior"),
                            labels = c("TOP<sub>eff</sub> + PP<sub>tox</sub>", "BOP<sub>eff</sub> + PP<sub>tox</sub>", "Simon + PP<sub>tox</sub>", "TOP<sub>eff/tox</sub> with 2 analyses", "TOP<sub>eff/tox</sub> with close monitoring of toxicity", "TOP<sub>eff/tox</sub> with update of the prior of toxicity")))
  
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

# I/ Figures en barplot ----

plot_accept <- plot_opchar(tableau_results, metadata_eng, "accept_ttt")
ggsave(filename = "figures/eval_carac/accept_ttt_30_Rpos.png",
       plot = plot_accept, device = "png",
       height = 14, width = 20)

plot_pts <- plot_opchar(tableau_results, metadata_eng, "moy_pts")
ggsave(filename = "figures/eval_carac/moy_pts_30_Rpos.png",
       plot = plot_pts, device = "png",
       height = 14, width = 20)

plot_duree <- plot_opchar(tableau_results, metadata_eng, "moy_duree")
ggsave(filename = "figures/eval_carac/moy_duree_30_Rpos.png",
       plot = plot_duree, device = "png",
       height = 14, width = 20)

plot_arretprec <- plot_opchar(tableau_results, metadata_eng, "arret_prema")
ggsave(filename = "figures/eval_carac/arret_prema_30_Rpos.png",
       plot = plot_arretprec, device = "png",
       height = 14, width = 20)

rm(plot_accept, plot_arretprec, plot_pts, plot_duree)

# II/ Risque alpha et puissance selon la corrélation ----

load("data/probas_h0h1_20220131.Rdata")

plot_alpha_pos <- tab_probas %>%
  pivot_longer(cols = c(rejet_topiva:rejet_toxprior)) %>%
  mutate(hypothese = factor(hypothese, levels = c("H0", "H1"), labels = c("Scenario 1", "Scenario 2")),
         name = factor(name,
                       levels = c("rejet_topiva", "rejet_bopiva", "rejet_simon", "rejet_classique", "rejet_toxrapp", "rejet_toxprior"),
                       labels = c("TOP<sub>eff</sub> + PP<sub>tox</sub>", "BOP<sub>eff</sub> + PP<sub>tox</sub>", "Simon + PP<sub>tox</sub>", "TOP<sub>eff/tox</sub> with 2 interim analyses", "TOP<sub>eff/tox</sub> with close monitoring of toxicity", "TOP<sub>eff/tox</sub> with update of the prior of toxicity"))) %>%
  ggplot(aes(x = correlation, y = value, color = name)) +
  geom_point() +
  geom_line() +
  geom_hline(data = data.frame(hypothese = c("Scenario 1", "Scenario 2"), intercept = c(.05, .9)), aes(yintercept = intercept), linetype = "dashed") +
  facet_wrap(vars(hypothese), scales = "free_y") +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(color = "Design",
        x = "Correlation between efficacy and toxicity",
        y = "% of acceptation of the treatment")
ggsave(plot = plot_alpha_pos,
       filename = "figures/eval_corr/alpha_puiss_pos.png",
       device = "png", height = 8, width = 12)

load("data/probas_h0h1_Rind_20220131.Rdata")

plot_alpha_ind <- tab_probas %>%
  pivot_longer(cols = c(rejet_topiva:rejet_toxprior)) %>%
  mutate(hypothese = factor(hypothese, levels = c("H0", "H1"), labels = c("Scenario 1", "Scenario 2")),
         name = factor(name,
                       levels = c("rejet_topiva", "rejet_bopiva", "rejet_simon", "rejet_classique", "rejet_toxrapp", "rejet_toxprior"),
                       labels = c("TOP<sub>eff</sub> + PP<sub>tox</sub>", "BOP<sub>eff</sub> + PP<sub>tox</sub>", "Simon + PP<sub>tox</sub>", "TOP<sub>eff/tox</sub> with 2 interim analyses", "TOP<sub>eff/tox</sub> with close monitoring of toxicity", "TOP<sub>eff/tox</sub> with update of the prior of toxicity"))) %>%
  ggplot(aes(x = correlation, y = value, color = name)) +
  geom_point() +
  geom_line() +
  geom_hline(data = data.frame(hypothese = c("Scenario 1", "Scenario 2"), intercept = c(.05, .9)), aes(yintercept = intercept), linetype = "dashed") +
  facet_wrap(vars(hypothese), scales = "free_y") +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(color = "Design",
       x = "Correlation between efficacy and toxicity",
       y = "% of acceptation of the treatment")
ggsave(plot = plot_alpha_ind,
       filename = "figures/eval_corr/alpha_puiss_ind.png",
       device = "png", height = 8, width = 12)
rm(plot_alpha_pos, plot_alpha_ind, tab_probas)

# III/ Figures avec les différentes corrélations ----

## A/ Pour une corrélation positive ----

load("data/resultats_tox30_rposvar_20220201.Rdata")
tableau_results_Rvar <- select(tableau_results_Rvar, -c(analyse_med_bopiva, analyse_med_topiva)) %>% 
  mutate(arret_both_topiva = arret_fut_topiva + arret_tox_topiva - (1 - accept_ttt_topiva),
         arret_both_bopiva = arret_fut_bopiva + arret_tox_bopiva - (1 - accept_ttt_bopiva),
         arret_both_simon = arret_fut_simon + arret_tox_simon - (1 - accept_ttt_simon),
         arret_both_classique = arret_fut_classique + arret_tox_classique - (1 - rejet_h0_classique),
         arret_both_toxrapp = arret_fut_toxrapp + arret_tox_toxrapp - (1 - rejet_h0_toxrapp),
         arret_both_toxprior = arret_fut_toxprior + arret_tox_toxprior - (1 - rejet_h0_toxprior))
names(tableau_results_Rvar) <- str_replace_all(names(tableau_results_Rvar), "nb_pts", "moy_pts")
names(tableau_results_Rvar) <- str_replace_all(names(tableau_results_Rvar), "temps_etude", "moy_duree")
names(tableau_results_Rvar) <- str_replace_all(names(tableau_results_Rvar), "rejet_h0", "accept_ttt")
tableau_results_Rvar$nom_corr <- factor(tableau_results_Rvar$nom_corr,
                                        levels = c("Rmin", "Rneg", "Rind", "Rpos1", "Rpos2", "Rmax"),
                                        labels = c("R<sub>min</sub>", "R<sub>neg</sub>", "R<sub>ind</sub>", "R<sub>pos,1</sub>", "R<sub>pos,2</sub>", "R<sub>max</sub>"))
tableau_results_Rvar$scenar_p <- factor(tableau_results_Rvar$scenar,
                                        levels = rev(c("H0", "H1", "Int1", "Int2", "Paseff1", "Paseff2", "Effint1", "Effint2", "Tox1", "Tox2")),
                                        labels = rev(c("Sc 1 : (0.15;0.30)", "Sc 2 : (0.30;0.20)", "Sc 3: (0.20;0.25)", "Sc 4 : (0.25;0.25)", "Sc 5 : (0.15;0.20)",
                                                       "Sc 6 : (0.10;0.15)", "Sc 7 : (0.20;0.20)", "Sc 8 : (0.20;0.30)", "Sc 9 : (0.30;0.30)", "Sc 10 : (0.40;0.35)")))
tableau_results_Rvar$scenarios <- factor(tableau_results_Rvar$scenar,
                                         levels = rev(c("H0", "H1", "Int1", "Int2", "Paseff1", "Paseff2", "Effint1", "Effint2", "Tox1", "Tox2")),
                                         labels = rev(c("Sc 1 :", "Sc 2 :", "Sc 3 :", "Sc 4 :", "Sc 5 :",
                                                        "Sc 6 :", "Sc 7 :", "Sc 8 :", "Sc 9 :", "Sc 10 :")))
tableau_results_Rvar$probas_plot <- factor(tableau_results_Rvar$scenar,
                                           levels = rev(c("H0", "H1", "Int1", "Int2", "Paseff1", "Paseff2", "Effint1", "Effint2", "Tox1", "Tox2")),
                                           labels = rev(c("(0.15;0.30)", "(0.30;0.20)", "(0.20;0.25)", "(0.25;0.25)", "(0.15;0.20)",
                                                          "(0.10;0.15)", "(0.20;0.20)", "(0.20;0.30)", "(0.30;0.30)", "(0.40;0.35)")))

plot_accept <- plot_opchar(tableau_results_Rvar %>% filter(nom_corr == "R<sub>pos,1</sub>") %>% select(-nom_corr), metadata_eng, "accept_ttt")
ggsave(filename = "figures/eval_corr/accept_ttt_Rpos.png",
       plot = plot_accept, device = "png",
       height = 14, width = 20)

plot_arret_prec <- plot_opchar(tableau_results_Rvar %>% filter(nom_corr == "R<sub>pos,1</sub>") %>% select(-nom_corr), metadata_eng, "arret_prema")
ggsave(filename = "figures/eval_corr/arret_prema_Rpos.png",
       plot = plot_arret_prec, device = "png",
       height = 14, width = 20)

plot_pts <- plot_opchar(tableau_results_Rvar %>% filter(nom_corr == "R<sub>pos,1</sub>") %>% select(-nom_corr), metadata_eng, "moy_pts")
ggsave(filename = "figures/eval_corr/moy_pts_Rpos.png",
       plot = plot_pts, device = "png",
       height = 14, width = 20)

plot_duree <- plot_opchar(tableau_results_Rvar %>% filter(nom_corr == "R<sub>pos,1</sub>") %>% select(-nom_corr), metadata_eng, "moy_duree")
ggsave(filename = "figures/eval_corr/moy_duree_Rpos.png",
       plot = plot_duree, device = "png",
       height = 14, width = 20)
rm(plot_accept, plot_arret_prec, plot_pts, plot_duree)

## B/ En faisant varier la corrélation entre efficacité et toxicité ----

plot_accept <- plot_opchar_facet(tableau_results_Rvar, metadata_eng, "accept_ttt")
ggsave(filename = "figures/eval_corr/accept_ttt_Rvar.png",
       plot = plot_accept, device = "png",
       height = 16, width = 20)

plot_arret_prec <- plot_opchar_facet(tableau_results_Rvar, metadata_eng, "arret_prema")
ggsave(filename = "figures/eval_corr/arret_prema_Rvar.png",
       plot = plot_arret_prec, device = "png",
       height = 16, width = 20)

plot_pts <- plot_opchar_facet(tableau_results_Rvar, metadata_eng, "moy_pts")
ggsave(filename = "figures/eval_corr/moy_pts_Rvar.png",
       plot = plot_pts, device = "png",
       height = 16, width = 20)

plot_duree <- plot_opchar_facet(tableau_results_Rvar, metadata_eng, "moy_duree")
ggsave(filename = "figures/eval_corr/moy_duree_Rvar.png",
       plot = plot_duree, device = "png",
       height = 16, width = 20)
rm(plot_accept, plot_arret_prec, plot_pts, plot_duree)

## C/ Risque alpha et puissance évoluant avec les 6 corrélations choisies ----

plot_alpha_pos <- tableau_results_Rvar %>% 
  filter(scenar %in% c("H0", "H1")) %>% 
  select(R, probas, scenar, starts_with("accept_ttt")) %>%
  pivot_longer(cols = c(accept_ttt_topiva:accept_ttt_toxprior)) %>%
  mutate(scenar = factor(scenar, levels = c("H0", "H1"), labels = c("Scenario 1", "Scenario 2")),
         name = factor(name,
                       levels = c("accept_ttt_topiva", "accept_ttt_bopiva", "accept_ttt_simon", "accept_ttt_classique", "accept_ttt_toxrapp", "accept_ttt_toxprior"),
                       labels = c("TOP<sub>eff</sub> + PP<sub>tox</sub>", 
                                  "BOP<sub>eff</sub> + PP<sub>tox</sub>", 
                                  "Simon + PP<sub>tox</sub>", 
                                  "TOP<sub>eff/tox</sub> with 2 interim analyses", 
                                  "TOP<sub>eff/tox</sub> with close monitoring of toxicity", 
                                  "TOP<sub>eff/tox</sub> with update of the prior of toxicity"))) %>% 
  ggplot(aes(x = R, y = value, color = name)) +
  geom_point() +
  geom_line() +
  geom_hline(data = data.frame(scenar = c("Scenario 1", "Scenario 2"), intercept = c(.05, .9)), aes(yintercept = intercept), linetype = "dashed") +
  facet_wrap(vars(scenar), scales = "free_y") +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(color = "Design",
       x = "Correlation between efficacy and toxicity",
       y = "% of acceptation of the treatment")
ggsave(plot = plot_alpha_pos,
       filename = "figures/eval_corr/alpha_puiss_pos.png",
       device = "png", height = 8, width = 12)




