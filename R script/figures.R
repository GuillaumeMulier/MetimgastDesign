# -------------------------------------- #
# Rscript for the figures of the results #
# -------------------------------------- #

library(tidyverse)
library(ggtext)
library(cowplot)
library(patchwork)

source("R script/metadata.R")
source("R script/functions.R")
load("data/resultats_tox30_rposvar_20240524.Rdata")

theme_set(theme_light() +
            theme(strip.background = element_rect(fill = NA),
                  strip.text = element_textbox(
                    size = 12, 
                    color = "white", fill = "#7888C0", box.color = "#000066",
                    halign = 0.5, linetype = 1, r = unit(3, "pt"), width = unit(0.75, "npc"),
                    padding = margin(2, 0, 1, 0), margin = margin(3, 3, 3, 3))))

# Clean the table and add some new variables
tableau_results_Rvar <- select(tableau_results_Rvar, -starts_with("analyse_med")) %>% 
  mutate(arret_both_topiva = arret_fut_topiva + arret_tox_topiva - (1 - accept_ttt_topiva),
         arret_both_bopiva = arret_fut_bopiva + arret_tox_bopiva - (1 - accept_ttt_bopiva),
         arret_both_simon = arret_fut_simon + arret_tox_simon - (1 - accept_ttt_simon),
         arret_both_classique = arret_fut_classique + arret_tox_classique - (1 - rejet_h0_classique),
         arret_both_toxrapp = arret_fut_toxrapp + arret_tox_toxrapp - (1 - rejet_h0_toxrapp),
         arret_both_toxprior = arret_fut_toxprior + arret_tox_toxprior - (1 - rejet_h0_toxprior))
names(tableau_results_Rvar) <- str_replace_all(names(tableau_results_Rvar), "nb_pts", "moy_pts")
names(tableau_results_Rvar) <- str_replace_all(names(tableau_results_Rvar), "temps_etude", "moy_duree")
names(tableau_results_Rvar) <- str_replace_all(names(tableau_results_Rvar), "rejet_h0", "accept_ttt")
tableau_results_Rvar <- left_join(tableau_results_Rvar %>% select(-nom_corr), liste_scenars)


# I/ With positive correlation between efficacy and toxicity ----

plot_accept <- plot_opchar(tableau_results_Rvar %>% filter(nom_corr == "R<sub>pos,1</sub>") %>% select(-nom_corr), metadata_eng, "accept_ttt")
plot_accept <- plot_opchar(tableau_results_Rvar %>% filter(nom_corr == "R<sub>pos,1</sub>") %>% select(-nom_corr), metadata_eng, "accept_ttt",
                           niveaux = c("topiva", "bopiva", "simon", "toxrapp", "toxprior"),
                           etiquettes = c("TOP<sub>eff</sub> + PP<sub>tox</sub>", "BOP<sub>eff</sub> + PP<sub>tox</sub>",
                                          "Simon + PP<sub>tox</sub>", 
                                          "TOP<sup>t</sup><sub>eff/tox</sub>",
                                          "iTOP<sub>eff/tox</sub>"))
ggsave(filename = "figures/results/accept_ttt_Rpos_rev2.png",
       plot = plot_accept, device = "png",
       height = 14, width = 20)
ggsave(filename = "figures/results/accept_ttt_Rpos_rev2.eps",
       plot = plot_accept, device = cairo_ps,
       height = 9000, width = 17000, units = "px", dpi = 800)

plot_arret_prec <- plot_opchar(tableau_results_Rvar %>% filter(nom_corr == "R<sub>pos,1</sub>") %>% select(-nom_corr), metadata_eng, "arret_prema")
plot_arret_prec <- plot_opchar(tableau_results_Rvar %>% filter(nom_corr == "R<sub>pos,1</sub>") %>% select(-nom_corr), metadata_eng, "arret_prema",
                               niveaux = c("topiva", "bopiva", "simon", "toxrapp", "toxprior"),
                               etiquettes = c("TOP<sub>eff</sub> + PP<sub>tox</sub>", "BOP<sub>eff</sub> + PP<sub>tox</sub>",
                                              "Simon + PP<sub>tox</sub>", 
                                              "TOP<sup>t</sup><sub>eff/tox</sub>",
                                              "iTOP<sub>eff/tox</sub>"))
ggsave(filename = "figures/results/arret_prema_Rpos_rev2.png",
       plot = plot_arret_prec, device = "png",
       height = 14, width = 20)
ggsave(filename = "figures/results/arret_prema_Rpos_rev2.eps",
       plot = plot_arret_prec, device = cairo_ps,
       height = 9000, width = 17000, units = "px", dpi = 800)

plot_pts <- plot_opchar(tableau_results_Rvar %>% filter(nom_corr == "R<sub>pos,1</sub>") %>% select(-nom_corr), metadata_eng, "moy_pts")
plot_pts <- plot_opchar(tableau_results_Rvar %>% filter(nom_corr == "R<sub>pos,1</sub>") %>% select(-nom_corr), metadata_eng, "moy_pts",
                        niveaux = c("topiva", "bopiva", "simon", "toxrapp", "toxprior"),
                        etiquettes = c("TOP<sub>eff</sub> + PP<sub>tox</sub>", "BOP<sub>eff</sub> + PP<sub>tox</sub>",
                                       "Simon + PP<sub>tox</sub>", 
                                       "TOP<sup>t</sup><sub>eff/tox</sub>",
                                       "iTOP<sub>eff/tox</sub>"))
ggsave(filename = "figures/results/moy_pts_Rpos_rev2.png",
       plot = plot_pts, device = "png",
       height = 14, width = 20)
ggsave(filename = "figures/results/moy_pts_Rpos_rev2.eps",
       plot = plot_pts, device = cairo_ps,
       height = 9000, width = 17000, units = "px", dpi = 800)

plot_duree <- plot_opchar(tableau_results_Rvar %>% filter(nom_corr == "R<sub>pos,1</sub>") %>% select(-nom_corr), metadata_eng, "moy_duree")
plot_duree <- plot_opchar(tableau_results_Rvar %>% filter(nom_corr == "R<sub>pos,1</sub>") %>% select(-nom_corr), metadata_eng, "moy_duree",
                          niveaux = c("topiva", "bopiva", "simon", "toxrapp", "toxprior"),
                          etiquettes = c("TOP<sub>eff</sub> + PP<sub>tox</sub>", "BOP<sub>eff</sub> + PP<sub>tox</sub>",
                                         "Simon + PP<sub>tox</sub>", 
                                         "TOP<sup>t</sup><sub>eff/tox</sub>",
                                         "iTOP<sub>eff/tox</sub>"))
ggsave(filename = "figures/results/moy_duree_Rpos_rev2.png",
       plot = plot_duree, device = "png",
       height = 14, width = 20)
ggsave(filename = "figures/results/moy_duree_Rpos_rev2.eps",
       plot = plot_duree, device = cairo_ps,
       height = 9000, width = 17000, units = "px", dpi = 800)
rm(plot_accept, plot_arret_prec, plot_pts, plot_duree)


# II/ Plot for all correlations ----

plot_accept <- plot_opchar_facet(tableau_results_Rvar, metadata_eng, "accept_ttt")
ggsave(filename = "figures/results/accept_ttt_Rvar_rev2.png",
       plot = plot_accept, device = "png",
       height = 16, width = 20)
ggsave(filename = "figures/results/accept_ttt_Rvar_rev2.eps",
       plot = plot_accept, device = cairo_ps,
       height = 9000, width = 17000, units = "px", dpi = 800)

plot_arret_prec <- plot_opchar_facet(tableau_results_Rvar, metadata_eng, "arret_prema")
ggsave(filename = "figures/results/arret_prema_Rvar_rev2.png",
       plot = plot_arret_prec, device = "png",
       height = 16, width = 20)
ggsave(filename = "figures/results/arret_prema_Rvar_rev2.eps",
       plot = plot_arret_prec, device = cairo_ps,
       height = 9000, width = 17000, units = "px", dpi = 800)

plot_pts <- plot_opchar_facet(tableau_results_Rvar, metadata_eng, "moy_pts")
ggsave(filename = "figures/results/moy_pts_Rvar_rev2.png",
       plot = plot_pts, device = "png",
       height = 16, width = 20)
ggsave(filename = "figures/results/moy_pts_Rvar_rev2.eps",
       plot = plot_pts, device = cairo_ps,
       height = 9000, width = 17000, units = "px", dpi = 800)

plot_duree <- plot_opchar_facet(tableau_results_Rvar, metadata_eng, "moy_duree")
ggsave(filename = "figures/results/moy_duree_Rvar_rev2.png",
       plot = plot_duree, device = "png",
       height = 16, width = 20)
ggsave(filename = "figures/results/moy_duree_Rvar_rev2.eps",
       plot = plot_duree, device = cairo_ps,
       height = 9000, width = 17000, units = "px", dpi = 800)
rm(plot_accept, plot_arret_prec, plot_pts, plot_duree)


# III/ Alpha risk varying with correlation ----

plot_alpha_pos <- tableau_results_Rvar %>% 
  filter(scenar %in% c("H0", "H1")) %>% 
  select(R, probas, scenar, starts_with("accept_ttt")) %>%
  pivot_longer(cols = c(accept_ttt_topiva:accept_ttt_toxprior)) %>%
  mutate(scenar = factor(scenar, levels = c("H0", "H1"), labels = c("Scenario 1", "Scenario 2")),
         name = factor(name,
                       levels = c("accept_ttt_topiva", "accept_ttt_bopiva", "accept_ttt_simon", "accept_ttt_toxrapp", "accept_ttt_toxprior"),
                       labels = c("TOP<sub>eff</sub> + PP<sub>tox</sub>", 
                                  "BOP<sub>eff</sub> + PP<sub>tox</sub>", 
                                  "Simon + PP<sub>tox</sub>",
                                  "TOP<sup>t</sup><sub>eff/tox</sub>", 
                                  "iTOP<sub>eff/tox</sub>"))) %>% 
  filter(!is.na(name)) %>% 
  ggplot(aes(x = R, y = value, color = name)) +
  geom_point() +
  geom_line() +
  geom_hline(data = data.frame(scenar = c("Scenario 1", "Scenario 2"), intercept = c(.05, .9)), aes(yintercept = intercept), linetype = "dashed") +
  facet_wrap(vars(scenar), scales = "free_y") +
  scale_y_continuous(labels = scales::percent_format()) +
  theme(legend.title = element_text(face = "bold"),
        legend.text = element_markdown()) +
  labs(color = "Design",
       x = "Correlation between efficacy and toxicity",
       y = "% of acceptation of the treatment")
# plot_alpha_pos <- tableau_results_Rvar %>% 
#   filter(scenar %in% c("H0", "H1")) %>% 
#   select(R, probas, scenar, starts_with("accept_ttt")) %>%
#   pivot_longer(cols = c(accept_ttt_topiva:accept_ttt_toxprior)) %>%
#   mutate(scenar = factor(scenar, levels = c("H0", "H1"), labels = c("Scenario 1", "Scenario 2")),
#          name = factor(name,
#                        levels = c("accept_ttt_topiva", "accept_ttt_bopiva", "accept_ttt_simon", "accept_ttt_classique", "accept_ttt_toxrapp", "accept_ttt_toxprior"),
#                        labels = c("TOP<sub>eff</sub> + PP<sub>tox</sub>", 
#                                   "BOP<sub>eff</sub> + PP<sub>tox</sub>", 
#                                   "Simon + PP<sub>tox</sub>", 
#                                   "TOP<sub>eff/tox</sub>", 
#                                   "TOP<sup>t</sup><sub>eff/tox</sub>", 
#                                   "iTOP<sub>eff/tox</sub>"))) %>% 
#   ggplot(aes(x = R, y = value, color = name)) +
#   geom_point() +
#   geom_line() +
#   geom_hline(data = data.frame(scenar = c("Scenario 1", "Scenario 2"), intercept = c(.05, .9)), aes(yintercept = intercept), linetype = "dashed") +
#   facet_wrap(vars(scenar), scales = "free_y") +
#   scale_y_continuous(labels = scales::percent_format()) +
#   theme(legend.title = element_text(face = "bold"),
#         legend.text = element_markdown()) +
#   labs(color = "Design",
#        x = "Correlation between efficacy and toxicity",
#        y = "% of acceptation of the treatment")
ggsave(plot = plot_alpha_pos,
       filename = "figures/results/alpha_puiss_pos_rev2.png",
       device = "png", height = 8, width = 12)
ggsave(filename = "figures/results/alpha_puiss_pos_rev2.eps",
       plot = plot_alpha_pos, device = cairo_ps,
       height = 6000, width = 10000, units = "px", dpi = 800)
