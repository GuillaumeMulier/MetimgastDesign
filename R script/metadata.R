# -------------------------- #
# Metadata for figures, etc. #
# -------------------------- #

library(tidyverse)
source("R script/functions.R")

theme_set(theme_light() +
            theme(strip.background = element_rect(fill = NA),
                  strip.text = element_textbox(
                    size = 12, 
                    color = "white", fill = "#7888C0", box.color = "#000066",
                    halign = 0.5, linetype = 1, r = unit(3, "pt"), width = unit(0.75, "npc"),
                    padding = margin(2, 0, 1, 0), margin = margin(3, 3, 3, 3))))

# Scenarios and probabilities ----

# Get all the scenarios with all correlations
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
liste_scenars <- map2_dfr(
  liste_scenars,
  noms_scenars,
  ~ generer_6corr_hyp(.x["Peff"], .x["Ptox"], .y)
)
liste_scenars$scenar_p <- factor(liste_scenars$scenar,
                                 levels = rev(c("H0", "H1", "Int1", "Int2", "Paseff1", "Paseff2", "Effint1", "Effint2", "Tox1", "Tox2")),
                                 labels = rev(c("Sc 1 : (0.15;0.30) [H<sub>0</sub>]", "Sc 2 : (0.30;0.20) [H<sub>1</sub>]", "Sc 3: (0.20;0.25)", "Sc 4 : (0.25;0.25)", "Sc 5 : (0.15;0.20)",
                                                "Sc 6 : (0.10;0.15)", "Sc 7 : (0.20;0.20)", "Sc 8 : (0.20;0.30)", "Sc 9 : (0.30;0.30)", "Sc 10 : (0.40;0.35)")))
liste_scenars$scenarios <- factor(liste_scenars$scenar,
                                  levels = rev(c("H0", "H1", "Int1", "Int2", "Paseff1", "Paseff2", "Effint1", "Effint2", "Tox1", "Tox2")),
                                  labels = rev(c("Sc 1 :", "Sc 2 :", "Sc 3 :", "Sc 4 :", "Sc 5 :",
                                                 "Sc 6 :", "Sc 7 :", "Sc 8 :", "Sc 9 :", "Sc 10 :")))
liste_scenars$probas_plot <- factor(liste_scenars$scenar,
                                    levels = rev(c("H0", "H1", "Int1", "Int2", "Paseff1", "Paseff2", "Effint1", "Effint2", "Tox1", "Tox2")),
                                    labels = rev(c("(0.15;0.30) [H<sub>0</sub>]", "(0.30;0.20) [H<sub>1</sub>]", "(0.20;0.25)", "(0.25;0.25)", "(0.15;0.20)",
                                                   "(0.10;0.15)", "(0.20;0.20)", "(0.20;0.30)", "(0.30;0.30)", "(0.40;0.35)")))
liste_scenars$nom_corr <- factor(liste_scenars$nom_corr,
                                 levels = c("Rmin", "Rneg", "Rind", "Rpos1", "Rpos2", "Rmax"),
                                 labels = c("R<sub>min</sub>", "R<sub>neg</sub>", "R<sub>ind</sub>", "R<sub>pos,1</sub>", "R<sub>pos,2</sub>", "R<sub>max</sub>"))

# Interim analyses for toxicity ----

vec_eval_tox <- c(5, 10, 15, 20, 30, 40, 50, 60, 70, 80, 90)

# Graphical parameters ----

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

