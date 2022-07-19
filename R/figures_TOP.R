# ----------------------------------------- #
# Script de figures pour comprendre le TOP  #
# Créé le 17/08/2021, modifié le 17/08/2021 #
# ----------------------------------------- #

library(tidyverse)
theme_set(theme_light())

# Influence du TESS
ggplot() +
  map(seq(5, 25, 2), function(x) geom_function(aes(color = paste0("TESS: ", x)), fun = ~ pbeta(.x, .3 + 3, .7 + x - 3))) +
  xlim(c(0, 1)) +
  geom_vline(xintercept = .3, linetype = "dashed", color = "red") +
  scale_color_discrete(name = "TESS", breaks = paste0("TESS: ", seq(5, 25, 2))) +
  labs(title = "Efficacité")

ggplot() +
  map(seq(5, 25, 2), function(x) geom_function(aes(color = paste0("TESS: ", x)), fun = ~ 1 - pbeta(.x, .3 + 3, .7 + x - 3))) +
  xlim(c(0, 1)) +
  geom_vline(xintercept = .3, linetype = "dashed", color = "red") +
  scale_color_discrete(name = "TESS", breaks = paste0("TESS: ", seq(5, 25, 2))) +
  labs(title = "Toxicité")

# Influence du nombre d'efficacité
ggplot() +
  map(seq(0, 10), function(x) geom_function(aes(color = paste0(x, " eff")), fun = ~ pbeta(.x, .3 + x, .7 + 10 - x))) +
  xlim(c(0, 1)) +
  geom_vline(xintercept = .3, linetype = "dashed", color = "red") +
  scale_color_discrete(name = "# d'efficacités", breaks = paste0(0:10, " eff")) +
  labs(title = "Efficacité")

# Influence du nombre de toxicités
ggplot() +
  map(seq(0, 10), function(x) geom_function(aes(color = paste0(x, " tox")), fun = ~ 1 - pbeta(.x, .3 + x, .7 + 10 - x))) +
  xlim(c(0, 1)) +
  geom_vline(xintercept = .3, linetype = "dashed", color = "red") +
  scale_color_discrete(name = "# d'efficacités", breaks = paste0(0:10, " tox")) +
  labs(title = "Toxicité")
