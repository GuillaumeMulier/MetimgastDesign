# --------------------------------------------------------------- #
# Script pour faire la figure du déroulement de l'essai METIMGAST #
# Créé le 01/10/2021, modifié le 05/10/2021                       #
# --------------------------------------------------------------- #

# Packages et helpers ----

library(tidyverse)

theme_set(theme_light())

# Fabriquer les données ----

donnee_deroulement <- tibble(
  patients = 1:90,
  screening = seq(0, 548, 6.1), # Date de screening des patients
  inclusion = screening + 14, # Date max d'inclusion
  tox_C2 = inclusion + 42, # Toxicité dans les 2 1ers cycles
  eff_M6 = inclusion + 183, # Efficacité à 6 mois
  fin_ttt = inclusion + 365.25, # Fin du traitement
  fin_suivi = inclusion + 730.5 # Fin du suivi
)

# GGPLOT ----

fig_essai <- ggplot(donnee_deroulement, aes(y = patients)) +
  geom_ribbon(aes(xmin = screening, xmax = inclusion, fill = "Screening")) +
  geom_ribbon(aes(xmin = inclusion, xmax = fin_ttt, fill = "Traitement pendant 1 an max")) +
  geom_ribbon(aes(xmin = fin_ttt, xmax = fin_suivi, fill = "Suivi jusqu'à 2 ans max")) +
  geom_polygon(data = data.frame(xx = c(14, 556.9, 739.9, 197), yy = c(1, 90, 90, 1)), aes(xx, yy), color = "blue", size = 2.5, fill = "transparent") +
  geom_polygon(data = data.frame(xx = c(14, 556.9, 598.9, 56), yy = c(1, 90, 90, 1)), aes(xx, yy), color = "red", size = 1.2, fill = "transparent") +
  annotate("text", x = 50, y = -5, label = "Toxicité au cours des 2 premiers cycles", color = "red", hjust = 0) +
  annotate("text", x = 50, y = -2, label = "Efficacité au cours des 8 premiers cycles", color = "blue", hjust = 0) +
  scale_y_reverse(breaks = c(1, 15, 30, 45, 60, 75, 90)) +
  scale_x_continuous(breaks = seq(0, 1278.375, 182.625), labels = function(x) paste0(x / 365.25, " an", ifelse(x / 365.25 > 1, "", "s"))) +
  scale_fill_manual(values = c("#beaed4", "#7fc97f", "#fdc086")) +
  labs(x = "Temps depuis le scrrening du premier patient",
       y = "Index du patient",
       fill = "Etape de l'essai")
  
ggsave(filename = "figures/brouillon/real_essai.png",
       plot = fig_essai, device = "png",
       width = 10, height = 7)  
  
