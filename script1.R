#### packages ####
library(ismev)
library(evd)
library(fExtremes)
library(stringr)
library(dplyr)

#### data ####
load("Data/donneesStations.RData")
load("Data/donneesVagues.RData")
source("funcs1.R")

str(buoysInfos)
str(donneesVague)

SA <- donneesVague$station6
head(SA)
m <- length(SA)

SB <- donneesVague$station14

############################### AJUSTELEMENT GEV ###############################

#### Extraction des maximas ####

df <- data.frame(date = donneesVague$date, station = SA)

# Extraction de l'année
df$year <- format(as.Date(df$date, format="%Y-%m-%d"),"%Y")
# Extraction du mois
df$month <- format(as.Date(df$date, format = "%Y-%m-%d"), "%m")

# Max par an
maxX1 <- df %>% group_by(year) %>% summarise(Max = max(station)) %>% pull()
length(maxX1)
# Nombre de mesures par an
eff1 <- df %>% group_by(year) %>% count() %>% pull()

# Max par mois
maxX2 <- df %>% group_by(year, month, .add = TRUE) %>%
  summarise(Max = max(station)) %>% pull()
length(maxX2)
# Nombre de mesures par mois
eff2 <- df %>% group_by(year, month, .add = TRUE) %>% count() %>% pull()

# # par bloc
# nb_bloc <- m / 1060 # 438 blocs (maximas)
# n <- m / nb_bloc # bloc de taille 1060
# 
# maxX2 <- extractMax(data =SA, nb_bloc, n, plot_bloc = FALSE,
#                    max.x = "Indice", max.y = "Hauteur",
#                    max.title = "Hauteurs de vagues maximales")

#### Ajustement ####

# ajout de tendence à discuter

maxX1_GEV1 <- gev.fit(maxX1)
maxX1_GEV1$mle
maxX1_GEV1$nllh

# Paramètres
mu1 <- maxX1_GEV1$mle[1]
sigma1 <- maxX1_GEV1$mle[2]
gamma1 <- maxX1_GEV1$mle[3] # gamma < 0 => Weibull

gev.diag(maxX1_GEV1) # Correct

# Intervalle de confiance
IC_GEV(maxX1_GEV1, alpha = 0.05) # IC_gamma contient 0, tester si Gumbel mieux

# Ajustement d'une Gumbel
maxX1_gum <- gum.fit(maxX1)
gum.diag(maxX1_gum) # Retrun level et Quantile plot moins convaincants

# Comparaison des modèles : Test gamma = 0
2 * (maxX1_GEV1$nllh - maxX1_gum$nllh) <= qchisq(0.95, 1)
# On rejette gamma = 0

#### Niveaux de retour ####

# Périodes et probas
T1 <- 100 ; q1 <- 1/T1 ; yq1 <- - log(1 - q1)
T2 <- 500 ; q2 <- 1/T2 ; yq2 <- - log(1 - q2)
T3 <- 1000 ; q3 <- 1/T3 ; yq3 <- - log(1 - q3)

# Quantiles
quant(maxX1_GEV1, 1/100, mu1, sigma1, gamma1, 0.05)
quant(maxX1_GEV1, 1/500, mu1, sigma1, gamma1, 0.05)
quant(maxX1_GEV1, 1/1000, mu1, sigma1, gamma1, 0.05)

############################### AJUSTELEMENT GPD ###############################

# Choix du seuil
mrlplot(SA)  # trop lourd à faire tourner, alternative nécessaire ?
