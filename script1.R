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

# Représentation graphique
plot(maxX1, xlab = "Années", ylab = "Hauteur mesurée",
     main = "Hauteur maximales des vagues par années")

# Max par mois
maxX2 <- df %>% group_by(year, month, .add = TRUE) %>%
  summarise(Max = max(station)) %>% pull()
length(maxX2)
# Nombre de mesures par mois
eff2 <- df %>% group_by(year, month, .add = TRUE) %>% count() %>% pull()

plot(maxX2, xlab = "Mois", ylab = "Hauteur mesurée",
     main = "Hauteur maximales des vagues par mois")

#### Ajustement ####

# ajout de tendance à discuter

# Max annuels

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

# Max mensuels

maxX2_GEV1 <- gev.fit(maxX2)
maxX2_GEV1$mle
maxX2_GEV1$nllh

# Paramètres
mu2 <- maxX2_GEV1$mle[1]
sigma2 <- maxX2_GEV1$mle[2]
gamma2 <- maxX2_GEV1$mle[3] # gamma < 0 => Weibull

gev.diag(maxX2_GEV1) # Correct

# Intervalle de confiance
IC_GEV(maxX2_GEV1, alpha = 0.05) # IC_gamma ne contient pas 0, test pas Gumbel

# # Ajustement d'une Gumbel
# maxX2_gum <- gum.fit(maxX2)
# gum.diag(maxX2_gum) # Retrun level et Quantile plot moins convaincants
# 
# # Comparaison des modèles : Test gamma = 0
# 2 * (maxX2_GEV1$nllh - maxX2_gum$nllh) <= qchisq(0.95, 1)
# # On rejette gamma = 0

#### Niveaux de retour ####

# Max annuels

# Périodes et probas
T1 <- 100 ; q1 <- 1/T1
T2 <- 500 ; q2 <- 1/T2
T3 <- 1000 ; q3 <- 1/T3

alpha <- 0.05

# Quantiles des max annuels

# 100 ans
quant(maxX1_GEV1, q1, mu1, sigma1, gamma1, alpha)

quant100 <- fgev(maxX1, prob = q1)

quant100$estimate[1]
c(quant100$estimate[1] - 1.96 * quant100$std.err[1],
  quant100$estimate[1] + 1.96 * quant100$std.err[1])

# 500 ans
quant(maxX1_GEV1, q2, mu1, sigma1, gamma1, alpha)

quant500 <- fgev(maxX1, prob = q2)

quant500$estimate[1]
c(quant500$estimate[1] - 1.96 * quant500$std.err[1],
  quant500$estimate[1] + 1.96 * quant500$std.err[1])

# 1000 ans
quant(maxX1_GEV1, q3, mu1, sigma1, gamma1, alpha)

quant1000 <- fgev(maxX1, prob = q3)

quant1000$estimate[1]
c(quant1000$estimate[1] - 1.96 * quant1000$std.err[1],
  quant1000$estimate[1] + 1.96 * quant1000$std.err[1])

# Max mensuels

# Périodes et probas
T1 <- 100 * 12 ; q1 <- 1/T1
T2 <- 500 * 12 ; q2 <- 1/T2
T3 <- 1000 * 12 ; q3 <- 1/T3

alpha <- 0.05

# Quantiles des max annuels

# 100 ans
quant(maxX2_GEV1, q1, mu2, sigma2, gamma2, alpha)

quant100 <- fgev(maxX2, prob = q1)

quant100$estimate[1]
c(quant100$estimate[1] - 1.96 * quant100$std.err[1],
  quant100$estimate[1] + 1.96 * quant100$std.err[1])

# 500 ans
quant(maxX2_GEV1, q2, mu2, sigma2, gamma2, alpha)

quant500 <- fgev(maxX2, prob = q2)

quant500$estimate[1]
c(quant500$estimate[1] - 1.96 * quant500$std.err[1],
  quant500$estimate[1] + 1.96 * quant500$std.err[1])

# 1000 ans
quant(maxX2_GEV1, q3, mu2, sigma2, gamma2, alpha)

quant1000 <- fgev(maxX2, prob = q3)

quant1000$estimate[1]
c(quant1000$estimate[1] - 1.96 * quant1000$std.err[1],
  quant1000$estimate[1] + 1.96 * quant1000$std.err[1])

# Les quantiles des max mensuels sont plus élevés que ceux des max annuels

############################### AJUSTELEMENT GPD ###############################

# Choix du seuil
# mrlplot(SA)  # trop lourd à faire tourner, alternative nécessaire ?

# on choisit le nombre r d'excès souhaité
r1 <- 52 # autant que d'années
r2 <- 12 * r1 # autant que de mois

# seuils correspondant
u1 <- SA[order(SA, decreasing = TRUE)[1:(r1 + 1)]][(r1 + 1)]
u2 <- SA[order(SA, decreasing = TRUE)[1:(r2 + 1)]][(r2 + 1)]

# clusters à envisager puisqu'il y a dépendance entre les données
GPD1 <- fpot(SA, u1, std.err = TRUE, cmax = TRUE, r = 5)
GPD1
plot(GPD1)

GPD2 <- fpot(SA, u2, cmax = TRUE, r = 5)
GPD2
plot(GPD2)

# Diagnostiques graphique pas terrible