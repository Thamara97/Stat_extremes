#### packages ####
library(ismev)
library(evd)
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
maxA1 <- df %>% group_by(year) %>% summarise(Max = max(station)) %>% pull()
length(maxA1)
# Nombre de mesures par an
eff1 <- df %>% group_by(year) %>% count() %>% pull()
# Représentation graphique des max
plot(1961:2012, maxA1, xlab = "Années", ylab = "Hauteur mesurée",
     main = "Hauteur maximales des vagues par années")

# Max par mois
maxA2 <- df %>% group_by(year, month, .add = TRUE) %>%
  summarise(Max = max(station)) %>% pull()
length(maxA2)
# Nombre de mesures par mois
eff2 <- df %>% group_by(year, month, .add = TRUE) %>% count() %>% pull()
# Représentation graphique des max
plot(maxA2, xlab = "Mois", ylab = "Hauteur mesurée",
     main = "Hauteur maximales des vagues par mois")

#### Ajustement ####

# Max annuels

maxA1_GEV1 <- gev.fit(maxA1)
maxA1_GEV1$mle
maxA1_GEV1$nllh

# ajout de tendance à discuter
mat1 <- matrix(c(1:52, (1:52)^2), ncol = 2)
head(mat1)
maxA1_GEV0 <- gev.fit(maxA1, ydat = mat1, mul = 1)
maxA1_GEV0$mle
maxA1_GEV0$nllh
D1 <- - 2 * (maxA1_GEV0$nllh - maxA1_GEV1$nllh)
D1 <= qchisq(0.99, df = 1)

# Paramètres
mu1 <- maxA1_GEV1$mle[1]
sigma1 <- maxA1_GEV1$mle[2]
gamma1 <- maxA1_GEV1$mle[3] # gamma < 0 => Weibull

gev.diag(maxA1_GEV1) # Correct

# Intervalle de confiance
IC_GEV(maxA1_GEV1, alpha = 0.05) # IC_gamma contient 0, tester si Gumbel mieux

# Ajustement d'une Gumbel
maxA1_gum <- gum.fit(maxA1)
gum.diag(maxA1_gum) # Retrun level et Quantile plot moins convaincants

# Comparaison des modèles : Test gamma = 0
2 * (maxA1_GEV1$nllh - maxA1_gum$nllh) <= qchisq(0.95, 1)
# On rejette gamma = 0

# Max mensuels

maxA2_GEV1 <- gev.fit(maxA2)
maxA2_GEV1$mle
maxA2_GEV1$nllh

# ajout de tendance à discuter
mat2 <- matrix(c(1:624, (1:624)^2), ncol = 2)
head(mat2)
maxA2_GEV0 <- gev.fit(maxA2, ydat = mat2, mul = 1)
maxA2_GEV0$mle
maxA2_GEV0$nllh
D2 <- - 2 * (maxA2_GEV0$nllh - maxA2_GEV1$nllh)
D2 <= qchisq(0.99, df = 1)

# Paramètres
mu2 <- maxA2_GEV1$mle[1]
sigma2 <- maxA2_GEV1$mle[2]
gamma2 <- maxA2_GEV1$mle[3] # gamma < 0 => Weibull

gev.diag(maxA2_GEV1) # Correct

# Intervalle de confiance
IC_GEV(maxA2_GEV1, alpha = 0.05) # IC_gamma ne contient pas 0, test pas Gumbel

# # Ajustement d'une Gumbel
# maxA2_gum <- gum.fit(maxA2)
# gum.diag(maxA2_gum) # Retrun level et Quantile plot moins convaincants
# 
# # Comparaison des modèles : Test gamma = 0
# 2 * (maxA2_GEV1$nllh - maxA2_gum$nllh) <= qchisq(0.95, 1)
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
quant(maxA1_GEV1, q1, mu1, sigma1, gamma1, alpha)

quant100 <- fgev(maxA1, prob = q1)

quant100$estimate[1]
c(quant100$estimate[1] - 1.96 * quant100$std.err[1],
  quant100$estimate[1] + 1.96 * quant100$std.err[1])

# 500 ans
quant(maxA1_GEV1, q2, mu1, sigma1, gamma1, alpha)

quant500 <- fgev(maxA1, prob = q2)

quant500$estimate[1]
c(quant500$estimate[1] - 1.96 * quant500$std.err[1],
  quant500$estimate[1] + 1.96 * quant500$std.err[1])

# 1000 ans
quant(maxA1_GEV1, q3, mu1, sigma1, gamma1, alpha)

quant1000 <- fgev(maxA1, prob = q3)

quant1000$estimate[1]
c(quant1000$estimate[1] - 1.96 * quant1000$std.err[1],
  quant1000$estimate[1] + 1.96 * quant1000$std.err[1])

# Max mensuels

# Périodes et probas
T1 <- 100 * 12 ; q1 <- 1/T1
T2 <- 500 * 12 ; q2 <- 1/T2
T3 <- 1000 * 12 ; q3 <- 1/T3

alpha <- 0.05

# Quantiles des max mensuels

# 100 ans
quant(maxA2_GEV1, q1, mu2, sigma2, gamma2, alpha)

quant100 <- fgev(maxA2, prob = q1)

quant100$estimate[1]
c(quant100$estimate[1] - 1.96 * quant100$std.err[1],
  quant100$estimate[1] + 1.96 * quant100$std.err[1])

# 500 ans
quant(maxA2_GEV1, q2, mu2, sigma2, gamma2, alpha)

quant500 <- fgev(maxA2, prob = q2)

quant500$estimate[1]
c(quant500$estimate[1] - 1.96 * quant500$std.err[1],
  quant500$estimate[1] + 1.96 * quant500$std.err[1])

# 1000 ans
quant(maxA2_GEV1, q3, mu2, sigma2, gamma2, alpha)

quant1000 <- fgev(maxA2, prob = q3)

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

################################## Station SB ##################################

SB <- donneesVague$station12

plot(SA, SB)

#### Extraction des maximas pour la station SB ####

df2 <- data.frame(date = donneesVague$date, station = SB)

# Extraction de l'année
df2$year <- format(as.Date(df2$date, format="%Y-%m-%d"),"%Y")
# Extraction du mois
df2$month <- format(as.Date(df2$date, format = "%Y-%m-%d"), "%m")

# Max par an
maxB1 <- df2 %>% group_by(year) %>% summarise(Max = max(station)) %>% pull()
length(maxB1)
# Nombre de mesures par an
eff3 <- df2 %>% group_by(year) %>% count() %>% pull()

# Représentation graphique
plot(1961:2012, maxB1, xlab = "Années", ylab = "Hauteur mesurée",
     main = "Hauteur maximales des vagues par années")

# Max par mois
maxB2 <- df2 %>% group_by(year, month, .add = TRUE) %>%
  summarise(Max = max(station)) %>% pull()
length(maxB2)
# Nombre de mesures par mois
eff4 <- df %>% group_by(year, month, .add = TRUE) %>% count() %>% pull()

plot(maxB2, xlab = "Mois", ylab = "Hauteur mesurée",
     main = "Hauteur maximales des vagues par mois")

plot(maxA1, maxB1)
plot(maxA2, maxB2)

#### Ajustement ####

# Max annuels

maxB1_GEV1 <- gev.fit(maxB1)
maxB1_GEV1$mle
maxB1_GEV1$nllh

# ajout de tendance à discuter
maxB1_GEV0 <- gev.fit(maxB1, ydat = mat1, mul = 1:2)
maxB1_GEV0$mle
maxB1_GEV0$nllh
D3 <- - 2 * (maxB1_GEV0$nllh - maxB1_GEV1$nllh)
D3 <= qchisq(0.99, df = 1)
# Rejet de la tendance

# Paramètres
mu3 <- maxB1_GEV1$mle[1]
sigma3 <- maxB1_GEV1$mle[2]
gamma3 <- maxB1_GEV1$mle[3] # gamma > 0 => Fréchet

gev.diag(maxB1_GEV1) # Correct

# Intervalle de confiance
IC_GEV(maxB1_GEV1, alpha = 0.05) # IC_gamma contient 0, tester si Gumbel mieux

# Ajustement d'une Gumbel
maxB1_gum <- gum.fit(maxB1)
gum.diag(maxB1_gum) # Retrun level et Quantile plot moins convaincants

# Comparaison des modèles : Test gamma = 0
- 2 * (maxB1_GEV1$nllh - maxB1_gum$nllh) <= qchisq(0.95, 1)
# On rejette gamma = 0

# Max mensuels

maxB2_GEV1 <- gev.fit(maxB2)
maxB2_GEV1$mle
maxB2_GEV1$nllh

# ajout de tendance à discuter
maxB2_GEV0 <- gev.fit(maxB2, ydat = mat2, mul = 1)
maxB2_GEV0$mle
maxB2_GEV0$nllh
D4 <- - 2 * (maxB2_GEV0$nllh - maxB2_GEV1$nllh)
D4 <= qchisq(0.95, df = 1)
# Rejet de la tendance

# Paramètres
mu4 <- maxB2_GEV1$mle[1]
sigma4 <- maxB2_GEV1$mle[2]
gamma4 <- maxB2_GEV1$mle[3] # gamma < 0 => Weibull

gev.diag(maxB2_GEV1) # Correct

# Intervalle de confiance
IC_GEV(maxB2_GEV1, alpha = 0.05) # IC_gamma ne contient pas 0, test pas Gumbel

# # Ajustement d'une Gumbel
maxB2_gum <- gum.fit(maxB2)
gum.diag(maxB2_gum) # Retrun level et Quantile plot moins convaincants

# Comparaison des modèles : Test gamma = 0
- 2 * (maxB2_GEV1$nllh - maxB2_gum$nllh) <= qchisq(0.95, 1)
# On rejette gamma = 0

#### Niveaux de retour ####

# Max annuels

# Quantiles des max annuels

# 100 ans
quant(maxB1_GEV1, q1, mu3, sigma3, gamma3, alpha)

quant100_B <- fgev(maxB1, prob = q1)

quant100_B$estimate[1]
c(quant100_B$estimate[1] - 1.96 * quant100_B$std.err[1],
  quant100_B$estimate[1] + 1.96 * quant100_B$std.err[1])

# 500 ans
quant(maxB1_GEV1, q2, mu3, sigma3, gamma3, alpha)

quant500_B <- fgev(maxB1, prob = q2)

quant500_B$estimate[1]
c(quant500_B$estimate[1] - 1.96 * quant500_B$std.err[1],
  quant500_B$estimate[1] + 1.96 * quant500_B$std.err[1])

# 1000 ans
quant(maxB1_GEV1, q3, mu3, sigma3, gamma3, alpha)

quant1000_B <- fgev(maxB1, prob = q3)

quant1000_B$estimate[1]
c(quant1000_B$estimate[1] - 1.96 * quant1000_B$std.err[1],
  quant1000_B$estimate[1] + 1.96 * quant1000_B$std.err[1])

# Max mensuels

# Quantiles des max mensuels

# 100 ans
quant(maxB2_GEV1, q1, mu4, sigma4, gamma4, alpha)

quant100_B <- fgev(maxB2, prob = q1)

quant100_B$estimate[1]
c(quant100_B$estimate[1] - 1.96 * quant100_B$std.err[1],
  quant100_B$estimate[1] + 1.96 * quant100_B$std.err[1])

# 500 ans
quant(maxB2_GEV1, q2, mu4, sigma4, gamma4, alpha)

quant500_B <- fgev(maxB2, prob = q2)

quant500_B$estimate[1]
c(quant500_B$estimate[1] - 1.96 * quant500_B$std.err[1],
  quant500_B$estimate[1] + 1.96 * quant500_B$std.err[1])

# 1000 ans
quant(maxB2_GEV1, q3, mu4, sigma4, gamma4, alpha)

quant1000_B <- fgev(maxB2, prob = q3)

quant1000_B$estimate[1]
c(quant1000_B$estimate[1] - 1.96 * quant1000_B$std.err[1],
  quant1000_B$estimate[1] + 1.96 * quant1000_B$std.err[1])

# Les quantiles des max mensuels sont plus élevés que ceux des max annuels


chiplot(cbind(SA, SB))
# plot 1 chi : semble tendre vers 0 (indique une indépendance asymptotique)
# plot 2 chi_bar : plus petit que 1 (confirme l'indépendance asymptotique)

#### Tranformation des données en Fréchet 1 ####
SA_F <- qgev(pgev(SA,loc = 3.9868599, scale = 1.1734179 ,
                            shape = -0.1526229),
                  loc = 1, shape = 1, scale = 1)
SB_F <- qgev(pgev(SB, loc = 1.43213305, scale = 0.71594773,
                            shape = -0.02667831),
                  loc = 1, shape = 1, scale = 1)
plot(SA_F, SB_F)
# On voit très nettement un indépendance asymptotique.

# Etude plus poussée de la dépendance avec A
abvnonpar(data = cbind(maxA2, maxB2), plot = TRUE, col = "red")

abvnonpar(data = cbind(maxA2, maxB2), plot = FALSE) * 2
# theta est proche de 2 donc on à bien de l'indépendance asymptotique.
# Indépendance assez logique puisque stations distantes d'environ 264km

# Fit d'un modèle
fbvevd(cbind(maxA2, maxB2), model = "alog")
MGEV1 <- fbvevd(cbind(maxA2, maxB2), model = "log") # meilleur AIC
MGEV1

#### Recherche de z_p tq P(X > z_p | Y > 8) = p
# X=SB et Y=SA

# On suppose qu'on a l'indépendance asymptotique, donc p ~ G^-1_X((1-p)^n)
quant_cond(MGEV1, 8, 10^-4, 624, dep = FALSE, cond.mar2 = FALSE)

