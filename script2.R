#### packages ####
library(ismev)
library(evd)
library(stringr)
library(dplyr)
library(ggplot2)

#### data ####
load("Data/donneesStations.RData")
load("Data/donneesVagues.RData")
source("funcs1.R")

str(buoysInfos)
str(donneesVague)

SA <- donneesVague$station6
head(SA)
m <- length(SA)

SB <- donneesVague$station12

############################### AJUSTELEMENT GEV ###############################

#### Extraction des maximas ####

# Extraction des données et des dates
df <- data.frame(date = donneesVague$date, station = SA)

# Extraction de l'année
df$year <- format(as.Date(df$date, format="%Y-%m-%d"),"%Y")
# Extraction du mois
df$month <- format(as.Date(df$date, format = "%Y-%m-%d"), "%m")

# Max par an
maxA1 <- df %>% group_by(year) %>% summarise(Max = max(station)) %>% pull()
n1 <- length(maxA1)
# Nombre de mesures par an
eff1 <- df %>% group_by(year) %>% count() %>% pull()
# Représentation graphique des max
plot(1961:2012, maxA1, xlab = "Années", ylab = "Hauteur mesurée",
     main = "Hauteur maximales des vagues par années à la station 6 (SA)")

# Max par mois
maxA2 <- df %>% group_by(year, month, .add = TRUE) %>%
  summarise(Max = max(station)) %>% pull()
n2 <- length(maxA2)
# Nombre de mesures par mois
eff2 <- df %>% group_by(year, month, .add = TRUE) %>% count() %>% pull()
# Représentation graphique des max
plot(maxA2, xlab = "Mois", ylab = "Hauteur mesurée",
     main = "Hauteur maximales des vagues par mois")

#### Ajustement ####

# Max annuels

maxA1_GEV1 <- gev.fit(maxA1)

# ajout de tendance à discuter
mat1 <- matrix(c(1:52, (1:52)^2), ncol = 2)
maxA1_GEV0 <- gev.fit(maxA1, ydat = mat1, mul = 1)
IC_GEV(maxA1_GEV0, alpha = 0.05)
D1 <- - 2 * (maxA1_GEV0$nllh - maxA1_GEV1$nllh)
D1 <= qchisq(0.95, df = 1) # on ne rejette pas que la tendance est nulle

# Paramètres
mu1 <- maxA1_GEV1$mle[1]
sigma1 <- maxA1_GEV1$mle[2]
gamma1 <- maxA1_GEV1$mle[3] # gamma < 0 => Weibull

gev.diag(maxA1_GEV1) # Correct

# Intervalle de confiance
round(IC_GEV(maxA1_GEV1, alpha = 0.05), 3) # IC_gamma contient 0, tester si Gumbel mieux

# Ajustement d'une Gumbel
maxA1_gum <- gum.fit(maxA1)
round(IC_GEV(maxA1_gum, alpha = 0.95, gum = TRUE), 3)
gum.diag(maxA1_gum) # Retrun level et Quantile plot moins convaincants

# Comparaison des modèles : Test gamma = 0
- 2 * (maxA1_GEV1$nllh - maxA1_gum$nllh) <= qchisq(0.95, 1)
# On rejette gamma = 0

# Max mensuels

maxA2_GEV1 <- gev.fit(maxA2)

# ajout de tendance à discuter
mat2 <- matrix(c(1:624, (1:624)^2), ncol = 2)
maxA2_GEV0 <- gev.fit(maxA2, ydat = mat2, mul = 1)
D2 <- - 2 * (maxA2_GEV0$nllh - maxA2_GEV1$nllh)
D2 <= qchisq(0.99, df = 1)

# Paramètres
mu2 <- maxA2_GEV1$mle[1]
sigma2 <- maxA2_GEV1$mle[2]
gamma2 <- maxA2_GEV1$mle[3] # gamma < 0 => Weibull

gev.diag(maxA2_GEV1) # Correct

# Intervalle de confiance
IC_GEV(maxA2_GEV1, alpha = 0.05) # IC_gamma ne contient pas 0, test pas Gumbel

#### Niveaux de retour ####

# Max annuels

# Périodes et probas
T1 <- 100 ; q1 <- 1/T1
T2 <- 500 ; q2 <- 1/T2
T3 <- 1000 ; q3 <- 1/T3

alpha <- 0.05

# Quantiles des max annuels

# 100 ans
quant(maxA1_GEV1, q1, alpha)

quant100 <- fgev(maxA1, prob = q1)

quant100$estimate[1]
c(quant100$estimate[1] - 1.96 * quant100$std.err[1],
  quant100$estimate[1] + 1.96 * quant100$std.err[1])

# 500 ans
quant(maxA1_GEV1, q2, alpha)

quant500 <- fgev(maxA1, prob = q2)

quant500$estimate[1]
c(quant500$estimate[1] - 1.96 * quant500$std.err[1],
  quant500$estimate[1] + 1.96 * quant500$std.err[1])

# 1000 ans
quant(maxA1_GEV1, q3, alpha)

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

################################## Station SB ##################################

# plot(SA, SB) # long a produire
plot(SA[1:5000], SB[1:5000])

#### Extraction des maximas pour la station SB ####

df2 <- data.frame(date = donneesVague$date, station = SB)

# Extraction de l'année
df2$year <- format(as.Date(df2$date, format="%Y-%m-%d"),"%Y")
# Extraction du mois
df2$month <- format(as.Date(df2$date, format = "%Y-%m-%d"), "%m")

# Max par an
maxB1 <- df2 %>% group_by(year) %>% summarise(Max = max(station)) %>% pull()
length(maxB1)

# Représentation graphique
plot(1961:2012, maxB1, xlab = "Années", ylab = "Hauteur mesurée",
     main = "Hauteur maximales des vagues par années")

# Max par mois
maxB2 <- df2 %>% group_by(year, month, .add = TRUE) %>%
  summarise(Max = max(station)) %>% pull()
length(maxB2)

# Représentation graphique
plot(maxB2, xlab = "Mois", ylab = "Hauteur mesurée",
     main = "Hauteur maximales des vagues par mois")

plot(maxA1, maxB1, xlab = "SA", ylab = "SB",
     main = "Hauteurs maximales des vagues mesurée aux stations SA et SB")
plot(maxA2, maxB2)

#### Ajustement ####

# Max annuels

maxB1_GEV1 <- gev.fit(maxB1)

# ajout de tendance à discuter
maxB1_GEV0 <- gev.fit(maxB1, ydat = mat1, mul = 1)
D3 <- - 2 * (maxB1_GEV0$nllh - maxB1_GEV1$nllh)
D3 <= qchisq(0.95, df = 1)
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

# ajout de tendance à discuter
maxB2_GEV0 <- gev.fit(maxB2, ydat = mat2, mul = 1)
D4 <- - 2 * (maxB2_GEV0$nllh - maxB2_GEV1$nllh)
D4 <= qchisq(0.95, df = 1)
# Rejet de la tendance

# Paramètres
mu4 <- maxB2_GEV1$mle[1]
sigma4 <- maxB2_GEV1$mle[2]
gamma4 <- maxB2_GEV1$mle[3] # gamma < 0 => Weibull

gev.diag(maxB2_GEV1) # Correct

# Intervalle de confiance
IC_GEV(maxB1_GEV1, alpha = 0.05) # IC_gamma contient 0, test Gumbel

# # Ajustement d'une Gumbel
maxB2_gum <- gum.fit(maxB2)
gum.diag(maxB2_gum) # Retrun level et Quantile plot moins convaincants

# Comparaison des modèles : Test gamma = 0
- 2 * (maxB1_GEV1$nllh - maxB2_gum$nllh) <= qchisq(0.95, 1)
# On rejette gamma = 0

#### Etude de la dépendance ####

par(mfrow = c(1,2))
chiplot(cbind(SA, SB), xlab = "u", ylab1 = expression(chi(u)),
        ylab2 = expression(bar(chi)(u)))
# plot 1 chi : semble tendre vers 0 (indique une indépendance asymptotique)
# plot 2 chi_bar : plus petit que 1 (confirme l'indépendance asymptotique)

#### Transformation des données en Fréchet 1 ####
maxA1_F <- qgev(pgev(maxA1,loc = mu1, scale = sigma1,
                     shape = gamma1),
                loc = 1, shape = 1, scale = 1)
maxB1_F <- qgev(pgev(maxB1, loc = mu3, scale = sigma3,
                     shape = gamma3),
                loc = 1, shape = 1, scale = 1)

# Représentation graphique
par(mfrow = c(1,1))
plot(maxA1_F, maxB1_F)
# On voit plus ou moins un indépendance asymptotique.

# Etude plus poussée de la dépendance avec A
abvnonpar(data = cbind(maxA1, maxB1), plot = TRUE, col = "red",
          main = "Estimateur de la fonction de dépendance de Pickands")

abvnonpar(data = cbind(maxA1, maxB1), plot = FALSE) * 2
# theta est proche de 2 donc on à bien de l'indépendance asymptotique.
# Indépendance assez logique puisque stations distantes d'environ 264km

# Fit d'un modèle
fbvevd(cbind(maxA1, maxB1), model = "alog", std.err = FALSE) # AIC = 221.7819
fbvevd(cbind(maxA1, maxB1), model = "log") # AIC = 217.7394
MGEV1 <- fbvevd(cbind(maxA1, maxB1), model = "hr") # AIC = 217.724

#### Quantile conditionnel ####

quant_cond(MGEV1, 4, 10^(-5), n1, dep = FALSE, cond.mar2 = FALSE)

# On suppose qu'on a l'indépendance asymptotique, donc p ~ G^-1_X((1-p)^n)
p <- 10^-seq(1,10, by = 0.1)
z_p <- c()
k <- 0
for (j in p) {
  k <- k + 1
  z_p[k] <- quant_cond(MGEV1, 4, j, n1, dep = FALSE, cond.mar2 = FALSE)
}

plot(seq(1,10, by = 0.1), z_p, type = "l",
     xlab = expression(10^-x),
     ylab = "Quantile",
     main = "Quantile conditionnel en cas d'indépendance asymptotique")

################################## Station SC ##################################

# SC (station 1) plus distante de SA que SB

SC <- donneesVague$station1

# Etude de la dépendance

plot(SA[1:5000], SC[1:5000])
# On voit une dépendance forte

#### Extraction des maximas pour la station SC ####

df3 <- data.frame(date = donneesVague$date, station = SC)

# Extraction de l'année
df3$year <- format(as.Date(df3$date, format="%Y-%m-%d"),"%Y")
# Extraction du mois
df3$month <- format(as.Date(df3$date, format = "%Y-%m-%d"), "%m")

# Max par an
maxC1 <- df3 %>% group_by(year) %>% summarise(Max = max(station)) %>% pull()

# Représentation graphique
plot(1961:2012, maxC1, xlab = "Années", ylab = "Hauteur mesurée",
     main = "Hauteur maximales des vagues par années")

# Max par mois
maxC2 <- df3 %>% group_by(year, month, .add = TRUE) %>%
  summarise(Max = max(station)) %>% pull()

# Représentation graphique
plot(maxC2, xlab = "Mois", ylab = "Hauteur mesurée",
     main = "Hauteur maximales des vagues par mois")

plot(maxA1, maxC1, xlab = "SA", ylab = "SC",
     main = "Hauteurs maximales des vagues mesurée aux stations SA et SC")
plot(maxA2, maxC2)

#### Etude de la dépendance ####

chiplot(cbind(SA, SC))
# plot 1 chi = 1 (indique une dépendance asymptotique)
# plot 2 chi_bar = 1 (confirme la dépendance asymptotique)

#### Transformation des données en Fréchet 1 ####
maxC1_GEV1 <- gev.fit(maxC1)
mu5 <- maxC1_GEV1$mle[1]
sigma5 <- maxC1_GEV1$mle[2]
gamma5 <- maxC1_GEV1$mle[3]

IC_GEV(maxC1_GEV1, alpha = 0.05) # 0 n'est pas de l'IC de gamma

# gumbel ? Possible mais non (à cause des graphes)
maxC1_gum <- gum.fit(maxC1)
- 2 * (maxC1_GEV1$nllh - maxC1_gum$nllh) <= qchisq(0.95, 1)

gum.diag(maxC1_gum)
gev.diag(maxC1_GEV1)

# tendance ? Non
maxC1_GEV0 <- gev.fit(maxC1, ydat = mat1, mul = 1:2)
-2*(maxC1_GEV0$nllh - maxC1_GEV1$nllh) <= qchisq(p = 0.99, df = 1)

maxC1_F <- qgev(pgev(maxC1, loc = mu5, scale = sigma5,
                     shape = gamma5),
                loc = 1, shape = 1, scale = 1)

# Représentation graphique
plot(maxA1_F, maxC1_F)
# On voit une dépendance claire

# Etude plus poussée de la dépendance avec A
abvnonpar(data = cbind(maxA1, maxC1), plot = TRUE, col = "red")

abvnonpar(data = cbind(maxA1, maxC1), plot = FALSE) * 2
# theta est proche de 1 donc on à bien de la dépendance asymptotique.

# Fit d'un modèle
fbvevd(cbind(maxA1, maxC1), model = "alog", std.err = FALSE) # AIC = 31.54714
fbvevd(cbind(maxA1, maxC1), model = "log") # AIC = 34.66942
MGEV2 <- fbvevd(cbind(maxA1, maxC1), model = "hr") # AIC = 3.649479
MGEV2

#### Quantile conditionnel ####

quant_cond(MGEV2, 8, 10^-10, n1, dep = TRUE, cond.mar2 = FALSE)

y <- 4
p <- 10^-seq(1,10, by = 0.1)
l <- 0
for (i in y){
  l <- l + 1
  z_p <- c()
  k <- 0
  for (j in p) {
    k <- k + 1
    z_p[k] <- quant_cond(MGEV2, i, 10^(-j), n1, dep = TRUE, cond.mar2 = FALSE)
  }
  plot(seq(1,10, by = 0.1), z_p, xlab = expression(10^-x), ylab = "Quantile",
       main = paste("Quantile conditionnel en quand Y >",y[l]), type = "l")
}

################################## Station SD ################################## 
# SD (station 4) plus proche de SA que SB

SD <- donneesVague$station4

# Etude de la dépendance

plot(SA[1:5000], SD[1:5000])

#### Extraction des maximas pour la station SD ####

df4 <- data.frame(date = donneesVague$date, station = SD)

# Extraction de l'année
df4$year <- format(as.Date(df4$date, format="%Y-%m-%d"),"%Y")
# Extraction du mois
df4$month <- format(as.Date(df4$date, format = "%Y-%m-%d"), "%m")

# Max par an
maxD1 <- df4 %>% group_by(year) %>% summarise(Max = max(station)) %>% pull()

# Représentation graphique
plot(1961:2012, maxD1, xlab = "Années", ylab = "Hauteur mesurée",
     main = "Hauteur maximales des vagues par années")

# Max par mois
maxD2 <- df4 %>% group_by(year, month, .add = TRUE) %>%
  summarise(Max = max(station)) %>% pull()

# Représentation graphique
plot(maxD2, xlab = "Mois", ylab = "Hauteur mesurée",
     main = "Hauteur maximales des vagues par mois")

plot(maxA1, maxD1)
plot(maxA2, maxD2)

#### Etude de la dépendance ####

chiplot(cbind(SA, SD))
# plot 1 chi semble tendre vers 0 (indique une indépendance asymptotique)
# plot 2 chi_bar < 1 (confirme l'indépendance asymptotique)

#### Transformation des données en Fréchet 1 ####
maxD1_GEV1 <- gev.fit(maxD1)
mu6 <- maxD1_GEV1$mle[1]
sigma6 <- maxD1_GEV1$mle[2]
gamma6 <- maxD1_GEV1$mle[3]

IC_GEV(maxD1_GEV1, alpha = 0.05) # 0 est dans l'IC de gamma

# gumbel ? Non
maxD1_gum <- gum.fit(maxD1)
-2 * (maxD1_GEV1$nllh - maxD1_gum$nllh)
# rejet de gamma = 0


# tendance ? Non
maxD1_GEV0 <- gev.fit(maxD1, ydat = mat1, mul = 1)
-2*(maxD1_GEV0$nllh - maxD1_GEV1$nllh) <= qchisq(p = 0.99, df = 1)

maxD1_F <- qgev(pgev(maxD1, loc = mu6, scale = sigma6,
                     shape = gamma6),
                loc = 1, shape = 1, scale = 1)

# Représentation graphique
plot(maxA1_F, maxD1_F)
# On voit pas grand chose

# Etude plus poussée de la dépendance avec A
abvnonpar(data = cbind(maxA1, maxD1), plot = TRUE, col = "red")

abvnonpar(data = cbind(maxA1, maxD1), plot = FALSE) * 2
# theta est proche de 2 donc par sur de l'indépendance asymptotique.

# Fit d'un modèle
fbvevd(cbind(maxA1, maxD1), model = "alog", std.err = FALSE) # AIC = 219.3239
fbvevd(cbind(maxA1, maxD1), model = "log") # AIC = 215.073
MGEV3 <- fbvevd(cbind(maxA1, maxD1), model = "hr") # AIC = 214.4105
MGEV3

#### Quantile conditionnel ####

quant_cond(MGEV3, 4, 10^-10, n1, dep = FALSE, cond.mar2 = FALSE)

p <- 10^-seq(1,10, by = 0.1)
z_p <- c()
k <- 0
for (j in p) {
  k <- k + 1
  z_p[k] <- quant_cond(MGEV3, 4, j, n1, dep = FALSE, cond.mar2 = FALSE)
}
plot(seq(1,10, by = 0.1), z_p, xlab = expression(10^-x),
     ylab = "Quantile",
     main = "Quantile conditionnel dans le cas d'indépendance asymptotique",
     type = "l")
