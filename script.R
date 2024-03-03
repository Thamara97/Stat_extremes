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

################################## Station SA ##################################

#### Extraction des maximas annuels ####

# Extraction des données et des dates
df <- data.frame(date = donneesVague$date, station = SA)

# Extraction de l'année
df$year <- format(as.Date(df$date, format="%Y-%m-%d"),"%Y")

# Max par an
maxA <- df %>% group_by(year) %>% summarise(Max = max(station)) %>% pull()
n <- length(maxA)
# Nombre de mesures par an
eff <- df %>% group_by(year) %>% count() %>% pull()
# Représentation graphique des max
maxima_annuels <- data.frame(maxima=maxA,
                             date=1961:2012)

ggplot(maxima_annuels, aes(x=date, y=maxima))+
  geom_point(aes(y=maxima), color="black", pch=3, size=0.25)+
  theme(strip.text.x=element_text(size=8),
        plot.title=element_text(size=10, hjust=0.5, face="bold"))+
  xlab("Année") +
  ylab("Hauteur de vagues") +
  labs(title="")

#### Ajustement d'une GEV ####

maxA_GEV1 <- gev.fit(maxA)
mu1 <- maxA_GEV1$mle[1]
sigma1 <- maxA_GEV1$mle[2]
gamma1 <- maxA_GEV1$mle[3]

# Intervalle de confiance des paramètres
IC_GEV(maxA_GEV1, alpha = 0.05)

# Représentations graphiques
gev.diag(maxA_GEV1)

# Test de l'hypothèse de stationnarité
mat <- matrix(c(1:52, (1:52)^2), ncol = 2)
maxA_GEV0 <- gev.fit(maxA, ydat = mat, mul = 1)
IC_GEV(maxA_GEV0, alpha = 0.05)
# Test du rapport de vraisemblance
- 2 * (maxA_GEV0$nllh - maxA_GEV1$nllh) <= qchisq(0.95, df = 1)

# Ajustement d'une Gumbel
maxA_gum <- gum.fit(maxA)
IC_GEV(maxA_gum, alpha = 0.95, gum = TRUE)
gum.diag(maxA_gum) 
# Test du rapport de vraisemblance
- 2 * (maxA_GEV1$nllh - maxA_gum$nllh) <= qchisq(0.95, 1)

#### Niveaux de retour ####

# Périodes et probas
T1 <- 100 ; q1 <- 1/T1
T2 <- 500 ; q2 <- 1/T2
T3 <- 1000 ; q3 <- 1/T3

alpha <- 0.05

# Quantiles des max annuels

# 100 ans
quant(maxA_GEV1, q1, alpha)
quant100 <- fgev(maxA, prob = q1)

# 500 ans
quant(maxA_GEV1, q2, alpha)
quant500 <- fgev(maxA, prob = q2)

# 1000 ans
quant(maxA_GEV1, q3, alpha)
quant1000 <- fgev(maxA, prob = q3)

################################## Station SB ##################################

SB <- donneesVague$station12

# plot(SA, SB) # long a produire
plot(SA[1:5000], SB[1:5000])

#### Extraction des maximas pour la station SB ####

df2 <- data.frame(date = donneesVague$date, station = SB)

# Extraction de l'année
df2$year <- format(as.Date(df2$date, format="%Y-%m-%d"),"%Y")

# Extraction des maxima
maxB <- df2 %>% group_by(year) %>% summarise(Max = max(station)) %>% pull()
length(maxB)

# Représentation graphique
plot(1961:2012, maxB, xlab = "Années", ylab = "Hauteur mesurée",
     main = "Hauteur maximales des vagues par années")

#### Ajustement ####

maxB_GEV1 <- gev.fit(maxB)

mu2 <- maxB_GEV1$mle[1]
sigma2 <- maxB_GEV1$mle[2]
gamma2 <- maxB_GEV1$mle[3]

IC_GEV(maxB_GEV1, alpha = 0.05)

gev.diag(maxB_GEV1)

# Test de la stationnarité
maxB_GEV0 <- gev.fit(maxB, ydat = mat, mul = 1)
IC_GEV(maxB_GEV0, alpha = 0.05)
- 2 * (maxB_GEV0$nllh - maxB_GEV1$nllh) <= qchisq(0.95, df = 1)

# Ajustement d'une Gumbel
maxB_gum <- gum.fit(maxB)
IC_GEV(maxB_gum, alpha = 0.05)
gum.diag(maxB_gum)
# Test du rapport de vraisemblance
- 2 * (maxB_GEV1$nllh - maxB_gum$nllh) <= qchisq(0.95, 1)

#### Etude de la dépendance ####

# Chi plot et Chi bar plot
par(mfrow = c(1,2))
chiplot(cbind(SA, SB), xlab = "u", ylab1 = expression(chi(u)),
        ylab2 = expression(bar(chi)(u)))

# Transformation des données en Fréchet 1
maxA_F <- qgev(pgev(maxA,loc = mu1, scale = sigma1,
                     shape = gamma1),
                loc = 1, shape = 1, scale = 1)
maxB_F <- qgev(pgev(maxB, loc = mu2, scale = sigma2,
                     shape = gamma2),
                loc = 1, shape = 1, scale = 1)
# Représentation graphique
par(mfrow = c(1,1))
plot(maxA_F, maxB_F)

# Estimateur de la fonction de dépendance de Pickands
abvnonpar(data = cbind(maxA, maxB), plot = TRUE, col = "red",
          main = "Estimateur de la fonction de dépendance de Pickands")
abvnonpar(data = cbind(maxA, maxB), plot = FALSE) * 2

#### Ajustement d'une MGEV ####

fbvevd(cbind(maxA, maxB), model = "alog", std.err = FALSE) # AIC = 221.7819
fbvevd(cbind(maxA, maxB), model = "log") # AIC = 217.7394
MGEV1 <- fbvevd(cbind(maxA, maxB), model = "hr") # AIC = 217.724

#### Quantile conditionnel ####

quant_cond(MGEV1, 4, 10^(-5), n, dep = FALSE, cond.mar2 = FALSE)

p <- 10^-seq(1,10, by = 0.1)
z_p <- c()
k <- 0
for (j in p) {
  k <- k + 1
  z_p[k] <- quant_cond(MGEV1, 4, j, n, dep = FALSE, cond.mar2 = FALSE)
}

plot(seq(1,10, by = 0.1), z_p, type = "l",
     xlab = expression(10^-x),
     ylab = "Quantile",
     main = "Quantile conditionnel en cas d'indépendance asymptotique")

################################## Station SC ##################################

SC <- donneesVague$station1

plot(SA[1:5000], SC[1:5000])

#### Extraction des maximas pour la station SC ####

df3 <- data.frame(date = donneesVague$date, station = SC)
# Extraction de l'année
df3$year <- format(as.Date(df3$date, format="%Y-%m-%d"),"%Y")

# Maxima
maxC <- df3 %>% group_by(year) %>% summarise(Max = max(station)) %>% pull()

# Représentation graphique
plot(1961:2012, maxC, xlab = "Années", ylab = "Hauteur mesurée",
     main = "Hauteur maximales des vagues par années")

#### Etude de la dépendance ####

# Chi et chi bar plot
chiplot(cbind(SA, SC))

# Transformation des données en Fréchet 1
maxC_GEV1 <- gev.fit(maxC)
mu3 <- maxC_GEV1$mle[1]
sigma3 <- maxC_GEV1$mle[2]
gamma3 <- maxC_GEV1$mle[3]

IC_GEV(maxC_GEV1, alpha = 0.05)
gev.diag(maxC_GEV1)

# Hypothèse de stationnarité
maxC_GEV0 <- gev.fit(maxC, ydat = mat, mul = 1)
-2*(maxC_GEV0$nllh - maxC_GEV1$nllh) <= qchisq(p = 0.95, df = 1)

# Ajustement d'une Gumbel
maxC_gum <- gum.fit(maxC)
gum.diag(maxC_gum)
- 2 * (maxC_GEV1$nllh - maxC_gum$nllh) <= qchisq(0.95, 1)

# Transformation des marginales en Fréchet 1
maxC_F <- qgev(pgev(maxC, loc = mu3, scale = sigma3,
                     shape = gamma3),
                loc = 1, shape = 1, scale = 1)
# Représentation graphique
plot(maxA_F, maxC_F)

# Estimation de la fonction de dépendance de Pickands
abvnonpar(data = cbind(maxA, maxC), plot = TRUE, col = "red")
abvnonpar(data = cbind(maxA, maxC), plot = FALSE) * 2

#### Ajustement d'une MGEV ####

fbvevd(cbind(maxA, maxC), model = "alog", std.err = FALSE) # AIC = 31.54714
fbvevd(cbind(maxA, maxC), model = "log") # AIC = 34.66942
MGEV2 <- fbvevd(cbind(maxA, maxC), model = "hr") # AIC = 3.649479
MGEV2

#### Quantile conditionnel ####

quant_cond(MGEV2, 8, 10^-10, n, dep = TRUE, cond.mar2 = FALSE)

y <- 4
p <- 10^-seq(1,10, by = 0.1)
l <- 0
for (i in y){
  l <- l + 1
  z_p <- c()
  k <- 0
  for (j in p) {
    k <- k + 1
    z_p[k] <- quant_cond(MGEV2, i, 10^(-j), n, dep = TRUE, cond.mar2 = FALSE)
  }
  plot(seq(1,10, by = 0.1), z_p, xlab = expression(10^-x), ylab = "Quantile",
       main = paste("Quantile conditionnel en quand Y >",y[l]), type = "l")
}

################################## Station SD ################################## 

SD <- donneesVague$station4

plot(SA[1:5000], SD[1:5000])

#### Extraction des maximas pour la station SD ####

df4 <- data.frame(date = donneesVague$date, station = SD)

# Extraction de l'année
df4$year <- format(as.Date(df4$date, format="%Y-%m-%d"),"%Y")

# Max par an
maxD <- df4 %>% group_by(year) %>% summarise(Max = max(station)) %>% pull()

# Représentation graphique
plot(1961:2012, maxD, xlab = "Années", ylab = "Hauteur mesurée",
     main = "Hauteur maximales des vagues par années")

plot(maxA, maxD)

#### Etude de la dépendance ####

# Chi et chi bar plot
chiplot(cbind(SA, SD))

# Transformation des données en Fréchet 1
maxD_GEV1 <- gev.fit(maxD)
mu4 <- maxD_GEV1$mle[1]
sigma4 <- maxD_GEV1$mle[2]
gamma4 <- maxD_GEV1$mle[3]

IC_GEV(maxD_GEV1, alpha = 0.05)

# Ajutemetn d'une Gumbel
maxD_gum <- gum.fit(maxD)
-2 * (maxD_GEV1$nllh - maxD_gum$nllh)

# Test de l'hypothèse de stationnarité
maxD_GEV0 <- gev.fit(maxD, ydat = mat, mul = 1)
-2*(maxD_GEV0$nllh - maxD_GEV1$nllh) <= qchisq(p = 0.99, df = 1)

# Transformation des marginales
maxD_F <- qgev(pgev(maxD, loc = mu4, scale = sigma4, shape = gamma4),
                loc = 1, shape = 1, scale = 1)
# Représentation graphique
plot(maxA_F, maxD_F)

# Estimation de la fonction de dépendance de Pickands
abvnonpar(data = cbind(maxA, maxD), plot = TRUE, col = "red")
abvnonpar(data = cbind(maxA, maxD), plot = FALSE) * 2

#### Ajustement d'une MGEV ####
fbvevd(cbind(maxA, maxD), model = "alog", std.err = FALSE) # AIC = 219.3239
fbvevd(cbind(maxA, maxD), model = "log") # AIC = 215.073
MGEV3 <- fbvevd(cbind(maxA, maxD), model = "hr") # AIC = 214.4105
MGEV3

#### Quantile conditionnel ####

quant_cond(MGEV3, 4, 10^-10, n, dep = FALSE, cond.mar2 = FALSE)

p <- 10^-seq(1,10, by = 0.1)
z_p <- c()
k <- 0
for (j in p) {
  k <- k + 1
  z_p[k] <- quant_cond(MGEV3, 4, j, n, dep = FALSE, cond.mar2 = FALSE)
}
plot(seq(1,10, by = 0.1), z_p, xlab = expression(10^-x),
     ylab = "Quantile",
     main = "Quantile conditionnel dans le cas d'indépendance asymptotique",
     type = "l")