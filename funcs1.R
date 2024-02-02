extractMax <- function(data, nb_bloc, n, plot_bloc = TRUE, plot_max = TRUE,
                       bloc.x = "x", bloc.y = "data[x]",
                       bloc.title = paste(data, "par bloc"),
                       max.x = "x", max.y = "max[x]",
                       max.title = "Maximas") {
  ## Extraction des maximas par bloc
  
  maxX <- rep(0, nb_bloc)
  ind <- rep(0, nb_bloc)
  
  for (i in 1:nb_bloc) {
    maxX[i] <- max(data[((i-1) * n + 1):(i*n)])
    ind[i] <- which.max(data[((i-1) * n + 1):(i*n)]) + ((i - 1) * n)
  }
  
  if (plot_bloc) {plot(data, xlab = bloc.x, ylab = bloc.y, main = bloc.title)
                  points(x = ind, y = maxX, col = "red")
                  }
  
  if (plot_max) {plot(maxX, xlab = max.x, ylab = max.y, main = max.title)}
  
  return(maxX)
}

IC_GEV <- function(data, alpha) {
  ## Renvoie les intervalles de confiance des paramètres d'un objet gev.fit
  ## Calculés par maximum de vraisemblance
  
  # data : objet gev.fit
  # alpha : niveau de contrôle
  
  IC <- cbind(data$mle - qnorm(p = 1 - alpha/2) * data$se,
              data$mle + qnorm(p = 1 - alpha/2) * data$se)
  rownames(IC) <- c("mu", "sigma", "gamma")
  colnames(IC) <- c("Borne Inf", "Borne Sup")
  IC
  return(IC)
}

quant <- function(gev, q, mu, sigma, gama, alpha) {
  ## Calcul des niveaux de retour et leur intervalle de confiance
  
  # gev : objet gev.fit
  # q : probabilité qu'un max dépasse le niveau de retour
  #     associé à la période de retour 1/q
  # mu : paramètre de moyenne (location)
  # sigma : paramètre de dispersion
  # gama : paramètre de forme
  # alpha : niveaux de l'intervalle de confiance
  
  yq <- -log(1-q)
  xq <- mu - (sigma / gama) * (1 - yq^(-gama)) 
  
  gradxq <- c(1, - gama^(-1) * (1 - yq^(-gama)),
              sigma * gama^(-2) * (1 - yq^(-gama)) - (sigma / gama) * 
                yq^(-gama) * log(yq))
  varxq <- t(gradxq) %*% gev$cov %*% gradxq
  IC <- c(xq - qnorm(1 - alpha/2) * sqrt(varxq),
          xq + qnorm(1 - alpha/2) * sqrt(varxq))
  
  return(list("xq" = xq, "IC" = IC))
}
