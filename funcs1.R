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

quant_cond <- function(model, y, proba, n, dep = TRUE, cond.mar2 = TRUE){
  ## Calcule du quantile conditionnel z_p tel que P(X > z_p | Y > y) = p
  
  # Entrée :
  # model : objet de type MGEV
  # y : quantile variable de conditionnement
  # proba : probabilité de dépassement conditionnelle
  # n : taille de l'échantillon
  # dep : dépendance (TRUE) ou d'indépendance (FALSE) des marginales
  # cond.mar2 : conditionnement par mar2 (TRUE) ou mar1 (FALSE) du modèle
  
  # Sortie
  # z_p : quantile conditionnel
  
  if (cond.mar2) {
    marX <- model$estimate[1:3]
    marY <- model$estimate[4:6] 
  } else {
    marY <- model$estimate[1:3]
    marX <- model$estimate[4:6]
  }
  
  # Cas de dépendance asymptotique et en supposant z_p > y
  if (dep) {
    g_y <- pgev(y, loc = marY[1], scale = marY[2], shape = marY[3])^(1/n)
    u <- (1-(1-g_y)*proba)^n
  } else {
    # Cas d'indépendance asymptotique
    u <- (1-proba)^n
  }
  
  z_p <- qgev(u, loc = marX[1], scale = marX[2], shape = marX[3])
  names(z_p) <- "quant"
  
  return(z_p)
}
