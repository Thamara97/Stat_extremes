IC_GEV <- function(gev, alpha, gum = FALSE) {
  ## Calcul des intervalles de confiance des paramètres par maximum de vraisemblance
  
  ## Entrée
  # gev : objet gev.fit
  # alpha : niveau de contrôle
  
  ## Sortie
  # IC : intervalles de confiances
  
  IC <- cbind(data$mle - qnorm(p = 1 - alpha/2) * data$se,
              data$mle + qnorm(p = 1 - alpha/2) * data$se)
  print(IC)
  return(IC)
}

quant <- function(gev, q, alpha) {
  ## Calcul des niveaux de retour et leur intervalle de confiance
  
  ## Entrée
  # gev : objet gev.fit
  # q : probabilité qu'un max dépasse le niveau de retour
  #     associé à la période de retour 1/q
  # alpha : niveaux de l'intervalle de confiance
  
  ## Sortie
  # x_q : niveau de retour
  # IC : intervalles de confiances
  
  mu <- gev$mle[1]
  sigma <- gev$mle[2]
  gamma <- gev$mle[3]
  
  yq <- -log(1-q)
  xq <- mu - (sigma / gamma) * (1 - yq^(-gamma)) 
  
  gradxq <- c(1, - gamma^(-1) * (1 - yq^(-gamma)),
              sigma * gamma^(-2) * (1 - yq^(-gamma)) - (sigma / gamma) * 
                yq^(-gamma) * log(yq))
  varxq <- t(gradxq) %*% gev$cov %*% gradxq
  IC <- c(xq - qnorm(1 - alpha/2) * sqrt(varxq),
          xq + qnorm(1 - alpha/2) * sqrt(varxq))
  
  return(list("xq" = xq, "IC" = IC))
}

quant_cond <- function(model, y, proba, n, dep = TRUE, cond.mar2 = FALSE){
  ## Calcule du quantile conditionnel z_p tel que P(X > z_p | Y > y) = p
  
  ## Entrée :
  # model : objet de type MGEV
  # y : quantile variable de conditionnement
  # proba : probabilité de dépassement conditionnelle
  # n : taille de l'échantillon
  # dep : dépendance (TRUE) ou d'indépendance (FALSE) des marginales
  # cond.mar2 : conditionnement par mar2 (TRUE) ou mar1 (FALSE) du modèle
  
  ## Sortie
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
