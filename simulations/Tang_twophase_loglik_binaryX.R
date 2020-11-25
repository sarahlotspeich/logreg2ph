Tang_twophase_loglik <- function(params, Val, Y_unval, Y_val, X_unval, X_val, C, data)
{
  theta0 <- params[1]
  theta1 <- params[2]
  theta2 <- params[3]
  theta3 <- params[4]
  theta4 <- params[5]
  delta0 <- params[6] 
  delta1 <- params[7] 
  delta3 <- params[8]
  beta0 <- params[9] 
  beta1 <- params[10] 
  beta2 <- params[11]
  gamma0 <- params[12]
  
  data_unval <- data[which(data[,Val] == 0),]
  data_val <- data[which(data[,Val] == 1),]
  
  # Main study 
  ## P(Y*|y, x, X*, C)
  eta1_y0x0 <- theta0 + theta1*0 + theta2*0 + theta3*data_unval[,X_unval] + theta4*data_unval[,C]
  eta1_y0x1 <- theta0 + theta1*0 + theta2*1 + theta3*data_unval[,X_unval] + theta4*data_unval[,C]
  eta1_y1x0 <- theta0 + theta1*1 + theta2*0 + theta3*data_unval[,X_unval] + theta4*data_unval[,C]
  eta1_y1x1 <- theta0 + theta1*1 + theta2*1 + theta3*data_unval[,X_unval] + theta4*data_unval[,C]
  
  pYstar_y0x0 <- ifelse(data_unval[,Y_unval] == 1, (1 + exp(-eta1_y0x0))^(-1), 1-(1 + exp(-eta1_y0x0))^(-1))
  pYstar_y0x1 <- ifelse(data_unval[,Y_unval] == 1, (1 + exp(-eta1_y0x1))^(-1), 1-(1 + exp(-eta1_y0x1))^(-1))
  pYstar_y1x0 <- ifelse(data_unval[,Y_unval] == 1, (1 + exp(-eta1_y1x0))^(-1), 1-(1 + exp(-eta1_y1x0))^(-1))
  pYstar_y1x1 <- ifelse(data_unval[,Y_unval] == 1, (1 + exp(-eta1_y1x1))^(-1), 1-(1 + exp(-eta1_y1x1))^(-1))
  
  pYstar <- data.frame(pYstar_y0x0, pYstar_y1x0, pYstar_y0x1, pYstar_y1x1)
  ### -----------------------------------------------------------------
  
  ## P(X*|x, C)
  eta2_x0y0 <- delta0 + delta1*0 + delta3*data_unval[,C]
  eta2_x0y1 <- delta0 + delta1*0 + delta3*data_unval[,C]
  eta2_x1y0 <- delta0 + delta1*1 + delta3*data_unval[,C]
  eta2_x1y1 <- delta0 + delta1*1 + delta3*data_unval[,C]
  
  pXstar_x0y0 <- ifelse(data_unval[,X_unval] == 1, (1 + exp(-eta2_x0y0))^(-1), 1-(1 + exp(-eta2_x0y0))^(-1))
  pXstar_x0y1 <- ifelse(data_unval[,X_unval] == 1, (1 + exp(-eta2_x0y1))^(-1), 1-(1 + exp(-eta2_x0y1))^(-1))
  pXstar_x1y0 <- ifelse(data_unval[,X_unval] == 1, (1 + exp(-eta2_x1y0))^(-1), 1-(1 + exp(-eta2_x1y0))^(-1))
  pXstar_x1y1 <- ifelse(data_unval[,X_unval] == 1, (1 + exp(-eta2_x1y1))^(-1), 1-(1 + exp(-eta2_x1y1))^(-1))
  
  pXstar <- data.frame(pXstar_x0y0, pXstar_x0y1, pXstar_x1y0, pXstar_x1y1)
  ### ------------------------------------------------------------------
  
  ## P(y|x)
  eta3_x0 <- beta0 + beta1*0 + beta2*data_unval[,C]
  eta3_x1 <- beta0 + beta1*1 + beta2*data_unval[,C]
  
  Py0_x0 <- 1-(1 + exp(-eta3_x0))^(-1)
  Py1_x0 <- (1 + exp(-eta3_x0))^(-1)
  Py0_x1 <- 1-(1 + exp(-eta3_x1))^(-1)
  Py1_x1 <- (1 + exp(-eta3_x1))^(-1)
  
  #pY <- data.frame(Py0_x0 = rep(Py0_x0, nrow(data_unval)), Py1_x0 = rep(Py1_x0, nrow(data_unval)), Py0_x1 = rep(Py0_x1, nrow(data_unval)), Py1_x1 = rep(Py1_x1, nrow(data_unval)))
  pY <- data.frame(Py0_x0, Py1_x0, Py0_x1, Py1_x1)
  
  ## P(x)
  Px1 <- (1 + exp(-gamma0))^(-1)
  Px0 <- 1 - (1 + exp(-gamma0))^(-1)
  
  pX <- data.frame(Px0 = rep(Px0, nrow(data_unval)), Px0 = rep(Px0, nrow(data_unval)), 
                   Px1 = rep(Px1, nrow(data_unval)), Px1 = rep(Px1, nrow(data_unval)))
  ## Sum over x = 0,1/ y = 0/1
  
  # Multiply P(Y*|y,x,X*) x P(X*|x,y) x P(y|x) x P(x)
  like_main <- rowSums(pYstar * pXstar * pY * pX)
  log_like_main <- log(like_main)
  sum_log_like_main <- sum(log_like_main)
  
  # Validation subsample (Phase II)
  ## P(Y*|Y)
  eta1 <- theta0 + theta1*data_val[,Y_val] + theta2*data_val[,X_val] + theta3*data_val[,X_unval] + theta4*data_val[,C]
  pYstar <- ifelse(data_val[,Y_unval] == 1, (1 + exp(-eta1))^(-1), 1-(1 + exp(-eta1))^(-1))
  
  ## P(X*|X)
  eta2 <- delta0 + delta1*data_val[,X_val] + delta3*data_val[,C]
  pXstar <- ifelse(data_val[,X_unval] == 1, (1 + exp(-eta2))^(-1), 1-(1 + exp(-eta2))^(-1))
  
  ## P(Y|X)
  eta3 <- beta0 + beta1*data_val[,X_val] + beta2*data_val[,C]
  pY <- ifelse(data_val[,Y_val] == 1, (1 + exp(-eta3))^(-1), 1-(1 + exp(-eta3))^(-1))
  
  ## P(X)
  pX <- ifelse(data_val[,X_val] == 1, (1 + exp(-gamma0))^(-1), 1 - (1 + exp(-gamma0))^(-1))
  
  # Multiply P(Y*|y,x,X*) x P(X*|x,y) x P(y|x) x P(x)
  like_val <- pYstar * pXstar * pY * pX
  log_like_val <- log(like_val)
  sum_log_like_val <- sum(log_like_val)
  
  return(-(sum_log_like_main + sum_log_like_val))
}