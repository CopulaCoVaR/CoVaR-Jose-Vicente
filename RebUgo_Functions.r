#-- t-asimetrica de Hansen, NO SE USA, es utiliza la de Rugarch ---#
f_Skewed_t_Hansen <- function(x, eta=0.5, nu=5){
  c = gamma((nu+1)/2) / (sqrt(pi*(nu-2)) * gamma(nu/2))
  a = 4 * eta* c * ( (nu-2)/(nu-1) )
  b = sqrt(1 + 3*eta^2 - a^2)
  #f = (b*c* ((1+(1/(nu-2))*((b*x+a)/(1-eta))^2)^(-(nu+1)/2)))*(x < -a/b) +
  #    (b*c* ((1+(1/(nu-2))*((b*x+a)/(1+eta))^2)^(-(nu+1)/2)))*(x >= -a/b)
  log.f = (log(b)+log(c)-0.5*(nu+1)*log(1+(1/(nu-2))*((b*x+a)/(1-eta))^2) )*(x < -a/b) + 
          (log(b)+log(c)-0.5*(nu+1)*log(1+(1/(nu-2))*((b*x+a)/(1+eta))^2) )*(x >= -a/b) 
  return(log.f)
}
#integrate(f_Skewed_t_Hansen, -Inf, Inf)

#---- ML de un modelo ARMA(0,0)-TGARCH(1,1)-Swt
ARMA_GARCH_Skt_ML <- function(param){
    #serie <<- serie.to.ML
    mu     = param[1]
    omega  = param[2]
    alpha1 = param[3]
    beta1  = param[4]
    lambda = param[5]
    eta    = param[6]
    nu     = param[7]
    T = length(serie.to.ML)
    
    #Modelo ARMA#
    epsilon = serie.to.ML - mu
    prop.ep.neg = sum(epsilon<0)/T
    #Modelo TGARCH#

    sigma2 = matrix(omega/(1-alpha1-beta1-lambda*prop.ep.neg),1,T)
    for(t in 2:T)
      sigma2[t] = omega +  beta1*sigma2[t-1] + (epsilon[t-1]^2)*(alpha1+lambda*(epsilon[t-1]<0))
    error.std = as.matrix(epsilon)/ t(as.matrix(sqrt(sigma2)))  
    LogLik = (f_Skewed_t_Hansen(error.std, eta, nu))
    LogLik[LogLik>1000]=0; LogLik[LogLik< -1000]=0 
    res = -sum(LogLik)
    cat('\n', res)
    return (res)
}
#---- t.Copula_pdf
T_pdf=function(u,v,k1,k2){
  cop=BiCop(family = 2, par = k1,par2 = k2)
  pdf=BiCopPDF(u, v, cop)
  return(pdf)
}
#--- Funciones de Woraphon Yamaka del paquete <DynamicCOP> en chrome-extension://efaidnbmnnnibpcajpcglclefindmkaj/viewer.html?pdfurl=https%3A%2F%2Fwyamaka.files.wordpress.com%2F2020%2F07%2Fpresentation.pdf&clen=2980059&chunk=true
dynamicnormal.LF <- function (data, plot, printout=FALSE) 
{
  lower = rep(-5, 3)
  upper = rep(5, 3)
  theta = c(w = 0.01, beta1 = 0.002, beta2 = -0.002)
  rhobar = cor(data[,1],data[,2], method='kendall')
  
  model <- optim(theta, normtvtp.LF, data = data, rhobar = rhobar, 
                 Call = "optim", control = list(maxit = 1e+05, fnscale = -1), 
                 method = "L-BFGS-B", lower = lower, upper = upper, hessian = TRUE,
                 printout=printout)
  n = nrow(data)
  coef <- model$par
  model$se <- sqrt(abs(-diag(solve(model$hessian))))
  S.E. = model$se
  stat = coef/S.E.
  pvalue <- 2 * (1 - pnorm(abs(stat)))
  result <- cbind(coef, S.E., stat, pvalue)
  BIC = -2 * model$value + (log(n) * length(coef))
  AIC = -2 * model$value + 2 * length(coef)
  rhot = normtvtp.LF(coef, data = data, rhobar = 0.5,printout=printout, Call = "filtering")$TVTP
  if (plot == TRUE) {
    plot(ts(rhot[5:n]))
  }
  output = list(result = result, AIC = AIC, BIC = BIC, Loglikelihood = model$value, 
                tvtpdep = rhot)
  output
}
normtvtp.LF <- function (theta, data, rhobar, printout=FALSE, Call = "optim") 
{
  # Correccion de pseudo-observaciones iguales a 0 o 1.
  data[(data[,1]<1e-10),1] = 0.000001
  data[(data[,2]<1e-10),2] = 0.000001
  data[(data[,1]>0.99999999),1] = 0.999999
  data[(data[,2]>0.99999999),2] = 0.999999
  x = qnorm(data[, 1])
  y = qnorm(data[, 2])
  T = nrow(data)
  u = data[, 1]
  v = data[, 2]
  w = theta[1]
  a = theta[2]
  b = theta[3]
  kappa = rep(-0.99, T)
  kappa[1] = rhobar
  psi1 = rep(0, T)
  for (jj in 2:T) {
    if (jj <= 10) {
      psi1[jj] = w + a * kappa[jj-1] + b * (mean(abs(x[1:(jj-1)] * y[1:(jj-1)])))
    }
    else {
      psi1[jj] = w + (a * kappa[jj-1]) + b * (mean(abs(x[(jj-10):(jj-1)] * y[(jj-10):(jj-1)])))
    }
    #kappa[jj] = 1.998/(1 + exp(-psi1[jj])) - 0.999
    kappa[jj] = (1 - exp(-psi1[jj])) / (1 + exp(-psi1[jj])) #LFM
  }
  rhohat = kappa
  # Log-Lik de la copula
  #CL = -1 * (2 * (1 - kappa^2))^(-1) * ((x^2) + (y^2) - 2 * kappa * x * y)
  #CL = CL + 0.5 * ((x^2) + (y^2))
  #CL = CL - 0.5 * log(1 - (kappa^2))
  #CL = sum((CL))
  CL = rep(0, T)
  for (i in 1:T){ #Ajuste por incompatibilidad de la Copula Gausiana. 
    if (rhohat[i]==1) rhohat[i]=0.99999999
    if (rhohat[i]==(-1)) rhohat[i]=-0.99999999
    cop   = BiCop(family = 1, par = rhohat[i])
    CL[i] = BiCopPDF(u1=u[i], u2=v[i], obj=cop) #PDF de una copula normal
    #CL[i] = BiCopPDF(u1=u[i], u2=v[i], family=1, par=rhohat[i]) #PDF de una copula normal
  }
  CL = sum(log(CL))
  
  n = T
  if (is.infinite(CL)) 
    CL <- -n * 100
  if (is.nan(CL)) 
    CL <- -n * 100
  if (printout) cat("Sum of log Likelihood for normdynamic ->", sprintf("%4.4f", 
                                                          c(CL, kappa[n])), "\n")
  if (Call == "optim") 
    return(CL)
  if (Call == "filtering") {
    specOut <- list(LL = CL, beta = c(theta), TVTP = rhohat)
    return(specOut)
  }
}

dynamicT.LF <-function (data, plot=FALSE, printout=FALSE) 
{
  # ML #
  lower = c(rep(-5, 3), 2)
  upper = c(rep( 5, 3), 50)
  theta = c(w = 0.01, beta1 = 0.002, beta2 = -0.002, nu=5)
  model <- optim(theta, Ttvtp.LF, data = data, Call = "optim", 
                 control = list(maxit = 1e+05, fnscale = -1, trace=1), method = "L-BFGS-B", 
                 lower = lower, upper = upper, hessian = TRUE,
                 printout=TRUE)
  # Impresion de resultados #
  n = nrow(data)
  coef <- model$par
  model$se <- sqrt(abs(-diag(solve(model$hessian))))
  S.E. = model$se
  stat = coef/S.E.
  pvalue <- 2 * (1 - pnorm(abs(stat)))
  result <- cbind(coef, S.E., stat, pvalue) # La ult fila de <result> es <NU> i.e. d.f 
  # Criterios de info #
  BIC = -2 * model$value + (log(n) * length(coef))
  AIC = -2 * model$value + 2 * length(coef)
  # Coef de correlacion cambiante en el tiempo <rhot>#
  rhot = Ttvtp.LF(coef, data = data, printout=printout, Call = "filtering")$TVTP
  if (plot == TRUE) {
    plot(ts(rhot[5:n]))
  }
  output = list(result = result, AIC = AIC, BIC = BIC, Loglikelihood = model$value, 
                tvtpdep = rhot)
  output
}
Ttvtp.LF <- function (theta, data, printout=FALSE, Call = "optim") 
{
  #cat('\nw:',theta[1],'Beta1:',theta[2],'Beta2:',theta[3],'nu:',theta[4])
  # Correccion de pseudo-observaciones iguales a 0 o 1.
  data[(data[,1]<1e-10),1] = 0.000001
  data[(data[,2]<1e-10),2] = 0.000001
  data[(data[,1]>0.99999999),1] = 0.999999
  data[(data[,2]>0.99999999),2] = 0.999999
  d = theta[4] #LF
  x = qt(data[, 1], df=d)
  y = qt(data[, 2], df=d)

  n = length(x)
  T = nrow(data)
  u = data[, 1]
  v = data[, 2]
  w = theta[1]
  a = theta[2]
  b = theta[3]

  kappa = rep(-0.99, n)
  #kappa[1] = 0.5
  kappa[1] = cor(u,v, method='kendall')
  psi1 = rep(0, n)
  for (jj in 2:n) {
    if (jj <= 10) 
      psi1[jj] = w + a * kappa[jj - 1] + b * ( mean(abs(x[1:(jj-1)] * y[1:(jj-1)])) )
    else
      psi1[jj] = w + a * kappa[jj - 1] + b * (mean(abs(x[(jj - 10):(jj - 1)] * y[(jj - 10):(jj - 1)])))
    #kappa[jj] = 1.998/(1 + exp(-psi1[jj])) - 0.999
    kappa[jj] = (1 - exp(-psi1[jj])) / (1 + exp(-psi1[jj])) #LFM
  }
  RHO = kappa
  NU = d
  CL = rep(0, T)
  for (i in 1:T) {
    if (RHO[i]==1) RHO[i] =0.99999999
    if (RHO[i]==(-1)) RHO[i]=-0.99999999
    CL[i] = T_pdf(u[i], v[i], RHO[i], NU)
  }
  CL = sum(log(CL))
  n = T
  if (is.infinite(CL)) 
    CL <- -n * 100
  if (is.nan(CL)) 
    CL <- -n * 100
  if (printout) cat("Sum of log Likelihood for Tdynamic ->", sprintf("%4.4f", 
                                                       c(CL, kappa[n])), "\n")
  if (Call == "optim") 
    return(CL)
  if (Call == "filtering") {
    specOut <- list(LL = CL, beta = c(theta), TVTP = RHO)
    return(specOut)
  }
}

dynamicGum.LF <-function (data, plot=FALSE, rotated=FALSE, printout=FALSE) 
{
  lower = c(rep(-4, 3))
  upper = c(rep( 4, 3))
  theta = c(w = 0.01, beta1 = 0.002, beta2 = -0.002) #Cambirlo???
  model = optim(theta, Gumtvtp.LF, data = data, Call = "optim", 
                control = list(maxit = 1e+05, fnscale = -1), method = "L-BFGS-B", 
                lower=lower, upper=upper, hessian=TRUE, rotated=rotated, printout=printout)
  n        = nrow(data)
  coef     = model$par
  model$se = sqrt(abs(-diag(solve(model$hessian))))
  S.E.     = model$se
  stat     = coef/S.E.
  pvalue = 2 * (1 - pnorm(abs(stat)))
  result = cbind(coef, S.E., stat, pvalue)
  BIC    = -2 * model$value + (log(n) * length(coef))
  AIC    = -2 * model$value + 2 * length(coef)
  rhot   = Gumtvtp.LF(coef, data=data, rotated=rotated, printout=printout, Call="filtering")$TVTP
  if (plot == TRUE) {
    plot(ts(rhot[5:n]))
  }
  output = list(result = result, AIC = AIC, BIC = BIC, Loglikelihood = model$value, 
                tvtpdep = rhot)
  output
}
Gumtvtp.LF <- function (theta, data, rotated=FALSE, printout=FALSE, Call = "optim") 
{
  T = nrow(data); n=T
  u = data[, 1]
  v = data[, 2]
  w = theta[1]
  a = theta[2]
  b = theta[3]
  
  kappa = rep(-0.99, n)
  #kappa[1] = 0.5
  #kappa[1] = 1/(1-cor(u,v,method='kendall'))
  kappa[1] = 1.2
  psi1 = rep(0, n)
  for (jj in 2:n) {
    if (jj <= 10) 
      psi1[jj] = w + a * kappa[jj-1] + b * (mean(abs(u[1:(jj-1)] - v[1:(jj-1)])))
    else 
      psi1[jj] = w + a * kappa[jj-1] + b * (mean(abs(u[(jj-10):(jj-1)] - v[(jj-10):(jj-1)])))
    kappa[jj] = 1 + sqrt(abs(psi1[jj])) + 1e-3
  }
  RHO = kappa
  CL = rep(0, T)
  for (i in 1:T) {
    #CL[i] = T_pdf(u[i], v[i], RHO[i], NU)
    if (RHO[i]>= 17) cat('\npar:', RHO[i])
    if (rotated==FALSE) CL[i] = BiCopPDF(u[i], v[i], family=4, par=RHO[i]) #Gumbel density copula
    else                CL[i] = BiCopPDF(u[i], v[i], family=14,par=RHO[i]) #rotated density Gumbel copula (180 degrees)
  }
  CL = sum(log(CL))
  n = T
  if (is.infinite(CL)) 
    CL <- -n * 100
  if (is.nan(CL)) 
    CL <- -n * 100
  if (printout) cat("Sum of log Likelihood for Tdynamic ->", 
                    sprintf("%4.4f", c(CL, kappa[n])), "\n")
  if (Call == "optim") 
    return(CL)
  if (Call == "filtering") {
    specOut <- list(LL = CL, beta = c(theta), TVTP = RHO)
    return(specOut)
  }
}

dynamicBB7.LF <-function (data, plot=FALSE, printout=FALSE) 
{
  lower = rep(-10, 6)
  upper = rep(10, 6)
  theta = c(w  = 0.01, beta1 = 0.002, beta2 = -0.002, 
            w2 = 0.01, beta3 = 0.002, beta4 = -0.002)
  model <- optim(theta, BB7vtp.LF, data = data, Call = "optim", 
                 control = list(maxit = 1e+05, fnscale = -1), method = "L-BFGS-B", 
                 lower = lower, upper = upper, hessian = TRUE, printout=FALSE)
  n        = nrow(data)
  coef     = model$par
  model$se = sqrt(abs(-diag(solve(model$hessian))))
  S.E.     = model$se
  stat     = coef/S.E.
  pvalue   = 2 * (1 - pnorm(abs(stat)))
  result   = cbind(coef, S.E., stat, pvalue)
  BIC = -2 * model$value + (log(n) * length(coef))
  AIC = -2 * model$value + 2 * length(coef)
  model.est = BB7vtp.LF(coef, data = data, Call = "filtering", printout=FALSE)
  par1.t    = model.est$par1.t
  par2.t    = model.est$par2.t
  
  if (plot) {
    par(mfrow=c(2,1))
    plot(ts(par1.t[5:n]))
    plot(ts(par2.t[5:n]))
    par(mfrow=c(1,1))
  }
  output = list(result=result, AIC=AIC, BIC=BIC, Loglikelihood=model$value, 
                tvtpdep1=par1.t, tvtpdep2=par2.t)
  output
}
BB7vtp.LF <- function (theta, data, printout=FALSE, Call = "optim") 
{
    T = nrow(data); n=T
    u = data[, 1]
    v = data[, 2]
    w = theta[1]
    a = theta[2]
    b = theta[3]
    w1 = theta[4]
    a1 = theta[5]
    b1 = theta[6]
    kappa = cbind(rep(1.2, n), rep(0.2, n))
    kappa[, 1] = 1.2
    kappa[, 2] = 0.2
    psi1 = matrix(0, n, 2)
    for (jj in 2:n) {
      if (jj <= 10) {
        psi1[jj, 1] = w  + a *kappa[jj-1, 1] + b *(mean(abs(u[1:(jj-1)] - v[1:(jj-1)]))) #Verificar esta funcion!!!
        psi1[jj, 2] = w1 + a1*kappa[jj-1, 2] + b1*(mean(abs(u[1:(jj-1)] - v[1:(jj-1)]))) #Verificar esta funcion!!!
       }
      else {
        psi1[jj, 1] = w  + a *kappa[jj-1, 1] + b *(mean(abs(u[(jj-10):(jj-1)] * v[(jj-10):(jj-1)])))
        psi1[jj, 2] = w1 + a1*kappa[jj-1, 2] + b1*(mean(abs(u[(jj-10):(jj-1)] * v[(jj-10):(jj-1)])))
      }
      #kappa[jj, 1] = 1 + (abs(psi1[jj, 1])^(2)) + 1e-3
      kappa[jj, 1] = 2*(1/(1 + exp(-psi1[jj,1])) ) + 1 #in [1,3]
      #kappa[jj, 2] = 0 + (psi1[jj, 2])^2 + 1e-3
      kappa[jj, 2] = 3*(1/(1 + exp(-psi1[jj,2])) ) + 0  #in [0,3]
    }
    par1.t = kappa[, 1]
    par2.t = kappa[, 2]
    CL = rep(0, T)
    for (i in 1:T) 
      CL[i] = BiCopPDF(u[i], v[i], family=9, par=par1.t[i], par2=par2.t[i]) # density of a BB7 copula 
    CL = sum(log(CL))
    n = T
    if (is.infinite(CL)) CL <- -n * 100
    if (is.nan(CL))      CL <- -n * 100
    if(printout) cat("Sum of log Likelihood for Tdynamic ->", sprintf("%4.4f", c(CL, kappa[n])), "\n")
    
    if (Call == "optim") return(CL)
    if (Call == "filtering") {
      specOut <- list(LL=CL, beta=c(theta), par1.t=par1.t, par2.t=par2.t)
      return(specOut)
    }
}  
  
  
#-- https://drawar.github.io/posts/dynamic-copula/ --#
### Estimates a time varying Gumbel copula
if(0){
 # Computes the negative log likelihood of a time varying Gumbel copula  
 GumbelTVLogL <- function(psi, data) {
  u = data[, 1]
  v = data[, 2]
  theta <- GumbelTV(psi, data)
  tu <- (-log(u))^theta
  tv <- (-log(v))^theta
  out1 <- exp(-(tu+tv)^(1/theta))
  out2 <- (u*v)^(-1)
  out3 <- (tu+tv)^(-2+2/theta)
  out4 <- (log(u)*log(v))^(theta-1)
  out5 <- 1+(theta-1)*(tu+tv)^(-1/theta)
  out <- out1*out2*out3*out4*out5
  LL = sum(log(out))
  LL = -LL
  return(LL)
}
 # Patton' methodology
 GumbelTV <- function(theta,data){
  u <- data[,1]
  v <- data[,2]
  t <- dim(data)[1]
  tau <- rep(1,t)
  psi <- rep(0,t)
  tau[1]     = cor(data,method="kendall")[1,2]
  for (i in 2:t){
    if(i <= 10){ #parece que no esta como en Patton, la funcion Lambda debe estar aqui en RHS como en la cop.t dinamica de arriva
      psi[i] <- theta[1]+theta[2]*psi[i-1]+theta[3]*mean(abs(u[1:i-1]-v[1:i-1]))
    }
    else{
      psi[i] <- theta[1]+theta[2]*psi[i-1]+theta[3]*mean(abs(u[(i-10):i-1]-v[(i-10):i-1]))
    }
    tau[i] <- 0.0001+0.75/(1+exp(-psi[i]))
  }
  psi <- 1/(1-tau)
  return(psi)
}
 
 # Optimization algorithm
 data=u.s
 out = nloptr(x0 = c(0.1, -0.3, -0.5), eval_f = GumbelTVLogL, lb = c(-10, -10, -10), ub = c(10, 10, 10),
             opts = list(algorithm = "NLOPT_LN_COBYLA", xtol_rel = 1e-05, maxeval=10000), data = data)
 theta = GumbelTV(out$solution, data = data)
 sol = out$solution
 aic = 2 * length(sol) - 2 * (-GumbelTVLogL(sol, data))
 se = diag(sqrt(abs(solve(optimHess(sol, GumbelTVLogL, data = data)))))
}

#---- Ejemplos Lambda function ----#
if(0){
  #kappa[jj] = 1.998/(1 + exp(-psi1[jj])) - 0.999
  #kappa[jj] = (1 - exp(-psi1[jj])) / (1 + exp(-psi1[jj])) #LFM
  #kappa[jj, 1] = 2*(1/(1 + exp(-psi1[jj,1])) ) + 1 #in [1,3]
  #kappa[jj, 2] = 3*(1/(1 + exp(-psi1[jj,2])) ) + 0  #in [0,3]
 Lambda <- function(X, a,b,c){    #Max: <c> + <b> #Min: <b>
  f = c*(a/(1 + exp(-X)) ) + b
  #f  = a/(1 + exp(-X)) - b
  #f = (1 - exp(-X)) / (1 + exp(-X)) 
  return(f)
 }

 X   = 3*rnorm(100)
 L.X = Lambda(X,a=1, b=0, c=3)
 #L.X = Lambda(X,a=0, b=0, c=0)
 #x11()
 plot(X,L.X)
}


