# 10/4/2022 ---------------------------------------------------------------
library(fANCOVA)
arma.seleccion=function(ts_object, AR.m=5, MA.m=0, d=0, drift=TRUE, metodo='ML'){
  index = 1
  data = data.frame(p=integer(), d=integer(), q=integer(), AIC=double(), BIC=double())
  for (p in 0:AR.m) {
    for (q in 0:MA.m)  {
      fitp <- try(arima(ts_object, order=c(p, d, q), include.mean=drift,
                    method=metodo, optim.control=list(maxit=1000)))
      if (class(fitp)=='try-error') fitp = arima(ts_object, order=c(p, d, q), include.mean=drift,
                                         method="CSS-ML", optim.control=list(maxit=1000))
      data[index,] = c(p, d, q, AIC(fitp), BIC(fitp))
      index = index + 1
    }
  }  
  Criterios     = matrix(NA, nrow=2, ncol=3, dimnames=list(c('AIC','BIC'), c('p','q', 'Min.IC')))
  model.min.AIC = which.min(data[,'AIC'])
  Criterios['AIC',] =  as.matrix(data[model.min.AIC, c('p','q','AIC')])
  model.min.BIC = which.min(data[,'BIC'])
  Criterios['BIC',] =  as.matrix(data[model.min.BIC, c('p','q','BIC')])
  return(Criterios)
}
if (0) {
  N=100
  r.boot = 10
  set.seed(131313)
  y   = as.matrix(rnorm(N)) 
  x   = as.matrix(rnorm(N)) 
  
  ecdf.VaR   = ecdf(VaR)
  ecdf.CoVaR = ecdf(CoVaR)
  Ecdf.VaR   = ecdf.VaR(seq(min(VaR),max(VaR),length.out=N))
  Ecdf.CoVaR = ecdf.CoVaR(seq(min(CoVaR),max(CoVaR),length.out=N))
  ks         = max(abs(Ecdf.VaR - Ecdf.CoVaR))  * (sqrt((N^2)/(2*N)))
  
  p = 5 
  data.VaR = embed(VaR,p+1)
  AR.VaR   = lm(data.VaR[,1]~data.VaR[,-1])
  res.boot = wild.boot(AR.VaR$res, nboot=r.boot)
  VaR.boot = matrix(rep(AR.VaR$fit, time=r.boot), ncol=r.boot) + res.boot
  
  data.CoVaR = embed(CoVaR,p+1)
  AR.CoVaR   = lm(data.CoVaR[,1]~data.CoVaR[,-1])
  res.boot   = wild.boot(AR.CoVaR$res, nboot=r.boot)
  CoVaR.boot = matrix(rep(AR.CoVaR$fit, time=r.boot), ncol=r.boot) + res.boot
  
  edcf.Col <- function(x, y, x.boot, y.boot){
    n.rep = ncol(VaR.boot)
    n     = nrow(VaR.boot)
    KS    = matrix(NA,n.rep,1)
    for (i in 1:n.rep){
      ecdf.boot.x   = ecdf(x.boot[,i])
      Ecdf.boot.x   = ecdf.boot.x(seq(min(x),max(x),length.out=n))
      ecdf.boot.y   = ecdf(y.boot[,i])
      Ecdf.boot.y   = ecdf.boot.y(seq(min(y),max(y),length.out=n))
      KS[i]         = max(abs(Ecdf.boot.x - Ecdf.boot.y))  * (sqrt((n^2)/(2*n)))
    }
    return(KS)
  }
  
  KS = edcf.Col(VaR, CoVaR, VaR.boot, CoVaR.boot)
  
  
  plot(ecdf.VaR)
  plot(Ecdf.VaR)
}
# --------------------------- KS.bootstrap --------------------------------
# Función que permite realizar la prueba de Kolmogorov-Smirnov con wild
# bootstrap para evaluar la hiótesis nula de igualdad entre el CoVaR y el VaR
# x: CoVaR
# y: VaR
# --------------------------- Argumentos de entrada -----------------------
# <x>: Vector de CoVaR estimado.
# <y>: Vector de VaR estimado.
# <n.boot>: Número de repeticiones para el bootstrapping. 
# <density.plot>: <TRUE> grafica la densidad de los estadísticos KS calculados.
#
# --------------------------- Argumentos de salida ------------------------
# <KS.1>: Valor del estadístico de KS para la muestra original. 
# <KS>: Vector con los valores del etadístico KS para las muestras simuladas.
# <p.value>: Valor P de la prueba asociado a H0: x=y 
#
# -------------------------------------------------------------------------
KS.bootstrap=function(x, y, n.boot=500, density.plot=TRUE, alternative=c("two.sided", "less", "greater")[3]){
  
  # Assumption of nrow(x)=nrow(y) 
  # T statistic for original series.
  if (0) {
    ecdf.x        = ecdf(x)
    Ecdf.x        = ecdf.x(seq(min(x,y),max(x,y), length.out=nrow(x)))
    ecdf.y        = ecdf(y)
    Ecdf.y        = ecdf.y(seq(min(x,y), max(x,y), length.out=nrow(y)))
    KS.1          = max(abs(Ecdf.x - Ecdf.y))*(sqrt((nrow(x)^2)/(2*nrow(x)))) # Estadístico de la serie original.
  }
  KS.1          = ks.test(x=x,y=y, alternative=alternative)
  KS.1          = KS.1$statistic    
  # AR order for series
  p.x           = arma.seleccion(x, AR.m=5, MA.m=0)
  p.y           = arma.seleccion(y, AR.m=5, MA.m=0)
  p             = max(p.x['AIC','p'], p.y['AIC','p']); if(p==0) p=1
  
  # Resampling for AR.x residuals
  data.x        = embed(x, p+1)
  AR.x          = lm(data.x[,1]~data.x[,-1])
  res.boot.x    = wild.boot(AR.x$res, nboot=n.boot)  # Residuals resampling
  Coef.AR.x     = coefficients(AR.x)
  x.boot        = matrix(NA, nrow(res.boot.x), n.boot)
  x.boot        = rbind(as.matrix(x[1:p])%x%matrix(1,1,n.boot), x.boot)
  
  # resampling for AR.y residuals
  data.y        = embed(y, p+1)
  AR.y          = lm(data.y[,1]~data.y[,-1])
  res.boot.y    = wild.boot(AR.y$res, nboot=n.boot)  # Residuals resampling
  Coef.AR.y     = coefficients(AR.y)
  y.boot        = matrix(NA, nrow(res.boot.y), n.boot)
  y.boot        = rbind(as.matrix(y[1:p])%x%matrix(1,1,n.boot), y.boot)
  
  # resampling of time series (x and y)
  for (t in (p+1):nrow(x)){
    y.boot[t,]   = Coef.AR.y[1] + Coef.AR.y[-1][p:1]%*%y.boot[(t-p):(t-1),] + res.boot.y[t-p,]
    x.boot[t,]   = Coef.AR.x[1] + Coef.AR.x[-1][p:1]%*%x.boot[(t-p):(t-1),] + res.boot.x[t-p,]
  }  
  x.boot   = x.boot[-(1:p),]
  y.boot   = y.boot[-(1:p),]
  
  # KS's Wild bootstrap 
  KS    = matrix(NA,n.boot,1)
  for (i in 1:n.boot){
    if (0) {
      ecdf.boot.x   = ecdf(x.boot[,i])
      Ecdf.boot.x   = ecdf.boot.x(seq(min(x.boot[,i],y.boot[,i]), max(x.boot[,i],y.boot[,i]), length.out=nrow(x)))
      ecdf.boot.y   = ecdf(y.boot[,i])
      Ecdf.boot.y   = ecdf.boot.y(seq(min(x.boot[,i],y.boot[,i]), max(x.boot[,i],y.boot[,i]), length.out=nrow(x)))
      KS[i]         = max(abs(Ecdf.boot.x - Ecdf.boot.y))*(sqrt((nrow(x)^2)/(2*nrow(x))))
    } # Cambiar por función!!!
    test  = ks.test(x=x.boot[,i],y=y.boot[,i], alternative=alternative)
    KS[i] = test$statistic
  }
  # P-value 
  #p.value          = (sum(KS>KS.1)/n.boot)
  p.value           = mean(KS>KS.1)
  cat('p-value: ', p.value,'\n')
  if (density.plot==TRUE) {
    d=density(KS)
    x11()
    plot(d, main="Density")
    polygon(d, col="steelblue", border="blue")
  }
  return(list(KS.1=KS.1, KS=KS, p.value=p.value))
}

# # Prueba con series de VaR & CoVaR
# equal=as.matrix(rnorm(200))
# KS.Col=KS.bootstrap(x=equal, y=equal, n.boot=300)
# #Estructura para objeto CoVaR_data

if (0){
  n.rep    = 1000        
  Object   = CoVaR_pre_post
  x        = Object$CoVaR
  y        = Object$VaR
  P.values = matrix(NA, nrow=ncol(x), ncol=3, dimnames=list(colnames(x), c('VaR average', 'CoVaR average','P-value')))
  for (i in 1:ncol(x)) {
    Bootstrap.KS = KS.bootstrapped(x=as.matrix(x[,i]), y=as.matrix(y[,i]), boot=n.rep, density.plot=F)
    P.values[i,'P-value']       = Bootstrap.KS$p.value
    P.values[i,'CoVaR average'] = mean(x[,i])
    P.values[i,'VaR average']   = mean(y[,i])
  }
  print(P.values)
}
