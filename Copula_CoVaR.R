#-----
#-------------------------------------------------------------------------------------#
# Objetivo: Usando la función arma.selection crea una matriz que contiene los mejores
#           ordenes p y q para cada serie en la matriz. 
#-------------------------------------------------------------------------------------#
# Argumentos:
# <Datos>: Matriz de datos con las series de tiempo a las cuales se le quiere 
#          determinar el orden p,q óptimo.
# <AR.m>: Orden autoregresivo máximo a tener en cuenta.
# <MA.m>: Orden de media móvil máximo a tener en cuenta.
# <d>:    Orden de integración
# <drift>: TRUE incluye término de deriva en el modelo, FALSE lo excluye.
#-------------------------------------------------------------------------------------#
# Salida: matriz que contiene los mejores ordenes p y q para cada serie .
#-------------------------------------------------------------------------------------#
#-----
ARMA.ORDER.DF = function(Datos, AR.max=5, MA.max=0, d=0, drift=T){
  N = length(colnames(Datos))
  ARMA.Order = matrix(0, N, 2, dimnames=list(colnames(Datos), c('p','q')))
  for (i in 1:N) 
    ARMA.Order[i,] = (arma.seleccion(Datos[,i], AR.m=AR.max, MA.m=MA.max, d=d, drift=drift))['BIC',c('p','q')]  #<--- <AR.max>, <MA.max>, <drift>
  return(ARMA.Order)
}

#-----
#-------------------------------------------------------------------------------------#
# Objetivo: Calcular el CoVaR mediante el uso de copulas bajo la metodología de Reboredo & Ugolini
#           ó Girardi & Ergün (Ver manual).
#-------------------------------------------------------------------------------------#
# Argumentos:
# <Serie.1>: Serie uno del CoVaR. e.g Países.
# <Serie.2>: Serie dos del CoVaR, aquella que va a ser estresada, e.g Petróleo  
# <Datos>:   Matriz de datos que contiene las series. 
# <ALPHA>:   Nivel de significancia asociada al CoVaR.
# <BETA>:    Nivel de significancia asociada al  VaR (Serie estresada).
# <type>:    Identificador de mejor copula asociada al par de datos a modelar.
#            Se utiliza el objeto resultante de <Copula_seleccion>
# <forecast.type>: Tipo de pronóstico a llevar a cabo en la estimación del GARCH para obtener las pseudo-observaciones:
#                  1. <rolling>: Realiza un pronóstico rolling, cuya especificación dependerá
#                                de los parámetros <window.size>, <refit.every> y <n.ahead>.
#                  2. <in sample>: Realiza un pronóstico dentro de muestra extraído directamente
#                                  del modelo GARCH.
# <COND>: Metodología del CoVaR a calcular (Liu et al. ó Girardi & Ergün).
#         1- 'Less':  Utiliza la metodología de Girardi & Ergün
#         2- 'Equal': Utiliza la metodología de Liu et al.
#        Ver manual para las diferencias metodológicas. 
# <window.size>: Determina el tamaño de las ventanas de estimación cuando <forecast.type='rolling'>
# <refit.every>: Especifica cada cuantas ventanas (pronósticos) se desean reestimar los parámetros del modelo.
# <n.ahead>: Define el horizonte de pronóstico en cada ventana cuando <forecast.type='rolling'>, por defecto 1.
# <AR.max>:  Orden autoregresivo máximo a tener en cuenta para los modelos de la media.
# <MA.max>:  Orden de media móvil máximo a tener en cuenta para los modelos de la media.
# <garch.series.2.first>: Objeto con la primera estimación del GARCH de la Serie.2, 
#                         puede ser tipo 'rolling' o 'in sample'. 
#                         -Si es 'rolling', el objeto debe ser una salida de la función 'garch_roll'.
#                         -Si es 'in sample', el objeto debe ser una salida de la función 'ugarchfit'.
#<external.regressors>: Por defecto NULL. Si no es NULL se brindan los nombres de los regresores 
#                        externos que se incluyen en la ecuación de la media del modelo.                 
#-------------------------------------------------------------------------------------#
# Salida: Lista que contiene las series estimadas para el CoVaR, VaR, CoVaRUp, VaRUp, DeltaCoVaR, 
#         DeltaCoVaRUp y CoVaRMedian (El CoVaR con la variable estresada en su mediana) para el par de series especificado. 
#-------------------------------------------------------------------------------------#
#-----
CoVaR.Copula <- function(Serie.1='CDS_Brazil', Serie.2='WTI', Datos=G20_Crisis_WNA, ALPHA=0.05,
                         BETA=0.05, type='Dyn.Rot.Gumbel', forecast.type='rolling',
                         COND='Equal', window.size=150, refit.every=20, n.ahead=1, 
                         AR.max=5, MA.max=0, garch.series.2.first=NULL, 
                         external.regressors=NULL, ARMA.Order=NULL){
  #--- Identificacion del modelo ARMA de la <Serie.1> ---#
  if (is.null(ARMA.Order)==TRUE) {
    ARMA.Order = ARMA.ORDER.DF(Datos=Datos[,Serie.1], AR.max = AR.max, MA.max = MA.max) #Usando valores por defecto
  }
  #if (Serie.1=='CDS_Korea') ARMA.Order[,'p']=5
  #--- Estimacion del modelo GARCH de la <Serie.1> ---# 
  #--- Supestos: Modelo <gjrGARCH(1,1)> y distribucion <sstd>
  #--- 1step-ahead Rolling forecasts ----#
  if(forecast.type == 'rolling'){
    cat(paste0('AR order: ',ARMA.Order[Serie.1,1],'\nMA order: ',ARMA.Order[Serie.1,2],'\n' ))
    cat('Rolling forecast started...\n')
    forecast.Serie.1       = garch.roll(Data=Datos[,Serie.1], model="gjrGARCH", garchOrder=c(1,1), armaOrder=ARMA.Order[Serie.1,], 
                                        distribution.model="sstd", window.size=window.size, refit.every=refit.every, 
                                        show.progress=TRUE, n.ahead=n.ahead, external.regressors = external.regressors)
    forecast.Serie.1.mean  = forecast.Serie.1$mean
    forecast.Serie.1.sigma = forecast.Serie.1$sigma
    forecast.Serie.1.skew  = forecast.Serie.1$skew
    forecast.Serie.1.shape = forecast.Serie.1$shape
    VaR   = as.xts(forecast.Serie.1.mean + forecast.Serie.1.sigma*(rugarch::qdist(distribution='sstd', p=ALPHA,   shape=forecast.Serie.1.shape, skew=forecast.Serie.1.skew)))
    VaRUp = as.xts(forecast.Serie.1.mean + forecast.Serie.1.sigma*(rugarch::qdist(distribution='sstd', p=1-ALPHA, shape=forecast.Serie.1.shape, skew=forecast.Serie.1.skew)))
  }
  #--- In sample forecasts (for ARMA(0,0) models ) ----#
  if(forecast.type == 'in sample'){
    model.Serie.1     = ugarchspec(variance.model=list(model="gjrGARCH", garchOrder=c(1,1)), 
                                   mean.model=list(armaOrder=ARMA.Order[Serie.1,]), distribution.model="sstd")
    model.Serie.1.fit = ugarchfit(spec=model.Serie.1, data=Datos[,Serie.1])
    if (ARMA.Order[1] == 0 &  ARMA.Order[2] == 0) {
      forecast.Serie.1.mean  =  as.matrix(rep(coef(model.Serie.1.fit)['mu'], nrow(Datos),)) 
      forecast.Serie.1.sigma =  as.xts(sigma(model.Serie.1.fit))  
    } else {
      forecast.Serie.1.mean  =  Datos[,Serie.1] - residuals(model.Serie.1.fit)
      forecast.Serie.1.sigma =  as.xts(sigma(model.Serie.1.fit))  
    }
    
    VaR   = as.xts(forecast.Serie.1.mean + forecast.Serie.1.sigma*(rugarch::qdist(distribution='sstd', p=ALPHA,   shape=coef(model.Serie.1.fit)['shape'],skew=coef(model.Serie.1.fit)['skew'])))
    VaRUp = as.xts(forecast.Serie.1.mean + forecast.Serie.1.sigma*(rugarch::qdist(distribution='sstd', p=1-ALPHA, shape=coef(model.Serie.1.fit)['shape'],skew=coef(model.Serie.1.fit)['skew'])))
  }
  
  #--- <Serie.2> model ('WTI')---#
  #--- Identificacion del modelo ARMA de la <Serie.2> ---#
  if (is.null(ARMA.Order)==TRUE)ARMA.Order = ARMA.ORDER.DF(Datos=Datos[,Serie.2]) 
  #--- Estimacion rolling ---#
  if (forecast.type=='rolling'){
    if (is.null(garch.series.2.first)==TRUE) {
      #--- Modelo ARMA-GARCH. Sup: GARCH(1,1) con errores distribuidos <sstd>
      cat(paste0('AR order: ',ARMA.Order[Serie.2,1],'\nMA order: ',ARMA.Order[Serie.2,2],'\n' ))
      cat('Rolling forecast (Serie.2) started...\n')
      forecast.Serie.2       = garch.roll(Data=Datos[,Serie.2], model="gjrGARCH", garchOrder=c(1,1), armaOrder=ARMA.Order[Serie.2,], 
                                          distribution.model="sstd", window.size=window.size, refit.every=refit.every, 
                                          show.progress=TRUE, n.ahead=n.ahead, external.regressors=NULL)
    } else{
       cat('Using garch.series.2.first input for rolling forecast of Serie.2\n')
       #Cuando hay objeto previo
       forecast.Serie.2       = garch.series.2.first
    }
  } 
  #--- Estimacion estandar ---#
  if (forecast.type=='in sample'){
    if (is.null(garch.series.2.first)==TRUE) {
      #--- Modelo ARMA-GARCH. Sup: GARCH(1,1) con errores distribuidos <sstd>
       model.Serie.2     = ugarchspec(variance.model = list(model="sGARCH", garchOrder = c(1,1)), 
                                   mean.model = list(ARMA.Order[Serie.2,]), 
                                   distribution.model = "sstd")
       model.Serie.2.fit = ugarchfit(spec = model.Serie.2, data = Datos[,Serie.2])
    } else{
       cat('Using garch.series.2.first input for rolling forecast of Serie.2\n')
       #Cuando hay objeto previo
       model.Serie.2.fit = garch.series.2.first
    }
  }
  
  #--- Copula(Serie.1, Serie.2) (i.e. Country,Brend) ----#
  #Pseudo-observaciones de la Serie.1
  if (forecast.type=='rolling') {
    sd_errors1 = (forecast.Serie.1$residuals)/(forecast.Serie.1$sigma)
    u1         = rugarch::pdist(distribution='sstd', q=sd_errors1, skew=forecast.Serie.1.skew, shape=forecast.Serie.1.shape)
  }
  if (forecast.type=='in sample') {
    sd_errors1 = residuals(model.Serie.1.fit)/sigma(model.Serie.1.fit)
    u1         = rugarch::pdist(distribution='sstd', q=sd_errors1, skew=coef(model.Serie.1.fit)['skew'], shape=coef(model.Serie.1.fit)['shape'])
  }
  #Pseudo-observaciones de la Serie.2
  if (forecast.type == 'rolling') {
      sd_errors2 = (forecast.Serie.2$residuals)/(forecast.Serie.2$sigma)
      u2         = rugarch::pdist(distribution='sstd', q=sd_errors2, skew=forecast.Serie.2$skew, shape=forecast.Serie.2$shape)
  } 
  if (forecast.type == 'in sample'){#--- Identificacion del modelo ARMA de la <Serie.2> ---#
    sd_errors2 = residuals(model.Serie.2.fit)/sigma(model.Serie.2.fit)
    u2         = rugarch::pdist(distribution='sstd', q=sd_errors2, skew=coef(model.Serie.2.fit)['skew'], shape=coef(model.Serie.2.fit)['shape'])
  }
  
  u.s        = cbind(u1, u2)
  
  cat('Copula fitting started...\n')
  
  if (type=='Gaussian' | type=='Gaussian.Up'){
    fit.Copula = BiCopEst(u1,u2, family=1, method='mle')
    par.Cop.1  = fit.Copula$par
    par.Cop.2  = NULL
    print(fit.Copula)
  }
  if (type=='Gumbel' | type=='Gumbel.Up'){
    fit.Copula = BiCopEst(u1,u2, family=4, method='mle')
    par.Cop.1  = fit.Copula$par
    par.Cop.2  = NULL
    print(fit.Copula)
  }
  if (type=='Joe' | type=='Joe.Up'){
    fit.Copula = BiCopEst(u1,u2, family=6, method='mle')
    par.Cop.1  = fit.Copula$par
    par.Cop.2  = NULL
    print(fit.Copula)
  }
  if (type=='Frank' | type=='Frank.Up'){
    fit.Copula = BiCopEst(u1,u2, family=5, method='mle')
    par.Cop.1  = fit.Copula$par
    par.Cop.2  = NULL
    print(fit.Copula)
  }
  if (type=='Rot.Gumbel' | type=='Rot.Gumbel.Up'){
    fit.Copula = BiCopEst(u1,u2, family=14, method='mle') #rotated Gumbel copula (180 degrees)
    #fit.Copula = BiCopEst(u1,u2, family=24, method='mle') #rotated Gumbel copula (90 degrees)
    #fit.Copula = BiCopEst(u1,u2, family=34, method='mle') #rotated Gumbel copula (270 degrees)
    par.Cop.1  = fit.Copula$par
    par.Cop.2  = NULL
    print(fit.Copula)
  }
  if (type=='AMH' | type=='AMH.Up'){
    cop_model  = amhCopula(dim=2)
    fit.Copula = fitCopula(cop_model, u.s , method='ml')
    par.Cop.1  = fit.Copula@estimate
    par.Cop.2  = NULL
    print(fit.Copula)
  }
  if (type=='Plackett' | type=='Plackett.Up'){
    cop_model  = plackettCopula()
    fit.Copula = fitCopula(cop_model, u.s , method='ml')
    par.Cop.1  = fit.Copula@estimate
    par.Cop.2  = NULL
    print(fit.Copula)
  }
  if (type=='BB7' | type=='BB7.Up'){
    fit.Copula = BiCopEst(u1,u2, family=9, method='mle')
    par.Cop.1  = fit.Copula$par
    par.Cop.2  = fit.Copula$par2 
    print(fit.Copula)
  }
  if (type=='Clayton' | type=='Clayton.Up'){
    fit.Copula = BiCopEst(u1,u2, family=3, method='mle')
    par.Cop.1  = fit.Copula$par
    par.Cop.2  = NULL 
    print(fit.Copula)
  }
  if (type=='Student' | type=='Student.Up'){
    fit.Copula = BiCopEst(u1,u2, family=2, method='mle')
    par.Cop.1  = fit.Copula$par
    par.Cop.2  = fit.Copula$par2
    print(fit.Copula)
  }
  if (type=='Dyn.Student' | type=='Dyn.Student.Up'){
    fit.Copula = dynamicT.LF(data=u.s, plot=FALSE, printout=FALSE)
    par.Cop.1  = fit.Copula$tvtpdep; par.Cop.1[1:10] = mean(fit.Copula$tvtpdep[-c(1:10)])
    par.Cop.2  = fit.Copula$result['nu','coef']
    print(fit.Copula$result)
  }
  if (type=='Dyn.Gaussian' | type=='Dyn.Gaussian.Up'){
    fit.Copula = dynamicnormal.LF(data=u.s,plot=FALSE, printout=FALSE) 
    par.Cop.1  = fit.Copula$tvtpdep; par.Cop.1[1:10] = mean(fit.Copula$tvtpdep[-c(1:10)])
    par.Cop.2  = NULL  
    print(fit.Copula$result)
  }  
  if (type=='Dyn.Gumbel' | type=='Dyn.Gumbel.Up'){
    fit.Copula = dynamicGum.LF(data=u.s,plot=FALSE, rotated=FALSE, printout=FALSE) 
    par.Cop.1  = fit.Copula$tvtpdep; par.Cop.1[1:10] = mean(fit.Copula$tvtpdep[-c(1:10)])
    par.Cop.2  = NULL  
    print(fit.Copula$result)
  } 
  if (type=='Dyn.Rot.Gumbel' | type=='Dyn.Rot.Gumbel.Up'){
    fit.Copula = dynamicGum.LF(data=u.s,plot=FALSE, rotated=TRUE, printout=FALSE) 
    par.Cop.1  = fit.Copula$tvtpdep; par.Cop.1[1:10] = mean(fit.Copula$tvtpdep[-c(1:10)])
    par.Cop.2  = NULL  
    print(fit.Copula$result)
  } 
  if (type=='Dyn.BB7' | type=='Dyn.BB7.Up'){
    fit.Copula = dynamicBB7.LF(data=u.s,plot=FALSE, printout=FALSE) 
    par.Cop.1  = fit.Copula$par1.t; par.Cop.1[1:10] = mean(fit.Copula$par1.t[-c(1:10)])
    par.Cop.2  = fit.Copula$par2.t; par.Cop.2[1:10] = mean(fit.Copula$par2.t[-c(1:10)])  
    print(fit.Copula$result)
  } 
  # ALPHA = 0.05 #<--- Serie.1 (Country) Stock CoVaR quantile
  # BETA  = 0.05 #<--- Oil VaR-quantile
  if (forecast.type == 'rolling') {
    CoVaR.Cop        = CoVaR.LF(y=ALPHA,x=BETA, par=par.Cop.1, par2=par.Cop.2,
                                dof=forecast.Serie.1.shape, gamma=forecast.Serie.1.skew, 
                                cond.mean=forecast.Serie.1.mean, cond.sigma=forecast.Serie.1.sigma, dist="tskew", type=type, COND=COND)  
    CoVaR.Cop.median = CoVaR.LF(y=ALPHA,x=0.5, par=par.Cop.1, par2=par.Cop.2,
                                dof=forecast.Serie.1.shape, gamma=forecast.Serie.1.skew, 
                                cond.mean=forecast.Serie.1.mean, cond.sigma=forecast.Serie.1.sigma, dist="tskew", type=type, COND=COND)
  }else {
    CoVaR.Cop        = CoVaR.LF(y=ALPHA, x=BETA, par=par.Cop.1, par2=par.Cop.2,
                                dof=coef(model.Serie.1.fit)['shape'], gamma=coef(model.Serie.1.fit)['skew'], 
                                cond.mean=forecast.Serie.1.mean, cond.sigma=forecast.Serie.1.sigma, dist="tskew", type=type, COND=COND)
    CoVaR.Cop.median = CoVaR.LF(y=ALPHA, x=0.5, par=par.Cop.1, par2=par.Cop.2,
                                dof=coef(model.Serie.1.fit)['shape'], gamma=coef(model.Serie.1.fit)['skew'], 
                                cond.mean=forecast.Serie.1.mean, cond.sigma=forecast.Serie.1.sigma, dist="tskew", type=type, COND=COND)  
  }
  DeltaCoVaRDown  = as.xts(CoVaR.Cop$CoVaR)    - as.xts(CoVaR.Cop.median$CoVaR)
  DeltaCoVaRUp    = as.xts(CoVaR.Cop$CoVaR.up) - as.xts(CoVaR.Cop.median$CoVaR.up)
  CoVaRDown       = as.xts(CoVaR.Cop$CoVaR)
  CoVaRUp         = as.xts(CoVaR.Cop$CoVaR.up)
  CoVaRMedian     = as.xts(CoVaR.Cop.median$CoVaR)
  if (sum(is.na(CoVaRDown))!= 0)  CoVaRDown = na.fill(CoVaRDown,'extend')
  if (sum(is.na(CoVaRUp))!= 0)    CoVaRUp   = na.fill(CoVaRUp,  'extend')
  if (sum(is.na(VaR))!= 0)        VaR       = na.fill(VaR,      'extend')
  if (sum(is.na(VaRUp))!= 0)      VaRUp     = na.fill(VaRUp,    'extend')
  if (sum(is.na(DeltaCoVaRDown))!= 0) DeltaCoVaRDown = na.fill(DeltaCoVaRDown, 'extend')
  if (sum(is.na(DeltaCoVaRUp))!= 0)   DeltaCoVaRUp   = na.fill(DeltaCoVaRUp, 'extend')
  if (sum(is.na(DeltaCoVaRUp))!= 0)   CoVaRMedian    = na.fill(CoVaRMedian,  'extend')
  
  if (forecast.type=='in sample') 
     return(list(CoVaR=CoVaRDown, VaR=VaR, CoVaRUp=CoVaRUp, VaRUp=VaRUp, DeltaCoVaR=DeltaCoVaRDown, 
                DeltaCoVaRUp=DeltaCoVaRUp, CoVaRMedian=CoVaRMedian, 
                garch.serie.2=model.Serie.2.fit, Par1=par.Cop.1, Par2=par.Cop.2))
  if (forecast.type=='rolling')
     return(list(CoVaR=CoVaRDown, VaR=VaR, CoVaRUp=CoVaRUp, VaRUp=VaRUp, DeltaCoVaR=DeltaCoVaRDown, 
                  DeltaCoVaRUp=DeltaCoVaRUp, CoVaRMedian=CoVaRMedian, 
                 garch.serie.2=forecast.Serie.2, Par1=par.Cop.1, Par2=par.Cop.2))
}


#-----
#-------------------------------------------------------------------------------------#
# Objetivo: Obtener y graficar un objeto que contenga las series estimadas del CoVaR, VaR, CoVaRUp, VaRUp, DeltaCoVaR, 
#            DeltaCoVaRUp y CoVaRMedian para todas las series disponibles.
#            Esto mediante el uso de la función Copula_CoVar que realiza el procedimiento para pares de series. 
#-------------------------------------------------------------------------------------#
# Argumentos:
# <Data>: Matriz de datos completa, con todas las series a modelar. 
# <Serie.1>: Serie a la cual se le va a a calcular el CoVaR. Debe indiacarse el nombre de
#            la columna que contiene la serie (Vector con el nombre de las columnas para hacerlo con todas las posibles).
#            NO se debe incluir la serie a estresar en este argumento.
# <Serie.2>: La serie a estresar, se debe indicar el nombre de la colimna correspondiente.
# <copulas>: Vector que contiene las mejores copulas para cada par de datos siguiendo el mismo ordenamiento de <Serie.1>.
#            En este argumento se debe poner el objeto resultante de la sección "Selección de copula"
# <alpha>:Nivel del CoVaR
# <beta>: Nivel del VaR (cuantil de la variable estresada) 
# <plot>: Lógica. TRUE imprime en un gráfico para cada serie todas las medidas de riesgo mencionadas en <Objetivo>
# <COND>: Metodología del CoVaR a calcular (Liu et al. ó Girardi & Ergün).
#         1- 'Less':  Utiliza la metodología de Girardi & Ergün
#         2- 'Equal': Utiliza la metodología de Liu et al.
#         Ver manual para las diferencias metodológicas. 
# <forecast.type>: Tipo de pronóstico a llevar a cabo en la estimación del GARCH para obtener las pseudo-observaciones:
#                  1. <rolling>: Realiza un pronóstico rolling, cuya especificación dependerá
#                                de los parámetros <window.size>, <refit.every> y <n.ahead>.
#                  2. <in sample>: Realiza un pronóstico dentro de muestra extraído directamente
#                                  del modelo GARCH.
# <window.size>: Determina el tamaño de las ventanas de estimación cuando <forecast.type='rolling'>
# <refit.every>: Especifica cada cuantas ventanas (pronósticos) se desean reestimar los parámetros del modelo.
# <n.ahead>:     Define el horizonte de pronóstico en cada ventana cuando <forecast.type='rolling'>, por defecto 1.
# <external.regressors>: Por defecto NULL. Sí se desea incluir regresores externos, deben introducirse en formato
#                        formato 'xts' previamente rezagados. 
#-------------------------------------------------------------------------------------#
# Salida: Lista que contiene en cada elemento matrices con las series estimadas del CoVaR, VaR, CoVaRUp, VaRUp, DeltaCoVaR, 
#         DeltaCoVaRUp y CoVaRMedian para el par todas las series disponibles.
#-------------------------------------------------------------------------------------#
#-----
CoVaR_DF = function(Data, Serie.1=Name.Series1, Serie.2='WTI', copulas, alpha, beta,
                    plot=TRUE,COND='Equal', forecast.type='in sample', window.size=150,
                    refit.every=20, n.ahead=1, external.regressors=NULL, ARMA.Order=NULL){
  
  time.ini                      = Sys.time() 
  #--- Inicializacion de los objetos de salida de la funcion ---#
  CoVaR.series                  = xts(matrix(0, nrow=nrow(Data[(window.size+1):nrow(Data),]) , 
                                             ncol=length(Serie.1)), order.by=ichimoku::index(Data[(window.size+1):nrow(Data),])) 
  colnames(CoVaR.series)        = Serie.1
  VaR.series  = CoVaRUp.series  = VaRUp.series = DeltaCoVaR.series = DeltaCoVaRUp.series = CoVaRMedian.series = CoVaR.series 
  
  i0 = 0 
  for (i in Serie.1){
    #print(i)
    #print(colnames(Data)[i], quote = FALSE)
    i0 = i0 + 1
    cat('\n',i0,'-',i,'\n')
    if(i==Serie.1[1]){ 
      CoVaR.1=CoVaR.Copula(Serie.1=i, Serie.2=Serie.2, Datos=Data, ALPHA=alpha,
                           BETA=beta, type=copulas[i], COND=COND, 
                           forecast.type=forecast.type, refit.every=refit.every,
                           window.size=window.size, n.ahead=n.ahead, 
                           external.regressors=external.regressors,
                           ARMA.Order=ARMA.Order)
      garch.model.2 = CoVaR.1$garch.serie.2 # Guardado del primer GARCH de <Serie.2>
    } else{ # Se usa el GARCH de <Serie.2> de la linea anterior
      CoVaR.1=CoVaR.Copula(Serie.1=i, Serie.2=Serie.2, Datos=Data, ALPHA=alpha,
                           BETA=beta, type=copulas[i], COND=COND, 
                           forecast.type=forecast.type, refit.every=refit.every,
                           window.size=window.size, n.ahead=n.ahead, 
                           garch.series.2.first=garch.model.2,
                           external.regressors=external.regressors, 
                           ARMA.Order=ARMA.Order)
    }  
    CoVaR.series[,i]        = xts(CoVaR.1$CoVaR)
    CoVaRUp.series[,i]      = xts(CoVaR.1$CoVaRUp) 
    VaR.series[,i]          = xts(CoVaR.1$VaR)
    VaRUp.series[,i]        = xts(CoVaR.1$VaRUp)
    DeltaCoVaR.series[,i]   = xts(CoVaR.1$DeltaCoVaR)
    DeltaCoVaRUp.series[,i] = xts(CoVaR.1$DeltaCoVaRUp)
    CoVaRMedian.series[,i]  = xts(CoVaR.1$CoVaRMedian)    
    if (plot==TRUE) {
      pdf(file = paste0(Resultados,'/Graficas_CoVaR_',i,'.pdf'), onefile=FALSE)
      # print(plot.xts(CoVaR.series[,i],type="l",col="red", ylim=range(CoVaR.1),xlab="Time",
      #                ylab="", lwd=3, main=i, format.labels="%m-%Y", major.ticks = 'months', 
      #                yaxis.left=TRUE, yaxis.right=TRUE, lty='dotted'))
      # En el <range> se quitan el 'CoVaRMedian' y 'garch.series.2' (obj. 7 y 8)
      print(plot.xts(CoVaR.series[,i],type="l",col="red", ylim=range(CoVaR.1[-c(7:8)]),xlab="Time",
                     ylab="", lwd=1.5, main=i, format.labels="%Y", major.ticks = 'years', 
                     yaxis.left=TRUE, yaxis.right=TRUE, lty='dotted'))
      print(lines(VaR.series[,i],         col="black", lwd=1.5,lty='dashed'))
      print(lines(CoVaRUp.series[,i],     col="red",   lwd=1.5,lty='dotted'))
      print(lines(VaRUp.series[,i],       col="black", lwd=1.5,lty='dashed'))
      print(lines(DeltaCoVaR.series[,i],  col="blue"  ,lwd=1.5,lty='solid'))
      print(lines(DeltaCoVaRUp.series[,i],col="blue",  lwd=1.5,lty='solid'))
      print(abline(h=0,        col="black", lwd=1.5))
      print(addLegend("topright", lwd=2,
                      legend.names = c('CoVaR', 'VaR',   'DeltaCoVaR'), 
                      lty = c('dotted','dashed','solid'), col = c('red',   'black', 'blue')))
      dev.off()
    }
    graphics.off()
  }
  
  #-- Time of execution --#
  time.diff = (Sys.time() - time.ini) 
  time.diff.over.60 = (Sys.time() - time.ini)/60 
  
  #-- Printout of final results --#
  print(time.diff)
  print(time.diff.over.60)
  
  # objetos de salida
  out = list(CoVaR=CoVaR.series, CoVaRUp=CoVaRUp.series, VaR=VaR.series, VaRUp=VaRUp.series,
             DeltaCoVaR=DeltaCoVaR.series, DeltaCoVaRUp=DeltaCoVaRUp.series,
             CoVaRMedian=CoVaRMedian.series)
  return(out)
}

#plotxts <- fix("plot.xts") 
# In the editor that opens, replace the lines below:
###    text.exp <- c(expression(text(xlim[1], 0.5, main, font = 2, 
###        col = theme$labels, offset = 0, cex = 1.1, pos = 4)), 
###        expression(text(xlim[2], 0.5, paste(start(xdata[xsubset]), 
###            end(xdata[xsubset]), sep = " / "), col = theme$labels, 
###            adj = c(0, 0), pos = 2)))
###    cs$add(text.exp, env = cs$Env, expr = TRUE)
# with these lines:
###    text.exp <- expression(text(xlim[1], 0.5, main, font = 2,
###        col = theme$labels, offset = 0, cex = 1.1, pos = 4))

# Finally, you need to ensure your copy's environment is the xts namespace:
#environment(plotxts) <- asNamespace("xts")

# Fix of <plot.xts> 
plotxts <-function (x, y = NULL, ..., subset = "", panels = NULL, multi.panel = FALSE, 
                    col = 1:8, up.col = NULL, dn.col = NULL, bg = "#FFFFFF", 
                    type = "l", lty = 1, lwd = 2, lend = 1, main = deparse(substitute(x)), 
                    observation.based = FALSE, ylim = NULL, yaxis.same = TRUE, 
                    yaxis.left = TRUE, yaxis.right = TRUE, major.ticks = "auto", 
                    minor.ticks = NULL, grid.ticks.on = "auto", grid.ticks.lwd = 1, 
                    grid.ticks.lty = 1, grid.col = "darkgray", labels.col = "#333333", 
                    format.labels = TRUE, grid2 = "#F5F5F5", legend.loc = NULL) 
{
  if (is.numeric(multi.panel)) {
    multi.panel <- min(NCOL(x), multi.panel)
    idx <- seq.int(1L, NCOL(x), 1L)
    chunks <- split(idx, ceiling(seq_along(idx)/multi.panel))
    if (length(lty) < ncol(x)) 
      lty <- rep(lty, length.out = ncol(x))
    if (length(lwd) < ncol(x)) 
      lwd <- rep(lwd, length.out = ncol(x))
    if (length(col) < ncol(x)) 
      col <- rep(col, length.out = ncol(x))
    if (!is.null(panels) && nchar(panels) > 0) {
      multi.panel <- FALSE
    }
    else {
      multi.panel <- TRUE
      panels <- NULL
      if (yaxis.same) 
        ylim <- range(x[subset], na.rm = TRUE)
    }
    for (i in 1:length(chunks)) {
      tmp <- chunks[[i]]
      p <- plot.xts(x = x[, tmp], y = y, ... = ..., subset = subset, 
                    panels = panels, multi.panel = multi.panel, col = col[tmp], 
                    up.col = up.col, dn.col = dn.col, bg = bg, type = type, 
                    lty = lty[tmp], lwd = lwd[tmp], lend = lend, 
                    main = main, observation.based = observation.based, 
                    ylim = ylim, yaxis.same = yaxis.same, yaxis.left = yaxis.left, 
                    yaxis.right = yaxis.right, major.ticks = major.ticks, 
                    minor.ticks = minor.ticks, grid.ticks.on = grid.ticks.on, 
                    grid.ticks.lwd = grid.ticks.lwd, grid.ticks.lty = grid.ticks.lty, 
                    grid.col = grid.col, labels.col = labels.col, 
                    format.labels = format.labels, grid2 = grid2, 
                    legend.loc = legend.loc)
      if (i < length(chunks)) 
        print(p)
    }
    return(p)
  }
  cs <- new.replot_xts()
  if (is.null(major.ticks)) {
    xs <- x[subset]
    mt <- c(years = nyears(xs), months = nmonths(xs), days = ndays(xs))
    major.ticks <- names(mt)[rev(which(mt < 30))[1]]
  }
  plot.call <- match.call(expand.dots = TRUE)
  if (isTRUE(multi.panel)) {
    if (NCOL(x) == 1) 
      cs$set_asp(3)
    else cs$set_asp(NCOL(x))
  }
  else {
    cs$set_asp(3)
  }
  cs$Env$cex <- if (hasArg("cex")) 
    eval.parent(plot.call$cex)
  else 0.6
  cs$Env$mar <- if (hasArg("mar")) 
    eval.parent(plot.call$mar)
  else c(3, 2, 0, 2)
  cs$Env$theme$up.col <- up.col
  cs$Env$theme$dn.col <- dn.col
  if (hasArg("colorset")) 
    col <- eval.parent(plot.call$colorset)
  if (length(col) < ncol(x)) 
    col <- rep(col, length.out = ncol(x))
  cs$Env$theme$col <- col
  cs$Env$theme$rylab <- yaxis.right
  cs$Env$theme$lylab <- yaxis.left
  cs$Env$theme$bg <- bg
  cs$Env$theme$grid <- grid.col
  cs$Env$theme$grid2 <- grid2
  cs$Env$theme$labels <- labels.col
  cs$Env$theme$srt <- if (hasArg("srt")) 
    eval.parent(plot.call$srt)
  else 0
  cs$Env$theme$las <- if (hasArg("las")) 
    eval.parent(plot.call$las)
  else 0
  cs$Env$theme$cex.axis <- if (hasArg("cex.axis")) 
    eval.parent(plot.call$cex.axis)
  else 0.9
  cs$Env$format.labels <- format.labels
  cs$Env$major.ticks <- if (isTRUE(major.ticks)) 
    "auto"
  else major.ticks
  cs$Env$minor.ticks <- if (isTRUE(minor.ticks)) 
    "auto"
  else minor.ticks
  cs$Env$grid.ticks.on <- if (isTRUE(grid.ticks.on)) 
    "auto"
  else grid.ticks.on
  cs$Env$grid.ticks.lwd <- grid.ticks.lwd
  cs$Env$grid.ticks.lty <- grid.ticks.lty
  cs$Env$type <- type
  if (length(lty) < ncol(x)) 
    lty <- rep(lty, length.out = ncol(x))
  if (length(lwd) < ncol(x)) 
    lwd <- rep(lwd, length.out = ncol(x))
  cs$Env$lty <- lty
  cs$Env$lwd <- lwd
  cs$Env$lend <- lend
  cs$Env$legend.loc <- legend.loc
  cs$Env$call_list <- list()
  cs$Env$call_list[[1]] <- plot.call
  cs$Env$observation.based <- observation.based
  if (is.character(x)) 
    stop("'x' must be a time-series object")
  cs$Env$xdata <- x
  cs$Env$xsubset <- subset
  cs$Env$column_names <- colnames(x)
  cs$Env$nobs <- NROW(cs$Env$xdata)
  cs$Env$main <- main
  if (cs$Env$observation.based) {
    cs$Env$xycoords <- xy.coords(1:NROW(cs$Env$xdata[subset]))
    cs$set_xlim(c(1, NROW(cs$Env$xdata[subset])))
    cs$Env$xstep <- 1
  }
  else {
    xycoords <- xy.coords(.index(cs$Env$xdata[cs$Env$xsubset]), 
                          cs$Env$xdata[cs$Env$xsubset][, 1])
    cs$Env$xycoords <- xycoords
    cs$Env$xlim <- range(xycoords$x, na.rm = TRUE)
    cs$Env$xstep <- diff(xycoords$x[1:2])
    cs$set_xlim(cs$Env$xlim)
  }
  if (is.null(ylim)) {
    if (isTRUE(multi.panel)) {
      if (yaxis.same) {
        yrange <- range(cs$Env$xdata[subset], na.rm = TRUE)
      }
      else {
        yrange <- range(cs$Env$xdata[, 1][subset], na.rm = TRUE)
      }
    }
    else {
      yrange <- range(cs$Env$xdata[subset], na.rm = TRUE)
    }
    if (yrange[1L] == yrange[2L]) {
      if (yrange[1L] == 0) {
        yrange <- yrange + c(-1, 1)
      }
      else {
        yrange <- c(0.8, 1.2) * yrange[1L]
      }
    }
    cs$set_ylim(list(structure(yrange, fixed = TRUE)))
    cs$Env$constant_ylim <- range(cs$Env$xdata[subset], na.rm = TRUE)
  }
  else {
    cs$set_ylim(list(structure(ylim, fixed = TRUE)))
    cs$Env$constant_ylim <- ylim
  }
  cs$set_frame(1, FALSE)
  if (!isNullOrFalse(grid.ticks.on)) {
    cs$add(expression(atbt <- axTicksByTime(xdata[xsubset], 
                                            ticks.on = grid.ticks.on), segments(xycoords$x[atbt], 
                                                                                get_ylim()[[2]][1], xycoords$x[atbt], get_ylim()[[2]][2], 
                                                                                col = theme$grid, lwd = grid.ticks.lwd, lty = grid.ticks.lty)), 
           clip = FALSE, expr = TRUE)
  }
  cs$add_frame(0, ylim = c(0, 1), asp = 0.5)
  cs$set_frame(1)
  cs$add(expression(if (NROW(xdata[xsubset]) < 400) {
    axis(1, at = xycoords$x, labels = FALSE, col = theme$grid2, 
         col.axis = theme$grid2, tcl = 0.3)
  }), expr = TRUE)
  if (!isNullOrFalse(major.ticks)) {
    cs$add(expression(axt <- axTicksByTime(xdata[xsubset], 
                                           ticks.on = major.ticks, format.labels = format.labels), 
                      axis(1, at = xycoords$x[axt], labels = names(axt), 
                           las = theme$las, lwd.ticks = 1.5, mgp = c(3, 
                                                                     1.5, 0), tcl = -0.4, cex.axis = theme$cex.axis, 
                           col = theme$labels, col.axis = theme$labels)), 
           expr = TRUE)
  }
  if (!isNullOrFalse(minor.ticks)) {
    cs$add(expression(axt <- axTicksByTime(xdata[xsubset], 
                                           ticks.on = minor.ticks, format.labels = format.labels), 
                      axis(1, at = xycoords$x[axt], labels = FALSE, las = theme$las, 
                           lwd.ticks = 0.75, mgp = c(3, 1.5, 0), tcl = -0.4, 
                           cex.axis = theme$cex.axis, col = theme$labels, 
                           col.axis = theme$labels)), expr = TRUE)
  }
  #  text.exp <- c(expression(text(xlim[1], 0.5, main, font = 2, 
  #               col = theme$labels, offset = 0, cex = 1.1, pos = 4)), 
  #                expression(text(xlim[2], 0.5, paste(start(xdata[xsubset]), 
  #                           end(xdata[xsubset]), sep = " / "), col = theme$labels, 
  #                           adj = c(0, 0), pos = 2)))
  #  cs$add(text.exp, env = cs$Env, expr = TRUE)
  text.exp <- expression(text(xlim[1], 0.5, main, font = 2,
                              col = theme$labels, offset = 0, cex = 1.05, pos = 4)) #Replacement line LFM
  #  text.exp.1 <- c(  expression(text(xlim[2], 0.5, paste(start(xdata[xsubset]), 
  #                     end(xdata[xsubset]), sep = " / "), col = theme$labels, 
  #                     adj = c(0, 0), cex=0.5, pos = 2)))
  text.exp.1 <- c(  expression(text(xlim[2], 0.5, paste(start(xdata[xsubset]), 
                                                        end(xdata[xsubset]), sep = " / "), col = theme$labels, 
                                    adj = c(0, 0), cex=0.8, pos = 2)))  #Replacement line LFM
  cs$add(text.exp, env = cs$Env, expr = TRUE)
  cs$add(text.exp.1, env = cs$Env, expr = TRUE)
  cs$set_frame(2)
  exp <- expression(segments(xlim[1], y_grid_lines(get_ylim()[[2]]), 
                             xlim[2], y_grid_lines(get_ylim()[[2]]), col = theme$grid, 
                             lwd = grid.ticks.lwd, lty = grid.ticks.lty))
  if (yaxis.left) {
    exp <- c(exp, expression(text(xlim[1], y_grid_lines(get_ylim()[[2]]), 
                                  noquote(format(y_grid_lines(get_ylim()[[2]]), justify = "right")), 
                                  col = theme$labels, srt = theme$srt, offset = 1, 
                                  pos = 2, cex = theme$cex.axis, xpd = TRUE)))
  }
  if (yaxis.right) {
    exp <- c(exp, expression(text(xlim[2], y_grid_lines(get_ylim()[[2]]), 
                                  noquote(format(y_grid_lines(get_ylim()[[2]]), justify = "right")), 
                                  col = theme$labels, srt = theme$srt, offset = 1, 
                                  pos = 4, cex = theme$cex.axis, xpd = TRUE)))
  }
  cs$add(exp, env = cs$Env, expr = TRUE)
  cs$set_frame(2)
  if (isTRUE(multi.panel)) {
    lenv <- cs$new_environment()
    lenv$xdata <- cs$Env$xdata[subset, 1]
    lenv$label <- colnames(cs$Env$xdata[, 1])
    lenv$type <- cs$Env$type
    if (yaxis.same) {
      lenv$ylim <- cs$Env$constant_ylim
    }
    else {
      lenv$ylim <- range(cs$Env$xdata[subset, 1], na.rm = TRUE)
    }
    exp <- quote(chart.lines(xdata, type = type, lty = lty, 
                             lwd = lwd, lend = lend, col = theme$col, up.col = theme$up.col, 
                             dn.col = theme$dn.col, legend.loc = legend.loc))
    exp <- as.expression(add.par.from.dots(exp, ...))
    cs$add(exp, env = lenv, expr = TRUE)
    text.exp <- expression(text(x = xycoords$x[2], y = ylim[2] * 
                                  0.9, labels = label, col = theme$labels, adj = c(0, 
                                                                                   0), cex = 1, offset = 0, pos = 4))
    cs$add(text.exp, env = lenv, expr = TRUE)
    if (NCOL(cs$Env$xdata) > 1) {
      for (i in 2:NCOL(cs$Env$xdata)) {
        lenv <- cs$new_environment()
        lenv$xdata <- cs$Env$xdata[subset, i]
        lenv$label <- cs$Env$column_names[i]
        if (yaxis.same) {
          lenv$ylim <- cs$Env$constant_ylim
        }
        else {
          yrange <- range(cs$Env$xdata[subset, i], na.rm = TRUE)
          if (all(yrange == 0)) 
            yrange <- yrange + c(-1, 1)
          lenv$ylim <- yrange
        }
        lenv$type <- cs$Env$type
        lenv$lty <- cs$Env$lty[i]
        lenv$lwd <- cs$Env$lwd[i]
        lenv$col <- cs$Env$theme$col[i]
        cs$add_frame(ylim = c(0, 1), asp = 0.25)
        cs$next_frame()
        text.exp <- expression(text(x = xlim[1], y = 0.5, 
                                    labels = "", adj = c(0, 0), cex = 0.9, offset = 0, 
                                    pos = 4))
        cs$add(text.exp, env = lenv, expr = TRUE)
        cs$add_frame(ylim = lenv$ylim, asp = NCOL(cs$Env$xdata), 
                     fixed = TRUE)
        cs$next_frame()
        exp <- quote(chart.lines(xdata[xsubset], type = type, 
                                 lty = lty, lwd = lwd, lend = lend, col = col, 
                                 up.col = theme$up.col, dn.col = theme$dn.col, 
                                 legend.loc = legend.loc))
        exp <- as.expression(add.par.from.dots(exp, ...))
        exp <- c(exp, expression(segments(xlim[1], y_grid_lines(ylim), 
                                          xlim[2], y_grid_lines(ylim), col = theme$grid, 
                                          lwd = grid.ticks.lwd, lty = grid.ticks.lty)), 
                 expression(x_grid_lines(xdata[xsubset], grid.ticks.on, 
                                         ylim)))
        if (yaxis.left) {
          exp <- c(exp, expression(text(xlim[1], y_grid_lines(ylim), 
                                        noquote(format(y_grid_lines(ylim), justify = "right")), 
                                        col = theme$labels, srt = theme$srt, offset = 1, 
                                        pos = 2, cex = theme$cex.axis, xpd = TRUE)))
        }
        if (yaxis.right) {
          exp <- c(exp, expression(text(xlim[2], y_grid_lines(ylim), 
                                        noquote(format(y_grid_lines(ylim), justify = "right")), 
                                        col = theme$labels, srt = theme$srt, offset = 1, 
                                        pos = 4, cex = theme$cex.axis, xpd = TRUE)))
        }
        cs$add(exp, env = lenv, expr = TRUE, no.update = TRUE)
        text.exp <- expression(text(x = xycoords$x[2], 
                                    y = ylim[2] * 0.9, labels = label, col = theme$labels, 
                                    adj = c(0, 0), cex = 1, offset = 0, pos = 4))
        cs$add(text.exp, env = lenv, expr = TRUE)
      }
    }
  }
  else {
    if (type == "h" & NCOL(x) > 1) 
      warning("only the univariate series will be plotted")
    exp <- quote(chart.lines(xdata[xsubset], type = type, 
                             lty = lty, lwd = lwd, lend = lend, col = theme$col, 
                             up.col = theme$up.col, dn.col = theme$dn.col, legend.loc = legend.loc))
    exp <- as.expression(add.par.from.dots(exp, ...))
    cs$add(exp, expr = TRUE)
    assign(".xts_chob", cs, .plotxtsEnv)
  }
  if (!is.null(panels) && nchar(panels) > 0) {
    panels <- parse(text = panels, srcfile = NULL)
    for (p in 1:length(panels)) {
      if (length(panels[p][[1]][-1]) > 0) {
        cs <- eval(panels[p])
      }
      else {
        cs <- eval(panels[p])
      }
    }
  }
  assign(".xts_chob", cs, .plotxtsEnv)
  cs
}
environment(plotxts) <- asNamespace("xts")
