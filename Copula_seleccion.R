#-------------------------------------------------------------------------------#
# Objetivo: Permite determinar el orden p y q para una serie de tiempo dada. 
#-------------------------------------------------------------------------------#
# Argumentos de la función: 
#-------------------------------------------------------------------------------#
# <ts_object>: Serie de tiempo a la cuál se le quiere determinar el orden p,q
# <AR.m>: Orden autoregresivo máximo a tener en cuenta.
# <MA.m>: Orden de media móvil máximo a tener en cuenta.
# <d>: Orden de integración
# <drift>: TRUE incluye término de deriva en el modelo, FALSE lo excluye. 
# <metodo>: Método de estimación, por defecto máxima verosimilitud. 
#-------------------------------------------------------------------------------#
# Salida: 
# Matriz con los ordenes <p>, <q> y criterios optimos bajo AIC (fila 1) y BIC( fila 2)
#-------------------------------------------------------------------------------#
arma.seleccion = function(ts_object, AR.m=5, MA.m=0, d=0, drift=TRUE, metodo="ML"){
  index = 1
  data = data.frame(p=integer(), d=integer(), q=integer(), AIC=double(), BIC=double())
  for (p in 0:AR.m) {
    for (q in 0:MA.m)  {
      fitp <- arima(ts_object, order=c(p, d, q), include.mean=drift,
                    method=metodo, optim.control=list(maxit=1000))
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
#-------------------------------------------------------------------------------#
# Copula seleccion: Esta funcion se encarga de escoger la mejor copula que modele cada pareja
# de datos entre una de las acciones del G20 y el petroleo, se estiman copulas t, Normal, Clayton,
# Frank, Gumbel, AMH, Joe, BB7, BB7 rotada, Gumbel Rotada, Clayton rotada, Normal dinamica,
# t dinamica, Gumbel dinamica y Gumbel rotada dinamica, luego se obtiene el respectivo AIC para
# cada una de las copulas obtenidas y se determina cual es aquella que modela de mejor manera
# la estructura de datos e.g accion/petroleo.
#-------------------------------------------------------------------------------#
# Argumentos de la funcion:
#-------------------------------------------------------------------------------#
# <muestra>: Datos, se espera dividir la muestra entre periodos con diferentes caracteristicas
# de volatilidad que pueden dar origen a diferentes relaciones entre las series.
# <arma.order>: 1. Vector de 4 elementos dónde los primeros dos son los ordenes <p> y <q> de la serie
#                  1, y los elementos 3 y 4 son los ordenes <p> y <q> de la serie 2. 
#               2. Si es <NULL>, la función encuentra los ordenes <p> y <q> óptimos para las series 1 y 2 
#                  basado en la función <arma.seleccion>.
# <GARCH.model>: Tipo de modelo garch a usar para modelar la varianza, por defecto gjrGARCH.
# <serie.1>: Primera serie
# <serie.2>: Segunda serie (que se choca), se supone que siempre es la misma, e.g WTI. 
#-------------------------------------------------------------------------------#
# Salida:
# Lista que contiene los criterios de información asociados a cada cópula y aquella que mejor 
# se ajusta a la pareja de datos. 
#-------------------------------------------------------------------------------#
Copula_seleccion = function(Datos=G20_Crisis_WNA, arma.order=NULL, GARCH.model="gjrGARCH", serie.1='Argentina', serie.2='WTI', CoVaR.type='Equal'){
  #Seleccion de los ordenes de los modelos ARMA. Nota: <AR.max>, <MA.max> y <drift> estan fijos 
  #Seleccion para la <serie.1>
  ARMA.ORDER.1 <- matrix(0,1,2, dimnames = list(serie.1, c('p','q')))
  if (is.null(arma.order)==FALSE) ARMA.ORDER.1[1,] = arma.order[1:2]
  else{
    p.q = arma.seleccion(Datos[,serie.1], AR.m=5, MA.m=0, d=0, drift=TRUE)  #<--- <AR.max>, <MA.max>, <drift>
    ARMA.ORDER.1[1,] = p.q['BIC', c('p', 'q')]
  } 
  #Seleccion para la <serie.2>
  ARMA.ORDER.2 <- matrix(0,1,2, dimnames = list(serie.2, c('p','q')))
  if (is.null(arma.order)==FALSE) ARMA.ORDER.2[2,] = arma.order[3:4]
  else{
    p.q = arma.seleccion(Datos[,serie.2], AR.m=5, MA.m=0, d=0, drift=TRUE)  #<--- <AR.max>, <MA.max>, <drift>
    ARMA.ORDER.2[1,] = p.q['BIC', c('p', 'q')]
  } 
  
  #Estimacion los GARCH para <serie.1> y <serie.2>, se supone distribucion <sstd>
  model.1     = ugarchspec(variance.model=list(model=GARCH.model, garchOrder=c(1,1)), 
                           mean.model=list(armaOrder=ARMA.ORDER.1), distribution.model="sstd")
  model.1.fit = ugarchfit(spec=model.1, data=Datos[,serie.1], solver='hybrid' )
  model.2     = ugarchspec(variance.model=list(model=GARCH.model, garchOrder=c(1,1)), 
                           mean.model=list(armaOrder=ARMA.ORDER.2), distribution.model="sstd")
  model.2.fit = ugarchfit(spec=model.2, data=Datos[,serie.2], solver='hybrid')
  
  #pobs (pseudo muestras suponiendo distribuciones marginales <sstd>)
  sd_errors1 = residuals(model.1.fit)/sigma(model.1.fit)
  u1         = rugarch::pdist(distribution='sstd', q=sd_errors1, skew=coef(model.1.fit)['skew'], shape=coef(model.1.fit)['shape'])
  sd_errors2 = residuals(model.2.fit)/sigma(model.2.fit)
  u2         = rugarch::pdist(distribution='sstd', q=sd_errors2, skew=coef(model.2.fit)['skew'], shape=coef(model.2.fit)['shape'])
  u.s        = cbind(u1, u2)
  
  #Estimacion todas las copulas consideradas, incluyendo AIC y BIC
  #Copulas consideradas para la metodologia de Liu (2022)  
  if (CoVaR.type == 'Equal') {
    CoPuLaS    = c('Student',         'Gaussian',  'Clayton',       'Frank',     'Gumbel', 'Joe',             'BB7',       'Rot.BB7',       'Rot.Gumbel','Rot.Clayton',
                   'Dyn.Gaussian',    'Dyn.Student',     'Dyn.Gumbel',    'Dyn.Rot.Gumbel')
    CoPuLaS.Fam= c( 2,                 1,            3,               5 ,          4,       6,                 9,           19,               14,          13,
                   -1,                -2,           -3,              -4)
  }
  #Copulas consideradas para la metodologia de RU (2016) 
  if (CoVaR.type == 'Less') {
    CoPuLaS    = c('Student',         'Gaussian',  'Clayton',       'Frank',       'Gumbel',    'AMH',    
                   'Joe',             'BB7',       'Rot.BB7',       'Rot.Gumbel',  'Rot.Clayton',
                   'Dyn.Gaussian',    'Dyn.Student',     'Dyn.Gumbel',    'Dyn.Rot.Gumbel', 'Plackett')
    CoPuLaS.Fam= c( 2,                 1,            3,                5 ,          4,          -9,
                    6,                 9,           19,               14,          13,
                   -1,                -2,           -3,               -4,          10)
  }
  #Estimacion de las copulas
  N.Cop = length(CoPuLaS)
  IC = matrix(NA,N.Cop,2, dimnames=list(CoPuLaS, c('AIC','BIC')))
  cat('\nCop:')
  for ( i in 1:N.Cop ) {
    if(CoPuLaS.Fam[i]>=0) {Cop.Model = try(BiCopEst(u1, u2, family = CoPuLaS.Fam[i]))
       if(class(Cop.Model)=='try-error') Cop.Model$AIC=Cop.Model$BIC=99999
    }
    if(CoPuLaS.Fam[i]==-1){ Cop.Model = try(dynamicnormal.LF(u.s, plot=FALSE))
      if(class(Cop.Model)=='try-error') Cop.Model$AIC=Cop.Model$BIC=99999
    }
    if(CoPuLaS.Fam[i]==-2){ Cop.Model = try(dynamicT.LF(data=u.s, plot=FALSE, printout=FALSE))
      if(class(Cop.Model)=='try-error') Cop.Model$AIC=Cop.Model$BIC=99999
    }
    if(CoPuLaS.Fam[i]==-3){ Cop.Model = try(dynamicGum.LF(u.s))
      if(class(Cop.Model)=='try-error') Cop.Model$AIC=Cop.Model$BIC=99999
    }
    if(CoPuLaS.Fam[i]==-4){ Cop.Model = try(dynamicGum.LF(u.s, rotated=TRUE))   
      if(class(Cop.Model)=='try-error') Cop.Model$AIC=Cop.Model$BIC=99999
    }
    if(CoPuLaS.Fam[i]==-9){ 
        cop_model= amhCopula(dim=2)
        Cop.Model= try(fitCopula(cop_model, u.s , method='ml'))
        if(class(Cop.Model)=='try-error'){IC[i,'AIC']=IC[i,'BIC']=99999}
        else{      
          IC[i,'AIC'] = AIC(Cop.Model)
          IC[i,'BIC'] = BIC(Cop.Model)
        }
    }
    if(CoPuLaS.Fam[i]==10){ 
        cop_model= plackettCopula()
        Cop.Model= try(fitCopula(cop_model, u.s , method='ml'))
        if(class(Cop.Model)=='try-error'){IC[i,'AIC']=IC[i,'BIC']=99999}
        else{      
          IC[i,'AIC'] = AIC(Cop.Model)
          IC[i,'BIC'] = BIC(Cop.Model)
        }
    }
    if(CoPuLaS.Fam[i] != -9 & CoPuLaS.Fam[i] != 10){ 
      IC[i,'AIC'] = Cop.Model$AIC
      IC[i,'BIC'] = Cop.Model$BIC
    }
    #cat('\nCop:',CoPuLaS[i],', AIC:',  IC[i,'AIC'], ', BIC:', IC[i,'BIC'])
    cat(' ',CoPuLaS[i])
  }
  #Nota: no se incluyeron las copulas <Dynamic_BB7> y <Dynamic_Rotated_BB7> ya que no convergian
  if(0){
   #Dynamic BB7 <Luis Melo>
   #cop_model_Dynamic_BB7 <- dynamicBB7.LF(u.s)         
   #AIC_model_Dynamic_BB7 <- cop_model_Dynamic_BB7$AIC                  #Error in solve.default(model$hessian) : 
   #Dynamic rotated BB7 <Luis Melo>
   #cop_model_Dynamic_Rotated_BB7 <- dynamicBB7.LF(u.s)                 
   #AIC_model_Dynamic_Rotated_BB7 <- cop_model_Dynamic_Rotated_BB7$AIC  #Error in solve.default(model$hessian) : 
  }

  #----Presentación de resultados----#
  bestfit        = c(names(which.min(IC[,'AIC'])), names(which.min(IC[, 'BIC'])))
  names(bestfit) = c('min.AIC', 'min.BIC')
  list( IC.Copulas=IC, Best.Copula=bestfit) 
}