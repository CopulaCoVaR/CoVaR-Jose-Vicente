#----
#-------------------------------------------------------------------------------------#
# Descripción: La función permite hacer un rolling forecast basada en  un modelo
#              ARMA-GARCH, asumiendo una distribución t asimétrica por defecto. 
#-------------------------------------------------------------------------------------#
# Argumentos:
# <Data>               : Vector de datos sobre el cuál se quiere realizar el rolling forecast
# <model>              : Tipo de modelo garch a tener en cuenta, por defecto 'gjrGARCH', se 
#                        pueden usar todos los modelos asociados al paquete <rugarch>.
# <garchOrder>         : vector que determina el orden del GARCH, por defecto c(1,1)
# <armaOrder>          : vector que determina el orden ARMA para el modelo de la media 
#                        en el GARCH
# <distribution.model> : Por defecto 'sstd' (skewed t), se puede usar cualquier 
#                        distribución disponible dentro del paquete rugarch. 
# <window.size>        : valor numérico entero que determina el tamaño de la ventana de 
#                        estimación. 
# <refit.every>        : Valor numérico entero que determina cada cuantos pronósticos se
#                        debe reestimar el modelo. 
# <n.ahead>            : Valor numérico entero que determina el horizonte de pronóstico. 
# 
# <show.progress>      : Variable lógica que especifica sí se debería guardar la gráfica asociada
#                        a cada una de las medidas de riesgo calculadas con el rolling. 
# <external.regressors>: Por defecto NULL. Si no es NULL se brindan los nombres de los regresores 
#                        externos que se incluyen en la ecuación de la media del modelo.                 
        
#-------------------------------------------------------------------------------------#
# Salida: La función retorna una matriz que contiene los valores pronosticados de
#         mean, sigma necesarios para el calculo del VaR, así como sus valores 
#         estimados de Skewness y Shape asociados a la t asimétrica. 
#-------------------------------------------------------------------------------------#
#----
garch.roll = function(Data=Datos[,serie], model="gjrGARCH", garchOrder=c(1,1), armaOrder=NULL, 
                      distribution.model="sstd", window.size=150, refit.every=20, 
                      n.ahead=1, show.progress=TRUE, external.regressors=NULL){
  if (is.null(armaOrder)==TRUE) {
    armaOrder = ARMA.ORDER.DF(Datos=Datos[,serie]) #Usando valores por defecto
  }
  #Vector de ceros con el índice temporal teniendo en cuenta la pérdida de datos del rolling forecast.
  forecast_index     = Data[(window.size+n.ahead):nrow(Data),] - Data[(window.size+n.ahead):nrow(Data),]  
  #Definicion de la matriz de salida
  forecast           = cbind(forecast_index, forecast_index, forecast_index, 
                             forecast_index, forecast_index, forecast_index)
  colnames(forecast) = c('original','mean', 'residuals', 'sigma','skew', 'shape')
  total.rolling  = nrow(Data)- window.size  
  
  if (Sys.info()["sysname"] == 'Darwin') #Compatible con cualqueir sistema operativo (Escrito originalmente para MacOs)
    pb <- txtProgressBar(  min   = 0,
                           max   = total.rolling,
                           style = 3,
                           width = total.rolling, # Necesario para evitar múltiples prints
                           char  = "=") 
  if (Sys.info()["sysname"] == 'Windows') # Del paquete 'utils', instalado por defecto en Windows. 
    pb <- winProgressBar(title = "Barra de progreso de Windows", # Título de la ventana
                         label = "Porcentaje completado",        # Texto de la ventana
                         min = 0,                                # Valor mínimo de la barra
                         max = total.rolling,                    # Valor máximo de la barra
                         initial = 0,                            # Valor inicial de la barra
                         width = 300L)                           # Ancho de la ventana 
  inicio = total.rolling 
  fin    = total.rolling
  
  error.estim    = 0
  first.estim    = 0
  
  # Especificación del modelo.
  model.spec         = ugarchspec(variance.model = list(model=model, garchOrder = garchOrder), 
                                  mean.model = list(armaOrder=armaOrder,
                                                    external.regressors=external.regressors), 
                                  distribution.model=distribution.model)
  for (i in 1:total.rolling){
    #cat(' i:',i)
    if(show.progress==TRUE | Sys.info()["sysname"] == 'Darwin') inicio[i]=Sys.time()
    window=Data[i:(window.size+i-1),]
    if (i==1 | i%%refit.every==0){ #Permite reestimar el modelo cada <refit.every> ventanas
      # if(show.progress == TRUE) cat('[',i,']','','refitting model...','\n')
      model.fit           = try(ugarchfit(spec=model.spec, data=window, solver='hybrid'), TRUE)
      if(class(model.fit) == 'try-error' | is.null(coef(model.fit))) error.estim=1 
      else {
             #print(coef(model.fit)) #Quitar!!!!
        if (armaOrder[1]>0){ 
          Mod.roots=Mod(roots(c(1,(-coef(model.fit)[2:(armaOrder[1]+1)]))))
          #Modelo no estacionario
          if (sum(Mod.roots<=1 & Mod.roots>=(-1))==TRUE) error.estim=1 
          else error.estim=0
        }
        else error.estim=0
      }
      if(error.estim      == 0){  
         #cat('\nT1') #Quitar!!!
        window.spec              = getspec(model.fit)
        setfixed(window.spec)    = as.list(coef(model.fit))
        model.fit.0              = model.fit
        window.spec.0            = window.spec
        first.estim              = 1
      }
      # Cuando no se puede estimar el modelo
      if (error.estim == 1){
        if (first.estim == 0){
          warning('No se ha podido estimar ningún modelo, error en ventana ',i,'\n')
          #cat('\nT13') #Quitar
          model.spec         = ugarchspec(variance.model = list(model=model, garchOrder = garchOrder), 
                                          mean.model = list(armaOrder=c(1,0),
                                          external.regressors = external.regressors), 
                                          distribution.model = distribution.model)
          model.fit           = try(ugarchfit(spec=model.spec, data=window, solver='hybrid'), TRUE)
          window.spec         = getspec(model.fit)
          setfixed(window.spec)    = as.list(coef(model.fit))
        }else{
          model.fit                = model.fit.0
          window.spec              = window.spec.0
          setfixed(window.spec)    = as.list(coef(model.fit))
        }
      }
    }
    #cat('\nT2, error.estim:', error.estim,"first.estim:",first.estim ); print(window.spec) #quitar!!
    model.forecast            = ugarchforecast(window.spec, n.ahead=n.ahead, data=window)
    forecast[i, 'original']   = Data[window.size+i-1+n.ahead,] 
    forecast[i, 'mean']       = model.forecast@forecast$seriesFor
    # if(forecast[i, 'mean']>(15) ){ #| forecast[i, 'mean']> 20#
    #   cat('\n ========================================= mean:', forecast[i, 'mean'])
    #   print(coef(model.fit))
    #   print(model.fit)
    # }
    forecast[i, 'residuals']  = Data[window.size+i-1+n.ahead,] - model.forecast@forecast$seriesFor
    forecast[i, 'sigma']      = sigma(model.forecast)
    forecast[i, 'skew']       = coef(model.fit)['skew']
    forecast[i, 'shape']      = coef(model.fit)['shape']
    
    if (show.progress == TRUE) {
      if (Sys.info()["sysname"] == 'Darwin'){
        fin[i] <- Sys.time()
        setTxtProgressBar(pb, i)
        tiempo <- round(seconds_to_period(sum(fin - inicio)), 0)
        
        # Tiempo restante estimado basado en el tiempo medio que tardaron en ejecutarse las iteraciones previas
        est <- total.rolling * (mean(fin[fin != 0] - inicio[inicio != 0])) - tiempo
        restante <- round(seconds_to_period(est), 0)
        
        cat(paste0(" // Tiempo ejecución:", tiempo," // Tiempo restante estimado:", restante))
      }
      if (Sys.info()["sysname"] == 'Windows'){
        pctg <- paste(round(i/total.rolling *100, 0),"% completado")
        setWinProgressBar(pb, i, label = pctg) # Al pasar pctg al argumento label se sobreescribirá
                                               # el texto del argumento label de la función winProgressBar
      }
    }
  }
  if (sum(is.na(forecast))!= 0) {
    forecast = na.fill(forecast, 'extend')
  }
  return(forecast)
}

# ------------------------------------------------------------------------

#------------------------------------- Ejemplo ----------------------------
if (0) {
  Prueba_ARG           = garch.roll(Data=Datos[,'ARG'], refit.every = 20)
  Prueba_ARG_Ext       = garch.roll(Data=Datos[,'ARG'], refit.every = 20, external.regressors = Datos[,2])
  VaR                  = as.xts(Prueba_IND$mean + Prueba_IND$sigma*(rugarch::qdist(distribution='sstd', p=ALPHA,   shape=Prueba_IND$shape,skew=Prueba_IND$skew)))
  VaRUp                = as.xts(Prueba_IND$mean + Prueba_IND$sigma*(rugarch::qdist(distribution='sstd', p=1-ALPHA, shape=Prueba_IND$shape,skew=Prueba_IND$skew)))
  x11()
  plot(VaR)
  abline(VarUp)
  save.image(file = paste0(today(),'.RData'))
  load(paste0('2022-08-04','.RData'))
  # 
  # model.spec         = ugarchspec(variance.model = list(model=model, garchOrder = garchOrder), 
  #                                 mean.model = list(armaOrder=armaOrder), 
  #                                 distribution.model = distribution.model)
  # model_fit=ugarchfit(spec=model.spec, data=Datos[,1], solver='hybrid'); model_fit
  # 
  # model.spec_ext         = ugarchspec(variance.model = list(model=model, garchOrder = garchOrder), 
  #                                 mean.model = list(armaOrder=armaOrder, external.regressors = Datos[,2:3]), 
  #                                 distribution.model = distribution.model)
  # model_fit_ext=ugarchfit(spec=model.spec, data=Datos[,1], solver='hybrid')
  # 
}


