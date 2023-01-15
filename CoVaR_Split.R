# Preparación -------------------------------------------------------------
setwd("C:/Users/maico/OneDrive - Universidad Nacional de Colombia/BanRep/Value at Risk/Highdimensional CoVaR network connectedness/Github-R/Copula-CoVaR-")
#setwd('C:\\Respaldo DLO\\Copula_CoVaR\\R')   #<<<<<--- Carpeta general
#setwd('/Users/lumelo/archivos/Copula_CoVaR/Github-R/Copula-CoVaR-')
wd = getwd()                               # Carpeta de trabajo adaptable al proyecto, se modifica automáticamente en cada computador.
Resultados <<- paste0(wd,'/Resultados')    #<<<<<--- Capeta de resultados qué depende deldirectorio de trabajo. 
CoVaR.type  = c('Equal','Less')[2]         #<<<--- Seleccion del tipo de metodologia para calcular el CoVaR <Equal> para Liu(2022) y <Less> para RU(2016)
data.file   = c('CDS_Data_VIX.xlsx','CDS_Data_usfcon.xlsx')[1]          #<<<--- String con el nombre del archivo de datos en la carpeta de trabajo.
log.Serie.1 = c(TRUE, FALSE)[2]            #<<<--- <T> se aplican logaritmos a las <Serie.1>, <F> No se aplican.
log.Serie.2 = c(TRUE, FALSE)[2]            #<<<--- <T> se aplican logaritmos a la <Serie.2>, <F> No se aplica.
dif.Serie.1 = c(TRUE, FALSE)[1]            #<<<--- <T> diferencia las <Serie.1>, <F> deja los datos en nivel. 
dif.Serie.2 = c(TRUE, FALSE)[1]            #<<<--- <T> diferencia la <Serie.2>, <F> deja los datos en nivel. 
crisis.break= c(2008)[1]                   
CoVaR.plot  = c('Up', 'Down', 'Both')[1]   #<<<--- 'Up' grafica <VaR.Up>, <CoVaR.Up> y <DeltaCoVaR.Up>. 'Down' grafica <VaR>, <CoVaR> y <DeltaCoVaR>. 'Both' grafica todo. 
VaR.plot    = c(TRUE, FALSE)[1]            #<<<--- <T>  grafica <VaR.Up> y <VaR>, <F> no los grafica. 
BigVAR      = c(TRUE, FALSE)[2]            #<<<--- <F> si solo se estima el CoVaR, <T> si ademas se realiza con VAR y conectividad
n.ahead.connectedness=20                   #<<<--- Número entero que determina el horizonte de pronóstico para las GIRF, GFEVD y conectividad.     
# Paquetes ----------------------------------------------------------------
library(ichimoku) # Para transformar de xts a df
library(xts)
library(circlize)
library(car)
library(readxl)
library(fpp3)
library(vars)
library(forecast)
library(copula)
library(ggplot2)
library(tseries)
library(rugarch)
library(moments)
library(aTSA)
library(FinTS)
library(VineCopula)
library(VC2copula)
library(sgt)
library(pracma)
library(quantmod)
library(PerformanceAnalytics)
library(nloptr)
library(expm)
library(vars)
library(lubridate)
library(BigVAR)
library(readxl)
library(tidyquant)
library(reshape2)
library(igraph)
# Funciones ---------------------------------------------------------------
source('RebUgo_Functions.R')
source('Copula_seleccion.R')
source('Copula_CoVaR.R')
source('DynCopulaCoVaR.R')
source('DynCopulaCoVaRUpper.R')
source('GIRF_BigVAR.R')
source('GFEVD_BigVAR.R')
source('Connectedness.R')
source('CoVaR_LF.R')
source('Rolling.R')
source('garch_roll.R')
source('Optimal_model.R')
# Datos -------------------------------------------------------------------
#--Base de datos de RETORNOS para el periodo de crisis. 
#---- Supuesto_1: la primera col es la FECHA de la base de datos
#---- Supuesto_2: Las col.2 hasta la penultima col. son las series que se modelan w.r.t la serie de la ultima columna (Ej. Ind Acc de varios paises)
#---- Supuesto_3: La ultima col. es la serie que se mantiene en todas copulas (Ej: WTI)

DATOS = read_xlsx(data.file)  #<<<---- Base de datos

N.Series1=ncol(DATOS)- 1 -1                  #Se resta la primera col. (FECHA) y la ultima columna (Ej: WTI)
Name.Serie2=colnames(DATOS)[ncol(DATOS)]  #Nombre de la ultima col. es la serie que se mantiene en todas copulas (Ej: WTI)
Name.Series1=colnames(DATOS[,-1])[1:N.Series1]  #Se quitan la col de  fechas  y se definen las series.1 (i.e. country' Stock returns)

# As XTS
DATOS_xts = xts(DATOS[,-1], order.by = as.POSIXct(DATOS[[1]])) #Se resta la primera col. (FECHA) 
# Interpolacion lineal para completar los datos faltantes (para datos interiores, en los extremos toma el valor mas cercano)
NAs=sum(is.na(DATOS_xts)) # Numero datos faltantes
if (NAs!=0){
  cat('\nNo. de datos faltantes: ',NAs,'\n') # Numero datos faltantes
  DATOS_WNA = na.fill(DATOS_xts, 'extend')
}else DATOS_WNA=DATOS_xts
#Se transforman las series si se requiere
#Logaritmos
if(log.Serie.1==TRUE) DATOS_WNA[,-ncol(DATOS_WNA)]=log(DATOS_WNA[,-ncol(DATOS_WNA)])  # Logaritmo de las Serie.1
if(log.Serie.2==TRUE) DATOS_WNA[, ncol(DATOS_WNA)]=log(DATOS_WNA[, ncol(DATOS_WNA)])  # Logaritmo de la Serie.2                                  # Logaritmos a todas las series.
#Diferencias
if(dif.Serie.1==TRUE) DATOS_WNA[,-ncol(DATOS_WNA)]=diff(DATOS_WNA[,-ncol(DATOS_WNA)]) #Diferencia de las  Serie.1
if(dif.Serie.2==TRUE) DATOS_WNA[, ncol(DATOS_WNA)]=diff(DATOS_WNA[, ncol(DATOS_WNA)]) #Diferencia de la  Serie.2
if(dif.Serie.1==TRUE|dif.Serie.2==TRUE) DATOS_WNA =DATOS_WNA[-1,]                     #Eliminamos la primera fila de datos (perdidos en la diferencia)

# Subperiodos -------------------------------------------------------------
if (crisis.break==2008){
  DATOS_PRE = window(DATOS_WNA, end='2008-09-15'); NAs=sum(is.na(DATOS_PRE))     # Numero datos faltantes
  if (NAs!=0) DATOS_PRE  = na.fill(DATOS_xts, 'extend')
  DATOS_POST= window(DATOS_WNA, start='2008-09-16'); NAs=sum(is.na(DATOS_POST)) # Numero datos faltantes
  if (NAs!=0) DATOS_POST = na.fill(DATOS_xts, 'extend')
}
# Estadísticas descriptivas -----------------------------------------------
if(1){
  Series = colnames(DATOS_WNA)
  N      = length(Series)
  stats  = matrix(0,10,N,dimnames=list(c('Mean','Max','Min','SD','Skew','Kurt','J-B','ARCH' ,paste0('Corr.',Name.Serie2),'No.NA'),Series))
  ARMA.Order = ARMA.ORDER.DF(Datos=DATOS_WNA)
  for (i in 1:N){
    arima      = arima(DATOS_WNA[,Series[i]], method = 'ML', order=c(ARMA.Order[Series[i],],0))
    ARCH       = arch.test(arima)
    stats[1,i] = mean(DATOS_WNA[,Series[i]])
    stats[2,i] = max(DATOS_WNA[,Series[i]])
    stats[3,i] = min(DATOS_WNA[,Series[i]])
    stats[4,i] = sd(DATOS_WNA[,Series[i]])
    stats[5,i] = skewness(DATOS_WNA[,Series[i]])
    stats[6,i] = kurtosis(DATOS_WNA[,Series[i]])
    stats[7,i] = jarque.bera.test(DATOS_WNA[,Series[i]])$p.value
    stats[8,i] = ARCH[5,'LM']
    stats[9,i] = cor(DATOS_WNA[,Series[i]],DATOS_WNA[,Name.Serie2], method = 'kendall')
    stats[10,i] = sum(is.na(DATOS_WNA[,Series[i]]))
  }
  print(t(stats), digits=3)
  write.csv(stats, file='Descriptive.csv')
}
# Graficación -------------------------------------------------------------
if(0){
  #Data.graph = DATOS_xts
  Data.graph = DATOS_WNA
  # Gráfica conjunta
  x11()
  plot(Data.graph, format.labels="%Y", major.ticks = 'years', grid.col = NA)
  print(addLegend("topright", legend.names = colnames(Data.graph), 
                  lty=1, lwd=1))#Gráficas individual
  x11()
  par(mfrow=c(4,3))
  for (i in 1:(ncol(Data.graph))) {
    print(plot(Data.graph[,i], main=paste(colnames(Data.graph)[i]), 
               col='steelblue', format.labels="%Y",
               major.ticks = 'years', grid.col = NA))
  }
}
# Selección de Copula -----------------------------------------------------
if(1){
  time.ini  = Sys.time()
  Serie2    = Name.Serie2 #Ej: WTI
  Res.Tot   = array(NA,dim=c(N.Series1,2), 
                    dimnames=list(Name.Series1, c('Min.AIC', 'Min.BIC')))
  ii.n      = 0 
  for (ii in Name.Series1){
    print(ii)
    ii.n            = ii.n + 1
    res             = Copula_seleccion(Datos=DATOS_WNA, arma.order=NULL, GARCH.model="gjrGARCH", serie.1=ii, serie.2=Serie2, CoVaR.type=CoVaR.type)
    Res.Tot[ii.n, ] = res$Best.Copula
    print(res)
  }
  
  #-- Time of execution --#
  time.diff = (Sys.time() - time.ini) 
  time.diff.over.60 = (Sys.time() - time.ini)/60 
  
  #-- Printout of final results --#
  print(Res.Tot)
  print(time.diff)
  print(time.diff.over.60)
  Copulas.seleccionadas = Res.Tot
  save(Copulas.seleccionadas, file=paste0('Cop_seleccionadas_',today()))
} # Toda la muestra
if (crisis.break==2008) {
  if(1){
    time.ini  = Sys.time()
    Serie2    = Name.Serie2 #Ej: WTI
    Res.Tot   = array(NA,dim=c(N.Series1,2), 
                      dimnames=list(Name.Series1, c('Min.AIC', 'Min.BIC')))
    ii.n      = 0 
    for (ii in Name.Series1){
      print(ii)
      ii.n            = ii.n + 1
      res             = Copula_seleccion(Datos=DATOS_PRE, arma.order=NULL, GARCH.model="gjrGARCH", serie.1=ii, serie.2=Serie2, CoVaR.type=CoVaR.type)
      Res.Tot[ii.n, ] = res$Best.Copula
      print(res)
    }
    
    #-- Time of execution --#
    time.diff = (Sys.time() - time.ini) 
    time.diff.over.60 = (Sys.time() - time.ini)/60 
    
    #-- Printout of final results --#
    print(Res.Tot)
    print(time.diff)
    print(time.diff.over.60)
    Copulas.seleccionadas = Res.Tot
    save(Copulas.seleccionadas, file=paste0('Cop_seleccionadas_Pre_Crisis_Vix',today()))
  } # Pre_Crisis (2008)
  if(1){
    time.ini  = Sys.time()
    Serie2    = Name.Serie2 #Ej: WTI
    Res.Tot   = array(NA,dim=c(N.Series1,2), 
                      dimnames=list(Name.Series1, c('Min.AIC', 'Min.BIC')))
    ii.n      = 0 
    for (ii in Name.Series1){
      print(ii)
      ii.n            = ii.n + 1
      res             = Copula_seleccion(Datos=DATOS_POST, arma.order=NULL, GARCH.model="gjrGARCH", serie.1=ii, serie.2=Serie2, CoVaR.type=CoVaR.type)
      Res.Tot[ii.n, ] = res$Best.Copula
      print(res)
    }
    
    #-- Time of execution --#
    time.diff = (Sys.time() - time.ini) 
    time.diff.over.60 = (Sys.time() - time.ini)/60 
    
    #-- Printout of final results --#
    print(Res.Tot)
    print(time.diff)
    print(time.diff.over.60)
    Copulas.seleccionadas = Res.Tot
    save(Copulas.seleccionadas, file=paste0('Cop_seleccionadas_Post_Crisis_Vix',today()))
  } # Post_Crisis (2008)
} # Con quiebre en 2008
# Calculo del CoVaR -----------------------------------
#---- Orden ARMA automático
ARMA.Order = ARMA.ORDER.DF(Datos=DATOS_WNA) #Full Sample
if (crisis.break==2008) {
  ARMA.Order.pre  = ARMA.ORDER.DF(Datos=DATOS_PRE)   #Pre Crisis
  ARMA.Order.post = ARMA.ORDER.DF(Datos=DATOS_POST) #Post Crisis
}
#---- Orden ARMA manual para cambios particulares
if (0) {
  N=ncol(DATOS_WNA)
  ARMA.Order = matrix(0, N, 2, dimnames=list(colnames(DATOS_WNA),c('p','q')))
  ARMA.Order['CDS_Brazil',]     = c(0,0)
  ARMA.Order['CDS_Chile',]      = c(0,0)
  ARMA.Order['CDS_China',]      = c(0,0)
  ARMA.Order['CDS_Colombia',]   = c(0,0)
  ARMA.Order['CDS_Indonesia',]  = c(5,0)
  ARMA.Order['CDS_Korea',]      = c(5,0)
  ARMA.Order['CDS_Malaysia',]   = c(5,0)
  ARMA.Order['CDS_Mexico',]     = c(0,0)
  ARMA.Order['CDS_Peru',]       = c(0,0)
  ARMA.Order['CDS_SouthAfrica',]= c(0,0)
  ARMA.Order['CDS_Turkey',]     = c(0,0)
  ARMA.Order['vix',]            = c(0,0)
  
}
#----Estimacion del VaR, CoVaR, DeltaCoVaR
#---- <alpha>: Nivel del VaR (Codicionamiento del CoVaR),  
#---- <beta>: Nivel del CoVaR
#-----<refit.every>: Cada cuantos periodos se reestima el modelo para el rolling
# Serie.1=Name.Series1[6]; copulas=copulas.Sel[6]  
if (1) {
  # CoVaR_data_Vix       = CoVaR_DF(Data=DATOS_WNA, Serie.1=Name.Series1, Serie.2='vix', 
  #                                 copulas=copulas.Sel.vix, alpha=0.05, beta=0.05, plot=TRUE, 
  #                                 COND=CoVaR.type, forecast.type='rolling', refit.every=10, 
  #                                 ARMA.Order=ARMA.Order)
  CoVaR_data_usfcon    = CoVaR_DF(Data=DATOS_WNA, Serie.1=Name.Series1, Serie.2='usfcon', 
                                  copulas=copulas.Sel.usfcon, alpha=0.05, beta=0.05, plot=TRUE, 
                                  COND=CoVaR.type, forecast.type='rolling', refit.every=10, 
                                  ARMA.Order=ARMA.Order)
  #save(CoVaR_data_Vix, file=paste0('CoVaR_Data_Vix',today()))
  save(CoVaR_data_usfcon, file=paste0('CoVaR_data_usfcon',today()))
  load('Cop_seleccionadas_2022-09-17')            #<<<--- USFCON-Archivo donde se guardaron las cop. seleccionadas
  copulas.Sel.usfcon=Copulas.seleccionadas[,'Min.AIC'] 
  # load('Cop_seleccionadas_2022-09-21')            #<<<--- VIX-Archivo donde se guardaron las cop. seleccionadas
  # copulas.Sel.vix = Copulas.seleccionadas[,'Min.AIC'] 
}
if (crisis.break==2008) {
  load('Cop_seleccionadas_Pre_Crisis_Vix2022-10-04')            #<<<--- USFCON-Archivo donde se guardaron las cop. seleccionadas
  copulas.Sel.vix_Pre=Copulas.seleccionadas[,'Min.AIC'] 
  load('Cop_seleccionadas_Post_Crisis_Vix2022-10-04')            #<<<--- VIX-Archivo donde se guardaron las cop. seleccionadas
  copulas.Sel.vix_Post=Copulas.seleccionadas[,'Min.AIC'] 
  #CoVaR Pre_Crisis
  CoVaR_data_vix_pre    = CoVaR_DF(Data=DATOS_PRE, Serie.1=Name.Series1, Serie.2='vix', 
                                  copulas=copulas.Sel.vix_Pre, alpha=0.05, beta=0.05, plot=TRUE, 
                                  COND=CoVaR.type, forecast.type='in sample', ARMA.Order=ARMA.Order.pre)
  
  save(CoVaR_data_vix_pre, file=paste0('CoVaR_data_vix_pre_',today()))
  
  #CoVaR Post_Crisis
  #Ajuste temporal
  CoVaR_data_vix_post    = CoVaR_DF(Data=DATOS_POST, Serie.1=Name.Series1, Serie.2='vix', 
                                      copulas=copulas.Sel.vix_Post, alpha=0.05, beta=0.05, plot=TRUE, 
                                      COND=CoVaR.type, forecast.type='in sample', ARMA.Order=ARMA.Order.post)
  
  save(CoVaR_data_usfcon_post, file=paste0('CoVaR_data_vix_post_',today()))
  
  if (1) {
    CoVaR_pre_post=CoVaR_data_usfcon_pre
    CoVaR_pre_post$CoVaR=rbind(CoVaR_data_usfcon_pre$CoVaR, CoVaR_data_usfcon_post$CoVaR)
    CoVaR_pre_post$CoVaRUp=rbind(CoVaR_data_usfcon_pre$CoVaRUp, CoVaR_data_usfcon_post$CoVaRUp)
    CoVaR_pre_post$VaR=rbind(CoVaR_data_usfcon_pre$VaR, CoVaR_data_usfcon_post$VaR)
    CoVaR_pre_post$VaRUp=rbind(CoVaR_data_usfcon_pre$VaRUp, CoVaR_data_usfcon_post$VaRUp)
    CoVaR_pre_post$DeltaCoVaR=rbind(CoVaR_data_usfcon_pre$DeltaCoVaR, CoVaR_data_usfcon_post$DeltaCoVaR)
    CoVaR_pre_post$DeltaCoVaRUp=rbind(CoVaR_data_usfcon_pre$DeltaCoVaRUp, CoVaR_data_usfcon_post$DeltaCoVaRUp)
    CoVaR_pre_post$CoVaRMedian=rbind(CoVaR_data_usfcon_pre$CoVaRMedian, CoVaR_data_usfcon_post$CoVaRMedian)
    
  } #Se unen las muestras
}
# BigVAR ------------------------------------------------------------------
if (BigVAR==TRUE){
  #----
  # #############################################################################
  #  Para información más detallada ver 'http://www.wbnicholson.com/BigVAR.html'  
  # #############################################################################  
  # Objetivos: 
  # (1) Estimar un modelo VAR penalizado con un numero grande de variables (ej. 20).
  #     En particular cuando <struct='BasicEN'> en la función <constructModel> del paquete BigVAR
  #     Se estima un VAR-ElasticNetwork(Lambda, alpha). alpha=0 (Ridge) y alpha=1 (Lasso)
  #    (1-a) Especificacion del modelo: <constructModel()> 
  #    (1-b) Estimacion del modelo:     <cv.BigVAR()> 
  # (2) Estimar las funciones de impulso-respuesta y GFEVD asociados al VAR de (1).
  #
  # ¿Como se estima Lambda en (1)? -----------------------------------------------------
  # La función <cv.BigVAR> estima Lambda mediante rolling windows (de tamaño <window.size>) minimizando 
  # el MSFE un paso adelante durante el período de entrenamiento(Por defecto: T/3 hasta 2T/3, donde T es el tamaño total de la muestra).
  #
  # Nota: También se puede usar el MSFE <h> pasos adelante.
  #
  # ¿Como se determina alpha en (1)?
  # Se especifica como argumento de la función <constructModel> así: <model.controls = list(alpha = alpha)> 
  # Sí es un escalar se toma dicho valor. 
  # Sí es un vector, la funcion selecciona de ese vector el <alpha> optimo que minimice el MSFE junto a Lambda. 
  # Nota: El valor predeterminado es 1/(k+1).x|
  # ----------------------------------------------------------------------------
  # 
  # (1.a) Los argumentos de entrada DE <constructModel> son:
  #   
  #  <Y>: serie de tiempo multivariada
  #  <p>: orden de rezago máximo
  #  <struct>:     Estructura de la penalización. Las opciones son
  #    Basic:     Penalización por Lazo
  #    BasicEN:   Lasso con Red Elástica (Lasso o Ridge)
  #    lag:       penalización por retraso
  #    OwnOther:  Propia/Otra Penalización
  #    SparseLag: penalización por retraso disperso
  #    SparseOwnOther: Escasa propia/otra penalización por retraso
  #    SCAD:      Desviación absoluta recortada suavemente
  #    MCP:       Penalización cóncava minimax
  #    EFX:       Endógeno-Primero VARX
  #    HLAGC:     HLAG por componentes
  #    HLAGOO:    Propia/Otra HLAG
  #    HLAGELEM:  HLAG elemental
  #    Tapered:   Lasso con penalización por retraso
  #    BGR:       VAR bayesiano como se detalla en Bańbura, Giannone y Reichlin (2010).
  #  Los primeros 7 se pueden aplicar a los modelos VAR y VARX, EFX solo se puede aplicar a los modelos VARX, los 5 restantes solo se aplican a los modelos VAR.
  # 
  #  <gran>: dos opciones para la cuadrícula de parámetros de penalización λ.
  #   - La primera opción controla la profundidad de la cuadrícula lambda (una buena opción predeterminada es 50). 
  #   - La segunda opción controla el número de valores de cuadrícula (un buen valor predeterminado es 10).
  #      Si la cuadrícula no es lo suficientemente profunda, los resultados de los pronósticos pueden ser subóptimos, pero si es demasiado profundo, la rutina puede tardar una cantidad considerable de tiempo en ejecutarse. 
  #      El índice del parámetro de penalización óptimo es monitoreado por <cv.BigVAR()>. Si está en el borde de la cuadrícula, se recomienda volver a ejecutar la función con un parámetro de granularidad mayor. 
  #      Si establece la <ownlambdas=TRUE>, <gran> se usa para proporcionar una cuadrícula definida por el usuario.
  # 
  #   <h>: Horizonte de pronóstico en el que minimizar MSFE (por defecto 1).
  #   <verbose>: Lógico, si es VERDADERO, mostrará una barra de progreso para los procedimientos de validación y evaluación.
  #   <window.size>: Cuanco es > 0: Tamaño de la ventana para el rolling cross validation. 
  #                  Cuando es 0 :  Se utiliza una ventana recursiva(expansiva) para el cross validation. (Este es el valor por defecto).
  #                  
  #   <cv>: Tipo de validación que se utiliza para seleccionar el parámetro de penalización 
  #         Las opciones son:
  #         1. "rolling" (predeterminado) 
  #         2. "LOO", un pseudo procedimiento cv de "Leave One Out Validation" que respeta la dependencia en el tiempo.
  # 
  #   <ownlambdas>: <FALSE> si se utilizan los dos parámetros de <gran> discutidos anteriormente. 
  #                  <TRUE> si el usuario desea sumnistrar un conjunto de valores de <Lambda> para minimizar el MSFE, en cuyo caso los indica en <gran>
  #                 
  # 
  #   <separate_lambdas>: <TRUE>  para usar parámetros de penalización separados para cada serie. Esta opción solo es válida para las estructuras Basic,BasicEN,HLAG,HLAGOO,HLAGELEM,SCAD,MCP.
  #                       <FALSE> Usa el mismo parámetro de penalización para todas las series (opción por defecto)
  #
  #   <model.controls>: Se especifican los parámetros en la lista <model.controls>. Estos parámetros incluyen:
  #     - intercept:  Tenga en cuenta que el intercepto se ajusta por separado y no está sujeto a penalización. ver. "Nicholson, W. B., Matteson, D. S., & Bien, J. (2014). Structured Regularization for Large Vector Autoregressions. 14853, 1–40." Apendix: 9.1
  #                  <TRUE>  se incluye en el modelo (por defecto)
  #                  <FALSE> no se incluye en el modelo
  #     - alpha: ver (¿Como se determina alpha en (1)?)
  #     - linear: <TRUE>  construye una cuadrícula de Lambda lineal(por defecto).
  #               <FALSE> construye una cuadrícula logarítmica-lineal.
  #
  #-------------------------------------#
  #----
  #--- Estimacion del BigVAR para la muestra completa ---# 
  if (0) {
    BigVAR.Model            = Optimal.model(Data=CoVaR_data$CoVaR , alpha=seq(0,1,0.05), auto.lambdas=TRUE, auto.length=30, grid=c(10000,10)) # Estimación del modelo óptimo bajo la opción
    save(BigVAR.Model, file = paste0("BigVAR-Optimo-",today(),".R"))
  }
  load("BigVAR-Optimo-2022-08-31.R") # ¡¡CUIDADO CON LA FECHA!! 
  x11()
  plot(BigVAR.Model$Model) 
  
  # Impulso respuesta y GFEVD ---------------------------------------------------
  GIRF.res              = GIRF_BigVAR(results=BigVAR.Model$Model, data=CoVaR_Data_BV, n.ahead=n.ahead.connectedness,
                                      plot.out='GIRF', ortho=FALSE, std_shock=TRUE, magnitude=1, 
                                      vars.to.shock=colnames(CoVaR_Data_BV[,-1]), vars.to.resp=colnames(CoVaR_Data_BV[,-1]),
                                      x11=FALSE, pdf=FALSE, grid.dims=c(1,1)) #para el calc. adecuado de GIRF, <std_shock> debe ser TRUE.
  GFEVD.res             = GFEVD(Girf=GIRF.res, N.ahead=n.ahead.connectedness, plot.out='GFEVD', vars.to.shock=colnames(CoVaR_Data_BV[,-1]), 
                                vars.to.resp=colnames(CoVaR_Data_BV[,-1]) )
  Static_Network        = Connectedness(GFEVD=GFEVD.res, Title='Net Pairwise directional connectedness', node.size=0.5, 
                                        Q1='95%', Q2='90%', Q3='85%', arrow.size=0.6, layout='circle', n.ahead = n.ahead.connectedness) 
  Dynamic_Network       = Rolling_GFEVD(Data=CoVaR_data$CoVaR, model=BigVAR.Model$Model, structure='BasicEN', window.size=100,
                                        gran=c(10000,10), n.ahead=n.ahead.connectedness, pdf=TRUE, x11=TRUE, alpha=model@alpha, plot=FALSE, BigVAR.window=24) 
  Dynamic_Delta_Network = Rolling_GFEVD(Data=CoVaR_data$DeltaCoVaR, model=BigVAR.Model$Model, structure='BasicEN', window.size=150, 
                                        gran=c(10000,10), n.ahead=n.ahead.connectedness, pdf=TRUE, x11=TRUE, alpha=model@alpha, plot=FALSE, BigVAR.window=24) 
  
}

# BackUp ----------------------------------------------------------------------
save.image(file = paste0(today(),'.RData'))
load(paste0('2022-0-20','.RData'))





