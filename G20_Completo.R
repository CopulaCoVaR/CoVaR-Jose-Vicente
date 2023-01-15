#setwd("C:/Users/maico/OneDrive - Universidad Nacional de Colombia/BanRep/Value at Risk/Highdimensional CoVaR network connectedness/Github-R/Copula-CoVaR-")
#setwd('C:\\Respaldo DLO\\Copula_CoVaR\\R')   #<<<<<--- Carpeta general
setwd('/Users/lumelo/archivos/Copula_CoVaR/Github-R/Copula-CoVaR-')
wd = getwd()                                # Carpeta de trabajo adaptable al proyecto, se modifica automáticamente en cada computador.
Resultados <<- paste0(wd,'/Resultados')     #<<<<<--- Capeta de resultados qué depende deldirectorio de trabajo. 
CoVaR.type  = c('Equal','Less')[1]          #<<<--- Seleccion del tipo de metodologia para calcular el CoVaR <Equal> para Liu(2022) y <Less> para RU(2016)
BigVAR      = c(TRUE, FALSE)[2]             #<<<--- <F> si solo se estima el CoVaR, <T> si ademas se realiza con VAR y conectividad
# Paquetes ----------------------------------------------------------------
library(ichimoku) # Para transformar de xts a df
library(xts)
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

 G20_Crisis = read_xlsx('G20_crisis.xlsx')  #<<<---- Base de datos

 sum(is.na(G20_Crisis)) # Numero datos faltantes
 N.Series1    = ncol(G20_Crisis)- 1 -1                  #Se resta la primera col. (FECHA) y la ultima columna (Ej: WTI)
 Name.Serie2  = colnames(G20_Crisis)[ncol(G20_Crisis)]  #Nombre de la ultima col. es la serie que se mantiene en todas copulas (Ej: WTI)
 Name.Series1 = colnames(G20_Crisis[,-1])[1:N.Series1]  #Se quitan la col de  fechas  y se definen las series.1 (i.e. country' Stock returns)

# As XTS
G20_Crisis_xts = xts(G20_Crisis[,-1], order.by = as.POSIXct(G20_Crisis[[1]])) #Se resta la primera col. (FECHA) 
# Interpolacion lineal para completar los datos faltantes (para datos interiores, en los extremos toma el valor mas cercano)
G20_Crisis_WNA = na.fill(G20_Crisis_xts, 'extend')
sum(is.na(G20_Crisis_WNA)) # Numero datos faltantes
G20_Crisis_WNA = G20_Crisis_WNA*100 #Multiplica los retornos por cien. Quedan en terminos porcentuales
# Estadísticas descriptivas -----------------------------------------------
if(1){
  Series = colnames(G20_Crisis_WNA)
  N      = length(Series)
  stats  = matrix(0,7,N,dimnames=list(c('Mean','Max','Min','SD','Skew','Kurt','J-B'),Series))
  for (i in 1:N){
    stats[1,i] = mean(G20_Crisis_WNA[,Series[i]])
    stats[2,i] = max(G20_Crisis_WNA[,Series[i]])
    stats[3,i] = min(G20_Crisis_WNA[,Series[i]])
    stats[4,i] = sd(G20_Crisis_WNA[,Series[i]])
    stats[5,i] = skewness(G20_Crisis_WNA[,Series[i]])
    stats[6,i] = kurtosis(G20_Crisis_WNA[,Series[i]])
    stats[7,i] = jarque.bera.test(G20_Crisis_WNA[,Series[i]])$statistic
  }
  print(t(stats), digits=3)
}
# Selección de Copula -----------------------------------------------------
if(0){
  time.ini  = Sys.time()
  Serie2    = Name.Serie2 #Ej: WTI
  Res.Tot   = array(NA,dim=c(N.Series1,2), 
                    dimnames=list(Name.Series1, c('Min.AIC', 'Min.BIC')))
  ii.n      = 0 
  for (ii in Name.Series1){
    print(ii)
    ii.n            = ii.n + 1
    res             = Copula_seleccion(Datos=G20_Crisis_WNA, arma.order=NULL, GARCH.model="gjrGARCH", serie.1=ii, serie.2=Serie2, CoVaR.type=CoVaR.type)
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
  
}
# Creación del objeto con las series de CoVaR -----------------------------------
if (1) {
  #---Seleccion de las mejores copulas 
  load('Cop_seleccionadas_2022-09-15')  #<<<--- Archivo donde se guardaron las cop. seleccionadas
  copulas.Sel = Copulas.seleccionadas[,'Min.AIC'] #Falta corregir los nombres de la copulas!!!
  #for(i12 in 1:length(copulas.Sel)) if (copulas.Sel[i12]=='Dyn.t') copulas.Sel[i12]='Dyn.Student'

  #----Estimacion del VaR, CoVaR, DeltaCoVaR
  #---- <alpha>: Nivel del VaR (Codicionamiento del CoVaR),  
  #---- <beta>: Nivel del CoVaR
  #-----<refit.every>: Cada cuantos periodos se reestima el modelo para el rolling
  time.ini   = Sys.time() 
  CoVaR_data = CoVaR_DF(Data=G20_Crisis_WNA, Serie.1=Name.Series1, Serie.2=Name.Serie2, copulas=copulas.Sel, alpha=0.05, beta=0.05, plot=TRUE, COND=CoVaR.type, forecast.type = 'rolling', refit.every = 20) 
  time.diff  = Sys.time() - time.ini
  time.diff.over.60 = (Sys.time() - time.ini)/60 
  print(time.diff)
  print(time.diff.over.60)
  save(CoVaR_data, file=paste0('CoVaR_Data',today()))
  
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
  if (0) {
    BigVAR.Model            = Optimal.model(Data=CoVaR_data$CoVaR , alpha=seq(0,1,0.05), auto.lambdas=TRUE, auto.length=30, grid=c(10000,10)) # Estimación del modelo óptimo bajo la opción
    save(BigVAR.Model, file = paste0("BigVAR-Optimo-",today(),".R"))
  }
  load("BigVAR-Optimo-2022-08-31.R") # ¡¡CUIDADO CON LA FECHA!! 
  x11()
  plot(BigVAR.Model$Model)

# Impulso respuesta y GFEVD ---------------------------------------------------
n.ahead = 20  #<<<--- Horizonte asociado al GIRF y al GFEVD  
GIRF.res              = GIRF_BigVAR(results=BigVAR.Model$Model, data=CoVaR_Data_BV, n.ahead=n.ahead,
                                    plot.out='GIRF', ortho=FALSE, std_shock=TRUE, magnitude=1, 
                                    vars.to.shock=colnames(CoVaR_Data_BV[,-1]), vars.to.resp=colnames(CoVaR_Data_BV[,-1]),
                                    x11=FALSE, pdf=FALSE, grid.dims=c(1,1))
GFEVD.res             = GFEVD(Girf=GIRF.res, N.ahead=n.ahead, plot.out='GFEVD', vars.to.shock=colnames(CoVaR_Data_BV[,-1]), 
                                    vars.to.resp=colnames(CoVaR_Data_BV[,-1]) )
Static_Network        = Connectedness(GFEVD=GFEVD.res, Title='Net Pairwise directional connectedness', node.size=0.5, 
                                      Q1='95%', Q2='90%', Q3='85%', arrow.size=0.6, layout='circle') 
Dynamic_Network       = Rolling_GFEVD(Data=CoVaR_data$CoVaR, model=BigVAR.Model$Model, structure='BasicEN', window.size=100,
                                      gran=c(10000,10), n.ahead=n.ahead, pdf=TRUE, x11=TRUE, alpha=model@alpha, plot=FALSE, BigVAR.window=24) 
Dynamic_Delta_Network = Rolling_GFEVD(Data=CoVaR_data$DeltaCoVaR, model=BigVAR.Model$Model, structure='BasicEN', window.size=150, 
                                      gran=c(10000,10), n.ahead=n.ahead, pdf=TRUE, x11=TRUE, alpha=model@alpha, plot=FALSE, BigVAR.window=24) 

}

# BackUp ----------------------------------------------------------------------

save.image(file = paste0(today(),'.RData'))
load(paste0('2022-07-19','.RData'))
