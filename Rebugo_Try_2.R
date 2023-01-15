# Preparación -------------------------------------------------------------
setwd("C:/Users/maico/OneDrive - Universidad Nacional de Colombia/BanRep/Value at Risk/Highdimensional CoVaR network connectedness/Github-R/Copula-CoVaR-")
#setwd('C:\\Respaldo DLO\\Copula_CoVaR\\R')   #<<<<<--- Carpeta general
#setwd('/Users/lumelo/archivos/Copula_CoVaR/Github-R/Copula-CoVaR-')
wd = getwd()                                # Carpeta de trabajo adaptable al proyecto, se modifica automáticamente en cada computador.
# ----------------------
Resultados <<- paste0(wd,'/Resultados')     #<<<<<--- Capeta de resultados qué depende deldirectorio de trabajo. 
# ----------------------
Start.date   = Sys.Date()                   # Fecha de inicio, a tener en cuenta para la carga de archivos.
# ----------------------
data.file    = 'Data_Rebugo.xlsx'               #<<<--- String con el nombre del archivo de datos en la carpeta de trabajo.
# ----------------------
Series.1     = c("Brazil", "Russia", "India", "China", "South.Africa", "UK", "US", "EMU") #<<<--- Columnas de series de análisis
# ----------------------
Serie.2      = 'Brent'      #<<<--- Columna de serie de referencia.
# ----------------------
Sample.name  = c('Pre_Crisis', 'Post_Crisis')[2]   #<<<--- Periodo a estimar
# ----------------------
estimate.copula = c(TRUE, FALSE)[2]
if (Sample.name=='Pre_Crisis'){
  Copulas_Pre     = matrix(NA, nrow = length(Series.1[-9]), ncol = 1, dimnames = list(Series.1[-9],'Min.AIC' ))
  Copulas_Pre[,1] = c('Student','Rot.Gumbel', 'Gaussian', 'Gumbel', 'Dyn.Student', 'Dyn.Gaussian', 'Dyn.Gaussian', 'Dyn.Gaussian')
}
if (Sample.name=='Post_Crisis'){
  Copulas_Post    = matrix(NA, nrow = length(Series.1[-9]), ncol = 1, dimnames = list(Series.1[-9],'Min.AIC' ))
  Copulas_Post[,1]= c('Student', 'BB7', 'Dyn.Rot.Gumbel', 'Rot.Gumbel', 'Dyn.Rot.Gumbel', 'Student', 'Rot.Gumbel', 'BB7' )  
}
# ----------------------
Samples.date = list(pre.crisis  = c(0, 451), post.crisis = c(452))   #<<<--- Lista para determinar las muestras.  
# ----------------------
#Series.reg   = c("bsoft", "benergy", "bindustrial", "EMBI_Global", "bcmpi", 'bprec') #<<<--- Columnas de series regresoras externas. 
# ----------------------
CoVaR.type   = c('Equal','Less')[2]         #<<<--- Seleccion del tipo de metodologia para calcular el CoVaR <Equal> para Liu(2022) y <Less> para RU(2016)
# ----------------------
Forecast     = c('in sample', 'rolling')[1] #<<<--- Determina el tipo de pronóstico a llevar a cabo.
if(Forecast=='rolling') refit=10            #<<<--- Sí el pronóstico es rolling, se define cada cuanto reestimar el modelo.       
# ----------------------
log.Serie.1  = c(TRUE, FALSE)[2]            #<<<--- <T> se aplican logaritmos a las <Serie.1>, <F> No se aplican.
log.Serie.2  = c(TRUE, FALSE)[2]            #<<<--- <T> se aplican logaritmos a la <Serie.2>, <F> No se aplica.
dif.Serie.1  = c(TRUE, FALSE)[1]            #<<<--- <T> diferencia las <Serie.1>, <F> deja los datos en nivel. 
dif.Serie.2  = c(TRUE, FALSE)[1]            #<<<--- <T> diferencia la <Serie.2>, <F> deja los datos en nivel. 
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

DATOS        = read_xlsx(data.file)                               #<<<---- Base de datos
N.Series1    = ncol(DATOS[,Series.1])                             #Se resta la primera col. (FECHA) y la ultima columna (Ej: WTI)
Name.Serie2  = Serie.2                                            # Nombre de la ultima col. es la serie que se mantiene en todas copulas (Ej: WTI)
Name.Series1 = Series.1                                           # Se quitan la col de  fechas  y se definen las series.1 (i.e. country' Stock returns)
DATOS_xts    = xts(DATOS[,-1], order.by = as.Date(ichimoku::index(DATOS))) # (objeto XTS) Se resta la primera col. (FECHA) 
NAs=sum(is.na(DATOS_xts))                                         # Numero datos faltantes

# Interpolacion lineal para completar los datos faltantes
if (NAs!=0){
  cat('\nNo. de datos faltantes: ',NAs,'\n') # Numero datos faltantes
  DATOS_WNA = na.fill(DATOS_xts, 'extend')
}else DATOS_WNA=DATOS_xts
#Se transforman las series si se requiere
#Logaritmos
if(log.Serie.1==TRUE) DATOS_WNA[,Name.Series1]=log(DATOS_WNA[,Name.Series1])  # Logaritmo de las Serie.1
if(log.Serie.2==TRUE) DATOS_WNA[,Name.Serie2]=log(DATOS_WNA[,Name.Serie2])  # Logaritmo de la Serie.2                                  # Logaritmos a todas las series.
#Diferencias
if(dif.Serie.1==TRUE) DATOS_WNA[,Name.Series1]=diff(DATOS_WNA[,Name.Series1]) #Diferencia de las  Serie.1
if(dif.Serie.2==TRUE) DATOS_WNA[,Name.Serie2]=diff(DATOS_WNA[,Name.Serie2]) #Diferencia de la  Serie.2
if(dif.Serie.1==TRUE|dif.Serie.2==TRUE) DATOS_WNA =DATOS_WNA[-1,]                     #Eliminamos la primera fila de datos (perdidos en la diferencia)

# Muestras -------------------------------------------------------------
Pre_Crisis     = try(window(cbind(DATOS_WNA[,Series.1], DATOS_WNA[,Serie.2]), end=as.Date(Samples.date$pre.crisis[2])));    NAs=sum(is.na(Pre_Crisis)); if (NAs!=0) Pre_Crisis  =na.fill(Pre_Crisis, 'extend')
Post_Crisis    = try(window(cbind(DATOS_WNA[,Series.1], DATOS_WNA[,Serie.2]), start=as.Date(Samples.date$post.crisis[1]))); NAs=sum(is.na(Post_Crisis));if (NAs!=0) Post_Crisis =na.fill(Post_Crisis, 'extend')

# ---- Se elige la muestra a trabajar en función del argumento inicial Sample ---- #
Sample  = list(Full_Sample=DATOS_WNA, Pre_Crisis=Pre_Crisis, Post_Crisis=Post_Crisis)[Sample.name] 
# Estadísticas descriptivas -----------------------------------------------
if(0){
  data.stats = Sample[[1]]
  Series     = colnames(data.stats)
  N          = length(Series)
  stats      = matrix(0,10,N,dimnames=list(c('Mean','Max','Min','SD','Skew','Kurt','J-B','ARCH' ,paste0('Corr.',Name.Serie2),'No.NA'),Series))
  ARMA.Order = ARMA.ORDER.DF(Datos=DATOS_WNA)
  for (i in 1:N){
    arima       = arima(DATOS_WNA[,Series[i]], method = 'ML', order=c(ARMA.Order[Series[i],],0))
    ARCH        = arch.test(arima)
    stats[1,i]  = mean(DATOS_WNA[,Series[i]])
    stats[2,i]  = max(DATOS_WNA[,Series[i]])
    stats[3,i]  = min(DATOS_WNA[,Series[i]])
    stats[4,i]  = sd(DATOS_WNA[,Series[i]])
    stats[5,i]  = skewness(DATOS_WNA[,Series[i]])
    stats[6,i]  = kurtosis(DATOS_WNA[,Series[i]])
    stats[7,i]  = jarque.bera.test(DATOS_WNA[,Series[i]])$p.value
    stats[8,i]  = ARCH[5,'LM']
    stats[9,i]  = cor(DATOS_WNA[,Series[i]],DATOS_WNA[,Name.Serie2], method = 'kendall')
    stats[10,i] = sum(is.na(DATOS_WNA[,Series[i]]))
  }
  print(t(stats), digits=3)
  write.csv(stats, file=paste0('Descriptive_',names(Sample),'.csv')) # Ajusta el objeto guardado para qué tenga el nombre de la muestra.
}
# Graficación -------------------------------------------------------------
if(0){
  Data.graph= DATOS_WNA # Gráfica de las series en nivel. 
  #Data.graph=list(Full_Sample=DATOS_WNA, Pre_Crisis=Pre_Crisis, Crisis=Crisis, Post_Crisis=Post_Crisis)[names(Sample)] # Gráfica de series transformadas y submuestras. 
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
# Selección de Copula----------------
if(estimate.copula==TRUE){
  Data.copula = Sample[[1]]
  time.ini    = Sys.time()
  Serie2      = Name.Serie2 #Ej: WTI
  Res.Tot     = array(NA,dim=c(N.Series1,2), 
                      dimnames=list(Name.Series1, c('Min.AIC', 'Min.BIC')))
  ii.n      = 0 
  for (ii in Name.Series1){
    print(ii)
    ii.n            = ii.n + 1
    res             = Copula_seleccion(Datos=Data.copula, arma.order=NULL, GARCH.model="gjrGARCH", serie.1=ii, serie.2=Serie2, CoVaR.type=CoVaR.type)
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
  Copulas.seleccionadas = Res.Tot; Copulas.date=Sys.Date()
  save(Copulas.seleccionadas, file=paste0('Cop_seleccionadas_',names(Sample),'_',today()))
} 
# Calculo del CoVaR -----------------------------------
ARMA.Order = matrix(0, length(Name.Series1)+length(Name.Serie2), 2, dimnames=list(c(Name.Series1,Name.Serie2),c('p','q')))
ARMA.Order['Brazil',]  = c(0,0)
ARMA.Order['Russia',]  = c(0,0)
ARMA.Order['India',]   = c(0,0)
ARMA.Order['China',]   = c(0,0)
ARMA.Order['South.Africa',]  = c(0,0)
ARMA.Order['UK',] = c(0,0)
ARMA.Order['US',] = c(0,0)
ARMA.Order['EMU',]  = c(0,0)
ARMA.Order['Brent',]  = c(0,0)
#----Estimacion del VaR, CoVaR, DeltaCoVaR y CoVaRmedian
if (names(Sample)=='Pre_Crisis')copulas.Sel=as.vector(Copulas_Pre)
if (names(Sample)=='Post_Crisis')copulas.Sel=as.vector(Copulas_Post)
names(copulas.Sel)=Name.Series1
#EXT.REG[[1]]
if (Sample.name=='Pre_Crisis') {
  CoVaR_data_Pre_Crisis=CoVaR_DF(Data=Sample[[1]], Serie.1=Name.Series1, Serie.2=Name.Serie2, 
                                 copulas=copulas.Sel, alpha=0.05, beta=0.05, plot=TRUE, 
                                 COND=CoVaR.type, forecast.type=Forecast, refit.every=refit, 
                                 ARMA.Order=ARMA.Order, external.regressors=NULL)
  save(CoVaR_data_Pre_Crisis, file=paste0('CoVaR_data_',Name.Serie2,'_', names(Sample),'_',today()))
}
if (Sample.name=='Post_Crisis') {
  CoVaR_data_Post_Crisis=CoVaR_DF(Data=Sample[[1]], Serie.1=Name.Series1, Serie.2=Name.Serie2, 
                                  copulas=copulas.Sel, alpha=0.05, beta=0.05, plot=TRUE,
                                  COND=CoVaR.type, forecast.type=Forecast, refit.every=refit, 
                                  ARMA.Order=ARMA.Order, external.regressors=NULL)
  save(CoVaR_data_Post_Crisis, file=paste0('CoVaR_data_',Name.Serie2,'_', names(Sample),'_',today()))
}

if (1){
  n.rep    = 1000        
  x        = CoVaR_data$CoVaR
  y        = CoVaR_data$VaR
  P.values = matrix(NA, nrow=ncol(x), ncol=3, dimnames=list(colnames(x), c('VaR average', 'CoVaR average','P-value')))
  for (i in 1:ncol(x)) {
    Bootstrap.KS                = KS.bootstrap(x=as.matrix(x[,i]), y=as.matrix(y[,i]), 
                                               n.boot=n.rep, density.plot=F)
    P.values[i,'P-value']       = Bootstrap.KS$p.value
    P.values[i,'CoVaR average'] = mean(x[,i])
    P.values[i,'VaR average']   = mean(y[,i])
  }
  print(P.values)
  save(P.values, file=paste0('KS_Test_',Name.Serie2,'_', names(Sample),'_',today()))
}

# Se unen las muestras (Manual)
if (0) {
  load("C:\\Users\\maico\\OneDrive - Universidad Nacional de Colombia\\BanRep\\Value at Risk\\Highdimensional CoVaR network connectedness\\Github-R\\Copula-CoVaR-\\CoVaR_data_Brent_Pre_Crisis_2022-10-19")
  load("C:\\Users\\maico\\OneDrive - Universidad Nacional de Colombia\\BanRep\\Value at Risk\\Highdimensional CoVaR network connectedness\\Github-R\\Copula-CoVaR-\\CoVaR_data_Brent_Post_Crisis_2022-10-19")
  #load("C:\\Users\\maico\\OneDrive - Universidad Nacional de Colombia\\BanRep\\Value at Risk\\Highdimensional CoVaR network connectedness\\Github-R\\Copula-CoVaR-\\CoVaR_data_Brent_Pre_Crisis_2022-10-20")
  #load("C:\\Users\\maico\\OneDrive - Universidad Nacional de Colombia\\BanRep\\Value at Risk\\Highdimensional CoVaR network connectedness\\Github-R\\Copula-CoVaR-\\CoVaR_data_Brent_Post_Crisis_2022-10-20")
  CoVaR_DATA=list()
  CoVaR_pre=CoVaR_data_Pre_Crisis 
  CoVaR_Post=CoVaR_data_Post_Crisis 
  CoVaR_DATA$CoVaR=rbind(CoVaR_pre$CoVaR, CoVaR_Post$CoVaR)
  CoVaR_DATA$CoVaRUp=rbind(CoVaR_pre$CoVaRUp, CoVaR_Post$CoVaRUp)
  CoVaR_DATA$VaR=rbind(CoVaR_pre$VaR, CoVaR_Post$VaR)
  CoVaR_DATA$VaRUp=rbind(CoVaR_pre$VaRUp, CoVaR_Post$VaRUp)
  CoVaR_DATA$DeltaCoVaR=rbind(CoVaR_pre$DeltaCoVaR, CoVaR_Post$DeltaCoVaR)
  CoVaR_DATA$DeltaCoVaRUp=rbind(CoVaR_pre$DeltaCoVaRUp, CoVaR_Post$DeltaCoVaRUp)
  CoVaR_DATA$CoVaRMedian=rbind(CoVaR_pre$CoVaRMedian, CoVaR_Post$CoVaRMedian)
  breaks=xts('', as.Date(451))
  CoVaR_DATA$CoVaR$horizontal_line = 0
  #Gráficas
  for (i in Name.Series1){
    pdf(file = paste0(Resultados,'/Graficas_CoVaR_',i,'.pdf'), onefile=FALSE)
    print(plot.xts(CoVaR_DATA$CoVaR[,i],type="l",col="blue", grid.col = NA, ylim=c(min(CoVaR_DATA$CoVaR[,i]),max(CoVaR_DATA$CoVaRUp[,i])),xlab="Time",
                   ylab="", lwd=1.5, main=i, format.labels="%Y", major.ticks = 'years', 
                   yaxis.left=TRUE, yaxis.right=TRUE, lty='solid'))
    print(lines(CoVaR_DATA$VaR[,i],         col="black", lwd=1.5,lty='dashed'))
    print(lines(CoVaR_DATA$CoVaRUp[,i],     col="red",   lwd=1.5,lty='dotted'))
    print(lines(CoVaR_DATA$VaRUp[,i],       col="green", lwd=1.5,lty='dashed'))
    print(lines(CoVaR_DATA$CoVaR$horizontal_line, col='black'))
    addEventLines(breaks, col='orange', lwd=1.5)
    print(addLegend("topright", lwd=2,legend.names = c('CoVaR', 'VaR', 'CoVaRUp', 'VaRUp'), 
                    lty = c('solid','dashed','dotted', 'dashed'), col = c('blue',   'black', 'red', 'green')))
    dev.off()
  }
}

x11()
par(mfrow=c(5,2))
for (i in colnames(DATOS_WNA)) {
  print(plot.xts(DATOS_WNA[,i], main=i, col=rand_color(n=1)))
}
