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
data.file    = 'Data_Rebugo.xlsx'              #<<<--- String con el nombre del archivo de datos en la carpeta de trabajo.
# ----------------------
Series.1     = c("Brazil", "Russia", "India", "China", "South Africa", "UK", "US",   
                 "EMU", "BRICS") #<<<--- Columnas de series de análisis
# ----------------------
Serie.2      = 'Brent'     #<<<--- Columna de serie de referencia.
# ----------------------
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
# Cargamos y modificamos datos (Weekly data) --------------------------------------------
# Weekly data from 7Jan2000 - 19Dec2014
Full_Sample = ts(read_xlsx(data.file))
#Cols: c("date", "Brazil","Russia", "India", "China", "South Africa", "UK", "US", "EMU", "BRICS", "Brent") 
Series = colnames(Full_Sample)[-1] #w/o <date>
Series = Series[-9]
Pre_Crisis    = window(ts(Full_Sample), start = 1, end=451) #2000 - 2008
Post_Crisis   = window(ts(Full_Sample), start = 452, end=780) #2008-2014
Post_Crisis_Ch= window(ts(Full_Sample), start = 452, end=735) #2008-2014Ene31, Solo para China
N = length(Series); 
T = matrix(0,4,1,dimnames=list(c('Full_Sample','Pre_Crisis','Post_Crisis','Post_Crisis_Ch'), 'Sample Size'))
T['Full_Sample',]= nrow(Full_Sample); T['Pre_Crisis',]    =nrow(Pre_Crisis);
T['Post_Crisis',]=nrow(Post_Crisis);  T['Post_Crisis_Ch',]=nrow(Post_Crisis_Ch)
#Full_Sample_BR_BRENT = ts(Full_Sample[,c("Brazil", "Brent")])
#Pre_Crisis_BR_BRENT  = ts(Pre_Crisis[,c("Brazil", "Brent")])

ARMA.Order = matrix(0, N, 2, dimnames=list(Series,c('p','q')))
ARMA.Order['Brazil',]  = c(0,0)
ARMA.Order['Russia',]  = c(1,1)
ARMA.Order['India',]   = c(1,1)
ARMA.Order['China',]   = c(0,0)
ARMA.Order['South Africa',]  = c(0,0)
ARMA.Order['UK',] = c(0,0)
ARMA.Order['US',] = c(1,0)
ARMA.Order['EMU',]  = c(0,0)
ARMA.Order['Brent',]  = c(0,0)


# Selección de Copula----------------
  Copulas_Pre     = matrix(NA, nrow = length(Series[-9]), ncol = 1, dimnames = list(Series[-9],'Min.AIC' ))
  Copulas_Post    = matrix(NA, nrow = length(Series[-9]), ncol = 1, dimnames = list(Series[-9],'Min.AIC' ))
  Copulas_Pre[,1] = c('Student','Rot.Gumbel', 'Gaussian', 'Gumbel', 'Dyn.Student', 'Dyn.Gaussian', 'Dyn.Gaussian', 'Dyn.Gaussian')
  Copulas_Post[,1]= c('Student', 'BB7', 'Dyn.Rot.Gumbel', 'Rot.Gumbel', 'Dyn.Rot.Gumbel', 'Student', 'Rot.Gumbel', 'BB7' )  
# Calculo del CoVaR -----------------------------------
#----Estimacion del VaR, CoVaR, DeltaCoVaR y CoVaRmedian
#EXT.REG[[1]]
if (Sample.name=='Pre_Crisis') {
  CoVaR_data_Pre_Crisis=CoVaR_DF(Data=xts(Pre_Crisis[,-1], order.by = as.Date(Pre_Crisis[,'date'])), Serie.1=Series[-9], Serie.2='Brent', 
                                 copulas=Copulas_Pre, alpha=0.05, beta=0.05, plot=TRUE, 
                                 COND='Less', forecast.type='in sample', refit.every=refit, 
                                 ARMA.Order=ARMA.Order, external.regressors=NULL)
  save(CoVaR_data_Pre_Crisis, file=paste0('CoVaR_data_',Name.Serie2,'_', names(Sample),'_',today()))
}
if (Sample.name=='Crisis') {
  CoVaR_data_Crisis=CoVaR_DF(Data=Sample[[1]], Serie.1=Name.Series1, Serie.2=Name.Serie2, 
                             copulas=copulas.Sel, alpha=0.05, beta=0.05, plot=TRUE, 
                             COND=CoVaR.type, forecast.type=Forecast, refit.every=refit, 
                             ARMA.Order=ARMA.Order, external.regressors=NULL)
  save(CoVaR_data_Crisis, file=paste0('CoVaR_data_',Name.Serie2,'_', names(Sample),'_',today()))
  
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
  # Con regresores externos de la reunión del 15 de octubre
  if (0) {
    load("C:\\Users\\maico\\OneDrive - Universidad Nacional de Colombia\\BanRep\\Value at Risk\\Highdimensional CoVaR network connectedness\\Github-R\\Copula-CoVaR-\\Resultados\\usfcon\\10-15-2022-ext-reg\\CoVaR_data_usfcon_Pre_Crisis_2022-10-14")
    load("C:\\Users\\maico\\OneDrive - Universidad Nacional de Colombia\\BanRep\\Value at Risk\\Highdimensional CoVaR network connectedness\\Github-R\\Copula-CoVaR-\\Resultados\\usfcon\\10-15-2022-ext-reg\\CoVaR_data_usfcon_Crisis_2022-10-14")
    load("C:\\Users\\maico\\OneDrive - Universidad Nacional de Colombia\\BanRep\\Value at Risk\\Highdimensional CoVaR network connectedness\\Github-R\\Copula-CoVaR-\\Resultados\\usfcon\\10-15-2022-ext-reg\\CoVaR_data_usfcon_Post_Crisis_2022-10-14")
    
  }
  # Con regresores externos y modificación del 17 de octubre
  if (0) {
    load("C:\\Users\\maico\\OneDrive - Universidad Nacional de Colombia\\BanRep\\Value at Risk\\Highdimensional CoVaR network connectedness\\Github-R\\Copula-CoVaR-\\Resultados\\usfcon\\10-17-2022-ext-reg - Changed\\CoVaR_data_usfcon_Pre_Crisis_2022-10-17")
    load("C:\\Users\\maico\\OneDrive - Universidad Nacional de Colombia\\BanRep\\Value at Risk\\Highdimensional CoVaR network connectedness\\Github-R\\Copula-CoVaR-\\Resultados\\usfcon\\10-17-2022-ext-reg - Changed\\CoVaR_data_usfcon_Crisis_2022-10-17")
    load("C:\\Users\\maico\\OneDrive - Universidad Nacional de Colombia\\BanRep\\Value at Risk\\Highdimensional CoVaR network connectedness\\Github-R\\Copula-CoVaR-\\Resultados\\usfcon\\10-17-2022-ext-reg - Changed\\CoVaR_data_usfcon_Post_Crisis_2022-10-17")
    
  }
  # Sin regresores externos 18 de octubre
  if (0) {
    load("C:\\Users\\maico\\OneDrive - Universidad Nacional de Colombia\\BanRep\\Value at Risk\\Highdimensional CoVaR network connectedness\\Github-R\\Copula-CoVaR-\\Resultados\\usfcon\\10-18-2022-no-ext-reg\\CoVaR_data_usfcon_Pre_Crisis_2022-10-18")
    load("C:\\Users\\maico\\OneDrive - Universidad Nacional de Colombia\\BanRep\\Value at Risk\\Highdimensional CoVaR network connectedness\\Github-R\\Copula-CoVaR-\\Resultados\\usfcon\\10-18-2022-no-ext-reg\\CoVaR_data_usfcon_crisis_2022-10-18")
    load("C:\\Users\\maico\\OneDrive - Universidad Nacional de Colombia\\BanRep\\Value at Risk\\Highdimensional CoVaR network connectedness\\Github-R\\Copula-CoVaR-\\Resultados\\usfcon\\10-18-2022-no-ext-reg\\CoVaR_data_usfcon_Post_Crisis_2022-10-18")
  }
  CoVaR_DATA=list()
  CoVaR_pre=CoVaR_data_Pre_Crisis 
  CoVaR_Crisis=CoVaR_data_Crisis 
  CoVaR_Post=CoVaR_data_Post_Crisis 
  CoVaR_DATA$CoVaR=rbind(CoVaR_pre$CoVaR, CoVaR_Crisis$CoVaR, CoVaR_Post$CoVaR)
  CoVaR_DATA$CoVaRUp=rbind(CoVaR_pre$CoVaRUp, CoVaR_Crisis$CoVaRUp, CoVaR_Post$CoVaRUp)
  CoVaR_DATA$VaR=rbind(CoVaR_pre$VaR, CoVaR_Crisis$VaR, CoVaR_Post$VaR)
  CoVaR_DATA$VaRUp=rbind(CoVaR_pre$VaRUp, CoVaR_Crisis$VaRUp, CoVaR_Post$VaRUp)
  CoVaR_DATA$DeltaCoVaR=rbind(CoVaR_pre$DeltaCoVaR, CoVaR_Crisis$DeltaCoVaR, CoVaR_Post$DeltaCoVaR)
  CoVaR_DATA$DeltaCoVaRUp=rbind(CoVaR_pre$DeltaCoVaRUp, CoVaR_Crisis$DeltaCoVaRUp, CoVaR_Post$DeltaCoVaRUp)
  CoVaR_DATA$CoVaRMedian=rbind(CoVaR_pre$CoVaRMedian, CoVaR_Crisis$CoVaRMedian, CoVaR_Post$CoVaRMedian)
  breaks=xts(letters[1:2], as.Date(c('2007-07-31','2016-11-30')))
  #Gráficas
  for (i in Name.Series1){
    pdf(file = paste0(Resultados,'/Graficas_CoVaR_',i,'.pdf'), onefile=FALSE)
    print(plot.xts(CoVaR_DATA$CoVaR[,i],type="l",col="red", grid.col = NA, ylim=c(min(CoVaR_DATA$CoVaR[,i]),max(CoVaR_DATA$CoVaRUp[,i])),xlab="Time",
                   ylab="", lwd=1, main=i, format.labels="%Y", major.ticks = 'years', 
                   yaxis.left=TRUE, yaxis.right=TRUE, lty='dotted'))
    print(lines(CoVaR_DATA$VaR[,i],         col="black", lwd=1,lty='dashed'))
    print(lines(CoVaR_DATA$CoVaRUp[,i],     col="red",   lwd=1,lty='solid'))
    print(lines(CoVaR_DATA$VaRUp[,i],       col="black", lwd=1,lty='dashed'))
    addEventLines(breaks, col='blue', lwd=3)
    print(addLegend("topright", lwd=2,legend.names = c('CoVaR', 'VaR'), 
                    lty = c('dashed','solid'), col = c('red',   'black')))
    dev.off()
  }
}
