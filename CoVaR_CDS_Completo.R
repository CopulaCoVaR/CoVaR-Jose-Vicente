# Preparación -------------------------------------------------------------
#setwd("C:/Users/maico/OneDrive - Universidad Nacional de Colombia/BanRep/Value at Risk/Highdimensional CoVaR network connectedness/Github-R/Copula-CoVaR-")
#setwd('C:\\Respaldo DLO\\Copula_CoVaR\\R')   #<<<<<--- Carpeta general
#setwd('/Users/lumelo/archivos/Copula_CoVaR/Github-R/CoVaR-Jose-Vicente')
wd = getwd()                                # Carpeta de trabajo adaptable al proyecto, se modifica automáticamente en cada computador.
# ----------------------
Resultados <<- paste0(wd,'/Resultados')     #<<<<<--- Capeta de resultados qué depende deldirectorio de trabajo. 
# ----------------------
data.file    = c('CDS_Data_VIX_EMBI.xlsx','CDS_Data_usfcon.xlsx')[1]              #<<<--- String con el nombre del archivo de datos en la carpeta de trabajo.
# ----------------------
Series.1     = c("CDS_Brazil", "CDS_Chile", "CDS_China", "CDS_Colombia", "CDS_Indonesia", "CDS_Korea", "CDS_Malaysia",   
                 "CDS_Mexico", "CDS_Peru", "CDS_SouthAfrica", "CDS_Turkey")  #<<<--- Columnas de series de análisis
# ----------------------
#Serie.2      = if (data.file=='CDS_Data_VIX.xlsx') 'vix' else 'usfcon'      #<<<--- Columna de serie de referencia.
Serie.2      = c("usfcon", 'vix', 'EMBI_Global')[2]
# ----------------------
Sample.name  = c('Full_Sample', 'Pre_Crisis', 'Crisis', 'Post_Crisis')[4]    #<<<--- Periodo a estimar
# ----------------------
Samples.date = list(pre.crisis  = c('2004-10-08', '2007-07-31'), 
                    crisis      = c('2007-08-01', '2016-11-30'),  #<<<--- Lista para determinar las muestras. 
                    post.crisis = c('2016-12-01'))  
# ----------------------
estimate.copula = c(TRUE, FALSE)[1]
if (Serie.2=='vix') {
  if (Sample.name=='Pre_Crisis')  copula.file = 'Cop_seleccionadas_vix_Pre_Crisis_2022-10-23'
  if (Sample.name=='Crisis')      copula.file = 'Cop_seleccionadas_vix_Crisis_2023-02-28'
  if (Sample.name=='Post_Crisis') copula.file = 'Cop_seleccionadas_vix_Post_Crisis_2022-10-23'
  
}
if (Serie.2=="EMBI_Global") {
  if (Sample.name=='Pre_Crisis')  copula.file = 'Cop_seleccionadas_EMBI_Global_Pre_Crisis_2022-10-25'
  if (Sample.name=='Crisis')      copula.file = 'Cop_seleccionadas_EMBI_Global_Crisis_2022-10-25'
  if (Sample.name=='Post_Crisis') copula.file = 'Cop_seleccionadas_EMBI_Global_Post_Crisis_2022-10-26'
}
# ----------------------
Ext.Reg      = c(TRUE, FALSE)[2]
Series.reg   = c("bsoft", "benergy", "bindustrial", "EMBI_Global", "bcmpi", 'bprec') #<<<--- Columnas de series regresoras externas. 
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
library(stargazer)
library(circlize)
library(car)
library(readxl)
library(fpp3)
library(kldtools)
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
library(plotly)
library(rgl)
library(processx)
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
source('CoVaR_LF.R')
source('garch_roll.R')







# Datos -------------------------------------------------------------------
#--Base de datos de RETORNOS para el periodo de crisis. 
#---- Supuesto_1: la primera col es la FECHA de la base de datos
#---- Supuesto_2: Las col.2 hasta la penultima col. son las series que se modelan w.r.t la serie de la ultima columna (Ej. Ind Acc de varios paises)
#---- Supuesto_3: La ultima col. es la serie que se mantiene en todas copulas (Ej: WTI)

DATOS        = read_xlsx(data.file)                               #<<<---- Base de datos
N.Series1    = ncol(DATOS[,Series.1])                             #Se resta la primera col. (FECHA) y la ultima columna (Ej: WTI)
Name.Serie2  = Serie.2                                            # Nombre de la ultima col. es la serie que se mantiene en todas copulas (Ej: WTI)
Name.Series1 = Series.1                                           # Se quitan la col de  fechas  y se definen las series.1 (i.e. country' Stock returns)
DATOS_xts    = xts(DATOS[,-1], order.by = as.POSIXct(DATOS[[1]])) # (objeto XTS) Se resta la primera col. (FECHA) 
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
if(dif.Serie.1==TRUE| dif.Serie.2==TRUE) DATOS_WNA =DATOS_WNA[-1,]                     #Eliminamos la primera fila de datos (perdidos en la diferencia)

# Regresores externos -----------------------------------------------------
ext.reg.data      = DATOS_xts[,Series.reg]
ext.reg.data.ext1 = diff(ext.reg.data[,Series.reg[4]])
ext.reg.data.ext2 = diff(log(ext.reg.data[,Series.reg[-4]])); 
reg.ext           = cbind(ext.reg.data.ext2[-1,], ext.reg.data.ext1[-1,])
ext.reg           = array(NA, dim=list(nrow(reg.ext), 2, length(Name.Series1)), dimnames = list(c(), c('EMBI','commodity_index'), Name.Series1))
ext.reg[,'EMBI',] = reg.ext[,'EMBI_Global']; sum(is.na(reg.ext[,'EMBI_Global']))

# Arreglo de segundo regresor externo para cada país. 
ext.reg[, 2, "CDS_Brazil"]                                                   = reg.ext[,"bsoft"]
ext.reg[, 2, c("CDS_Chile", "CDS_Korea", "CDS_Peru")]                        = reg.ext[,"bindustrial"]
ext.reg[, 2, c("CDS_China", "CDS_Indonesia", "CDS_Malaysia", "CDS_Turkey")]  = reg.ext[,"bcmpi"]  
ext.reg[, 2, c("CDS_Colombia", "CDS_Mexico")]                                = reg.ext[,"benergy"]
ext.reg[, 2, "CDS_SouthAfrica"]                                              = reg.ext[,"bprec"]

# Muestras -------------------------------------------------------------
Pre_Crisis     = try(window(cbind(DATOS_WNA[,c(Series.1, 'vix')], DATOS_WNA[,Serie.2]), end=Samples.date$pre.crisis[2]));    NAs=sum(is.na(Pre_Crisis)); if (NAs!=0) Pre_Crisis  =na.fill(Pre_Crisis, 'extend')
Crisis         = try(window(cbind(DATOS_WNA[,c(Series.1, 'vix')], DATOS_WNA[,Serie.2]), start=Samples.date$crisis[1], end=Samples.date$crisis[2])); NAs=sum(is.na(Crisis)); if (NAs!=0) Crisis=na.fill(Crisis, 'extend') 
Post_Crisis    = try(window(cbind(DATOS_WNA[,c(Series.1, 'vix')], DATOS_WNA[,Serie.2]), start=Samples.date$post.crisis[1])); NAs=sum(is.na(Post_Crisis));if (NAs!=0) Post_Crisis =na.fill(Post_Crisis, 'extend')

# while(class(Pre_Crisis)=='try-error' | class(Crisis)=='try-error' | class(Post_Crisis)=='try-error') {
#   if (class(Pre_Crisis)=='try-error') Pre_Crisis  = try(window(DATOS_WNA, end=as.Date(Samples.date$pre.crisis[2])+1))
#   if (class(Pre_Crisis)=='try-error') Crisis      = try(window(DATOS_WNA, start=as.Date(Samples.date$crisis[1]+1)))
#   if (class(Pre_Crisis)=='try-error') Post_Crisis = try(window(DATOS_WNA, start=as.Date(Samples.date$post.crisis[1]+1)))
# } #<<<---- Evita fines de semana sí el índice temporal no existe.

ext.reg.pre    = ext.reg[1:nrow(Pre_Crisis), , ]  
ext.reg.crisis = ext.reg[(nrow(Pre_Crisis)+1):(nrow(Pre_Crisis)+nrow(Crisis)),,]
ext.reg.post   = ext.reg[(nrow(Pre_Crisis)+nrow(Crisis)+1):nrow(DATOS_WNA),,]

# ---- Se elige la muestra a trabajar en función del argumento inicial Sample ---- #
EXT.REG = list(Pre_Crisis=ext.reg.pre, Crisis=ext.reg.crisis, Post_Crisis=ext.reg.post)[Sample.name]
Sample  = list(Full_Sample=DATOS_WNA, Pre_Crisis=Pre_Crisis, Crisis=Crisis, Post_Crisis=Post_Crisis)[Sample.name] 
# Estadísticas descriptivas -----------------------------------------------
if(1){
  data.stats = Sample[[1]]
  #data.stats = Pre_Crisis
  Series     = colnames(data.stats) 
  N          = length(Series)
  stats      = matrix(0,12,N,dimnames=list(c('Mean','Max','Min','SD','Skew','Kurt','J-B', 'p-value' ,'ARCH','p-value' ,paste0('Corr.',Name.Serie2),'No.NA'),Series))
  ARMA.Order = ARMA.ORDER.DF(Datos=data.stats)
  for (i in 1:N){
    arima       = try(arima(data.stats[,Series[i]], method = 'ML', order=c(ARMA.Order[Series[i],],0)))
    ARCH        = try(arch.test(arima))
    stats[1,i]  = mean(data.stats[,Series[i]])
    stats[2,i]  = max(data.stats[,Series[i]])
    stats[3,i]  = min(data.stats[,Series[i]])
    stats[4,i]  = sd(data.stats[,Series[i]])
    stats[5,i]  = skewness(data.stats[,Series[i]])
    stats[6,i]  = kurtosis(data.stats[,Series[i]])
    stats[7,i]  = jarque.bera.test(data.stats[,Series[i]])$statistic
    stats[8,i]  = jarque.bera.test(data.stats[,Series[i]])$p.value
    stats[9,i]  = try(ARCH[5,'LM'])
    stats[10,i] = try(ARCH[5,5])
    stats[11,i] = cor(data.stats[,Series[i]],data.stats[,Name.Serie2], method = 'kendall')
    stats[12,i] = sum(is.na(data.stats[,Series[i]]))
  }
  print(t(stats), digits=3)
  write.csv(stats, file=paste0('Descriptive_',names(Sample),'_',Serie.2,'.csv')) # Ajusta el objeto guardado para qué tenga el nombre de la muestra.
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
  save(Copulas.seleccionadas, file=paste0('Cop_seleccionadas_',Name.Serie2,'_',names(Sample),'_',today()))
} 
# Calculo del CoVaR -----------------------------------
#---- Orden ARMA automático
ARMA.data  = list(Full_Sample=DATOS_WNA, Pre_Crisis=Pre_Crisis, Crisis=Crisis, Post_Crisis=Post_Crisis)[names(Sample)]
ARMA.Order = ARMA.ORDER.DF(Datos=ARMA.data[[1]]) 
#---- Carga de las mejores Copulas seleccionadas según fecha.
if(estimate.copula==FALSE) load(copula.file)
copulas.Sel=Copulas.seleccionadas[,'Min.AIC']
#----Estimacion del VaR, CoVaR, DeltaCoVaR y CoVaRmedian
if (Sample.name=='Pre_Crisis') {
  CoVaR_data_Pre_Crisis=CoVaR_DF(Data=Sample[[1]], Serie.1=Name.Series1, Serie.2=Name.Serie2, 
                                 copulas=copulas.Sel, alpha=0.05, beta=0.05, plot=TRUE, 
                                 COND=CoVaR.type, forecast.type=Forecast, refit.every=refit, 
                                 ARMA.Order=ARMA.Order)
  save(CoVaR_data_Pre_Crisis, file=paste0('CoVaR_data_',Name.Serie2,'_', names(Sample),'_',today()))
}
if (Sample.name=='Crisis')     {
  CoVaR_data_Crisis=CoVaR_DF(Data=Sample[[1]], Serie.1=Name.Series1, Serie.2=Name.Serie2, 
                             copulas=copulas.Sel, alpha=0.05, beta=0.05, plot=TRUE, 
                             COND=CoVaR.type, forecast.type=Forecast, refit.every=refit, 
                             ARMA.Order=ARMA.Order)
  save(CoVaR_data_Crisis, file=paste0('CoVaR_data_',Name.Serie2,'_', names(Sample),'_',today()))
  
}
if (Sample.name=='Post_Crisis'){
  CoVaR_data_Post_Crisis=CoVaR_DF(Data=Sample[[1]], Serie.1=Name.Series1, Serie.2=Name.Serie2, 
                                  copulas=copulas.Sel, alpha=0.05, beta=0.05, plot=TRUE,
                                  COND=CoVaR.type, forecast.type=Forecast, refit.every=refit, 
                                  ARMA.Order=ARMA.Order, external.regressors=ifelse(Ext.Reg==TRUE,EXT.REG[[1]], NULL))
  save(CoVaR_data_Post_Crisis, file=paste0('CoVaR_data_',Name.Serie2,'_', names(Sample),'_',today()))
}
# Post-Estimación ---------------------------------------------------------

# Prueba KS
if (1){
  
  if (1) {
    load("CoVaR_data_vix_Pre_Crisis_2023-02-27", verbose = TRUE)
    load("CoVaR_data_vix_Crisis_2023-02-28", verbose = TRUE)
    load("CoVaR_data_vix_Post_Crisis_2023-03-02", verbose = TRUE)
    
  } #Carga de datos VIX sin regresores externos.
  
  if (0) {
    load('CoVaR_data_EMBI_Global_Pre_Crisis_2022-10-25', verbose = TRUE)    
    load('CoVaR_data_EMBI_Global_Crisis_2022-10-26', verbose = TRUE)    
    load('CoVaR_data_EMBI_Global_Post_Crisis_2022-10-26', verbose = TRUE)    
    
  } #Carga de datos EMBI sin regresores externos.
  
  if (0) {
    load('CoVaR_data_usfcon_Pre_Crisis_2022-11-07', verbose = TRUE)    
    load('CoVaR_data_usfcon_Crisis_2022-11-08', verbose = TRUE)    
    load('CoVaR_data_usfcon_Post_Crisis_2022-11-08', verbose = TRUE)    
    
  } #Carga de datos USFCON sin regresores externos.
  test     = c('ksboot', 'created.function')[1]
  samp     = c('Pre_Crisis', 'Crisis', 'Post_Crisis')[3]
  
  if (1) {
    if (samp=='Pre_Crisis') sample = CoVaR_data_Pre_Crisis
    if (samp=='Crisis')     sample = CoVaR_data_Crisis
    if (samp=='Post_Crisis')sample = CoVaR_data_Post_Crisis
  } # Arreglo para la muestra en función de samp
  
  n.rep    = 1000       
  x        = sample$CoVaRUp
  y        = sample$VaRUp
  alt      = c("greater", "two.sided", "less")[3]
  P.values = matrix(NA, nrow=ncol(x), ncol=6, dimnames=list(colnames(x), c('VaR average', 'CoVaR average','SD-VaR', 'SD-CoVaR', 'Statistic','P-value')))
  
  for (i in 1:ncol(x)) {
    Bootstrap.KS                = ksboot(x=as.matrix(x[,i]), y=as.matrix(y[,i]), alternative=alt)
    P.values[i,'P-value']       = Bootstrap.KS$ksboot.pvalue 
    KS.1                        = ks.test(x=as.matrix(x[,i]),y=as.matrix(y[,i]), alternative=alt)
    P.values[i,'Statistic']     = KS.1$statistic
    P.values[i,'CoVaR average'] = mean(x[,i])
    P.values[i,'VaR average']   = mean(y[,i])
    P.values[i,'SD-VaR']        = sd(y[,i])
    P.values[i,'SD-CoVaR']      = sd(x[,i])
  }
  print(P.values, digits=3)
  if (samp=='Pre_Crisis') write.csv(P.values, file=paste('KS_Test_','Pre_Crisis','.csv'))
  if (samp=='Crisis') write.csv(P.values, file=paste('KS_Test_','Crisis','.csv'))
  if (samp=='Post_Crisis') write.csv(P.values, file=paste('KS_Test_','Post_Crisis','.csv'))
} 
# Graficación ECDF
if (0){
  Serie.2=c('VIX', 'EMBI')[2]
  if (Serie.2=='VIX') {
    #<<<<<<< HEAD
    load("CoVaR_data_vix_Pre_Crisis_2022-10-25", verbose=TRUE)
    load("CoVaR_data_vix_Crisis_2022-10-25", verbose=TRUE)
    load("CoVaR_data_vix_Post_Crisis_2022-10-25", verbose=TRUE)

  } #Carga de datos VIX sin regresores externos.
  if (Serie.2=='EMBI'){
    load('CoVaR_data_EMBI_Global_Pre_Crisis_2022-10-25', verbose=TRUE)    
    load('CoVaR_data_EMBI_Global_Crisis_2022-10-26', verbose=TRUE)    
    load('CoVaR_data_EMBI_Global_Post_Crisis_2022-10-26', verbose=TRUE)    
    
  } #Carga de datos EMBI sin regresores externos.
  samp     = c('Pre_Crisis', 'Crisis', 'Post_Crisis')[3]
  titles=c('Brazil', 'Chile', 'China', 'Colombia', 'Indonesia',
           'Korea', 'Malaysia', 'Mexico', 'Peru', 'South Africa', 'Turkey')
  class    = c('pdf','x11')[1]
  if (1) {
    if (samp=='Pre_Crisis') sample = CoVaR_data_Pre_Crisis
    if (samp=='Crisis')     sample = CoVaR_data_Crisis
    if (samp=='Post_Crisis')sample = CoVaR_data_Post_Crisis
  } # Arreglo para la muestra en función de samp
  x        = sample$CoVaRUp
  y        = sample$VaRUp
  if (class=='x11'){
    x11()
    par(mfrow=c(4,3))
  }
  for (i in 1:ncol(x)) {
    ecdx=ecdf(as.matrix(x)[,i])
    ecdy=ecdf(as.matrix(y)[,i])
    if (class=='pdf') pdf(file = paste0(Resultados,'/CDF_',colnames(x)[i],'_',samp,'.pdf'), onefile=FALSE)
    print(plot(ecdx, main=titles[i], col='red', lwd=1))
    print(lines(ecdy, lwd=1))
    if(class=='pdf')dev.off()
  }
}
# Graficación CoVaR
if (0){
  serie2= c('VIX', 'EMBI')[1]
  if (serie2=='VIX') {
    load("CoVaR_data_vix_Pre_Crisis_2022-10-25", verbose=TRUE)
    load("CoVaR_data_vix_Crisis_2022-10-25", verbose=TRUE)
    load("CoVaR_data_vix_Post_Crisis_2022-10-25", verbose=TRUE)
    
  } #Carga de datos VIX sin regresores externos.
  if (serie2=='EMBI') {
    load('CoVaR_data_EMBI_Global_Pre_Crisis_2022-10-25', verbose=TRUE)    
    load('CoVaR_data_EMBI_Global_Crisis_2022-10-26', verbose=TRUE)    
    load('CoVaR_data_EMBI_Global_Post_Crisis_2022-10-26', verbose=TRUE)    
    
  } #Carga de datos EMBI sin regresores externos.
  if (1) {
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
  } # Unión de datos en la misma lista. 
  plot.class = c('Up', 'Down', 'Both')[1]
  titles=c('Brazil', 'Chile', 'China', 'Colombia', 'Indonesia',
           'Korea', 'Malaysia', 'Mexico', 'Peru', 'South Africa', 'Turkey')
  names(titles)=Series.1
  breaks=xts(c('',''), as.Date(c('2007-07-31','2016-11-30')))
  CoVaR_DATA$CoVaR$horizontal_line = 0
  nr    = nrow(CoVaR_DATA$CoVaR)
  shade = cbind(upper = rep(4000, nr), lower = rep(-1, nr))
  shade = xts(shade, ichimoku::index(CoVaR_DATA$CoVaR))
  #Gráficas
  for (i in Series.1){
    if(plot.class=='Up')   range=c(0,max(CoVaR_DATA$CoVaRUp[,i]))
    if(plot.class=='Down') range=c(min(CoVaR_DATA$CoVaR[,i]),0)
    if(plot.class=='Both') range=c(min(CoVaR_DATA$CoVaR[,i]),max(CoVaR_DATA$CoVaRUp[,i]))
    pdf(file = paste0(Resultados,'/Graficas_CoVaR_',serie2,'_',i,'.pdf'), onefile=FALSE)
    print(plotxts(if (plot.class=='Down'|plot.class=='Both')
                { CoVaR_DATA$CoVaR[,i]} else{
                  CoVaR_DATA$CoVaRUp[,i] },
          type="l",col="red", grid.col = NA, 
          ylim=range, xlab="Time", ylab="", lwd=1, main=titles[i],format.labels="%Y", major.ticks = 'years', 
          yaxis.left=TRUE, yaxis.right=TRUE, lty='solid', cex=1.20, cex.axis=0.8)) 
          #, cex.main=2.5, cex.sub=1.5, cex.axis=1.5
    if (plot.class=='Down'|plot.class=='Both') print(lines(CoVaR_DATA$VaR[,i],     col="black", lwd=1,lty='dashed'))
    if (plot.class=='Both')                    print(lines(CoVaR_DATA$CoVaRUp[,i], col="red",   lwd=1,lty='solid'))
    if (plot.class=='Up'|plot.class=='Both')   print(lines(CoVaR_DATA$VaRUp[,i],   col="black", lwd=1,lty='dashed'))
    print(lines(CoVaR_DATA$CoVaR$horizontal_line, col='darkgrey'))
    addEventLines(breaks, col=grey(0.7), lwd=1)
    addPolygon(shade[c('2007-07-31','2016-11-30')], col = grey(0.9), on = -1)
    print(addLegend("top", lwd=2,legend.names = c('CoVaR', 'VaR'), ncol=2, 
                    lty = c('solid','dashed'), col = c('red',   'black')))
    dev.off()
  }
}
# Graficación 3D
if (0){
  Sample  = list(Full_Sample=DATOS_WNA, Pre_Crisis=Pre_Crisis, Crisis=Crisis, Post_Crisis=Post_Crisis)
  
  EMBI_Pre_Crisis_3D=GRAPH.3D(Data=Sample$Pre_Crisis, Serie.1=Series.1, Serie.2=Serie.2, copulas=copulas.Sel,
                              alpha=c(0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9), beta=c(0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9), ARMA.Order=ARMA.Order)
  
  EMBI_Crisis_3D=GRAPH.3D(Data=Sample$Crisis, Serie.1=Series.1, Serie.2=Serie.2, copulas=copulas.Sel,
                              alpha=seq(0.00000001, 0.99999999,0.05), beta=seq(0.00000001, 0.99999999,0.05), ARMA.Order=ARMA.Order, CoVaR.tail='Up')
  
  EMBI_Post_Crisis_3D=GRAPH.3D(Data=Sample$Post_Crisis, Serie.1=Series.1, Serie.2=Serie.2, copulas=copulas.Sel,
                              alpha=seq(0.00000001, 0.99999999,0.05), beta=seq(0.00000001, 0.99999999,0.05), ARMA.Order=ARMA.Order, CoVaR.tail='Up')
}
# Garch individual
if (0){
  sample=c('Full_Sample', 'Pre_Crisis', 'Crisis', 'Post_Crisis')[4] # Selección de la muestra.
  if (1) {
    Sample=list(Full_Sample=DATOS_WNA, Pre_Crisis=Pre_Crisis, Crisis=Crisis, Post_Crisis=Post_Crisis)
    if (sample=='Full_Sample') Data.GARCH=Sample$Full_Sample
    if (sample=='Pre_Crisis') Data.GARCH=Sample$Pre_Crisis
    if (sample=='Crisis') Data.GARCH=Sample$Crisis
    if (sample=='Post_Crisis') Data.GARCH=Sample$Post_Crisis
  } # Se crea el objeto Data.GARCH que contiene los datos para la muestra elegida.
  Series.names = colnames(Data.GARCH)
  ARMA.Order   = ARMA.ORDER.DF(Datos=Data.GARCH) #Se obtiene el orden ARMA de la serie.
  models.table = matrix(NA, 24, length(Series.names), 
                        dimnames = list(c('mu','','ar1','','ar2','',
                                          'ar3','','ar4','','ar5','','omega','','alpha1','','beta1', '','gamma1','','skew','', 'shape',''), Series.names))

  for (s in Series.names) {
    model                   = ugarchspec(variance.model=list(model="gjrGARCH", garchOrder=c(1,1)), 
                                         mean.model=list(armaOrder=ARMA.Order[s,]), distribution.model="sstd")
    model.fit               = ugarchfit(spec=model, data=Data.GARCH[,s], solver='hybrid')
    models.table['mu',s]    = model.fit@fit$matcoef['mu',1]  
    models.table[2,s ]      = model.fit@fit$matcoef['mu',2]
    models.table['ar1',s]   = try(model.fit@fit$matcoef['ar1',1]) 
    models.table[4,s]       = try(model.fit@fit$matcoef['ar1',2])
    models.table['ar2',s]   = try(model.fit@fit$matcoef['ar2',1])  
    models.table[6,s]       = try(model.fit@fit$matcoef['ar2',2])
    models.table['ar3',s]   = try(model.fit@fit$matcoef['ar3',1])  
    models.table[8,s]       = try(model.fit@fit$matcoef['ar3',2])
    models.table['ar4',s]   = try(model.fit@fit$matcoef['ar4',1])  
    models.table[10,s]      = try(model.fit@fit$matcoef['ar4',2])
    models.table['ar5',s]   = try(model.fit@fit$matcoef['ar5',1])  
    models.table[12,s]      = try(model.fit@fit$matcoef['ar5',2])
    models.table['omega',s] = model.fit@fit$matcoef['omega',1]  
    models.table[14,s]      = model.fit@fit$matcoef['omega',2]
    models.table['alpha1',s]= model.fit@fit$matcoef['alpha1',1]  
    models.table[16,s]      = model.fit@fit$matcoef['alpha1',2]
    models.table['beta1',s] = model.fit@fit$matcoef['beta1',1]  
    models.table[18,s]      = model.fit@fit$matcoef['beta1',2]
    models.table['gamma1',s]= model.fit@fit$matcoef['gamma1',1]  
    models.table[20,s]      = model.fit@fit$matcoef['gamma1',2]
    models.table['skew',s]  = model.fit@fit$matcoef['skew',1]  
    models.table[22,s]      = model.fit@fit$matcoef['skew',2]
    models.table['shape',s] = model.fit@fit$matcoef['shape',1]  
    models.table[24,s]      = model.fit@fit$matcoef['shape',2]
    print(models.table, digits=3)
    write.csv(models.table, file=paste0(sample,'_GARCH.csv'))
  }
  
  # Comparación con modelo de la función
  # Prueba.garch=CoVaR.Copula(Serie.1=Series.names[1], Serie.2=Series.names[12], Datos=Data.GARCH, ALPHA=0.05,
  #              BETA=0.05, type=copulas.Sel[1], COND='Less', 
  #              forecast.type='in sample', ARMA.Order=ARMA.Order)
  # Prueba.garch$model.garch.fit@fit$matcoef
}
# Copulas result
Copulas    = matrix(NA,12,11, dimnames=list(rep(c('Type', 'par1', 'par2', 'AIC'),3),Series.1))
if (0) {
  samp = c('Full_Sample', 'Pre_Crisis', 'Crisis', 'Post_Crisis')[3]
  Sample  = list(Full_Sample=DATOS_WNA, Pre_Crisis=Pre_Crisis, Crisis=Crisis, Post_Crisis=Post_Crisis)[samp]
  Sample  = Sample[[1]]
  ARMA.Order = ARMA.ORDER.DF(Datos=Sample) 
  Serie.2=c('vix', 'EMBI_Global')[2]
  if (Serie.2=='vix') {
    if (samp=='Pre_Crisis')  copula.file = 'Cop_seleccionadas_vix_Pre_Crisis_2022-10-23'
    if (samp=='Crisis')      copula.file = 'Cop_seleccionadas_vix_Crisis_2022-10-23'
    if (samp=='Post_Crisis') copula.file = 'Cop_seleccionadas_vix_Post_Crisis_2022-10-23'
  }
  if (Serie.2=="EMBI_Global") {
    if (samp=='Pre_Crisis')  copula.file = 'Cop_seleccionadas_EMBI_Global_Pre_Crisis_2022-10-25'
    if (samp=='Crisis')      copula.file = 'Cop_seleccionadas_EMBI_Global_Crisis_2022-10-26'
    if (samp=='Post_Crisis') copula.file = 'Cop_seleccionadas_EMBI_Global_Post_Crisis_2022-10-26'
  }
  load(copula.file, verbose=TRUE)
  copulas.Sel=Copulas.seleccionadas[,'Min.AIC']
  model.Serie.2     = ugarchspec(variance.model = list(model="sGARCH", garchOrder = c(1,1)), 
                                 mean.model = list(ARMA.Order[Serie.2,]), 
                                 distribution.model = "sstd")
  model.Serie.2.fit = ugarchfit(spec = model.Serie.2, data = Sample[,Serie.2])
  sd_errors2 = residuals(model.Serie.2.fit)/sigma(model.Serie.2.fit)
  u2         = rugarch::pdist(distribution='sstd', q=sd_errors2, skew=coef(model.Serie.2.fit)['skew'], shape=coef(model.Serie.2.fit)['shape'])
  for (i in Series.1) {
    # Models Serie.1
    model.Serie.1     = ugarchspec(variance.model=list(model="gjrGARCH", garchOrder=c(1,1)), 
                                   mean.model=list(armaOrder=ARMA.Order[i,]), distribution.model="sstd")
    model.Serie.1.fit = ugarchfit(spec=model.Serie.1, data=Sample[,i], solver='hybrid')
    sd_errors1 = residuals(model.Serie.1.fit)/sigma(model.Serie.1.fit)
    u1         = rugarch::pdist(distribution='sstd', q=sd_errors1, skew=coef(model.Serie.1.fit)['skew'], shape=coef(model.Serie.1.fit)['shape'])
    u.s        = cbind(u1, u2)
    if (1) {
      cat('Copula fitting started...\n')
      type=copulas.Sel[i]
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
        fit.Copula = BiCopEst(u1,u2, family=14, method='mle') 
        par.Cop.1  = fit.Copula$par
        par.Cop.2  = NULL
        print(fit.Copula)
      }
      if (type=='Rot.Clayton' | type=='Rot.Clayton.Up'){
        fit.Copula = BiCopEst(u1,u2, family=13, method='mle')
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
      if (type=='Rot.BB7' | type=='Rot.BB7.Up'){
        fit.Copula = BiCopEst(u1,u2, family=19, method='mle') 
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
    }
    if (samp=='Pre_Crisis') {
      Copulas[1,i] = copulas.Sel[i]
      Copulas[2,i] = mean(par.Cop.1)
      Copulas[3,i] = mean(par.Cop.2)
      Copulas[4,i] = try(fit.Copula$AIC)

    }
    if (samp=='Crisis')     {
      Copulas[5,i] = copulas.Sel[i]
      Copulas[6,i] = mean(par.Cop.1)
      Copulas[7,i] = mean(par.Cop.2)
      Copulas[8,i] = try(fit.Copula$AIC)
    }
    if (samp=='Post_Crisis'){
      Copulas[9,i] = copulas.Sel[i]
      Copulas[10,i] = mean(par.Cop.1)
      Copulas[11,i] = mean(par.Cop.2)
      Copulas[12,i] = try(fit.Copula$AIC)
    }
  }
}

