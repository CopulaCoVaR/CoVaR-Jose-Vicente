#rango de las gráficas
if (0){
  if(CoVaR.plot=='Up'  & VaR.plot==TRUE) rango=range(CoVaRUp.series[,i],VaRUp.series[,i],DeltaCoVaRUp.series[,i])
  if(CoVaR.plot=='Up'  & VaR.plot==FALSE)rango=range(CoVaRUp.series[,i],DeltaCoVaRUp.series[,i])
  if(CoVaR.plot=='Down'& VaR.plot==TRUE) rango=range(CoVaR.series[,i],VaR.series[,i],DeltaCoVaR.series[,i])
  if(CoVaR.plot=='Down'& VaR.plot==FALSE)rango=range(CoVaR.series[,i],DeltaCoVaR.series[,i])
  if(CoVaR.plot=='Both'& VaR.plot==TRUE) rango=range(CoVaRUp.series[,i],CoVaR.series[,i],VaR.series[,i],VaRUp.series[,i],DeltaCoVaR.series[,i],DeltaCoVaRUp.series[,i])
  if(CoVaR.plot=='Both'& VaR.plot==FALSE)rango=range(CoVaRUp.series[,i],CoVaR.series[,i],DeltaCoVaR.series[,i],DeltaCoVaRUp.series[,i])
  
}
# --------------------------- Graficas --------------------------------
# Función que permite graficar los resultados de <CoVaR_DF> en distintos formatos
# de forma agrupada o individual. 
# --------------------------- Argumentos de entrada -----------------------
# <CoVaR_data>: Objeto de salida de la función CoVaR_DF
# <Name.Series1>: Series (Países) a gráficar.
# <grid.dims>: c(1,1) imprime gráficas individuales, cualquier ptrp valor determina
#              el número de filas y columnas para ordenar las gráficas conjuntas.
# <grid.col>: Color del grid de fondo (Por defecto nulo)
# <format.labels>: Determina el formato del label del eje x, por defecto '%Y' indica
#                  un formato anual.
# <major.ticks>: Determina bajo qué periodicidad agregar líneas de gúia en el eje
#                y, por defecto 'years'.
# <lwd>: Indica el tamaño de las líneas a gráficar en órden <c('CoVaR','VaR','DeltaCoVaR')>
# <lty>: Indica el tipo de línea en órden c('CoVaR','CoVaRUp','VaR', 'VaRUp', 'DeltaCoVaR', 'DeltaCoVaRUp')>
# <col>: Indica el color de las líneas a graficar
# <plot.type>: Tres opciones:
#             <up>: Se grafican el VaRUp, CoVaRUp y DeltaCoVaRUp
#             <down>: Se grafican el VaR, CoVaR y DeltaCoVaR
#             <both>: Se grafican los dos anteriores.
# <plot> Puede ser 'pdf' o 'x11' para elegir el método de graficación.
#
# --------------------------- Argumentos de salida ------------------------
# Sí <plot=='pdf'> La función guarda las gráficas en el directorio de trabajo. 
# <CDS>: Array con las series estimadas ordenadas por país. 
# -------------------------------------------------------------------------
Graficas=function(CoVaR_data=CoVaR_data, Name.Series1=Name.Series1, grid.dims=c(2,2), 
                  grid.col = NULL, format.labels="%Y", major.ticks='years', 
                  lwd=c(rep(0.5,length(CoVaR_data)-3)), lty=c(rep('dotted',2),rep('dashed',2)),
                  col=c(rep("red",2),rep('black',2)), plot.type=c('up','down','both')[3], plot='pdf'){
  # Creamos el Array que contiene las series
  columnas=names(CoVaR_data)[1:6]
  filas=rep(0,nrow(CoVaR_data$CoVaR)) 
  series=Name.Series1
  CDS=array(NA, dim=c(length(filas), length(columnas), length(series)), dimnames=list(c(),columnas,series))
  
  #Llenamos el Array con los objetos de salida de CoVaR_data
  for (i in Name.Series1) {
    CDS[,"CoVaR",i]        = xts(CoVaR_data$CoVaR[,i],        order.by = ichimoku::index(CoVaR_data$CoVaR[,i]))
    CDS[,"CoVaRUp",i]      = xts(CoVaR_data$CoVaRUp[,i],      order.by = ichimoku::index(CoVaR_data$CoVaRUp[,i]))
    CDS[,"VaR",i]          = xts(CoVaR_data$VaR[,i],          order.by = ichimoku::index(CoVaR_data$VaR[,i]))
    CDS[,"VaRUp",i]        = xts(CoVaR_data$VaRUp[,i],        order.by = ichimoku::index(CoVaR_data$VaRUp[,i]))
    CDS[,"DeltaCoVaR",i]   = xts(CoVaR_data$DeltaCoVaR[,i],   order.by = ichimoku::index(CoVaR_data$DeltaCoVaR[,i]))
    CDS[,"DeltaCoVaRUp",i] = xts(CoVaR_data$DeltaCoVaRUp[,i], order.by = ichimoku::index(CoVaR_data$DeltaCoVaRUp[,i]))
  }
  
  #Establecemos los parámetros gráficos iniciales.
  grid.dims=grid.dims
  max.graph=grid.dims[1]*grid.dims[2]
  ngraphs=0
  times=0
  par(mfrow=grid.dims)
  for (i in Name.Series1){
    ngraphs=ngraphs+1
    #Abrimos un nuevo dispositivo gráfico si se llena el actual dado <grid.dims>
    if (((ngraphs - 1)%% max.graph/(grid.dims[1]*grid.dims[2]))==0){
      times=times+1
      if (plot=='x11')x11()
      else {
        if (dev.cur() > 1) dev.off() #Cierra todos menos el último
        pdf(file = paste0(Resultados,'/Graficas_CoVaR_',times,'.pdf'), onefile=FALSE)
      }
      par(mfrow=grid.dims)
    }
    
    if (plot.type=='both'){
      print(plot.xts(xts(CDS[,c('CoVaR','VaR', 'CoVaRUp','VaRUp'),i], order.by = ichimoku::index(CoVaR_data$DeltaCoVaRUp[,i])), main=paste(i,'-',Serie2),lwd = lwd, lty=lty, col=c('red', 'blue', 'red', 'blue'), format.labels=format.labels, major.ticks=major.ticks, grid.col=grid.col))
      print(addLegend("topright", lwd=2, legend.names = c('CoVaR', 'VaR'), 
                lty = c('dotted','dashed'), col = c('red',   'blue')))
    }
    else if (plot.type=='down'){ 
      print(plot.xts(xts(CDS[,c('CoVaR','VaR'),i], order.by = ichimoku::index(CoVaR_data$DeltaCoVaRUp[,i])), main=paste(i,'-',Serie2),lwd = lwd, lty=
                       
                       , col=col, format.labels=format.labels, major.ticks=major.ticks, grid.col=grid.col))
      print(addLegend("topright", lwd=2, legend.names = c('CoVaR', 'VaR',  ), 
                      lty = c('dotted','dashed'), col = c('red',   'black')))
    }
    else if (plot.type=='up'){
      print(plot.xts(xts(CDS[,c('CoVaRUp','VaRUp'),i], order.by = ichimoku::index(CoVaR_data$DeltaCoVaRUp[,i])), main=paste(i,'-',Serie2),lwd = lwd, lty=lty, col=col, format.labels=format.labels, major.ticks=major.ticks, grid.col=grid.col))
      print(addLegend("topright", lwd=2, legend.names = c('CoVaR', 'VaR'), 
                      lty = c('dotted','dashed'), col = c('red',   'black')))
    }
  }
  if (plot=='pdf') {
    graphics.off()
  }
  return(CDS)
}



