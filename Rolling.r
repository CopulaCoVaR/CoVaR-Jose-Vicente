library(readxl)
library(expm)
library(BigVAR)
library(ggplot2)
#-------------------------------------------------------------------------------
#No usada, ayuda a determinar el mejor argumento para gran[1] en constructModel en BigVAR
grid_test = function(Data=CoVaR_Data_BV, values=c(50, 100, 500, 1000, 2500, 5000, 10000, 20000),
                     p=2, verbose=FALSE, structure='BasicEN', alpha=seq(0,1,0.1), window.size = 150, linear = FALSE){
  
  Optimal_lambdas = matrix(0, nrow = 4, ncol = length(values), dimnames = list(c('Lambda Value', 'OOSMSFE','seoosmsfe', 'alpha'), as.character(values)))
  
  for (i in values) {
    
    Model             = constructModel(Y = as.matrix(Data), p = p , struct=structure, gran=c(i,10),
                                       cv = 'Rolling', linear = linear, window.size = window.size, 
                                       verbose=verbose, model.controls = list(alpha = alpha))
    
    Model_Results     = cv.BigVAR(Model)
    
    Optimal_lambdas[1,as.character(i)] = Model_Results@OptimalLambda
    Optimal_lambdas[2,as.character(i)] = mean(Model_Results@OOSMSFE)
    Optimal_lambdas[3,as.character(i)] = Model_Results@seoosmsfe
    Optimal_lambdas[4,as.character(i)] = Model_Results@alpha

  }
  # Extraemos los valores óptimos de <gran[1]> y alpha que minimizan el MSFE
  
  bestgrid        = as.numeric(names(which.min(Optimal_lambdas[2,])))
  bestalpha       = Optimal_lambdas[4,as.character(bestgrid)]
    
  # Estimamos el modelo con los parámetros óptimos encontrados anteriormente:
  
  Optimal_Model   = constructModel(Y = as.matrix(Data), p = p , struct=structure, gran=c(as.numeric(bestgrid),10),
                                     cv = 'Rolling', linear = linear, window.size = window.size, 
                                     verbose=verbose, model.controls = list(alpha = bestalpha))
  
  Optimal_Results = cv.BigVAR(Optimal_Model)
  
  x11()
  plot(Optimal_Results)
  
  list            = list(Optimal_lambdas=Optimal_lambdas, grid=bestgrid, alpha = bestalpha, Optimal_model = Optimal_Results)
  return(list)
}

#Auto_Demirer = Auto_BigVAR(Data = Demirer_21, p = 10, structure = 'BasicEN', gran = c(50,10), alpha = seq(0,1,0.01), window = ceiling(0.25*nrow(Data)))


#---------------------------------------- Función Rolling_GFEVD --------------------------------------------
# Supuestos:
# - Se supone una muestra de tamaño T para el rolling del BigVAR en cada ventana.
# - Se supone que ya se estimó un modelo para toda la muestra del cuál se toma el alpha.
# Objetivo:       
# a. Obtener la serie de conectividad total a lo largo del tiempo. (Conectividad dinámica total)
# b. Obtener las series de conectividad dinámica hacia otros de cada una de las variables del VAR. 
# c. Obtener las series de conectividad dinámica desde otros de cada una de las variables del VAR. 
# d. Obtener las series de conectividad dinámica neta (c-b) de cada una de las viarbles del VAR.
#
# Procedimiento: En esta función, tomando como input las series de CoVaR estimadas y determinando un tamaño de 
#                ventana <window.size> se realiza el siguiente procedimiento rolling para cada ventana:
#                (Dónde la 1ra ventana es 1:<window.size>, 2da ventana es 2:(<window.size>+1), y así sucesivamente) 
#                
#                 1- Se estima el BigVAR con los datos de las series de CoVaR tomando <window.size> observaciones.
#                 2- Calculamos los impulso respuesta (<GIRF_BigVAR>) y la descomposición de varianza (<GFEVD_BigVAR>) para el horizonte desde 1 hasta <n.ahead>.
#                 3- Calculamos la matriz de conectividad a partir de los valores del (<GFEVD_BigVAR)> para el horizonte <n.ahead>.
#                    Esta matriz es de dimensión mxn (Dónde m: , n: ). Por ejemplo, el componente 2,1 de esta matriz es...
#                 4- Se extrae el valor de la conectividad total de la matriz calculada en el paso anterior.
#                    La conectividad total se define como...  
#                Haciendo esto en cada ventana obtenemos una serie de tiempo que nos brinda información acerca
#                de las dinámicas de conectividad total a lo largo del tiempo.
#
# ----------------------------------- Parámetros de entrada de la función  ----------------------------------
# <Data> :         Matriz de datos de formato válido para BigVAR
# <structure> :    Estructura de penalización del VAR
# <window.size>:   Tamaño de ventana de estimación para el rolling window
# <BigVAR.window>: Tamaño de ventana para el rolling del BigVAR.
#                   Notas: 
#                         * Dado que el periodo de evaluación del BigVAR es entre T/3 y 2T/3, es decir de una longitud de T/3
#                           este parámetro debe ser inferior a T/3. En este caso se exige que sea al menos
#                           inferior a 1/4 de T/3 para poder tener suficientes números de ventanas rolling.
#                         * El programa se detiene cuando este parámetro es mayor a (T/3)/4.
# <n.ahead> :      Número de  pasos adelante para las impulso respuesta y la descomposición de varianza. 
#                  Determina también el horizonte de tiempo en el cual se calcula la matriz de conectividad. 
# <alpha> :        Valor de alpha para la penalización, se supone constante para todas las ventanas. 
# <pdf> :          Variable lógica para guardar archivos pdf en la carpeta de trabajo.
# <x11> :          Variable lógica para abrir una ventana con gráficas. ELIMIAR ARGUMENTO DE X11() Y SÍ PDF=FALSE QUE SOLO IMPRIMA GRÁFICA
# <show.progress>: Variable lógica. <TRUE> muestra el progreso de la función (Adaptado para Mac y Windows)
# <plot>: Variable lógica. <TRUE> grafica las series estimadas. <FALSE> no grafica y sólo se guarda el objeto con las series correspondientes de cada país. 
# <Num.alphas.rolling>: Entero que determina cuantas veces reestimar alpha. 
#
#---------------------------------- Objetos de salida de la función ---------------------------------------
# Salida :  La función retorna una lista con 4 objetos, y los grafica como serie de tiempo. (graficación opcional)
#
#               1- <$Total.Dynamic> vector columna que contiene el total de la matriz de conectividad
#                  estimado en cada ventana, permite ver como cambia la conectividad total del sistema 
#                  a lo largo del tiempo.
#               2- <$From.Degree.Dynamic> Contiene una matriz con número de columnas igual al número de 
#                  de variables del sistema, contiene en cada fila la conectividad total "desde otros"
#                  hacia la variable de interés, permite ver los cambios a lo largo del tiempo de la conectividad
#                  'desde otros', lo que permite ver en qué momento del tiempo, esta variable fué más vulnerable a
#                  ser receptora de riesgo (Si usamos como series medidas de riesgo).
#                  Se define como... (Rango: ...)
#               3- <$To.Degree.Dynamic> Contiene una matriz con número de columnas igual al número de 
#                  de variables del sistema, y en cada fila la conectividad total "hacia otros"
#                  desde la variable de interés, permite ver los cambios a lo largo del tiempo de la conectividad
#                  'hacia otros', lo que permite analizar en qué momento del tiempo esta variable fué transmisora de riesgo.
#                  Se define como... (Rango: ...)
#               4- <$Net.Dynamic> Contiene una matriz con número de columnas igual al número de 
#                  de variables del sistema, y contiene en cada fila la conectividad total neta
#                  de la variable de interés, es decir, <To.Degree.Dynamic> - <From.Degree.Dynamic>, lo que permite ver en qué momento del
#                  tiempo esta variable fué transmisora o receptora neta de riesgo.
#--------------------------------------------------------------------------------------------------------------
Rolling_GFEVD = function(Data, structure='BasicEN', window.size, BigVAR.window=NULL,
                         n.ahead=10, alpha=0.3, pdf=FALSE, x11=FALSE,
                         show.progress=TRUE, plot=FALSE, Num.alphas.rolling=5){
  
  New_length         = nrow(Data)-window.size+1 # Número de observaciones totales del ejercicio rolling.
  
  # En qué puntos de la estimación se va a reestimar alpha
  alpha.seq          = seq(1, New_lenght, length.out=(Num.alphas.rolling+1))
  alpha.seq          = ceiling(alpha.seq[-(Num.alphas.rolling+1)])
  
  # Transformación de datos al formato adecuado para BigVAR en caso de ser necesario (Formato apropiado matrix())
  if (is.zoo(Data)==TRUE) {
    alphas.mat     = xts(matrix(0, nrow = (nrow(Data)-window.size+1), ncol = 1),          order.by = ichimoku::index(Data[(window.size):nrow(Data),])) # Matriz que contendrá la conectividad total estimada en cada ventana
    Total          = xts(matrix(0, nrow = (nrow(Data)-window.size+1), ncol = 1),          order.by = ichimoku::index(Data[(window.size):nrow(Data),])) # Matriz que contendrá la conectividad total estimada en cada ventana
    From_Degree    = xts(matrix(0, nrow = (nrow(Data)-window.size+1), ncol = ncol(Data)), order.by = ichimoku::index(Data[(window.size):nrow(Data),])) # Matriz que contendrá la conectividad 'desde otros' estimada en cada ventana
    To_Degree      = xts(matrix(0, nrow = (nrow(Data)-window.size+1), ncol = ncol(Data)), order.by = ichimoku::index(Data[(window.size):nrow(Data),])) # Matriz que contendrá la conectividad 'hacia otros' estimada en cada ventana
    Net_Dynamic    = xts(matrix(0, nrow = (nrow(Data)-window.size+1), ncol = ncol(Data)), order.by = ichimoku::index(Data[(window.size):nrow(Data),])) # Matriz que contendrá la conectividad total neta estimada en cada ventana
    Data.mat       = as.matrix(fortify.zoo(Data)[,-1])
  }else{
    alphas.mat     = matrix(0, nrow = (nrow(Data)-window.size+1), ncol = 1)          # Matriz que contendrá la conectividad total estimada en cada ventana
    Total          = matrix(0, nrow = (nrow(Data)-window.size+1), ncol = 1)          # Matriz que contendrá la conectividad total estimada en cada ventana
    From_Degree    = matrix(0, nrow = (nrow(Data)-window.size+1), ncol = ncol(Data)) # Matriz que contendrá la conectividad 'desde otros' estimada en cada ventana
    To_Degree      = matrix(0, nrow = (nrow(Data)-window.size+1), ncol = ncol(Data)) # Matriz que contendrá la conectividad 'hacia otros' estimada en cada ventana
    Net_Dynamic    = matrix(0, nrow = (nrow(Data)-window.size+1), ncol = ncol(Data)) # Matriz que contendrá la conectividad total neta estimada en cada ventana
    Data.mat       = as.matrix(Data)
  }
  colnames(Total)       = 'Total'                                                       
  colnames(From_Degree) = colnames(Data)                                         
  colnames(To_Degree)   = colnames(Data)                                         
  colnames(Net_Dynamic) = colnames(Data)                                         
  
  if (is.null(BigVAR.window)) {
    BigVAR.window = floor((window.size/3)/4)  #  El periodo de evaluación tiene tamaño T/3, por lo que se elige 1/4 de este para el rolling
    cat('\nWindow size for BigVAR cv (<BigVAR.window>) set to:', BigVAR.window,'\n') 
  }
  if (BigVAR.window > floor(window.size/4)){
    stop(paste('<BigVAR.window> must be smaller than',floor(window.size/4)))
  }
  
  # Especificación de las barras de progreso
  if (Sys.info()["sysname"] == 'Darwin'){ # Barra de progreso compatible con cualqueir sistema operativo (Escrito originalmente para MacOs)
    
    pb <- txtProgressBar(  min   = 0,
                           max   = New_length,
                           style = 3,
                           width = New_length, # Necesario para evitar múltiples prints
                           char  = "=") 
  }
  if (Sys.info()["sysname"] == 'Windows'){ # Del paquete 'utils', instalado por defecto en Windows. 
    
    pb <- winProgressBar(title = "Barra de progreso de Windows", # Título de la ventana
                         label = "Porcentaje completado",        # Texto de la ventana
                         min = 0,                                # Valor mínimo de la barra
                         max = New_length,                       # Valor máximo de la barra
                         initial = 0,                            # Valor inicial de la barra
                         width = 300L)                           # Ancho de la ventana 
  }
  
  inicio = New_length #Parámetro de la ventana de carga de Mac
  fin    = New_length #Parámetro de la ventana de carga de Mac
  
  #----- Loop del rolling window -----#
  for (i in 1:New_length){
    if(show.progress==TRUE | Sys.info()["sysname"]=='Darwin') inicio[i]=Sys.time()
    # Datos de cada ventana rollling.  
    DF              = Data.mat[i:(window.size+i-1),] # Matriz de datos de número de columnas igual al número de variables y número de filas igual al tamaño de ventana, esta matriz cambia en cada ventana.
    #Etapas de estimación
    #Nota: Se asume que el <alpha> de toda la muestra es el mismo para cada ventana rolling.
    if (i %in% alpha.seq) { #Obj en la cual se re-estima <alpha>
      BigVaR.spec   = Optimal.model(Data=DF, alpha=seq(0,1,0.05), auto.lambdas=T, auto.length=20, plot=plot, window.size=BigVAR.window)
      BigVAR        = BigVaR.spec$Model 
      alpha.actual  = BigVAR@alpha
    }else{ #Se usa el <alpha> estimado en las lineas anteriores
      BigVaR.spec   = Optimal.model(Data=DF, alpha=alpha.actual, auto.lambdas=T, auto.length=20, plot=plot, window.size=BigVAR.window)
      BigVAR        = BigVaR.spec$Model 
    }
    # Calculo de GIRF
    GIRF     = GIRF_BigVAR(pdf=FALSE, x11=FALSE, results=BigVAR,          # Con la salida de la <Auto_BigVAR> se calculan las impulso respuesta <n.ahead> pasos adelante en cada ventana.
                           data=DF, n.ahead=n.ahead, ortho=FALSE, 
                           vars.to.shock=colnames(Data), vars.to.resp=colnames(Data),
                           grid.dims=c(2,2))
    # Calculo de GFEVD
    GFEVD_BV = GFEVD(Girf=GIRF, N.ahead=n.ahead, 
                     vars.to.shock=colnames(Data), vars.to.resp=colnames(Data))   # Con la salida de <GIRF_BigVAR> se calcula la descomposición de varianza generalizada del VAR estimado.
    # Calculo de la tabla de conectividad
    Dynamic  = Connectedness(GFEVD_BV, n.ahead=n.ahead, type='GFEVD',       # Usando la función <Connectedness> Obtenemos las respectivas matrices de conectividad en cada ventana.
                             plot.network=FALSE)                            # <type> puede ser 'GFEVD' o 'GFEVD.A'
    Matrix   = Dynamic$gephi # Matrix que se genera cada ventana            # Del objeto de salida de <Connectedness> extraemos <$gephi>, que corresponde a la matriz de conectividad apropiada para graficación donde la diagonal se convierte en 0's para no tener en cuenta conectividad consigo mismo.
    
    #Estadísticas de conectividad en cada ventana rolling. 
    alphas.mat[i,]   = BigVAR@alpha
    Total[i,]        = sum(Matrix)/ncol(DF)                  #De la matriz calculada en cada ventana se extrae la conectivdad total y se llena la matriz correspondiente
    From_Degree[i,]  = apply(Matrix, MARGIN=1, FUN=sum)      # From-Degree. De la matriz calculada en cada ventana se extrae la conectivdad 'Desde otros'  y se llena la matriz correspondiente
    To_Degree[i,]    = apply(Matrix, MARGIN=2, FUN=sum)      # To-Degree. De la matriz calculada en cada ventana se extrae la conectivdad 'Hacia otros' y se llena la matriz correspondiente
    Net_Dynamic[i,]  = To_Degree[i,] - From_Degree[i,]       # Se restan las dos matrices obtenidas anteriormente para obtener el efecto neto de conectividad dinámica
    # Cuerpo de la ventana de progreso para cada OS.
    if (show.progress == TRUE) {
      if (Sys.info()["sysname"]=='Darwin'){
        
        fin[i] <- Sys.time()
        setTxtProgressBar(pb, i)
        tiempo <- round(seconds_to_period(sum(fin - inicio)), 0)
        
        # Tiempo restante estimado basado en el tiempo
        # medio que tardaron en ejecutarse las iteraciones previas
        est <- New_length * (mean(fin[fin != 0] - inicio[inicio != 0])) - tiempo
        restante <- round(seconds_to_period(est), 0)
        
        cat(paste(" // Tiempo ejecución:", tiempo,
                  " // Tiempo restante estimado:", restante), "")
      }
      if (Sys.info()["sysname"]=='Windows'){
        pctg <- paste(round(i/New_length *100, 0),"% completado")
        setWinProgressBar(pb, i, label=pctg) # Al pasar pctg al argumento label se sobreescribirá
        # el texto del argumento label de la función winProgressBar
      }
    }
  }
  
  # Gráficas de conectividad dinámica total. 
  if (x11==TRUE) x11()
  if (pdf==TRUE) pdf(file=paste0(Resultados,'/Conectividad_total.pdf'), onefile=FALSE)
  print(plot.xts(Total, main='Conectividad Total dinámica', xlab='Time',ylab='%',type='l', col='steelblue', lwd=2)) # Se grafica la conectividad total
  
  # Graficas de conectividad dinámicas (From, To y Net) de cada variable. 
  # Se grafican, en un ordenamiento de 3*1, tanto la conectividad desde otros, hacia otros y neta. El usuario puede decidir entre guardarlas en PDF o imprimirlas en el dispositivo gráfico.
  if (pdf==TRUE & dev.cur() > 1) dev.off() #Cierra todos menos el último
  for (kk in colnames(Data)){
      if (x11==TRUE) x11()
      if (pdf==TRUE) pdf(file = paste0(Resultados,'/Conectividad_total_de_', kk,'.pdf'), onefile=FALSE)
      par(mfrow=c(3,1))
      print(plot.xts(From_Degree[, kk], main=paste('Conectividad desde otros para',kk), xlab='Time',ylab='Porcentaje',type='l', col='steelblue', lwd=2))
      print(plot.xts(To_Degree[, kk],main=paste('Conectividad hacia', kk,  'deste otros'), xlab='Time',ylab='Porcentaje',type='l', col='steelblue', lwd=2))
      print(plot.xts(Net_Dynamic[, kk],main=paste('Conectividad neta de',kk), xlab='Time',ylab='Porcentaje',type='l', col='steelblue', lwd=2))
      print(abline(h = 0, col='tomato', lty = 'dotted', lwd = 2))
      print(mtext(paste('Conectividad dinámica de', kk), side=3, cex =1.3, line = 0.5, outer=TRUE))
  }
  if (pdf==TRUE) graphics.off()
  
  # Argumentos de salida 
  list = list(Total.Dynamic = Total, From.Degree.Dynamic = From_Degree, To.Degree.Dynamic = To_Degree, Net.Dynamic = Net_Dynamic) 
  return(list)
}