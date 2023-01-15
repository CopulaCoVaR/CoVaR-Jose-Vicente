#----------------------------- values --------------------------------------#
# Indica en la consola los valores de la secuencia de Lambdas (lower.bound, upper.bound, length.out)
# Se utiliza en el <caso c> de la funcion <optimal.model> (explicado más adelante) para 
# solicitar los valores requeridos de la secuencia. 
#------------------------- Argumentos de entrada ---------------------------#
# Ningugo
#------------------------- Argumentos de salida ---------------------------#
# lista con los valores ingresados de la secuencia de Lambdas
#--------------------------------------------------------------------------#
values = function(){
  lower.bound = readline(('Insert value for the lower bound of lambda: '))
  lower.bound = as.numeric(lower.bound)
  upper.bound = readline(('Insert value for the upper bound of lambda: '))
  upper.bound = as.numeric(upper.bound)
  length.out  = readline(('Insert value for the length of lambda grid: '))
  length.out  = as.numeric(length.out)
  out = list(lower = lower.bound, upper = upper.bound, length.out = length.out)
  return(out)
}
#-------------------------------------- Optimal.model function ----------------------------------------#
# Supuestos:
# En el modelo base (Ver etapa 1), el valor de <depth>, es decir gran[1] asociado al Lambda grid está predefinido (en 10000). 
# Sí este no es lo suficientemente grande, se incrementa automáticamente en la función para alcanzar un mínimo global. 
# 
# OBJETIVO:
# Estimar los parámetros óptimos de penalizacion(alpha y Lambda) para un modelo <BigVAR> bajo la estructura "BasicEN" (Elastic Network) en dos etapas.
# -Etapa 1: Se estima el modelo con depth = <grid[1]> i.e 10000 por defecto y se extrae el Lambda y alpha óptimos para generar un intervalo 
#           más pequeño alrededor del posible  valor óptimo de  Lambda (<lower>, <upper>)  que será el input de <gran()> en la etapa 2. 
# -Etapa 2: Con el intervalo anterior se genera una secuencia lineal de Lambdas. Se usa la opción <ownlambdas=TRUE> y <gran> corresponde 
#           a la secuencia del paso 1 (seq(lower=<lower>, upper=<upper>, length.out=<seq.length>)).
#
# Para llevar a cabo este procedimiento de dos etapas se cuenta con tres opciones:
# a- Opción automática: En este caso el intervalo para evaluar Lambdas se establece de forma automática dentro de la función. 
#    * Para esto se utilizan los argumentos <auto.lambdas=TRUE> y <auto.length="integer">, "integer" es el número de valores de la secuencia de
#      Lambdas de la etapa 2.
# b- Opción con intervalo conocido: Se calcula el intervalo de Lambdas usando los argumentos <lower="numeric">, <upper="numeric"> y <seq.length="integer">.
#    * En este caso se supone que el intervalo óptimo de Lambdas de la etapa 2 es conocido, la función se ejecuta sin solicitar valores adicionales.  
# c- Opción manual(Para análisis gráfico): Si ninguno de los argumentos de las opciones a y b (<auto.lambdas>, <auto.length>, <lower>, <upper> y <seq.length>) 
#    son usados, la función estima el modelo base y muestra el gráfico asociado a los Lambdas para determinar manualmente el intervalo deseado. 
#    * La función estima el modelo base y se detiene una vez se presenta la gráfica para solicitar (en la consola) los valores asociados al 
#      límite inferior, superior y dimensión de la secuencia de Lambdas a evaluar en la segunda etapa. 
# Nota:
#  + En la etapa 1, sí <alpha> es un vector se lleva a acabo una optimización doble(Para <alpha> y <Lambda>, minimizando el MMSFE un paso adelante para la muestra
#    entre T/3 y 2T/3. En este caso se extrae el <alpha> óptimo encontrado para usarlo en la etapa 2 (En esta etapa <alpha> ya NO es un vector). 
#    Esto con la finalidad de poder llevar a cabo un análisis gráfico válido.
# ----------------------------------------------------------------------------------------------------------#
# Argumentos de entrada:
# <Data>:      Matriz de datos de formato válido para BigVAR 
# <grid>:      Dos opciones para la cuadrícula de parámetros de penalización λ.
#                 - La primera opción controla la profundidad de la cuadrícula lambda (una buena opción predeterminada es 50). 
#                 - La segunda opción controla el número de valores de cuadrícula (un buen valor predeterminado es 10).
#                   Si la cuadrícula no es lo suficientemente profunda, los resultados de los pronósticos pueden ser subóptimos, pero si es demasiado profundo, la rutina puede tardar una cantidad considerable de tiempo en ejecutarse. 
#                   El índice del parámetro de penalización óptimo es monitoreado por <cv.BigVAR()>. Si está en el borde de la cuadrícula, se recomienda volver a ejecutar la función con un parámetro de granularidad mayor. 
#                   Si establece la <ownlambdas=TRUE>, <gran> se usa para proporcionar una cuadrícula definida por el usuario.
# <structure>:  Estructura de penalización del VAR. i.e 'BasicEN'
# <alpha>:      Puede ser un entero o una secuencia, en el último caso se lleva a cabo optimización doble en la etapa 1 (Ver Nota). 
# <window.size: Tamaño de la ventana para el rolling window. 
# <linear>:     Argumento lógico, si TRUE se construye una secuencia lineal de Lambdas, FALSE construye una secuencia con escala logarítmica.
# <lower>:      Indica el límite inferior del vector de Lambdas nuevo a evaluar según el Lambda óptimo optenido previamente, la función no se interrumpe para solicitar valores. 
# <upper>:      Indica el límite superior del vector de Lambdas nuevo a evaluar según el Lambda óptimo optenido previamente, la función no se interrumpe para solicitar valores.
# <seq.length>: Indica el número de elementos qué tendrá la nueva matriz de Lambdas.
# <auto.Lambdas>: <TRUE>: se calcula el intervalo de Lambdas para la etapa 2, dónde el límite inferior es
#                   el Lambda previo al óptimo del modelo base, y el límite superior será el Lambda siguiente al óptimo en el modelo base.
#                 <FALSE>: 
#                   - Cuando lower=NULL, upper=NULL, seq.length=NULL( Se hace análisis gráfico), la función solicita en la consola los valores requeridos para 
#                     construir el vector de Lambdas para la etapa 2. 
#                   - Sí <lower>, <upper> y <seq.length> no son NULL se lleva a cabo el proceso con intervalo conocido (ver obción b, línea 28).
# <auto.length>: Indica el número de elementos de la secuencia de Lambdas para la etapa 2. Esta opción sólo es válida sí <auto.Lambdas=TRUE> (ver obción a, línea 28).
# <plot>:        Determina sí imprimir gráficas
#-------------------------------------------------------------------------------
# Agumentos de salida:
# Model  = objeto <BigVAR.results> con el modelo bajo parámetros óptimos. 
# Lambda = Lambda óptimo encontrado y usado en la etapa 2. 
# Grid   = Matriz de Lambdas óptimos estimada en la etapa 1. 
#-------------------------------------------------------------------------------
Optimal.model = function(Data=CoVaR_data$CoVaR, grid=c(10000,10), structure='BasicEN', alpha=seq(0,1,0.1),
                         window.size=150,verbose=FALSE, linear=FALSE, lower=NULL, upper=NULL, seq.length=NULL, auto.lambdas = FALSE, auto.length=NULL, plot=TRUE){
  if (auto.lambdas==TRUE & is.null(lower)==FALSE & is.null(upper)==FALSE){
    lower      = NULL
    upper      = NULL
    seq.length = NULL
    warning('auto.lambdas=TRUE; lower, upper and seq.length parameters will be disregarded.')
  }
  if(is.null(auto.length)==TRUE) auto.length=20
  
  if (is.zoo(Data)==TRUE) Data=as.matrix(fortify.zoo(Data)[,-1]) else Data=as.matrix(Data)
  cat("Stage one in progress...")
  # Inicio de la primera etapa
  if (is.null(window.size)) {
    window.size = floor(floor(nrow(Data)/3)/4)   #  El periodo de evaluación tiene tamaño T/3, por lo que se elige 0.25(T/3)
    cat('\nWindow size for BigVAR cv set to:', window.size,'\n') 
  }
  if (window.size >= floor(nrow(Data)/4)){
    stop(paste('window.size must be smaller than',floor(nrow(Data)/4)))
  }
  # Seleccion del rezago optimo
  VAR_order           = VARselect(Data)
  # Estimación del Modelo de la etapa 1.  Optimizando <Lambda> y <alpha>
  Model.base          = constructModel(Data, p=VAR_order$selection[1], struct=structure, gran=grid,
                                         cv='Rolling', window.size=window.size, linear=linear, 
                                         model.controls=list(alpha=alpha), verbose=verbose)
  Model.base.results  = cv.BigVAR(Model.base)
  index               = Model.base.results@index[1]
  counter             = 0
  #------------------------------------------------
  # Determinacion del <depth> apropiado. i.e. el primer parámetro de <grid>, en el caso de que 
  # el valor por defecto (10.000) no sea suficiente para obtener una grafica adecuada de Lambdas vs. MSFE,
  # El siguiente <while> determina el valor de <depth> aumentandolo en 500 hasta llegar a un optimo.
  # El valor optimo de Lambda se selecciona minimizando el [meanMSFE]. 
  #------------------------------------------------

  while((index==1 | index==grid[2]) & counter<30){
    warning(paste("Depth is not deep enough to reach a convex function","\n","Increasing gran[1] value by 500","\n","New value: ", (grid[1]+500)))
    grid                = c((grid[1]+500),grid[2])
   
   #--- Estimación del modelo de la etapa 1 optimizando Lambda y alpha, verificando que la grafica de los Lambdas sea convexa. 
    Model.base          = constructModel(Data, p=VAR_order$selection[1], struct=structure, gran=grid,
                                         cv='Rolling', window.size=window.size, linear=linear, 
                                         model.controls=list(alpha = alpha),verbose=verbose)
    Model.base.results  = cv.BigVAR(Model.base)
    index               = Model.base.results@index[1]
    counter             = counter+1 
    if (counter==30) {
      warning("No convex function reached for MSFE(Lambda)")
    }
  }
  if(index==1 | index==grid[2]){
    Model.base          = constructModel(Data, p=VAR_order$selection[1], struct=structure, gran=grid,
                                         cv='Rolling', window.size=window.size, linear=linear,
                                         model.controls=list(alpha = seq(0,1,0.05)),verbose=verbose)
    index               = Model.base.results@index[1]
    if(index==1 | index==grid[2]) stop("No convex function reached for MSFE(Lambda)")
    else Model.base.results  = cv.BigVAR(Model.base)
  }
  
  #--- Estimacion de la Etapa 2 ---#
  if (is.null(lower) & is.null(upper) & is.null(seq.length)){
    # Procedimiento bajo opción (a)
    if (auto.lambdas==TRUE) {
      if (length(alpha)!=1) {
        # Secuecia de Lambdas usada en la etapa 2 opción a, cuando alpha es vector 
        # Nota: En este caso <LambdaGrid> es una matriz, las filas es la seq ini. de <lambdas> y las cols son <lambdas> asociados
        #       al  vector de <alphas> considerado.
        Lambda.range.seq = try(seq(from=Model.base.results@LambdaGrid[ Model.base.results@index[1]+1, Model.base.results@index[2] ],
                                   to=Model.base.results@LambdaGrid[   Model.base.results@index[1]-1, Model.base.results@index[2] ], 
                                   length.out=auto.length))
        if (class(Lambda.range.seq)=='try-error') Lambda.range.seq = Model.base.results@LambdaGrid[,Model.base.results@index[2]]
        # Modelo de la etapa 1 - opción a, con alpha óptimo(escalar) obtenido en la inicializacicón del modelo en la línea 93
        Model.base          = constructModel(Data, p=VAR_order$selection[1], struct=structure, gran=grid,
                                             cv='Rolling', window.size=window.size, linear=linear, model.controls=list(alpha = Model.base.results@alpha))
        Model.base.results  = cv.BigVAR(Model.base)
        if (plot==TRUE) {
          x11()
          plot(Model.base.results)
          title(paste('Optimal Lambda for initial model: ', round(Model.base.results@OptimalLambda, 5) ))
        }
      }else{
        # Secuencia de Lambdas usada en la etapa 2 - opción a cuando alpha es un escalar. 
        Lambda.range.seq = try(seq(from=Model.base.results@LambdaGrid[((Model.base.results@index[1])+1)], 
                                   to=Model.base.results@LambdaGrid[((Model.base.results@index[1])-1)], 
                                   length.out=auto.length))
        if (class(Lambda.range.seq)=='try-error') Lambda.range.seq = Model.base.results@LambdaGrid[,Model.base.results@index[2]]
      }
    }else{
      # Procedimiento bajo opción (c)
      if (length(alpha)!=1) {
        # Modelo de la etapa 2 - opción (c), con alpha vector obtenido en la inicializacicón del modelo en la línea 93
        Model.base          = constructModel(Data, p=VAR_order$selection[1], struct=structure, gran=grid,
                                             cv='Rolling', window.size=window.size, linear=linear, 
                                             model.controls=list(alpha = Model.base.results@alpha))
        Model.base.results  = cv.BigVAR(Model.base)
        x11()
        plot(Model.base.results)
          title(paste('Optimal Lambda for initial model: ', round(Model.base.results@OptimalLambda, 5) ))
          cat("Lamda grid\n",Model.base.results@LambdaGrid,'\n') #Verificar que salga nombre impreso
          cat("Indice del Lambda óptimo:", Model.base.results@index)
        sequence         = values()
        Lambda.range.seq = seq(sequence$lower,sequence$upper, length.out = sequence$length.out)
        
      }else{ # Modelo de la etapa 2 - opción (c), con alpha escalar
        x11()
        plot(Model.base.results)
          title(paste('Optimal Lambda for initial model: ', round(Model.base.results@OptimalLambda, 5) ))
          cat("Lamda grid\n",Model.base.results@LambdaGrid,'\n') #Verificar que salga nombre impreso
          cat("Indice del Lambda óptimo:", Model.base.results@index)
        sequence         = values()
        Lambda.range.seq = seq(sequence$lower,sequence$upper, length.out = sequence$length.out)
       }
    }
  # Procedimiento bajo opción (b)
  }else{
    # Secuecia de Lambdas usada en la etapa 2 bajo la opción b. 
    Lambda.range.seq = seq(lower, upper, length.out=seq.length)
  }
  cat("Stage two in progress...")
  # Estimacion de la etapa 2 (Para las opciones a, b ó c)
  Model.optimal         = constructModel(Data, p=VAR_order$selection[1], struct=structure, gran=Lambda.range.seq,
                                         cv='Rolling', ownlambdas = TRUE, window.size = window.size, 
                                         linear = TRUE, verbose=verbose, model.controls = list(alpha = Model.base.results@alpha))
  Model.optimal.results = cv.BigVAR(Model.optimal)
  
  if (plot==TRUE) {
    x11()
    plot(Model.optimal.results)
     title(paste('Optimal Lambda for optimal model: ', round(Model.optimal.results@OptimalLambda,5) ))
  }  
  out = list(Model = Model.optimal.results, Lambda = Model.optimal.results@OptimalLambda, Grid = Model.optimal.results@LambdaGrid)
  
  return(out)
}

if (0) {
  Optimal.auto   = Optimal.model(Data=Y, window.size=25, p=4, grid = c(10,10), auto.lambdas = TRUE)                                   # Ejemplo de <optimal_model()> bajo caso a
  Optimal.intcon = Optimal.model(Data=Y, window.size=25, p=4, grid = c(10000,10), lower=0.005 , upper=0.015 , seq.length=20)  # Ejemplo de <optimal_model()> bajo caso b
  Optimal.manual = Optimal.model(Data=Y, window.size=25, p=4, grid = c(10000,10))                                                        # Ejemplo de <optimal_model()> bajo caso c
}
