# Prueba de modificación
#-------------------------------- Funcion Connectedness  -------------------------------
# Permite obtener la tabla y grafico de conectividad propuesta por Diebold & Yilmaz (2014)
# Además realiza el gráfico de red asociado a la tabla, por lo que crea una 
# tabla adicional con el formato adecuado para ello. 
# ---------------- Parámetros de entrada de la función  ------------------------------
# <GFEVD>          : Descomposición de varianza extraída de la función GFEVD
# <n.ahead>        : Paso adelante de la descomposición de varianza que se quiere analizar (default = 10)
# <type>           : Dado que la función <GFEVD> extrae diferentes tipos de descomposición de varianza,
#                    se puede elegir cuál usar ("FEVD", "GFEVD", "GFEVD.A", "GFEVD.A")
# <plot.network>   : Variable lógica, si <TRUE> se realiza la gráfica de conectividad, si <FALSE>
#                    sólo se obtienen las matrices de conectividad. 
# <layout>         :   FALTA!!!  
# <Title>          : String que determina el título de la gráfica. 
# <edge.width.scale>: Determina el escalado del grosor de las conexiones en función del peso de la misma. 
# <node.size>      : Determina el escalado de los nodos, por defecto 2. 
# <arrow.soze>     :   FALTA!!!
# <node.label.size>: Determina el escalado del label de los nodos.
# <Q1>             : Cuantil superior a graficar, formato de porcentaje entre comillas ej.   ''95%''
# <Q2>             : Cuantil intermedio a graficar, formato de porcentaje entre comillas ej. ''85%''
# <Q3>             : Cuantil inferior a graficar, formato de porcentaje entre comillas ej.   ''80%''
# <Q1.scale>       : Por defecto 1.5, permite escalar las conexiones para el cuartil superior.
# <Q2.scale>       : Por defecto 1.5, permite escalar las conexiones para el cuartil intermedio
# <Q3.scale>       : Por defecto 1.5, permite escalar las conexiones para el cuartil inferior.
#---------------------- Objetos de salida de la función ----------------------------
# <table>:   FALTA !!!
# <gephi>:   FALTA !!!
#-----------------------------------------------------------------------------------
Connectedness = function(GFEVD, n.ahead=10, type='GFEVD', plot.network=TRUE, layout = NULL, 
                         Title='Insert title', edge.width.scale=30, node.size=2, arrow.size = 1,
                         node.label.size=1, Q1 = '95%', Q2 = '90%', Q3 = '85%', Q1.scale = 1.5, 
                         Q2.scale = 1, Q3.scale = 1 ){
  
  #Verificacion de que los cuantiles enten on orden descendente
  if (Q1 < Q2 | Q1 < Q3 | Q2 < Q3) {
    stop('The Q1, Q2 and Q3 quantiles must satisfy: Q1 > Q2 > Q3')
  }
  
  #Definimos el objeto correspondiente segun la opción de <type>
  n.vars = ncol(GFEVD$GFEVD)
  if (type == 'GFEVD')   FE = GFEVD$GFEVD
  if (type == 'GFEVD.A') FE = GFEVD$GFEVD.A
  if (type == 'FEVD')    FE = FEVD$FEVD
  if (type == 'FEVD.A')  FE = FEVD.A$FEVD.A
  
  #Definicion y calculo de la tabla de conectividad. Table 3 de Diebold y Yilmaz (2014)
  table  = matrix(nrow = n.vars, ncol = n.vars)
  From   = rep(0, n.vars)
  To     = rep(0, n.vars+1)
  
  for (i in 1:n.vars) {
    table[i,] = FE[n.ahead, ,i]
  }
  table            = t(table)*100
  diag(table)      = 0
  
  for (i in 1:n.vars) {
    From[i] = sum(table[i,])
    To[i]   = sum(table[,i])   
  } 
  To[n.vars+1] = sum(table)/n.vars
  table = cbind(table, From)
  table = rbind(table, To)
  colnames(table)   = c(colnames(FE), 'From')
  rownames(table)   = c(colnames(FE), 'To')
  
  #Definicion de la tabla de conectividad usada para graficar. 
  #Se quita la ultima fila y columna ('From' y 'to') y se ignoran los valores menores
  table_gephi       = (table[1:n.vars,1:n.vars])/100 #Se retorna a los valores originales
  diag(table_gephi) = 0 # Creamos la tabla para el formato Gephi
  # <table_gephi_original>: Objeto gephi para graficacion de conectividad en R (igraph)
  # <table_gephi>:  Objeto gephi para graficacion de conectividad fuera de R (Gephi)
  table_gephi_original = (table[1:n.vars,1:n.vars])/100 #Se retorna a los valores originales
  diag(table_gephi_original) = 0 # Creamos la tabla para el formato Gephi
  quanti = quantile(table_gephi, probs=seq(0,1,0.1))
  for (i in 1:(nrow(table_gephi)*ncol(table_gephi))){
    if (table_gephi[i] < quanti['70%'] ) table_gephi[i] = 0
  }
  
  # Grafica de conectividad 
  if (plot.network == TRUE) {
    #Creacion del objeto <melt>, que es adecuado para el paquete <igraph>
    melt1                      <- melt(table_gephi)
    melt                       <- melt1 %>% filter( value != 0 ) #Solo para valores dif. de cero
    colnames(melt)[3]          <- "weight"
    melt_original              <- melt(table_gephi_original)
    melt_original              <- melt_original %>% filter( value != 0 ) #Solo para valores dif. de cero
    colnames(melt_original)[3] <- "weight"
    ### calculate each percentile of the net pairwise connectedness values
    ### and choose only the top 10%
    quant    <- quantile(melt_original$weight,prob=seq(0, 1, by=0.01))
    new.net1 <- melt[melt$weight >= quant['80%'],]      #Solo se dejan cuantiles > 80%
    new.net  <- new.net1[new.net1[,1] != new.net1[,2],] #Se quitan conexiones consigo mismo
    ### create igraph graph
    network  <- graph.data.frame(new.net,direct=T)
    #network <- set_edge_attr(network, "weight", value= new.net$weight)
    E(network)$weight <- as.numeric(new.net$weight) #Ponderacion de los Edges 
    ### set graph nodes
    #V(network)
    ### set edge colors
    #Deficion de los cuantiles con base en lo parametros de entrada
    Q1 = Q1
    q1.u = 100 - as.numeric(strsplit(Q1, '%')[1])
    Q2 = Q2
    q2.u = 100 - as.numeric(strsplit(Q2, '%')[1])
    Q3 = Q3
    q3.u = 100 - as.numeric(strsplit(Q3, '%')[1])
    #Definicion de los colores de los Edges segun los tres cuantiles indicados 
    E(network)$color <- ifelse(E(network)$weight >= quant[Q1],"black",
                               ifelse(E(network)$weight < quant[Q1] & E(network)$weight >= quant[Q2],"red",
                                      ifelse(E(network)$weight < quant[Q2] & E(network)$weight >= quant[Q3],"orange", 'orange')))
    ### scaling edges by their quantile
    for (i in 1:length(E(network)$weight)) {

      if (E(network)$weight[i] >= quant[Q1]) { 
        E(network)$weight[i] = E(network)$weight[i]*Q1.scale
      }
      if (E(network)$weight[i] < quant[Q1] & E(network)$weight[i] >= quant[Q2]) { 
        E(network)$weight[i] = E(network)$weight[i]*Q2.scale
        
      }
      if (E(network)$weight[i] <= quant[Q2] & E(network)$weight[i] > quant[Q3]) { 
        E(network)$weight[i] = E(network)$weight[i]*Q3.scale
      }
      
    }
    ### set node size
    V(network)$size <- node.size*degree(network) 
    # Grafica de conectivad
    if (is.null(layout) == TRUE) { #Estilo por defecto (teleranha)
      x11()
      plot(network, edge.width = ((E(network)$weight)*edge.width.scale),
           main = Title, vertex.label.cex = node.label.size, edge.arrow.size = arrow.size,
           xlab = paste0("black: Upper " ,  q1.u ,"th percentile \n red: Upper ", q2.u ,"th percentile \n orange: Upper ", q3.u , "th percentile \n node size: number of edges connected to the node"))
      
    }
    else if (layout == 'circle'){ #Estilo poliedro
      x11()
      plot(network,layout=layout.circle(network),
           edge.width = ((E(network)$weight)*edge.width.scale),
           main = Title,vertex.label.cex = node.label.size,edge.arrow.size = arrow.size,
           xlab = paste0("black: Upper " ,  q1.u ,"th percentile \n red: Upper ", q2.u ,"th percentile \n orange: Upper ", q3.u , "th percentile \n node size: number of edges connected to the node"))
      
    }
    
  }
  list(table=table, gephi=table_gephi)
}

#for (i in 1:nrow(Network.mat)) {
 # Network.mat[i,(ncol(Network.mat)+1)] = sum(Network.mat[1,])
#}

