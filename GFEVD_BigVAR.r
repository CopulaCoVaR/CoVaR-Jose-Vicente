#-------------------- Objetos de salida de la función GIRF -------------------------
# IRF:    Arreglo de tres dimensiones con los IRF (<n.ahead> x <n.var> x <n.var> )
#         [ (horizonte de los IRF) x (variables rta) x (variables chocadas)  ] 
# GIRF:   Arreglo de tres dimensiones con los IRF Generalizados
# IRF.A:  Arreglo de tres dimensiones con los IRF acumulados
# GIRF.A: Arreglo de tres dimensiones con los GIRF acumulados
#-----------------------------------------------------------------------------------
library(ggplot2)


#-------------------------------- Funcion GFEVD  --------------------------------------
# Calcula los FEVD y GFEVD, ambos en niveles y acumulados basados en los objetos de 
# salida de la funcion <GIRF>, explicados en las lineas 1 a 7. 
# Los FEVD y GFEVD se calculan de acuerdo a Lanne & Nymberg (2016)
# ---------------- Parámetros de entrada de la función  ------------------------------
# <Girf>:      Nombre del objeto que contiene los parametros de salida de la
#              funcion <GIRF>
# <n.ahead>:   Número de pasos adelante de la FEVD
# <vars.to.shock>:  Labels de las variables que se desean chocar
# <vars.to.resp> :  Labels de las variables sobre las cuales se desean ver las
#                   descomposiciones
# <plot.out>:       Indica que tipo de ka FEVD que se grafican. Valores:
#                 <'FEVD'>     :  Para la FEVD en niveles
#                 <'GFEVD'>    :  Para la GFEVD en niveles
#                 <'cum.FEVD'> :  Para la FEVD acumulada
#                 <'cum.GFEVD'>:  Para la GFEVD acumulada
#---------------------- Objetos de salida de la función ----------------------------
# FEVD:    Arreglo de tres dimensiones con los FEVD (<n.ahead> x <n.var> x <n.var> )
#         [ (horizonte de los FEVD) x (variables rta) x (variables chocadas)  ] 
# GFEVD:   Arreglo de tres dimensiones con los FEVD Generalizados
# FEVD.A:  Arreglo de tres dimensiones con los FEVD acumulados
# GFEVD.A: Arreglo de tres dimensiones con los GFEVD acumulados
#-----------------------------------------------------------------------------------
GFEVD <- function(Girf, N.ahead=NULL, vars.to.shock=NULL,vars.to.resp=NULL, plot.out=NULL ){
  N.Ahead = dim(Girf$IRF)[1]
  n.ahead = min(N.ahead, N.Ahead)
  n.var   = dim(Girf$IRF)[2]
  y.names = names(Girf$IRF[1,1,])
  #----------------- FEVD y  FEVD.A (A: acumulados) -------------------#
  FEVD     = Girf$IRF - Girf$IRF #Ceros
  FEVD.A   = FEVD
  # <FEVD.den>: denominator eqn (4) in Lanne & Nymberg (2016)
  FEVD.den   = array(data=0, dim = c(n.ahead, n.var), dimnames = list(NULL,y.names))   
  FEVD.A.den = FEVD.den   
  #<temp>: Sum_(j=1,...,k){ IRF^2(h,i,j) }
  temp     = matrix(0,n.ahead, n.var) 
  temp.A   = temp
  for (i in 1:n.var)
    for (h in 1:n.ahead){
      temp[h,i]   =  sum(Girf$IRF[h,i,]^2)
      temp.A[h,i] =  sum(Girf$IRF.A[h,i,]^2)
    }  
  for (i in 1:n.var){
    FEVD.den[,i]   = cumsum(temp[,i])    
    FEVD.A.den[,i] = cumsum(temp.A[,i])    
  }  
  for (i in 1:n.var)
    for (j in 1:n.var){# Eq 4 de Lanne & Nymberg (2016)
      FEVD[,i,j]   = cumsum(Girf$IRF[,i,j]^2) / FEVD.den[,i] 
      FEVD.A[,i,j] = cumsum(Girf$IRF.A[,i,j]^2) / FEVD.A.den[,i] 
    }  
  
  #---------------- GFEVD y GFEVD.A (A: acumulados)  ------------------#
  GFEVD     = Girf$IRF - Girf$IRF #Ceros
  GFEVD.A   = GFEVD
  # <GFEVD.den>: denominator eqn (9) in Lanne & Nymberg (2016)
  GFEVD.den   = array(data=0, dim = c(n.ahead, n.var), dimnames = list(NULL,y.names))   
  GFEVD.A.den = GFEVD.den 
  #<temp>: Sum_(j=1...k){ GIRF^2(h,i,j) }
  temp   = matrix(0,n.ahead, n.var) 
  temp.A = temp
  for (i in 1:n.var)
    for (h in 1:n.ahead){
      temp[h,i]   =  sum(Girf$GIRF[h,i,]^2)
      temp.A[h,i] =  sum(Girf$GIRF.A[h,i,]^2)
    }  
  for (i in 1:n.var){
    GFEVD.den[,i]   = cumsum(temp[,i])   
    GFEVD.A.den[,i] = cumsum(temp.A[,i])   
  }  
  for (i in 1:n.var)
    for (j in 1:n.var){ # Eq 9 de Lanne & Nymberg (2016)
      GFEVD[,i,j]   = cumsum(Girf$GIRF[,i,j]^2) / GFEVD.den[,i] 
      GFEVD.A[,i,j] = cumsum(Girf$GIRF.A[,i,j]^2) / GFEVD.A.den[,i] 
    }
  
  #--- Plot ---#
  if(!is.null(plot.out)){
    if (plot.out=='FEVD')      { Plot.FEVD = FEVD;    Title = 'FEVD' }
    if (plot.out=='GFEVD')     { Plot.FEVD = GFEVD;   Title = 'GFEVD'}
    if (plot.out=='cum.FEVD')  { Plot.FEVD = FEVD.A;  Title = 'FEVD acum.'  }
    if (plot.out=='cum.GFEVD') { Plot.FEVD = GFEVD.A; Title = 'GFEVD acum.' }
    data      = Plot.FEVD
    name.vars = vars.to.resp
    for (variable in name.vars){
      x11()
      # par(mfrow=c(1,1))
      barplot(t(data[,variable,vars.to.shock]), main=paste0(Title, " for ",variable),  xlab="Horizon",
              ylab="Proportion", col=gray(seq(0,1,length.out=(dim(data)[2]+1))[-dim(data)[2]-1]),
              beside=FALSE, args.legend = list(xjust=0, cex = 0.5), legend = vars.to.shock)
      
      # ggplot(data=t(data[,variable,vars.to.shock]), aes(x=Horizon, y=Proportion, fill=Serie)) +
      #   geom_bar(position="fill", stat="identity")
      
    }  
  }
  
  #---- Objetos de salida ----#
  list(FEVD=FEVD, GFEVD=GFEVD, FEVD.A=FEVD.A, GFEVD.A=GFEVD.A) 
  
}

#Ej de comparacion de FEVD con estimacion del paquete <vars> (estimacion en <vars> y FEVD en <vars> y en nuestro pgm)
if(0){
  require(vars)
  #---- 1 Cholesky ----#
  if(0){
    model.vars  = VAR(y=Y, p=4, type='const')
    FEVD.vars   = fevd(model.vars, impulse=colnames(Y)[1], response=colnames(Y), n.ahead=20, ortho=TRUE, cumulative=FALSE)
    x11(); plot(FEVD.vars)
    irf.1 = GIRF(results=model.vars, data = Y, n.ahead = 20, plot.out='IRF', ortho = TRUE,  std_shock = FALSE, magnitude = 1,
                      vars.to.shock=colnames(Y)[1],vars.to.resp=colnames(Y))
    FEVD.vars.1 = GFEVD(Girf=irf.1, vars.to.shock=colnames(Y), vars.to.resp=colnames(Y)[1], plot.out='FEVD')
    FEVD.vars.2 = GFEVD(Girf=irf.1, vars.to.shock=colnames(Y), vars.to.resp=colnames(Y)[1], plot.out='GFEVD')
  }
}  