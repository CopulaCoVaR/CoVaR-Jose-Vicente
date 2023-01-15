#-------------------------------- Funcion GIRF  --------------------------------------
# Calcula los IRF y GIRF. Ambos en niveles y acumulados basados en objetos <BigVAR> o <VARS>
# Los GIRF se calculan de acuerdo a Pesaran y Shin (1998)
# ---------------- Parámetros de entrada de la función  ------------------------------
# <results>:   Objeto class <BigVAR.results> del paquete <BigVAR> o del paquete <vars>
# <data> :     Matriz de datos 
# <p> :        orden p del VAR. Si es <null> se obtiene de <results>
# <n.ahead>:   Número de pasos adelante del impulso respuesta
# <ortho>:     <TRUE>:  Se calculan las impulso respuesta suponiendo el ordenamiento 
#                       de exogenidad contemporanea dado por Cholesky. 
#              <FALSE>: Se calculan las impulso respuesta SIN realizar la transformacion 
#                       de ortogonalidad.
#              Esta opcion solo afecta a las IRF y no tiene incidencia sobre las GIRF.
# <std_shock>: <TRUE>:  calcula las impulso respuesta con un choque de una desviación 
#                       estándar.
#              <FALSE>: Las calcula con un choque determinado por <magnitude>, 
#                       por defecto 1. 
# <magnitude>: Variable numérica que indica la magnitud del choque para calcular 
#              las impulso respuesta. 
# <plot.out>:       Indica que tipo de impulsos-respuesta se grafican. Valores:
#              <'IRF'>     :  Para IRF en niveles
#              <'GIRF'>    :  Para GIRF en niveles
#              <'cum.IRF'> :  Para IRF acumulada
#              <'cum.GIRF'>:  Para GIRF acumulada
# <vars.to.shock>:  Labels de las variables que se desean chocar
# <vars.to.resp> :  Labels de las variables sobre las cuales se desean ver las rtas
# <pdf>       : FALTA!!!!!   
# <x11>       : FALTA!!!!!
# <grid.dims> : FALTA!!!!!
#-----------------------------------------------------------------------------------
#---------------------- Objetos de salida de la función ----------------------------
# IRF:    Arreglo de tres dimensiones con los IRF (<n.ahead> x <n.var> x <n.var> )
#         [ (horizonte de los IRF) x (variables rta) x (variables chocadas)  ] 
# GIRF:   Arreglo de tres dimensiones con los IRF Generalizados
# IRF.A:  Arreglo de tres dimensiones con los IRF acumulados
# GIRF.A: Arreglo de tres dimensiones con los GIRF acumulados
#-----------------------------------------------------------------------------------
GIRF_BigVAR <- function(results, data, p=NULL, n.ahead = 20, ortho = TRUE,  std_shock = TRUE, 
                        magnitude = 1, plot.out=NULL, vars.to.shock=NULL,vars.to.resp=NULL,
                        pdf = FALSE, x11 = TRUE, grid.dims = c(1,1)){
  
  y.names  = colnames(data)
  impulse  = y.names
  response = y.names
  n.ahead = abs(as.integer(n.ahead))
  n.var   = ncol(data)
  # Crear la estructura del <IRF> y <GIRF> para los calculos     
  # [1:nlags, resp. vars, shocked vars] 
  IRF  = array(data = 0, dim = c(n.ahead, n.var, n.var), dimnames = list(1:n.ahead,y.names,y.names))       
  GIRF = IRF
  
  if (class(results) == 'BigVAR.results'){
    Phi = results@betaPred
    if (results@intercept) Phi = Phi[,-1] #Se quita la columna del intercepto 
    if (is.null(p))        p   = results@lagmax
    sigma.u = crossprod(results@resids)/ nrow(results@resids) #ML
    #sigma.u = crossprod(results@resids)/(nrow(results@resids) - (sum(Phi!=0)/n.var) - (results@intercept)) #LS
  }  
  if (class(results) == 'varest'){
    Phi = NULL
    for (ii in 1:n.var)   Phi = rbind(Phi, coef(results$varresult[[ii]]))
    if (results$type=='const') Phi = Phi[,-ncol(Phi)] #Se quita la columna del intercepto 
    if (is.null(p))            p   = results$p
    sigma.u = crossprod(residuals(results))/nrow(residuals(results)) #ML
    #sigma.u = crossprod(residuals(results))/(nrow(residuals(results)) - n.var*p - (results$type=='const')) #LS
    #sigma.u = summary(results)$covres
  }  
  # Estimación de "orthogonalized" y "generalized" IRFs
  # Representacion VAR(1) de un VAR(p)
  if(p>1)  Phi.VAR1 <- VarptoVar1MC(Phi, p, n.var)
  else     Phi.VAR1 <- Phi
  
  #params  = ncol(data) #datamat son los datos, K de variables del var
  J       = matrix(0, nrow=n.var, ncol=n.var*p)
  diag(J) = 1
  P       = t(chol(sigma.u))
  sig     = diag(sigma.u)
  if (std_shock == TRUE)   shock = sqrt(sig) else shock = magnitude
  if (ortho     == TRUE)   PP    = P         else PP    = diag(n.var)   #<ortho==TRUE>  asume orden de exo. contemporanea dado por Cholesky <ortho==FALSE> no hace transformaciones
  if (ortho     == TRUE)   shock = 1 # Cuando se usan choques ortogonales-Cholesky, el choque es aprox. una SD, 
  # => <shock = 1> para no multiplicarlo por nada mas (por lo tanto la opc. de <std_shock> no aplica)
  # Loop para cada variable 
  for (jj in 1:n.var){ # Choque en las vars jj
    # Vector índice
    indx_       = matrix(0,n.var,1)
    indx_[jj,1] = 1   
    # Loop para cada horizonte del IRF 
    for (kk in 0:(n.ahead-1)){  ## kk counts the horizon
      phi1              = J %*% (Phi.VAR1%^%kk) %*% t(J) # Coef MA en el rezago kk
      IRF [(kk+1), ,jj] = phi1 %*% PP %*% indx_ * shock  # Pesaran-Shin (1998) eq 7 (OIRF)
      GIRF[(kk+1), ,jj] = phi1 %*% sigma.u %*% indx_ * (sig[jj]^(-1)) * shock     # delta_j = sig[jj]^(0.5) # Pesaran-Shin (1998) eq 10 (GIRF)
    }
  }
  
  #Impulsos acumulados
  IRF.A  = IRF-IRF
  GIRF.A = GIRF-GIRF
  for(jj in 1:n.var)
    for(kk in 1:n.var){
      IRF.A[,jj,kk]  = cumsum(IRF[,jj,kk])
      GIRF.A[,jj,kk] = cumsum(GIRF[,jj,kk])
    }  
  
  # Grafica del GIRF asociado a <vars.to.shock>  y <vars.to.resp>  
  if(!is.null(plot.out))
  {  
    if (length(vars.to.shock)*length(vars.to.resp)<9) wind = c(4,2)  else wind = grid.dims #Falta generalizar!!
    if (plot.out=='IRF')       { Plot.IRF = IRF;    Title = 'IRF en niveles' }
    if (plot.out=='GIRF')      { Plot.IRF = GIRF;   Title = 'GIRF en niveles'}
    if (plot.out=='cum.IRF')   { Plot.IRF = IRF.A;  Title = 'IRF acumulada'  }
    if (plot.out=='cum.GIRF')  { Plot.IRF = GIRF.A; Title = 'GIRF acumulada' }
    ngraphs = 0
    max.graph = wind[1]*wind[2]
    for(jj in vars.to.shock)
      for(kk in vars.to.resp){
        ngraphs = ngraphs + 1
        if (((ngraphs - 1)%% max.graph/(grid.dims[1]*grid.dims[2]))==0) {
          if (pdf == TRUE) {
            if (dev.cur() > 1) dev.off() #Cierra todos menos el último
            pdf(file = paste('Impulso',kk,jj,'.pdf'))
          }
          if (x11 == TRUE) x11()
          par(mfrow=wind, mar=c(5-3, 4, 4-2, 2) + 0.1, oma=c(0,0,3,0))
        }
        plot(1:n.ahead,Plot.IRF[,kk,jj], main=paste('Choque en ',jj,'Rta en ',kk), xlab='',ylab='IRF',type='l', col='steelblue', lwd=2) #Falta verif. choque y rta
        abline(h = 0, col='tomato', lty = 'dotted', lwd = 2)
        mtext(Title, side=3, cex =1.3, line = 0.5, outer=TRUE)
      }
    if (pdf ==TRUE ) graphics.off() # Cierra cualquier device que haya quedado abierto
  }
  #---- Objetos de salida ----#
  list(IRF=IRF, GIRF=GIRF, IRF.A=IRF.A, GIRF.A=GIRF.A)
}
# Ej de comparacion IRF vs GIRF para estimacion por paquete <BigVAR>
if(0){
IRF1 = GIRF(results=Model1Results, data=Y, n.ahead=20, plot.out='IRF', ortho=TRUE,  std_shock=TRUE, magnitude=1,
            vars.to.shock=colnames(Y), vars.to.resp=colnames(Y))

IRF2 = GIRF(results=Model1Results, data=Y, n.ahead=20, plot.out='GIRF', ortho=FALSE,  std_shock=FALSE, magnitude=1,
            vars.to.shock=colnames(Y), vars.to.resp=colnames(Y))
}

# Ej de comparacion de IRFs con estimacion del paquete <vars> (estimacion en <vars> e IRF en <vars> y en nuestro pgm)
if(0){
  require(vars)
  #---- 1 Cholesky ----#
  if(0){
   model.vars = VAR(y=Y, p=4, type='const')
   irf.vars   = irf(model.vars, impulse=colnames(Y)[1], response=colnames(Y), n.ahead=20, ortho=TRUE, cumulative=FALSE, boot=FALSE, ci=50)
   x11(); plot(irf.vars)
   IRF.vars   = GIRF(results=model.vars, data = Y, n.ahead = 20, plot.out='IRF', ortho = TRUE,  std_shock = FALSE, magnitude = 1,
                    vars.to.shock=colnames(Y)[1],vars.to.resp=colnames(Y))
  }
  #--- 2 sin Cholesky, choque=1 ----#
  if(0){
    model.vars = VAR(y=Y, p=4, type='const')
    irf.vars   = irf(model.vars, impulse=colnames(Y)[1], response=colnames(Y), n.ahead=20, ortho=FALSE, cumulative=FALSE, boot=FALSE, ci=50)
    x11(); plot(irf.vars)
    IRF.vars   = GIRF(results=model.vars, data = Y, n.ahead = 20, plot.out='IRF', ortho = FALSE,  std_shock = FALSE, magnitude = 1,
                      vars.to.shock=colnames(Y)[1],vars.to.resp=colnames(Y))
  }
  #---- 3 Acum-Cholesky ----#
  if(0){
    model.vars = VAR(y=Y, p=4, type='const')
    irf.vars   = irf(model.vars, impulse=colnames(Y)[1], response=colnames(Y), n.ahead=20, ortho=TRUE, cumulative=TRUE, boot=FALSE, ci=50)
    x11(); plot(irf.vars)
    IRF.vars   = GIRF(results=model.vars, data = Y, n.ahead = 20, plot.out='cum.IRF', ortho = TRUE,  std_shock = FALSE, magnitude = 1,
                      vars.to.shock=colnames(Y)[1],vars.to.resp=colnames(Y))
  }
  
}

# Ej de comparacion de IRFs usando paquetes <vars> y <BigVAR> (tambien se puede comparar con el PRM <Ej_BigVAR.R>)
if(0){
  model.vars = VAR(y=Y, p=4, type='const')
  IRF1 = GIRF(results=model.vars,    data=Y, n.ahead=20, plot.out='IRF', ortho=TRUE,  std_shock=FALSE, magnitude=1,
              vars.to.shock=colnames(Y)[1], vars.to.resp=colnames(Y))
  
  IRF2 = GIRF(results=Model1Results, data=Y, n.ahead=20, plot.out='IRF', ortho=TRUE,  std_shock=FALSE, magnitude=1,
              vars.to.shock=colnames(Y)[1], vars.to.resp=colnames(Y))
}
