#--------------------- Calculo del CoVaR -------------------#
# 1. Se calcula Fy(CoVaR_z) o F1(CoVaR_z) Eqs. (12) de RU (2016) y (22) de Liu (2022)
# 2. Se calcula el CoVaR estandarizado:  Fy^(-1)[Fy(Covar_z)] 
# 3. Se calcula EL CoVaR.
#------------------------------------ Argumentos de entrada ------------------------------------------#
# y:          Nivel de significancia asociada al CoVaR (ALPHA)
# x:          Nivel de significancia asociada al  VaR (BETA)
# par:        Primer parámetro asociado a la cópula.
# par2:       Segundo parámetro asociado a la cópula, si esta lo requiere.  
# dof:        Grados de libertad de la distribucion marginal t (de <dist>)
# gamma:      Coeficiente de asimetría de la distribucion marginal t (de <dist>)
# cond.mean:  Serie de la media condicional marginal sobre la serie analizada (country' stock index)
#             Se obtiene del modelo ARMA 
# cond.sigma: Serie de la STD  condicional marginal sobre la serie analizada (country' stock index)
#             Se obtiene del modelo GARCH
# dist:       Nombre de la distribucion marginal, se asume que es t-asimetrica
# type:       Tipo de copula usada.
# COND:       Asociado al evento condicional en la definicion del CoVaR. Toma dos valores:
#             <'Less'>:  Si la variable estresada es menor o igual que el cuantil <ALPHA>-esimo
#             <'Equal'>: Si la variable estresada es igual al cuantil <ALPHA>-esimo
#------------------------------------------------------------------------------------------------------#
# ARGUMENTOS DE SALIDA: <CoVaR>:    CoVaR,  <alpha>:   Fy(CoVaR_z) y <perc>:    Covar_z  [For the left tail]
#                       <CoVaR.up>: CoVaR, <alpha.up>: Fy(CoVaR_z) y <perc.up>: Covar_z  [For the right tail]
#------------------------------------------------------------------------------------------------------#
CoVaR.LF <- function(y, x, par, par2, dof, gamma, cond.mean, cond.sigma, dist,type, COND='Equal')
{
  level.CoVaR = y
  level.VaR   = x
  
  #--- Fy(CoVaR) o F1(CoVaR) ---#
  cond = NULL; #cond.up = NULL
  par  = as.matrix(par)
  N    = length(par)
  fam  <- switch (type,
            Gumbel        = 4,
            Dyn.Gumbel    = 4,
            Rot.Gumbel    = 14,
            Dyn.Rot.Gumbel= 14,
            Clayton       = 3,
            Rot.Clayton   = 13,
            Frank         = 5,
            Joe           = 6,
            Gaussian      = 1,
            Dyn.Gaussian  = 1,
            Student       = 2,
            Dyn.Student   = 2,
            BB1           = 7,
            BB7           = 9,
            Rot.BB7       = 19,
            BB6           = 8,
            BB8           = 10, 
            Dyn.BB7       = 9,
            Plackett      = 20,
            AMH           = 21,
    )       
    #Plackett     = Plackett_CoVaR(level.CoVaR = y , level.VaR = x, par = par),
  #SCJ          = SCJ_CoVaR(level.CoVaR = y , level.VaR = x, par = par, par2 = par2)
  if (COND == 'Less'){#Paper Reboredo & Ugolini 2016
    if (fam==20) cond = Plackett_CoVaR(level.CoVaR = y , level.VaR = y, par = par)
    if (fam==21) cond = amh_CoVaR(level.CoVaR = y , level.VaR = x, par = par)
    if (fam < 20){
      for (i in 1:N){
        CoVaR=function(x){ 
          #pCopula(c(x,level.VaR), normalCopula(par[i]))-(level.CoVaR*level.VaR)
          # Supuesto: El parametro <par2> nunca cambia en el tiempo
          if (is.null(par2)) COPula = BiCop(family=fam, par=par[i]) else
            COPula = BiCop(family=fam, par=par[i], par2=par2)
          BiCopCDF(u1=x, u2=level.VaR,  COPula) - (level.CoVaR*level.VaR) 
        }  
        cond[i] = fzero(CoVaR,0, tol = 10^-16)$x
      } 
    }
  }
  if (COND == 'Equal'){ #Paper High-Dimensional - Liu et al. 2022
    for (i in 1:N){
      # Supuesto: El parametro <par2> nunca cambia en el tiempo
      if (is.null(par2)) COPula = BiCop(family=fam, par=par[i]) else
        COPula    = BiCop(family=fam, par=par[i], par2=par2)
      cond[i]     = BiCopHinv2(level.CoVaR, level.VaR, COPula) #Inverse Conditional Distr. Function of a Bivar. Copula
      #cond.up[i] = BiCopHinv2(1-level.CoVaR, 1-level.VaR, COPula ) #Inverse Conditional Distr. Function of a Bivar. Copula
    }
  }
  alpha = cond; #alpha.up = cond.up
  alpha.up = 1 - alpha
  
  #-- Standarized CoVaR ( Fy^(-1)[Fy(Covar)]) --#
  perc = switch(dist,
                t      = as.matrix(qt(alpha,df=as.numeric(dof))),
                gauss  = as.matrix(qnorm(alpha)),
                #tskew = skewtdis_inv(alpha,nu=dof,lambda=gamma) #Lambda: Grado de asimetría de la t
                tskew  = rugarch::qdist(distribution='sstd', p=alpha, shape=dof, skew=gamma)
         )
  perc.up = switch(dist,
                t     = as.matrix(qt(alpha.up,df=as.numeric(dof))),
                gauss = as.matrix(qnorm(alpha.up)),
                tskew = rugarch::qdist(distribution='sstd', p=alpha.up, shape=dof, skew=gamma)
  )
  #-- CoVaR  ---#
  CoVaR = NULL; CoVaR.up = NULL
  #---- 
  CoVaR               = cond.mean + cond.sigma * as.numeric(perc)
  VaR                 = as.xts(cond.mean + cond.sigma*(rugarch::qdist(distribution='sstd', p=y,   shape=dof,skew=gamma)))
  cond.mean.corrected = cond.mean
  for (i in which(CoVaR>0)) {
    a=i-1
    b=i-1
    if (length(perc.up)!=1) CoVaR.alt = as.numeric(cond.mean[a]) + cond.sigma[i] * as.numeric(perc[i])
    if (length(perc.up)==1) CoVaR.alt = as.numeric(cond.mean[a]) + cond.sigma[i] * as.numeric(perc)
    while (CoVaR.alt>0) {
      a=a-1
      if (a<=0) {
        CoVaR.alt = CoVaR[b]
        VaR.alt   = VaR[b]
        b=b-1
      }else{
        if (length(perc.up)!=1)CoVaR.alt=as.numeric(cond.mean[a]) + cond.sigma[i] * as.numeric(perc[i])
        if (length(perc.up)==1)CoVaR.alt=as.numeric(cond.mean[a]) + cond.sigma[i] * as.numeric(perc)
      }
    }
    if (a>0) cond.mean.corrected[i]=as.numeric(cond.mean[a])
    CoVaR[i]         = CoVaR.alt
    if (a<=0) VaR[i] = VaR.alt
  }
  #CoVaR           = cond.mean.corrected + cond.sigma * as.numeric(perc)
  VaR             = as.xts(cond.mean.corrected + cond.sigma*(rugarch::qdist(distribution='sstd', p=y,   shape=dof,skew=gamma)))
  colnames(CoVaR) = "CoVaR"
  colnames(VaR)   = "VaR"
 #----
  CoVaR.up               = cond.mean.corrected + cond.sigma * as.numeric(perc.up)
  VaR.up                 = as.xts(cond.mean.corrected + cond.sigma*(rugarch::qdist(distribution='sstd', p=1-y, shape=dof,skew=gamma)))
  cond.mean.corrected.up = cond.mean.corrected
  for (i in which(CoVaR.up<0)){
    a=i-1
    b=i-1
    if (length(perc.up)!=1) CoVaR.up.alt = as.numeric(cond.mean.corrected[a]) + cond.sigma[i] * as.numeric(perc.up[i])
    if (length(perc.up)==1) CoVaR.up.alt = as.numeric(cond.mean.corrected[a]) + cond.sigma[i] * as.numeric(perc.up)
    while (CoVaR.up.alt<0) {
      a=a-1
      if (a<=0) {
        CoVaR.up.alt = CoVaR.up[b]
        VaR.up.alt   = VaR.up[b]
        b=b-1
      }else{
        if (length(perc.up)!=1)CoVaR.up.alt = as.numeric(cond.mean.corrected[a]) + cond.sigma[i] * as.numeric(perc.up[i])
        if (length(perc.up)==1)CoVaR.up.alt = as.numeric(cond.mean.corrected[a]) + cond.sigma[i] * as.numeric(perc.up)
      }
    }
    if (a>0) cond.mean.corrected.up[i]=as.numeric(cond.mean[a])
    CoVaR.up[i]         = CoVaR.up.alt
    if (a<=0) VaR.up[i] = VaR.up.alt
  }
  CoVaR.up           = cond.mean.corrected.up + cond.sigma * as.numeric(perc.up)
  VaR.up             = as.xts(cond.mean.corrected.up + cond.sigma*(rugarch::qdist(distribution='sstd', p=1-y, shape=dof,skew=gamma)))
  colnames(CoVaR.up) = "CoVaR.up"
  colnames(VaR.up)   = "VaR.up"
  #-- Output ---#
  list(CoVaR=CoVaR,       perc   =perc,    alpha   =alpha, 
       CoVaR.up=CoVaR.up, perc.up=perc.up, alpha.up=alpha.up,
       VaR=VaR, VaR.up=VaR.up)
}




