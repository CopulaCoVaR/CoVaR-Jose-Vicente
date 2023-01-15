# ----------------------- Función outliersReplace ---------------------
# Función que permite reemplazar datos atípicos de una base de datos. 
#
# ----------------------- Argumnetos de entrada -----------------------
# <data>:      Matriz de datos
# <lowLimit>:  Límite inferior. Un dato menor a <lowLimit> es atipico.
# <highLimit>:  Límite superior. Un dato superior a <highLimit> es atipico.
#
# ----------------------- Obejetos de salida --------------------------
# Matriz de datos con corrección de valores atípicos. 
#----------------------------------------------------------------------
outliersReplace <- function(data, lowLimit, highLimit){
  data[data<lowLimit] = mean(data)
  data[data>highLimit]= median(data)
  # data[data<lowLimit] = NA
  # data[data>highLimit]= NA
  # data = na.fill(data, 'extend')
  return(data)     #devolvemos el dato       
}
#--------------------- Ejemplo ----------------------------------------
Datos.corregir=CoVaR_data_Vix$CoVaR
# Grafica inicial para determinar un límite inferior y superior válido. 
x11()
plot(Datos.corregir)
low=-500
high=500
# Ejecutamos la función,
Prueba_Outliers=outliersReplace(data=Datos.corregir, lowLimit=low, highLimit=high) 
# Vemos diferencias gráficas. 
x11() 
par(mfrow=c(2,1))
plot(Datos.corregir)
plot(Prueba_Outliers)


