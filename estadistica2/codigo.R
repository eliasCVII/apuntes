# DATOS DE ENTRADA:
#
#   n TAMAÑO DE LA MUESTRA
#   prom PROMEDIO MUESTRAL
#   sigma DESVIACIÓN ESTÁNDAR POBLACIONAL
#   nc NIVEL DE CONFIANZA

# IC para mu con sigma, sigma2 conocido
ICMU1=function(n,prom,sigma,nc=0.95) 
  {
    z=qnorm(1-(1-nc)/2)
    EE=z*sigma/sqrt(n)
    LI=prom-EE
    LS=prom+EE
    cat("Un intervalo de confianza para MU","\n",
         "al nivel de confianza",nc*100,"%", "es:","\n",
    "Limite inferior:",LI,"Limite superior:",LS,"\n")
  }

# DATOS DE ENTRADA:
#
#   n TAMAÑO DE LA MUESTRA
#   prom PROMEDIO MUESTRAL
#   s DESVIACIÓN ESTÁNDAR MUESTRAL
#   nc NIVEL DE CONFIANZA

# IC para mu con sigma, sigma2 desconocido
ICMU2=function(n,prom,s,nc=0.95)
  {
    t=qt(1-(1-nc)/2,n-1)
    EE=t*s/sqrt(n)
    LI=prom-EE
    LS=prom+EE
    cat("Un intervalo de confianza para MU","\n",
         "al nivel de confianza",nc*100,"%", "es:","\n",
    "Limite inferior:",LI,"Limite superior:",LS,"\n")
  }

# DATOS DE ENTRADA:
#
#   n TAMAÑO DE LA MUESTRA
#   s DESVIACIÓN ESTÁNDAR MUESTRAL
#   nc NIVEL DE CONFIANZA


# IC para sigma2 (varianza poblacional)
ICSIGMA2=function(n,s,nc=0.95)
  {
    v1=qchisq((1-nc)/2,n-1)
    v2=qchisq(1-(1-nc)/2,n-1)
    LI=(n-1)*s^2/v2
    LS=(n-1)*s^2/v1
    cat("Un intervalo de confianza para SIGMA^2","\n",
         "al nivel de confianza",nc*100,"%", "es:","\n",
    "Limite inferior:",LI,"Limite superior:",LS,"\n")
  }

# DATOS DE ENTRADA:
#
#   n TAMAnO DE LA MUESTRA
#   p PROPORCIoN MUESTRAL
#   nc NIVEL DE CONFIANZA

# IC para pi (proporcion poblacional)
ICPI=function(n,p,nc=0.95)
  {
    z=qnorm(1-(1-nc)/2)
    EE=z*sqrt(p*(1-p)/n)
    LI=p-EE
    LS=p+EE
    cat("Un intervalo de confianza para PI","\n",
         "al nivel de confianza",nc*100,"%", "es:","\n",
    "Limite inferior:",LI,"Limite superior:",LS,"\n")
  }

# SE NECESITA:
#     EE ERROR DE ESTIMACIoN
#     nc NIVEL DE CONFIANZA
#     sigma DESVIACIoN ESTaNDAR
#

# tamano de la muestra para media poblacional mu
NMU=function(EE,nc,sigma)
  {
    z=qnorm(1-(1-nc)/2)
    n=round((z/EE)^2*sigma^2,0)
    cat("El tamano de la muestra para estimar MU","\n",
         "al nivel de confianza",nc*100,"%","con un error
 de estimacion de",EE,", es:",n,"\n")
  }

# SE NECESITA:
#     EE ERROR DE ESTIMACIoN
#     nc NIVEL DE CONFIANZA
#     p PROPORCIoN CONOCIDA
#

# tamano de muestra para proporcion poblacional
NPI=function(EE,nc,p)
  {
    z=qnorm(1-(1-nc)/2)
    n1=round((z/EE)^2*p*(1-p),0)
    n2=round(1/4*(z/EE)^2,0)
    cat("El tamano de la muestra para estimar PI","\n",
         "al nivel de confianza",nc*100,"%","con un error
 de estimacion de",EE,", es:",n1,"\n")
    cat("El tamano de la muestra para estimar PI","\n",
         "al nivel de confianza",nc*100,"%","con un error
 de estimacion de",EE,",en el peor de los casos, es:",n2,"\n")
  }

# DATOS:
# n1, n2 TAMAÑOS DE LAS MUESTRAS
# prom1, prom2 PROMEDIOS MUESTRALES
# sigma1, sigma2 DESVIACIONES POBLACIONALES
# nc NIVEL DE CONFIANZA (PROPORCIÓN)

# Diferencia de medias
## Caso1. IC para dos poblaciones, con varianzas poblacionales conocidas
difmed1=function(n1,n2,prom1,prom2,sigma1,sigma2,nc=0.95)
{
  z=qnorm(1-(1-nc)/2)
  EE=z*sqrt(sigma1^2/n1+sigma2^2/n2)
  LI= prom1-prom2-EE
  LS=prom1-prom2+EE
  cat("Un intervalo de confianza para MU1-MU2","\n",
         "al nivel de confianza",nc*100,"%", "es:","\n",
  "Limite inferior:",LI,"Limite superior:",LS,"\n")
}

# DATOS:
# n1, n2 TAMAÑOS DE LAS MUESTRAS
# prom1, prom2 PROMEDIOS MUESTRALES
# s1, s2 DESVIACIONES ESTANDAR MUESTRALES
# nc NIVEL DE CONFIANZA (PROPORCIÓN)
#

# Diferencia de medias
## Caso2. IC para dos poblaciones, con varianzas poblacionales desconocidas pero iguales
difmed2=function(n1,n2,prom1,prom2,s1,s2,nc=0.95)
     {
        t=qt(1-(1-nc)/2,n1+n2-2)
  sp=sqrt(((n1-1)*s1^2+(n2-1)*s2^2)/(n1+n2-2))
  EE=t*sp*sqrt(1/n1+1/n2)
  LI= prom1-prom2-EE
  LS=prom1-prom2+EE
  cat("Un intervalo de confianza para MU1-MU2","\n",
         "al nivel de confianza",nc*100,"%", "es:","\n",
  "Limite inferior:",LI,"Limite superior:",LS,"\n")
     }

# DATOS:
# n1, n2 TAMAÑOS DE LAS MUESTRAS
# prom1, prom2 PROMEDIOS MUESTRALES
# s1, s2 DESVIACIONES ESTANDAR MUESTRALES
# nc NIVEL DE CONFIANZA (PROPORCIÓN)
#

# Diferencia de medias
## Caso3. IC para dos poblaciones, con varianzas poblacionales desconocidas pero y diferentes
difmed3=function(n1,n2,prom1,prom2,s1,s2,nc=0.95)
     {
  eta=(s1^2/n1+s2^2/n2)^2/((s1^2/n1)^2/(n1-1)+
             (s2^2/n2)^2/(n2-1))
        t=qt(1-(1-nc)/2,eta)
  EE=t*sqrt(s1^2/n1+s2^2/n2)
  LI= prom1-prom2-EE
  LS=prom1-prom2+EE
  cat("Un intervalo de confianza para MU1-MU2","\n",
         "al nivel de confianza",nc*100,"%", "es:","\n",
  "Limite inferior:",LI,"Limite superior:",LS,"\n")
     }

# DATOS:
# s1, s2 DESVIACION ESTANDAR MUESTRALES
# n1, n2 TAMAÑOS DE LAS MUESTRAS
# nc NIVEL DE CONFIANZA (PROPORCION)
#

# Diferencia de medias
## Cuociente de varianzas poblacionales
cuovar=function(n1,n2,s1,s2,nc=0.95)
     {
  f1=qf((1-nc)/2,n1-1,n2-1)
  f2=qf(1-(1-nc)/2,n1-1,n2-1)
  cuo=s1^2/s2^2
  LI=cuo/f2
  LS=cuo/f1
        cat("Un intervalo de confianza para SIGMA1^2/SIGMA2^2",
            "\n","al nivel de confianza",nc*100,"%", "es:","\n",
  "Limite inferior:",LI,"Limite superior:",LS,"\n")
     }

# DATOS:
#   n1, n2 TAMAÑOS DE LAS MUESTRAS
#   p1, p2 PROPORCIONES MUESTRALES
#   nc NIVEL DE CONFIANZA (PROPORCIÓN)

# IC para diferencia de proporciones pi1 - pi2
difprop=function(n1,n2,p1,p2,nc=0.95){
  z=qnorm(1-(1-nc)/2)
  v1=p1*(1-p1)/n1
  v2=p2*(1-p2)/n2
  EE = z*sqrt(v1+v2)
  LI=p1-p2-EE
  LS=p1-p2+EE
  cat("Un intervalo de confianza para PI1-PI2","\n",
         "al nivel de confianza",nc*100,"%", "es:","\n",
      "Limite inferior:",LI,"Limite superior:",LS,"\n")
}
