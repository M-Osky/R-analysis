#Ver el directorio de trabajo
getwd()
#Cambiarlo al que contenga el input file
setwd("C:/Users/Mosky/Documents/AdegenetADE4/Files")
#comprobamos
getwd()


#Algunas librerías que podrían ser necesarias
library("ape")
library("genetics")
library("Rgenetics")
library("pegas")
library("Seqinr")
library("ggplot2")
library("adegenet")
library("spdep")
library ("poppr")

#Podemos Crear un objeto (archivo conteniendo los datos) en el formato que suelen usar
#las funciones de Adegenet a partir de un input de Genpop
#Genind13 <- import2genind("13loci13popMixID.gen", missing=NA, quiet=FALSE)

#Luego habría que poner coordenadas en un Genind
#Genind13$other$xy <- XYunic
#Genind13$other$xy
#Genind13

#Pero es más facil hacerlo a través de un Genalex con coordenadas y nombres de poblaciones
Genind13 <- read.genalex("genalexfin.csv", ploidy = 2, geo = TRUE, region = FALSE, genclone = FALSE, sep = ";")

#Y es exportable: Transformar Genind a Genalex
#genind2genalex(Genind13, filename = "genalex1.csv", quiet = FALSE, geo=TRUE)

#Vemos la Info
Genind13

Genind13$loc.names
Genind13$ind.names
Genind13$pop.names
Genind13$pop
plot(Genind13$other$xy)

#Creamos un archivo resumido, parece que algunas funciones lo usan
ResumGenind <- summary(Genind13)
#Podemos ver los descriptores
ResumGenind



#Podemos transformarlo a otro formato que puede ser también necesario
Genpformat13 <- genind2genpop (Genind13, quiet=FALSE)
#Comprobar si un "objeto" tiene el formato requerido
is.genpop(Genpformat13)
#Vemos el archivo/objeto
Genpformat13

#Creamos tambien un resumen de este archivo
ResumGenpop <- summary(Genpformat13)
#Podemos ver ahora la info del objeto
ResumGenpop

#Podemos generar un conjunto que contenga los datos separados por loci
perLoci <- seploc(Genpformat13)
#Comprobamos algunas etiquetas
class(perLoci)
names(perLoci)

#Generamos otro conjunto con los datos separados por población
perPop <- seppop(Genind13)
#Lo mismo, comprobaciones
class(perPop)
names(perPop)
perPop$Almirez

#Una vez separados por poblacion podemos separarlos por cada loci
# Tendremos así un conjunto que tiene la Info de cada población para cada marcador por separado
perUnit <- lapply(perPop, seploc)
#Comprobaciones
names(perUnit)
class(perUnit$Almirez)
names(perUnit$Almirez)
perUnit$Almirez$Pan16

#Hacemos un sencillo gráfico
#Primero hay que ver como están configurados por defecto el tamaño de los ejes
par("mar")
#Ponemos unos valores más adecuados
par(mar=c(3,3,3,3))
#Le pedimos un gráfico con los descriptores contenidos en el objeto resumen del objeto Genind
#Numero de alelos relacionado con el número de muestras
#Habrá que establecer los límites de los ejes de acuerdo con los valores
#decirle que ponga puntos (type=p) http://www.statmethods.net/graphs/line.html 
#el tipo de puntos que queremos (pch=1) http://www.statmethods.net/advgraphs/parameters.html, el color, etc. 
plot(ResumGenind$pop.eff, ResumGenind$pop.nall, ylim=c(40,100), xlab="Localities sample size", ylab="Number of alleles", main="Alleles numbers and Sample sizes", type="p", pch=1, col="blue")
#Opciones de las etiquetas, Queremos que las saque de la info de pop.eff
#Queremos que la posición (pos) y el tamaño (cex) sean diferentes a los que pone por defecto
text(ResumGenind$pop.eff, ResumGenind$pop.nall, lab=names(ResumGenind$pop.eff), pos=1, cex=0.5)
#Como quisquilloso y desconfiado, compruebo que los valores usados están bien
ResumGenind$pop.nall
ResumGenind$pop.eff

#Algunos tests que se pueden hacer con la Heterozigosidad almacenada en los "resúmenes"
bartlett.test(list(ResumGenind$Hexp, ResumGenind$Hobs))
t.test(ResumGenind$Hexp,ResumGenind$Hobs,pair=T,var.equal=TRUE,alter="greater")

#HWE
HWEqTest <- HWE.test.genind(Genind13, res="matrix")
warnings()
#comprobaciones
dim(HWEqTest)
colnames(HWEqTest)
#Ver cuales han dado significativos
idx <- which(HWEqTest<0.001, TRUE)
#Si las filas son las poblaciones y las columnas los marcadores
idx

#Calculamos las distáncias genéticas (Pairwise Fst) y las guardamosen un objeto
pairFst <- pairwise.fst(Genind13, pop=NULL, truenames=TRUE, res.type="dist")
#Vemos el objeto
pairFst
#Comprobamos que es euclidea
is.euclid(pairFst)

#Calculamos los coeficientes de inbreeding por individuo
inbreedcoef <- inbreeding(Genind13)
#Generamos un gráfico
Fbar <- sapply(inbreedcoef, mean)
hist(Fbar, col="firebrick", main="Average inbreeding")
#Le pedimos que nos muestre los individuos con un coef > 0.4
which(Fbar>0.4)

#Vamos a los PCAs!!
data(Genind13)
#Vemos cuantos missing values hay (para cada alelo)
sum(is.na(Genind13$tab))
#Sustituimos los missing por la frecuencia media y formamos una matriz escalada
PCAmatrix <- scaleGen(Genind13, missing="mean")
#Comprobamos que tipo de fichero es
class(PCAmatrix)
#Comprovamos sus dimensiones
dim(PCAmatrix)
#Quitamos la escala y comenzamos con el PCA
pca1 <- dudi.pca(PCAmatrix,cent=FALSE,scale=FALSE,scannf=FALSE,nf=3)
#dibujamos un histograma que nos muestre la importancia de los diferentes eigenvalues
barplot(pca1$eig[1:50],main="PCA eigenvalues", col=heat.colors(50))
pca1
#dibujamos un gráfico con la posición de las muestras según los prales eigenvalues
s.label(pca1$li)
title("PCA regulero del input con 13 micros \ ejes 1 y 2")
add.scatter.eig(pca1$eig[1:20], 3,1,2)
#Uno un poco más detallado
s.class(pca1$li, pop(Genind13))
title("PCA regulero del input con 13 micros II  \ ejes 1 y 2")
add.scatter.eig(pca1$eig[1:20], 3,1,2)
#Seguimos variando la representación gráfica
s.class(pca1$li,pop(Genind13),xax=1,yax=3,sub="PCA 1-3",csub=2)
title("PCA regulero del input con 13 micros III  \ ejes 1 y 3")
add.scatter.eig(pca1$eig[1:20],nf=3,xax=1,yax=3)
#vamos a ver si dándole color y transparencia es más fácil detectar patrones
col <- funky(15)
s.class(pca1$li, pop(Genind13),xax=1,yax=2, col=transp(col,.6), axesell=FALSE,
        cstar=0, cpoint=3, grid=FALSE)
#Ahora sin poblaciones, sólo coloreado según los grupos formados
colorplot(pca1$li, pca1$li, transp=TRUE, cex=3, xlab="PC 1", ylab="PC 2")
title("PCA sin etiquetas con colores")
abline(v=0,h=0,col="grey", lty=2)

#Correspondance Analisis
#Vamos a reemplazar los missing por su Chi^2
data(Genind13)
Genpop13chi <- genind2genpop(Genind13,missing="chi2")
#Gráfico de los eigenvalues para el CA
ca1 <- dudi.coa(as.data.frame(Genpop13chi$tab),scannf=FALSE,nf=3)
barplot(ca1$eig,main="Correspondance Analysis eigenvalues",
        col=heat.colors(length(ca1$eig)))
#El CA en un gráfico de dispersión, ejes 1 y 2
s.label(ca1$li,lab=Genpop13chi$pop.names,sub="CA 1-2",csub=2)
add.scatter.eig(ca1$eig,nf=3,xax=1,yax=2,posi="bottomright")
#Ejes 1 y 3
s.label(ca1$li,xax=1,yax=3,lab=Genpop13chi$pop.names,sub="CA 1-3",csub=2)
add.scatter.eig(ca1$eig,nf=3,xax=2,yax=3,posi="bottomright")

#Analisis espacial
#Vamos a intertar añadiendo las coordenadas de las poblaciones
XYpop <- read.table("popUTM.txt")
# Comprobamos que está bien
head(XYpop)
dim(XYpop)
#Vemos la distribución de los puntos
plot(XYpop)
text(XYpop,label=Genpformat13$pop.names, pos=1)

#IsolationbyDistance
#Convertimos los missing values en "0" en un archivo formato genpop
Genpmiss0 <- genind2genpop(Genind13, miss="0")
#distancias geneticas
Dgen <- dist.genpop(Genpmiss0, method=2)
#distancias geograficas
Dgeo <- dist(XYpop)
#IbD
ibd <- mantel.randtest(Dgen, Dgeo)
ibd
#Vemos el valor real contra los simulados
plot(ibd)
#Si el valor está claramente fuera de la distribución podemos afirmar que hay estructura
#Pero es una clina continua o tenemos distintor parches?
#Representamos las distancias
plot(Dgeo, Dgen)
abline(lm(Dgen~Dgeo), col="red",lty=2)
#Si no vemos una nube más o menos continua de puntos
# agrupados alrededor de un nucleo es que no es claramente una clina
#Puede quedar mas claro con colores
dens <- kde2d(Dgeo,Dgen, n=150)
myPal <- colorRampPalette(c("white","blue","gold", "orange", "red"))
plot(Dgeo, Dgen, pch=20,cex=.5)
image(dens, col=transp(myPal(150)), add=TRUE)
abline(lm(Dgen~Dgeo))
title("Isolation by distance plot")

#Ahora una aproximación parecida a la Voronoi Tesselation:
# Monmonier's algorithm

#El tutorial hace antes un gráfico no entiendo muy bien para qué.
#Se supone que lo hace con coordenadas pero los valores van de 0 a 100 en ambos ejes
#Probamos con las coordenadas únicas

XYunic <- Genind13$other$xy
# Comprobamos que está bien
head(XYunic)
dim(XYunic)
#Vemos la distribución de los puntos
plot(XYunic)

#Ahora dibujamos las muestras en un su distribucion
Genind13$pop
poptab <- Genind13$pop
levels(poptab) <- c(1,2,3,4,5,6,15,16,17,18,21,24,0)
poptab <- as.numeric(as.character(poptab))
plot(Genind13$other$xy,pch=poptab,cex=1.5,xlab='x',ylab='y')
Genind13$pop.names
#legend("bottomright",leg=c(Genind13$pop.names),pch=c(1,2,3,4,5,6,0,16,17,18,21,24,22))
#No ponemos leyenda, porque tapa todo el mapa
#Habría que configurar otro gráfico en blanco y poner ahi la leyenda

args(monmonier)
#Nos dirá lo que necesita:
#function (xy, dist, cn, threshold = NULL, bd.length = NULL, nrun = 1, 
#skip.local.diff = rep(0, nrun), scanthres = is.null(threshold), 
#allowLoop = TRUE) 
#NULL
#xy -> Son las coordenadas
#dist -> Son las distáncias genéticas
D <- dist(Genind13$tab)
#cn -> Connection network. Lo que serían la red de conexión entre los diferentes puntos
#parece algo muy complicado, así que en el tutorial lo simplifican
gab <- chooseCN(XYunic,ask=FALSE,type=2)
#Una vez está la network hecha procedemos al analisis
mon1 <- monmonier(XYunic, D, gab)
#Vemos las diferencias locales ordenadas en orden decreciente.
#Si las diferencias no caen abruptamente después del límite quiere decir que
#o no hay estructura o que no puede sobreponerse al "ruido" de las distancias genéticas

#Observamos las distancias geneticas
pairFst <- pairwise.fst(Genind13, pop=NULL, truenames=TRUE, res.type="dist")

#deberíamos considerar las Fst como el valor de una diferenciación débil
#pero para ver si esa diferenciacion es significativa
#deberíamos realizar unas permutaciones para ver si es mayor o menor que ellas

#replicate(200, pairwise.fst(Genind13, pop=sample(pop(Genind13))))
#ponemos 10 por ser rápido, pero con unas 200 tendremos un mejor p-valor

#suponemos que si (eso ya lo hemos hecho en Genalex y otros softwares)
#Y vamos a ver si la estructura está enmascarada por el ruido o si no hay

#Una forma de librarse del ruido sería mediante un PCA

#La forma de hacerlo sería así
#pco1 <- dudi.pco(D,scannf=FALSE,nf=1)
#pero da un error de "non euclidean distances" así que probamos otro método
#aunque transforma un poco los valores

pco1 <- dudi.pco(cailliez(as.dist(D)), scannf=FALSE)

barplot(pco1$eig,main="Eigenvalues")
pco1
#Redefinimos las distancias geneticas sólo con los componentes del principal eigenvalues

D <- dist(pco1$li)

#corremos ahora el analisis con ese valor

mon1 <- monmonier(XYunic, D, gab)
d
mon1
#Los resultados contienen información sobre los límites entre los diferentes puntos
names(mon1)
names(mon1$run1)
#Encontramos las coordenadas de los puntos con barreras y los valores
mon1$run1$dir1
#podemos ver a qué puntos afectan las barreras
coords.monmonier(mon1)
#Nos dará las coordenadas de los puntos en la barrera y el número que identifica
# Los dos puntos cuyo baricentro está en la barrera

#para ver esto
plot(mon1)
plot(mon1,add.arrows=FALSE,bwd=10,col="black")
points(XYunic, cex=2, pch=20,col=fac2col(pop(Genind13), col.pal=spectral))
#legend("bottomright",leg=c("Pop A", "Pop B"), pch=c(20),col=spectral(2), pt.cex=2)
#dejamos la leyenda, que hay que configurarla mucho y luego no va a caber bien

#hay que aprender a interpretar bien este output

#Por lo demás el 1er tutorial (general) de Adegenet ya está =·)
