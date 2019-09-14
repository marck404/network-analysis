### Proyecto Final Inferencia II
### Marcos Torres Vivanco

rm(list=ls())

## usare esta libreria para usar graficas
library(igraph)

#### funci\'on que crea una grafica de n vertices con probabilidad p
#### con una comunidad de tamano m y probabilidad q
graph_comunidad <- function(n,p,m,q){
  a <- matrix(nrow = n, ncol=n)
  for (i in 1:n) {
    a[i,i] <- 0 
  }
  e <- ((n*(n-1))/2)
  e2 <- ((m*(m-1))/2)
  e1 <- e-e2
  edges1 <- rbinom(e1,1,p)
  edges2 <- rbinom(e2,1,q)
  edges <- c(edges1,edges2)
  c <- 1
  for (i in 1:(n-1)) {
    for (j in 1:(n-i)) {
      a[i,i+j] <- edges[c] 
      a[i+j,i] <-  edges[c]
      c <- c+1
    }
  }
  g <- graph_from_adjacency_matrix(a,mode="undirected")
  return(g)
}

## funci\'on que devuelve el triangulo superior de una matriz
triangulo <- function(m){
  n <- length(m[1,])
  t <- c()
  for (i in 1:n-1) {
    t <- c(t,m[i,(i+1):n])
  }
  return(t)
}

## funcion que cuenta la cantidad de aristas en una grafica
aristas <- function(g){
  mg <- as_adjacency_matrix(g)
  mg <- as.matrix(mg)
  aris <- sum(triangulo(mg))
  return(aris)
}

### funcion para estimar p0, en una grafica aleatoria con el modelo G(N,p0)
emvp0 <- function(g){
  mg <- as_adjacency_matrix(g)
  mg <- as.matrix(mg)
  n <- length(mg[1,])
  n2 <- n*(n-1)/2
  emv <- aristas(g)/n2
  return(emv)
}

#### -------- Veamos la prueba para detectar clique con p0 conocido

## establecemos la cantidad de v\'ertices
N <- 30
## establecemos la probabilidad de la hipotesis nula
p0 <- 0.8

## calibramos la prueba con estad\'istica clique number usando bootstrap
clique_boot <- c()
for (i in 1:10000) {
  a <- sample_gnp(N,p0)
  clique.number(a)
  clique_boot[i] <- clique.number(a)
}

hist(clique_boot, freq = F,main="Histograma No. Clique, N=30, p_0=0.8", xlab = "No. clique", ylab = "", col="tomato")
quantile(clique_boot,0.95)

## verifiquemos la potencia cuando p0 y m son conocidos
M <- 1000
potencia0 <- 0
for (i in 1:M) {
  g <- sample_gnp(N,p0)
  if(clique.number(g)>=14){
    potencia0 <- potencia0+1
  }
}
pot0 <- potencia0/M


potencia1 <- 0
for (i in 1:M) {
  g <- graph_comunidad(N,p0,10,1)
  if(clique.number(g)>=13){
    potencia1 <- potencia1+1
  }
}
pot1 <- potencia1/M



### Prueba deteccion de clique cuando p0 es desconocido

gra_w_clique <- graph_comunidad(N,p0,10,1)
emvp <- emvp0(gra_w_clique)
  
  boot_emv_clique <- c()
  for (i in 1:1000) {
    a <- sample_gnp(N,emvp)
    clique.number(a)
    boot_emv_clique[i] <- clique.number(a)
  }
  
  q <- as.integer(quantile(boot_emv_clique,probs = 0.95))
  c <- clique.number(gra_w_clique)

## realic\'e una prueba de potencia y solo rechaza 0.16
## esto se debe a que NCn p0^(nC2) =1308.645
## cuando este valor va a cero la prueba es potente

## W suma de aristas

esta01 <- function(g){
  aristas(g)
}

## Ws Scan

esta02 <- function(g,m){
  mg <- as_adjacency_matrix(g)
  mg <- as.matrix(mg)
  l <- length(mg[1,])
  combin <- combn(l,m)
  scan <- c()
  for (i in 1:length(combin[1,])) {
    subgr <- combin[,i]
    scan[i] <- sum(mg[subgr,subgr])
  }
  return(max(scan)/2)
}  

## T cantidad de triangulos

esta03 <- function(g){
  mg <- as_adjacency_matrix(g)
  mg <- as.matrix(mg)
  mg3 <- mg%*%mg%*%mg
  a <- sum(diag(mg3))/6
  return(a)
}


### calculando la potencia de las estad\'isticas

## cambiamos la p0 a una menos densa
p0 <- 0.4

## calibramos la prueba W con bootstrap

bootW <- c()
for (i in 1:5000) {
  gra <- sample_gnp(N,p0)
  bootW[i] <- esta01(gra)
}

hist(bootW,main = "Densidad total aristas N=30, p_0=0.4",xlab = "Aristas",ylab="",freq = F)
w1 <- as.integer(quantile(bootW,0.95))
ww1 <- qbinom(0.95,N*(N-1)/2,p0)

hori <- 500
test <- 500

pes <- seq(p0,1,length=hori)
po1 <- c()
for (i in 1:hori) {
  pt1 <- 0
  for (j in 1:test) {
    gp <- graph_comunidad(N,p0,10,pes[i])
    if(esta01(gp)>=w1){
      pt1 <- pt1+1
    }  
  }
  po1[i] <- pt1/test
}

## calibramos la prueba de triangulos con bootstrap

bootT <- c()
for (i in 1:5000) {
  gra <- sample_gnp(N,p0)
  bootT[i] <- esta03(gra)
}

hist(bootT, main = "Densidad Total de tri\'angulos N=30, p_0=0.4", xlab = "Tri\'angulos",ylab = "",freq = T)
w2 <- as.integer(quantile(bootT,0.95))

po2 <- c()
for (i in 1:hori) {
  pt2 <- 0
  for (j in 1:test) {
    gp <- graph_comunidad(N,p0,10,pes[i])
    if(esta03(gp)>=w2){
      pt2 <- pt2+1
    }  
  }
  po2[i] <- pt2/test
}
plot(pes,po2,type = 'l',col='blue',main = "Funciones potencia de las pruebas",xlab = "Probabilidad conexi\'on comunidad" , ylab = "Prob. rechazo")
lines(pes,po1,col='tomato')
legend("topleft", inset=.02, legend=c("Triangulos", "Aristas"),
       col=c("blue", "tomato"), lty=c(1,1), cex=1)

