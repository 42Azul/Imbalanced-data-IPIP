#Ejecutar todos los casos:


quant = 0.75

nmin <- nrow(muertos)
np <- ceiling(nmin*quant)
p <- ceilinglog(alpha)/(log(1-1/nmin)*np)
alpha <- .01


dfs <- list()
# Seleccionamos todos los  sanos que usaremos, y los reordenamos

for(k in 1:p){
  id.muertos <- sample(x = 1:nrow(muertos), size = np) #Índices de clase minoritaria para cada subconjunto
  id.sanos <- sample(x= 1:nrow(sanos), size = round(np*prop.mayoritaria/(1-prop.mayoritaria))) #Índices de la clase mayoritaria para cada subconjunto
  
  dfs[[k]] <- rbind(muertos[id.muertos,],sanos[id.sanos,])
}

alpha_b <- 0.01
b <- ceiling(log(alpha_b)/(log(1-1/np)*np))