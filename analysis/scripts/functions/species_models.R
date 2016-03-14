# Returns the probability of observing the species i independently 
# of the environment
C0 <- function(data,selection) {
   with(as.list(data),{
        model = list()
        model$i <- glm(Xi ~ 1, family = "binomial")
        model$j <- glm(Xj ~ 1, family = "binomial")
        if(selection) {return(list(i=step(model$i,trace=0),j=step(model$j,trace=0)))}
        else {return(model)}
   })
}

# Returns the probability of observing the pair ij independently 
# of the environment
C1 <- function(data, selection) {
   with(as.list(data),{
        model = list()    
        model$ij <- glm(Xij ~ 1, family = "binomial")
        if(selection) {return(list(ij = step(model$ij,trace=0)))}
        else {return(list(ij = model$ij))}
   })
}

# Returns the probability of observing the pair ij as a function 
# of the environment
C2 <- function(data, selection) {
   with(as.list(data),{
      Enames = expression("T","T2","PP","PP2")
      fmla <- as.formula(paste("Xij ~ ", paste(Enames, collapse= "+")))     
      model = list()    
      model$ij <- glm(fmla, family = "binomial", data = E)
      if(selection) {return(list(ij = step(model$ij,trace=0)))}
      else {return(list(ij = model$ij))}
   })
}

# Returns the probability of observing the species i as a function 
# of the environment
C3 <- function(data, selection) {
   with(as.list(data),{
 #     Ei = E[match(unique(data$sites[data$IDi == i]),data$sites[data$IDi == i]),]
 #     Ej = E[match(unique(data$sites[data$IDj == j]),data$sites[data$IDj == j]),]
      Enames = expression("T","T2","PP","PP2")
      fmlai <- as.formula(paste("Xi ~ ", paste(Enames, collapse= "+")))  
      fmlaj <- as.formula(paste("Xj ~ ", paste(Enames, collapse= "+"))) 
      model = list()
      model$i <- glm(fmlai, family = "binomial", data=E)
      model$j <- glm(fmlaj, family = "binomial", data=E)
      if(selection) {return(list(i = step(model$i,trace = 0), j= step(model$j,trace=0)))}
      else {return(model)}
   })
}

