# Returns the probability of observing the interaction ij 
# independently of co-occurrence and of the environment
L0 <- function(data, selection) {
   with(as.list(data),{
        model <- glm(Lij ~ 1, family = "binomial")
        if(selection) {return(step(model,trace=0))}
        else {return(model)}
   })
}

# Returns the probability of observing the interaction ij when i 
# and j are present independently of the environment
L1 <- function(data, selection) {
   with(as.list(data),{
        model <- glm(Lij[Xij==1] ~ 1, family = "binomial")
        if(selection) {return(step(model,trace=0))}
        else {return(model)}
   })
}

# Returns the probability of observing the interaction ij when i 
# and j are present for each level of the environment
L2 <- function(data, selection) {
   with(as.list(data),{
    subLij = Lij[Xij==1] 
    subE = E[Xij==1,]  
    Enames = expression("T","T2","PP","PP2")
    fmla <- as.formula(paste("subLij ~ ", paste(Enames, collapse= "+")))     
    model <- glm(fmla, family = "binomial",data=subE)
    if(selection) {return(step(model,trace=0,direction = "forward"))}
    else {return(model)}
   })
}
