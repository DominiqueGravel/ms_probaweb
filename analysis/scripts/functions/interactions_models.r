# Returns the probability of observing the interaction ij 
# independently of co-occurrence and of the environment
# Havens type of interaction
# We have to trick the glm function so that it returns a probability 
# of 1 if the species do interact at least once, and 0 if they never
# interact
L0 <- function(data, selection) {
   with(as.list(data),{

        if(sum(Lij)>0) Lij_trick = numeric(length(Lij)) + 1
          else Lij_trick = numeric(length(Lij))
        model <- glm(Lij_trick ~ 1, family = "binomial")

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
