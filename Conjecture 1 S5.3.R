#This script is to provide numerical evidence that the equality I(d)= L(d)= D(d) holds only on the 
#region of constant relative potency. Otherwise, we have I(d)>D(d) and L(d)> D(d).

#Here, I(d) refers to the NJR equation (NJR), L(d) refers to Berenbaum's linear isobole method (Berenbaum), and
#D(d) refers to Zaider's equation (Zaider).

#These three helper functions return a vector of Dose-Effec for each method.
#x goes from 0.1 - 10 with 0.1 step, and hence there are 100 x-coordinates.
#For example, to plot the output - type x <- NJR(2, 0, 0, 3, 0.4, 0.6), y <- Berenbaum(2, 0, 0, 3, 0.4, 0.6)
#You can do plot(x) to see how Dose varies with Effect for this x.
#You can do lines((1:100)/10, y[,2]) to graph y on the same plot for comparison.

#NJR returns a vector of Dose-Effect under the solution of NJR.

NJR <- function(alpha1, beta1, alpha2, beta2, r1, r2){
  y.initial <- 1*10^(-30)
  de <- function(t, y, parms)
    list(r1*(sqrt(alpha1^2 + 4*beta1*y)) + r2*(sqrt(alpha2^2 + 4*beta2*y))) #Function form of the Differenial Equation.
  
  library(deSolve)
  times <- seq(from = 0, to = 10, by = 10^(-1))
  out <- ode(y = y.initial, times = times, func = de, parms = NULL, rtol=10^-10, atol=10^-10)
  return(out[-1,]) #Don't return the first value so that there are 100 data points.
}

#Zaider returns a vector of Dose-Effect under Zaider's method.

Zaider <- function(alpha1, beta1, alpha2, beta2, r1, r2){
  out <- matrix(0, 100, 2)
  for (i in 1:100) {
    D <- i/10
    out[i,1] <- D
    out[i,2] <- (alpha1*r1 + alpha2*r2)*D + (beta1*(r1^2) + beta2*(r2^2) + 2*r1*r2*sqrt(beta1*beta2))*D^2 
    #Zaider's Equation
  }
  return(out)
}

#Berenbaum returns a vector of Dose-Effect under Berenbaum's method.

Berenbaum <- function(alpha1, beta1, alpha2, beta2, r1, r2){
  Effect <-function(d,alpha,beta){  
    e = alpha*d+beta*(d^2)
  }
  
  Dose <-function(e,alpha,beta){
    d=(-alpha+sqrt(alpha^2+4*e*beta))/(2*beta) 
  }
  out <- matrix (0, 100, 2)
  for (i in 1:100){
    D = i/10
    
    #Special case where Beta = 0, where the Dose function just becomes Dose = Effect/alpha 
    if (beta2 == 0){
      iso=function(D1){ 
        D1*(1 - (r2*D*alpha2)/Effect(D1,alpha1,beta1)) -r1*D 
      } 
    }
    else {
      iso=function(D1){ 
        D1*(1 - (r2*D)/Dose(Effect(D1,alpha1,beta1),alpha2,beta2)) -r1*D 
      } 
    }
    
    D1givenD = uniroot(iso,c(r1*D+0.001,5*r1*D),extendInt="upX",tol=1.0e-12) 
    out[i,1] = D
    out[i,2]= Effect(D1givenD$root,alpha1,beta1)
  }
  return(out)
}

#NJRZaiderBerenbaum compares all 3 methods in one plot.
#In the output, Column 2 is NJR, Column 3 is Zaider and Column 4 is Berenbaum.

NJRZaiderBerenbaum <- function(alpha1, beta1, alpha2, beta2, r1, r2){
  #This part just repeats code from above to graph all methods on one plot. Except in one case
  #the differences are so small it is hard to see them or (for 2 constant potency cases) are zero.
  y.initial <- 1*10^(-30)
  de <- function(t, y, parms)
    list(r1*(sqrt(alpha1^2 + 4*beta1*y)) + r2*(sqrt(alpha2^2 + 4*beta2*y)))
  library(deSolve)
  times <- seq(from = 0, to = 10, by = 10^(-1)) 
  out <- ode(y = y.initial, times = times, func = de, parms = NULL, rtol=10^-10, atol=10^-10)
  plot(out, main = "Comparison of LQ IDER Methods", xlim = c(0,10), ylim = c(0,200), 
       xlab = "Total Dose D", ylab = "Effect E") 
  curve ((alpha1*r1 + alpha2*r2)*x + (beta1*(r1^2) + beta2*(r2^2) + 2*r1*r2*sqrt(beta1*beta2))*x^2,
         add = TRUE, col = "red")
  
  Effect <-function(d,alpha,beta){  
    e = alpha*d+beta*(d^2)
  }
  Dose <-function(e,alpha,beta){
    d=(-alpha+sqrt(alpha^2+4*e*beta))/(2*beta) 
  }
  Effect_B = rep(0, 100)
  
  for (i in 1:100){
    D = i/10
    if (beta2 == 0){
      iso=function(D1){ 
        D1*(1 - (r2*D*alpha2)/Effect(D1,alpha1,beta1)) -r1*D 
      } 
    }
    else {
      iso=function(D1){ 
        D1*(1 - (r2*D)/Dose(Effect(D1,alpha1,beta1),alpha2,beta2)) -r1*D 
      } 
    }
    D1givenD = uniroot(iso,c(r1*D+0.001,5*r1*D),extendInt="upX",tol=1.0e-12) 
    Effect_B[i]= Effect(D1givenD$root,alpha1,beta1)
  }
  lines((1:100)/10, Effect_B, col = 'green')
  legend(0.5, 170, c("I(d)","D(d)", "L(d)"), lty=c(1,1,1), lwd=c(2.5,2.5,2.5),
         col=c("black","red", "green"))
  
  #This part calls the helper functions to get vectors of Dose-Effect to output along with the graph.
  all3 <- matrix(0,100,4)
  a <- NJR(alpha1, beta1, alpha2, beta2, r1, r2)
  b <- Zaider(alpha1, beta1, alpha2, beta2, r1, r2)
  c <- Berenbaum(alpha1, beta1, alpha2, beta2, r1, r2)
  all3[,1] <- a[,1]
  all3[,2] <- a[,2]
  all3[,3] <- b[,2]
  all3[,4] <- c[,2]
  return(all3)
}

# #a) Checking to see all 3 methods give the same values under constant potency.
# for (i in 1:3){
#   for (j in 1:3) {
#     for (k in 1:3){
#       print(NJRZaiderBerenbaum (i, ((i/j)^2)*k, j, k, 0.4, 0.6))
#     }
#   }
# }
#b) Checking to see I(d)>D(d) and L(d)> D(d) for different values of parameters.
NJRZaiderBerenbaum(2, 0, 0, 3, 0.4, 0.6) #Alpha2 = Beta1 = 0
NJRZaiderBerenbaum(0, 2, 0, 3, 0.4, 0.6) #Alpha1 = Alpha2 = 0 constant potency
NJRZaiderBerenbaum(5, 3, 3, 2, 0.4, 0.6) #Alpha1 > Alpha2, Beta1 > Beta2
NJRZaiderBerenbaum(3, 3, 5, 2, 0.4, 0.6) #Alpha1 < Alpha2, Beta1 > Beta2
NJRZaiderBerenbaum(5, 2, 3, 3, 0.4, 0.6) #Alpha1 > Alpha2, Beta1 < Beta2
NJRZaiderBerenbaum(3, 2, 5, 3, 0.4, 0.6) #Alpha1 < Alpha2, Beta1 < Beta2
NJRZaiderBerenbaum(6, 8, 2, 8/9, 0.4, 0.6) #Constant potency
