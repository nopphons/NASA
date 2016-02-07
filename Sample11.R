#This is fully free, open-source software. It comes with no warranty. Used for Siranart et al, submitted to Rad. Res. 2016.
#Compares E_I=I(d) and E_A=S(d) for a mixture of 2 ions or a mixture of 10 ions. Uses targeted effect (TE) individual
#dose-effect relations (IDER) and parameters from Cucinotta. Kim and Chappell (2013).  
rm(list=ls()); library(deSolve)
#next 2 lines:Cucinotta et al. (2013) table 5.4 TE parameter central values; & alph=(IDER slope near zero)
parameter = list(alph0=7.53, alph1=1.261, alph2=0.0037, beta_p = 6.3, lamd0=.25, lamd1=.0051, lamd2=.0034) 
alph=function(LET, param) {with(param,alph0+alph1*LET*exp(-alph2*LET))} #uses component LET (keV/micron)
lamd=function(LET, param) {with(param,lamd0+lamd1*LET*exp(-lamd2*LET))} # used to modulate LQ IDER below

#R function calculates IDER (effect E for each mixture component). proton=1 (H+ or He++) or proton=0 (Z>2)
E_single = function(dose,LET,proton){ #dose in Gy, not equivalent dose in Sv and not flux 
  param = parameter
  beta = 0 #only H+ and He++ are assigned non-zero term quadratic in dose for low doses
  if(proton == 1){ beta= param$beta_p} # next is linear-quadratic (LQ) modulated by exponential decrease
  P_E = (alph(LET,param)*dose + beta*(dose^2))*exp(-lamd(LET,param)*dose) #lamd is ostensibly a "cell killing parameter"
  return(P_E)
}

#function that will be used to find dose when effect E is the independent and dose the dependent variable
solve_d=function(dose,LET,proton,E) {E_single(dose,LET,proton)-E}

E_A <- function(r, doses, LET, proton){#R function that can calculate simple effect additivity E_A=S(d) for mixture
  #r will be a vector of dose proportions; length(r)=nL; sum(r)=100%; dose of jth component=r[j]*(total mixture doses)
  output = rep(0,cols); k0=rep(0,nL)#cols length(doses);LET will be a vector of LET values with length(L)=NL
  for(l in 1:nL){#doses will be a vector of total mixture doses; cols=length(doses)=number of different dose points
    output = output + E_single(r[l]*doses,LET[l],proton[l])#will add contribution from each component to get E_A
  }
  return(output)#E_A=S(d) at various doses 
}

bound_Ed = function(LET,proton){# can calculate, far each IDER, the maximum as dose varies & the dose at the max
  E = c(); D = c(); ##
  for(ll in 1:length(LET)){#LET vector of LETs of mixture components, as above
    p = proton[ll]; l = LET[ll]; param = parameter
    b = param$beta; a = alph(l,param); ld = lamd(l,param)
    if(p == 1){ b= param$beta_p}#now set the derivative=0 and calculate the dose 
    if(b !=0){
      d = (sqrt((a*ld)^2+4*b^2)-a*ld+2*b)/(2*b*ld);#the derivative of the IDER is zero at maximum
    }else{d = 1/(ld)} #
    D = c(D,d)#values of d at the maximum, I'm not sure.
    E = c(E,solve_d(d,l,p,0))#not needed here; was used earlier, retained in case needed again
  }
  upper_E = min(E)#not needed here; was used earlier, retained in case needed again
  max_d = D[which(E==upper_E)]# not needed here; was used earlier, retained in case needed again
  return(list(upper_E,max_d,D))#
}

dose_E=function(E, LET, proton){#not really needed here; was used earlier, retained in case needed again
    result0 = bound_Ed(LET,proton); upper_E = result0[[1]]; d_max = result0[[3]]
  output=c()
  for(i in 1:length(E)){
    e =min(E[i], upper_E-1e-4)#
    dose = (uniroot(solve_d,c(0,.001), upper=d_max, extendInt="upX",tol = 1e-6, LET = LET, 
                    proton=proton, E = e))$root
    output= c(output,dose)
  }
  return(output)
}

first_deriv = function(f,E){#calculate IDER slopes used to get E_I=I(d) slope
  delta = 1e-6
  output = (f(E) - f(E-delta))/delta
  return(output)
}

sum_rE = function(E,parameter){#used to calculate slope for incremental effect additivity default prediction E_I=I(d) 
  output = c()
  for(j in 1:length(E)){
    cumrE = 0
    for(i in 1:nL){
      f = function(x) E_single(x,LET[i],proton[i])
      d = dose_E(E[j],LET[i],proton[i])
      cumrE = cumrE + r[i]*first_deriv(f,d)
    }
    output = c(output,cumrE)
  }
  return(output)
}

solve_ode <- function(t, state, parameters) {#R function used when integrating ordinary differential equation for I(d)
  with(as.list(c(state, parameters)), {
    dE <- sum_rE(E, parameters) 
    list(c(dE))
  })
}

E_I <- function(state,doses,LET, proton){#R function used to actually integrate once user adjusted values are specified
  output = ode(y = state, times = doses, func = solve_ode, parms = parameter)
  result = output[,2]
  return(result)
}

#User adjusted values: dose range; IDER LETs, dose ratios r with sum(r)=100%, and proton (proton=1 if Z<3; else = 0)
end=1.5; by = 0.01; cols = end/by+1; doses <- seq(0, end, by = by); 
LET=25*(1:10); proton = rep(0,10);r=.02*c(9,8,7,6,5,5,4,3,2,1);nL = length(LET)# 10 component HZE mixture
#LET = c(.4,190); r = c(.8,.2); proton = c(1,0) #more realistic 2 component mixture with one component low LET

#Plot one component dose-effect relations (IDER); to increase smoothness decrease "by"
plot(c(0,end),c(0,60),bty="l",ylab ='Effect (%)', xlab = 'dose (Gy') #set size of the plot
for(kk in 1:nL){#
  lines(doses,E_single(doses,LET[kk],proton[kk]),type="l",col = colors()[20*kk+2]) #plot each IDER 
}

E_A_mean = E_A(r, doses, LET ,proton)#find E_A=S(d) with central values of the parameters
state<- c(E = 0);E_I_mean=E_I(state,doses, LET, proton)#solution of differential equation for E_I=I(d); require I(0)=0

plot(c(0,max(doses)),c(0,max(c(E_A_mean,E_I_mean))),bty="l", ylab = 'Effect (%)', xlab = 'dose (Gy)')
lines(doses,E_I_mean,col='red'); lines(doses,E_A_mean,col='blue')