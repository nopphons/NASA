rm(list=ls())
library(deSolve)

#NJR returns a vector of Dose-Effect under the solution of NJR.

alden.NJR = function(alpha1, beta1, alpha2, beta2, r1, r2, Dvec = 1:100){
  y.initial = 1*10^(-30)
  de = function(t, y, parms)
    list(r1*(sqrt(alpha1^2 + 4*beta1*y)) + r2*(sqrt(alpha2^2 + 4*beta2*y)))
  times = c(0,Dvec)
  out = ode(y = y.initial, times = times, func = de, parms = NULL, rtol=10^-10, atol=10^-10)
  return(out[-1,])
}

#Zaider returns a vector of Dose-Effect under Zaider's method.
alden.Zaider = function(alpha1, beta1, alpha2, beta2, r1, r2, Dvec = 1:100){
  out = matrix(0, length(Dvec), 2)
  for (i in 1:length(Dvec)) {
    D = Dvec[i]
    d1 = r1*D
    d2 = r2*D
    out[i,1] = D
    out[i,2] = (alpha1*d1 + alpha2*d2) + (beta1*d1^2 + beta2*d2^2) + (2*sqrt(beta1*beta2)*d1*d2)
  }
  return(out)
}

#Berenbaum returns a vector of Dose-Effect under Berenbaum's method.

alden.Berenbaum = function(alpha1, beta1, alpha2, beta2, r1, r2, Dvec = 1:100){
  out <- matrix (0, length(Dvec), 2)
  for (i in 1:length(Dvec)){
    D = Dvec[i]
    d1 = r1*D
    d2 = r2*D
    if (beta2 == 0){
      iso = function(D1){
        D1 - d1 - alpha2*d2*D1/(alpha1*D1+beta1*D1^2)
      }
    }else{
      iso=function(D1){
        D1 - d1 - 2*beta2*d2*D1/(-alpha2+sqrt(alpha2^2+4*beta2*(alpha1*D1+beta1*D1^2)))
      }
    }
    D1e=uniroot(iso,c(d1+.0001,5*d1),extendInt="upX",tol=1.0e-12)
    out[i,2]= alpha1*D1e$root+beta1*D1e$root^2
  }
  out[,1] = Dvec
  return(out)
}

# alden.Berenbaum(2,1,3,4,0.2,0.8, 1:100)
# alden.Zaider(2,1,3,4,0.2,0.8, 1:100)
# alden.NJR(2,1,3,4,0.2,0.8, 1:100)

combined = function(alpha1,beta1,alpha2,beta2,r1,r2, Dvec = 1:100, tol=10^-7){
  E_B_mat = alden.Berenbaum(alpha1,beta1,alpha2,beta2,r1,r2, Dvec)
  E_I_mat = alden.NJR(alpha1,beta1,alpha2,beta2,r1,r2, Dvec)
  E_Z_mat = alden.Zaider(alpha1,beta1,alpha2,beta2,r1,r2, Dvec)
  
  dif1 = (E_B_mat[,2] - E_I_mat[,2])/E_I_mat[,2]
  dif2 = (E_I_mat[,2] - E_Z_mat[,2])/E_Z_mat[,2]
  dif3 = (E_B_mat[,2] - E_Z_mat[,2])/E_Z_mat[,2]
  
  Eb_geq_Ei = dif1 >= -tol
  Ei_geq_Ez = dif2 >= -tol
  Eb_geq_Ez = dif3 >= -tol
  
  Eb_eq_Ei = abs(dif1)<=tol
  Ei_eq_Ez = abs(dif2)<=tol
  Eb_eq_Ez = abs(dif3)<=tol
  
  df = data.frame(alpha1 = alpha1, beta1=beta1, alpha2 = alpha2, beta2 = beta2,
                  Dose1 = r1*Dvec, Dose2 = r2*Dvec,
                  E_B = E_B_mat[,2], E_I = E_I_mat[,2], E_Z = E_Z_mat[,2],
                  Relative_Dif_EbEi = dif1 , Relative_Dif_EiEz = dif2, Relative_Dif_EbEz = dif3,
                  Eb_geq_Ei = Eb_geq_Ei, Ei_geq_Ez = Ei_geq_Ez, Eb_geq_Ez = Eb_geq_Ez,
                  Eb_eq_Ei = Eb_eq_Ei, Ei_eq_Ez=Ei_eq_Ez, Eb_eq_Ez=Eb_eq_Ez,
                  ConstRelPot = (beta2/beta1) * (alpha1/alpha2)^2)
  
  return(df)
}

# alpha1: mean=9, var=16
# beta1: mean=6, var=9
# alpha2: mean=110, var=1600

# rate = mean/var
# shape = (mean^2)/var
mean.alpha1 = 9; var.alpha1 = 16;
rate.alpha1 = mean.alpha1/var.alpha1; shape.alpha1 = (mean.alpha1^2)/var.alpha1

mean.beta1 = 6; var.beta1 = 9
rate.beta1 = mean.beta1/var.beta1; shape.beta1 = (mean.beta1^2)/var.beta1

mean.alpha2 = 110; var.alpha2 = 1600
rate.alpha2 = mean.alpha2/var.alpha2; shape.alpha2 = (mean.alpha2^2)/var.alpha2

N=100
set.seed(1234567)
sim.alpha1 = rgamma(N, shape = shape.alpha1, rate = rate.alpha1)
set.seed(1234567)
sim.beta1 = rgamma(N, shape = shape.beta1, rate = rate.beta1)
set.seed(1234567)
sim.alpha2 = rgamma(N, shape = shape.alpha2, rate = rate.alpha2)

##### 3 Cases of interest: r1 = 0.1, 0.5, 0.9

doses = seq(from = 0.01, to = 1.50, length=150)

### r1 = 0.1
r1 = 0.1; r2 = 1 - r1
df1 = combined(sim.alpha1[1], sim.beta1[1], sim.alpha2[1], 0,
               r1, r2, Dvec = doses, tol=10^-7)
for (i in 2:N){
  temp.df = combined(sim.alpha1[i], sim.beta1[i], sim.alpha2[i], 0,
                     r1, r2, Dvec = doses, tol=10^-7)
  df1 = rbind(df1,temp.df)
  if (i%%100==0) print(i)
}

### r1 = 0.5
r1 = 0.5; r2 = 1 - r1
df2 = combined(sim.alpha1[1], sim.beta1[1], sim.alpha2[1], 0,
               r1, r2, Dvec = doses, tol=10^-7)
for (i in 2:N){
  temp.df = combined(sim.alpha1[i], sim.beta1[i], sim.alpha2[i], 0,
                     r1, r2, Dvec = doses, tol=10^-7)
  df2 = rbind(df2,temp.df)
  if (i%%100==0) print(i)
}

### r1 = 0.9
r1 = 0.9; r2 = 1 - r1
df3 = combined(sim.alpha1[1], sim.beta1[1], sim.alpha2[1], 0,
               r1, r2, Dvec = doses, tol=10^-7)
for (i in 2:N){
  temp.df = combined(sim.alpha1[i], sim.beta1[i], sim.alpha2[i], 0,
                     r1, r2, Dvec = doses, tol=10^-7)
  df3 = rbind(df3,temp.df)
  if (i%%100==0) print(i)
}

write.csv(df1, "~/rewriteupfornasapaper/SimulationResults/Jan3/r1_1.csv")#my R is in folder "Documents" and I had to create a file
#create a folder path "rewriteupfornasapaper/SimulationResults/Jan3/" in "Documents"
write.csv(df2, "~/rewriteupfornasapaper/SimulationResults/Jan3/r1_5.csv")
write.csv(df3, "~/rewriteupfornasapaper/SimulationResults/Jan3/r1_9.csv")

save.image("~/rewriteupfornasapaper/SimulationResults/Jan3/temp.Rdata")

