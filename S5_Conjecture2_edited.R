### This script reproduces the density plots in S5 of the supplement. On some platforms it may be necessary
### to source the whole file first and then use Ctrl-Enter to run the plots by hand afterwards
### User can produce density plots for other values of the parameters by changing the parameters in the
###section "Simulations for beta2 = 0". The plots will be OK, but the plot labels will not track the changes

rm(list=ls())
### Uncomment the following line if you never installed the "deSolve" packages
# install.packages("deSolve")
library(deSolve)

################################## Functions for Predicted Dose Effects ####################################

# Linear-Quadratic IDER: (Effect i) = alpha_i*d_i + beta_i*d_i^2,
# d_i = r_i*D, where D is the total dose.

#  NJR=I(d) is a function that takes as inputs, parameter values for two Linear-Quadratic (LQ) IDERS
# as well as a vector of total doses, and outputs a 2-column dose-effect matrix.
# The first column of the matrix contains the input doses, and the second column 
# contains predicted effects under the hypothesis of incremental additivity.

NJR = function(alpha1, beta1, alpha2, beta2, r1, r2, Dvec = 1:100){
  y.initial = 1*10^(-30)
  de = function(t, y, parms)
    list(r1*(sqrt(alpha1^2 + 4*beta1*y)) + r2*(sqrt(alpha2^2 + 4*beta2*y)))
  times = c(0,Dvec)
  out = ode(y = y.initial, times = times, func = de, parms = NULL, rtol=10^-10, atol=10^-10)
  return(out[-1,])
}

# Zaider=D(d) is a function that takes as inputs, parameter values for two Linear-Quadratic (LQ) IDERS
# as well as a vector of doses, and outputs a 2-column dose-effect matrix.
# The first column of the matrix contains the input doses, and the second column 
# contains predicted effects under Zaider's hypothesis.

Zaider = function(alpha1, beta1, alpha2, beta2, r1, r2, Dvec = 1:100){
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

# Berenbaum=L(d) is a function that takes as inputs, parameter values for two Linear-Quadratic (LQ) IDERS
# as well as a vector of total doses, and outputs a 2-column dose-effect matrix.
# The first column of the matrix contains the input doses, and the second column 
# contains predicted effects under Berenbaum's linear isobole hypothesis.

Berenbaum = function(alpha1, beta1, alpha2, beta2, r1, r2, Dvec = 1:100){
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



# The following function takes as inputs, parameter values for two Linear-Quadratic (LQ) IDERS
# as well as a vector of total doses, and returns the 3 predictions as well as the 
# relative differences between the predicted effects.

combined = function(alpha1,beta1,alpha2,beta2,r1,r2, Dvec = 1:100){
  E_B_mat = Berenbaum(alpha1,beta1,alpha2,beta2,r1,r2, Dvec)
  E_I_mat = NJR(alpha1,beta1,alpha2,beta2,r1,r2, Dvec)
  E_Z_mat = Zaider(alpha1,beta1,alpha2,beta2,r1,r2, Dvec)
  
  dif1 = (E_B_mat[,2] - E_I_mat[,2])/E_I_mat[,2]
  dif2 = (E_I_mat[,2] - E_Z_mat[,2])/E_Z_mat[,2]
  dif3 = (E_B_mat[,2] - E_Z_mat[,2])/E_Z_mat[,2]
  
  df = data.frame(alpha1 = alpha1, beta1=beta1, alpha2 = alpha2, beta2 = beta2,
                  Dose1 = r1*Dvec, Dose2 = r2*Dvec,
                  E_B = E_B_mat[,2], E_I = E_I_mat[,2], E_Z = E_Z_mat[,2],
                  Relative_Dif_EbEi = dif1 , Relative_Dif_EiEz = dif2, Relative_Dif_EbEz = dif3)
  return(df)
}


# Next, we run the simulations described in S5 of the supplement, using the functions defined above.

############################# Simulations for beta2 = 0 ############################

###### Draw from distributions that yield parameter values of interest

### As described in detail in S5 of the supplement, we draw parameter values from gamma distributions,
### with means and variances of interest.

### Parameter means and variances of interest:
# alpha1: mean=9, var=16  dimension=Gy^-1
# beta1: mean=6, var=9  dimension=Gy^-2
# alpha2: mean=110, var=1600

### Relationship between mean and variance of a gamma distribution, and its rate and shape parameters:
# rate = mean/var
# shape = (mean^2)/var

mean.alpha1 = 9; var.alpha1 = 16;
rate.alpha1 = mean.alpha1/var.alpha1; shape.alpha1 = (mean.alpha1^2)/var.alpha1

mean.beta1 = 6; var.beta1 = 9
rate.beta1 = mean.beta1/var.beta1; shape.beta1 = (mean.beta1^2)/var.beta1

mean.alpha2 = 110; var.alpha2 = 1600
rate.alpha2 = mean.alpha2/var.alpha2; shape.alpha2 = (mean.alpha2^2)/var.alpha2

### Here, we obtain 10,000 random combinations of parameter values. Use 100 while debugging
N=10000
set.seed(1234567)
sim.alpha1 = rgamma(N, shape = shape.alpha1, rate = rate.alpha1)
set.seed(1234567)
sim.beta1 = rgamma(N, shape = shape.beta1, rate = rate.beta1)
set.seed(1234567)
sim.alpha2 = rgamma(N, shape = shape.alpha2, rate = rate.alpha2)

doses = seq(from = 0.01, to = 1.5, length=150)# doses are in Gy, not cGy, not Sv

### As described in S5 of the supplement, we consider 3 values of r1: 0.1, 0.5, 0.9.

### r1 = 0.1
r1 = 0.1; r2 = 1 - r1
df1 = combined(sim.alpha1[1], sim.beta1[1], sim.alpha2[1], 0,
               r1, r2, Dvec = doses)
for (i in 2:N){
  temp.df = combined(sim.alpha1[i], sim.beta1[i], sim.alpha2[i], 0,
                     r1, r2, Dvec = doses)
  df1 = rbind(df1,temp.df)
  if (i%%100==0) print(i)
}

### r1 = 0.5
r1 = 0.5; r2 = 1 - r1
df2 = combined(sim.alpha1[1], sim.beta1[1], sim.alpha2[1], 0,
               r1, r2, Dvec = doses)
for (i in 2:N){
  temp.df = combined(sim.alpha1[i], sim.beta1[i], sim.alpha2[i], 0,
                     r1, r2, Dvec = doses)
  df2 = rbind(df2,temp.df)
  if (i%%100==0) print(i)
}

### r1 = 0.9
r1 = 0.9; r2 = 1 - r1
df3 = combined(sim.alpha1[1], sim.beta1[1], sim.alpha2[1], 0,
               r1, r2, Dvec = doses)
for (i in 2:N){
  temp.df = combined(sim.alpha1[i], sim.beta1[i], sim.alpha2[i], 0,
                     r1, r2, Dvec = doses)
  df3 = rbind(df3,temp.df)
  if (i%%100==0) print(i)
}

### Uncomment these lines and substitute in your desired output file destination 
### if you want to save the simulation results as csv files.

#write.csv(df1, "~/Desktop/rewriteupfornasapaper/SimulationResults/r1_1.csv")
#write.csv(df2, "~/Desktop/rewriteupfornasapaper/SimulationResults/r1_5.csv")
#write.csv(df3, "~/Desktop/rewriteupfornasapaper/SimulationResults/r1_9.csv")


######################### Produce the Density Plots in S5 of the Supplement ##########################

### Uncomment the following lines if you do not have the packages installed in R
# install.packages('ggplot2')
library('ggplot2')

df1$id = "r1 = 0.1"; df2$id="r1 = 0.5"; df3$id= "r1 = 0.9"
df = rbind(df1,df2,df3); df$id = as.factor(df$id)

### Plot 1: Relative difference between predictions from linear isobole, and incremental additivity. If needed
### rerun lines 189-196 again after sourcing the whole file.
density_EbEi = ggplot(df, aes(x=Relative_Dif_EbEi, fill = id)) + geom_density(alpha=0.2)
# The user may change the plot title, x-labels and file location of the saved plot
# by modifying the following lines of code accordingly.
density_EbEi + ggtitle(expression(atop("Density plot of (Eb-Ei)/Ei",
                                       atop(italic("1,500,000 simulations for each value of r1, doses from 0.01 to 1.5"), "")))) + 
  labs(x = "Relative Difference") +
  theme(plot.title = element_text(size = 15, face = "bold", colour = "black", vjust = -1))
ggsave(file='EbEi_DensityPlot.pdf', width=7, height=7)

### Plot 2: Relative difference between predictions from linear isobole, and Zaider; see comments on plot 1
density_EbEz = ggplot(df, aes(x=Relative_Dif_EbEz, fill = id)) + geom_density(alpha=0.2)
# The user may change the plot title, x-labels and file location of the saved plot
# by modifying the following lines of code accordingly.
density_EbEz + ggtitle(expression(atop("Density plot of (Eb-Ez)/Ez",
                                       atop(italic("1,500,000 simulations for each value of r1, doses from 0.01 to 1.5"), "")))) + 
  labs(x = "Relative Difference") +
  theme(plot.title = element_text(size = 15, face = "bold", colour = "black", vjust = -1))
ggsave(file='EbEz_DensityPlot.pdf', width=7, height=7)

# 
# # # The code below (commented out) allows the user to create the same graphs in .eps format.
# # # However, this code was tested on a MacOSX system, and is not guaranteed to work on other platforms.
# # 
# ### Uncomment the following lines if you do not have the packages installed in R
# # install.packages('ggplot2')
# # install.packages('Cairo')
# library('ggplot2')
# library(Cairo)
# 
# df1$id = "r1 = 0.1"; df2$id="r1 = 0.5"; df3$id= "r1 = 0.9"
# df = rbind(df1,df2,df3); df$id = as.factor(df$id)
# 
# options(bitmapType="cairo")
# ### Plot 1: Relative difference between predictions from linear isobole, and incremental additivity
# #cairo_ps(file='~/Desktop/rewriteupfornasapaper/EbEi_DensityPlot.eps', width=7, height=7)
# cairo_ps(file='EbEi_DensityPlot.eps', width=7, height=7)# more stable though less convenient location
# density_EbEi = ggplot(df, aes(x=Relative_Dif_EbEi, fill = id)) + geom_density(alpha=0.2)
# # The user may change the plot title, x-labels and file location of the saved plot
# # by modifying the following lines of code accordingly.
# density_EbEi + ggtitle(expression(atop("Density plot of (Eb-Ei)/Ei",
#                                        atop(italic("1,500,000 simulations for each value of r1, doses from 0.01 to 1.5"), "")))) + 
#   labs(x = "Relative Difference") +
#   theme(plot.title = element_text(size = 15, face = "bold", colour = "black", vjust = -1))
# dev.off()
# 
# ### Plot 2: Relative difference between predictions from linear isobole, and Zaider
# cairo_ps(file='~/Desktop/rewriteupfornasapaper/EbEz_DensityPlot.eps', width=7, height=7)
# density_EbEz = ggplot(df, aes(x=Relative_Dif_EbEz, fill = id)) + geom_density(alpha=0.2)
# # The user may change the plot title, x-labels and file location of the saved plot
# # by modifying the following lines of code accordingly.
# density_EbEz + ggtitle(expression(atop("Density plot of (Eb-Ez)/Ez",
#                                        atop(italic("1,500,000 simulations for each value of r1, doses from 0.01 to 1.5"), "")))) + 
#   labs(x = "Relative Difference") +
#   theme(plot.title = element_text(size = 15, face = "bold", colour = "black", vjust = -1))
# 
# dev.off()