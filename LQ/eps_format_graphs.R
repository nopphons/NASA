rm(list=ls())
library('ggplot2')
load("~/rewriteupfornasapaper/SimulationResults/Jan3/temp.Rdata")

df1$id = "r1 = 0.1"; df2$id="r1 = 0.5"; df3$id= "r1 = 0.9"
df = rbind(df1,df2,df3); df$id = as.factor(df$id)

### Plot 1
library(Cairo)
options(bitmapType="cairo")

cairo_ps(file='~/rewriteupfornasapaper/EbEi_DensityPlot.eps', width=7, height=7)
density_EbEi = ggplot(df, aes(x=Relative_Dif_EbEi, fill = id)) + geom_density(alpha=0.2)
density_EbEi + ggtitle(expression(atop("Density plot of (Eb-Ei)/Ei",
                                       atop(italic("1,500,000 simulations for each value of r1, doses from 0.01 to 1.5"), "")))) + 
  labs(x = "Relative Difference") +
  theme(plot.title = element_text(size = 15, face = "bold", colour = "black", vjust = -1))
dev.off()

cairo_ps(file='~/rewriteupfornasapaper/EbEz_DensityPlot.eps', width=7, height=7)#on Window machine these lines have to be run manually
density_EbEz = ggplot(df, aes(x=Relative_Dif_EbEz, fill = id)) + geom_density(alpha=0.2)
density_EbEz + ggtitle(expression(atop("Density plot of (Eb-Ez)/Ez",
                                       atop(italic("1,500,000 simulations for each value of r1, doses from 0.01 to 1.5"), "")))) + 
  labs(x = "Relative Difference") +
  theme(plot.title = element_text(size = 15, face = "bold", colour = "black", vjust = -1))

dev.off()#Result should be in folder rewriteupfornasapaper

