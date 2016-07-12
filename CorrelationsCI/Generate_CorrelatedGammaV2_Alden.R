rm(list=ls())
require('MASS')
require('clusterGeneration')

N=1000
TE = list(alph0=7.65, alph1=1.25, alph2=0.0038, beta = 0, beta_p = 6.02, lamd0=.243, lamd1=.006, lamd2=.0043,
          e_alph0=3.94, e_alph1=0.14, e_alph2=0.0004, e_beta = 0, e_beta_p =3.51, e_lamd0=.07, e_lamd1=.0036, e_lamd2=.0027)
my.seed = 1234567

myCovMat = matrix(c(1, -0.6, 0.1, 0.8, -0.3, -0.2,
                    -0.6, 1, 0.6, -0.8, 0.65, 0.2,
                    0.1, 0.6, 1, -0.3, 0.63, 0.03,
                    0.8, -0.8, -0.3, 1, -0.75, -0.57,
                    -0.3, -0.65, 0.63, -0.75, 1, 0.63,
                    -0.2, 0.2, 0.03, -0.57, 0.63, 1), ncol = 6)


set.seed(my.seed)
gauss.sim = mvrnorm(n=N, mu = rep(0,6), Sigma = myCovMat)

alph0 = pnorm(gauss.sim[,1])
alph1 = pnorm(gauss.sim[,2])
alph2 = pnorm(gauss.sim[,3])
lamd0 = pnorm(gauss.sim[,4])
lamd1 = pnorm(gauss.sim[,5])
lamd2 = pnorm(gauss.sim[,6])

set.seed(my.seed)
alph0 = qgamma(alph0, (TE$alph0/TE$e_alph0)^2, TE$alph0/(TE$e_alph0)^2)
alph1 = qgamma(alph1, (TE$alph1/TE$e_alph1)^2, TE$alph1/(TE$e_alph1)^2)
alph2 = qgamma(alph2, (TE$alph2/TE$e_alph2)^2, TE$alph2/(TE$e_alph2)^2)
lamd0 = qgamma(lamd0, (TE$lamd0/TE$e_lamd0)^2, TE$lamd0/(TE$e_lamd0)^2)
lamd1 = qgamma(lamd1, (TE$lamd1/TE$e_lamd1)^2, TE$lamd1/(TE$e_lamd1)^2)
lamd2 = qgamma(lamd2, (TE$lamd2/TE$e_lamd2)^2, TE$lamd2/(TE$e_lamd2)^2)

cor(alph0, alph1)
cor(alph1, alph2)
cor(lamd0, lamd1)
cor(lamd1, lamd2)

L = c(76,174)
cor(alph0+alph1*L[1],alph2)
cor(alph0+alph1*L[2],alph2)
cor(lamd0+lamd1*L[1],lamd2)
cor(lamd0+lamd1*L[2],lamd2)

alph = alph0 + alph1*exp(-alph2)
lamd = lamd0 + lamd1*exp(-lamd2)

cor(alph, lamd)

alph0.vec = alph0
alph1.vec = alph1
alph2.vec = alph2
lamd0.vec = lamd0
lamd1.vec = lamd1
lamd2.vec = lamd2

set.seed(my.seed)
beta_p.vec=rgamma(N,(TE$beta_p/TE$e_beta_p)^2,TE$beta_p/(TE$e_beta_p)^2)

save(alph0.vec, alph1.vec, alph2.vec, lamd0.vec, lamd1.vec, lamd2.vec, beta_p.vec,
     file = "~/Desktop/RadiationPaperReview/gamma_realisticV2.RData")