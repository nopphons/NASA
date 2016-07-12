	##This is free, open-source software under GNU GPLv3. It comes with no warranty. Concerns radiogenic mouse HG gland tumorigenesis.
	##Uses individual dose-effect relations (IDER) & parameters from Cucinotta1 et al. (2013). Space Radiation Cancer Risk Projections and
	##Uncertainties - 2012. NASA Center for AeroSpace Information, Hanover, MD; http://ston.jsc.nasa.gov/collections/TRS
	##User adjustable parameters are those near line 160 and also limits on x and y axes for some of the graphs#147->153RKS
	##Gives main calculations for Siranart, Cheng, Handa, Sachs, manuscript submitted to Rad. Res.
	##Main outputs: predictions E_A=S(d) and E_I=I(d) with 95% CI, with graphs, for effect of mixtures of ions (fully ionized atomic nuclei)
	library(ggplot2); library(deSolve)
	rm(list=ls())
	load("~/Desktop/RadiationPaperReview/gamma_realisticV2.RData")
	outputfile = "~/Desktop/RadiationPaperReview/fig8B_RealisticCorV2.eps"
	###FIXED PARAMETERS###next line has Cucinotta et al. (2013) table 5.4 IDER parameters for ions. TE=targeted effects, NTE = non-TE.
# 	parameter = list(TE = list(alph0=7.53, alph1=1.261, alph2=0.0037, beta = 0, beta_p = 6.3, lamd0=.25, lamd1=.0051, lamd2=.0034,
# 	e_alph0=3.96, e_alph1=0.213, e_alph2=0.00058, e_beta = 0, e_beta_p = 3.41, e_lamd0=.065, e_lamd1=.0029, e_lamd2=.0027),
# 	NTE=list(alph0=10.02,alph1=0.679,alph2=0.0033,beta=0,beta_p=5.08,lamd0=.231,lamd1=.0033,lamd2=.005,k1=0.12,k2=0.0053,
# 	e_alph0=2.07,e_alph1=0.187,e_alph2=0.0006,e_beta=0,e_beta_p=3,e_lamd0=.016,e_lamd1=.0042,e_lamd2=.0064,e_k1=0.06,e_k2=0.002))
#   
  ##next lines: values in Blakely et al. "Harderian Gland Tumorigenesis: Low-Dose and LET Response." Radiat Res. Epub 04/19/2016.  
	parameter = list(TE = list(alph0=7.65, alph1=1.25, alph2=0.0038, beta = 0, beta_p = 6.02, lamd0=.243, lamd1=.006, lamd2=.0043,
	e_alph0=3.94, e_alph1=0.14, e_alph2=0.0004, e_beta = 0, e_beta_p =3.51, e_lamd0=.07, e_lamd1=.0036, e_lamd2=.0027),
  NTE=list(alph0=10.05,alph1=0.90,alph2=0.0039,beta=0,beta_p=4.61,lamd0=.219,lamd1=.0047,lamd2=.0051,k1=0.048,k2=0.0028,#for NTE1
	e_alph0=3.56,e_alph1=0.21,e_alph2=0.0009,e_beta=0,e_beta_p=3.33,e_lamd0=.078,e_lamd1=.0059,e_lamd2=.0059,e_k1=0.023,e_k2=0.0019))
  
  alph=function(L, param) {with(param,alph0+alph1*L*exp(-alph2*L))} #slope near zero of the IDER; depends on LET L (keV/micron)
	lamd=function(L, param) {with(param,lamd0+lamd1*L*exp(-lamd2*L))}
	#lamd is controversially interpreted as giving a cell killing modulation of LQ IDER. It leads to IDER which have maxima.
	kappa=function(L,param){with(param, (k1*L)*exp(-k2*L))}#kappa interpreted as NTE, important near d=0; approximated as initial value.
	
	###R FUNCTIONS### Used later after user-defined parameters have been set (near line 147)
	
	#Next R function calculates an IDER "targeted part" T (Cucinotta et al. 2013); uses proton=0 for ion charge Z>2; proton=1 for H, He ions
	T_single = function(dose, L, Eff, proton, param = parameter[[Eff]]){ #if Eff=: 1=TE=targeted effect model; =2=NTE=non-targeted (bystander) effect.
	beta = param$beta #initialize to be beta=0 with Z>2, i.e near d=0 TE effect is linear no threshhold (LNT) for Z>2.
	if(proton == 1){ beta= param$beta_p} #ions with Z<3 are approximated as LQ (linear quadratic) near d=0.
	P_E = (alph(L,param)*dose + beta*(dose^2))*exp(-lamd(L,param)*dose) #LQ or LNT modulated by exponentially decreasing factor.
	return(P_E)
	}
	
	E_single = function(dose, L, Eff, proton, param = parameter[[Eff]]){#if Eff=1, E_single=T; #if Eff=2, E_single=T+NT, where NT=kappa*exp(-lamd*dose).  
	P_E = T_single(dose,L,Eff,proton)
	if (Eff==2){
	P_E = P_E + kappa(L,param)*exp(-lamd(L,param)*dose)#P_E here is sum of targeted and non-targeted effects
	}
	return(P_E)
	}
	
	#Each T_single has a maximum, often around dose=1 Gy. Next is a function to find the maxima and the doses where they occur.
	bound_Ed = function(L,Eff,proton, param){##
	E = c(); D = c(); ##
	for(ll in 1:length(L)){#when this function is used, L will be a vector of the LETs of the individual components of a mixed-ion beam
	p = proton[ll]; l = L[ll];
	b = param[[ll]]$beta; a = alph(l,param[[ll]]); ld = lamd(l,param[[ll]])
	if(p == 1){ b= param[[ll]]$beta_p}
	if(b !=0){
	d = (sqrt((a*ld)^2+4*b^2)-a*ld+2*b)/(2*b*ld);#the derivative of the IDER is zero at maximum, occuring at this dose.
	}else{d = 1/(ld)}
	D = c(D,d)
	E = c(E,solve_d(d,l,Eff,p,0, param[[ll]]))
	}
	upper_E = min(E)#smallest of the maxima in a mixed beam
	return(list(upper_E,D))
	}
	
	#R function which uniroot can use to get dose from L,Eff,proton, and T_single;  treats E=T_single as independent and
	#d as dependent variable instead of vice versa as in RBE calculations.
	solve_d=function(dose,L,Eff,proton,E, param = parameter[[Eff]]) {
	eq = T_single(dose,L,Eff,proton,param)-E
	return(eq)
	}
	
	gamma_param = function(Eff,iter){#When calculating 95% CI for mixture, IDER parameters are taken to have gamma distributions (and stay >0)
	  param = parameter[[Eff]]
	  param_error=with(param,list(alph0=alph0.vec[iter], alph1=alph1.vec[iter],
	                              alph2=alph2.vec[iter], beta=0, beta_p=beta_p.vec[iter],
	                              lamd0=lamd0.vec[iter], lamd1=lamd1.vec[iter],
	                              lamd2=lamd2.vec[iter]))
	  if(Eff ==2){
	    param_error1 = with(param,list(k1 = rgamma(1,(k1/e_k1)^2,k1/(e_k1)^2), k2 = rgamma(1,(k2/e_k2)^2,k2/(e_k2)^2)))
	    param_error = c(param_error,param_error1)
	  }
	  return(param_error)
	}
	gen_param = function(Eff,N,random,iter){
	  param = list()
	  p = gamma_param(Eff,iter)
	  for(ii in 1:N){
	    param[[ii]] = p #random parameters
	    if (random ==0) param[[ii]] = parameter[[Eff]]
	  }
	  return(param)
	}
	
	#R function for d = D_j(E), where D_j will be the compositional inverse function to the jth IDER for the jth ion in a mixture.
	dose_E=function(E, L, Eff,proton, param){#sometimes used with E a vector of effects, not just a single effect
	result0 = bound_Ed(L,Eff,proton, param); upper_E = result0[[1]]; d_max = result0[[2]]
	output=c()
	e =min(E, upper_E)
	dose = (uniroot(solve_d,c(0,.001), upper=d_max, extendInt="upX",tol = 1e-6, L = L,
	Eff = Eff, proton=proton, E = e, param = param[[1]]))$root
	output= c(output,dose)
	return(output)
	}
	
	#Next is function to get simple effect additivity default prediction E_A=S(d) for Eff=1 or E_A=S*(d) for Eff=2
	E_A <- function(doses, L, Eff, proton, param = parameter[[Eff]]){#use T (i.e subtract out kappa term), manipulate, later add biggest back in.
	output = rep(0,cols); k0=rep(0,nL)
	for(l in 1:nL){
	output = output + T_single(r[l]*doses,L[l],Eff,proton[l],param[[l]])
	}
	if (Eff==2){
	maxx=which.max(kappa(L,param[[l]]))
	output = output + kappa(L[maxx],param[[l]])*exp(-lamd(L[maxx],param[[l]])*doses)
	}
	return(output)
	}
	
	E_Anm <- function(doses, L, Eff, proton, param = parameter[[Eff]]){#whole function is RKS
	  output = rep(0,cols); k0=rep(0,nL)
	  for(l in 1:nL){
	    output = output + E_single(r[l]*doses,L[l],Eff,proton[l],param[[l]])
	  }
	  return(output)
	}	
  
	#R function to get first derivative of function f evaluated at E, using effect rather than dose as independent variable
	first_deriv = function(f,E){
	delta = 1e-6
	output = (f(E) - f(E-delta))/delta
	return(output)
	}
	
	sum_rE = function(E, param){#used to calculate slopes for integrating E_I=I(d); apply to each T_single; then add in kappa term later
	output = c()
	for(j in 1:length(E)){
	cumrE = 0
	for(i in 1:nL){
	f = function(x) T_single(x,L[i],Eff,proton[i],param[[i]])
	temp = list(); temp[[1]] = param[[i]]
	d = dose_E(E[j],L[i],Eff,proton[i],param=temp)
	cumrE = cumrE + r[i]*first_deriv(f,d)
	}
	output = c(output,cumrE)
	}
	return(output)
	}
	
	solve_ode <- function(t, state, parameters) {#R function used when integrating the ordinary differential eq. for I(d)
	with(as.list(c(state, parameters)), {
	dE <- sum_rE(E, parameters)
	list(c(dE))
	})
	}
	
	E_I <- function(state,doses,L, Eff, proton, param){#R function used to actually integrate
	output = ode(y = state, times = doses, func = solve_ode, parms = param)
	result = output[,2]
	if(Eff ==2){
	maxx=which.max(kappa(L,param[[1]]))
	result = result + kappa(L[maxx],param[[1]])*exp(-lamd(L[maxx],param[[1]])*doses)}
	return(result)
	}
	
	#####EXAMPLES OF USER-DEFINED PARAMETERS FOR MIXTURE CALCULATIONS#####
	#r=mixture proportions (must sum to 1); L=LETs; proton=1 for H, He, and proton =0 for heavier ions (0 means LQ beta is 0). Examples follow
	#r=1; L=25;proton=0;Eff=2;N=50#Example of 1-ion "mixture"
	#r=.01*c(50,20,rep(5,6));L=c(.4,1.4,25,21,76,107,174,464);proton=c(1,1,rep(0,6)); Eff=1; N=1000
	#r=c(.5,.5);L=c(76,175);proton=c(0,0);Eff=2;N=50 #Ellie 60-40; next is a mixture of 10 different heavy ions.
	#r=c(.02*c(9,8,7,6,5,5,4,3,2,1)); L=c(25*(1:10)); proton =rep(0,10); Eff=2; N = 50 #use N >= 50 to get 95% CI
	#	r=c(.7,.006*c(9,8,7,6,5,5,4,3,2,1)); L=c(.4,25*(1:10)); proton =c(1,rep(0,10)); Eff=1; N = 50# another example RKS
	#r=c(.5,.5); L=c(50,200); proton = c(0,0); Eff=1; N = 100 #Use even bigger N  (e.g. N=10000)for more accurate CI
	#r=c(.4,.4,.2); L=c(50,50,200); proton = c(0,0,0); Eff=2; N = 100 #test obedience to sham mixture principle
	# r=c(.9,.1);L=c(.4,190);proton=c(1,0);Eff=1;N=50# row A in new table 4
	# r=c(.65,.35);L=c(.4,190);proton=c(1,0);Eff=1;N=50# row B in table 4
	# r=c(.80,.20);L=c(.4,60);proton=c(1,0);Eff=1;N=50# row C in table 4
	# r=c(.8,.12,.08); L=c(.4,25,60); proton = c(1,0,0); Eff=1; N = 50#row D
	#r=c(.76,.18,.06); L=c(.4,60,190); proton = c(1,0,0); Eff=1; N = 50#row E
	# r=c(.70,.30);L=c(60,190);proton=c(0,0);Eff=1;N=50#row G
	# r=c(.6,.3,.1); L=c(25,110,190); proton = c(0,0,0); Eff=1; N = 50 #row H 
	r=c(.5,.5); L=c(76,174); proton=c(0,0); Eff=1; N=1000
	nL = length(L); end =1; by = 0.01; cols = end/by+1; doses <- seq(0, end, by = by); #to increase smoothness decrease "by"
	#####END OF MAIN USER DEFINED PARAMETERS######
	
	############E_A with central values of the parameters##############
	E_A_mean = E_A(doses, L ,Eff,proton, gen_param(Eff,nL,0,iter=1))
  E_Anm=E_Anm(doses, L ,Eff,proton, gen_param(Eff,nL,0,iter=1))#RKS
	############E_I with central values of the parameters###############
	state <- c(E = 0);
	E_I_mean = E_I(state,doses, L, Eff, proton,gen_param(Eff,nL,0,iter=1))
  
	##########95% confidence intervals;  Most CPU intensive part of the calculation; comment out if only want central parameter value results.
	E_A_mat = matrix(0, N, cols); ##now do Monte Carlo for error structure; at first without any kappa term even if Eff=2
	set.seed(1234567)
	for(i in 1:N){
	  E_A_mat[i,] = E_A(doses, L, Eff,proton, gen_param(Eff,nL,1,iter=i))
	}
	E_A_sort = apply(E_A_mat,2,sort,decreasing=F);
	E_A_low = E_A_sort[floor(N*0.025),]; E_A_high = E_A_sort[ceiling(N*0.975),]; print(E_A_high[cols])
	
	E_I_mat = matrix(0, N, cols);
	for(i in 1:N){
	  print(i)
	  param = gen_param(Eff,nL,1,iter=i)
	  E_I_mat[i,] = E_I(state,doses,L,Eff,proton, param)
	}
	E_I_sort = apply(E_I_mat,2,sort,decreasing=F);
	E_I_low = E_I_sort[floor(N*0.025),]; E_I_high = E_I_sort[ceiling(N*0.975),];
	
	####################GRAPHS####################Comment out unneeded ones
	#IDER Plot. Allows adjustment of upper dose limit to less than "end" under user adjsted parameters.
	xstart = 0; xend = 1.5; ystart= 0; yend = 60
	plot(c(xstart,xend),c(ystart,yend),type='n',bty="l", ylab = 'E', xlab = 'd')
	for(l in 1:length(L)){
	lines(doses,E_single(doses,L[l],Eff,proton[l]),type="l", col = colors()[20*l+2]) #automatic choice of colors 
	}
	
	##Use next plot only for Eff=1
  ## plot E_I and E_A together in one graph with error bars, allowing the user to adjust the upper dose limit (<2 Gy; if no protons <1)
	df= data.frame(md = doses, EI=E_I_mean, EA=E_Anm)##dev.off()
	######Elegant ribbon plots###### 
	ggplot(df,aes(md)) + geom_line(aes(y = EI),color ='red') + geom_ribbon(aes(ymin=E_I_low, ymax=E_I_high), alpha=0.2,fill='red')+
	geom_line(aes(y = EA),color ='blue') + geom_ribbon(aes(ymin=E_A_low, ymax=E_A_high), alpha=0.2,fill='blue')+ggtitle("95% CI")
	ggsave(outputfile, device=cairo_ps)# Only applies to Eff=1 
	
	#####plot E_I, and E_A together in one graph for the central parameters, without error bars
	#end=1; by=.1# adjust here RKS
	plot(c(0,end),c(0,40),type='n',bty="l",yaxs='i',xaxs='i',xaxp=c(0,end,2),
       yaxp=c(0,40,4),xlab="dose (cGy)", ylab="Effect (%)",ann=F,las=1)
	lines(doses,E_I_mean,col='red',lw=3,lty=2)
  lines(doses,E_Anm,col='blue',lw=3)
  abline(h=c(10,20,30))
# 	if (Eff==2){
# 	  lines(doses,E_A_mean,lty=3,lw=2)#RKS
# 	  lines(c(0,.002),c(0,E_Anm[1]),lw=3,col='blue')
# 	  lines(c(0,.001),c(0,E_A_mean[1]),lw=3,col='red',lty=3)
# 	}
  
  #Next section estimates the minimum dose where our I(d) is valid for Eff=2 under construction.
#   Eff=2
# 	LE_single = function(dose, L, Eff, proton, param = parameter[[Eff]]){#if Eff=1, E_single=T; #if Eff=2, E_single=T+NT, where NT=kappa*exp(-lamd*dose).  
# 	  if (dose > .001){
#       P_E = T_single(dose,L,Eff,proton)
# 	 
# 	    P_E = P_E + kappa(L,param)*exp(-lamd(L,param)*dose)#P_E here is sum of targeted and non-targeted effects
# 	  }
#	  return(P_E)
#	}
  