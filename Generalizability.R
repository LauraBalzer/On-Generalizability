#============================
# R Code to implement simulations in the Invited Commentary:
# `All generalizations are dangerous, even this one.' - Dumas
# Written by Laura Balzer for Epidemiology 2017 
#
# Programmer: Laura Balzer
# lbbalzer@hsph.harvard.edu
#
# Full reference to Lesko et al: 
# C.R. Lesko, A.L. Buchanan, D. Westreich, J.K. Edwards, M.G. Hudgens, and S.R. Cole. 
# Generalizing study results: a potential outcomes perspective. Epidemiology, 2017.
#============================


rm(list=ls()) 

set.seed(1)

#============================
# get.data: function to simulate the full data 
#		(covariates, enrollment, txt, counterfactual outcomes, observed outcome)
#	input: number of units, 
#		biased.sample (if true, generate according to biased scheme in Lesko et al.; 
#			if false, take a simple random sample)
#		equal.allocation (if true, ensure n/2 units receive the treatment & n/2 receive the control;
#			if false, completely randomize)
# output: full data 
#============================
get.data<- function(n, biased.sample=F, equal.allocation=F){
	
	# generate baseline covariates
	W1 <- rbinom(n, 1, .15)
	W2 <- rbinom(n, 1, .20)
	
	# if a biased sample, then selection probability depends on covariates
	if(biased.sample){
		strata1 <- sample(which(W1==0 & W2==0), 320)
		strata2 <- sample(which(W1==1 & W2==0), 480)
		strata3 <- sample(which(W1==0 & W2==1), 480)
		strata4 <- sample(which(W1==1 & W2==1), 720)

		S<- rep(0, n)
		S[c(strata1, strata2, strata3, strata4)]<- 1
	}else{
		# otherwise select all units
		S<- rep(1,n)
	}
	
	A<- rep(NA, n)
	# if randomize treatment with equal allocation
	#		ie ensure n/2 A=1 and n/2 with A=0
	if(equal.allocation){
		Atemp<- rbinom(sum(S)/2, 1, .5)
		Atemp2<- ifelse(Atemp==1, 0, 1)
		Atemp3<- sample( c(Atemp, Atemp2))
	}else{
		# randomize completely as done by Lesko et al.
		Atemp3<- rbinom(sum(S), 1, .5)
	}
	A[S==1]<- Atemp3
	
	# generate the unmeasured factor U_Y that deterimnes the counterfactual & observed outcome
	UY <- runif(n, 0, 1)
	
	# generate the counterfactual outcomes
	Y1<- get.Y(W1, W2, A=1, UY)
	Y0<- get.Y(W1, W2, A=0, UY)
	
	# the observed outcome is equal to the counterfactual outcome under the observed txt
	Y<- ifelse(A, Y1, Y0)
	
	# set observed exposure and outcome to 0 if not selected/enrolled 
	A[S==0] <- 0 
	Y[S==0] <- 0 
	data.frame(W1, W2, S, A, Y1, Y0, Y)
}

# get.Y: function to generate the outcome
get.Y <- function(W1, W2, A, UY){
	as.numeric( UY < (0.1073 - 0.05*A + 0.2*W1 + 0.2*W2 - 0.15*A*W1*W2))
}

#============================
# Calculate the true value of the population effect with Monte Carlo simulations
out<- get.data(n= 750000)
PATE<- mean(out$Y1 - out$Y0)


#============================
# For varying sample sizes, repeatly generate the full data, calculate the SATE,
#	and for the biased sample, implement the unadjusted, Gcomp & IPW estimator 
nReps<- 5000
SATE.100 <- SATE.500 <- SATE.2000 <- SATE.biased <- Gcomp <- IPW <- unadj <- rep(NA, nReps)
for(r in 1:nReps){
	
	#------------------------------
	# Table 1 - first three rows
	# getting the true value of the sample effect with simple random sampling
	#------------------------------
	temp<- get.data(n=100)
	SATE.100[r]<- mean(temp$Y1 - out$Y0)
	temp<- get.data(n=500)
	SATE.500[r]<- mean(temp$Y1 - out$Y0)	
	temp<- get.data(n=2000)
	SATE.2000[r]<- mean(temp$Y1 - out$Y0)
	
	#--------------------------------
	# Under biased sampling - threat to external validity
	#----------------------------------
	# draw the population of 50000 units
	full.data <- get.data(n=50000, biased.sample=T)
	
	# calculate the SATE for the enrolled units (S==1)
	# Table 1 final row
	SATE.biased[r]<- mean(full.data[full.data$S==1, 'Y1'] - full.data[full.data$S==1, 'Y0'])
	
	# subset on the observed data
	obs.data<- subset(full.data, select=c(W1,W2,S,A,Y))
	
	#-----------------------------------------
	# unadjusted estimator (ave difference in mean outcomes)
	unadj[r]<- mean(obs.data[obs.data$A==1 & obs.data$S==1, 'Y']) - 
		 	mean(obs.data[obs.data$A==0 & obs.data$S==1, 'Y']) 	
		 	
	#-----------------------------------------
	# Gcomputation with fully saturated parametric regression for the outcome 
	outcome.regression <- glm(Y ~ A*W1*W2, family='binomial', data=obs.data[obs.data$S==1, ])
	X.11<- X.10<- obs.data
	X.11$A<- 1; X.11$S<- 1
	X.10$A<- 0; X.10$S<- 1
	predict.outcome.exp<- predict(outcome.regression, newdata=X.11, type='response')
	predict.outcome.unexp<- predict(outcome.regression, newdata=X.10, type='response')
	Gcomp[r]<- mean(predict.outcome.exp - predict.outcome.unexp)
	
	#-----------------------------------------
	# IPW with fully saturated parametric regression for the enrollment (selection) mechanism	
	# and with fully saturated parametric regression for the exposure mechanism
	# (alternative we know that in a randomized trial P(A=1|S=1, W) = 0.5)
	selection.regression<- glm(S~ W1*W2, family='binomial', data=obs.data)
	predict.prob.select <- predict(selection.regression, type='response')	
	exp.regression<- glm(A~ W1*W2, family='binomial', data=obs.data[obs.data$S==1,])
	predict.prob.exp <- predict(exp.regression, newdata=X.11, type='response')
	
	# calculate the denominator for the weights as the product of the probability of being selected
	#	and conditional probability of the obs. exposure
	den.11<- predict.prob.select*predict.prob.exp
	den.10<- predict.prob.select*(1-predict.prob.exp)
	IPW[r] <- mean( 
		( as.numeric(obs.data$A==1 & obs.data$S==1)/den.11 - 
		  as.numeric(obs.data$A==0 & obs.data$S==1)/den.10 )*obs.data$Y )	
	print(r)

}
# save the results
save(PATE, SATE.100, SATE.500, SATE.2000, SATE.biased, unadj, Gcomp, IPW, file='Generalizability_v7Mar.Rdata')

#=====================================
# Create the Table 1 with parameter values
make.table1<- function(est, digits.summary=2, digits.var=3){
	print( round(summary(est)*100, digits.summary))
	print(round(var(est)*100, digits.var))
}
make.table1(SATE.100, digits.summary=1, digits.var=2)
make.table1(SATE.500, digits.summary=1, digits.var=2)
make.table1(SATE.2000, digits.summary=1, digits.var=3)
make.table1(SATE.biased, digits.summary=1, digits.var=3)

# Examine estimator performance as a function of  the causal parameter
performance <- function(estimates, truth) {
	ave <- mean(estimates)
	bias <- mean(estimates - truth)
	var <- var(estimates)
	mse <- mean((estimates - truth) ^ 2)
	z <- bias / sqrt(var)
	cover <- pnorm(1.96 - z) - pnorm(-1.96 - z)
	data.frame(ave,bias,var, mse,cover) * 100 # Return in percent
}

# Evaluating the performance 
performance(estimates=unadj, truth=SATE.biased)
performance(estimates=unadj, truth=PATE)
performance(estimates=Gcomp, truth=PATE)
performance(estimates=IPW, truth=PATE)

# Looking at equivalence of Gcomp and IPW
summary(Gcomp- IPW)