############################################################################################################
#Object: Perform a Progressive-Change BACIPS analysis using step, linear, asymptotic and sigmoid models	#
#Authors: L. Thiault, L. Kernaléguen, C.W. Osenberg & J. Claudet												#
############################################################################################################

# Load packages
require(minpack.lm) # Fitting non-linear models
require(nls2) # Fitting non-linear models
require(AICcmodavg) # calculate second order AIC (AICc)

### The function requires 4 inputs of class integer with the same length (STEP 1) :
	# control: includes response variable measured in the Control site at each sampling time.model;
	# impact: includes response variable measured in the Impact site at each sampling time.model;
	# time.true: time corresponding to each sample;
	# time.model: surveys from the Before period are assigned time.model=0, and surveys from the After period are assigned time.models that represent the time since the intervention started.

### Create the ProgressiveChangeBACIPS function
ProgressiveChangeBACIPS<-function(control, impact, time.true, time.model) 
{
	### STEP 2 - Calculate the delta at each sampling date
	delta <- impact - control
	
	# Plot delta against time.true
	dev.new(width=10, height=5)
	par(mfrow=c(1,2))
	plot(delta~time.true, type="n")
	time.model.of.impact=max(which(time.model==0))
	rect(time.model.of.impact, min(delta)-10, max(time.model)+10, max(delta)+10, col = "grey")
	points(delta~time.true, pch=24, bg="white", cex=2)
	
	### STEP 3 - Fit and compete models
	## Create a 'period' variable
	period <- ifelse(time.model==0, "Before","After")
	
	## Fit a step model
	step.Model<-aov(delta ~ period)

	## Fit a linear model
	linear.Model<-lm(delta ~ time.model)
	
	## Fit an asymptotic model
	# Create an asymptotic function
	myASYfun<-function(delta, time.model)
	{
		funAsy<-function(parS, time.model)	(parS$M * time.model) / (parS$L + time.model) + parS$B
		residFun<-function(p, observed, time.model) observed + funAsy(p,time.model)
		parStart <- list(M=mean(delta[time.model.of.impact:length(time.true)]), B=mean(delta[1:time.model.of.impact]), L=1)
		nls_ASY_out <- nls.lm(par=parStart, fn= residFun, observed=delta, time.model=time.model, control = nls.lm.control(maxfev = integer(), maxiter = 1000))
		foAsy<-delta~(M * time.model) / (L + time.model) + B
		startPar<-c(-coef(nls_ASY_out)[1], coef(nls_ASY_out)[2], coef(nls_ASY_out)[3])
		asyFit<-nls2(foAsy, start=startPar, algorithm="brute-force") # nls2 enables to calculate AICc
		asyFit
	}
	# Fit the asymptotic model
	asymptotic.Model<-myASYfun(delta=delta,time.model=time.model)
	
	
	## Fit a sigmoid model
	## Create a sigmoid function
	mySIGfun<-function(delta, time.model)
	{
		funSIG<-function(parS, time.model)	(parS$M * (time.model/parS$L)^parS$K) / (1 + (time.model/parS$L) ^ parS$K) + parS$B
		residFun<-function(p, observed, time.model) observed + funSIG(p,time.model)
		parStart <- list(M=mean(delta[time.model.of.impact:length(time.true)]), B=mean(delta[1:time.model.of.impact]), L=mean(time.model), K=5)
		nls_SIG_out <- nls.lm(par=parStart, fn= residFun, observed=delta, time.model=time.model, control = nls.lm.control(maxfev = integer(), maxiter = 1000))
		foSIG<-delta~(M * (time.model/L) ^ K) / (1 + (time.model/L) ^ K) + B
		startPar<-c(-coef(nls_SIG_out)[1],-coef(nls_SIG_out)[2],coef(nls_SIG_out)[3],coef(nls_SIG_out)[4])
		sigFit<-nls2(foSIG, start=startPar, algorithm="brute-force") # nls2 enables to calculate AICc
		sigFit
	}
	# Fit the sigmoid model
	sigmoid.Model<-mySIGfun(delta=delta,time.model=time.model)


	## Compete models
	# Perform AIC tests
	AIC.test<-AIC(step.Model, linear.Model, asymptotic.Model, sigmoid.Model)
	AICc.test<-as.data.frame(cbind(AIC.test[,1], c(AICc(step.Model), AICc(linear.Model), AICc(asymptotic.Model), AICc(sigmoid.Model))))
	rownames(AICc.test)<-rownames(AIC.test)
	names(AICc.test)<-names(AIC.test)
	
	# Calculate AICc weight and selected the best model
	for(i in 1:dim(AICc.test)[1])
	{
		AICc.test$diff[i]<-AICc.test$AIC[i]-min(AICc.test$AIC)
	}
	AICc.test$RL<-exp(-0.5* AICc.test$diff)
	RL_sum<-sum(AICc.test$RL)
	AICc.test$aicWeights<-(AICc.test$RL/RL_sum)*100
	w<-AICc.test$aicWeights
	names(w)<-rownames(AICc.test)

	# Display AICc weights
	print(w)
	barplot(w, col="white", ylab="Relative likelihood (%)", cex.names = 0.9, names.arg =c("Step","Linear","Asymptotic","Sigmoid"))
	best.Model<-which(w==max(w))
	
	### STEP 4 - Derive inference based on the best model (i.e., with the higher AICc weight)
	if(best.Model==1) {writeLines(paste("\n\nSTEP MODEL SELECTED - Likelihood = ", round(w[1],1), "%\n\n", sep=""))
		print(summary(step.Model))}
	if(best.Model==2) {writeLines(paste("\n\nLINEAR MODEL SELECTED - Likelihood = ", round(w[2],1), "%\n\n", sep=""))
		print(summary(linear.Model))}
	if(best.Model==3) {writeLines(paste("\n\nASYMPTOTIC MODEL SELECTED - Likelihood = ", round(w[3],1), "%\n\n", sep=""))
		print(asymptotic.Model)}
	if(best.Model==4) {writeLines(paste("\n\nSIGMOID MODEL SELECTED - Likelihood = ", round(w[4],1), "%\n\n", sep=""))
		print(sigmoid.Model)}
	
	w <<- w
	asymptotic.Model <<- asymptotic.Model
	step.Model <<- step.Model
	linear.Model <<- linear.Model
	sigmoid.Model <<- sigmoid.Model
	
	}





### Implement the ProgressiveChangeBACIPS function using exemples on simulated datasets (R=1)

# # Exemple with a positive step change	
# control<-rpois(20, 100) # Nbefore=100
# impact<-rpois(20, 100) # Nbefore=100
# impact<-c(impact[1:5], impact[6:20]+rep(100, 15))
# time.true<-seq(1:20)
# time.model<-c(rep(0,5), seq(1,20, length.out=15)) # Sbefore=5; Safter=15
# ProgressiveChangeBACIPS(control, impact, time.true, time.model)
# 
# 
# # Exemple with a positive sigmoid change
# control<-rpois(15, 100) # Nbefore=100
# impact<-rpois(15, 100) # Nbefore=100
# time.true<-seq(1:15)
# time.model<-c(rep(0,5), seq(1,15, length.out=10)) # Sbefore=5; Safter=10
# M=100 # Effect magnitude (Etrue=100)
# K=5 # Steepness of the slope
# L=5 # Half-life time (based on time.model, not time.true)
# impact <- c(impact[1:5], M * (time.model[6:15]/L) ^ K / (1 + (time.model[6:15]/L) ^ K) + impact[6:15])
# ProgressiveChangeBACIPS(control, impact, time.true, time.model)

# End of the script