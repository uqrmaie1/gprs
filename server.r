# server.r

# note: for upload boxes, the default file size limit is 5 MB
# options(shiny.maxRequestSize = 9*1024^2)
# will raise the limit to 9 MB


# load the necessary libraries
library(rms)
library(pROC)
library(MBESS)
library(MASS)
library(ggplot2)
#library(shinyRGL)
#library(rgl)
library(reshape2)

# source dudbridge code containing functions, i.e. AEGP code
source("polygenescore.r")


# function to disable ActionButtons
disableActionButton <- function(id,session) {
  session$sendCustomMessage(type="jsCode",
                            list(code= paste("$('#",id,"').prop('disabled',true)"
                                             ,sep="")))
}

# function to enable ActionButtons
enableActionButton <- function(id,session) {
  session$sendCustomMessage(type="jsCode",
                            list(code= paste("$('#",id,"').prop('disabled',false)"
                                             ,sep="")))
}

clean = function() {
  for(file in c('GPRS.csv', 'GPRS_deciles.csv', 'Vg2FromP.csv', 'CorrFromP.csv', 'www/plot.pdf')) {
  if(file.exists(file)) {
    file.remove(file)
  }
  }
}

getExDir = function() {'./www/examples/'}

getVarDescription = function(vars) {
  desc = c(scoreFile   = 'name of profile score file',
           pThreshold  = 'p-value Threshold',
           NKr2        = 'Nagelkerke R2',
           pval        = 'Nagelkerke R2 p-value',
           h2l_r2n     = 'Heritability of liability explained by profile score',
           se_h2l_r2   = 'standard error of h2l_r2n',
           AUC         = 'area under the ROC curve of the profile score (not adjusted for covariates)',
           OR10decile  = 'odds ratio between highest and lowest decile',
           ORL95       = 'lower end of 95% confidence interval',
           ORH95       = 'higher end of 95% confidence interval',
           vg          = 'estimate of the genetic variance explained by all SNPs',
           vgLo        = 'lower end of 95% confidence interval',
           vgHi        = 'higher end of 95% confidence interval',
           rG          = 'estimate of genetic correlation',
           rGLo        = 'lower end of 95% confidence interval',
           rGHi        = 'higher end of 95% confidence interval'
           )
  mat = matrix(NA, 0, 2)
  for(var in vars) {
    if(var %in% names(desc)) {
      mat = rbind(mat, c(var, desc[var]))
    }
  }
  colnames(mat) = c('variable', 'description')
  as.data.frame(mat)
}

# main code
shinyServer(function(input, output, session) {
  
  clean()  

	# function that does the main calculations, creates "GPRS.csv"
	originalScriptBQ <- function(numFiles,						# number of files
	                           scoreFilenames,				# vector of uploaded score file filenames
	                           covariatesFilenames,	# vector of uploaded covariate file filenames
	                           dataTables,						# list containing tables (one table per file/file pair)
	                           covariatesSelected,		# vector of selected covariate columns' names
	                           targetK								# baseline risk of the target trait, e.g. 0.01 for schizophrenia, N/A for quantitative data (K = 0 for these)
	){
	  # Sunny's comments:
	  # added annotations to most lines, trying to explain what they do
	  #
	  # removed all code not related to the production of GPRS.csv
	  # genericized inputs and glm/lm formulas
	  # this code is now a function
	  #
	  #
	  #
	  # np currently hard-coded to 1 during linear model R2 calculation
	  
	  ###########################################################
	  #PGC_GPRSincORse.R 
	  #Genomic Profile Risk Score analysis
	  #Considering covariates
	  #Hong Lee & Naomi Wray December 2013 updated January 2014
	  ###########################################################
	  
	  ### functions
	  h2l_R2 <- function(k, r2, p) {
	    #K baseline disease risk
	    #r2 from a linear regression model attributable to genomic profile risk score
	    #P proportion of sample that are cases
	    #calculates proportion of variance explained on the liability scale
	    #from ABC at http://www.complextraitgenomics.com/software/
	    #Lee SH, Goddard ME, Wray NR, Visscher PM. (2012) A better coefficient of determination for genetic profile analysis. Genet Epidemiol. 2012 Apr;36(3):214-24.
	    x= qnorm(1-K)
	    z= dnorm(x)
	    i=z/K
	    C= k*(1-k)*k*(1-k)/(z^2*p*(1-p))
	    theta= i*((p-k)/(1-k))*(i*((p-k)/(1-k))-x)
	    h2l_R2 = C*r2 / (1 + C*theta*r2)
	  } 
	  
	  se_h2l_R2 <- function(k,h2,se, p) {
	    #K baseline disease risk
	    #r2 from a linear regression model attributable to genomic profile risk score
	    #P proportion of sample that are cases
	    #calculates proportion of variance explained on the liability scale
	    #from ABC at http://www.complextraitgenomics.com/software/
	    #Lee SH, Goddard ME, Wray NR, Visscher PM. (2012) A better coefficient of determination for genetic profile analysis. Genet Epidemiol. 2012 Apr;36(3):214-24.
	    
	    #SE on the liability (From a Taylor series expansion)
	    #var(h2l_r2) = [d(h2l_r2)/d(R2v)]^2*var(R2v) with d being calculus differentiation
	    x= qnorm(1-K)
	    z= dnorm(x)
	    i=z/K
	    C= k*(1-k)*k*(1-k)/(z^2*p*(1-p))
	    theta= i*((p-k)/(1-k))*(i*((p-k)/(1-k))-x)
	    se_h2l_R2 = C*(1-h2*theta)*se
	  }
	  
	  h2l_R2N <- function(k, r2n, p) {
	    #k baseline disease risk
	    #r2n Nagelkerke's attributable to genomic profile risk score
	    #proportion of sample that are cases
	    #calculates proportion of variance explained on the liability scale
	    #from ABC at http://www.complextraitgenomics.com/software/
	    #Lee SH, Goddard ME, Wray NR, Visscher PM. (2012) A better coefficient of determination for genetic profile analysis. Genet Epidemiol. 2012 Apr;36(3):214-24.
	    x <- qnorm(1 - k)
	    z <- dnorm(x)
	    i <- z / k
	    cc <- k * (1 - k) * k * (1 - k) / (z^2 * p * (1 - p))
	    theta <- i * ((p - k)/(1 - k)) * (i * ((p - k) / ( 1 - k)) - x)
	    e <- 1 - p^(2 * p) * (1 - p)^(2 * (1 - p))
	    h2l_R2N <- cc * e * r2n / (1 + cc * e * theta * r2n) 
	  }
	  
	  ### main program
	  
	  nj = numFiles # number of files
	  
	  # initialize the data frame and output file it will be writing to at the end of each iteration of the main for loop
	  # This output file gives the measures that should be used
	  O=data.frame("scoreFile",        # scoreFilename = score input filename (FID IID PHENO CNT CNT2 SCORE)
	               "covariatesFile",  # covariatesFilename = covariate input filename (FID IID covariates..)
	               "N",               # N = total sample size
	               "Propcase",        # P = proportion of target sample that are cases and not controls ("Propcase" column)
	               "NKr2",            # NKv = Nagelkerke's R2 (R2 measure proposed in the original paper by Purcell, Wray et al in Nature where they first described the genetic profile score method)
	               "pval",            # pval = p-value for the significance of the Nagelkerkes R2
	               "PopRisk",         # K = population risk (prevalence) of disease used in converting NKv to liability scale
	               "h2l_r2n",         # h2l_r2n = proportion of variance explained (R2) by the score on the liability scale, calculated from NKv (this is the NKr2 converted to the liability scale)
	               "se_h2l_r2",       # se_h2l_r2 = standard error of the R2 on the liability scale
	               "AUC",             # aucvS = what we think is the most appropriate estimate of AUC (area under the Receiver Operator Curve) attributed to the score
	               "OR10decile",      # ORD = the odds ratio when comparing top to bottom decile
	               "ORL95",           # ORDL = lower boundary of the CI of the ORD
	               "ORH95",           # ORDH = upper boundary of the CI of the ORD
	               "pval_r2",         # pval_r2 = p-value for the significance of the R2 (##AV change)
	               "R2v",
                 "vr"
                 #"r2",              # r2 = proportion of variance explained (R2) by the score (##AV added AV)
	               #"se_r2"            # se_r2 = standard error of the R2 (##AV added AV)
	  )
	  write.table(O,"GPRS.csv",row.names=F,col.names=F,sep=",")
	  
	  # initialize data frame and output csv for decile data
	  OO=data.frame("scoreFile", "decile", "estimate", "lowerCI", "upperCI")
	  
	  write.table(OO, "GPRS_deciles.csv", row.names=F, col.names=F, sep=",")
	  
	  # the main for loop, each j is for another file
	  for (j in (1:nj)){
	    K = targetK[j]
	    scoreFilename = scoreFilenames[j]
	    covariatesFilename = covariatesFilenames[j]
	    ### these filenames will be in the output to let the user know which of their files generated the output
	    
	    ri <- dataTables[[j]] # ri is the main working frame, assigned to the jth item in dataTables (a function parameter)
      binary = length(unique(ri$PHENO)) == 2 # determine if this file has a binary phenotype
	    
	    originalLength <- length(ri)
	    # save original length to know which column is the end of covariates, 
	    # because additional columns will be added in the analysis later, increasing the length
	    
	    # main calculations begin
	    ###normalize the score
	    (ri$SCORE-mean(ri$SCORE))/sd(ri$SCORE)->ri$NSCORE # NSCORE is normalized score (i.e. how many SDs away)
	    
      ri$PHENO1=ri$PHENO-1 # new column PHENO1 from PHENO, converting plink's 1-2 coding to 0-1 coding **** only for binary data
	    P=sum(ri$PHENO1)/length(ri$PHENO1) # proportion of target sample that are cases **** only for binary data
	    # P used in output. inputs used: PHENO
	    
	    # define N, missing in Anna's code
	    N=length(ri$PHENO1) # number of entries, essentially
      
	    # genericize formulas and adapt as covariates are selected/unselected
	    tstFchar <- paste("PHENO1 ~ NSCORE + ", paste(covariatesSelected, collapse="+"), sep="")
	    tstFform <- as.formula(tstFchar)
	    
	    tstRchar <- paste("PHENO1 ~ ", paste(covariatesSelected, collapse="+"), sep="")
	    tstRform <- as.formula(tstRchar)
	    
	    if(!binary) {
  	    tstF = glm(tstFform,        data = ri) #linear model
  	    tstS = glm(PHENO1 ~ NSCORE, data = ri) #linear model
  	    tstR = glm(tstRform,        data = ri) #linear model
	    }
      
	    if(binary) {
  	    ###logistic models
  	    tstF = glm(tstFform,        data = ri, family = binomial(logit)) #logit model
  	    tstS = glm(PHENO1 ~ NSCORE, data = ri, family = binomial(logit)) #logit model
  	    tstR = glm(tstRform,        data = ri, family = binomial(logit)) #logit model
  	    # F fits NSCORE + covariates, S just NSCORE, R just covariates
	    
      
  	    #library(pROC)
  	    aucvS = auc(ri$PHENO1,tstS$linear.predictors)
  	    # area under curves for S (NSCORE). 
  	    # auc(response, predictors)
  	    # aucvS used in output. inputs used: PHENO, SCORE
  	    
  	    #Cox&Snell R2 
  	    LLF=logLik(tstF) # "extract log-likelihood" from logistic model F (everything)
  	    LLR=logLik(tstR) # from logistic model R (covariates)
  	    CSv=1 - exp(  (2/N) * (LLR[1]-LLF[1])  ) # some math calculates Cox & Snell R2 involving model R (cov) and F (all)
  	    
  	    #Nagelkerke's R2
  	    NKv= CSv/ (1 - exp(  (2/N) * LLR[1])  ) 
  	    # NKv used in output. inputs used: PHENO, SCORE, covariates
  	    
  	    #pvalue
  	    devdiff=tstR$deviance-tstF$deviance
  	    df=tstR$df.residual-tstF$df.residual
  	    pval=pchisq(devdiff,df,lower.tail=F) # chi-squared distribution. devdiff = quantity, df = degrees of freedom
  	    # pval used in output. inputs used: PHENO, SCORE, covariates
      }
      
      
	    #linear model R2 *********************************************
	    std_y=ri$PHENO1 # an alias for PHENO1 to be used in the next line
	    ri$std_y=(std_y-mean(std_y))/sd(std_y) # new column "std_y" in ri, how many SDs away is the PHENO1
	    
	    # genericize formulas and adapt to covariates selected
	    lmfChar <- paste("std_y ~ NSCORE + ", paste(covariatesSelected, collapse="+"), sep="")
	    lmfForm <- as.formula(lmfChar)
	    
	    lmrChar <- paste("std_y ~ ", paste(covariatesSelected, collapse="+"), sep="")
	    lmrForm <- as.formula(lmrChar)
	    
	    lmf=lm(lmfForm, data = ri)
	    lmr=lm(lmrForm, data = ri)
	    # linear regression models fitting everything (f), covariates (r), or 1 (0) to ri$std_y
	    
	    R2v=1-exp((2/N)*(logLik(lmr)[1]-logLik(lmf)[1]))
	    # math for R2 values of the linear regression
	    
	    #standard error of R2v
	    #from Olkin and Finn (Psychological Bulletin, 1995, Vol. 118, No. 1, 155-164)
	    np=1    #number of parameters
	    # np is hard-coded = 1, only used in vr calculation (below), which is in turn only used in se_h2l_r2 calculation, which is in output
	    vr= 4/length(std_y) * R2v * (1-R2v)^2 * (  1-(2*np+3) / length(std_y)  )
	    test <- anova (tstR, tstF, test="Chisq" ) ##AV added
	    pval_r2 <- test[2,5]                         ##AV added
	    # p-value
      
	    if(binary) {
  	    #calculate liability R2
  	    h2l_r2 = h2l_R2(K,R2v,P) #linear model
  	    se_h2l_r2=se_h2l_R2(K,h2l_r2,vr^.5,P)
  	    # se_h2l_r2 used in output. inputs used: hard-coded (K, np), PHENO, SCORE, covariates
  	    
  	    h2l_r2n = h2l_R2N(K,NKv,P) #nagelkerke's
  	    # h2l_r2n used in output. inputs used: hard-coded (K), PHENO, SCORE, covariates
  	    
  	    #make deciles
  	    oNSCORE=ri$NSCORE[order(ri$NSCORE)] # copy of NSCORE column, ordered by NSCORE (default ascending order)
  	    oPHENO1=ri$PHENO1[order(ri$NSCORE)] # copy of PHENO1 column, ordered by NSCORE
  	    rio=ri[order(ri$NSCORE),] # copy of working data frame, ordered by NSCORE
  	    N10=round(N/10) # bin size for deciles
  	    dumv=matrix(0,length(oNSCORE),9) #dummy varaible
  	    # a matrix of length rows x 9 columns, initialized to 0
  	    for (zi in 1:9) {
  	      fst=length(oNSCORE)-zi*N10+1 # imagine length = 100. when zi = 1, fst = 100-10+1 = 91
  	      lst=length(oNSCORE)-zi*N10+N10 # and lst = 100. then the next expression will populate rows 91-100/column 1 with "1"
  	      # and when zi = 9, fst = 100-90+1 = 11 and lst = 20, populating these rows in column 9 with "1"
  	      dumv[fst:lst,zi]=1
  	      # each iteration in the loop is another column, where "1" is assigned to cells of rows fst to lst
  	      # this creates a matrix where each decile has its own column, with its members having value "1" and everyone else "0"
  	      # 1st column has rows farthest from 0th row
  	      # since our row data will be ordered ascendingly, this would make column 1 the top decile
  	      # the bottom decile is not represented in this matrix
  	    }
  	    
  	    # genericize formula and adapt to covariates selected
  	    tstFchar <- paste("PHENO1 ~ dumv + ", paste(covariatesSelected, collapse="+"), sep="")
  	    tstFform <- as.formula(tstFchar)
  	    
  	    # binomial logistic regressions using data ordered by NSCORE
  	    tstF = glm(tstFform, data = rio,family = binomial(logit)) #logit model
  	    # tests whether membership in deciles (column number lower = top decile) correlate with PHENO1 ordered ascendingly by NSCORE
  	    
  	    # when glm takes dumv (a matrix) in its fit, it treats every column in dumv as a new variable
  	    # tstF then has coefficients and std errors, etc. for dumv1, dumv2, ..., dumv9
  	    # these represent the first column (highest decile or 10th decile) to the 9th column (2nd lowest decile, or 2nd decile, since the lowest decile is not represented)
  	    
  	    # tstF is a glm object that includes many other objects like $coefficients (a numeric vector)
  	    # tstF$coefficients looks like this: [1] (Intercept) [2] dumv1 [3] dumv2 ... [10] dumv9 [11] C1 [12] C2 ... [15] C5
  	    ORD=exp(tstF$coefficients[2])
  	    # accesses the coefficient of dumv1 (the second item in $coefficients), and raises e to its power (exp function)
  	    
  	    # summary() is used to summarize model fitting functions like glm, and in turn produces an object containing other objects like $coefficients
  	    # summary(tstF)$coefficients is a matrix that has 
  	    #		rows [1] (Intercept), [2-9] dumv1-9, [10-15] C1-5 (of course this depends on how many coefficients there are)
  	    #		columns [1] Estimate, [2] Std. Error, [3] z value, and [4] Pr(>|z|)
  	    ORDL=exp(tstF$coefficients[2]-1.96*summary(tstF)$coefficients[2,2])
  	    # accesses the Std. Error of dumv1 (summary(tstF)$coefficients[2,2])
  	    # multiply std error (which is one std dev) by 1.96 (in a normal distribution, you go above and below the mean 1.96 std deviations to capture 95% of the data)
  	    # subtract that from tstF$coefficients[2], or the estimate of dumv1, which gives you the lower boundary of the 95% confidence interval
  	    # and raise e to its power
  	    ORDH=exp(tstF$coefficients[2]+1.96*summary(tstF)$coefficients[2,2])
  	    # same as above, but add instead of subtract, to get the higher boundary of the 95% confidence interval
  	    
  	    thDecile <- 10
  	    for(i in 2:10){
  	      assign(paste0("ORD",thDecile,collapse=""), exp(tstF$coefficients[i]) )
  	      assign(paste0("ORD",thDecile,"L",collapse=""), exp(tstF$coefficients[i]-1.96*summary(tstF)$coefficients[i,2]) )
  	      assign(paste0("ORD",thDecile,"H",collapse=""), exp(tstF$coefficients[i]+1.96*summary(tstF)$coefficients[i,2]) )
  	      thDecile <- thDecile - 1
  	    }
	    }
	    
	    #if(binary) {pval_r2=R2v=vr=NA}
      if(!binary) {NKv=pval=K=h2l_r2n=se_h2l_r2=aucvS=ORD=ORDL=ORDH=NA}
      
	    O=data.frame(scoreFilename,covariatesFilename,N,P,NKv,pval,K,h2l_r2n,se_h2l_r2,aucvS,ORD,ORDL,ORDH,pval_r2,R2v,vr)
	    write.table(O,"GPRS.csv",row.names=F,col.names=F,sep=",",append=T)
	    
      if(binary) {
  	    # OO has headings scoreFile, decile, estimate, upperCI, lowerCI
  	    OO=data.frame(scoreFilename, 1, 1, 1, 1)
  	    write.table(OO, "GPRS_deciles.csv", row.names=F, col.names=F, sep=",", append=T)
  	    for(i in 2:10){
  	      OO=data.frame( scoreFilename, i, get(paste0("ORD",i)), get(paste0("ORD",i,"L")), get(paste0("ORD",i,"H")) )
  	      write.table(OO, "GPRS_deciles.csv", row.names=F, col.names=F, sep=",", append=T)
  	    }
      }
	  } # end for loop
	}
  
  # function that does the main calculations, creates "GPRS.csv"
  simplifiedScriptB <- function(numFiles,						# number of files
                               scoreFilenames,				# vector of uploaded score file filenames
                               covariatesFilenames,	# vector of uploaded covariate file filenames
                               dataTables,						# list containing tables (one table per file/file pair)
                               covariatesSelected,		# vector of selected covariate columns' names
                               targetK								# baseline risk of the target trait, e.g. 0.01 for schizophrenia, N/A for quantitative data (K = 0 for these)
  ){
    
  }

	# determines if single (just score file)/two files (score file + covariates file) uploaded, calls the appropriate function (processOneUpload/processTwoUploads), returns "data"
	# it's a reactive expression instead of a function to take advantage of lazy evaluation
	# also determines whether "D" (coded 0,1, used in Anna's example files) or "PHENO" (coded 1,2) is used as column name for disease status, and standardize to "PHENO" (1,2)
	# if "D.ls" used (Anna's example Q files), rename to "PHENO" (quantitative)
	processUpload <- reactive({
    # return NULL if input$file1 is empty
    if(is.null(input$file1)) {
      if(input$exampleButton == 0) {
        return(NULL)
    	} else {
    	  data <- processExample()
    	}
      
    }
		# if neither upload slots are empty, we must be using score+covariates files
		if (!is.null(input$file1) & !is.null(input$file2)) {
			# combine score and covariate files, assign return to variable data
			data <- processTwoUploads(upScore = input$file1, upCovariate = input$file2, upLabel = input$file3)

			# if Anna's "D" column instead of "PHENO", adjust 0,1 coding to 1,2 and rename back to "PHENO"
			if ( length( names(data$tables[[1]])[names(data$tables[[1]])=="D"] ) > 0 ) {
				# for every table...
				for(i in 1:length(data$tables)) {
					# adjust phenotype coding 1,0 to 2,1
					data$tables[[i]]$D <- data$tables[[i]]$D + 1
					# rename "D" back to "PHENO"
					names(data$tables[[i]])[names(data$tables[[i]])=="D"] <- "PHENO"
				}
			}

			# if Anna's "D.ls" column (quantitative trait) instead of "PHENO", rename back to "PHENO"
			if("D.ls" %in% names(data$tables[[1]])){
				# for every table...
				for(i in 1:length(data$tables)){
					#rename "D.ls" back to "PHENO"
					names(data$tables[[i]])[names(data$tables[[i]])=="D.ls"] <- "PHENO"
				}
			}
		} 


		# if first slot filled but second empty, we must be using a single file
		else if (!is.null(input$file1)) {
			data <- processOneUpload(upScore = input$file1, upLabel = input$file3)

			# if Anna's "D" column instead of "PHENO", adjust 0,1 coding to 1,2 and rename back to "PHENO"
			if ( length( names(data$tables[[1]])[names(data$tables[[1]])=="D"] ) > 0 ) {
				# for every table...
				for(i in 1:length(data$tables)) {
					# adjust phenotype coding 1,0 to 2,1
					data$tables[[i]]$D <- data$tables[[i]]$D + 1
					# rename "D" back to "PHENO"
					names(data$tables[[i]])[names(data$tables[[i]])=="D"] <- "PHENO"
				}
			}

			# if Anna's "D.ls" column (quantitative trait) instead of "PHENO", rename back to "PHENO"
			if("D.ls" %in% names(data$tables[[1]])){
				# for every table...
				for(i in 1:length(data$tables)){
					#rename "D.ls" back to "PHENO"
					names(data$tables[[i]])[names(data$tables[[i]])=="D.ls"] <- "PHENO"
				}
			}
		}
    
    data$check = checkScoreFiles(data)
		
    data$quant <- NULL
		q.count = b.count = 0
    
		for(df in data$tables) {
		  if(length(unique(df$PHENO)) == 2)
		    b.count = b.count+1
		  else
		    q.count = q.count+1
		}
		type = 'both binary and quantitative'
		if(b.count == 0) {data$quant <- TRUE}
		if(q.count == 0) {data$quant <- FALSE}
		return(data)
	})
  
  # checks input score files; returns string with error message
  checkScoreFiles = function(data) {
    msg = ''
    dflist = data$tables
    req.colnames = c('FID', 'IID', 'PHENO', 'SCORE')
    for(i in 1:length(dflist)) {
      missing = setdiff(req.colnames, names(dflist[[i]]))
      if(length(missing) > 0) {
        msg = paste(msg, 'column', missing[1], 'missing in score file', data$scoreFilenames[i],'!')
        return(msg)
      }
    }
    if(!is.null(input$file3)) {
      req.colnames = c('scoreFile','discoverySample','targetSample','discoveryTrait','discoveryK','targetTrait','targetK','dtTraitCorrelation','Vg1','Vg2','discoveryCaseN','discoveryControlN','pThreshold','nIndependentSNPs')
      missing = setdiff(req.colnames, names(getLabelFile()))
      if(length(missing) > 0)
        msg = paste(msg, 'column', missing[1], 'missing in label file!')
    }
    return(msg)
  }

	# function to process the upload data (merge score and covariate files, reorder columns, extract K info) for when both score and covariate files are used
	# called by processUpload
	# returns a single list object containing: $tables (list), $scoreFilenames (vector), $covariatesFilenames (vector), $K (vector)
	processTwoUploads <- function(upScore, upCovariate, upLabel) {
		# use dataframes generated by upload boxes as parameters, in this case input$file1 and input$file2

		# declare mergedReorderedColumns first for later use in the for loop
		mergedReorderedColumns <- list()

		for (i in (1:length(upScore$name))){
			# for every file uploaded...
			
			# read the ith file in upScore/upCovariate and assign that dataframe to variables x and y
			x <- read.table(upScore$datapath[i], header = T)
			y <- read.table(upCovariate$datapath[i], header = T) 

			# merge x and y and assign to dataframe xy
			xy <- merge(x,y)
			
			# reorder the columns in xy and add this dataframe to the ith position in mergedReorderedColumns
			# column 1 = FID, 2 = IID, 3 = PHENO, 6 = SCORE, 7+ = covariates, 4 = CNT, 5 = CNT2
			mergedReorderedColumns[[i]] <- xy[c(1,2,3,6:length(xy),4,5)]
			# this arranges columns thus: FID IID PHENO SCORE ..covariates.. CNT CNT2
			# this is because CNT and CNT2 will be processed as covariates later on
			# from now on, positions 5 to length() will all be covariates
		}
		# now we end up with mergedReorderedColumns, a list of dataframes in the format we want
		# upScore and upCovariate still have filename information that we want
		# and we want target K information from the label file

		# read label file (for K information)
		labelTable <- read.csv(upLabel$datapath, header = T, sep = ",")
		# only those entries where uploaded score file name matches that in the label file will be merged
		matchedLabels <- merge(upScore$name,labelTable)
		# now we have a row for each uploaded score file, including a column "K"



		data <- list("tables" = mergedReorderedColumns,
									"scoreFilenames" = upScore$name,
									"covariatesFilenames" = upCovariate$name,
									"targetK" = matchedLabels$targetK
		) 

		return(data)
	}

	# function to process the upload data (only score files uploaded)
	# called by processUpload
	# returns a single list in same format as above (but covariatesFilenames <- "")
	processOneUpload <- function(upScore, upLabel) {
		# for every file, read as table and write to a dataframe in list "tables"
		tables <- list()
		for (i in 1:nrow(upScore)){
			tables[[i]] <- read.table(upScore$datapath[i], header = TRUE)
			# for read.table(), default sep="", which means any whitespace (tab or series of spaces)
			# this is more compatible than using read.csv(sep="\t"), which only reads tab-delimited and not space-delimited
		}

		# read label file (purpose is to extract K information)
    if(!is.null(upLabel))
  		labelTable <- read.csv(upLabel$datapath, header = T, sep = ",")
		# only those entries where uploaded score file name matches that in the label file will be merged
		# but first, we need the column names to match
		# the first element of upScore used to be "filename", but here we rename it to "scoreFile"
		x <- upScore
		names(x)[1] <- "scoreFile"
		# now both x and labelTable have a column named "scoreFile" which they can use as an ID column in merge()
		if(!is.null(upLabel))
      matchedLabels <- merge(x,labelTable)
    else
      matchedLabels = list(targetK=0)
		# now we have a row for each uploaded score file, including a column "targetK"

		# standard, for input into originalScript function
		data <- list("tables" = tables,
									"scoreFilenames" = upScore$name,
									"covariatesFilenames" = "",
									"targetK" = matchedLabels$targetK
		)

		return(data)
	}

  getLabelFile = reactive({
    if(!is.null(input$file3)) {
      labelFile <- input$file3$datapath
      labelTable <- read.csv(labelFile, header=TRUE, sep=",")
    } else {
      dir = getExDir()
      labelTable <- read.csv(paste0(dir, 'disease1_label.csv'), header = T, sep = ",")
    }
    labelTable
  })
  
	# take uploaded label file and apply label columns to GPRS.csv and GPRS_deciles.csv
	# also formats for 3 decimal places in both GPRS.csv and GPRS_deciles.csv
	# called by "run" observer
	processLabelFile <- function() {
	  labelTable = getLabelFile()
    x <- read.csv("GPRS.csv", sep=",", header=TRUE)
		y <- read.csv("GPRS_deciles.csv", sep=",", header=TRUE)

		a <- merge(labelTable,x)
		# b <- merge(labelTable,y) # not necessary to append label file to deciles data

		# do the formatting and write to .csv
		a <- format(a, digits=3)
		write.table(a, "GPRS.csv", row.names=F, sep=",")
		y <- format(y, digits=3)
		write.table(y, "GPRS_deciles.csv", row.names=F, sep=",")
	}

	# first used in the decile graph to select data, made into a function so it's easier to do the same for other graphs
	# takes 5 strings to use as names for the 5 selectInput widgets, enclosed in a wellPanel
	select5 <- function(name1,name2,name3,name4,name5) {
		x <- read.csv("GPRS.csv", header = T, sep = ",")

		wellPanel(
			h5("Data Selection"),
			selectInput(name1, 
				"Include score files that use the following discovery sample(s)...",
				choices = levels(factor(x$discoverySample)),
				selected = levels(factor(x$discoverySample))[1],
				width = "100%",	selectize = TRUE,	multiple = TRUE
			),
			selectInput(name2, 
				"...if they describe the following trait(s)...",
				choices = levels(factor(x$discoveryTrait)),
				selected = levels(factor(x$discoveryTrait))[1],
				width = "100%", selectize = TRUE, multiple = TRUE
			),
			selectInput(name3, 
				"...And use the following target sample(s)...",
				choices = levels(factor(x$targetSample)),
				selected = levels(factor(x$targetSample))[1],
				width = "100%",	selectize = TRUE,	multiple = TRUE
			),
			selectInput(name4, 
				"...if they describe the following trait(s)...",
				choices = levels(factor(x$targetTrait)),
				selected = levels(factor(x$targetTrait))[1],
				width = "100%",	selectize = TRUE,	multiple = TRUE
			),
			selectInput(name5, 
				"...All with SNPs that pass the following p-value threshold(s):",
				choices = levels(factor(x$pThreshold)),
				selected = levels(factor(x$pThreshold))[length(levels(factor(x$pThreshold)))],
				width = "100%",	selectize = TRUE,	multiple = TRUE
			)
		) # end wellPanel
	}
  
  # function to read in example files
  # called by processUpload
  # returns a single list in same format as above (but covariatesFilenames <- "")
  processExampleFiles <- function() {
    dir = getExDir()
    # for every file, read as table and write to a dataframe in list "tables"
    files = list.files(dir, pattern='^example_*')
    # limit to no more than 10 example files, just so it runs faster
    files = na.omit(files[1:min(10, length(files))])
    tables <- list()
    for (i in 1:length(files)){
      tables[[i]] <- read.table(paste0(dir, files[i]), header = TRUE)
      # for read.table(), default sep="", which means any whitespace (tab or series of spaces)
      # this is more compatible than using read.csv(sep="\t"), which only reads tab-delimited and not space-delimited
    }
    
    # read label file (purpose is to extract K information)
      labelTable <- read.csv(paste0(dir, 'disease1_label.csv'), header = T, sep = ",")
    # only those entries where uploaded score file name matches that in the label file will be merged
    # but first, we need the column names to match
    # the first element of upScore used to be "filename", but here we rename it to "scoreFile"
    x <- data.frame(files)
    names(x)[1] <- "scoreFile"
    # now both x and labelTable have a column named "scoreFile" which they can use as an ID column in merge()
    matchedLabels <- merge(x,labelTable)

    # now we have a row for each uploaded score file, including a column "targetK"
    
    # standard, for input into originalScript function
    data <- list("tables" = tables,
                 "scoreFilenames" = files,
                 "covariatesFilenames" = "",
                 "targetK" = matchedLabels$targetK
    )
    
    return(data)
  }

	# take a dataframe 'y' that includes column "scoreFile", subset it based on intersection of 5 criteria (sel1-5), return the resulting dataframe
	# meant to be used in combination with select5 above. select5 makes the UI elements that collect input, and selectDataFrom5 processes that input
	selectDataFrom5 <- function(y, sel1, sel2, sel3, sel4, sel5) {
		x <- read.csv("GPRS.csv", header = T, sep = ",")

		# subset x based on intersection of five traits selected
		xx <- subset(x, discoverySample %in% sel1 &
										discoveryTrait %in% sel2 &
										targetSample %in% sel3 &
										targetTrait %in% sel4 &
										pThreshold %in% sel5
		)

		# subset decile data by filenames in selected subset of x (i.e. xx)
		d <- subset(y, scoreFile %in% xx$scoreFile)

		return(d)
	}

	# used in decile graph, returns a data frame that is used to plot a theory line on top of decile graph
	hongsTheoryCode <- function(bins, k, h2) {
		theory <- data.frame() # initialize data frame to collect output

		for (zi in 1:bins) {

			thd=-qnorm(k)   # thd = threshold, qnorm is the inverse of pnorm, qnorm(quantile) = how many SD away to find such a quantile
				# here, -qnorm(k) gives -2.326, or 2 and a bit SDs away on the left of the distribution
			h2=h2 # SNP heritability on liability scale, from Dudbridge or user-specified
				# actually R2 on the liability scale

			se=sqrt(1/h2-1) # standard error
			kv=1-pnorm(thd) # what probability of finding a value above thd
			bv=dnorm(thd)   # what density at thd
			iv=bv/kv				# density at thd/probability above thd

			#ckv=0.1
			ckv=1-(1/bins)*zi  # ckv = 1 - (size of each quantile; if decile, 0.1, if quartile, 0.25, if percentile, 0.01)*zi
												# since bins = 30, quantile size is 0.033...
												# = the proportion leftover not counting zi quantiles. e.g. if zi = 3, ckv = 1 - 0.033*3, or 1 - 3 quantiles
			cthd=-qnorm(ckv)	# how many SDs to the left do you have to go to have above proportion to the right of you
			czv=dnorm(cthd)		# what's the density at above SD?
			civ=czv/ckv				# density at cthd / proportion to the right of ckv

			var_l_s=(1-civ*(civ-cthd))*h2+(1-h2)
			p_case_s_int=(thd-civ*h2^.5)/var_l_s^.5 # uses above
			p_case_s11=(1-pnorm(p_case_s_int))*ckv # uses above

			var_l_s_=(1+civ*(-civ+cthd))*h2+(1-h2) # name has an extra _ at the end compared to var_l_s
			p_case_s_int_=(thd+civ*h2^.5)/var_l_s_^.5 # uses above
			p_case_s11_=(1-pnorm(p_case_s_int_))*ckv	# uses above

			p_case_s11_=(1-pnorm((thd-(-qnorm(p_case_s11/ckv)*var_l_s^.5-thd))/var_l_s^.5))*ckv # uses 4above, 6above, seems to overwrite above


			#this is for top and bottom x % ************************************
				p_case_s12=ckv-p_case_s11
				p_case_s12_=ckv-p_case_s11_

				iv2=-iv*kv/(1-kv)
				v103=sqrt(h2*((iv2-iv)^2+(iv2*(iv2-thd)))/(2-h2*iv*(iv-thd)))
				v103=pnorm(v103)           #AUC
				#cat ("aucv: ",v103,"\n")

				estim_OR2=(p_case_s11/p_case_s12)/(p_case_s11_/p_case_s12_)
				#cat("OR2 estimation : ",h2,v103,ckv,p_case_s11/kv,estim_OR2,"\n")  #************
			
			#this is for top x % vs. rest ***************************************
				p_case_s21=kv-p_case_s11
				p_case_s22=(1-ckv)-p_case_s21

				# estim_OR2=(p_case_s11/p_case_s12)/(p_case_s21/p_case_s22)
				#cat("OR2 estimation : ",h2,v103,ckv,p_case_s11/kv,estim_OR2,"\n")  #************

			#this is for top x % vs. pop ***************************************
				# estim_OR2=(p_case_s11/p_case_s12)/(kv/(1-kv))
				# cat("OR2 estimation : ",h2,v103,ckv,p_case_s11/kv,estim_OR2,"\n")  #************
				# estim_OR2 is what you want if comparing to population, which is what robert's "cumulative" plot is doing
				x <- data.frame("ckv" = ckv, "estim_OR2" = estim_OR2)
				theory <- rbind(theory,x)
		} # end for loop

		# fill in last bin
		x <- data.frame("ckv" = 1, "estim_OR2" = 1)	# create last bin
		theory <- rbind(theory,x) 									# add last bin created above
		theory$decile <- bins-theory$ckv*bins + 1 	# convert "ckv" to quantile
		theory[complete.cases(theory),] 						# drop missing values, return the rest
	}

  output$labelh4 <- renderUI({ h4(ifelse(is.null(input$file3) & !is.null(input$file1), '', 'Label file')) })
  tab.input.preview = tabPanel("Input Preview",
                               
                               h4("Risk profile score files"),
                               #uiOutput("inputPreviewUI"),
                               tableOutput("inputPreviewFiles"),
                               br(),
                               h4("First lines of first risk profile score file"),
                               tableOutput("inputPreview"),
                               br(),
                               uiOutput("labelh4"),
                               tableOutput("inputPreviewLabel")
  )
  
  tab.output.summary.b = tabPanel("Output - Summary", 
                                fluidRow(column(11, p("This tab shows various statistics related to the profile score prediction accuracy.")),
                                         column(1, actionButton("outputSummaryBHelp", label = "?"))),
                                uiOutput("outputSummaryBUI"),
                                downloadButton("osdGPRS", "Download GPRS.csv"),
                                dataTableOutput('outputSummary'),
                                br(),
                                br(),
                                p("Shows various statistics produced by the analysis relating to deciles"),
                                downloadButton("osdGPRSdeciles", "Download GPRS_deciles.csv"),
                                dataTableOutput("decileSummary"))
  tab.output.summary.b.r = tabPanel("Output - Summary", 
                                  fluidRow(column(11, p("This tab shows various statistics related to the profile score prediction accuracy.")),
                                             column(1, actionButton("outputSummaryBRHelp", label = "?"))),                                 
                                  p("You can also download these tables as .csv files."),
                                  uiOutput("outputSummaryBRUI"),
                                  downloadButton("osdGPRS", "Download GPRS.csv"),
                                  dataTableOutput('outputSummaryR'),
                                  br(),
                                  br(),
                                  p("Shows various statistics produced by the analysis relating to deciles"),
                                  downloadButton("osdGPRSdeciles", "Download GPRS_deciles.csv"),
                                  dataTableOutput("decileSummary"))
  tab.output.summary.q = tabPanel("Output - Summary", 
                                  fluidRow(column(11, p("This tab shows various statistics related to the profile score prediction accuracy.")),
                                           column(1, actionButton("outputSummaryBQHelp", label = "?"))),  
                                  p("You can also download these tables as .csv files."),
                                  uiOutput("outputSummaryBQUI"),
                                  downloadButton("osdGPRS", "Download GPRS.csv"),
                                  dataTableOutput('outputSummaryQ'),
                                  br())
  tab.output.bar = tabPanel("Barplot",
                            fluidRow(column(11, p("")),
                                     column(1, actionButton("barplotHelp", label = "?"))),
                            uiOutput("barplotHelpUI"),
                            fluidRow(
                              column(3,
                                     #p('This is a clustered bar graph. Select the x (consisting outer and inner clusters) and y axis variables to generate the graph.'),
                                     #p('This graph is modeled on Figures 2 and 4 of Purcell, Wray et al. (Nature 2009,', 
                                    #   a('PMID: 19571811', href='http://www.ncbi.nlm.nih.gov/pubmed/19571811', target='_blank'),
                                    #   '). '
                                    # ),
                                     wellPanel(
                                       #p('For the drop-down boxes, you can also type variable names directly into the box if they are not listed. The graph will still work if such a column exists in the table in Output - Summary'),
                                       
                                       # barY, a dropdown selection box to select Y axis for bar graph
                                       selectizeInput("barY", 
                                                      "Y variable",
                                                      choices = list(  "Nagelkerke's R2 (NK2r)"="NKr2", 
                                                                       "R2 on the liability scale (h2l_r2n)"="h2l_r2n", 
                                                                       "Area Under the Receiver Operator Curve (AUC)"="AUC", 
                                                                       "Odds Ratio comparing 10th and 1st decile (OR10decile)"="OR10decile"
                                                      ),
                                                      selected = "NKr2",
                                                      options = list(create = TRUE)
                                       ),
                                       
                                       # barAcrossClusters, to select what the clusters represent
                                       selectizeInput("barAcrossClusters",
                                                      "Across-cluster variable",
                                                      choices = c("discoverySample", "discoveryTrait", "targetSample", "targetTrait", "pThreshold"),
                                                      selected = "targetSample",
                                                      options = list(create = TRUE)
                                       ),
                                       
                                       # barWithinClusters, to select what's within clusters
                                       selectizeInput("barWithinClusters",
                                                      "Within-cluster variable",
                                                      choices = c("discoverySample", "discoveryTrait", "targetSample", "targetTrait", "pThreshold"),
                                                      selected = "pThreshold",
                                                      options = list(create = TRUE)
                                       ),
                                       
                                       # set title, axis labels
                                       textInput("barMain", "Title"),
                                       uiOutput("barXlabelUI"),
                                       uiOutput("barYlabelUI"),
                                       uiOutput("barLegendTitleUI")
                                     )
                              ),
                              
                              column(9, 
                                     plotOutput("bar")
                              )
                            )
  )
  tab.output.decile = tabPanel("OR deciles",
                               fluidRow(column(11, p("")),
                                        column(1, actionButton("decileHelp", label = "?"))),
                               uiOutput("decileHelpUI"),
                               fluidRow(
                                 column(3,
                                        #p("In this graph, x is score decile and y is odds of being a case in this decile relative to the first decile. Lines represent 95% confidence intervals."),
                                        #p("This graph is modeled on Figure 3 of PGC-SCZ (Nature 2014,", a('PMID: 25056061', href='http://www.ncbi.nlm.nih.gov/pubmed/25056061', target='_blank'), "). This graph should be interpreted recognizing the degree of ascertainment of cases to the Target sample"),
                                        #wellPanel(h5("Refresh Plot"),
                                        #          p("Sometimes the plot does not update to changes in uploaded data. In those cases, click this button."),
                                        #          actionButton("dnRefresh", "Refresh")
                                        #),
                                        uiOutput("decileNewSelectUI"),
                                        uiOutput("decileNewLabelUI"),
                                        uiOutput("decileNewScaleUI"),
                                        uiOutput("decileNewColourUI"),
                                        uiOutput("decileNewSizeUI"),
                                        wellPanel(h5("Theory Line"),
                                                  p("The theory line is based on expectation using the liability threshold model (Hong Lee, unpublished), assuming GPRS are normally distributed. The theory line has 3 input parameters which can be varied to visualize the fit:"),
                                                  checkboxInput("dnTheoryCheck", "Check to compare against theory line"),
                                                  numericInput("dnTheoryBins", "Number of bins (default 10 for decile)", 10),
                                                  numericInput("dnTheoryK", "K for the Target trait (default 0.01 for schizophrenia)", 0.01),
                                                  numericInput("dnTheoryH2", "R2 on the liability scale for Target sample (default 0.08 is a ballpark estimate)", 0.08)
                                        )
                                 ),
                                 column(9,
                                        plotOutput("decileNew"),
                                        wellPanel(
                                          downloadButton("dndPlot", "Download current plot"),
                                          downloadButton("dndData", "Download selected data as .csv"),
                                          br(),
                                          p("If you know how to use ggplot2 in R, download the selected data (above button) and use the following code as a template to perform more in-depth customization."),
                                          uiOutput("decileNewCodeUI")
                                          #br(),
                                          #p("For more information on how to customize graphs made using the ggplot2 package, visit", a(href="http://www.cookbook-r.com/Graphs/", target="_blank", "http://www.cookbook-r.com/Graphs/"))
                                        )
                                 )
                               )
             )
  
  tab.output.decile.r = tabPanel("OR deciles",
                      fluidRow(column(11, p("")),
                               column(1, actionButton("decileHelp", label = "?"))),
                      uiOutput("decileHelpUI"),
                      fluidRow(
                        column(12,
                               plotOutput("decileNew"),
                               wellPanel(
                                 downloadButton("dndPlot", "Download current plot"),
                                 downloadButton("dndData", "Download selected data as .csv")
                                 #br(),
                                 #p("For more information on how to customize graphs made using the ggplot2 package, visit", a(href="http://www.cookbook-r.com/Graphs/", target="_blank", "http://www.cookbook-r.com/Graphs/"))
                               )
                        )
                        )
                      )

  tab.output.cumulative = tabPanel("OR cumulative",
                                   fluidRow(column(11, p("")),
                                            column(1, actionButton("cumulativeHelp", label = "?"))),
                                   uiOutput("cumulativeHelpUI"),
                                   fluidRow(
                                     column(5,
                                            wellPanel(
                                              #p("This graph is provided as an alternative to the decile graph. The decile graph compares odds of being a case in each decile to the odds of the first decile. Here, the odds of being a case are compared to the odds of being a case for the total sample."),
                                              #p("In this graph, GPRS is normalized within the sample and mapped to the x axis, and the y axis is cumulative odds of disease for all individuals with the GPRS score and higher, relative to the odds of disease for the whole sample. This graph is more representative of the utility of GPRS when applied to a population sample"), 
                                              #p("We recommend that the results from top 100 scores are hidden, due to instability caused by using decreasingly small samples to calculate the odds ratio."),
                                              checkboxInput("cumulativeShowTop", "Hide top 100 scores", 1),
                                              # The purpose of this checkbox is to quickly add all profile scores if the intent is to get the cumulative
                                              # OR data, not the plots
                                              checkboxInput("cumulativeAddAll", "Add all profile scores", 0),
                                              uiOutput("fileForCumulativeUI")
                                            ),
                                            downloadButton("cumulData", "Download data as .csv")
                                     ),
                                     column(7,
                                            plotOutput("cumulative")
                                     )
                                   )
  )
  tab.output.auc = tabPanel("AUC plot",
                               fluidRow(
                                 column(12,
                                        plotOutput("aucPlot")
                                 )
                               )
  )
  tab.output.aegp = tabPanel("Output - AEGP", 
                             fluidRow(column(11, p("")),
                                      column(1, actionButton("outputAEGPHelp", label = "?"))),  
                             #p("If the server recognizes you, and you did not run a fresh analysis, there might be files left over from your last session that show up below. If you are not sure, use the 'Cleanup' button at the top and perform a fresh analysis."),
                             uiOutput("outputAEGPUI"),
                             br(),
                             p("This table shows target set genetic variances estimated from profile scores."),
                             br(),
                             downloadButton("aegpdVg2", "Download Vg2FromP.csv"),
                             dataTableOutput("aegpVg2"),
                             br(),
                             br(),
                             p("This table shows estimates of genetic correlation between the discovery and the target set."),
                             br(),
                             downloadButton("aegpdCorr", "Download CorrFromP.csv"),
                             dataTableOutput("aegpCorr")
  )
  tab.wrong.format = tabPanel("Wrong input format!", 
                             textOutput("inputError")
  )
  tab.citations = tabPanel("Citations",
                           p('Dudbridge, F. (2013). "Power and predictive accuracy of polygenic risk scores." ',
                             tags$u('Plos Genetics'),
                             strong('9'),
                             '(3): e1003348. ', 
                             a('PMID: 23555274', href='http://www.ncbi.nlm.nih.gov/pubmed/23555274', target='_blank')
                           ),
                           p('Lee, S. H., M. E. Goddard, et al. (2012). "A Better Coefficient of Determination for Genetic Profile Analysis."',
                             tags$u('Genetic Epidemiology'),
                             strong('36'),
                             '(3): 214-224. ', 
                             a('PMID: 22714935', href='http://www.ncbi.nlm.nih.gov/pubmed/22714935', target='_blank')
                           ),
                           p('Purcell, S. M., N. R. Wray, et al. (2009). "Common polygenic variation contributes to risk of schizophrenia and bipolar disorder."',
                             tags$u('Nature'),
                             strong('460'),
                             '(7256): 748-752. ', 
                             a('PMID: 19571811', href='http://www.ncbi.nlm.nih.gov/pubmed/19571811', target='_blank')
                           ),
                           p('Schizophrenia Working Group of the Psychiatric Genomic Consortium. (2014). "Biological insights from 108 schizophrenia-associated genetic loci."',
                             tags$u('Nature'),
                             strong('511'),
                             '(7510): 421-7. ', 
                             a('PMID: 25056061', href='http://www.ncbi.nlm.nih.gov/pubmed/25056061', target='_blank')
                           )
  )
  tab.gxe = tabPanel("GxE", 
                                  fluidRow(column(11, strong("This tab shows the results of the GxE analysis")),
                                           column(1, actionButton("GxEHelp", label = "?"))),
                                  uiOutput("GxEHelpUI"),
                                  fluidRow(column(4, p("Choose the environmental variables that should be tested:")), column(6, uiOutput("envChecklistUI"))),
                                  downloadButton("gxeDownload", "Download GxE results as csv"),
                                  dataTableOutput('gxeTable'),
                                  plotOutput("gxeHeatmap"),
                                  fluidRow(column(4, uiOutput("GxEselectSCOREUI")), column(4, uiOutput("GxEselectEUI")), column(4, fluidRow(br(), uiOutput("adjustScatterUI")))),
                                  plotOutput("gxeScatter"),
                                  #webGLOutput("gxe3d"),
                                  br())
  
  # select score for gxe scatterplot
  output$GxEselectSCOREUI <- renderUI({
    
    input$GxEaction
    selectizeInput("GxEselectSCORE", "Choose score",
                   choices = processUpload()[[2]],
                   selected = processUpload()[[2]][[1]],
                   width = "100%", multiple = FALSE,
                   options = list(create = TRUE)
    )
  })
  
  # select e for gxe scatterplot
  output$GxEselectEUI <- renderUI({
    
    selectizeInput("GxEselectE", "Choose E",
                   choices = input$envChecklist,
                   selected = input$envChecklist[1],
                   width = "100%", multiple = FALSE,
                   options = list(create = TRUE)
    )
  })
  
  # checkbox for GxE scatterplot
  output$adjustScatterUI <- renderUI({
    checkboxInput('adjustForScatter', 'adjust for covariates', value=FALSE)
  })
  
  # select e for gxe scatterplot
#   output$GxEselectEUI2 <- renderUI({
#     
#     selectizeInput("GxEselectE2", "Choose E2",
#                    choices = input$covariatesChecklist,
#                    selected = input$covariatesChecklist[1],
#                    width = "100%", multiple = FALSE,
#                    options = list(create = TRUE)
#     )
#   })
  
  
  getTabs = reactive({

    
    if(is.null(input$file1) & input$exampleButton == 0) {
      return(list())
    }
    data = processUpload()
    if (is.null(data$quant)) {
      myTabs = list()
    }
    else if(data$check != '') {
      myTabs = list(tab.wrong.format)
    }
    else if(input$GxEaction > 0) {
      myTabs = list(tab.gxe)
    }
    else if(input$action == 0) {
      myTabs = list(tab.input.preview)
    }
    # tabs for quantitative input
    else if(data$quant) {
      myTabs = list(tab.input.preview,
                    tab.output.summary.q,
                    tab.output.aegp,
                    tab.citations  
    )}
    # tabs for reduced binary input
    else if(is.null(input$file3) & !is.null(input$file1)) {
      myTabs = list(tab.input.preview,
                    tab.output.summary.b.r,
                    tab.output.bar,
                    tab.output.decile.r,
                    #tab.output.auc,
                    tab.output.cumulative)
    }
    # tabs for binary input
    else {
      myTabs = list(tab.input.preview,
                    tab.output.summary.b,
                    tab.output.bar,
                    tab.output.decile,
                    tab.output.cumulative,
                    tab.output.aegp,
                    tab.citations)
    }
    myTabs
  })

  output$paneltest <- renderUI({
    myTabs = getTabs()
    myTabs$type='tabs'
    myTabs$id='outputTabs'
    do.call(tabsetPanel, myTabs)
  })
  



	# when UI elements need to change dynamically based on other UI elements, they need to be defined as custom UI elements in the server.r file as opposed to the ui.r file

	# tutorialUI, an HTML-based tutorial for the app, show/hide based on the button 'showTutorial'
	output$tutorialUI <- renderUI({
		if (input$showTutorial%%2 != 0)
			wellPanel(
# 				h3("Tutorial"),
# 				# Introduction
# 					h5("Introduction"),
# 					p("This is a flexible online software-tool that implements two related methods that predict genetic risk and estimate genetic (co)variance from aggregate SNP effects."),
# 					p("The first method is a ", em("Genetic Profile Risk Score"), "(GPRS) regression method, originally implemented in the International Schizophrenia Consortium study (Purcell, Wray et al. 2009). It utilizes identified risk alleles and their effect sizes derived from a genome-wide association study (GWAS), denoted the 'Discovery' sample. These risk alleles and their effect sizes are used to generate GPRSs for individuals in an independent 'Target' sample. The GPRS are then used to predict the disease status or trait level in the Target sample."),
# 					p("The second method, which we call the ", em("Approximate Estimation of Genetic Parameters"), "(AEGP) method estimates genetic parameters attributable to set of risk alleles from the GPRS regression analysis results. The parameters estimated, dependent on the combination of traits in the Discovery and Target samples are the genetic variance (V", tags$sub("g"), ") and genetic correlation (r", tags$sub("g"), ") (Dudbridge 2013), which approximate estimates made directly from SNPs using a GREML approach (Yang, Benyamin et al. 2010; Yang, Lee et al. 2011), sometimes called the SNP-heritability and SNP-correlation, respectively."),
# 					p("The online tool consists of three steps: (1.) Upload Input; (2.) Choose Covariates; (3.) Run GPRS Regression Script."),
# 					p("The required input is a GPRS for each individual; inclusion of covariates is optional. The default input is compatible GPRS files generated by the –score function of PLINK (Purcell, Neale et al. 2007), but other methods, such as Best Linear Unbiased Prediction can be used to generate GPRS and users can upload GRPS files and covariate files in different formats."),
# 					p("The tool takes both disease traits and quantitative traits. We note that if the Discovery sample and the Target sample traits are neither ", em("both"), " disease nor ", em("both"), " quantitative traits, then the AEGP method cannot be applied."),
# 					p("The output consists of a tabular and a graphical part for each method. The specific statistics listed in the table and the types of graphs (e.g., bar graph, decile plot) generated depend on the measurement scale of the traits (i.e. binary or quantitative). Users can customize graphs using the tool, and can download the output table and R-code to customize graphs locally in R."),
# 				# 1. Upload Input
# 					h5("1. Upload Input"),
# 					p("Download the example data", a("here", href="exampleData.zip", target="_blank"), " and unzip to the directory of your choice. It contains 60 'score files'. Also download the accompanying label file", a("here", href="label.csv", target="_blank"), "."),
# 					p("First, we need to upload our score files. Score files are tab-delimited text files containing profile score data generated by PLINK. It has columns 'FID', 'IID', 'SCORE', and 'PHENO' (or 'D', for 'disease', as in the example data) with any additional number of covariate columns. You can upload multiple files at once."),
# 					em("Each of our example files has 1000 individuals with columns FID, IID, SCORE, D, and 5 covariates. Go ahead and upload all 60 example files. The upload process might take a few minutes."),
# 						br(),
# 						br(),
# 					p("Optionally, you can put additional covariates in a separate file with columns 'FID', 'IID', and any number of covariate columns. There should be an equal number of these files as the score files, if you choose to use them."),
# 					em("Our example files do not use separate covariate files."),
# 						br(),
# 						br(),
# 					p("Finally, upload the label file, which is a comma-separated values file. The label file variables provide parameters that define the nature of the GPRS Discovery and Target samples that are used in the GPRS and AEGP scripts. The label file also provides labels for each of the uploaded score files for grouping purposes in graphs. You must construct this file manually."),
# 					# table here
# 					tableOutput("labelExplanation"),
# 					p("Note: To obtain approximate estimates of genetic parameters of variance explained by SNPs or genetic correlation using the AEGP method, the GPRS should be generated from ", em("independent"), " SNPs. To construct a set of quasi-independent (clumped) SNPs to use in construction of the GPRS, use, for example, PLINK with the following thresholds: --clump-p1 1, --clump-p2 1, --clump-r2 0.2, --clump-kb 500 to the GWAS SNPs to select SNPs that contribution to the calculation of the GPRS."),
# 					em("Our example label file is already filled out, so open it in excel to see what the format looks like, then go ahead and upload it. You must make your own label file when using your own data. An easy way to do this (in Windows) is to open command prompt at the directory where all the score files are (look for the option when you shift-right click in File Explorer), enter the command 'dir /w', right click and select 'Mark', select all the score file names, then press Enter to copy. Now you can paste the file names into Excel (use the example label file as a template), and fill out the rest of the columns manually. When you save, make sure to save as a comma-separated values file."),
# 						br(),
# 						br(),
# 				# 2. Choose Covariates
# 					h5("2. Choose Covariates"),
# 					p("After uploading the input files, the app detects any covariate columns in your first uploaded file (this assumes all other files uploaded share the same format, which is true in our example and should be true in yours too), produces a checklist, and allows you to select which ones to include in the upcoming analysis. Look in the 'Input Preview' tab to see an example of what will be fed into the upcoming analysis."),
# 					em("After uploading our example score and label files, you should see five checkboxes named C1 - C5, and the table in the 'Input Preview' tab should change as you check/uncheck each covariate. Try it out, but make sure they are all checked before moving on for our tutorial."),
# 						br(),
# 						br(),
# 				# 3. Run GPRS Regression Script
# 					h5("3. Run GPRS Regression Script"),
# 					p("After uploading input files and selecting covariates to include, click the button to run the regression script. The script takes information from the score and label files, analyzes them, and stores the output in two files called 'GPRS.csv' and 'GPRS_deciles.csv'. The summary and graphs below depend on these files, so if you change anything in the previous fields, click the button again to ensure the changes are reflected in the output. Note that it may take a minute for the analysis to complete, if you have a lot of data."),
# 					em("Go ahead and press the button. A loading message should appear near the top of the window while the app is busy. Our example data should only take a 5-10 seconds to analyze. For reference, we have 60 files with 1000 rows each. After the script finishes, the 'Output - Summary' tab should display 'GPRS.csv' and 'GPRS_deciles.csv' in tabular format."),
# 						br(),
# 						br(),
# 				# 4. Run AEGP Script (optional)
# 					h5("4. Run AEGP Script (optional)"),
# 					p("This script calculates approximate estimates of genetic parameters, as described by Frank Dudbridge (Dudbridge, 2013)."),
# 					p('Dudbridge, F. (2013). "Power and predictive accuracy of polygenic risk scores." ',
# 						tags$u('Plos Genetics'),
# 						strong('9'),
# 						'(3): e1003348. ', 
# 						a('PMID: 23555274', href='http://www.ncbi.nlm.nih.gov/pubmed/23555274', target='_blank')
# 					),
# 					p("Note: The current implantation of the AEGP method assumes GPRS calculated using the PLINK –score function. To obtain approximate estimates of genetic parameters of variance explained by SNPs or genetic correlation using the AEGP method, the GPRS should be generated from ", em("independent"), " SNPs. To construct a set of quasi-independent (clumped) SNPs to use in construction of the GPRS, use, for example, PLINK with the following thresholds: --clump-p1 1, --clump-p2 1, --clump-r2 0.2, --clump-kb 500 to the GWAS SNPs to select SNPs that contribution to the calculation of the GPRS."),
# 					p("The two buttons each run one of two versions of the AEGP script, which takes mainly summary statistics from the label file and 'GPRS.csv', and stores the output in either 'Vg2FromP.csv' or 'CorrFromP.csv' depending on which version is run. You can view these files in a tabular format in the tab 'Output - AEGP'. Currently, none of the graphs depend on the AEGP code."),
# 						br(),
# 				# 5. Output Presentation
# 					h5("5. Output Presentation"),
# 					# 5.1. Tables
# 					strong("5.1. Tables"),
# 					p("After running the GPRS script, click on the 'Output - Summary' tab to view 'GPRS.csv' and 'GPRS_deciles.csv' in tabular formats. You should be able to see labels from your label file as some of the columns. An explanation of various headings in the tables can be found below."),
# 					# table here
# 					tableOutput("summaryExplanation"),
# 						br(),
# 					# 5.2. Bar Graph
# 					strong("5.2. Bar Graph"),
# 					p("To generate bar graphs of h2l_r2, the proportion of variance explained, go to the 'Variance explained' tab. At the left, select how you want to arrange the data, as well as customize the axis labels, etc."),
# 					em("Our example data was the result of PLINK analyses that saw ..."),
# 						br(),
# 					em("... 5 Target samples, with samples 1 and 2 measuring trait 'a', while samples 3 through 5 measured traits 'b', 'c', and 'd' respectively."),
# 						br(),
# 					em("... each had GPRS derived from 2 Discovery samples, both measuring trait 'a'."),
# 						br(),
# 					em("... with each of those 5 x 2 = 10 comparisons using SNPs that pass 6 progressively less stringent p-value thresholds."),
# 						br(),
# 					em("... which brings us to 60 combinations. These 60 PLINK output files are then used as input files for this app."),
# 						br(),
# 						br(),
# 					em("This provides us with some options about how to arrange our bar graph. For example, you can compare the five Target samples as clusters, with the six p-value thresholds within each cluster. This is the default arrangement."),
# 					em("Note that this obviously includes results from both Discovery samples. If you want to restrict the data to one Discovery sample only, you can restrict your upload in the first step to only those files that use the desired Discovery sample."),
# 					em("In doing this, you won't need to adjust your label file, as the algorithm matches labels to files based on the file name. In other words, one label file with information for all your files can suffice for analyses using any subset of those files."),
# 						br(),
# 						br(),
# 					# 5.3. Decile Plot
# 					strong("5.3. Decile Plot"),
# 					p("To generate decile plots similar to Figure 3 of PGC-SCZ (Nature 2013, ", a("PMID 25056061", href="http://www.ncbi.nlm.nih.gov/pubmed/25056061", target="_blank"), "), go to the 'Output - Decile' tab, where each decile based on ranked GPRS is compared to the first decile (lowest risk score). Select the data you want to include in the graphs by setting inclusion criteria, and customize the graph using the options provided."),
# 					em("This app uses special 'selectize' boxes that allows more operations than regular select boxes. Select one or more items simply by clicking on them, and deselect them by clicking on the items, then pressing delete. Depending on the context, you may even be able to create your own items by typing directly into the box."),
# 						br(),
# 						br(),
# 					# 5.4. Cumulative Odds Ratio Plot
# 					strong("5.4. Cumulative Odds Ratio Plot"),
# 					p("To generate a graph similar to decile graphs that compares each individual's (instead of each decile's) odds ratio to the overall odds ratio, go to the 'Output - Cumulative' tab. Additional information is provided in that tab."),
# 					# 5.5. AEGP Tables
# 					strong("5.5. AEGP Tables"),
# 					p("After running one of the two versions of the AEGP script, view the corresponding table in the 'Output - AEGP' tab. The column headings include scoreFile, pThreshold, and pval from GPRS.csv, and also the calculated value (either V", tags$sub("g"), "or r", tags$sub("g"), ") and its 95% confidence interval.")
			  h3("Tutorial"),
			  # Introduction
			  p(em("For the impatient: To get an idea of what to expect from this tool, press the button 'Load examples', and after that 'Run Main Script'.")),
			  h5("Introduction"),
			  

p("GPRS is a flexible online software tool for the simultaneous evaluation of multiple genomic profile risk scores. Using programs like PLINK, SNP effects can be calculated in a discovery set, which can then be used for individual risk prediction in a target set. The target set may consist of a different trait or just a different cohort."),
p("The comparison of the predicted trait values to the measured trait values in the target set allows the estimation of several parameters, such as the accuracy of the profile score, the genetic variance explained by the profile score, the SNP-heritability of the trait in the target set or the genetic correlation between the traits in the discovery and in the target set. The method for estimating the latter two parameters is here referred to as 'Approximate Estimation of Genetic Parameters' (AEGP) and is an integration of the Dudbridge (2013) polygenescore R script."),
			  
h5("Input data"),

p("The main input are profile score files, which in each line contain predicted scores and observed phenotypes for an individual. The profile score input files follow the format of the profile score files produced by PLINK (Purcell, Neale et al. 2007) with the --score option, but the scores themselves can be calculated by other methods, such as Best Linear Unbiased Prediction. Covariates can be supplied in additional columns in the profile score files or in separate covariate files."),
p("The profile score input files require at least the following columns, and have to include a header line:"),
tableOutput("profileScoreFileExplanation"),

p("In addition to the profile score files, a label file is required, which has to contain a header file and the following data:"),
			  
			  tableOutput("labelExplanation"),
			  
p("Note: To obtain approximate estimates of genetic parameters of variance explained by SNPs or genetic correlation using the AEGP method, the GPRS should be generated from ", em("independent"), " SNPs. To construct a set of quasi-independent (clumped) SNPs to use in construction of the GPRS, use, for example, PLINK with the following thresholds: --clump-p1 1, --clump-p2 1, --clump-r2 0.2, --clump-kb 500 to the GWAS SNPs to select SNPs that contribution to the calculation of the GPRS."),
			  p("Note that only 'dtTraitCorrelation' or 'Vg2' has to be provided, the other can be set to NA. The AEGP method can estimate each of these parameters if the other is provided."),
p("If no label file is provided, only some basic analyses can be performed."),

			  br(),
			  h5("Running the analyses"),
			  p("After files have been uploaded, available covariates will appear on the left hand side. After selecting the covariates to be included in the analysis (by default all are selected), the analysis can be started using the button 'Run Main Script'."),
			  
			  br(),

			  h5("Output"),
			  p("As soon as the analysis is finished, several tabs will appear, each of which contains tables or plots which summarise either the prediction accuracy of the profile scores, or parameters estimated from them. For more information on the output tabs, press the corresponding question mark buttons."),
			  br(),

h5("Example files"),
			  p("We have included several example files which can either be downloaded or directly loaded into the app. The files contain randomly (but not uniformly) generated profile scores for three discovery and target set pairs with binary disease phenotypes. In each of the three, the profile scores were calculated using 10 different p-value thresholds."),
			  p("The filenames indicate the discovery set (e.g. 'nocohort1'), the target set (e.g. 'cohort1') and the index of the p-value Threshold (S1 - S10). Further information on the example profile score files can be found in the corresponding label file ('disease1_label.csv').")
	)
})

  # annotation for inputPreview
#   output$inputPreviewUI <- renderUI({
#     if (input$inputPreviewHelp%%2 != 0)
#       wellPanel(
#         h3("Tutorial"),
#         h5("Introduction"),
#         p("This is a flexible online software-tool that implements two related methods that predict genetic risk and estimate genetic (co)variance from aggregate SNP effects."),
#         p(paste0(names(processUpload()$tables[[1]]), collapse=' ' ))
#       )
#   })
  
  # annotation for outputSummaryB
  output$outputSummaryBUI <- renderUI({
    df = getVarDescription(c('scoreFile', 'pThreshold', 'NKr2', 'pval', 'h2l_r2n', 'se_h2l_r2', 'AUC', 'OR10decile', 'ORL95', 'ORH95'))
    if (input$outputSummaryBHelp%%2 != 0)
      wellPanel(
        renderTable(df)
      )
  })
  
  # annotation for outputSummaryBR
  output$outputSummaryBRUI <- renderUI({
    df = getVarDescription(c('scoreFile', 'pThreshold', 'NKr2', 'pval', 'h2l_r2n', 'se_h2l_r2', 'AUC', 'OR10decile', 'ORL95', 'ORH95'))
    if (input$outputSummaryBRHelp%%2 != 0)
      wellPanel(
        renderTable(df)
      )
  })
  
  # annotation for outputSummaryQ
  output$outputSummaryQUI <- renderUI({
    df = getVarDescription(c('scoreFile', 'pThreshold', 'NKr2', 'pval', 'h2l_r2n', 'se_h2l_r2', 'AUC', 'OR10decile', 'ORL95', 'ORH95'))
    if (input$outputSummaryQHelp%%2 != 0)
      wellPanel(
        renderTable(df)
      )
  })

  # annotation for outputSummaryQ
  output$barplotHelpUI <- renderUI({
    if (input$barplotHelp%%2 != 0)
      wellPanel(
        p('This is a clustered bar graph of the data in the first table of the "Output - Summary" tab.'),
        p('The quantity displayed and the clustering along the x-axis can be changed in the drop-down menus on the left hand side.'),
        p('This graph is modeled after Figures 2 and 4 of Purcell, Wray et al. (Nature 2009,', 
          a('PMID: 19571811', href='http://www.ncbi.nlm.nih.gov/pubmed/19571811', target='_blank'),
          '). '
        )
      )
  })

  # annotation for outputSummaryQ
  output$decileHelpUI <- renderUI({
    if (input$decileHelp%%2 != 0)
      wellPanel(
        p("In this graph, x is score decile and y is odds of being a case in this decile relative to the first decile. Lines represent 95% confidence intervals."),
        p('The graph is based on the data shown in the second table of the "Output - Summary" tab.'),
        p("This graph is modeled after Figure 3 of PGC-SCZ (Nature 2014,", a('PMID: 25056061', href='http://www.ncbi.nlm.nih.gov/pubmed/25056061', target='_blank'), "). This graph should be interpreted recognizing the degree of ascertainment of cases to the Target sample.")
      )
  })

  # annotation for cumulative tab
  output$cumulativeHelpUI <- renderUI({
    if (input$cumulativeHelp%%2 != 0)
      wellPanel(
        p("This graph is provided as an alternative to the decile graph. The decile graph compares odds of being a case in each decile to the odds of the first decile. Here, the odds of being a case are compared to the odds of being a case for the total sample."),
        p("In this graph, GPRS is normalized within the sample and mapped to the x axis, and the y axis is cumulative odds of disease for all individuals with the GPRS score and higher, relative to the odds of disease for the whole sample. This graph is more representative of the utility of GPRS when applied to a population sample"), 
        p("We recommend that the results from top 100 scores are hidden, due to instability caused by using decreasingly small samples to calculate the odds ratio.")
      )
  })

  # annotation for GxE
  output$GxEHelpUI <- renderUI({
    if (input$GxEHelp%%2 != 0)
      wellPanel(
        p("Interactions between the profile scores and the selected covariates are tested."),
        p("For each environmental variable, the following model is fitted in a linear model for quantitative phenotypes and in a logistic model for binary phenotypes:"), 
        #p("PHENO ~ SCORE + C_i + E + SCORE*C_i + SCORE*E + C_i*E"),
        #p("where PHENO is the phenotype, SCORE the profile score, E the environmental variable for which an interaction is tested and C_i the ith covariate"),
        p("PHENO ~ SCORE.adj + E.adj + SCORE.adj*E.adj"),
        p("where SCORE.adj is the SCORE adjusted for all selected covariates and E.adj is E adjusted for all selected covariates."),
        p("The numbers represent uncorrected p-values of the estimate of the GxE interaction term."),
        p("Below the heatmap, a scatterplot of a given profile score file and environmental variable is plotted for binary phenotypes. This does, however, not take other covariates into account.")
      )
  })
  
  # annotation for outputAEGP
  output$outputAEGPUI <- renderUI({
    df = getVarDescription(c('scoreFile', 'pThreshold', 'pval', 'vg', 'vgLo', 'vgHi', 'rG', 'rGLo', 'rGHi'))
    if (input$outputAEGPHelp%%2 != 0)
      wellPanel(
        #p("Only profile scores with a p-value Threshold of 1 are being used here to get an unbiased estimate."),
        p("The AEGP method uses the R2 p-value (3rd column) and several parameters from the label file to estimate the genetic variance explained by all SNPs (vg)
          or the genetic correlation between the traits in the discovery and target set (rG). The parameter which should be estimated can be set to NA in the label file, the other parameters have to be present."),
        p("The estimated genetic variances and correlations are quite sensitive to some of the input parameters, in particular Vg1, the SNP-heritability in the discovery set. The 95% confidence intervals assume that input parameters are all correct."),
        renderTable(df)
      )
  })
  
	# labelExplanation, a table to explain headings in the label file
	output$labelExplanation <- renderTable({
		df = read.table("labelExplanation.txt", sep='\t', header=TRUE, quote="\"")
    names(df) = c('Column name', 'Explanation')
    df
	})

output$profileScoreFileExplanation = renderTable({
  df = as.data.frame(matrix(c('FID',  	'family ID',
                            'IID',		'individual ID',
                            'PHENO',	'measured phenotype',
                            'SCORE',	'predicted score'), 4, byrow=T))
  names(df) = c('Column name', 'Explanation')
  df
})

	# table for summaryExplanation
	output$summaryExplanation <- renderTable({
		read.table("summaryExplanation.txt", sep='\t', header=TRUE, quote="\"")
	})


  output$sf <- renderUI({
    if(is.null(input$file1) & input$exampleButton == 0) {
      return(NULL)
    }
    data = processUpload()
    txt = ''
    if(input$exampleButton > 0 & is.null(input$file1)) txt = 'Example files loaded.'
    if(input$exampleButton == 0 & is.null(input$file3)) txt = span('Reduced output - no label file.', style = "color:red")
    h4(txt)
    #h4(paste('Showing output for', ifelse(data$quant, 'quantitative', 'binary'), 'phenotypes.'), span(ifelse(is.null(input$file3) & input$exampleButton == 0, 'Reduced output - no label file.' ,''), style = "color:red"))
    
  })
  
  

	# covariatesChecklist, a dynamic UI checkbox group to select columns
	output$covariatesChecklistUI <- renderUI({
		# don't run until at least score and label files are uploaded
		if (is.null(input$file1) & input$exampleButton == 0) return(NULL)

		# process the upload and get data in the right format for originalScript
		data <- processUpload()
    

		table <- data$tables[[1]]
		columnNames <- names(table) # names of columns
		options <- columnNames[columnNames!="FID" & 
														columnNames!= "IID" & 
														columnNames!= "PHENO" & 
														columnNames!= "SCORE" &
														columnNames!= "COUNT"
													]

		# checkbox group for covariate columns, returns a vector with separate checked items
		checkboxGroupInput("covariatesChecklist", 
											"", 
											choices = options,
											selected = options,
                      inline=TRUE
		)
	})


  # envChecklist, a dynamic UI checkbox group to select env variables
  output$envChecklistUI <- renderUI({
    # don't run until at least score and label files are uploaded
    if (is.null(input$file1) & input$exampleButton == 0) return(NULL)
    
    # process the upload and get data in the right format for originalScript
    data <- processUpload()

    table <- data$tables[[1]]
    columnNames <- names(table) # names of columns
    
    options <- columnNames[columnNames!="FID" & 
                             columnNames!= "IID" & 
                             columnNames!= "PHENO" & 
                             columnNames!= "SCORE" &
                             columnNames!= "COUNT"
                           ]
    if(is.null(input$envChecklist))
      sel = options
    else sel <- input$envChecklist
    
    checkboxGroupInput("envChecklist", 
                       "", 
                       choices = options,
                       selected = sel,
                       inline=TRUE
    )
  })






	# run, a reactive observer, watches button, runs (isolate) the following:
	# processUpload, originalScript, and processLabelFile
	run <- observe({
    # deactivate button, if files not uploaded
	  if(is.null(input$file1) & input$exampleButton == 0) {
      disableActionButton("action",session)
      disableActionButton("GxEaction",session)
      return(NULL) } else {
	    enableActionButton("action",session)
	    enableActionButton("GxEaction",session)
	  }
	  if(!is.null(input$file1) | input$exampleButton > 0) {
	    disableActionButton("exampleButton",session)
	  }
    
		# depends on action button, but not at the default value
		if (input$action==0) return(NULL)

		isolate({
		  disableActionButton("action",session)
			# select Output - Summary tab
			updateTabsetPanel(session, "outputTabs", selected = "Output - Summary")
			data <- processUpload()
      
      if(is.null(data$targetK)) data$targetK = 0
			# run the calculation! originalScript() takes 6 parameters and produces GPRS.csv
			originalScriptBQ(numFiles = length(data$scoreFilenames),
										 scoreFilenames = data$scoreFilenames, 
										 covariatesFilenames = data$covariatesFilenames, 
										 dataTables = data$tables,
										 covariatesSelected = input$covariatesChecklist,
										 targetK = data$targetK
			)

			# process label file if uploaded
			if (!is.null(input$file3) | input$exampleButton > 0){
				processLabelFile()
			}
      
			df <- read.csv("GPRS.csv", header=T, sep=",")
			# select files where pThreshold is 1
			#df = df[df$pThreshold == 1, ]
			# select files where discovery and target are on the same scale
			df = df[(df$discoveryK == 0 & df$targetK == 0) | (df$discoveryK != 0 & df$targetK != 0), ]
			
			if(nrow(df) == 0) {enableActionButton("action",session); return(NULL)}
      
			runVg2FromP(df)
			runCorrFromP(df)
      
			enableActionButton("action",session)
		})
	})

runGxE <- reactive({
  
  # depends on GxEaction button, but not at the default value
  if (input$GxEaction==0) return(NULL)
  isolate({ 
    disableActionButton("GxEaction",session)
    data = processUpload()
    pvals = data.frame(matrix(NA, 0, 3))
    if(is.null(input$envChecklist))
      envvars = input$covariatesChecklist
    else
      envvars = input$envChecklist
    
    # choose gxe method
#     if(input$gxeMethod)
#       fun = gxeFor1File
#     else fun = gxeFor1File2
    fun = gxeFor1File2
      
    for(i in 1:length(data$scoreFilenames)) {
    #for(i in 1:10) {
      pvals = rbind(pvals, data.frame(data$scoreFilenames[i], fun(data[[1]][[i]], input$covariatesChecklist, envvars, data$quant)))
    }
    names(pvals) = c('scoreFilename', 'E', 'p')
    enableActionButton("GxEaction",session)
  })
  pvals
})

# analyses GxE in one profile score file, for each E in covars
# returns data.frame of interaction p-values (one for each E in covars)
gxeFor1File = function(dat, covars, envvars, quant=F) {
  
  #browser()
  #covars = setdiff(covars, envvars)
  pvals = c()
  for(co in 1:length(envvars)) {
    allcov = covars
    thiscov = envvars[co]
    othercov = setdiff(allcov, thiscov)
    # controls for covariates properly
    formula1 = as.formula(paste0('PHENO-1 ~ SCORE + ', thiscov, ' + ', paste0(othercov, collapse=' + '), ' + ', thiscov, '*SCORE + ', paste0(paste0(othercov, '*SCORE'), collapse=' + '), ' + ', paste0(paste0(othercov, '*', thiscov), collapse=' + ')))
    # does not control for covariates properly
    # formula2 = as.formula(paste0('PHENO-1 ~ SCORE + ', paste0(allcov, collapse=' + '), ' + SCORE*', thiscov))
    gxeterm = paste0('SCORE:', thiscov)
    model = glm( formula1, data=dat, family=if(quant) gaussian else binomial(logit))
    res = summary(model)$coefficients
    pvals[thiscov] = res[rownames(res)==gxeterm, 4]
    
  }
  data.frame(names(pvals), pvals)
  
}

# analyses GxE in one profile score file, for each E in covars
# returns data.frame of interaction p-values (one for each E in covars)
# adjusts SCORE and E for covariates
gxeFor1File2 = function(dat, covars, envvars, quant=F) {
  
  #browser()
  #covars = setdiff(covars, envvars)
  pvals = c()
  for(co in 1:length(envvars)) {
    allcov = covars
    thiscov = envvars[co]
    othercov = setdiff(allcov, thiscov)
    d2 = dat
    d2[,thiscov] = lm(as.formula(paste0(thiscov, ' ~ ', paste0(othercov, collapse=' + ') )), d2)$residuals
    d2[,'SCORE'] = lm(as.formula(paste0('SCORE ~ ', paste0(othercov, collapse=' + ') )), d2)$residuals
    formula1 = as.formula(paste0('PHENO-1 ~ SCORE + ', thiscov, ' + SCORE*', thiscov))

    gxeterm = paste0('SCORE:', thiscov)
    model = glm( formula1, data=d2, family=if(quant) gaussian else binomial(logit))
    res = summary(model)$coefficients
    pvals[thiscov] = res[rownames(res)==gxeterm, 4]
    
  }
  data.frame(names(pvals), pvals)
  
}


# reads GPRS.csv, writes Vg2FromP.csv
runVg2FromP = function(df) {
  # initialize dataframe and output file
  Output=data.frame("scoreFile",	 # identifier
                    "pThreshold",	 # i.e. pupper, from GPRS.csv
                    "pval",				 # i.e. p, or pval from GPRS.csv
                    "vg",					 # output from function estimateVg2FromP
                    "vgLo",				 # .., lower boundary for 95% CI
                    "vgHi"#,				 # .., higher boundary for 95% CI
                    #"theoreticalP" # theoretical p-value
  )
  write.table(Output,"Vg2FromP.csv",row.names=F,col.names=F,sep=",")
  
  # aegpVg2 for loop
  for(i in 1:nrow(df)){
    x <- df[i,]	# one row at a time
    binary = !is.na(x$AUC) # determine if this file has a binary phenotype
    
    TraitD_ncase=x$discoveryCaseN 		# number of cases in discovery sample, from label file
    TraitD_ncont=x$discoveryControlN 	# number of cases in target sample, from label file
    TraitT_ncase=x$Propcase*x$N 			# number of cases in target sample, from GPRS.csv
    TraitT_ncont=x$N-TraitT_ncase 		# number of controls in target sample, from GPRS.csv
    
    n1=TraitD_ncase+TraitD_ncont 			# N in discovery sample. "1" refers to discovery
    n2=TraitT_ncase+TraitT_ncont 			# N in target sample. "2" refers to target
    
    sampling1=TraitD_ncase/n1 				# proportion of cases
    sampling2=TraitT_ncase/n2 				# proportion of cases
    
    prevalence1=x$discoveryK 					# from label file
    prevalence2=x$targetK 						# from label file
    
    corr=x$dtTraitCorrelation 				
    # correlation between genetic effect sizes in the two populations, from label file. 
    # also known as rG (r = correlation coefficient, G = genetic effects)
    
    nsnp=x$nIndependentSNPs
    # number of independent SNPs used to create SCORE, extract from label file for each threshold 
    # emphasize in tutorial that SNPs have to be independent
    p=x$pval
    # p value from polygenic score model, value of "pval" column in GPRS.csv
    plower=0
    # lower bound on p-value for selection from discovery sample, value defaults to 0
    pupper=x$pThreshold
    # upper bound on p-value for selection from discovery sample, value of "pThreshold" (from label file) column for the same entry
    
    # estimate Vg 2 FromP = estimate genetic variance in the target sample from pval, requires source("polygenescore.R")
    Vg2FromP <- estimateVg2FromP( n1=n1, n2=n2, 
                                  sampling1=sampling1, sampling2=sampling2,
                                  prevalence1=prevalence1, prevalence2=prevalence2,
                                  corr=corr, 
                                  vg1=0, # problem here? ask anna
                                  nsnp=nsnp, p=p, plower=plower, pupper=pupper,
                                  weighted=T, binary=binary, # RIGHT HERE 
                                  lambdaS1=NA,lambdaS2=NA, nullfraction=0, shrinkage=F, logrisk=F
    )
    if(is.na(corr)) Vg2FromP[1:3] = NA 
    Vg2FromP ## output corresponds with Dudbridge et al. (estimate + upper/lower CI)
    
#     vg1 = x$Vg1
#     vg2 = x$Vg2
#     theop = polygenescore(n1=n1, n2=n2,
#                   sampling1=sampling1, sampling2=sampling2,
#                   prevalence1=prevalence1, prevalence2=prevalence2,
#                   corr=corr,
#                   vg1=vg1, vg2=vg2,
#                   nsnp=nsnp, plower=plower, pupper=pupper,
#                   weighted=T, binary=binary,
#                   lambdaS1=NA,lambdaS2=NA,nullfraction=0,shrinkage=F,logrisk=F)$p
    
    # output
    Output=data.frame(x$scoreFile, pupper, p, Vg2FromP$vg, Vg2FromP$vgLo, Vg2FromP$vgHi)
    write.table(Output,"Vg2FromP.csv",row.names=F,col.names=F,sep=",",append=T)
  } # end aegpVg2 for loop
}

# reads GPRS.csv, writes CorrFromP.csv
runCorrFromP = function(df) {
  
  # initialize dataframe and output file
  Output=data.frame("scoreFile",	# identifier
                    "pThreshold",	# i.e. pupper, from GPRS.csv
                    "pval",				# i.e. p, or pval from GPRS.csv
                    "rG",					# output from function estimateCorrFromP
                    "rGLo",				# .., lower boundary for 95% CI
                    "rGHi"				# .., higher boundary for 95% CI
  )
  write.table(Output,"CorrFromP.csv",row.names=F,col.names=F,sep=",")
  
  # aegpCorr for loop
  for(i in 1:nrow(df)){
    x <- df[i,]	# one row at a time
    binary = !is.na(x$AUC) # determine if this file has a binary phenotype
    
    
    TraitD_ncase=x$discoveryCaseN 		# number of cases in discovery sample, from label file
    TraitD_ncont=x$discoveryControlN 	# number of cases in target sample, from label file
    TraitT_ncase=x$Propcase*x$N 			# number of cases in target sample, from GPRS.csv
    TraitT_ncont=x$N-TraitT_ncase 		# number of controls in target sample, from GPRS.csv
    
    n1=TraitD_ncase+TraitD_ncont 			# N in discovery sample. "1" refers to discovery
    n2=TraitT_ncase+TraitT_ncont 			# N in target sample. "2" refers to target
    
    sampling1=TraitD_ncase/n1 				# proportion of cases
    sampling2=TraitT_ncase/n2 				# proportion of cases
    
    prevalence1=x$discoveryK 					# from label file
    prevalence2=x$targetK 						# from label file
    
    Vg1=x$Vg1         								# from label file
    Vg2=x$Vg2					      					# from label file
    
    nsnp=x$nIndependentSNPs
    # number of independent SNPs used to create SCORE, extract from label file for each threshold 
    # emphasize in tutorial that SNPs have to be independent
    p=x$pval
    # p value from polygenic score model, value of "pval" column in GPRS.csv
    plower=0
    # lower bound on p-value for selection from discovery sample, value defaults to 0
    pupper=x$pThreshold
    # upper bound on p-value for selection from discovery sample, value of "pThreshold" (from label file) column for the same entry
    
    # estimate Vg 2 FromP = estimate genetic variance in the target sample from pval, requires source("polygenescore.R")
    CorrFromP <- estimateCorrFromP(n1=n1, n2=n2, 
                                   sampling1=sampling1, sampling2=sampling2,
                                   prevalence1=prevalence1, prevalence2=prevalence2,
                                   vg1=Vg1, vg2=Vg2,
                                   nsnp=nsnp, p=p, plower=plower, pupper=pupper,
                                   weighted=T, binary=binary 
    )
    if(is.na(Vg2))  CorrFromP[1:3] = NA 
    CorrFromP ## output corresponds with Dudbridge et al. (estimate + upper/lower CI)
    
    # output
    Output=data.frame(x$scoreFile, pupper, p, CorrFromP$corr, CorrFromP$corrLo, CorrFromP$corrHi)
    write.table(Output,"CorrFromP.csv",row.names=F,col.names=F,sep=",",append=T)
  } # end aegpCorr for loop
}



	# cleanup, a reactive observer, watches button cleanupButton, deletes GPRS.csv, GPRS_deciles.csv, Vg2FromP.csv, CorrFromP.csv, www/plot.pdf
# 	cleanup <- observe({
# 		if(input$cleanupButton>0) {
#       clean()
# 		}
# 	})

  processExample <- reactive({
    if(input$exampleButton > 0) {
      disableActionButton("exampleButton",session)
      processExampleFiles()
    }
  })

  output$inputPreviewFiles <- renderTable(
  {
    # don't run until at least score and label files are uploaded
    if (is.null(input$file1) & input$exampleButton == 0) return(NULL)
    # process the upload
    data <- processUpload()
    
    mat = cbind(data$scoreFilenames, t(sapply(data$tables, dim)))
    colnames(mat) = c('filename', 'rows', 'columns')
    as.data.frame(mat)
  }, 
  options = list(bPaginate = FALSE, bFilter = FALSE)
  )

  output$inputPreviewLabel <- renderTable(
  {
    # don't run until at least score and label files are uploaded
    if (is.null(input$file3) & input$exampleButton == 0) return(NULL)
    # process the upload
    getLabelFile()
  }, 
  options = list(bPaginate = FALSE, bFilter = FALSE)
  )

	# tab inputPreview
		output$inputPreview <- renderTable(
			{
				# don't run until at least score and label files are uploaded
				if (is.null(input$file1) & input$exampleButton == 0) return(NULL)
				
				# depend on button so loading bar shows up
				input$action

				# process the upload
				data <- processUpload()
       
				# apply the column filter from "data" and "input$covariatesChecklist"
					# data is a list with 4 items: $tables (list of dataframes), $scoreFilenames (vector), $covariatesFilenames (vector), $targetK (vector)
					# covariatesChecklist is a vector containing column names of the columns selected
				defaultColumns <- names(data$tables[[1]])[names(data$tables[[1]]) =="FID" | 
																								names(data$tables[[1]])=="IID" | 
																								names(data$tables[[1]])=="PHENO" | 
																								names(data$tables[[1]])=="SCORE"
																							]
				for (i in 1:length(data$scoreFilenames)) {
					# for each table, i.e. for each uploaded file
					data$tables[[i]] <- data$tables[[i]][c(defaultColumns, input$covariatesChecklist)]
				}

				preview <- head(data$tables[[1]])
			}, 
			options = list(bPaginate = FALSE, bFilter = FALSE)
		)

	# tab outputSummary
		# download for GPRS.csv
		output$osdGPRS <- downloadHandler(
			filename = function() {
				"GPRS.csv"
			},
			content = function(file) {
				file.copy("GPRS.csv", file, overwrite=TRUE)
			}
		)

		# download for GPRS_deciles.csv
		output$osdGPRSdeciles <- downloadHandler(
			filename = function() {
				"GPRS_deciles.csv"
			},
			content = function(file) {
				file.copy("GPRS_deciles.csv", file, overwrite=TRUE)
			}
		)

# download for GPRS.csv
output$gxeDownload <- downloadHandler(
  filename = 'GxE_pvalues.csv',
  content = function(file) {
    write.csv(dcast(runGxE(), scoreFilename ~ E, value.var='p'), file)
  }
)

getFormattedGPRS = function() {
  x <- read.csv("GPRS.csv", header = T, sep = ",")
  for(col in c("NKr2", "h2l_r2n","se_h2l_r2",  "AUC",  "OR10decile","ORL95", "ORH95", 'pval_r2', 'R2v', 'vr')) {
    x[,col] = sprintf("%.2f", x[,col])
  }
  x$pval = prettyNum(x$pval, digits=2, scientific=T)
  x
}

		# table for GPRS.csv
		output$outputSummary <- renderDataTable(
			{
				# re-run if button pressed
				input$action

				# render "GPRS.csv" to the datatable widget
				x <- getFormattedGPRS()
				x = x[, names(x) %in% c('scoreFile', 'pThreshold', "NKr2", "pval", "h2l_r2n","se_h2l_r2",  "AUC",  "OR10decile","ORL95", "ORH95")]
				#x$pval = prettyNum(x$pval, digits=2, scientific=T); for(i in c(3,5:10)) {x[,i] = sprintf("%.2f", x[,i])}
        x
			},
			options = list(bPaginate = FALSE, bFilter = FALSE)
		)

  # table for GPRS.csv for Q data
  output$outputSummaryQ <- renderDataTable(
  {
    # re-run if button pressed
    input$action
    
    # render "GPRS.csv" to the datatable widget
    x <- read.csv("GPRS.csv", header = T, sep = ",")
    x[, names(x) %in% c('scoreFile', 'pThreshold', 'pval_r2', 'R2v', 'vr')]
    
  },
  options = list(bPaginate = FALSE, bFilter = FALSE)
  )
  # table for GPRS.csv for Q data
  output$outputSummaryR <- renderDataTable(
  {
    # re-run if button pressed
    input$action
    
    # render "GPRS.csv" to the datatable widget
    x <- getFormattedGPRS()
    x[, names(x) %in% c('scoreFile', "AUC",  "OR10decile","ORL95", "ORH95", "NKr2", "pval")]
    
  },
  options = list(bPaginate = FALSE, bFilter = FALSE)
  )

  # table for GPRS.csv
  output$gxeTable <- renderDataTable(
  {
    prettyNum(dcast(runGxE(), scoreFilename ~ E, value.var='p'), digits=2, scientific=T)
  },
  options = list(bPaginate = FALSE, bFilter = FALSE)
  )

		# table for GPRS_deciles.csv
		output$decileSummary <- renderDataTable(
			{
				# re-run if button pressed
				input$action

				# render "GPRS_deciles.csv" to the datatable widget
				x = read.csv("GPRS_deciles.csv", header=TRUE, sep=",")
        for(i in 3:5) {
          x[,i] = sprintf("%.2f", x[,i])
        }
        x
			},
			options = list(bPaginate = FALSE, bFilter = FALSE)
		)

	# tab bar
		# some UI elements
			output$barXlabelUI <- renderUI({
				textInput("barXlabel", "X Axis Label", value = input$barAcrossClusters)
			})
			output$barYlabelUI <- renderUI({
				textInput("barYlabel", "Y Axis Label", value = input$barY)
			})
			output$barLegendTitleUI <- renderUI({
				textInput("barLegendTitle", "Legend Title", value = input$barWithinClusters)
			})

		# bar plot
			output$bar <- renderPlot({
				# read data
				x <- read.csv("GPRS.csv", header = T, sep = ",")

        if(is.null(input$file3) & !is.null(input$file1)) {
          ggplot(data = x, 
                 environment = environment(), # required to make get() work later
                 aes( factor(scoreFile), get(input$barY) )
          ) + geom_bar(  aes(fill = factor(scoreFile)), 
                         position = "dodge", 
                         stat = "identity"
          ) + ggtitle(input$barMain) + xlab('') + ylab(input$barYlabel) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
            
          
        } else {
				# ggplot(data, environment, aes(x,y), position, stat)
				ggplot(data = x, 
								environment = environment(), # required to make get() work later
								aes( factor(get(input$barAcrossClusters)), get(input$barY) )
				) +

				# + title, x and y labels
				ggtitle(input$barMain) + xlab(input$barXlabel) + ylab(input$barYlabel) +

				# + geom_bar(aes(fill), position, stat)
				geom_bar(	aes(fill = factor(get(input$barWithinClusters))), 
									position = "dodge", 
									stat = "identity"
				) +

				# + scale_fill_discrete(name) legend title
				scale_fill_discrete(name = input$barLegendTitle) +

				# + geom text 
				geom_text(aes(y=get(input$barY)+0.01, 
											label=factor(get(input$barWithinClusters)), 
											group=factor(get(input$barWithinClusters))
									), 
									position = position_dodge(width=0.9), 
									size=4
				)
        }
			})

	# tab decile new
		# decileNewSelectUI, a dynamic UI dropdown selection box to select data source for new decile graph
		output$decileNewSelectUI <- renderUI({
			input$action
			select5("dns1", "dns2", "dns3", "dns4", "dns5")
		})

		# decileNewLabelUI, a dynamic UI wellPanel that customizes the chart, axes, and legend labels
		output$decileNewLabelUI <- renderUI({
			input$action
			y <- read.csv("GPRS_deciles.csv", header = T, sep = ",")
			d <- selectDataFrom5(y, input$dns1, input$dns2, input$dns3, input$dns4, input$dns5)

			wellPanel(
				h5("Labels & Legends"),
				textInput("decileNewTitle", "Enter title"),
				textInput("decileNewXlabel", "Enter x axis label", "Decile"),
				textInput("decileNewYlabel", "Enter y axis label", "OR"),
				hr(),
				textInput("decileNewLegendTitle", "Enter legend title", names(d)[1]),
				selectizeInput("decileNewLegendLabels", 
					"Choose labels from the list. Note: order matters, so make sure you choose the right label for the right data. You can type into the box to add custom labels to the list.",
					choices = d$scoreFile,
					selected = d$scoreFile,
					width = "100%", multiple = TRUE,
					options = list(create = TRUE)
				),
				radioButtons("decileNewLegendLocation", 
					label = "Legend location",
					choices = c("Inside", "Outside"), 
					selected = "Inside"
				)
			)
		})

		# decileNewScaleUI, a dynamic UI that customizes axis scales, ticks, and guidelines
		output$decileNewScaleUI <- renderUI({
			wellPanel(
				h5("Scales"),
				checkboxInput("dnYlogCheck", "Log scale on y-axis"),
				selectizeInput("dnYticks", 
					"Enter positions of y-axis tick marks", 
					choices = "0",
					multiple=TRUE, options = list(create=TRUE)
				),
				hr(),
				radioButtons("dnGuidelines",
					"Guidelines",
					choices = c("All", "Major Only", "Off"),
					selected = "All"
				),
				hr(),
				checkboxInput("dnAxisLines", "Axis Lines")
			)
		})

		# decileNewColourUI, change the colour scheme, bg colour
		output$decileNewColourUI <- renderUI({
			wellPanel(
				h5("Colours"),
				radioButtons("dnColour",
					"Either select one of the following colour schemes...",
					choices = c("Default ggplot2", "Colour-blind-friendly"),
					selected = "Default ggplot2"
				),
				selectizeInput("dnColourCustom",
					"...Or enter colours in hex code to create your own colour scheme",
					choices = "#000000",
					selected = NULL,
					multiple=TRUE, options = list(create=TRUE)
				),
				hr(),
				textInput("dnbgFill"	,	"Enter Background Fill Colour in Hex Code"		, "#"),
				textInput("dnbgColour",	"Enter Background Outline Colour in Hex Code"	, "#")
			)
		})

		# decileNewSizeUI, change the plot dimensions
		output$decileNewSizeUI <- renderUI({
			wellPanel(
				h5("Size (Only Affects Download)"),
				numericInput("dnSizeX", "Enter Width (cm), enter 0 for as-displayed"	, 0),
				numericInput("dnSizeY", "Enter Height (cm), enter 0 for as-displayed"	, 0)
			)
		})

		# dnPlotString, a reactive expression to create a string containing a ggplot() command
		dnPlotString <- reactive({
			# create string variable to store graph command
			g <- "ggplot(d)"
			# base
				g <- paste(g, "geom_pointrange(aes(decile, estimate, ymin=lowerCI, ymax=upperCI, colour=scoreFile), position=position_dodge(width=0.3))", sep="+")
			# labels
				g <- paste(g, "ggtitle(input$decileNewTitle)", "xlab(input$decileNewXlabel)",	"ylab(input$decileNewYlabel)", sep="+")
			# legend - labels
				g <- paste(g, "scale_colour_discrete(name = input$decileNewLegendTitle, labels = input$decileNewLegendLabels)", sep="+")
			# legend - location
# 				if(input$decileNewLegendLocation=="Inside") {
# 					g <- paste(g, "theme(legend.justification=c(0,1), legend.position=c(0,1))", sep="+")
# 				}
			# scale - y log
				if(input$dnYlogCheck){
					g <- paste(g, 'coord_trans(y="log2")', sep='+')
				}
			# scale - y ticks
				if(is.null(input$dnYticks)) {
				} else {
					g <- paste(g, 'scale_y_continuous(breaks=as.numeric(input$dnYticks))', sep='+')
				}
			# scale - guidelines
				if(input$dnGuidelines=="All") {
				}
				if(input$dnGuidelines=="Major Only") {
					g <- paste(g, 'theme(panel.grid.minor=element_blank())', sep='+')
				}
				if(input$dnGuidelines=="Off") {
					g <- paste(g, 'theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank())', sep='+')
				}
			# scale - axis lines
				if(input$dnAxisLines) {
					g <- paste(g, 'theme(axis.line = element_line(colour = "black"))', sep='+')
				}
			# colour - schemes
				if(!is.null(input$dnColourCustom)){
					palette <- input$dnColourCustom
					# super fucking complicated 4-layer-deep paste nest
					g <- paste(g, paste0(	'scale_colour_manual(values=', 
																paste0(	"c(", 
																				paste0(	"'", 
																								palette, 
																								"'", collapse=","
																				), 
																				")"
																), 
																')'
												), sep='+')
				} else {
					if(input$dnColour=="Colour-blind-friendly"){
						palette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
						g <- paste(g,paste0('scale_colour_manual(values=',paste0("c(",paste0("'",palette,"'",collapse=","),")"),')'),sep='+')
					} 
				}
			# colour - background
				if(!is.null(input$dnbgFill) & !is.null(input$dnbgColour) & input$dnbgFill!="#" & input$dnbgColour!="#"){
					g <- paste(g, 'theme(panel.background = element_rect(fill=input$dnbgFill, colour=input$dnbgColour,))', sep="+" )
				}
			# size - only affects ggsave, in ggsave call
			# theory line on/off
				if(input$dnTheoryCheck) {
					g <- paste(g, "geom_line(data = theory, aes(decile,estim_OR2))", sep="+")
				}
			# return the whole string in the end
			return(g)
		})

  # dnPlotString, a reactive expression to create a string containing a ggplot() command
  aucPlotString <- reactive({
    g = paste("ggplot(iris, aes(Sepal.Length, Sepal.Width) ) + geom_point()")
    return(g)
  })

		# decileNew plot
		output$decileNew <- renderPlot({
			# watch dnRefresh button
			input$dnRefresh

			# set up deciles data
			y <- read.csv("GPRS_deciles.csv", header = T, sep = ",")
			# change decile column to factor
			y$decile <- factor(y$decile)
      
      # plot with reduced functionality if no label file is provided
			if(is.null(input$file3) & !is.null(input$file1)) {
			  
			  shift = (as.numeric(factor(y$scoreFile))-1)/length(unique(y$scoreFile))
			  shift = (shift - median(shift))/2
        y$x = as.numeric(y$decile) + shift
       
			  g = ggplot(y) + geom_pointrange(aes(x, estimate, ymin=lowerCI, ymax=upperCI, colour=scoreFile), position=position_dodge(width=0.3)) + 
           ggtitle(input$decileNewTitle) + xlab(input$decileNewXlabel) + ylab(input$decileNewYlabel) + scale_x_continuous(breaks=1:10) #+ 
           #scale_colour_discrete(name = input$decileNewLegendTitle, labels = input$decileNewLegendLabels)
			  
			} else {
      
			# subset y based on intersection of 5 criteria
			d <- selectDataFrom5(y, input$dns1, input$dns2, input$dns3, input$dns4, input$dns5)

			# set up Hong's Theory data
			theory <- hongsTheoryCode(input$dnTheoryBins, input$dnTheoryK, input$dnTheoryH2)

			# create plot object from string containing command
			g <- eval(parse(text=dnPlotString()))
			# save plot, subject to size and resolution chosen by user
			para <- list(filename="www/plot.pdf", plot=g)
			if(input$dnSizeX!=0	){para$width 	<- input$dnSizeX*0.3937} # 0.3937 conversion factor for cm -> inches
			if(input$dnSizeY!=0	){para$height	<- input$dnSizeY*0.3937}
			
			#do.call(ggsave, para)
			# ggsave("www/plot.pdf", g, width=input$dnSizeX, height=input$dnSizeY, dpi=input$dnDPI)
			#g
			}
      return(g)
		})

  # decileNew plot
  output$aucPlot <- renderPlot({
    # watch dnRefresh button
    input$dnRefresh
    
    g <- eval(parse(text=aucPlotString()))
  
    return(g)
  })

  # plots roc curve (no ggplot)
  plot.roc = function(obs, pred, prec.recall=FALSE) {
    stopifnot(length(unique(union(c(0,1), obs))) == 2)
    obs.sort = obs[order(pred)]
    fn = cumsum(obs.sort)
    tn = cumsum(!obs.sort)
    tp = rev(cumsum(rev(obs.sort)))
    fp = rev(cumsum(rev(!obs.sort)))
    tpr = c(1,tp/(tp+fn),0) # sensitiyity, recall, tpr
    fpr = c(1,fp/(fp+tn),0)
    ppv = c(1,tp/(tp+fp),0) # precision, ppv, 1-fdr
    if(!prec.recall) {
      app = approx(fpr, tpr, n=1e3)
      auc = mean(app$y)
      plot(fpr, tpr, xlim=c(0,1), ylim=c(0,1), xlab='False positive rate', ylab='True positive rate', type='l', main=paste0('AUC: ', round(auc,3)))
      abline(0,1)
    }
    else {
      plot(tpr, ppv, xlim=c(0,1), ylim=c(0,1), xlab='True positive rate (recall, sensitivity)', ylab='Precision', type='l')
    }
  }


		# decileNew download plot
		output$dndPlot <- downloadHandler(
			filename = function() {
				"plot.pdf"
			},
			content = function(file) {
				file.copy("www/plot.pdf", file, overwrite=TRUE)
			}
		)

		# decileNew download selected data
		output$dndData <- downloadHandler(
			filename = function() {
				"selectedData.csv"
			},
			content = function(file) {
				y <- read.csv("GPRS_deciles.csv", header = T, sep = ",")
				d <- selectDataFrom5(y, input$dns1, input$dns2, input$dns3, input$dns4, input$dns5)
				write.table(d, file, sep = ",", row.names = FALSE)
			}
		)

  # download cumulative OR data
  output$cumulData <- downloadHandler(
    filename = function() {
      "cumulative_OR.csv"
    },
    content = function(file) {
      d <- getCumulDf()
      write.table(d, file, sep = ",", row.names = FALSE)
    }
  )

		# decileNewCodeUI
		output$decileNewCodeUI <- renderUI({
			x <- dnPlotString()
			y <- gsub(pattern="+", replacement=" +<br>&nbsp;&nbsp;", x=x, fixed=TRUE)
			z <- 'd <- read.csv("selectedData.csv", header=T, sep=",")<br>'
			a <- paste0(z,y)
			pre(HTML(a))
		})

	# tab cumulative
		# fileForDecile, a dynamic UI dropdown selection box to select data source file for decile graph
			output$fileForCumulativeUI <- renderUI({

				data <- read.csv("GPRS.csv", header = TRUE, sep = ",")
        

				selectInput("fileForCumulative", 
					"Add/remove files to draw data from",
					choices = levels(factor(data$scoreFile)),
					selected = if(input$cumulativeAddAll) levels(factor(data$scoreFile)) else levels(factor(data$scoreFile))[length(data$scoreFile)],
					width = "100%",
					selectize = TRUE,
					multiple = TRUE
				)
			})

    getCumulDf <- reactive({
      
      input$action
      data <- processUpload()
      if(is.null(input$fileForCumulative))
        index <- length(data$scoreFilenames)
      else
        index <- match(input$fileForCumulative, data$scoreFilenames)
      # index will be a vector of numerics
      d <- data.frame()
      # for every selected file...
      for (i in 1:length(index)) {
        
        df = data$tables[[index[i]]]
        
        df$NORM_SCORE = scale(df$SCORE)
        pheno.sorted = df$PHENO[order(df$SCORE)]-1 
        norm.score.sorted = sort(df$NORM_SCORE)
        total.odds = sum(pheno.sorted)/sum(!pheno.sorted)
        
        ##########################################################
        ################ continuous, all from top ################
        ##########################################################
        
        # odds for each subset 1..x, with x going from 1 to all individuals (individuals sorted by risk score)
        cumulative.odds = rev(
          cumsum(
            rev(pheno.sorted)
          )											/
            cumsum(
              rev(!pheno.sorted)
            ) 
        )
        # !pheno.sorted returns TRUE for 0 and FALSE for 1
        # rev(!pheno.sorted) is a reverse order version of the above (now ordered descendingly by score)
        # cumsum( rev(!pheno.sorted) ) replaces each element with the cumulative sum up until that point
        # 	i.e. the cumulative sum at each element of a descending version of pheno.sorted, where cases are 0 and controls are 1
        # the numerator is similar, but for cases = 1 and controls = 0
        #		therefore, 	the numerator 	is the cumulative number of cases, 		starting from top score going down
        #								the denominator is the cumulative number of controls, starting from top score going down
        # you take the odds of those, then reverse the order again (cumulative.odds is still ordered ascendingly by score)
        
        or = cumulative.odds / total.odds
            
        # norm.score.percentile can be used instead of norm.score.sorted in the dataframe "d" below; adjust the axis labels accordingly
        # norm.score.percentile = rank(norm.score.sorted)/length(norm.score.sorted)
        # rank(x) replaces each element with its rank, from 1 (i.e. the biggest value) to length(x) (the smallest value)
        
        e = data.frame(OR=or, risk.score=norm.score.sorted, scoreFile=data$scoreFilenames[index[i]])
        d <- rbind(d,e)
      }
      d
      
    })

		# cumulative plot
			output$cumulative <- renderPlot({

			  d <- getCumulDf()
        # set top 100 risk.score ORs to NA for each risk score file (they can be very unstable since they are based on few observations)
			  if(input$cumulativeShowTop) {
  			  top100 = c()
          for(f in unique(d$scoreFile)) {
            subs = d[d$scoreFile == f,]
            subs = subs[order(subs[,2], decreasing=TRUE),]
            top100 = c(top100, which(d$scoreFile == f & d$risk.score %in% subs$risk.score[1:min(100, nrow(subs))]))
          }
          d[top100, 1] = NA
			  }

        ggplot(d) + 
					geom_line( aes(risk.score, OR, colour=scoreFile) ) + 
					# theme(text=element_text(size=25)) + 
					geom_abline(slope=0, intercept=1) + 
					xlab('normalized risk score') + 
					ylab('odds ratio for all individuals scoring higher than x, relative to overall target sample odds')
			})

  output$gxeHeatmap = renderPlot({
    
    d = runGxE()
    d$scoreFilename = factor(d$scoreFilename, levels=rev(levels(d$scoreFilename)))
    siz = 16
    ggplot(d, aes(E, scoreFilename)) + 
      geom_tile(aes(fill = -log10(p)), colour = "white") + 
      scale_fill_gradient(low = "white", high = "steelblue") + 
      #scale_y_reverse() +
      xlab('') + 
      ylab('') + 
      theme(axis.text=element_text(size=siz), legend.text=element_text(size=siz), legend.title=element_text(size=siz))
    
  })

  
  # GxE scatterplot for selected score and E; only for binary phenotypes
  output$gxeScatter = renderPlot({
    
    data = processUpload()
    m = match(input$GxEselectSCORE, data$scoreFilenames)
    if(length(m) == 0) return(NULL)
    d = data[[1]][[m]]

    if(input$adjustForScatter) {
      allcov = input$covariatesChecklist
      thiscov = input$GxEselectE
      othercov = setdiff(allcov, thiscov)
      d[,thiscov] = lm(as.formula(paste0(thiscov, ' ~ ', paste0(othercov, collapse=' + ') )), d)$residuals
      d[,'SCORE'] = lm(as.formula(paste0('SCORE ~ ', paste0(othercov, collapse=' + ') )), d)$residuals
    }
    
    siz = 16
    ggplot(d, aes_string(input$GxEselectE, 'SCORE', color='as.factor(PHENO)')) + 
      geom_point() + 
      geom_smooth(method=lm, se=TRUE, fullrange=T) +
      theme(axis.text=element_text(size=siz), legend.text=element_text(size=siz), legend.title=element_blank(),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
            panel.background = element_blank(), axis.line = element_line(colour = "black"))
    
  })

#   output$gxe3d <- renderWebGL({
#     data = processUpload()
#     m = match(input$GxEselectSCORE, data$scoreFilenames)
#     if(length(m) == 0) return(NULL)
#     d = data[[1]][[m]]
#     plot3d(d$SCORE, d[,input$GxEselectE], d[,input$GxEselectE2], xlab='', ylab='', zlab='', col=ifelse(d==1, 'red', 'blue'))
#   })

	# tab AEGP
		# download for Vg2FromP.csv
		output$aegpdVg2 <- downloadHandler(
			filename = function() {
				"Vg2FromP.csv"
			},
			content = function(file) {
				file.copy("Vg2FromP.csv", file, overwrite=TRUE)
			}
		)

		# download for CorrFromP.csv
		output$aegpdCorr <- downloadHandler(
			filename = function() {
				"CorrFromP.csv"
			},
			content = function(file) {
				file.copy("CorrFromP.csv", file, overwrite=TRUE)
			}
		)

		output$aegpVg2 <- renderDataTable(
			{
				# re-run if button pressed
				input$action
				shiny::validate(
				  need(input$action != 0, "Please run the script!")
				)
				
				shiny::validate(
				  need(file.exists("Vg2FromP.csv"), "There are no input files which have a pThreshold of 1 and where discovery and target set are on the same scale!")
				)
				df = read.csv("Vg2FromP.csv", header=T, sep=",")
        df$pval = prettyNum(df$pval, digits=2, scientific=T); for(i in 4:6) {df[,i] = sprintf("%.2f", df[,i])}
        df
			},
			options = list(bPaginate = FALSE, bFilter = FALSE)
		)

		output$aegpCorr <- renderDataTable(
			{
				# re-run if button pressed
				input$action
				shiny::validate(
				  need(input$action != 0, "Please run the script!")
				)
        
				shiny::validate(
				  need(file.exists("CorrFromP.csv"), "There are no input files which have a pThreshold of 1 and where discovery and target set are on the same scale!")
				)
				df = read.csv("CorrFromP.csv", header=T, sep=",")
				df$pval = prettyNum(df$pval, digits=2, scientific=T); for(i in 4:6) {df[,i] = sprintf("%.2f", df[,i])}
        df
        
			},
			options = list(bPaginate = FALSE, bFilter = FALSE,
                     iDisplayLength=55, bSortClasses = TRUE, bAutoWidth=FALSE,
			     aoColumnDefs = list(list(sWidth=c("300px"), aTargets=c(list(3)))))
		)

  output$inputError = renderText({
    data = processUpload()
    data$check
  })

	# tab citations
		# entirely in ui.r for now
})