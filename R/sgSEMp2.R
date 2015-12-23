##' This function builds an gSEM model using gSEM principle 2. Principle 2 resembles the multiple regression principle in the way multiple predictors are considered simultaneously. Specifically, the first-level predictors to the system level variable, such as, Time and unit level variables, acted on the system level variable collectively by an additive model. This collective additive model can be found with a generalized stepwise variable selection (using the step() function in R, which performs variable selection on the basis of AIC) and this proceeds iteratively.
##' 
##' Data is analysed first using Principle 1 to find the best models. If needed, transformations based on the best models are applied to the predictors. Starting from the system response variable, each variable is regressed on all other variables except for the system response in an additive multiple regression model, which is reduced by a stepwise selection using stepAIC(). Then, for each selected variable, fitted regression for 6 selected functional forms and pick the best.
##'
##' @title Semi-supervised Generalized Structural Equation Modelling (gSEM) - Principle 2

##' @param x A dataframe, requiring at least 2 columns. By default its first column stores the main or primary influencing predictor, or exogenous variable e.g.., time, or a main predictor, the second column stores the response variable, and other columns store intermediate variables.
##' @param predictor A character string of the column name of the system predictor OR a numeric number indexing the column of the main predictor.
##' @param response A character string of the column name of the main response OR a numeric number indexing the column of the system response.
##' @return A list of the following items:
##'
##' \itemize{
##' \item "Graph": A network graph that contains the group and individual relationships between response and predictors determined by principle 2.
##' \item "res.print": A matrix. For each row, first column is the response variable, second column is the predictor, the other columns show corresponding summary information.
##'}

##' @seealso sgSEMp1() and plot.sgSEMp2()
##'
##' @export
##' @importFrom 'stats' 'coef'
##' @examples
##' data(acrylic)
##' ans <- sgSEMp2(acrylic)
##' ans$res.print
##' plot(ans)
##' 
##'
sgSEMp2 <- function(x, predictor = NULL, response = NULL){
  
	if(!missing(predictor)){
		if(!is.character(predictor)){
			if(is.wholenumber(predictor)){
				if(predictor < 1 | predictor > length(colnames(x)))
					stop("Predictor location out of range!")
				predictor.loc <- predictor
			}
		}else{ 
			if(!(predictor %in% colnames(x))){
				stop(paste0("Predictor '", predictor, "' does not exist!"))
			}
			predictor.loc <- which(colnames(x) == predictor)
		}
		neworder <- 1:length(colnames(x))
		neworder[1] <- predictor.loc
		neworder[predictor.loc] <- 1
		x <- x[neworder]
	}
	if(!missing(response)){
		response.loc <- which(colnames(x) == response)
		if(!is.character(response)){
			if(is.wholenumber(response)){
				if(response < 1 | response > length(colnames(x)))
					stop("Response location out of range!")
				response.loc <- response
			}
		}else{
			if(!(response %in% colnames(x))){
				stop(paste0("Response '", response, "' does not exist!"))
			}
			response.loc <- which(colnames(x) == response)
		}
		neworder <- 1:length(colnames(x))
		neworder[2] <- response.loc
		neworder[response.loc] <- 2
		x <- x[neworder]
	}

	###############################
	## The following function does multiple selection
	############################### 
	Multiple.relation <- function(iResp){

		Resp <- colnames(x)[iResp]   
		Var <- colnames(x)[c(-2, -iResp)]

		## get rid of NA's
		x1 <- x[c(Resp,Var)]
		x1 <- x1[apply(x1, 1, FUN = function(x) sum(is.na(x)) == 0),]
		
		################## START: APPLY TRANSFORMATION TO PREDICTORS ###################
		trans <- rep(NA,length=length(Var))
		trans1 <- rep(NA,length=length(Var))
		names(trans) <- Var
		for ( i in 1:length(Var) ) {
			# If "SL", x ---> x
			if ( p1.bestModels[Var[i],Resp] == "SL" ) {
				trans1[i] <- NA
				names(trans1)[i] <- Var[i]
				Quad_value <- coef(p1.allModels[Var[i],Resp,'SL'][[1]])
				if ( Quad_value[2]>0 ) {
					trans[i] <- paste(Resp,"=",format(Quad_value[1],scitific=T,digits=2),"+",format(Quad_value[2],scitific=T,digits=2),"*",Var[i],sep="")
				} else {
					trans[i] <- paste(Resp,"=",format(Quad_value[1],scitific=T,digits=2),format(Quad_value[2],scitific=T,digits=2),"*",Var[i],sep="")
				}
			}
			# If "SQuad", x ---> x^2
			if ( p1.bestModels[Var[i],Resp] == "SQuad" ) {
				x1[,Var[i]] <- (x1[,Var[i]])^2
				trans1[i] <- paste(Var[i],"_t=", Var[i],"^2",sep="")
				names(trans1)[i] <- Var[i]
				names(x1)[i+1] <- paste(names(x1)[i+1],"_t",sep="")
				Quad_value <- coef(p1.allModels[Var[i],Resp,'SQuad'][[1]])
				if ( Quad_value[2]>0 ) {
					trans[i] <- paste(Resp,"=",format(Quad_value[1],scitific=T,digits=2),"+",format(Quad_value[2],scitific=T,digits=2),"*",Var[i],"^2",sep="")
				} else {
					trans[i] <- paste(Resp,"=",format(Quad_value[1],scitific=T,digits=2),format(Quad_value[2],scitific=T,digits=2),"*",Var[i],"^2",sep="")
				}
			}
			# If "Exp", x ---> exp(x)
			if ( p1.bestModels[Var[i],Resp] == "Exp" ) {
				x1[,Var[i]] <- exp(x1[,Var[i]])
				trans1[i] <- paste(Var[i],"_t=e^", Var[i],sep="")
				names(trans1)[i] <- Var[i]
				names(x1)[i+1] <- paste(names(x1)[i+1],"_t",sep="")
				Quad_value <- coef(p1.allModels[Var[i],Resp,'Exp'][[1]])
				if ( Quad_value[2]>0 ) {
					trans[i] <- paste(Resp,"=",format(Quad_value[1],scitific=T,digits=2),"+",format(Quad_value[2],scitific=T,digits=2),"*e^",Var[i],sep="")
				} else {
					trans[i] <- paste(Resp,"=",format(Quad_value[1],scitific=T,digits=2),format(Quad_value[2],scitific=T,digits=2),"*e^",Var[i],sep="")
				}
			}
			# If "Log", x ---> log(x)
			if ( p1.bestModels[Var[i],Resp] == "Log" ) {
				x1[,Var[i]] <- log(x1[,Var[i]])
				trans1[i] <- paste(Var[i],"_t=log$", Var[i],"$",sep="")
				names(trans1)[i] <- Var[i]
				names(x1)[i+1] <- paste(names(x1)[i+1],"_t",sep="")
				Quad_value <- coef(p1.allModels[Var[i],Resp,'Log'][[1]])
				if ( Quad_value[2]>0 ) {
					trans[i] <- paste(Resp,"=",format(Quad_value[1],scitific=T,digits=2),"+",format(Quad_value[2],scitific=T,digits=2),"*log_",Var[i],sep="")
				} else {
					trans[i] <- paste(Resp,"=",format(Quad_value[1],scitific=T,digits=2),format(Quad_value[2],scitific=T,digits=2),"*log_",Var[i],sep="")
				}
			}
			# If "Quad", x ---> (x+b/2c)^2
			if ( p1.bestModels[Var[i],Resp] == "Quad" ) {
			    temp_column <- x1[,Var[i]]^2
				x1 <- cbind(x1,temp_column)
				names(x1)[ncol(x1)] <- paste(Var[i],"__2",sep="")
			
				Quad_value <- coef(p1.allModels[Var[i],Resp,'Quad'][[1]])
#				x1[,Var[i]] <- (x1[,Var[i]]+Quad_value[2]/(2*Quad_value[3]))^2
				names(trans1)[i] <- Var[i]
#				names(x1)[i+1] <- paste(names(x1)[i+1],"_t",sep="")
				if ( (Quad_value[2]>0)&(Quad_value[3]>0) ) {
					trans[i] <- paste(Resp,"=",format(Quad_value[1],scitific=T,digits=2),"+",format(Quad_value[2],scitific=T,digits=2),'*',Var[i],'+',format(Quad_value[3],scitific=T,digits=2),'*',Var[i],"^2",sep="")
					trans1[i] <- paste(Var[i],"_t=$", Var[i],"+",format(Quad_value[2]/(2*Quad_value[3]),scitific=T,digits=2),"$^2",sep="")
				}
				if ( (Quad_value[2]>0)&(Quad_value[3]<0) ) {
					trans[i] <- paste(Resp,"=",format(Quad_value[1],scitific=T,digits=2),"+",format(Quad_value[2],scitific=T,digits=2),'*',Var[i],format(Quad_value[3],scitific=T,digits=2),'*',Var[i],"^2",sep="")
					trans1[i] <- paste(Var[i],"_t=$", Var[i],"",format(Quad_value[2]/(2*Quad_value[3]),scitific=T,digits=2),"$^2",sep="")
				}
				if ( (Quad_value[2]<0)&(Quad_value[3]>0) ) {
					trans[i] <- paste(Resp,"=",format(Quad_value[1],scitific=T,digits=2),format(Quad_value[2],scitific=T,digits=2),'*',Var[i],'+',format(Quad_value[3],scitific=T,digits=2),'*',Var[i],"^2",sep="")
					trans1[i] <- paste(Var[i],"_t=$", Var[i],"",format(Quad_value[2]/(2*Quad_value[3]),scitific=T,digits=2),"$^2",sep="")
				}
				if ( (Quad_value[2]<0)&(Quad_value[3]<0) ) {
					trans[i] <- paste(Resp,"=",format(Quad_value[1],scitific=T,digits=2),format(Quad_value[2],scitific=T,digits=2),'*',Var[i],format(Quad_value[3],scitific=T,digits=2),'*',Var[i],"^2",sep="")
					trans1[i] <- paste(Var[i],"_t=$", Var[i],"+",format(Quad_value[2]/(2*Quad_value[3]),scitific=T,digits=2),"$^2",sep="")
				}
			}
			# If "nls", x ---> exp(cx)
			if ( p1.bestModels[Var[i],Resp] == "nls" ) {
				Quad_value <- coef(p1.allModels[Var[i],Resp,'nls'][[1]])
				x1[,Var[i]] <- exp(Quad_value[3]*x1[,Var[i]])
				trans1[i] <- paste(Var[i],"_t=exp$", format(Quad_value[3],scitific=T,digits=2),"*",Var[i],"$",sep="")
				names(trans1)[i] <- Var[i]
				names(x1)[i+1] <- paste(names(x1)[i+1],"_t",sep="")
				if ( Quad_value[2]>0 ) {
					trans[i] <- paste(Resp,"=",format(Quad_value[1],scitific=T,digits=2),"+",format(Quad_value[2],scitific=T,digits=2),"*e^",format(Quad_value[3],scitific=T,digits=2),Var[i],sep="")
					trans1[i] <- paste(Var[i],"_t=exp$", format(Quad_value[3],scitific=T,digits=2),"*",Var[i],"$",sep="")
				} else {
					trans[i] <- paste(Resp,"=",format(Quad_value[1],scitific=T,digits=2),format(Quad_value[2],scitific=T,digits=2),"*e^",format(Quad_value[3],scitific=T,digits=2),Var[i],sep="")
				}
			}
		}
		################## END: APPLY TRANSFORMATION TO PREDICTORS ###################
		Var_char <- paste(names(x1)[-1], collapse = "+")
		Rel <- paste0(Resp, "~", Var_char)
		lm.full.model <- do.call("lm", list(Rel, data=as.name("x1")))
		lm.best <- stepAIC(lm.full.model, direction = "backward", trace =FALSE)
		R2 <- summary(lm.best)$r.squared
		adj.R2 <- summary(lm.best)$adj.r.squared

		Var.relate <- names(coef(lm.best))[-1] 
		for ( m in 1:length(Var.relate) ) {
			Var.relate[m] <- gsub("_t","",Var.relate[m])
		}
		
		best.vars <- names(coef(lm.best))[-1] 
		for ( m in 1:length(best.vars) ) {
			best.vars[m] <- gsub("_t","",best.vars[m])
		}
				
		Var.relate1 <- subset(Var.relate, Var.relate != colnames(x)[1]) ## remove the predictor
#		tofit.flag[Var.relate1] <<- TRUE        ############## <<
#		fitted.flag[iResp-1] <<- TRUE       ############### <<
#        best.vars <- best.vars[-grep("__2",best.vars)]
		
		coeff_number <- format(lm.best$coeff,digits=2)
		names(coeff_number) <- NULL
		coeff_names <- names(lm.best$coeff)
		lm.best.text <- paste(Resp,"=",sep="")
		for ( ii in 1:length(coeff_names) ) {
			if (as.numeric(gsub(" ","",coeff_number[ii]))<0) {
			    if (coeff_names[ii]=='(Intercept)') {
					lm.best.text <- paste(lm.best.text,coeff_number[ii],sep="")
				} else {
					if (length(grep("_t",coeff_names[ii]))+length(grep("__2",coeff_names[ii]))==0) {
						lm.best.text <- paste(lm.best.text,coeff_number[ii],'*',coeff_names[ii],sep="")
					}
					if (length(grep("_t",coeff_names[ii]))!=0) {
						names_temp <- substr(coeff_names[ii],1,regexpr("_t",coeff_names[ii])[1]-1)
						names_temp2 <- substr(trans1[names_temp],regexpr("=",trans1[names_temp])[1]+1,nchar(trans1[names_temp]))
						lm.best.text <- paste(lm.best.text,coeff_number[ii],'*',names_temp2,sep="")
					}
					if (length(grep("__2",coeff_names[ii]))!=0) {
						names_temp <- gsub("__","^",coeff_names[ii])
						lm.best.text <- paste(lm.best.text,coeff_number[ii],'*',names_temp,sep="")
					}
				}
			} else {
				if (coeff_names[ii]=='(Intercept)') {
					lm.best.text <- paste(lm.best.text,gsub(" ","",coeff_number[ii]),sep="")
				} else {
					if (length(grep("_t",coeff_names[ii]))+length(grep("__2",coeff_names[ii]))==0) {
						lm.best.text <- paste(lm.best.text,'+',gsub(" ","",coeff_number[ii]),'*',coeff_names[ii],sep="")
					}
					if (length(grep("_t",coeff_names[ii]))!=0) {
						names_temp <- substr(coeff_names[ii],1,regexpr("_t",coeff_names[ii])[1]-1)
						names_temp2 <- substr(trans1[names_temp],regexpr("=",trans1[names_temp])[1]+1,nchar(trans1[names_temp]))
						lm.best.text <- paste(lm.best.text,'+',gsub(" ","",coeff_number[ii]),'*',names_temp2,sep="")
					}
					if (length(grep("__2",coeff_names[ii]))!=0) {
						names_temp <- gsub("__","^",coeff_names[ii])
						lm.best.text <- paste(lm.best.text,'+',gsub(" ","",coeff_number[ii]),'*',names_temp,sep="")
					}
				}
			}
		}
		
		# Make a new row in the Res.print table for each of the predictor variables in the additive model
		
		for (n in 1:length(best.vars)){
			if (length(grep("__2",best.vars[n]))==0) {
				Res.print.newrow <- matrix(NA, nrow = 1, ncol = nRes)
				Res.print.newrow <- as.data.frame(Res.print.newrow)
				colnames(Res.print.newrow) <- c("resp",          "var",          "adj-r2",      "Transformation",  "Group_Model", 
	                         "Group_p-value", "Group_R2_aR2", "Indv_Model",  "Indv_p-value",     "Indv_R2_aR2",
							 "Rank")
				# This fills the table with blanks for now while the real values are obtained
				# Need to break this into individual cell assignments to prevent factorization from combining datatypes
				Res.print.newrow[1,1] <- Resp
				Res.print.newrow[1,2]<- best.vars[n]
#			Res.print.newrow[1,3] <- as.character(lm.best$call[2])
				Res.print.newrow[1,5] <- lm.best.text
				Res.print.newrow[1,4] <- trans1[best.vars[n]]
				Res.print.newrow[1,6] <- paste(round(summary(lm.best)$coeff[,4],digits=3),collapse="  ")
				Res.print.newrow[1,7] <- paste(round(R2,digits=3),round(adj.R2,digits=3),sep="  ")
				Res.print.newrow[1,3] <- format(adj.R2,scitific=T,digits=2)
				Res.print.newrow[1,11] <- which(order(summary(lm.best)$coeff[-1,4])==n)
				Res.print.newrow[1,8] <- trans[best.vars[n]]
				Individual_summary <- summary(p1.allModels[best.vars[n],Resp,p1.bestModels[best.vars[n],Resp]][[1]])
				Res.print.newrow[1,9] <- paste(round(Individual_summary$coeff[,4],digits=3),collapse="  ")
				if (p1.bestModels[best.vars[n],Resp]!="nls") {
					Res.print.newrow[1,10] <- paste(round(Individual_summary$r.squared,digits=3),round(Individual_summary$adj.r.squared,digits=3),sep="  ")
				}
				Res.print <<- rbind(Res.print, Res.print.newrow)
			} else {
			    new_string <- substr(best.vars[n],1,regexpr("__",best.vars[n])[1]-1)
				Res.print.newrow <- matrix(NA, nrow = 1, ncol = nRes)
				Res.print.newrow <- as.data.frame(Res.print.newrow)
				colnames(Res.print.newrow) <- c("resp",          "var",          "adj-r2",      "Transformation",  "Group_Model", 
	                         "Group_p-value", "Group_R2_aR2", "Indv_Model",  "Indv_p-value",     "Indv_R2_aR2",
							 "Rank")
				# This fills the table with blanks for now while the real values are obtained
				# Need to break this into individual cell assignments to prevent factorization from combining datatypes
				Res.print.newrow[1,1] <- Resp
				Res.print.newrow[1,2]<- new_string
#			Res.print.newrow[1,3] <- as.character(lm.best$call[2])
				Res.print.newrow[1,5] <- lm.best.text
				Res.print.newrow[1,4] <- trans1[new_string]
				Res.print.newrow[1,6] <- paste(round(summary(lm.best)$coeff[,4],digits=3),collapse="  ")
				Res.print.newrow[1,7] <- paste(round(R2,digits=3),round(adj.R2,digits=3),sep="  ")
				Res.print.newrow[1,3] <- format(adj.R2,scitific=T,digits=2)
				Res.print.newrow[1,11] <- which(order(summary(lm.best)$coeff[-1,4])==n)
				Res.print.newrow[1,8] <- trans[new_string]
				Individual_summary <- summary(p1.allModels[new_string,Resp,p1.bestModels[new_string,Resp]][[1]])
				Res.print.newrow[1,9] <- paste(round(Individual_summary$coeff[,4],digits=3),collapse="  ")
				if (p1.bestModels[new_string,Resp]!="nls") {
					Res.print.newrow[1,10] <- paste(round(Individual_summary$r.squared,digits=3),round(Individual_summary$adj.r.squared,digits=3),sep="  ")
				}
				insert_ind <- 0
				Res.temp <- Res.print[-1,]
				for ( mm in 1:nrow(Res.temp) ) {
					if ( Res.temp[mm,1]==Resp ) {
						if ( Res.temp[mm,2]==new_string ) {
							insert_ind <- 1
						}
					}
				}
				if (insert_ind==0) {
					Res.print <<- rbind(Res.print, Res.print.newrow)
				}
			}	
		}
	}
	
	###############################
	## Main scripts; Above two functions are called
	############################### 
  
	# Apply principle 1 on the data to find the best models
	if ( !is.null(predictor) ) {
		predictor.p1 <- predictor
	} else {
		predictor.p1 <- 1
	}
	if ( !is.null(response) ) {
		response.p1 <- response
	} else {
		response.p1 <- 2
	}
	p1.result <- sgSEMp1(x, predictor = predictor.p1, response = response.p1)
	p1.bestModels <- rbind(p1.result$bestModels[1,],rep(NA,length(colnames(p1.result$bestModels))),p1.result$bestModels[-1,])
	p1.bestModels <- cbind(rep(NA,length(rownames(p1.bestModels))),p1.bestModels)
	rownames(p1.bestModels) <- c(rownames(p1.result$bestModels)[1],colnames(p1.result$bestModels)[1],rownames(p1.result$bestModels)[-1])
	colnames(p1.bestModels) <- rownames(p1.bestModels)
	
	p1.allModels <- vector( "list", length(colnames(p1.bestModels)) * length(colnames(p1.bestModels)) * length(names(p1.result$allModels[1,1,])) )
	dim(p1.allModels) <- c(length(colnames(p1.bestModels)),length(colnames(p1.bestModels)), length(names(p1.result$allModels[1,1,])))
    p1.allModels[,,] <- NA
    dimnames(p1.allModels) <- list(colnames(p1.bestModels), colnames(p1.bestModels), names(p1.result$allModels[1,1,]))
	
	temp <- rbind(p1.result$allModels[1,,"SL"],rep(NA,length(colnames(p1.result$allModels[,,"SL"]))),p1.result$allModels[-1,,"SL"])
	temp <- cbind(rep(NA,length(rownames(temp))),temp)
	p1.allModels[,,"SL"] <- temp

	temp <- rbind(p1.result$allModels[1,,"Quad"],rep(NA,length(colnames(p1.result$allModels[,,"Quad"]))),p1.result$allModels[-1,,"Quad"])
	temp <- cbind(rep(NA,length(rownames(temp))),temp)
	p1.allModels[,,"Quad"] <- temp
	
	temp <- rbind(p1.result$allModels[1,,"SQuad"],rep(NA,length(colnames(p1.result$allModels[,,"SQuad"]))),p1.result$allModels[-1,,"SQuad"])
	temp <- cbind(rep(NA,length(rownames(temp))),temp)
	p1.allModels[,,"SQuad"] <- temp
	
	temp <- rbind(p1.result$allModels[1,,"Exp"],rep(NA,length(colnames(p1.result$allModels[,,"Exp"]))),p1.result$allModels[-1,,"Exp"])
	temp <- cbind(rep(NA,length(rownames(temp))),temp)
	p1.allModels[,,"Exp"] <- temp

	temp <- rbind(p1.result$allModels[1,,"Log"],rep(NA,length(colnames(p1.result$allModels[,,"Log"]))),p1.result$allModels[-1,,"Log"])
	temp <- cbind(rep(NA,length(rownames(temp))),temp)
	p1.allModels[,,"Log"] <- temp

	temp <- rbind(p1.result$allModels[1,,"nls"],rep(NA,length(colnames(p1.result$allModels[,,"nls"]))),p1.result$allModels[-1,,"nls"])
	temp <- cbind(rep(NA,length(rownames(temp))),temp)
	p1.allModels[,,"nls"] <- temp

#	temp <- rbind(p1.result$allModels[1,,"CPSLSL"],rep(NA,length(colnames(p1.result$allModels[,,"CPSLSL"]))),p1.result$allModels[-1,,"CPSLSL"])
#	temp <- cbind(rep(NA,length(rownames(temp))),temp)
#	p1.allModels[,,"CPSLSL"] <- temp

	rownames(p1.allModels) <- colnames(p1.bestModels)
	# End apply principle 1
  
	nVar <- ncol(x) #total number of variables
	nRVar <- nVar - 1 #total number of possible response variables
  
#	fitted.flag <- rep(FALSE, nRVar)       #indicator of response variables have been tested, dim2 response
#	names(fitted.flag) <- colnames(x)[-1]
#	tofit.flag <- c(TRUE, rep(FALSE, nRVar-1)) #indicator of all reponse varibals need to be tested, dim2 response
#	names(tofit.flag) <- colnames(x)[-1]

	nRes <- 11 # Number of cells in print variable; see below
	Res.print <- matrix(NA, nrow = 1, ncol = nRes) 
#	colnames(Res.print) <- c("resp",	      "var",	  "Group_Model", 	"Transformation",	"Group_p-value",
#							 "Group_R2_aR2",  "Rank",     "Indv_Model",     "Indv_p-value",     "Indv_R2_aR2",
#                            "adj-r2")
	colnames(Res.print) <- c("resp",          "var",          "adj-r2",      "Transformation",  "Group_Model", 
	                         "Group_p-value", "Group_R2_aR2", "Indv_Model",  "Indv_p-value",     "Indv_R2_aR2",
							 "Rank")

#	while(sum(tofit.flag - fitted.flag) != 0){   
#		iResp <- which(tofit.flag != fitted.flag)[1] + 1
	for ( iResp in 2:ncol(x) ) {
		Multiple.relation(iResp)
	}
	Res.print <- Res.print[-1,]
	res <- list(res.print = Res.print) #outputs are three tables
	class(res) <- c("sgSEMp2","list")
	invisible(res)
}
