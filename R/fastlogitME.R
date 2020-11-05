fastlogitME<-function(model, at = NULL, vars = NULL, conf.band = .95){
  if(!(conf.band>0&conf.band<1)){
    stop("The bandwidth of the confidence interval should be between 0 and 1")
  }
  dataframe<-model.frame(model)
  
  regvars=names(dataframe)[2:length(names(dataframe))]
  base<-data.frame("Intercept"=1)
  for(var in regvars){
    if(is.factor(dataframe[,var])|length(unique(dataframe[,var]))==2&all(c(0,1)%in%dataframe[,var])){#Test for dummy variable
      if(is.factor(dataframe[,var])){#A dummy variable with levels (factor)
        add<-as.data.frame(t(rep(0, nlevels(dataframe[,var]))))
        if(var%in%names(at)){#If a specific value is specified in at then use this value. Else use the reference category
          add[which(levels(dataframe[,var])==at[var])]=1
        }else{
          add[1]=1
        }
        names(add)=paste0(var, levels(dataframe[,var]))
        add<-add[2:length(add)]#Reference category does not have a coefficient
        base<-cbind(base, add)
      }else{#A dummy variable with zeros and ones
        if(var%in%names(at)){#If a specific value is specified in at then use this value. Else use the reference category
          base[,var]=at[var]  
        }else{
          base[,var]=0
        }
      }
    }else{
      if(var%in%names(at)){#If a specific value is specified in at then use this value. Else use mean
        base[,var]=at[var]  
      }else{
        base[,var]=mean(dataframe[,var])
      }
    }
  }
  
  #Test if interaction variables exist and if so create variables for them
  if(any(grepl(":", names(model$coefficients)))){
    for(var in names(model$coefficients)[grepl(":", names(model$coefficients))]){
      if(length(gregexpr(":", var)[[1]])==1){
        var1=substr(var, 1, gregexpr(":", var)[[1]][1]-1)
        var2=substr(var, gregexpr(":", var)[[1]][1]+1, nchar(var))
        base[,var]=base[,var1]*base[,var2]
      }
      if(length(gregexpr(":", var)[[1]])==2){
        var1=substr(var, 1, gregexpr(":", var)[[1]][1]-1)
        var2=substr(var, gregexpr(":", var)[[1]][1]+1, gregexpr(":", var)[[1]][2]-1)
        var3=substr(var, gregexpr(":", var)[[1]][2]+1, nchar(var))
        base[,var]=base[,var1]*base[,var2]*base[,var3]
      }
      if(length(gregexpr(":", var)[[1]])>2){
        print("This function is designed to incorporate up to 3 interaction terms and not more.")
        break
      }
    }
  }
  
  if(is.null(vars)){#If the user does not specify which variables then all variables are used
    vars=regvars
  }
  for(var in vars){#Replace factor variables by their respective levels
    if(is.factor(dataframe[,var])){
      if(which(vars==var)==1){#Factor variable that needs to be replaced by its levels is the first in the vars list
        vars<-c(paste0(var,levels(dataframe[,var])[2:nlevels(dataframe[,var])]), vars[2:length(vars)])
        next
      }
      if(which(vars==var)>1&which(vars==var)<length(vars)){#Factor variable that needs to be replaced by its levels is between the first and the last in the vars list
        vars<-c(vars[1:(which(vars==var)-1)],paste0(var,levels(dataframe[,var])[2:nlevels(dataframe[,var])]), vars[(which(vars==var)+1):length(vars)])
        next
      }
      if(which(vars==var)>1&which(vars==var)==length(vars)){#Factor variable that needs to be replaced by its levels is the last in the vars list
        vars<-c(vars[1:(which(vars==var)-1)],paste0(var,levels(dataframe[,var])[2:nlevels(dataframe[,var])]))
        next
      }
    }
  }
  fit=as.matrix(base)%*%model$coef
  Results<-data.frame("Variable" = vars, "ME" = 0, "Confupper" = 0, "Conflower" = 0, "p" = 0)
  for(var in vars){
    dybase<-base
    if(!(var%in%names(dataframe))&var%in%names(model$coefficients)){#Test for dummy variable of the factor type
      dybase[,var]=1
      #Test if var is involved in interaction variables and if so update these values
      if(any(grepl(":", names(model$coefficients))&grepl(var, names(model$coefficients)))){
        for(intervar in names(model$coefficients)[grepl(":", names(model$coefficients))&grepl(var, names(model$coefficients))]){
          if(length(gregexpr(":", intervar)[[1]])==1){
            intervar1=substr(intervar, 1, gregexpr(":", intervar)[[1]][1]-1)
            intervar2=substr(intervar, gregexpr(":", intervar)[[1]][1]+1, nchar(intervar))
            dybase[,intervar]=dybase[,intervar1]*dybase[,intervar2]
          }
          if(length(gregexpr(":", intervar)[[1]])==2){
            intervar1=substr(intervar, 1, gregexpr(":", intervar)[[1]][1]-1)
            intervar2=substr(intervar, gregexpr(":", intervar)[[1]][1]+1, gregexpr(":", intervar)[[1]][2]-1)
            intervar3=substr(intervar, gregexpr(":", intervar)[[1]][2]+1, nchar(intervar))
            dybase[,intervar]=dybase[,intervar1]*dybase[,intervar2]*dybase[,intervar3]
          }
          if(length(gregexpr(":", var)[[1]])>2){
            print("This function is designed to incorporate up to 3 interaction terms and not more.")
            break
          }
        }
      }
      fit2=as.matrix(dybase)%*%model$coef
      Results$ME[Results$Variable==var]=exp(fit2)/(1+exp(fit2))-exp(fit)/(1+exp(fit))
      
      SEbase<-dybase-base
      SEbase[SEbase<0]=0
      se=diag(sqrt(as.matrix(SEbase) %*% vcov(model) %*% t(as.matrix(SEbase))))
      Results$Confupper[Results$Variable==var]=exp(fit2+se*(qnorm((1-(1-conf.band)/2))))/(1+exp(fit2+se*(qnorm((1-(1-conf.band)/2)))))-exp(fit)/(1+exp(fit))
      Results$Conflower[Results$Variable==var]=exp(fit2-se*(qnorm((1-(1-conf.band)/2))))/(1+exp(fit2-se*(qnorm((1-(1-conf.band)/2)))))-exp(fit)/(1+exp(fit))
    }else{
      if(length(unique(dataframe[,var]))==2&all(c(0,1)%in%dataframe[,var])){#Non factor but still dummy variable type
        dybase[,var]=1
        #Test if var is involved in interaction variables and if so update these values
        if(any(grepl(":", names(model$coefficients))&grepl(var, names(model$coefficients)))){
          for(intervar in names(model$coefficients)[grepl(":", names(model$coefficients))&grepl(var, names(model$coefficients))]){
            if(length(gregexpr(":", intervar)[[1]])==1){
              intervar1=substr(intervar, 1, gregexpr(":", intervar)[[1]][1]-1)
              intervar2=substr(intervar, gregexpr(":", intervar)[[1]][1]+1, nchar(intervar))
              dybase[,intervar]=dybase[,intervar1]*dybase[,intervar2]
            }
            if(length(gregexpr(":", intervar)[[1]])==2){
              intervar1=substr(intervar, 1, gregexpr(":", intervar)[[1]][1]-1)
              intervar2=substr(intervar, gregexpr(":", intervar)[[1]][1]+1, gregexpr(":", intervar)[[1]][2]-1)
              intervar3=substr(intervar, gregexpr(":", intervar)[[1]][2]+1, nchar(intervar))
              dybase[,intervar]=dybase[,intervar1]*dybase[,intervar2]*dybase[,intervar3]
            }
            if(length(gregexpr(":", var)[[1]])>2){
              print("This function is designed to incorporate up to 3 interaction terms and not more.")
              break
            }
          }
        }
        fit2=as.matrix(dybase)%*%model$coef
        Results$ME[Results$Variable==var]=exp(fit2)/(1+exp(fit2))-exp(fit)/(1+exp(fit))
        
        SEbase<-dybase-base
        SEbase[SEbase<0]=0
        se=diag(sqrt(as.matrix(SEbase) %*% vcov(model) %*% t(as.matrix(SEbase))))
        Results$Confupper[Results$Variable==var]=exp(fit2+se*(qnorm((1-(1-conf.band)/2))))/(1+exp(fit2+se*(qnorm((1-(1-conf.band)/2)))))-exp(fit)/(1+exp(fit))
        Results$Conflower[Results$Variable==var]=exp(fit2-se*(qnorm((1-(1-conf.band)/2))))/(1+exp(fit2-se*(qnorm((1-(1-conf.band)/2)))))-exp(fit)/(1+exp(fit))
        Results$p[Results$Variable==var]=pnorm(abs((fit2-fit)/se), lower.tail = FALSE)
      }else{#Continuous variable so continuous marginal effect
        dybase[,var]=dybase[,var]+1e-7
        #Test if var is involved in interaction variables and if so update these values
        if(any(grepl(":", names(model$coefficients))&grepl(var, names(model$coefficients)))){
          for(intervar in names(model$coefficients)[grepl(":", names(model$coefficients))&grepl(var, names(model$coefficients))]){
            if(length(gregexpr(":", intervar)[[1]])==1){
              intervar1=substr(intervar, 1, gregexpr(":", intervar)[[1]][1]-1)
              intervar2=substr(intervar, gregexpr(":", intervar)[[1]][1]+1, nchar(intervar))
              dybase[,intervar]=dybase[,intervar1]*dybase[,intervar2]
            }
            if(length(gregexpr(":", intervar)[[1]])==2){
              intervar1=substr(intervar, 1, gregexpr(":", intervar)[[1]][1]-1)
              intervar2=substr(intervar, gregexpr(":", intervar)[[1]][1]+1, gregexpr(":", intervar)[[1]][2]-1)
              intervar3=substr(intervar, gregexpr(":", intervar)[[1]][2]+1, nchar(intervar))
              dybase[,intervar]=dybase[,intervar1]*dybase[,intervar2]*dybase[,intervar3]
            }
            if(length(gregexpr(":", var)[[1]])>2){
              print("This function is designed to incorporate up to 3 interaction terms and not more.")
              break
            }
          }
        }
        fit2=as.matrix(dybase)%*%model$coef
        Results$ME[Results$Variable==var]=(exp(fit2)/(1+exp(fit2))-exp(fit)/(1+exp(fit)))/1e-7
        
        SEbase<-dybase-base
        SEbase[SEbase<0]=0
        se=diag(sqrt(as.matrix(SEbase) %*% vcov(model) %*% t(as.matrix(SEbase))))
        Results$Confupper[Results$Variable==var]=(exp(fit2+se*(qnorm((1-(1-conf.band)/2))))/(1+exp(fit2+se*(qnorm((1-(1-conf.band)/2)))))-exp(fit)/(1+exp(fit)))/1e-7
        Results$Conflower[Results$Variable==var]=(exp(fit2-se*(qnorm((1-(1-conf.band)/2))))/(1+exp(fit2-se*(qnorm((1-(1-conf.band)/2)))))-exp(fit)/(1+exp(fit)))/1e-7
        Results$p[Results$Variable==var]=pnorm(abs((fit2-fit)/se), lower.tail = FALSE)
      }
    }
  }
  Results<-Results[!(Results$ME==0&Results$Confupper==0&is.nan(Results$p)),]
  return(Results)
}
