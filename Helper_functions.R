
# Rank normalize function 
rankNorm <- function(x) {
  qnorm((rank(x, ties.method = "random", na.last = "keep") - 0.5) / sum(!is.na(x)))
}

########################     Check normality and missing value
check_NM <- function(data, id) {
  data_check<-data %>% remove_rownames() %>% column_to_rownames(id)
  cat("\n==============================================\n")
  cat("\n --- Checking for normality and missing values --- \n")
  res <- list()
  for (i in 1:ncol(data_check)){
    ks<-suppressWarnings(ks.test(data_check[,i], "pnorm", mean(data_check[,i]), sd(data_check[,i])))
    res[[names(data_check)[i]]] <- list(pval<-ks[[2]],
                                      mean<-round(mean(data_check[,i], na.rm = TRUE), 3),
                                      sd<-round(sd(data_check[,i], na.rm = TRUE), 3),
                                      missing<- sum(is.na(data_check[,i])) )
  }
  check_results<-as.data.frame(do.call(rbind, res))
  check_results$features<-row.names(check_results)
  check_results[,1:4]<-as.data.frame(lapply(check_results[,1:4], as.numeric))
  colnames(check_results)<-c('pval','mean','sd','n_missing','features')
  
  cat("Missing:", sum(check_results$n_missing>0), "out of", nrow(check_results),  "\n")
  
  check_nomissing<-subset(check_results,n_missing==0)
  cat("Failed normality %:", round((sum(check_nomissing$pval<0.05)/nrow(check_results))*100,2), "\n")
  
  cat("Range mean: [", round(min(check_results$mean, na.rm = TRUE), 3), ",", 
      round(max(check_results$mean, na.rm = TRUE), 3), "]\n")
  cat("Range sd: [", round(min(check_results$sd, na.rm = TRUE), 3), ",", 
      round(max(check_results$sd, na.rm = TRUE), 3), "]\n")
  cat("\n==============================================\n")
  
  return(check_results) 
}

########################    Rank-normalize and/or Half-minimum imputation
RN_IMP<-function(data, id, miss_prop, check_results, RN=TRUE, IMP=TRUE) {
   
   cat("\n==============================================\n")
   cat("Removing features with missing value proportion over", miss_prop, "\n")
   missing_prop <- data.frame(Compound_ID = check_results$features, 
                              missing_prop = round(check_results$n_missing/nrow(data),3) )
   
   inds_acceptable <- which(missing_prop$missing_prop < miss_prop)
   Compound_ID_keep <- missing_prop$Compound_ID[inds_acceptable]
   data_proc <- data %>% remove_rownames() %>% column_to_rownames(id)
   data_proc <- data_proc[, Compound_ID_keep]
   
   # Half-minimum imputation
   if (IMP==TRUE) {
     cat("Performing half-minimum imputation...\n")
     feature_mins <- apply(data_proc, 2, function(x){min(x, na.rm = TRUE)})
   
   for (i in 1:ncol(data_proc)){
     mis_inds_cur_chem <- which(is.na(data_proc[,colnames(data_proc)[i]]))
     if (length(mis_inds_cur_chem) > 0){
       data_proc[mis_inds_cur_chem,colnames(data_proc)[i]] <- feature_mins[colnames(data_proc)[i]]/2
     }
     } 
   } 
   
   # Rank-normalization
   if (RN==TRUE) {
     cat("Performing rank normalization...\n")
     cat("\n==============================================\n")
     for (i in 1:ncol(data_proc)) {
     data_proc[ ,colnames(data_proc)[i]] <- rankNorm(data_proc[ ,colnames(data_proc)[i]])
     }
   }
   return(data_proc)
}

########################    Filtering features if necessary 
# Filter features for association with outcome of interest 
# Covariates are provided as string, such as categorical_covariates<- c('sex','race')
# "id" should be the name of the id column in phefile
filter<-function(data_proc, id=NULL, phefile, outcome, categorical_covariates, numeric_covariates) {
  if (is.null(id)) {
    stop("Please provide the id column in phenotype file (phefile)")
  } 
  
  # Check for categorical and numeric covariates
  phefile[categorical_covariates]<-lapply(phefile[categorical_covariates], as.character)
  phefile[numeric_covariates]<-lapply(phefile[numeric_covariates], function(x) as.numeric(as.character(x)))
  
  # Filtering for association
  data_temp<- data_proc %>% rownames_to_column(var = id)
  dall<-merge(phefile[,c(id,outcome,categorical_covariates,numeric_covariates)],data_temp,by=id)
  nvar<-length(outcome)+length(categorical_covariates)+length(numeric_covariates)+1
  rall<-list()
  for (i in 1:length(outcome)) {
    if (all(dall[,outcome[i]] %in% c(0, 1))) {
      res<-data.frame()
      for (j in (nvar+1):ncol(dall)){
        form <- as.formula(paste(outcome[i],'~', colnames(dall)[j], '+',
                               paste(c(categorical_covariates,numeric_covariates), collapse='+'))) 
        fit<-glm( form,data=dall, family = binomial(link = "logit")) 
        coef<-as.data.frame(summary(fit)$coefficients)[2,]
        coef$Feature<-colnames(dall)[j]
        res<-rbind(res,coef)
      }
      res$outcome<-outcome[i]
      names(res)[4]<-'pval'
      rall[[as.character(outcome[i])]]<-res
      cat("----------------------------------------------\n\n")
      cat(sum(res$pval<0.05),"of",nrow(res),"features are associated with",outcome[i], "\n")
      cat("----------------------------------------------\n\n")
  } else {
    res<-data.frame()
    for (j in (nvar+1):ncol(dall)){
      form <- as.formula(paste(outcome[i],'~', colnames(dall)[j], '+',
                               paste(c(categorical_covariates,numeric_covariates), collapse='+')))
      fit<-lm( form,data=dall) #
      coef<-as.data.frame(summary(fit)$coefficients)[2,]
      coef$Feature<-colnames(dall)[j]
      res<-rbind(res,coef)
    }
    res$outcome<-outcome[i]
    names(res)[4]<-'pval'
    rall[[as.character(outcome[i])]]<-res
    cat("----------------------------------------------\n\n")
    cat(sum(res$pval<0.05),"of",nrow(res),"features are associated with",outcome[i], "\n")
    cat("----------------------------------------------\n\n")
  }
  }
  return(rall)
}

########################    Elastic net
ENselect<-function(data_proc,id=NULL, phefile, filtered_res=NULL, outcome, 
                   alpha, categorical_covariates, numeric_covariates) {
  if (is.null(id)) {
    stop("Please provide the id column in phenotype file (phefile)")
  } 
  
  # Check for categorical and numeric covariates
    phefile[categorical_covariates]<-lapply(phefile[categorical_covariates], as.character)

    phefile[numeric_covariates]<-lapply(phefile[numeric_covariates], function(x) as.numeric(as.character(x)))
  
  # Build matrix for EN
  row.names(phefile)<-phefile[,id]
  form<-as.formula(paste('~',paste(c(categorical_covariates,numeric_covariates,outcome), collapse='+')))
  phe <- as.data.frame(model.matrix(form, data = phefile)[, -1])  
  phe <- phe %>% rownames_to_column(var = id)
  
  alpha_grid <- alpha   # 1.0 = LASSO
  best_param<-data.frame(Outcome=outcome,alpha=0,lambda=0,performance=0)
  EN_feature<-data.frame(Feature=colnames(data_proc))
  
  for (i in 1:length(outcome)) {
    cat("----------------------------------------------\n\n")
    cat("Running Elastic Net on outcome", outcome[i], "\n")
    
    # If filtered:
    if (!is.null(filtered_res)){
      filtered<-subset(filtered_res[[outcome[i]]], pval<0.05)
      data_temp<- data_proc[,c(which(colnames(data_proc)%in%filtered$Feature))]
    } else {
      data_temp <- data_proc
    }

    data_temp<- data_temp %>% rownames_to_column(var = id)
    
    x<-merge(phe[, c(colnames(phe)[1:(ncol(phe)-length(outcome))],outcome[i] )],data_temp,by=id) 
    x<-x%>% remove_rownames() %>% column_to_rownames(id)
    x<-x[,c(which(colnames(x)==outcome[i]),
            1:(which(colnames(x)==outcome[i])-1),
            (which(colnames(x)==outcome[i])+1):ncol(x))]
    
    n_chem <- ncol(data_temp)-1
    p.fac = c(rep(0, ncol(x) - n_chem-1), rep(1, n_chem)) # intercept and covariates not included 
    x<-as.matrix(x[complete.cases(x),])
    y <- as.numeric(x[,1])

    perf <- data.frame()
    is_binomial <- all(y %in% c(0, 1))
    fam <- if(is_binomial) "binomial" else "gaussian"
    type <- if(is_binomial) "auc" else "mse"
      
      for (a in alpha_grid) {
        message("Running alpha = ", a)
        
        set.seed(1997) # Keep folds consistent across alphas
        cvfit <- cv.glmnet(
          x = x[,-1], y = y, family = fam, 
          alpha = a, nfolds = 5, type.measure = type,
          standardize = TRUE, penalty.factor = p.fac)
        
        lambda_best <- cvfit$lambda.min     # Best lambda
        
        # Calculate performance (AUC or CV-R2)
        perf_val <- if(is_binomial) max(cvfit$cvm) else (1 - (min(cvfit$cvm) / var(y)))
        
        # Save results
        perf <- rbind(perf, data.frame(alpha = a, lambda = lambda_best,
                                       performance = perf_val))
      }
      
      # Compare and pick best α
      best <- which.max(perf$performance)
      best_alpha <- perf$alpha[best]
      selected_lambda <- perf$lambda[best]
      best_param$alpha[i]<-best_alpha # 
      best_param$lambda[i]<-selected_lambda # 
      best_param$performance[i]<-perf$performance[best]# 
      
      # Run final model
      fit_full <- glmnet(x = x[,-1], y = y,  family = fam,
                         alpha = best_alpha, lambda = selected_lambda,
                         penalty.factor = p.fac, standardize = TRUE)
      all_coefs <- as.matrix(coef(fit_full))[,1]
      all_coefs <- all_coefs[names(all_coefs) != "(Intercept)"]
      # select only features
      chem_coefs <- tail(all_coefs, n_chem) 
      selected <- as.data.frame(chem_coefs[chem_coefs != 0])
      names(selected)<-outcome[i] 
      selected$Feature<-row.names(selected)
      EN_feature<-merge(EN_feature,selected,by='Feature',all=T)
      
  }
  return(EN_feature)
}




