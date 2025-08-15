# asymptotically anytime valid f test

BayesFactor <- function(Y,
                        X,
                        Z,
                        varhat = 1, # pre experimental estimate of the variance
                        mde = 1, # minimum detectable effect
                        tap = 0.5 # treatment assignment probability
){
  stopifnot(nrow(X) == nrow(Z))
  n <- nrow(Z)
  d <- ncol(Z)
  
  
  lambda <- varhat/ (mde^2 * tap * (1-tap)) #how to choose this?
  
  model1 <- lm(Y ~ X + Z)
  model2 <- lm(Y ~ X)
  fstat <- anova(model1, model2)
  fstat <- fstat[2, "F"]
  
  
  return( (lambda/(lambda + n)^(d/2)*exp(d/2 * n/(n+lambda) * fstat)) )
  
}

# E need to be a dummy variable for in or not in environment
# returns bayes factors / e values for every time point from start
avChowBayesFactors <- function(X,Y,E,varhat = 1, mde = 1, tap = 0.5,start = 10){
  
  
  X <- as.matrix(X)
  Z <- factor(levels = unique(E))
  
  # define dummies for each E except for the first one
  # The first E is baseline
  for(e in unique(E)[-1]){
    Z <- cbind(Z,E == e, sweep(X,1,E == e,"*"))
  }
  
  
  
  results <- c()
  for(i in start:nrow(X)){
    
    Ytest <- Y[1:i]
    Ztest <- Z[1:i,]
    
    if(ncol(X) > 1){
      Xtest <- X[1:i,]
    }else{
      Xtest <- X[1:i]
    }
    
    
    results <- c(results,BayesFactor(Ytest,Xtest,Ztest,varhat, mde, tap))
  }
  
  return(results)
}

#avchow using lindon 2024

avChowICP <- function(dat, 
                      target_name = "Y", 
                      env_name = "E", 
                      m = 10, 
                      alpha = 0.05 # alpha/|E|?
){
  
  n <- as.numeric(nrow(dat))  
  predictor_names <- colnames(dat)[! colnames(dat) %in% c(target_name,env_name)]
  
  
  Y <- dat[,target_name]
  E <- dat[,env_name]
  X <- dat[,predictor_names]
  
  # for making empty list with subset names
  all_subsets <- unlist(lapply(0:length(predictor_names), function(k) {
    combn(predictor_names, k, simplify = FALSE)
  }), recursive = FALSE)
  
  subset_keys <- sapply(all_subsets, function(vars) {
    if (length(vars) == 0) "(intercept)" else paste(vars, collapse = "_")
  })
  
  
  myresult <- c()
  
  rescolnamearr= c()
  # loop over all subsets
  for (k in 1:length(predictor_names)) {
    combos <- combn(predictor_names, k, simplify = FALSE)
    for (vars in combos) {
      
      key <- paste(vars, collapse = "_")
      
      print(key)
      
      
      rescolnamearr <- c(rescolnamearr, key)
      resultmat <- c()

      result <- avChowBayesFactors(X[,vars], Y, E,
                                   #varhat = ?, mde = ?, tap = ?,
                                   start = m) 
      
      # still need to figure out how to get 
      # pre experimental variance
      # minimum detectable effect
      # treatment assignment probability
      # or standardized effect size
      
      resultmat <- cbind(resultmat, result)
      

      
      # result is a matrix multiple e processes
      # cumulative maximum for every environment
      # then choose the maximum for each time (for this one should use a bonferroni)
      # WRONG do not take cumulative maximum
      
      myresult <- cbind(myresult , apply(resultmat, 1, max))
    }
  }
  
  colnames(myresult) <- rescolnamearr
  
  
  return(myresult)
  
}

#' Compute intersection of covariate subsets where null is NOT rejected
#'
#' @param pvals A named vector of p-values. Names are comma-separated variable indices (e.g., "1,2,3").
#' @param alpha Significance level for hypothesis testing.
#' @return A vector of variable indices in the intersection set.
compute_icp_intersection <- function(evals, alpha = 0.05) {
  # Filter to subsets where null is NOT rejected
  non_rejected <- names(evals)[evals >= 1/alpha]
  
  if (length(non_rejected) == 0) {
    return(character(0))  # empty set
  }
  
  # Convert subset names into list of integer vectors
  sets <- lapply(non_rejected, function(name) unlist(strsplit(name, "_")))
  
  # Compute intersection
  intersection <- Reduce(intersect, sets)
  
  return(sort(intersection))
}

#' get_icp_set
#'
#' Computes the intersection set of predictors accepted by Invariant Causal Prediction (ICP)
#' at a given evaluation time, based on minimum p-values across time steps.
#'
#' @param pvals A matrix of p-values where each column corresponds to a predictor
#'              and each row corresponds to a time step or evaluation iteration.
#' @param alpha Significance level for the ICP test (default is 0.05).
#' @param eval_at Optional integer. The number of rows (time steps) to consider when evaluating.
#'                If NULL, defaults to the last row of the matrix.
#'
#' @return A character vector or index vector of variables accepted at level `alpha`.

get_icp_set <- function(evals, alpha = 0.05, eval_at = NULL){
  # Default to the last row if time is not given
  if (is.null(eval_at)) eval_at <- nrow(evals)
  
  # Safety check
  if (eval_at > nrow(evals)) stop("`eval_at` exceeds number of rows in evals.")
  
  evec <- evals[1:eval_at, ]
  
  return(compute_icp_intersection(evec, alpha))
}


# intersection of strings for last step of icp
intersect_all_strings <- function(str_vec) {
  # Handle NA or empty input

  if (length(str_vec) == 0 || any(is.na(str_vec))) return("")
  
  # Split and intersect

  split_list <- strsplit(str_vec, "_")

  common <- split_list[[1]]
  # print(common)
  
  if(length(split_list) >= 2){
    for (i in 2:length(split_list)) {
      common <- intersect(common, split_list[[i]])
      if (length(common) == 0) return("")
    }
  }
  
  
  return(paste(common, collapse = "_"))
}

# I have a array of strings A = ("X1", "X2",...,"X5") and a string of the form S="X1_X3". 
# I want to return the subtraction of these two sets. e.g. A-S = ("X2","X4","X5")
subtract_sets <- function(A, S) {
  if (is.na(S) || S == "") return(A)  # nothing to subtract
  S_parts <- unlist(strsplit(S, "_"))
  return(setdiff(A, S_parts))
}


#method in lindon 2025

asympAvF <- function(Y,X,Z,g=1){
  stopifnot(nrow(X) == nrow(Z))
  
  X <- as.matrix(X)
  
  n <- nrow(Z)
  p <- ncol(X)
  d <- ncol(Z)
  nu <- n-p-d
  
  model1 <- lm(Y ~ X + Z)
  model2 <- lm(Y ~ X)
  fstat <- anova(model1, model2)
  fstat <- fstat[2, "F"]
  

  
  return( (g/(g+n))^(d/2) * ( (1+g/(g+n) * d/nu * fstat) / (1+ d/nu * fstat)  )^(-(nu+d)/2) )
}

# intercept model if X is empty
asympAvFempty <- function(Y,Z,g=1){
  stopifnot(length(Y) == nrow(Z))
  

  
  n <- length(Z)
  d <- 1
  nu <- n-d
  
  model1 <- lm(Y ~ Z)
  model2 <- lm(Y ~ 1)
  
  # print(summary(model1))
  # print(summary(model2))
  
  fstat <- anova(model1, model2)
  fstat <- fstat[2, "F"]
  
  
  
  return( (g/(g+n))^(d/2) * ( (1+g/(g+n) * d/nu * fstat) / (1+ d/nu * fstat)  )^(-(nu+d)/2) )
}

asympAvFChow <- function(X,Y,E,g = 1,start = 10){
  
  
  X <- as.matrix(X)
  Z <- factor(levels = unique(E))
  
  # define dummies for each E except for the first one
  # The first E is baseline
  for(e in unique(E)[-1]){
    Z <- cbind(Z,E == e, sweep(X,1,E == e,"*"))
  }
  
  

  results <- c()
  for(i in start:nrow(X)){
    
    Ytest <- Y[1:i]
    Ztest <- Z[1:i,]
    
    if(ncol(X) > 1){
      Xtest <- X[1:i,]
    }else{
      Xtest <- X[1:i]
    }
    
    
    results <- c(results,asympAvF(Ytest,Xtest,Ztest,g))
  }
  
  return(results)
}

# avChow if X is empty
asympAvFChowempty <- function(Y,E,g = 1,start = 10){
  

  Z <- as.factor(E)
  
  results <- c()
  for(i in start:length(Y)){
    
    Ytest <- Y[1:i]
    Ztest <- Z[1:i]
    
    
    results <- c(results,asympAvFempty(Ytest,Ztest,g))
  }
  
  return(results)
}



# avChowICP using g prior based method

avChowICP2 <- function(dat, 
                       target_name = "Y", 
                       env_name = "E", 
                       m = 10, # start
                       g = 1, # hyperparameter
                       alpha = 0.05 # alpha/|E|?
){
  
  n <- as.numeric(nrow(dat))  
  predictor_names <- colnames(dat)[! colnames(dat) %in% c(target_name,env_name)]
  
  
  Y <- dat[,target_name]
  E <- dat[,env_name]
  X <- dat[,predictor_names]
  
  # for making empty list with subset names
  all_subsets <- unlist(lapply(0:length(predictor_names), function(k) {
    combn(predictor_names, k, simplify = FALSE)
  }), recursive = FALSE)
  
  subset_keys <- sapply(all_subsets, function(vars) {
    if (length(vars) == 0) "" else paste(vars, collapse = "_")
  })
  
  
  myresult <- c()
  
  rescolnamearr= c()
  # loop over all subsets
  for (k in 0:length(predictor_names)) {
    combos <- combn(predictor_names, k, simplify = FALSE)
    for (vars in combos) {
      
      if(k == 0){
        key <- ""
        print(key)
        
        rescolnamearr <- c(rescolnamearr, key)
        resultmat <- c()
        
        result <- asympAvFChowempty(Y, E, g = g, start = m) 
        resultmat <- cbind(resultmat, result)
        myresult <- cbind(myresult , apply(resultmat, 1, max))
        
      }else{
        
        key <- paste(vars, collapse = "_")
        
        print(key)
        
        
        rescolnamearr <- c(rescolnamearr, key)
        resultmat <- c()
        
        
        result <- asympAvFChow(X[,vars], Y, E, g = g, start = m) 
        
        # still need to figure out how to get 
        # pre experimental variance
        # minimum detectable effect
        # treatment assignment probability
        # or standardized effect size
        
        resultmat <- cbind(resultmat, result)
        
        
        
        # result is a matrix multiple e processes
        # cumulative maximum for every environment
        # then choose the maximum for each time (for this one should use a bonferroni)
        # WRONG do not take cumulative maximum
        
        myresult <- cbind(myresult , apply(resultmat, 1, max))
        
      }
      
      
    }
  }
  
  colnames(myresult) <- rescolnamearr
  
  
  return(myresult)
  
}

#-----------------------------------------
# active learning

avChowActiveICP <- function(dat, 
                       target_name = "Y", 
                       env_name = "E", 
                       m = 10, # start
                       g = 1, # hyperparameter
                       alpha = 0.05 # alpha/|E|?
){
  
  n <- as.numeric(nrow(dat))  
  predictor_names <- colnames(dat)[! colnames(dat) %in% c(target_name,env_name)]
  
  
  Y <- dat[,target_name]
  E <- dat[,env_name]
  X <- dat[,predictor_names]
  
  # for making empty list with subset names
  all_subsets <- unlist(lapply(0:length(predictor_names), function(k) {
    combn(predictor_names, k, simplify = FALSE)
  }), recursive = FALSE)
  
  subset_keys <- sapply(all_subsets, function(vars) {
    if (length(vars) == 0) "" else paste(vars, collapse = "_")
  })
  
  
  myresult <- c()
  
  rescolnamearr= c()
  # loop over all subsets
  for (k in 0:length(predictor_names)) {
    combos <- combn(predictor_names, k, simplify = FALSE)
    for (vars in combos) {
      
      if(k == 0){
        key <- ""
        print(key)
        
        rescolnamearr <- c(rescolnamearr, key)
        resultmat <- c()
        
        result <- asympAvFChowempty(Y, E, g = g, start = m) 
        resultmat <- cbind(resultmat, result)
        myresult <- cbind(myresult , apply(resultmat, 1, max))
        
      }else{
        
        key <- paste(vars, collapse = "_")
        
        print(key)
        
        
        rescolnamearr <- c(rescolnamearr, key)
        resultmat <- c()
        
        
        result <- asympAvFChow(X[,vars], Y, E, g = g, start = m) 
        
        # still need to figure out how to get 
        # pre experimental variance
        # minimum detectable effect
        # treatment assignment probability
        # or standardized effect size
        
        resultmat <- cbind(resultmat, result)
        
        
        
        # result is a matrix multiple e processes
        # cumulative maximum for every environment
        # then choose the maximum for each time (for this one should use a bonferroni)
        # WRONG do not take cumulative maximum
        
        myresult <- cbind(myresult , apply(resultmat, 1, max))
        
      }
      
      
    }
  }
  
  colnames(myresult) <- rescolnamearr
  
  
  return(myresult)
  
}

jaccard_similarity <- function(result, parent){
  
  if(result == "" && parent == "") return(1)
  if(result != "" && parent == "" ) return(0)

  
  # Split and intersect
  
  res_split_list <- strsplit(result, "_")[[1]]
  par_split_list <- strsplit(parent, "_")[[1]]
  
  numerator <- length(intersect(res_split_list, par_split_list))
  #print(numerator)
  denominator <- length(union(res_split_list, par_split_list))
  #print(denominator)
  
  return(numerator/denominator)
  
}

is_subset <- function(X, S) {
  all(X %in% S | is.na(X))
}