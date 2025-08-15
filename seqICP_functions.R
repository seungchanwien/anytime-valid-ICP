#library(mgcv)
#library(xgboost)


# Robbins-Siegmund distribution
RS <- function(r){
  1 - 2 * (1 - pnorm(sqrt(r)) + sqrt(r) * dnorm(sqrt(r)))
}

# GCM by peters
GCM <- function(prodres){
  n <- length(prodres)
  sqrt(n)*mean(prodres)/sd(prodres)
}

seqGCM <- function(prodres){
  mean(prodres)/sd(prodres)
}

#' GCM_prodres
#'
#' Computes the product of residuals from two models—one predicting a target variable Y 
#' and the other predicting an environment variable E—across all subsets of predictors.
#' This is designed for use in generalized covariance measure (GCM) style tests for 
#' conditional independence between Y and E given subsets of predictor variables.
#'
#' @param dat.trn A data.frame containing the training data used to fit the models.
#' @param dat.tst A data.frame containing the test data used to compute residuals.
#' @param target_name Character string: the name of the target variable (default: "Y").
#' @param env_name Character string: the name of the environment variable (default: "E").
#'
#' @return A named list. Each element corresponds to a model (subset of predictors),
#'         and contains the element-wise product of residuals from the Y and E models
#'         fitted on the training set and evaluated on the test set.
#'
#' @examples
#' result <- GCM_prodres(train_data, test_data)
#' head(result[["X1_X2"]])  # product of residuals for predictors X1 and X2
#'
#' @details
#' - Models for Y are fitted using linear regression (`lm`).
#' - Models for E are fitted using logistic regression (`glm` with logit link).
#' - Probabilities are clipped to [1e-6, 1 - 1e-6] to avoid numerical issues.
#' - Residuals from E are standardized.
#' - All predictor combinations (up to full set) are considered.

GCM_prodres <- function(dat.trn, dat.tst, target_name = "Y", env_name = "E", E_continuous = FALSE){
  
  predictor_names <- colnames(data.obs)[! colnames(data.obs) %in% c(target_name,env_name)]
  residual_list_Y <- list()
  residual_list_E <- list()
  prod_res_list <- list()
  # --- Empty set (intercept-only model) ---
  # get residuals of for model on Y
  formulaY <- as.formula(paste(target_name,"~ 1"))
  modelY_0 <- lm(formulaY, data = dat.trn)

  residY_0 <- dat.tst[,target_name] - predict(modelY_0, newdata = dat.tst)
  
  
  # get residuals for model on E
  formulaE <- as.formula(paste(env_name,"~ 1"))
  
  if(E_continuous){
    modelE_0 <- lm(formulaE, data = dat.trn) # no gam needed
    residE_0 <- dat.tst[,env_name] - predict(modelE_0, newdata = dat.tst)

  }else{
    formulaE <- as.formula(paste(env_name,"~ 1"))
    modelE_0 <- gam(formulaE, family = binomial(link = "logit"), data = dat.trn) # should be improved with different classification method
    # modelE_0 <- multinom(formulaE, data = dat.trn) # intercept model for multiclass

    
    pred_probs_0 <- predict(modelE_0, newdata = dat.tst, type = "response")
    pred_probs_0 <- pmin(pmax(pred_probs_0, 1e-6), 1 - 1e-6)
    residE_0 <- dat.tst[,env_name] - pred_probs_0 / sqrt(pred_probs_0 * (1 - pred_probs_0))
  }

  prod_res_list[["(intercept)"]] <- residY_0 * residE_0

  
  # Loop over all subset sizes
  for (k in 1:length(predictor_names)) {
    combos <- combn(predictor_names, k, simplify = FALSE)
    for (vars in combos) {
      
      # get residuals for model on Y and model on E
      
      
      # for linear models use this:
      formulaY <- as.formula(paste0(target_name,"~", paste(vars, collapse = " + ")))
      formulaE <- as.formula(paste(env_name, "~", paste(vars, collapse = " + ")))
      
      # for smooth terms e.g. gam model, use this:
      # formulaY <- as.formula(paste0(target_name, " ~ ", paste0("s(", vars, ")", collapse = " + ")))
      # formulaE <- as.formula(paste0(env_name, " ~ ", paste0("s(", vars, ")", collapse = " + ")))
      
      # residuals for model on Y
      modelY <- lm(formulaY, data = dat.trn)
      residY <- dat.tst[,target_name] - predict(modelY, newdata = dat.tst)
      
      
      
      if(E_continuous){
        modelE <- lm(formulaE, data = dat.trn) # continuous environment predicted with GAM, can be improved/ made flexible 
        residE <- dat.tst[,env_name] - predict(modelE, newdata = dat.tst)
      }else{
        formulaE <- as.formula(paste(env_name, "~", paste(vars, collapse = " + ")))
        suppressWarnings({
          modelE <- gam(formulaE, family = binomial(link = "logit"), data = dat.trn) # should be improved with different classification method
        })
        
        # dtrain <- xgb.DMatrix(data = dat.trn[,vars], label = dat.trn[,env_name])
        # # For multiclass: use objective = "multi:softprob" and set num_class = N.
        # modelE <- xgboost(data = dtrain, objective = "binary:logistic", nrounds = 50, verbose = 0) 
        
        pred_probs <- predict(modelE, newdata = dat.tst)
        pred_probs <- pmin(pmax(pred_probs, 1e-6), 1 - 1e-6)
        residE <- dat.tst[,env_name] - pred_probs/ sqrt(pred_probs * (1 - pred_probs))
      }
      
      
      # calculate residuals
      
      
      key <- paste(vars, collapse = "_")
      
      residual_list_Y[[key]] <- residY
      residual_list_E[[key]] <- residE
      prod_res_list[[key]] <- residY * residE
      
    }
    
  }
  
  return(prod_res_list)
  
}

# @param residY residuals wrt target
# @param residY residuals wrt environment
# @param m start time 
# @param k is the current iteration

seqGCM_pval <- function(GCMval, m, k){
  pval <- 1 - RS(max(c(k * GCMval^2 - log(k/m), 0 )))
}


#' Compute intersection of covariate subsets where null is NOT rejected
#'
#' @param pvals A named vector of p-values. Names are comma-separated variable indices (e.g., "1,2,3").
#' @param alpha Significance level for hypothesis testing.
#' @return A vector of variable indices in the intersection set.
compute_icp_intersection <- function(pvals, alpha = 0.05) {
  # Filter to subsets where null is NOT rejected
  non_rejected <- names(pvals)[pvals >= alpha]
  
  if (length(non_rejected) == 0) {
    return("")  # empty set
  }
  
  # Convert subset names into list of integer vectors
  sets <- lapply(non_rejected, function(name) unlist(strsplit(name, "_")))

  # Compute intersection
  intersection <- Reduce(intersect, sets)

  if(length(intersection) == 0){
    return("")
  }
  
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

get_icp_set <- function(pvals, alpha = 0.05, eval_at = NULL){
  # Default to the last row if time is not given
  if (is.null(eval_at)){ eval_at <- nrow(pvals)}
  
  # Safety check
  if (eval_at > nrow(pvals)) stop("`eval_at` exceeds number of rows in pvals.")
  
  pvec <- pvals[eval_at, ]
  
  return(compute_icp_intersection(pvec, alpha))
}

GCMmat2Pvalmat <- function(GCMmat, m){

  
  pval_mat <- matrix(NA, nrow = nrow(GCMmat), ncol = ncol(GCMmat))
  
  for (k in seq_len(nrow(GCMmat))) {
    if (k >= m) {
      for (j in seq_len(ncol(GCMmat))) {
        GCMval <- GCMmat[k, j]
        pval_mat[k, j] <- seqGCM_pval(GCMval, m, k)
      }
    }
  }
  
  colnames(pval_mat) <- colnames(GCMmat)
  
  return(pval_mat)
}

#' seqGCM_ICP
#'
#' Performs asymptotically anytime valid Invariant Causal Prediction using the sequential Generalized Covariance Measure (GCM).
#' Iteratively computes GCM-based test statistics and p-values across time, comparing residuals
#' between models for the target and environment variables.
#'
#' @param dat A data.frame containing the full dataset, including target, environment, and predictors.
#' @param target_name Name of the target variable (default: "Y").
#' @param env_name Name of the environment variable (default: "E").
#' @param m Starting index for sequential evaluation (default: 10).
#' @param eval_at Optional stopping index (e.g. to evaluate at a specific time).
#' @param withplot Logical. Whether to plot the p-values over time (default: FALSE).
#'
#' @return A vector of accepted variables (intersection set) at level `alpha`, evaluated up to `eval_at`.


seqGCM_ICP <- function(dat, target_name = "Y", env_name = "E", m = 10, times , alpha = 0.05, eval_at = NULL, withplot = FALSE, E_continuous = FALSE){
  n <- as.numeric(nrow(dat)/2)
  predictor_names <- colnames(dat)[! colnames(dat) %in% c(target_name,env_name)]
  
  
  # for making empty list with subset names
  all_subsets <- unlist(lapply(0:length(predictor_names), function(k) {
    combn(predictor_names, k, simplify = FALSE)
  }), recursive = FALSE)
  
  
  subset_keys <- sapply(all_subsets, function(vars) {
    if (length(vars) == 0) "(intercept)" else paste(vars, collapse = "_")
  })
  

  emptylist <- setNames(vector("list", length(all_subsets)), subset_keys)
  
  prodres_list <- emptylist
  residY_list <- emptylist
  residE_list <- emptylist
  GCMvec_list <- lapply(emptylist, function(x) rep(0, m-1))
  pval_list <- lapply(emptylist, function(x) rep(1, m-1))
    


    for(i in m:n){
      

      
      if(i %% 200 == 0) print( paste0("iteration number: ", i))
      
      (idx <- 1:i * 2 ) # even: training; odd: test
      dat.trn <- dat[idx,]
      
      
      if(i == m){
        dat.tst <- dat[idx -1,]
      }else{
        dat.tst <- dat[i*2-1,]
      }

        prodres_list <- mapply(c, 
                               prodres_list, 
                               GCM_prodres(dat.trn, dat.tst,
                                           target_name = target_name, env_name = env_name,
                                           E_continuous = E_continuous), 
                               SIMPLIFY = FALSE)

      
      tmp_GCM_list <- lapply(prodres_list, seqGCM)
      GCMvec_list <-  mapply(c, GCMvec_list, tmp_GCM_list, SIMPLIFY = FALSE)
      
      
    }
  
  
  

  GCMmat <- do.call(cbind, GCMvec_list)
  
  pval_mat <- GCMmat2Pvalmat(GCMmat, m)
  
  if(withplot){
    # Basic line plot with matplot:
    matplot(
      pval_mat,          # all columns
      type = "l",        # line plot
      lty = 1,           # solid lines
      xlab = "Index",    # custom X-axis label
      ylab = "p-value",
      main = "P-values over Index"
    )
    
    abline(h = 0.05, col = "red")
    # Optional: add a legend
    legend(
      "topright",
      legend = colnames(pval_mat),
      col = seq_len(ncol(pval_mat)),
      lty = 1
    )
  }

  
  return(list(predicted_set = get_icp_set(pval_mat, alpha, eval_at),
              GCMmat = GCMmat,
              pval_mat = pval_mat))
  
  
}