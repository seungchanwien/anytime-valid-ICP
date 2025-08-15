library(InvariantCausalPrediction)



setwd("C:/ETH/Thesis/Code/lindon related")
source("anytimeValidFTest_utils3.R")

# -----------------------------
# Outer simulation (corrected)
# -----------------------------
alpha <- 0.05
outer <- 300
inner <- 200
start <- 50

parentlist <- character(outer)

final_result_is_subset        <- logical(outer)  # avICP (anytime) validity over time
final_fixed_result_is_subset  <- logical(outer)  # fixed ICP validity over time

icp_set_over_time_list        <- vector("list", outer)
fixed_icp_set_over_time_list  <- vector("list", outer)
result_is_subset_list         <- vector("list", outer)
fixed_result_is_subset_list   <- vector("list", outer)
pvals_over_time_list          <- vector("list", outer)
fixed_jaccard_similarity_over_time<- vector("list", outer)

for (o in 31:outer) {
  
  # data loading
  {
    filename <- paste0("C:/ETH/Thesis/Code/sempler/new_data_", o - 1, ".csv")
    dagname  <- paste0("C:/ETH/Thesis/Code/sempler/new_dag_",  o - 1, ".txt")
    
    scm_data <- read.csv(filename)
    scm_data <- as.data.frame(scm_data)
    names(scm_data)[1] <- "Y"
    names(scm_data)[length(scm_data)] <- "E"
    scm_data[,"E"] <- as.factor(scm_data[,"E"])
    
    dag <- read.csv(dagname, header = FALSE)
    rownames(dag) <- names(scm_data[-length(scm_data)])
    colnames(dag) <- names(scm_data[-length(scm_data)])
    parents <- rownames(dag)[dag[, 1] != 0]
    parents <- paste(parents, collapse = "_")
    parentlist[o] <- parents
    
    target_name    <- "Y"
    env_name       <- "E"
    predictor_names <- colnames(scm_data)[! colnames(scm_data) %in% c(target_name, env_name)]
    
    # Randomize order, but vary across outer iterations
    set.seed(14 + o)
    idx <- sample(1:nrow(scm_data), size = nrow(scm_data))
    scm_data <- scm_data[idx, ]
    scm_data <- scm_data[1:inner, ]
  }
  
  # fixed sample ICP
  {
    # Fixed-sample ICP over time (for comparison)
    icpset_fixed <- rep(NA,nrow(scm_data))
    
    # dynamic starting point to ensure functionality of ICP
    {
      E_full <- as.factor(scm_data$E)
      # position of the 3rd occurrence of each environment (in the randomized order)
      pos3 <- sapply(levels(E_full), function(lvl) {
        idx <- which(E_full == lvl)
        if (length(idx) >= 3) idx[3] else NA_integer_
      })
      if (any(is.na(pos3))) {
        stop("Not enough observations (≥3) for at least one environment in the first 'inner' rows.")
      }
      start_dyn <- max(pos3, na.rm = TRUE)
    }
    
   
    # --- NEW: preallocate a matrix to store p-values over time for this run ---
    Tseq <- 1:nrow(scm_data)
    p <- length(predictor_names)
    pval_mat <- matrix(NA_real_, nrow = nrow(scm_data), ncol = p,
                       dimnames = list(time = Tseq, var = predictor_names))
    
    for (t in start_dyn:nrow(scm_data)) {
      X <- as.matrix(scm_data[1:t, predictor_names])
      E <- as.factor(scm_data$E)[1:t]
      Y <- scm_data$Y[1:t]
      
      fit <- ICP(X, Y, E, alpha = alpha)
      
      # intersection set for reporting (unchanged)
      if (length(fit$acceptedSets)) {
        S_hat_idx_used <- Reduce(intersect, fit$acceptedSets)
      } else {
        S_hat_idx_used <- integer(0)
      }
      S_hat_idx_orig  <- if (!is.null(fit$usedvariables)) fit$usedvariables[S_hat_idx_used] else S_hat_idx_used
      S_hat_names     <- colnames(X)[S_hat_idx_orig]
      
      icpset_fixed[t] <-  paste(S_hat_names, collapse = "_")
      
      # --- NEW: extract per-variable p-values (mapped back to original columns) ---
     
      pval_mat[t, ] <- fit$pvalues
    }
    
  }
    
  # plot fixed sample ICP
  {
    # --- NEW: save matrix for this outer run ---
    pvals_over_time_list[[o]] <- pval_mat
    
    # --- NEW: plot progression of p-values for this run ---
    plot_name <- paste0("ICP_pvalues_sim_", o, ".jpeg")
    jpeg(file = plot_name, width = 1440, height = 1080)
    matplot(Tseq, pval_mat, type = "l", lty = 1, xlab = "Time (t)",
            ylab = "p-value", main = paste0("Fixed-n ICP p-values over time (run ", o, ")"))
    abline(h = alpha, col = "red", lwd = 2)   # significance threshold
    legend("topright", legend = colnames(pval_mat), col = seq_len(ncol(pval_mat)), lty = 1)
    dev.off()
    
    # (optional) also write the matrix to CSV for later analysis
    write.csv(pval_mat, file = paste0("ICP_pvalues_sim_", o, ".csv"), row.names = TRUE)
    
    print(o)
  }
    
  
  
  # avICP via avChow
  {
    result <- avChowICP2(
      scm_data,
      target_name = "Y",
      env_name    = "E",
      g           = 1,
      m           = start,
      alpha       = alpha
    )
    
  }
  
  
  # Plot e-values (log scale)
  {
    colors <- seq_len(ncol(result))
    plot_name <- paste0("avICP_avChow_sim_", o, ".jpeg")
    jpeg(file = plot_name, width = 1440, height = 1080)
    matplot(log(result),
            type = "l",
            lty  = 1,
            col  = colors,
            xlab = "Time",
            ylab = "log e-value",
            main = paste0("avICP-avChow: E-values over time (run ", o, ")"))
    abline(h = log(1/alpha), col = "red")
    legend("topright", legend = colnames(result), col = colors, lty = 1)
    dev.off()
  }
  
  
  # Anytime-valid ICP set over time from e-values
  {
    notreject <- vector("list", nrow(result))
    icpset_over_time <- vector("list", nrow(result))
    icpset_over_time_strsplit <- vector("list", nrow(result))
    fixed_icpset_over_time_strsplit <- vector("list", nrow(result))
    
    for (k in 1:nrow(result)) {
      notreject[[k]] <- colnames(result)[ result[k, ] < 1/alpha ]
      icpset_over_time[[k]] <- intersect_all_strings(notreject[[k]])
      icpset_over_time_strsplit[[k]] <- unlist(strsplit(icpset_over_time[[k]], "_"))
      fixed_icpset_over_time_strsplit[[k]] <- unlist(strsplit(icpset_fixed[[k]], "_"))
    }
    
    icp_set_over_time_list[[o]] <- unlist(icpset_over_time)
    
    Sstar <- unlist(strsplit(parents, split = "_"))
    result_is_subset <- sapply(icpset_over_time_strsplit, is_subset, S = c(Sstar, ""))
    fixed_result_is_subset <- sapply(fixed_icpset_over_time_strsplit, is_subset, S = c(Sstar, ""))
    
    result_is_subset_list[[o]] <- result_is_subset
    fixed_result_is_subset_list[[o]] <- fixed_result_is_subset
    
    final_result_is_subset[o]       <- all(result_is_subset)
    final_fixed_result_is_subset[o] <- all(fixed_result_is_subset)
  }
  
}

 # Post-analysis examples
table(final_result_is_subset)
last_elements <- sapply(icp_set_over_time_list, function(x) x[length(x)])
table(last_elements == parentlist)
table(final_fixed_result_is_subset)
which(final_fixed_result_is_subset == FALSE)
which(parentlist=="")
parentlist
below_alpha_list <- lapply(pvals_over_time_list, function(mat) {
  mat < 0.05
})

# get fixed_icp_set_over_time_list
for(i in 1:outer){
  fixed_icp_set_over_time_list[[i]]<- apply(below_alpha_list[[i]],
                                            MARGIN = 1,
                                             function(row){
                                               if(any(is.na(colnames(below_alpha_list[[i]])[row]))){
                                                 #print(is.na(colnames(below_alpha_list[[i]])[row]))
                                                 icpset <- ""
                                               }else{
                                                 icpset <- paste(colnames(below_alpha_list[[i]])[row],collapse= "_")
                                               }
                                                return(icpset)
                                                 
                                             }
                                            )
  
  fixed_jaccard_similarity_over_time[[i]] <-  sapply(fixed_icp_set_over_time_list[[i]], jaccard_similarity, parentlist[i])
}

fixed_icp_last_entries <- sapply(fixed_jaccard_similarity_over_time, function(x) x[length(x)])
# return jaccard similarity over time

table(fixed_icp_last_entries)
