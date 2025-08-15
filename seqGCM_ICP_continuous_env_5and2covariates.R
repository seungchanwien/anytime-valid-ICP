# experiment for 2 covariates, as well as for 5 covariates.

is_subset <- function(X, S) {
  all(X %in% S | is.na(X))
}


library("InvariantCausalPrediction")
source("C:/ETH/Thesis/Code/seqGCM related/seqICP_functions.R")
setwd("C:/ETH/Thesis/plots/1000sample_seqGCM_ICP_E_cont")


outer_sim <- 100

inner_sim <- 500
m <- 50
res <- c()

# Initialize result containers
predicted_sets <- vector("list", outer_sim)
GCMmat_list <- vector("list", outer_sim)
pval_mat_list <- vector("list", outer_sim)


target <- "Y"
myenv <- "E"
for(i in 1:outer_sim){
    # source("C:/ETH/Thesis/Code/seqICP_functions.R")
    set.seed(i*24)
    
    print(paste0("outer_sim: ", i))# update progress

    # E <- rnorm(inner_sim*2)
    # X1 <- E + E^2 / 2 + rnorm(inner_sim*2)
    # Y <- 3 * X1 + rnorm(inner_sim*2)
    # X2 <- 0.5 * E - 2* E ^ 2 - Y + rnorm(inner_sim*2)
    # data.obs <- data.frame(cbind(Y,E,X1,X2))
    # n <- as.numeric(nrow(data.obs)/2)
    # predictor_names <- colnames(data.obs)[! colnames(data.obs) %in% c("Y","E")]

    E <- rnorm(inner_sim*2)
    X1 <- E + rnorm(inner_sim*2)
    X2 <- E + X1 + rnorm(inner_sim*2)
    X3 <- E - X1 + rnorm(inner_sim*2)
    Y <- 2 * X2 - 3 * X3 + rnorm(inner_sim*2)
    X4 <- E + X3 - 3 * Y + rnorm(inner_sim*2)
    X5 <- E - Y + rnorm(inner_sim*2)
    data.obs <- data.frame(cbind(Y,E,X1,X2,X3,X4,X5))
    n <- as.numeric(nrow(data.obs)/2)
    predictor_names <- colnames(data.obs)[! colnames(data.obs) %in% c(target, myenv)]
    
    X <- as.matrix(data.obs[, predictor_names])
    
    # print(colnames(data.obs))
    seqgcm_icp <- seqGCM_ICP(data.obs,
                             target_name = target,
                             m = m,
                             withplot = FALSE,  # turn off plotting in parallel
                             E_continuous = TRUE)
    
    pred_set <- if (length(seqgcm_icp$predicted_set) == 0) {
      "(intercept)"
    } else {
      paste(seqgcm_icp$predicted_set, collapse = "_")
    }
    
    # list(
    #   pred_set = pred_set,
    #   GCMmat = seqgcm_icp$GCMmat,
    #   pval_mat = seqgcm_icp$pval_mat
    # )


  predicted_sets[[i]] <- pred_set
  GCMmat_list[[i]] <- seqgcm_icp$GCMmat
  pval_mat_list[[i]] <-seqgcm_icp$pval_mat



  
  # # now print your results
  # print(unlist(predicted_sets))
  
  
  # Define a color vector
  colors <- 1:ncol(pval_mat_list[[1]])  # or use rainbow(), topo.colors(), etc.

  
  matplot(
    pval_mat_list[[i]],          # all columns
    type = "l",        # line plot
    lty = 1,           # lines
    col = colors,
    xlab = "Index",    # custom X-axis label
    ylab = "p-value",
    main = paste0("P-values over Index, nr: ", i)
  )
  
  abline(h = 0.05, col = "red")
  # Optional: add a legend
  legend(
    "topright",
    legend = colnames(pval_mat_list[[i]]),
    col = colors,
    lty = 1
  )
  
  
  matplot(
    GCMmat_list[[i]],          # all columns
    type = "l",        # line plot
    lty = 1,           # solid lines
    col = colors,
    xlab = "Index",    # custom X-axis label
    ylab = "GCM values",
    main = paste0("GCM values, nr: ", i)
  )
  
  # Optional: add a legend
  legend(
    "topright",
    legend = colnames(GCMmat_list[[i]]),
    col = colors,
    lty = 1
  )
}

unlist(predicted_sets)
table(unlist(predicted_sets))
tail(pval_mat_list[[1]])

# check anytime validity
icpset_over_time <- vector("list", outer_sim)
result_is_subset <- vector("list", outer_sim)
final_result_is_subset <- rep(NA, outer_sim)

for(j in 1:outer_sim){
  
  result <- vector("list", outer_sim)
  
  for(i in m:inner_sim){
    
    result[[i]] <- get_icp_set(pval_mat_list[[j]], eval_at = i)
    
  }
  
  # 2 vars: S = c("X1", ""), 5 vars: S = c("X2","X3", "")
  result_is_subset[[j]] <- sapply(result, is_subset, S = c("X2","X3", ""))
  final_result_is_subset[[j]] <- all(result_is_subset[[j]])
  
  
  collapsed_icpset <- sapply(result, function(x) {
    if (length(x) == 0 || all(is.na(x))) {
      return(NA)  # or "" if you prefer
    } else {
      paste(x, collapse = "_")
    }
  })
  
  if(i == 1) print(collapsed_icpset)
  
  icpset_over_time[[j]] <- collapsed_icpset
  
  
}

which(icpset_over_time[[81]]== "X1")
result_is_subset
which(final_result_is_subset == FALSE)
table(final_result_is_subset)

# Replace NA with "NA" string and keep "" as is
x_clean <- as.character(icpset_over_time[[81]])
x_clean[is.na(x_clean)] <- "NA"
x_clean[x_clean == ""] <- "empty"  # Optional: name empty strings

# Create a data frame with index
df <- data.frame(index = 1:length(x_clean), value = x_clean)
library(ggplot2)

ggplot(df, aes(x = index, y = value)) +
  geom_point(aes(color = value), shape = 15, size = 0.5) +
  scale_y_discrete(limits = unique(df$value)) +
  theme_minimal() +
  labs(title = "Categorical Values Over Time", x = "Index", y = "Value")
















# saving the plots

for(i in 1:outer_sim){
  # Define a color vector
  colors <- 1:ncol(pval_mat_list[[1]])  # or use rainbow(), topo.colors(), etc.
  
  
  plot_name <- paste0("pval_nr", i, ".jpeg")
  jpeg(file=plot_name, width = 1440, height = 1080)
  matplot(
    pval_mat_list[[i]],          # all columns
    type = "l",        # line plot
    lty = 1,           # lines
    col = colors,
    xlab = "Index",    # custom X-axis label
    ylab = "p-value",
    main = paste0("P-values over Index, nr: ", i)
  )
  
  abline(h = 0.05, col = "red")
  # Optional: add a legend
  legend(
    "topright",
    legend = colnames(pval_mat_list[[i]]),
    col = colors,
    lty = 1
  )
  dev.off()
  
  plot_name <- paste0("gcmvals_nr", i, ".jpeg")
  jpeg(file=plot_name, width = 1440, height = 1080)
  matplot(
    GCMmat_list[[i]],          # all columns
    type = "l",        # line plot
    lty = 1,           # solid lines
    col = colors,
    xlab = "Index",    # custom X-axis label
    ylab = "GCM values",
    main = paste0("GCM values, nr: ", i)
  )
  
  # Optional: add a legend
  legend(
    "topright",
    legend = colnames(GCMmat_list[[i]]),
    col = colors,
    lty = 1
  )
  dev.off()
}

