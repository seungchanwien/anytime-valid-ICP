is_subset <- function(X, S) {
  all(X %in% S | is.na(X))
}

setwd("C:/ETH/Thesis/Code/lindon related")
source("anytimeValidFTest_utils3.R") # has chow test with all dummy variables included


# test area#############################
set.seed(12)
n <- 1000
E <- sample(0:5, n, replace = TRUE)
X1 <- 3 * (E == 1) + rnorm(n)
X2 <- 3 * (E == 2) + X1 + rnorm(n)
X3 <- 3 * (E == 3) - X1 + rnorm(n)
Y <- 2 * X2 - 3 * X3 + rnorm(n)
X4 <- 3 * (E == 4) + X3 - 3 * Y + rnorm(n)
X5 <- 3 * (E == 5) - Y + rnorm(n)
data.obs <- data.frame(cbind(Y,E,X1,X2,X3,X4,X5))

target_name <- "Y"
env_name <- "E"
predictor_names <- colnames(data.obs)[! colnames(data.obs) %in% c(target_name,env_name)]


X <- data.obs[,predictor_names]
E <- as.factor(E)

# res <- avChowBayesFactors(X,Y,E, start = 30)
m <- 30
result <- avChowICP2(data.obs,
                     target_name = "Y",
                     env_name = "E", 
                     m = m, 
                     alpha = alpha)

alpha <- 0.05
# plot results
colors <- 1:ncol(result)
matplot(log(result), type = "l")
abline(h = log(1/alpha), col = "red")
legend(
  "topright",
  legend = colnames(result),
  col = colors,
  lty = 1
)

# which subsets are not rejected
alpha <- 0.05
notreject <- colnames(result)[result[nrow(result),]< 1/alpha]
icp_results <- intersect_all_strings(notreject)
notreject

# outer_sim <- 1
# 
# # check anytime validity
# icpset_over_time <- vector("list", outer_sim)

# final_result_is_subset <- rep(NA, outer_sim)


notreject <- vector("list", nrow(result))
icpset_over_time <- vector("list", nrow(result))
icpset_over_time_strsplit <- vector("list", nrow(result))
result_is_subset <- rep(NA, nrow(result))

nrow(result)
for(k in 1:nrow(result)){
  
  notreject[[k]] <- colnames(result)[result[k,]< 1/alpha]
  icpset_over_time[[k]] <- intersect_all_strings(notreject[[k]])
  icpset_over_time_strsplit[[k]] <- unlist(strsplit(icpset_over_time[[k]], "_"))
  
  
}

result_is_subset <- sapply(icpset_over_time_strsplit, is_subset, S = c("X2","X3", ""))

final_result_is_subset <- all(result_is_subset)
final_result_is_subset


# collapsed_icpset <- sapply(result, function(x) {
#   if (length(x) == 0 || all(is.na(x))) {
#     return(NA)  # or "" if you prefer
#   } else {
#     paste(x, collapse = "_")
#   }
# })
# 
# if(i == 1) print(collapsed_icpset)
# 
# icpset_over_time[[j]] <- collapsed_icpset







########################################
# test successful
#




outer_sim <- 100
n <- 500
m <- 30
alpha <- 0.05

icp_results <- rep(NA,times = outer_sim)
final_result_is_subset <- rep(NA,times = outer_sim)
icp_set_over_time_list <- vector("list", outer_sim)
result_is_subset_list <- vector("list", outer_sim)

for(i in 1:outer_sim){
  set.seed(7474+3*i)
  
  E <- sample(0:5, n, replace = TRUE)
  X1 <- 3 * (E == 1) + rnorm(n)
  X2 <- 3 * (E == 2) + X1 + rnorm(n)
  X3 <- 3 * (E == 3) - X1 + rnorm(n)
  Y <- 2 * X2 - 3 * X3 + rnorm(n)
  X4 <- 3 * (E == 4) + X3 - 3 * Y + rnorm(n)
  X5 <- 3 * (E == 5) - Y + rnorm(n)
  data.obs <- data.frame(cbind(Y,E,X1,X2,X3,X4,X5))
  
  
  
  data.obs <- as.matrix(data.obs[])
  
  result <- avChowICP2(data.obs,
                       target_name = "Y",
                       env_name = "E", 
                       m = 30)
  
  
  colors <- 1:ncol(result)
  matplot(log(result), type = "l")
  abline(h = log(1/alpha), col = "red")
  legend(
    "topright",
    legend = colnames(result),
    col = colors,
    lty = 1
  )
  
  #############################
  # ab hier 
  ############################
  notreject <- vector("list", nrow(result))
  icpset_over_time <- vector("list", nrow(result))
  icpset_over_time_strsplit <- vector("list", nrow(result))
  result_is_subset <- rep(NA, nrow(result))
  
  
  for(k in 1:nrow(result)){
    
    notreject[[k]] <- colnames(result)[result[k,]< 1/alpha]
    icpset_over_time[[k]] <- intersect_all_strings(notreject[[k]])
    icpset_over_time_strsplit[[k]] <- unlist(strsplit(icpset_over_time[[k]], "_"))

  }
  
  icp_set_over_time_list[[i]] <- unlist(icpset_over_time)  
  result_is_subset <- sapply(icpset_over_time_strsplit, is_subset, S = c("X2","X3", ""))
  result_is_subset_list[[i]] <- result_is_subset
  
  final_result_is_subset[i] <- all(result_is_subset)
  
}

final_result_is_subset
result_is_subset_list
icp_set_over_time_list
notreject

table(final_result_is_subset)

last_elements <- sapply(icp_set_over_time_list, function(x) x[length(x)])
last_elements

