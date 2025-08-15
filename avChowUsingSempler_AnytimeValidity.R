is_subset <- function(X, S) {
  all(X %in% S | is.na(X))
}



setwd("C:/ETH/Thesis/Code/lindon related")
source("anytimeValidFTest_utils3.R")

# ------------------------------------
# test area

i <- 2
inner <- 300
#print(paste0("jaccard after next loop: ", jaccard))
#print(paste0("outer simulation nr: ",i))

filename <- paste0("C:/ETH/Thesis/Code/sempler/new_data_",i-1,".csv")
dagname <- paste0("C:/ETH/Thesis/Code/sempler/new_dag_",i-1,".txt")
scm_data <- read.csv(filename)
scm_data <- as.data.frame(scm_data)
names(scm_data)[1] <- "Y" 
names(scm_data)[length(scm_data)] <- "E"
scm_data[,"E"] <- as.factor(scm_data[,"E"])

dag <-read.csv(dagname, header = FALSE)
rownames(dag) <- names(scm_data[-length(scm_data)])
colnames(dag) <- names(scm_data[-length(scm_data)])
parents <- rownames(dag)[dag[,1] != 0]
parents <- paste(parents, collapse = "_")


target_name <- "Y"
env_name <- "E"
predictor_names <- colnames(scm_data)[! colnames(scm_data) %in% c(target_name,env_name)]


# X <- data.obs[,predictor_names]
# E <- as.factor(scm_data$E)

# Since the observations are blocked via environments, randomize the order of the observations
set.seed(14)
idx<- 1:nrow(scm_data) %% 300 < 50
#idx <- sample(1:nrow(scm_data), size = nrow(scm_data))

scm_data <- scm_data[idx,]
scm_data <- scm_data[1:inner,]

result <- avChowICP2(scm_data,
                     target_name = "Y",
                     env_name = "E", 
                     g = 1,
                     m = 51, 
                     alpha = alpha)



# plot results
alpha <- 0.05
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
notreject <- colnames(result)[result[nrow(result),]< 1/alpha]
icp_results <-intersect_all_strings(notreject)
notreject
icp_results
parents


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

Sstar <- unlist(strsplit(parents, split = "_"))
Sstar

result_is_subset <- sapply(icpset_over_time_strsplit, is_subset, S = c(Sstar, ""))

final_result_is_subset <- all(result_is_subset)
final_result_is_subset



length(icp_results)
jaccard_similarity(icp_results, parents)
#-------------------------------------



outer <- 300
inner <- 200
icp_results_list <- rep(NA,times = outer)
jaccard <- rep(NA,times = outer)
parentlist <- rep(NA,times = outer)


icp_results <- rep(NA,times = outer)
final_result_is_subset <- rep(NA,times = outer)
icp_set_over_time_list <- vector("list", outer)
result_is_subset_list <- vector("list", outer)


for(i in 1:outer){
  
  #print(paste0("jaccard after next loop: ", jaccard))
  #print(paste0("outer simulation nr: ",i))
  
  filename <- paste0("C:/ETH/Thesis/Code/sempler/new_data_",i-1,".csv")
  dagname <- paste0("C:/ETH/Thesis/Code/sempler/new_dag_",i-1,".txt")
  scm_data <- read.csv(filename)
  scm_data <- as.data.frame(scm_data)
  names(scm_data)[1] <- "Y" 
  names(scm_data)[length(scm_data)] <- "E"
  scm_data[,"E"] <- as.factor(scm_data[,"E"])
  
  dag <-read.csv(dagname, header = FALSE)
  rownames(dag) <- names(scm_data[-length(scm_data)])
  colnames(dag) <- names(scm_data[-length(scm_data)])
  parents <- rownames(dag)[dag[,1] != 0]
  parents <- paste(parents, collapse = "_")
  parentlist[i] <- parents 
  
  target_name <- "Y"
  env_name <- "E"
  predictor_names <- colnames(scm_data)[! colnames(scm_data) %in% c(target_name,env_name)]
  
  
  # X <- data.obs[,predictor_names]
  # E <- as.factor(scm_data$E)
  
  # Since the observations are blocked via environments, randomize the order of the observations
  set.seed(14)
  
  # use for blocked data by environments
  {
    idx<- 1:nrow(scm_data) %% 300 < 50
    scm_data <- scm_data[idx,]
  }
  
  # use for randomly sampled environments
  # {
  #   idx <- sample(1:nrow(scm_data), size = nrow(scm_data))
  #   scm_data <- scm_data[idx,]
  #   scm_data <- scm_data[1:inner,]
  # }

  
  # head(scm_data)
  result <- avChowICP2(scm_data,
                       target_name = "Y",
                       env_name = "E", 
                       g = 1,
                       m = 50, 
                       alpha = alpha)
  
  
  # plot results
  alpha <- 0.05
  colors <- 1:ncol(result)
  plot_name <- paste0("avICP_avChow simulation nr:", i, ".jpeg")
  jpeg(file=plot_name, width = 1440, height = 1080)
  matplot(log(result), 
          type = "l",        # line plot
          lty = 1,           # lines
          col = colors,
          xlab = "Time",    # custom X-axis label
          ylab = "e-value",
          main = paste0("avICP-avChow: E-values over time, nr: ", i))
  abline(h = log(1/alpha), col = "red")
  legend(
    "topright",
    legend = colnames(result),
    col = colors,
    lty = 1
  )
  dev.off()
  
  #############################
  # For anytime validity 
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
  
  Sstar <- unlist(strsplit(parents, split = "_"))
  
  result_is_subset <- sapply(icpset_over_time_strsplit, is_subset, S = c(Sstar, ""))
  result_is_subset_list[[i]] <- result_is_subset
  

  
  final_result_is_subset[i] <- all(result_is_subset)

}


#################################################
# analysis
final_result_is_subset #anytime validity for ICP
result_is_subset_list
icp_set_over_time_list
notreject

table(final_result_is_subset)

last_elements <- sapply(icp_set_over_time_list, function(x) x[length(x)])
last_elements
parentlist
table(last_elements == parentlist)


which(final_result_is_subset == FALSE)

icp_set_over_time_list[[27]]

parentlist[which(final_result_is_subset == FALSE)] 
# given that the anytime validity did not hold, the true set of parents was empty

which(parentlist == "") %in% which(final_result_is_subset == FALSE) 
# search for datasets, where the empty set were parents, avICP hold anytime validity
idx <- which(parentlist == "")[!which(parentlist == "") %in% which(final_result_is_subset == FALSE)]
icp_set_over_time_list[[56]]


