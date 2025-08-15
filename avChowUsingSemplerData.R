setwd("C:/ETH/Thesis/Code/lindon related")
source("anytimeValidFTest_utils3.R")

# ------------------------------------
# test area

i <- 1
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
# idx<- 1:nrow(scm_data) %% 300 < 50
idx <- sample(1:nrow(scm_data), size = nrow(scm_data))
sum(idx)

scm_data <- scm_data[idx,]
scm_data <- scm_data[1:inner,]

##################################################
# X <- as.matrix(scm_data[,c("X1","X2","X3","X4","X5")])
# Y <- scm_data[,"Y"]
# E <- scm_data[,"E"]
# 
# Z <- factor(levels = unique(E))
# 
# for(e in unique(E)[-1]){
#   Z <- cbind(Z,E == e, sweep(X,1,E == e,"*"))
# }


######## why does the f statstic work without minimal batch size
# 
# for(i in 1:nrow(scm_data)){
# 
#     Ytest <- Y[1:i]
#     Ztest <- Z[1:i,]
# 
# 
# 
#     if(ncol(X) > 1){
#       Xtest <- X[1:i,]
#     }else{
#       Xtest <- X[1:i]
#     }
# 
#     if(i == 80){
#       print(dim(Xtest))
#       print(dim(Ztest))
# 
#       model1 <- lm(Ytest ~ Xtest + Ztest)
#       model2 <- lm(Ytest ~ Xtest)
# 
#       print(summary(model1))
#       print(summary(model2))
#     }
# 
#   }


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



length(icp_results)
jaccard_similarity(icp_results, parents)
#-------------------------------------



outer <- 300
inner <- 200
icp_results_list <- rep(NA,times = outer)
jaccard <- rep(NA,times = outer)
parentlist <- rep(NA,times = outer)


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
  idx <- sample(1:nrow(scm_data), size = nrow(scm_data))
  scm_data <- scm_data[idx,]
  scm_data <- scm_data[1:inner,]
  
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
  
  
  icp_results_list[i] <- intersect_all_strings(notreject)
  #print(icp_results)
  #print(parents)
  #print(jaccard_similarity(icp_results, parents))
  jaccard[i] <-   jaccard_similarity(icp_results, parents)

  # print(paste0("jaccard before next loop: ", jaccard))
}

notreject
icp_results
icp_results_list
parentlist
mean(jaccard)  
jaccard
table(jaccard)
hist(jaccard)
