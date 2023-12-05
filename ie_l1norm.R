rm(list = ls())

getwd()
options(stringsAsFactors = F)

#library("xlsx")

#miRNA
miR_sim_matrix <- read.csv('./network/ie_toni_path.txt',header = F,sep = '\t')
#不断的替换这个矩阵
miR_sim_matrix[is.na(miR_sim_matrix)] <- 0
all_mir <- read.csv('./data/all_mirna.txt',sep = '\t',header = F)
row.names(miR_sim_matrix) <- all_mir$V2
colnames(miR_sim_matrix) <- all_mir$V2

#disease
disease_sim_matrix <- read.csv('./data/disease_semantic_similarity.txt',sep = '\t',header = F)
all_dis <- read.csv('./data/all_disease.txt',sep = '\t',header = F)
row.names(disease_sim_matrix) <- all_dis$V2
colnames(disease_sim_matrix) <- all_dis$V2

#association
miR_disease_matrix <- read.csv('./network/hmdd_causal_for_train_ie.txt',sep = '\t',header = T,row.names = 1)

miR_sim_matrix <- as.matrix(miR_sim_matrix)
disease_sim_matrix <- as.matrix(disease_sim_matrix)
miR_disease_matrix <- as.matrix(miR_disease_matrix )

#Cite from: Liang C, Yu S, Wong KC, Luo J: A novel semi-supervised model for miRNA-disease association prediction based on [Formula: see text]-norm graph. J Transl Med 2018, 16(1):357.
con_condition <- 1e-6
factor <- 0.5

LP_norm1 <- function(W, Y, con_condition, n){
  D <- diag(apply(W, 1, sum))
  
  U <- diag(1, nrow(W), ncol(W))
  if(!missing(n)){
    if(length(n) == 1) { 
      diag(U)[Y[, n, drop=T]==1] <- 1000000
    } else {
      diag(U)[apply(Y[, n]==1, 1, sum) >= 1] <- 1000000
    }		
  }
  
  L <- D - W
  Q <- solve(L + U, U %*% Y)
  
  distQ <- sqrt(abs(L2_distance(Q, Q)) + .Machine$double.eps)
  
  last_obj <- Inf
  cur_obj <- 0.5 * sum(distQ * W) + sum(diag(t(Q - Y) %*% U %*% (Q - Y)))
  #	print(cur_obj)
  
  i <- 1
  #	while(i <= 20){
  while(last_obj - cur_obj > con_condition){
    W1 <- W / (2 * distQ)
    D1 <- diag(apply(W1, 1, sum))
    L1 <- D1 - W1
    Q <- solve(L1 + U, U %*% Y)
    
    distQ <- sqrt(abs(L2_distance(Q, Q)) + .Machine$double.eps)
    last_obj <- cur_obj
    cur_obj <- 0.5 * sum(distQ * W) + sum(diag(t(Q - Y) %*% U %*% (Q - Y)))
    
    #		print(cur_obj)
    #		i <- i + 1
  }
  
  list(Q=Q, W1=W1)
}

L2_distance <- function(A, B){
  rowA <- apply(A*A, 1, sum)
  matrixA <- matrix(rep(rowA, each=length(rowA)), nrow=length(rowA), byrow=T)
  rowB <- apply(B*B, 1, sum)
  matrixB <- matrix(rep(rowB, each=length(rowB)), nrow=length(rowB), byrow=F)
  C <- 2 * A %*% t(B)
  dis <- Re(matrixA + matrixB - C)
}


miR_info <- LP_norm1(miR_sim_matrix, miR_disease_matrix, con_condition)
disease_info <- LP_norm1(disease_sim_matrix, t(miR_disease_matrix), con_condition)
finalQ <- miR_info$Q * factor + t(disease_info$Q) * (1 - factor)

write.table (finalQ, file ="./result/result_matrix_ie_toni_path.txt", sep ="\t", row.names =TRUE, col.names =TRUE, quote =F)
