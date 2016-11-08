k_aov <-
function(data,n, p=0.01, method){
  require(MASS)
  require(cclust)
  require(clusterSim)


  if(method == "kmeans") {
	  clust <- kmeans(data, n ,iter.max = 100,nstart = 100)
  } else if (method == "neuralgas") {
	  clust <-  cclust(data, centers=n, method=method)
  }
  ##print(paste("cluster number is set to be",n,sep=" "))
  dbs <- index.DB(data, as.integer(clust$cluster), centrotypes="centroids")[[1]]
  print(dbs)
  pValue<-sapply(1:ncol(data),function(i)summary(aov(gene~cl, data=data.frame(gene=data[,i],cl=as.character(clust$cluster))))[[1]][1,5])
  FDR <- p.adjust(pValue, method = "fdr", n = length(pValue))
  keep2 <- which(FDR<= p)
  data_table <- data[,keep2]
  print(paste("number of genes left is",ncol(data_table)))
  output <- list(dbScore=dbs,cluster=clust$cluster,table=data_table)
  output
}
