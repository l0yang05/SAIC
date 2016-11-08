## data: in put logg2fpkm/rpkm data matrix  :  row--sample names; columns-- gene names
## n: kmeans cluster number
## p: p value cut off for one-way anova
## nstart of kmeans, use 100 random initaiations and return the optimum

k_aov <- function(data,n, p=0.01, method){
  require(MASS)
  require(cclust)
  require(clusterSim)
  #print(method)

  if(method == "kmeans") {
	  clust <- kmeans(data, n ,iter.max = 100,nstart = 100)
  } else if (method == "neuralgas") {
	  clust <-  cclust(data, centers=n, method=method)
  }
  ##print(paste("cluster number is set to be",n,sep=" "))
  dbs <- index.DB(data, as.integer(clust$cluster), centrotypes="centroids")[[1]]
  print(dbs)
  pValue<-sapply(1:ncol(data),function(i)summary(aov(gene~cl, data=data.frame(gene=data[,i],cl=as.character(clust$cluster))))[[1]][1,5])
  keep2 <- which(pValue<= p)
  data_table <- data[,keep2]
  print(paste("number of genes left is",ncol(data_table)))
  output <- list(dbScore=dbs,cluster=clust$cluster,table=data_table)
  output
}

cluster_analysis <- function(data, cluster_number, pvalue, number=50, method){
  i =1
  result <- k_aov(data = data,n=cluster_number,p=pvalue, method=method)
  #print(paste("first gene set number is ",ncol(result$table)))
  
  ###--- check if deg number is less than 50 at this step, if so, change db_score to NA  ----
  if(ncol(as.matrix(result$table))<number){
    result$dbScore <- "NA"
	print(paste("Gene set number is less than ",number," dbScore is set to NA"))
	return(NULL)
  #return(result)
  }
  ns <- ncol(as.matrix(result$table))
  current_table <- result$table
  while(ns>number){
  #  browser()
	current_table <- result$table
    result <- k_aov(data = result$table,n=cluster_number,p=pvalue, method=method)
    if(ncol(as.matrix(result$table)) < ns){
      i=i+1
      ns <- ncol(as.matrix(result$table))
    }
    else{
		print(paste("done cluster", cluster_number, "and", pvalue))
      break
    }
  }
  #browser()
  result$table <- current_table
  print(paste("final number of genes kept is",ncol(result$table)))
  return(result)
}

para_select <- function(data,cluster,pvalue, number, method="neuralgas"){
  require(parallel)
  result_all <- mclapply(1:length(cluster), function(i) {
    print(paste("cluster number",cluster[i],"starts..."))
    result1 <- array(dim=length(pvalue))
	for (j in 1:length(pvalue)){
      result <- cluster_analysis_new(data = data,cluster_number = cluster[i],pvalue = pvalue[j], number, method=method)
      if(!is.null(result)) {
		result1[j] <- result$dbScore
      } else {
		result1[j] <- NA
    #result1[j] <- paste(NA,result$dbScore,sep=".")
	  }
    }
	return(result1)
  }, mc.cores=length(cluster))
  #}, mc.cores=1)
  result_all <- do.call("rbind",result_all)
  colnames(result_all) <- pvalue; row.names(result_all) <- cluster
  return(result_all)
  
}
