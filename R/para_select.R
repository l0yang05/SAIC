para_select <-
function(data,cluster,pvalue, number, method="neuralgas"){
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
	  #write.table(result$cluster,paste("cluster_ij",i,"_",j,".txt",sep=""),sep="\t",col.names=NA)
    }
	return(result1)
  }, mc.cores=length(cluster))
  #}, mc.cores=1)
  result_all <- do.call("rbind",result_all)
  colnames(result_all) <- pvalue; row.names(result_all) <- cluster
  return(result_all)
  
}
