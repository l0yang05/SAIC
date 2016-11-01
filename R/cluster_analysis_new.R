cluster_analysis_new <-
function(data, cluster_number, pvalue, number=50, method){
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
