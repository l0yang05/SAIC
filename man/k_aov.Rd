\name{k_aov}
\alias{k_aov}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
k_aov(data, n, p = 0.01, method)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
%%     ~~Describe \code{data} here~~
}
  \item{n}{
%%     ~~Describe \code{n} here~~
}
  \item{p}{
%%     ~~Describe \code{p} here~~
}
  \item{method}{
%%     ~~Describe \code{method} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (data, n, p = 0.01, method) 
{
    require(clusterCrit)
    require(MASS)
    require(edgeR)
    require(cclust)
    require(cluster)
    require(clusterSim)
    if (method == "kmeans") {
        clust <- kmeans(data, n, iter.max = 100, nstart = 100)
    }
    else if (method == "neuralgas") {
        clust <- cclust(data, centers = n, method = method)
    }
    dbs <- index.DB(data, as.integer(clust$cluster), centrotypes = "centroids")[[1]]
    print(dbs)
    pValue <- sapply(1:ncol(data), function(i) summary(aov(gene ~ 
        cl, data = data.frame(gene = data[, i], cl = as.character(clust$cluster))))[[1]][1, 
        5])
    FDR <- p.adjust(pValue, method = "fdr", n = length(pValue))
    keep2 <- which(FDR <= p)
    data_table <- data[, keep2]
    print(paste("number of genes left is", ncol(data_table)))
    output <- list(dbScore = dbs, cluster = clust$cluster, table = data_table)
    output
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
