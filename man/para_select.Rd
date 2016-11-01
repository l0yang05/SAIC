\name{para_select}
\alias{para_select}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
para_select(data, cluster, pvalue, number, method = "neuralgas")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
%%     ~~Describe \code{data} here~~
}
  \item{cluster}{
%%     ~~Describe \code{cluster} here~~
}
  \item{pvalue}{
%%     ~~Describe \code{pvalue} here~~
}
  \item{number}{
%%     ~~Describe \code{number} here~~
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
function (data, cluster, pvalue, number, method = "neuralgas") 
{
    require(parallel)
    result_all <- mclapply(1:length(cluster), function(i) {
        print(paste("cluster number", cluster[i], "starts..."))
        result1 <- array(dim = length(pvalue))
        for (j in 1:length(pvalue)) {
            result <- cluster_analysis_new(data = data, cluster_number = cluster[i], 
                pvalue = pvalue[j], number, method = method)
            if (!is.null(result)) {
                result1[j] <- result$dbScore
            }
            else {
                result1[j] <- NA
            }
        }
        return(result1)
    }, mc.cores = length(cluster))
    result_all <- do.call("rbind", result_all)
    colnames(result_all) <- pvalue
    row.names(result_all) <- cluster
    return(result_all)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
