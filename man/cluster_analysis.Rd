\name{cluster_analysis}
\alias{cluster_analysis}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
cluster_analysis(data, cluster_number, pvalue, number = 50, method)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
%%     ~~Describe \code{data} here~~
}
  \item{cluster_number}{
%%     ~~Describe \code{cluster_number} here~~
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
function (data, cluster_number, pvalue, number = 50, method) 
{
    i = 1
    result <- k_aov(data = data, n = cluster_number, p = pvalue, 
        method = method)
    if (ncol(as.matrix(result$table)) < number) {
        result$dbScore <- "NA"
        print(paste("Gene set number is less than ", number, 
            " dbScore is set to NA"))
        return(NULL)
    }
    ns <- ncol(as.matrix(result$table))
    current_table <- result$table
    while (ns > number) {
        current_table <- result$table
        result <- k_aov(data = result$table, n = cluster_number, 
            p = pvalue, method = method)
        if (ncol(as.matrix(result$table)) < ns) {
            i = i + 1
            ns <- ncol(as.matrix(result$table))
        }
        else {
            print(paste("done cluster", cluster_number, "and", 
                pvalue))
            break
        }
    }
    result$table <- current_table
    print(paste("final number of genes kept is", ncol(result$table)))
    return(result)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
