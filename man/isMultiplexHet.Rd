% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CheckObjects.R
\name{isMultiplexHet}
\alias{isMultiplexHet}
\title{Is this R object a Multiplex Heterogeneous object?}
\usage{
isMultiplexHet(x)
}
\arguments{
\item{x}{An R object}
}
\value{
A logical constant, \code{TRUE} if argument \code{x} is
a MultiplexHet object.
}
\description{
A Multiplex Heterogeneous object is an R object generated as the result
of calling the function \code{create.multiplexHet}
}
\details{
\code{isMultiplexHet(x)} checks whether an R object is MultiplexHet
}
\examples{
m1 <- igraph::graph(c(1,2,1,3,2,3), directed = FALSE)
m2 <- igraph::graph(c(1,3,2,3,3,4,1,4), directed = FALSE)
multiObject_1 <- create.multiplex(list(m1=m1,m2=m2))
h1 <- igraph::graph(c("A","C","B","E","E","D","E","C"), directed = FALSE)
multiObject_2 <- create.multiplex(list(h1=h1))
bipartite_relations <- data.frame(m=c(1,3),h=c("A","E"))
multiHetObject <- 
    create.multiplexHet(multiObject_1,multiObject_2,bipartite_relations)
isMultiplexHet(multiHetObject)
isMultiplexHet(h1)

}
\seealso{
\code{\link{create.multiplexHet}},
\code{\link{isMultiplex}}
}
\author{
Alberto Valdeolivas Urbelz \email{alvaldeolivas@gmail.com}
}
