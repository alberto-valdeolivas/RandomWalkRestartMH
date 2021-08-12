## R Code for the Random Walk with Restart Package (RandomWalkRestartMH).

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Functions to check the objects used within the package.
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

## Roxy Documentation comments
#' Is this R object a Multiplex object?
#'
#' A Multiplex object is an R object generated as the result of calling
#' the function \code{create.multiplex}
#'+
#' \code{isMultiplex(x)} checks whether an R object is Multiplex.
#'
#' @param x An R object
#'
#' @return A logical constant, \code{TRUE} if argument \code{x} is
#' a Mutiplex object.
#'
#' @seealso \code{\link{create.multiplex}}, \code{\link{isMultiplexHet}}
#'
#' @author Alberto Valdeolivas Urbelz \email{alvaldeolivas@@gmail.com}
#'
#' @examples
#' m1 <- igraph::graph(c(1,2,1,3,2,3), directed = FALSE)
#' m2 <- igraph::graph(c(1,3,2,3,3,4,1,4), directed = FALSE)
#' multiObject <- create.multiplex(list(m1=m1,m2=m2))
#' isMultiplex(multiObject)
#' isMultiplex(m1)
#'
#'@importFrom methods is
#'@export

isMultiplex <- function (x)
{
  is(x,"Multiplex")
}


## Roxy Documentation comments
#' Is this R object a Multiplex Heterogeneous object?
#'
#' A Multiplex Heterogeneous object is an R object generated as the result
#' of calling the function \code{create.multiplexHet}
#'
#' \code{isMultiplexHet(x)} checks whether an R object is MultiplexHet
#'
#' @param x An R object
#'
#' @return A logical constant, \code{TRUE} if argument \code{x} is
#' a MultiplexHet object.
#'
#' @seealso \code{\link{create.multiplexHet}},
#' \code{\link{isMultiplex}}
#'
#' @author Alberto Valdeolivas Urbelz \email{alvaldeolivas@@gmail.com}
#'
#' @examples
#' m1 <- igraph::graph(c(1,2,1,3,2,3), directed = FALSE)
#' m2 <- igraph::graph(c(1,3,2,3,3,4,1,4), directed = FALSE)
#' multiObject_1 <- create.multiplex(list(m1=m1,m2=m2))
#' h1 <- igraph::graph(c("A","C","B","E","E","D","E","C"), directed = FALSE)
#' multiObject_2 <- create.multiplex(list(h1=h1))
#' bipartite_relations <- data.frame(m=c(1,3),h=c("A","E"))
#' multiHetObject <- 
#'     create.multiplexHet(multiObject_1,multiObject_2,bipartite_relations)
#' isMultiplexHet(multiHetObject)
#' isMultiplexHet(h1)
#'
#'@importFrom methods is
#'@export
isMultiplexHet <- function (x)
{
  is(x,"MultiplexHet")
}

## Roxy Documentation comments
#' Is this R object a RWR on Multiplex object (Results of the RWR-M)?
#'
#' A RWR on Multiplex object is an R object generated as the result
#' of calling the function \code{Random.Walk.Restart.Multiplex}
#' (Results of the RWR-M)
#'
#' \code{isRWRM_Results(x)} checks whether an R object is RWRM_Results
#'
#' @param x An R object
#'
#' @return A logical constant, \code{TRUE} if argument \code{x} is
#' a RWRM_Results object.
#'
#' @seealso \code{\link{Random.Walk.Restart.Multiplex}},
#' \code{\link{isRWRMH_Results}}
#'
#' @author Alberto Valdeolivas Urbelz \email{alvaldeolivas@@gmail.com}
#'
#' @examples
#' m1 <- igraph::graph(c(1,2,1,3,2,3), directed = FALSE)
#' m2 <- igraph::graph(c(1,3,2,3,3,4,1,4), directed = FALSE)
#' multiObject <- create.multiplex(list(m1=m1,m2=m2))
#' AdjMatrix <- compute.adjacency.matrix(multiObject)
#' AdjMatrixNorm <- normalize.multiplex.adjacency(AdjMatrix)
#' Seed <- c(1)
#' RWR_MultiResults <- 
#'     Random.Walk.Restart.Multiplex(AdjMatrixNorm, multiObject,Seed)
#' isRWRM_Results(RWR_MultiResults)
#' isRWRM_Results(m1)
#'@importFrom methods is
#'@export
isRWRM_Results <- function (x)
{
  is(x,"RWRM_Results")
}


## Roxy Documentation comments
#' Is this R object a RWR on Multiplex-Heterogeneous object (Results of the
#' RWR-MH)?
#'
#' A RWR on Multiplex Heterogeneous object is an R object generated as the
#' result of calling the function \code{Random.Walk.Restart.MultiplexHet}
#' (Results of the RWR-MH)
#'
#' \code{isRWRMH_Results(x)} checks whether an R object is RWRMH_Results
#'
#' @param x An R object
#'
#' @return A logical constant, \code{TRUE} if argument \code{x} is
#' a RWRMH_Results object.
#'
#' @seealso \code{\link{Random.Walk.Restart.MultiplexHet}},
#' \code{\link{isRWRM_Results}}
#'
#' @author Alberto Valdeolivas Urbelz \email{alvaldeolivas@@gmail.com}
#'
#' @examples
#' m1 <- igraph::graph(c(1,2,1,3,2,3), directed = FALSE)
#' m2 <- igraph::graph(c(1,3,2,3,3,4,1,4), directed = FALSE)
#' multiObject_1 <- create.multiplex(list(m1=m1,m2=m2))
#' h1 <- igraph::graph(c("A","C","B","E","E","D","E","C"), directed = FALSE)
#' multiObject_2 <- create.multiplex(list(h1=h1))
#' bipartite_relations <- data.frame(m=c(1,3),h=c("A","E"))
#' multiHetObject <-
#'     create.multiplexHet(multiObject_1,multiObject_2,bipartite_relations)
#' MultiHetTranMatrix <- compute.transition.matrix(multiHetObject)
#' Multiplex1_Seeds <- c(1)
#' Multiplex2_Seeds <- c("E")
#' RWR_MultiHetResults <- 
#'     Random.Walk.Restart.MultiplexHet(MultiHetTranMatrix,multiHetObject,
#'     Multiplex1_Seeds ,Multiplex2_Seeds)
#' isRWRMH_Results(RWR_MultiHetResults)
#' isRWRMH_Results(m1)
#'
#'@importFrom methods is
#'@export
isRWRMH_Results <- function (x)
{
  is(x,"RWRMH_Results")
}