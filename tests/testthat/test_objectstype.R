library(RandomWalkRestartMH)

############################################################
context("Objects Comprobation")
############################################################

## Multiplex
m1 <- igraph::graph(c(1,2,1,3,2,3), directed = FALSE)
m2 <- igraph::graph(c(1,3,2,3,3,4,1,4), directed = FALSE)
multiObject <- create.multiplex(m1,m2)
AdjMatrix <- compute.adjacency.matrix(multiObject)
AdjMatrixNorm <- normalize.multiplex.adjacency(AdjMatrix)
Multiplex_Seeds <- c(1)

RWR_MultiResults <- Random.Walk.Restart.Multiplex(AdjMatrixNorm, multiObject,
    Multiplex_Seeds)


## Multiplex-Heterogeneous
h1 <- igraph::graph(c("A","C","B","E","E","D","E","C"), directed = FALSE)
bipartite_relations <- data.frame(m=c(1,3),h=c("A","E"))

multiHetObject <- create.multiplexHet(multiObject,
    h1,bipartite_relations)
MultiHetTranMatrix <- compute.transition.matrix(multiHetObject)

SecondNet_Seeds <- c("E")

RWR_MultiHetResults <- Random.Walk.Restart.MultiplexHet(MultiHetTranMatrix,
    multiHetObject,Multiplex_Seeds,SecondNet_Seeds)

test_that("Objects are of the expected type", {
    expect_equal(isMultiplex(multiObject), TRUE)
    expect_equal(isMultiplexHet(multiHetObject), TRUE)
    expect_equal(isRWRM_Results(RWR_MultiResults), TRUE)
    expect_equal(isRWRMH_Results(RWR_MultiHetResults), TRUE)
})

