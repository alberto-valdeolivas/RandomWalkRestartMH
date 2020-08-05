library(RandomWalkRestartMH)

############################################################
context("Objects Comprobation")
############################################################

## Multiplex
m1 <- igraph::graph(c(1,2,1,3,2,3), directed = FALSE)
m2 <- igraph::graph(c(1,3,2,3,3,4,1,4), directed = FALSE)
multiObject_1 <- create.multiplex(list(m1=m1,m2=m2))
AdjMatrix <- compute.adjacency.matrix(multiObject_1)
AdjMatrixNorm <- normalize.multiplex.adjacency(AdjMatrix)
Multiplex1_Seeds <- c(1)

RWR_MultiResults <- Random.Walk.Restart.Multiplex(AdjMatrixNorm, multiObject_1,
    Multiplex1_Seeds)

## Multiplex-Heterogeneous
h1 <- igraph::graph(c("A","C","B","E","E","D","E","C"), directed = FALSE)
multiObject_2 <- create.multiplex(list(h1=h1))
    
bipartite_relations <- data.frame(m=c(1,3),h=c("A","E"))

multiHetObject <- 
    create.multiplexHet(multiObject_1,multiObject_2,bipartite_relations)
MultiHetTranMatrix <- compute.transition.matrix(multiHetObject)

Multiplex2_Seeds <- c("E")

RWR_MultiHetResults <- Random.Walk.Restart.MultiplexHet(MultiHetTranMatrix,
    multiHetObject,Multiplex1_Seeds,Multiplex2_Seeds)

test_that("Objects are of the expected type", {
    expect_equal(isMultiplex(multiObject_1), TRUE)
    expect_equal(isMultiplexHet(multiHetObject), TRUE)
    expect_equal(isRWRM_Results(RWR_MultiResults), TRUE)
    expect_equal(isRWRMH_Results(RWR_MultiHetResults), TRUE)
})

