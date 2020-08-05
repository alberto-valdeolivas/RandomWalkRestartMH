library(RandomWalkRestartMH)

############################################################
context("Random Walk Computation")
############################################################

## Multiplex
m1 <- igraph::graph(c(1,2,1,3,2,3), directed = FALSE)
m2 <- igraph::graph(c(1,3,2,3,3,4,1,4), directed = FALSE)
multiObject_1 <- create.multiplex(list(m1=m1,m2=m2))
AdjMatrix <- compute.adjacency.matrix(multiObject_1)
AdjMatrixNorm <- normalize.multiplex.adjacency(AdjMatrix)
Multiplex1_Seeds <- c(1)

RWR_MultiResults <- 
    Random.Walk.Restart.Multiplex(AdjMatrixNorm, multiObject_1,Multiplex1_Seeds)

RWRM_ExpectedResults <- data.frame(NodeNames = c("3","2","4"),
                                    Score = c(0.04934363, 
                                    0.01963822, 0.01421403),
                                    stringsAsFactors = FALSE)

## Multiplex-Heterogeneous
h1 <- igraph::graph(c("A","C","B","E","E","D","E","C"), directed = FALSE)
multiObject_2 <- create.multiplex(list(h1=h1))
bipartite_relations <- data.frame(m=c(1,3),h=c("A","E"))

multiHetObject <- 
    create.multiplexHet(multiObject_1, multiObject_2, bipartite_relations)
MultiHetTranMatrix <- compute.transition.matrix(multiHetObject)

Multiplex2_Seeds <- c("E")

RWR_MultiHetResults <- Random.Walk.Restart.MultiplexHet(MultiHetTranMatrix,
    multiHetObject,Multiplex1_Seeds,Multiplex2_Seeds)


RWRMH_Multi1ExpectedResults <- data.frame(NodeNames = c("3","2","4"),
                                    Score = c(0.041135932,  
                                    0.005778811,0.003637140), 
                                    stringsAsFactors = FALSE)

RWRMH_Multi2ExpectedResults <- data.frame(NodeNames = c("A","C","B","D"),
                                    Score = c(0.06189589, 0.02817966,
                                    0.01889527, 0.01889527),
                                    stringsAsFactors = FALSE)


rownames(RWR_MultiHetResults$RWRMH_Multiplex1) <- NULL
rownames(RWRMH_Multi1ExpectedResults) <- NULL

rownames(RWR_MultiHetResults$RWRMH_Multiplex2) <- NULL
rownames(RWRMH_Multi2ExpectedResults) <- NULL

test_that("Random Walk Results are correct", {
    expect_equal(RWR_MultiResults$RWRM_Results, RWRM_ExpectedResults, 
                tolerance = 0.00001)
    expect_equal(RWR_MultiHetResults$RWRMH_Multiplex1,
                RWRMH_Multi1ExpectedResults,
                tolerance = 0.00001)
    expect_equal(RWR_MultiHetResults$RWRMH_Multiplex2, 
                RWRMH_Multi2ExpectedResults,
                tolerance = 0.00001)
})

