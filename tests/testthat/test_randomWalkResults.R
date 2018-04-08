library(RandomWalkRestartMH)

############################################################
context("Random Walk Computation")
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

RWRM_ExpectedResults <- data.frame(NodeNames = c("3","2","4"),
                                    Score = c(0.04934363, 
                                    0.01963822, 0.01421403),
                                    stringsAsFactors = FALSE)

## Multiplex-Heterogeneous
h1 <- igraph::graph(c("A","C","B","E","E","D","E","C"), directed = FALSE)
bipartite_relations <- data.frame(m=c(1,3),h=c("A","E"))

multiHetObject <- create.multiplexHet(multiObject,
    h1,bipartite_relations)
MultiHetTranMatrix <- compute.transition.matrix(multiHetObject)

SecondNet_Seeds <- c("E")

RWR_MultiHetResults <- Random.Walk.Restart.MultiplexHet(MultiHetTranMatrix,
    multiHetObject,Multiplex_Seeds,SecondNet_Seeds)


RWRMH_MultiExpectedResults <- data.frame(NodeNames = c("3","2","4"),
                                    Score = c(0.041135932,  
                                    0.005778811,0.003637140), 
                                    stringsAsFactors = FALSE)

RWRMH_SecondExpectedResults <- data.frame(SecondNet_node = c("A","C","B","D"),
                                    Score = c(0.06189589, 0.02817966,
                                    0.01889527, 0.01889527),
                                    stringsAsFactors = FALSE)


rownames(RWR_MultiHetResults$RWRMH_Results_MultiplexNodes) <- NULL
rownames(RWRMH_MultiExpectedResults) <- NULL

rownames(RWR_MultiHetResults$RWRMH_Results_SecondNetNodes) <- NULL
rownames(RWRMH_SecondExpectedResults) <- NULL

test_that("Random Walk Results are correct", {
    expect_equal(RWR_MultiResults$RWRM_Results, RWRM_ExpectedResults, 
                tolerance = 0.00001)
    expect_equal(RWR_MultiHetResults$RWRMH_Results_MultiplexNodes,
                RWRMH_MultiExpectedResults,
                tolerance = 0.00001)
    expect_equal(RWR_MultiHetResults$RWRMH_Results_SecondNetNodes, 
                RWRMH_SecondExpectedResults,
                tolerance = 0.00001)
})

