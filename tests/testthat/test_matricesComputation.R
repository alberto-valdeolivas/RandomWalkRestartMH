library(RandomWalkRestartMH)

############################################################
context("Matrices Computation")
############################################################

## Multiplex
m1 <- igraph::graph(c(1,2,1,3,2,3), directed = FALSE)
m2 <- igraph::graph(c(1,3,2,3,3,4,1,4), directed = FALSE)
multiObject_1 <- create.multiplex(list(m1=m1,m2=m2))
AdjMatrix <- compute.adjacency.matrix(multiObject_1)
AdjMatrixNorm <- normalize.multiplex.adjacency(AdjMatrix)

AdjMatrixExpected <- Matrix::Matrix(c(0,0.5,0.5,0,0.5,0,0,0,
                                0.5,0,0.5,0,0,0.5,0,0,
                                0.5,0.5,0,0,0,0,0.5,0,
                                0,0,0,0,0,0,0,0.5,
                                0.5,0,0,0,0,0,0.5,0.5,
                                0,0.5,0,0,0,0,0.5,0,
                                0,0,0.5,0,0.5,0.5,0,0.5,
                                0,0,0,0.5,0.5,0,0.5,0),
                                byrow = TRUE, nrow = 8, ncol = 8)

colnames(AdjMatrixExpected) <- c("1_1","2_1","3_1","4_1","1_2","2_2",
                                "3_2","4_2")
rownames(AdjMatrixExpected) <- c("1_1","2_1","3_1","4_1","1_2","2_2",
                                 "3_2","4_2")

AdjMatrixExpected <- as(AdjMatrixExpected, "dgCMatrix")


AdjMatrixNormExpected <- matrix(c(0,0.3333333,0.3333333,0,0.3333333,0,0,0,
                0.3333333,0,0.3333333,0,0,0.5,0,0,
				0.3333333,0.3333333,0,0,0,0,0.25,0,
				0,0,0,0,0,0,0,0.3333333,
				0.3333333,0,0,0,0,0,0.25,0.3333333,
				0,0.3333333,0,0,0,0,0.25,0,
				0,0,0.3333333,0,0.3333333,0.5,0,0.3333333,
				0,0,0,1,0.3333333,0,0.25,0),
				byrow = TRUE, nrow = 8, ncol = 8)

colnames(AdjMatrixNormExpected) <- c("1_1","2_1","3_1","4_1","1_2","2_2",
                                "3_2","4_2")
rownames(AdjMatrixNormExpected) <- c("1_1","2_1","3_1","4_1","1_2","2_2",
                                 "3_2","4_2")


AdjMatrixNormExpected <- as(AdjMatrixNormExpected, "dgCMatrix")


## Multiplex-Heterogeneous
h1 <- igraph::graph(c("A","C","B","E","E","D","E","C"), directed = FALSE)
multiObject_2 <- create.multiplex(list(h1=h1))

bipartite_relations <- data.frame(m=c(1,3),h=c("A","E"))

multiHetObject <- 
    create.multiplexHet(multiObject_1, multiObject_2, bipartite_relations)
MultiHetTranMatrix <- compute.transition.matrix(multiHetObject)

MultiHetTranMatrixExpected <- matrix(
        c(0,0.3333333,0.1666667,0,0.1666667,0,0,0,0.25,0,0,0,0,
        0.1666667,0,0.1666667,0,0,0.5,0,0,0,0,0,0,0,
        0.1666667,0.3333333,0,0,0,0,0.125,0,0,0,0,0,0.2500000,
        0,0,0,0,0,0,0,0.3333333,0,0,0,0,0,
        0.1666667,0,0,0,0,0,0.125,0.3333333,0.25,0,0,0,0,
        0,0.3333333,0,0,0,0,0.125,0,0,0,0,0,0,
        0,0,0.1666667,0,0.1666667,0.5,0,0.3333333,0,0,0,0,0.2500000,
        0,0,0,1,0.1666667,0,0.125,0,0,0,0,0,0,
        0.5000000,0,0,0,0.5000000,0,0,0,0,0,0.5,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0.1666667,
        0,0,0,0,0,0,0,0,0.50,0,0,0,0.1666667,
        0,0,0,0,0,0,0,0,0,0,0,0,0.1666667,
        0,0,0.5000000,0,0,0,0.500,0,0,1,0.5,1,0),
        byrow = TRUE, nrow = 13, ncol = 13)



colnames(MultiHetTranMatrixExpected) <- c("1_1","2_1","3_1","4_1","1_2","2_2",
                                 "3_2","4_2","A_1","B_1","C_1","D_1","E_1")
rownames(MultiHetTranMatrixExpected) <- c("1_1","2_1","3_1","4_1","1_2","2_2",
                                 "3_2","4_2","A_1","B_1","C_1","D_1","E_1")

MultiHetTranMatrixExpected <- as(MultiHetTranMatrixExpected, "dgCMatrix")


test_that("Check that Matrices computation is ok", {
    expect_equal(AdjMatrix, AdjMatrixExpected)
    expect_equal(AdjMatrixNorm, AdjMatrixNormExpected, tolerance = 0.00001)
    expect_equal(MultiHetTranMatrix, MultiHetTranMatrixExpected, tolerance = 0.00001)
})
