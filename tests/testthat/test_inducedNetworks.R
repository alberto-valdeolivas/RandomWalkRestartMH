library(RandomWalkRestartMH)

context("Top-result network construction")

test_that("multiplex top-result networks retain induced edges and components", {
  layer <- igraph::graph_from_edgelist(
    matrix(c("A", "B", "C", "D"), ncol = 2, byrow = TRUE),
    directed = FALSE)
  multiplex <- create.multiplex(list(interaction = layer))

  results <- list(
    RWRM_Results = data.frame(
      NodeNames = c("B", "C"), Score = c(2, 1),
      stringsAsFactors = FALSE),
    Seed_Nodes = "A")
  class(results) <- "RWRM_Results"

  induced <- create.multiplexNetwork.topResults(results, multiplex, k = 2)

  expect_s3_class(induced, "igraph")
  expect_setequal(igraph::V(induced)$name, c("A", "B", "C"))
  expect_equal(igraph::ecount(induced), 1)
  expect_equal(igraph::E(induced)$type, "interaction")
  expect_equal(igraph::V(induced)$comp[match(c("A", "B"),
                                             igraph::V(induced)$name)],
               c(1, 1))
  expect_equal(igraph::V(induced)$comp[match("C", igraph::V(induced)$name)], 2)
})

test_that("heterogeneous top-result networks retain all edge types", {
  layer1 <- igraph::graph_from_edgelist(
    matrix(c("A", "B", "C", "D"), ncol = 2, byrow = TRUE),
    directed = FALSE)
  layer2 <- igraph::graph_from_edgelist(
    matrix(c("X", "Y", "Z", "W"), ncol = 2, byrow = TRUE),
    directed = FALSE)
  multiplex1 <- create.multiplex(list(first = layer1))
  multiplex2 <- create.multiplex(list(second = layer2))
  relations <- data.frame(left = "B", right = "X")
  heterogeneous <- create.multiplexHet(
    multiplex1, multiplex2, relations)

  results <- list(
    RWRMH_GlobalResults = data.frame(),
    RWRMH_Multiplex1 = data.frame(NodeNames = "B", Score = 2),
    RWRMH_Multiplex2 = data.frame(NodeNames = "Y", Score = 1),
    Seed_Nodes = c("A", "X"))
  class(results) <- "RWRMH_Results"

  induced <- create.multiplexHetNetwork.topResults(
    results, heterogeneous, relations, bipartite_name = "between", k = 1)

  expect_s3_class(induced, "igraph")
  expect_setequal(igraph::V(induced)$name, c("A", "B", "X", "Y"))
  expect_equal(igraph::ecount(induced), 3)
  expect_setequal(igraph::E(induced)$type,
                  c("first", "second", "between"))
  expect_true(all(igraph::V(induced)$comp == 1))
})
