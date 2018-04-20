## R Code for the Random Walk with Restart Package (RandomWalkRestartMH).

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Functions to compute the matrices and perform the Random Walks.
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

## Roxy Documentation comments
#' Computes the adjacency matrix of a multiplex network
#'
#' \code{compute.adjacency.matrix} is a function to compute the adjacency
#' matrix of a multiplex network provided as a \code{Multiplex} object.
#'
#' @usage compute.adjacency.matrix(x,delta = 0.5)
#'
#' @details The parameter \code{delta} sets the probability to change between
#' layers at the next step. If delta = 0, the particle will always remain
#' in the same layer after a non-restart iteration. On the other hand, if
#' delta = 1, the particle will always change between layers, therefore
#' not following the specific edges of each layer.
#'
#' @param x A \code{Multiplex} object describing a multiplex network generated
#' by the function \code{create.multiplex}.
#' @param delta A numeric value between 0 and 1. It sets the probability
#' of performing inter-layer versus intra-layer transitions. It is set by
#' default to 0.5. See more details below.
#'
#' @return A square sparse adjacency matrix created with the \code{Matrix}
#' package.
#'
#' @seealso \code{\link{create.multiplex},\link{normalize.multiplex.adjacency},
#' \link{compute.transition.matrix}}
#'
#' @author Alberto Valdeolivas Urbelz \email{alvaldeolivas@@gmail.com}
#'
#' @examples
#' m1 <- igraph::graph(c(1,2,1,3,2,3), directed = FALSE)
#' m2 <- igraph::graph(c(1,3,2,3,3,4,1,4), directed = FALSE)
#' multiObject <- create.multiplex(m1,m2)
#' compute.adjacency.matrix(multiObject)
#'
#'@import igraph
#'@import Matrix
#'@importFrom methods as
#'@export
compute.adjacency.matrix <- function(x,delta = 0.5)
{
    if (!isMultiplex(x) & !isMultiplexHet(x)) {
        stop("Not a Multiplex or Multiplex Heterogeneous object")
    }
    if (delta > 1 || delta <= 0) {
        stop("Delta should be between 0 and 1")
    }
    
    N <- x$Number_of_Nodes_Multiplex
    L <- x$Number_of_Layers
    Layers_Names <- names(x)[seq(L)]
    
    ## IDEM_MATRIX.
    Idem_Matrix <- Matrix::Diagonal(N, x = 1)
    
    counter <- 0 
    Layers_List <- lapply(x[Layers_Names],function(x){
        
        counter <<- counter + 1;    
        Adjacency_Layer <-  as_adjacency_matrix(x,sparse = TRUE)
        Adjacency_Layer <- Adjacency_Layer[order(rownames(Adjacency_Layer)),
            order(colnames(Adjacency_Layer))]
        colnames(Adjacency_Layer) <- 
            paste0(colnames(Adjacency_Layer),"_",counter)
        rownames(Adjacency_Layer) <- 
            paste0(rownames(Adjacency_Layer),"_",counter)
        Adjacency_Layer
    })
    
    MyColNames <- unlist(lapply(Layers_List, function (x) unlist(colnames(x))))
    MyRowNames <- unlist(lapply(Layers_List, function (x) unlist(rownames(x))))
    names(MyColNames) <- c()
    names(MyRowNames) <- c()
    SupraAdjacencyMatrix <- (1-delta)*(bdiag(unlist(Layers_List)))
    colnames(SupraAdjacencyMatrix) <-MyColNames
    rownames(SupraAdjacencyMatrix) <-MyRowNames
    
    offdiag <- (delta/(L-1))*Idem_Matrix
    
    i <- seq_len(L)
    Position_ini_row <- 1 + (i-1)*N
    Position_end_row <- N + (i-1)*N
    j <- seq_len(L)
    Position_ini_col <- 1 + (j-1)*N
    Position_end_col <- N + (j-1)*N
    
    for (i in seq_len(L)){
        for (j in seq_len(L)){
            if (j != i){
                SupraAdjacencyMatrix[(Position_ini_row[i]:Position_end_row[i]),
                    (Position_ini_col[j]:Position_end_col[j])] <- offdiag
            }    
        }
    }
    
    SupraAdjacencyMatrix <- as(SupraAdjacencyMatrix, "dgCMatrix")
    return(SupraAdjacencyMatrix)
}

## Roxy Documentation comments
#' Computes column normalization of an adjacency matrix
#'
#' \code{normalize.multiplex.adjacency} is a function to compute the column
#' normalization of a sparse matrix of the package \code{Matrix}.
#'
#' @usage normalize.multiplex.adjacency(x)
#'
#' @param x A \code{Matrix} object describing an adjacency matrix of a network.
#'
#' @return A square sparse column normalized matrix created with the
#' \code{Matrix} package.
#'
#' @seealso \code{\link{compute.adjacency.matrix},
#' \link{Random.Walk.Restart.Multiplex}}
#'
#' @author Alberto Valdeolivas Urbelz \email{alvaldeolivas@@gmail.com}
#'
#' @examples
#' m1 <- igraph::graph(c(1,2,1,3,2,3), directed = FALSE)
#' m2 <- igraph::graph(c(1,3,2,3,3,4,1,4), directed = FALSE)
#' multiObject <- create.multiplex(m1,m2)
#' AdjMatrix <- compute.adjacency.matrix(multiObject)
#' normalize.multiplex.adjacency(AdjMatrix)
#'
#'@import Matrix
#'@export

normalize.multiplex.adjacency <- function(x)
{
    if (!is(x,"dgCMatrix")){
        stop("Not a dgCMatrix object of Matrix package")
    }
    
    Adj_Matrix_Norm <- t(t(x)/(Matrix::colSums(x, na.rm = FALSE, dims = 1,
        sparseResult = FALSE)))
    
    return(Adj_Matrix_Norm)
}

## Roxy Documentation comments
#' Performs Random Walk with Restart on a Multiplex Network
#'
#' \code{Random.Walk.Restart.Multiplex} is a function to perform a Random Walk
#' with Restart on a Multiplex network (on a \code{Multiplex} object). See
#' more details about the algorithm below.
#'
#' @usage Random.Walk.Restart.Multiplex(...)
#'
#' @details Random Walk with Restart simulates an imaginary particle that
#' starts on a seed(s) node(s) and follows randomly the edges of a network. At
#' each step, there is a restart probability, r, meaning that the particle comes
#' back to the seed(s). The extension to multiplex networks allows the particle
#' to explore different monoplex networks (layers). At each step, the particle
#' can also jump to the same node in a different layer.
#'
#'
#' \itemize{
#' \item \code{Seeds}: A vector containing the name of the different seed
#' node(s). It's mandatory to provide at least one seed. The seed(s) node(s)
#' should belong to any of the layers. The length of this vector should be
#' smaller than the total number of nodes in the multiplex network.
#' \item \code{r}: A numeric value representing the restart probability on the
#' seeds for the random walker. It must be between 0 and 1. It is set by default
#' to 0.7, which is the most common value in this kind of approaches. It means
#' that, at each step, the walker has a 70\% of probability of coming back to
#' one of the seeds.
#' \item \code{tau}: A numeric vector containing the probability of restarting
#' in the nodes of the different layers of the multiplex. In the example below,
#' we define the node 1 as the seed node. However, we can find this node in both
#' layers. Therefore, the walker can restart in any of these seed nodes. It is
#' a way to give different relevance (weight) to the different layers.
#' }
#'
#'
#' @param x An object of the \code{Matrix} package describing a column
#' normalized adjacency matrix of a multiplex network.
#' @param MultiplexObject A \code{Multiplex} object generated by the function
#' \code{create.multiplex} representing a multiplex network.
#' @param Seeds A vector containing the names of the seeds for the Random
#' Walk algorithm. See more details below.
#' @param r A numeric value between 0 and 1. It sets the probability of
#' restarting to a seed node after each step. See more details below.
#' @param tau A vector containing the probability of restart on the seeds
#' of the different layers (layers weights). It must have the same length than
#' the number of layers of the multpiplex network. The sum of its components
#' divided by the number of layers must be 1. See more details below.
#' @param ... Further arguments passed to \code{Random.Walk.Restart.Multiplex}
#'
#' @return A \code{RWRM_Results} object. It contains a sorted ranking of all
#' the nodes of the multiplex network, except the seeds, along with their score.
#' In addition, it contains in a different field the nodes used as seeds.
#'
#' @seealso \code{\link{create.multiplex},
#' \link{compute.adjacency.matrix}, \link{normalize.multiplex.adjacency},
#' \link{isRWRM_Results}, \link{Random.Walk.Restart.MultiplexHet}}
#'
#' @author Alberto Valdeolivas Urbelz \email{alvaldeolivas@@gmail.com}
#'
#' @examples
#' m1 <- igraph::graph(c(1,2,1,3,2,3), directed = FALSE)
#' m2 <- igraph::graph(c(1,3,2,3,3,4,1,4), directed = FALSE)
#' multiObject <- create.multiplex(m1,m2)
#' AdjMatrix <- compute.adjacency.matrix(multiObject)
#' AdjMatrixNorm <- normalize.multiplex.adjacency(AdjMatrix)
#' SeedNodes <- c(1)
#' Random.Walk.Restart.Multiplex(AdjMatrixNorm,multiObject,SeedNodes)
#'
#'@import igraph
#'@import Matrix
#'@rdname Random.Walk.Restart.Multiplex
#'@export
Random.Walk.Restart.Multiplex <- function(...) {
    UseMethod("Random.Walk.Restart.Multiplex")
}

#'@rdname Random.Walk.Restart.Multiplex
#'@export
Random.Walk.Restart.Multiplex.default <- 
    function(x, MultiplexObject, Seeds, r=0.7,tau,...){
        
    ### We control the different values.
    if (!is(x,"dgCMatrix")){
        stop("Not a dgCMatrix object of Matrix package")
    }
        
    if (!isMultiplex(MultiplexObject)) {
        stop("Not a Multiplex object")
    }
        
    L <- MultiplexObject$Number_of_Layers
    N <- MultiplexObject$Number_of_Nodes
        
    Seeds <- as.character(Seeds)
    if (length(Seeds) < 1 | length(Seeds) >= N){
        stop("The length of the vector containing the seed nodes is not 
            correct")
    } else {
        if (!all(Seeds %in% MultiplexObject$Pool_of_Nodes)){
            stop("Some of the seeds are not nodes of the network")
        }
    }
        
    if (r >= 1 || r <= 0) {
        stop("Restart partameter should be between 0 and 1")
    }
        
    if(missing(tau)){
        tau <- rep(1,L)/L
    } else {
        tau <- as.numeric(tau)
        if (sum(tau)/L != 1) {
            stop("The sum of the components of tau divided by the number of 
                layers should be 1")
        }
    }
        
    ## We define the threshold and the number maximum of iterations for
    ## the random walker.
    Threeshold <- 1e-10
    NetworkSize <- ncol(x)
        
    ## We initialize the variables to control the flux in the RW algo.
    residue <- 1
    iter <- 1
        
    ## We compute the scores for the different seeds.
    Seeds_Score <- get.seed.scoresMultiplex(Seeds,L,tau)
        
    ## We define the prox_vector(The vector we will move after the first RWR
    ## iteration. We start from The seed. We have to take in account
    ## that the walker with restart in some of the Seed genes, depending on
    ## the score we gave in that file).
    prox_vector <- matrix(0,nrow = NetworkSize,ncol=1)
        
    prox_vector[which(colnames(x) %in% Seeds_Score[,1])] <- (Seeds_Score[,2])
        
    prox_vector  <- prox_vector/sum(prox_vector)
    restart_vector <-  prox_vector
        
    while(residue >= Threeshold){
            
        old_prox_vector <- prox_vector
        prox_vector <- (1-r)*(x %*% prox_vector) + r*restart_vector
        residue <- sqrt(sum((prox_vector-old_prox_vector)^2))
        iter <- iter + 1;
    }
        
    NodeNames <- character(length = N)
    Score = numeric(length = N)
        
    rank_global <- data.frame(NodeNames = NodeNames, Score = Score)
    rank_global$NodeNames <- gsub("_1", "", row.names(prox_vector)[seq_len(N)])
    rank_global$Score <- geometric.mean(as.vector(prox_vector[,1]),L,N)
        
    ## We sort the genes according to their score.
    Global_results <- rank_global[with(rank_global, order(-Score, NodeNames)), ]
        
    ### We remove the seed genes from the Ranking and we write the results.
    Global_results_NoSeeds <- 
        Global_results[which(!Global_results$NodeNames %in% Seeds),]
        
    rownames(Global_results_NoSeeds) <- c()
        
    RWRM_ranking <- list(RWRM_Results = Global_results_NoSeeds,
        Seed_Nodes = Seeds)
        
    class(RWRM_ranking) <- "RWRM_Results"
    return(RWRM_ranking)
}

#' @method print RWRM_Results
#' @export
print.RWRM_Results <- function(x,...)
{
    cat("Top 10 ranked Nodes:\n")
    print(head(x$RWRM_Results,10))
    cat("\nSeed Nodes used:\n")
    print(x$Seed_Nodes)
}

## Roxy Documentation comments
#' Computes the transition matrix of a multiplex and heterogeneous network
#'
#' \code{compute.transition.matrix} is a function to compute the transition
#' matrix of a multiplex heterogeneous network provided as a \code{MultiplexHet}
#' object.
#'
#' @usage compute.transition.matrix(x,lambda = 0.5, delta=0.5)
#'
#' @details We clarify the role of the different parameters in this point:
#' \itemize{
#' \item \code{lambda}: For a given node, if a bipartite association exists,
#' the particle can either jump between networks or stay in the current
#' graph with a probability given by this parameter. The closer lambda is to
#' one, the higher is the probability of jumping between networks following
#' bipartite interactions.
#' \item \code{delta}: This parameter sets the probability to change between
#' layers at the next step. If delta = 0, the particle will always remain
#' in the same layer after a non-restart iteration. On the other hand, if
#' delta = 1, the particle will always change between layers, therefore
#' not following the specific edges of each layer.
#' }
#'
#' @param x A \code{MultiplexHet} object describing a multiplex and
#' heterogeneous network generated by the function
#' \code{create.multiplexHet}.
#'
#' @param lambda A numeric value between 0 and 1. It sets the probability of
#' jumping within a network or change to the other network of the heterogeneous
#' system. It is set by default to 0.5. See more details below.
#'
#' @param delta A numeric value between 0 and 1. It sets the probability
#' of performing inter-layer versus intra-layer transitions. It is set by
#' default to 0.5. See more details below.
#'
#' @return A square sparse transition matrix created with the \code{Matrix}
#' package. It is the transition matrix for the Random Walk with Restart on
#' Multiplex and Heterogeneous networks algorithm.
#'
#' @seealso \code{\link{create.multiplexHet},
#' \link{compute.adjacency.matrix}}
#'
#' @author Alberto Valdeolivas Urbelz \email{alvaldeolivas@@gmail.com}
#'
#' @examples
#' m1 <- igraph::graph(c(1,2,1,3,2,3), directed = FALSE)
#' m2 <- igraph::graph(c(1,3,2,3,3,4,1,4), directed = FALSE)
#' multiObject <- create.multiplex(m1,m2)
#' h1 <- igraph::graph(c("A","C","B","E","E","D","E","C"), directed = FALSE)
#' bipartite_relations <- data.frame(m=c(1,3),h=c("A","E"))
#' multiHetObject <- 
#'     create.multiplexHet(multiObject, h1,bipartite_relations)
#' compute.transition.matrix(multiHetObject)
#'
#'@import igraph
#'@import Matrix
#'@export
compute.transition.matrix <- function(x,lambda = 0.5, delta=0.5)
{
    if (!isMultiplexHet(x)) {
        stop("Not a Multiplex Heterogeneous object")
    }
    
    if (delta > 1 || delta <= 0) {
        stop("Delta should be between 0 and 1")
    }
    
    if (lambda > 1 || lambda <= 0) {
        stop("Lambda should be between 0 and 1")
    }
    
    N <- x$Number_of_Nodes_Multiplex
    L <- x$Number_of_Layers
    M <- x$Number_of_Nodes_Second_Network
    SupraBipartiteMatrix <- x$Bipartite_Matrix_Expanded
    
    message("Computing adjacency matrix of the Multiplex network...")
    AdjMatrix_MulNet <- compute.adjacency.matrix(x,delta)
    
    message("Computing adjacency matrix of the second network component...")
    ## We have to sort the adjacency matrix
    AdjMatrix_HetNet <- as_adjacency_matrix(x$Second_Network)
    AdjMatrix_HetNet <- 
        AdjMatrix_HetNet[order(rownames(AdjMatrix_HetNet)),
            order(colnames(AdjMatrix_HetNet))]
    
    
    ## Transition Matrix for the inter-subnetworks links
    message("Computing inter-subnetworks transitions...")
    Transition_Multiplex_SecondNet <- 
        get.transition.multiplex.secondNet(N,L,M, SupraBipartiteMatrix,lambda)
    Transition_SecondNet_Multiplex <- 
        get.transition.secondNet.multiplex(N,L,M, SupraBipartiteMatrix,lambda)
    
    ## Transition Matrix for the intra-subnetworks links
    message("Computing intra-subnetworks transitions...")
    Transition_Multiplex_Network <- 
        get.transition.multiplex(N,L, lambda,AdjMatrix_MulNet,
            SupraBipartiteMatrix)
    Transition_SecondNet_Network <- 
        get.transition.secondNet(M,lambda,AdjMatrix_HetNet,SupraBipartiteMatrix)
    
    ## We generate the global transiction matrix and we return it.
    message("Combining inter e intra layer probabilities into the global 
        Transition Matix")
    Transition_Multiplex_Heterogeneous_Matrix_1 <-
        cbind(Transition_Multiplex_Network, Transition_Multiplex_SecondNet)
    Transition_Multiplex_Heterogeneous_Matrix_2 <-
        cbind(Transition_SecondNet_Multiplex, Transition_SecondNet_Network)
    Transition_Multiplex_Heterogeneous_Matrix <-
        rbind(Transition_Multiplex_Heterogeneous_Matrix_1,
              Transition_Multiplex_Heterogeneous_Matrix_2)
    
    return(Transition_Multiplex_Heterogeneous_Matrix)
}

## Roxy Documentation comments
#' Performs Random Walk with Restart on a Multiplex and Heterogeneous
#' Network
#'
#' \code{Random.Walk.Restart.MultiplexHet} is a function to perform a Random
#' Walk with Restart on a Multiplex and Heterogeneous network (on a
#' \code{MultiplexHet} object). See more details about the algorithm below.
#'
#' @usage Random.Walk.Restart.MultiplexHet(...)
#'
#' @details Random Walk with Restart simulates an imaginary particle which
#' starts on a seed(s) node(s) and follows randomly the edges of a network. At
#' each step, there is a restart probability, r, meaning that the particle comes
#' back to the seed(s). The extension to multiplex networks allows the particle
#' to explore different monoplex networks (layers). At each step, the particle
#' can also jump to the same node in a different layer. The extension to
#' heterogeneous networks allows the particle to jump between nodes of different
#' nature thanks to bipartite relationships between them. We can combine both,
#' the multiplex and heterogeneous extension, by allowing the particle to jump
#' from a node in every layer of the multiplex network to the other network, and
#' the other way around.
#'
#' \itemize{
#' \item \code{Multiplex_Seed_Nodes}: A vector containing the name of the
#' different seed node(s) of the multiplex network. It's mandatory to provide at
#' least one seed (taking in account both types of seeds) The seed(s) node(s)
#' should belong to any of the layers of the multiplex network. The length of
#' this vector should be smaller than the total number of nodes in the multiplex
#' network.
#' \item \code{SecondNet_Seed_Nodes}: A vector containing the name of the
#' different seed node(s) of the second network. It's mandatory to provide at
#' least one seed (taking in account both types of seeds) The seed(s) node(s)
#' should belong to the second network. The length of this vector should be
#' smaller than the total number of nodes in the second network.
#' \item \code{r}: A numeric value representing the restart probability on the
#' seeds for the random walker. It must be between 0 and 1. It is set by default
#' to 0.7, which is the most common value in this kind of approaches. It means
#' that, at each step, the walker has a 70\% of probability of coming back to
#' one of the seeds.
#' \item \code{tau}: A numeric vector containing the probability of restarting
#' in the nodes of the different layers of the multiplex. In the example below,
#' we define the node 1 as the seed node. However, we can find this node in both
#' layers. Therefore, the walker can restart in any of these seed nodes. It is
#' a way to give different relevance (weight) to the different layers.
#' \item \code{eta}: A numeric value between 0 and 1 controlling the
#' probability of restarting in the nodes of each network. In the example below,
#' we define the node 1 as a multiplex seed node and "E" as a second network
#' seed node. Therefore, the walker can restart either in the seed 1 or in the
#' seed "E" with different probabilities (it is a way to give more relevance
#' to the different components of the heterogeneous system). If eta < 0.5
#' the particle will be more likely to restart in one of the multiplex seeds.
#' }
#'
#' @param x An object of the \code{Matrix} package describing the possible
#' transitions in a multiplex and heterogeneous network.
#' @param MultiplexHet_Object A \code{MultiplexHet} object generated by the
#' function \code{create.multiplexHet} representing a multiplex
#' and heterogeneous network.
#' @param Multiplex_Seed_Nodes A vector containing the names of the seeds of
#' the multiplex network for the Random Walk algorithm. See more details below.
#' @param SecondNet_Seed_Nodes A vector containing the names of the seeds of
#' the second network for the Random Walk algorithm. See more details below.
#' @param r A numeric value between 0 and 1. It sets the probability of
#' restarting to a seed node after each step. See more details below.
#' @param tau A vector containing the probability of restart on the seeds
#' of the different multiplex layers (layers weights).
#' It must have the same length than the number of layers of the multpiplex
#' network. The sum of its components divided by the number of layers must be 1.
#' See more details below.
#' @param eta A numeric value between 0 and 1. It controls the probability of
#' restarting in each network of the heterogeneous system (Multiplex or second
#' network). See more details below.
#' @param ... Further arguments passed to
#' \code{Random.Walk.Restart.MultiplexHet}
#'
#' @return A \code{RWRMH_Results} object. It contains two sorted rankings: The
#' first one contains the nodes of the multiplex network, except the seeds,
#' along with their score; The second one contains the nodes of the second
#' network, except the seeds, along with their score.
#' In addition, it contains two more fields describing the nodes of
#' different nature used as seeds.
#'
#' @seealso \code{\link{create.multiplexHet},
#' \link{compute.transition.matrix}, \link{Random.Walk.Restart.Multiplex}
#' \link{isRWRMH_Results}}
#'
#' @author Alberto Valdeolivas Urbelz \email{alvaldeolivas@@gmail.com}
#'
#' @examples
#' m1 <- igraph::graph(c(1,2,1,3,2,3), directed = FALSE)
#' m2 <- igraph::graph(c(1,3,2,3,3,4,1,4), directed = FALSE)
#' multiObject <- create.multiplex(m1,m2)
#' h1 <- igraph::graph(c("A","C","B","E","E","D","E","C"), directed = FALSE)
#' bipartite_relations <- data.frame(m=c(1,3),h=c("A","E"))
#' multiHetObject <- 
#'     create.multiplexHet(multiObject,h1,bipartite_relations)
#' MultiHetTranMatrix <- compute.transition.matrix(multiHetObject)
#' Multiplex_Seeds <- c(1)
#' SecondNet_Seeds <- c("E")
#' Random.Walk.Restart.MultiplexHet(MultiHetTranMatrix,
#'     multiHetObject,Multiplex_Seeds,SecondNet_Seeds)
#'
#'@import igraph
#'@import Matrix
#'@rdname Random.Walk.Restart.MultiplexHet
#'@export
Random.Walk.Restart.MultiplexHet <- function(...){
    UseMethod("Random.Walk.Restart.MultiplexHet")
}

#'@rdname Random.Walk.Restart.MultiplexHet
#'@export
Random.Walk.Restart.MultiplexHet.default <- 
    function(x, MultiplexHet_Object, Multiplex_Seed_Nodes,SecondNet_Seed_Nodes,
        r=0.7,tau,eta=0.5,...){
        
    ## We control the different values.
    if (!"dgCMatrix" %in% class(x)){
        stop("Not a dgCMatrix object of Matrix package")
    }
        
    if (!isMultiplexHet(MultiplexHet_Object)) {
        stop("Not a Multiplex Heterogeneous object")
    }
        
    L <- MultiplexHet_Object$Number_of_Layers
    N <- MultiplexHet_Object$Number_of_Nodes_Multiplex
    M <- MultiplexHet_Object$Number_of_Nodes_Second_Network
        
    All_nodes_Multiplex <- MultiplexHet_Object$Pool_of_Nodes
    All_nodes_SecondNet <- V(MultiplexHet_Object$Second_Network)$name
        
    MultiplexSeeds <- as.character(Multiplex_Seed_Nodes)
    SecondNet_Seed_Nodes <- as.character(SecondNet_Seed_Nodes)
        
    if (length(MultiplexSeeds) < 1 & length(SecondNet_Seed_Nodes) < 1){
        stop("You did not provided any seeds")
    } else {
        if (length(MultiplexSeeds) >= N | length(SecondNet_Seed_Nodes) >= M){
            stop("The length of some of the vectors containing the seed nodes 
                is not correct")
        }  else {
            if (!all(MultiplexSeeds %in% All_nodes_Multiplex)){
                stop("Some of the multiplex seeds are not nodes of the 
                    multiplex network")
            } else {
                if (!all(SecondNet_Seed_Nodes %in% All_nodes_SecondNet)){
                    stop("Some of the second network seeds are not nodes of 
                        this network")
                }
            }
        }
    }
        
    if (r >= 1 || r <= 0) {
        stop("Restart partameter should be between 0 and 1")
    }
        
    if (eta >= 1 || eta <= 0) {
        stop("Eta partameter should be between 0 and 1")
    }
        
    if(missing(tau)){
        tau <- rep(1,L)/L
    } else {
        tau <- as.numeric(tau)
        if (sum(tau)/L != 1) {
            stop("The sum of the components of tau divided by the number of 
                layers should be 1")}
        }
        
    ## We define the threshold and the number maximum of iterations
    ## for the random walker.
    Threeshold <- 1e-10
    NetworkSize <- ncol(x)
        
    ## We initialize the variables to control the flux in the RW algo.
    residue <- 1
    iter <- 1
        
    ## We compute the scores for the different seeds.
    Seeds_Score <- 
        get.seed.scores.multHet(Multiplex_Seed_Nodes, SecondNet_Seed_Nodes,eta,
            L,tau)
        
    ## We define the prox_vector(The vector we will move after the first
    ## RWR iteration. We start from The seed. We have to take in account
    ## that the walker with restart in some of the Seed genes,
    ## depending on the score we gave in that file).
    prox_vector <- matrix(0,nrow = NetworkSize,ncol=1)
        
    prox_vector[which(colnames(x) %in% Seeds_Score[,1])] <- (Seeds_Score[,2])
        
    prox_vector  <- prox_vector/sum(prox_vector)
    restart_vector <-  prox_vector
        
    while(residue >= Threeshold){
            
        old_prox_vector <- prox_vector
        prox_vector <- (1-r)*(x %*% prox_vector) + r*restart_vector
        residue <- sqrt(sum((prox_vector-old_prox_vector)^2))
        iter <- iter + 1;
    }
        
    ## We remove the seed genes from the ranking and we sort by score
    final_rank_multiplex_nodes <- 
        rank.nodes.multiplex(N,L,prox_vector,Multiplex_Seed_Nodes)
    final_rank_secondNet_nodes <- 
        rank.nodes.secondNet(N,L,M,prox_vector,SecondNet_Seed_Nodes)
        
    RWRMH_ranking <- list(RWRMH_Results_MultiplexNodes =
        final_rank_multiplex_nodes, Multiplex_Seed_Nodes = Multiplex_Seed_Nodes,
            RWRMH_Results_SecondNetNodes = final_rank_secondNet_nodes,
                SecondNet_Seed_Nodes = SecondNet_Seed_Nodes)
        
    class(RWRMH_ranking) <- "RWRMH_Results"
    return(RWRMH_ranking)
}

#' @method print RWRMH_Results
#' @export
print.RWRMH_Results <- function(x,...)
{
    cat("Top 10 ranked Multiplex nodes:\n")
    print(head(x$RWRMH_Results_MultiplexNodes,10))
    cat("\nMultiplex Seed Nodes used:\n")
    print(x$Multiplex_Seed_Nodes)
    cat("Top 10 ranked Second Network Nodes:\n")
    print(head(x$RWRMH_Results_SecondNetNodes,10))
    cat("\nSecond Network Seed Nodes used:\n")
    print(x$SecondNet_Seed_Nodes)
}
