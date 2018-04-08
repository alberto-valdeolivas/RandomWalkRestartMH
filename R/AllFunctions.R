## R Code for the Random Walk with Restart Package (RandomWalkRestartMH).

## Functions to check the objects used within the package.

## Roxy Documentaiton comments
#' Is this R object a Multiplex object?
#'
#' A Multiplex object is an R object generated as the result of calling
#' the function \code{create.multiplex}
#'
#' \code{is.multiplex(x)} checks whether an R object is Multiplex.
#'
#' @param x An R object
#'
#' @return A logical constant, \code{TRUE} if argument \code{x} is
#' a Mutiplex object.
#'
#' @seealso \code{\link{create.multiplex}}, \code{\link{is.multiplex.het}}
#'
#' @author Alberto Valdeolivas Urbelz \email{alvaldeolivas@@gmail.com}
#'
#' @examples
#' m1 <- igraph::graph(c(1,2,1,3,2,3), directed = FALSE)
#' m2 <- igraph::graph(c(1,3,2,3,3,4,1,4), directed = FALSE)
#' multiObject <- create.multiplex(m1,m2)
#' is.multiplex(multiObject)
#' is.multiplex(m1)
#'
#'@export

is.multiplex <- function (x)
{
    "Multiplex" %in% class(x)
}


## Roxy Documentaiton comments
#' Is this R object a Multiplex Heterogeneous object?
#'
#' A Multiplex Heterogeneous object is an R object generated as the result
#' of calling the function \code{create.multiplexHet}
#'
#' \code{is.multiplex.het(x)} checks whether an R object is MultiplexHet
#'
#' @param x An R object
#'
#' @return A logical constant, \code{TRUE} if argument \code{x} is
#' a MultiplexHet object.
#'
#' @seealso \code{\link{create.multiplexHet}},
#' \code{\link{is.multiplex}}
#'
#' @author Alberto Valdeolivas Urbelz \email{alvaldeolivas@@gmail.com}
#'
#' @examples
#' m1 <- igraph::graph(c(1,2,1,3,2,3), directed = FALSE)
#' m2 <- igraph::graph(c(1,3,2,3,3,4,1,4), directed = FALSE)
#' multiObject <- create.multiplex(m1,m2)
#' h1 <- igraph::graph(c("A","C","B","E","E","D","E","C"), directed = FALSE)
#' bipartite_relations <- data.frame(m=c(1,3),h=c("A","E"))
#' multiHetObject <- create.multiplexHet(multiObject,
#'     h1,bipartite_relations)
#' is.multiplex.het(multiHetObject)
#' is.multiplex.het(h1)
#'
#'@export
is.multiplex.het <- function (x)
{
    "MultiplexHet" %in% class(x)
}

## Roxy Documentaiton comments
#' Is this R object a RWR on Multiplex object (Results of the RWR-M)?
#'
#' A RWR on Multiplex object is an R object generated as the result
#' of calling the function \code{Random.Walk.Restart.Multiplex}
#' (Results of the RWR-M)
#'
#' \code{is.RWRM.Results(x)} checks whether an R object is RWRM_Results
#'
#' @param x An R object
#'
#' @return A logical constant, \code{TRUE} if argument \code{x} is
#' a RWRM_Results object.
#'
#' @seealso \code{\link{Random.Walk.Restart.Multiplex}},
#' \code{\link{is.RWRMH.Results}}
#'
#' @author Alberto Valdeolivas Urbelz \email{alvaldeolivas@@gmail.com}
#'
#' @examples
#' m1 <- igraph::graph(c(1,2,1,3,2,3), directed = FALSE)
#' m2 <- igraph::graph(c(1,3,2,3,3,4,1,4), directed = FALSE)
#' multiObject <- create.multiplex(m1,m2)
#' AdjMatrix <- compute.adjacency.matrix(multiObject)
#' AdjMatrixNorm <- normalize.multiplex.adjacency(AdjMatrix)
#' Seed <- c(1)
#' RWR_MultiResults <- Random.Walk.Restart.Multiplex(AdjMatrixNorm, multiObject,
#'     Seed)
#' is.RWRM.Results(RWR_MultiResults)
#' is.RWRM.Results(m1)
#'
#'@export
is.RWRM.Results <- function (x)
{
    "RWRM_Results" %in% class(x)
}


## Roxy Documentaiton comments
#' Is this R object a RWR on Multiplex-Heterogeneous object (Results of the
#' RWR-MH)?
#'
#' A RWR on Multiplex Heterogeneous object is an R object generated as the
#' result of calling the function \code{Random.Walk.Restart.MultiplexHet}
#' (Results of the RWR-MH)
#'
#' \code{is.RWRMH.Results(x)} checks whether an R object is RWRMH_Results
#'
#' @param x An R object
#'
#' @return A logical constant, \code{TRUE} if argument \code{x} is
#' a RWRMH_Results object.
#'
#' @seealso \code{\link{Random.Walk.Restart.MultiplexHet}},
#' \code{\link{is.RWRM.Results}}
#'
#' @author Alberto Valdeolivas Urbelz \email{alvaldeolivas@@gmail.com}
#'
#' @examples
#' m1 <- igraph::graph(c(1,2,1,3,2,3), directed = FALSE)
#' m2 <- igraph::graph(c(1,3,2,3,3,4,1,4), directed = FALSE)
#' multiObject <- create.multiplex(m1,m2)
#' h1 <- igraph::graph(c("A","C","B","E","E","D","E","C"), directed = FALSE)
#' bipartite_relations <- data.frame(m=c(1,3),h=c("A","E"))
#' multiHetObject <- create.multiplexHet(multiObject,
#'     h1,bipartite_relations)
#' MultiHetTranMatrix <- compute.transition.matrix(multiHetObject)
#' Multiplex_Seeds <- c(1)
#' SecondNet_Seeds <- c("E")
#' RWR_MultiHetResults <- Random.Walk.Restart.MultiplexHet(MultiHetTranMatrix,
#'     multiHetObject,Multiplex_Seeds,SecondNet_Seeds)
#' is.RWRMH.Results(RWR_MultiHetResults)
#' is.RWRMH.Results(m1)
#'
#'@export
is.RWRMH.Results <- function (x)
{
    "RWRMH_Results" %in% class(x)
}

## Geometric mean: R computation.

geometric.mean <- function(Scores, L, N) {

    FinalScore <- numeric(length = N)

    for (i in 1:N){
        FinalScore[i] <- prod(Scores[seq(from = i, to = N*L, by=N)])^(1/L)
    }

    return(FinalScore)
}



## Functions to perform Random Walk with Restart on Multiplex Networks.

simplify.layers <- function(Input_Layer){

    ## Undirected Graphs
    Layer <- as.undirected(Input_Layer, mode = c("collapse"),
                            edge.attr.comb = igraph_opt("edge.attr.comb"))

    ## Unweighted Graphs
    if (is_weighted(Layer)){
        Layer <- delete_edge_attr(Layer, "weight")
    }

    ## Simple Graphs
    Layer <- igraph::simplify(Layer, remove.multiple = TRUE,
                                remove.loops = TRUE)

    return(Layer)
}

add.missing.nodes <- function (Layers,Nr_Layers,NodeNames) {

    ## We generate a new list of layers.
    Layers_New <- vector("list", Nr_Layers)

    ## We add to each layer the missing nodes of the pool of nodes (all nodes).
    for (i in 1:Nr_Layers){
        Node_Names_Layer <- V(Layers[[i]])$name
        Missing_Nodes <- NodeNames[which(!NodeNames %in% Node_Names_Layer)]
        Layers_New[[i]] <- add_vertices(Layers[[i]] ,length(Missing_Nodes),
                                        name=Missing_Nodes)
    }
    return(Layers_New)
}

## Roxy Documentaiton comments
#' Create multiplex graphs from individual networks
#'
#' \code{create.multiplex} is a function to create a multiplex network
#' (\code{Multiplex} object) from up to 4 individual networks defined as
#' \code{igraph} objects. See more details about multiplex networks below.
#' If just one network is provided, a Multiplex object with one layer is
#' therefore created (A monoplex network).
#'
#' @usage create.multiplex(...)
#'
#' @details A multiplex network is a collection of layers (monoplex networks)
#' sharing the same nodes, but in which the edges represent relationships of
#' different nature.
#' The number of layers of a multiplex object can vary from 1 (monoplex
#' network) up to 4. Therefore, only the first layer is mandatory. We have
#' limited the number of layers to 4 in order to reduce computation times for
#' very large networks.
#' \itemize{
#'     \item \code{Layers_Name}: A vector contianing the name of the different
#'     layers. It's optional, but if provided the number of layers should match
#'     the length of this vector. Its elements should be in the same order than
#'     the graphs (\code{Layers_Name = c(L1_name, L2_name, L3_name, L4_name)}}
#'
#' @param L1 An igraph object describing a monoplex network. It will be
#' integrated as the first layer of the multiplex network.
#' @param L2 An igraph object describing a monoplex network. It will be
#' integrated as the second layer of the multiplex network. It's optional.
#' @param L3 An igraph object describing a monoplex network. It will be
#' integrated as the third layer of the multiplex network. It's optional.
#' @param L4 An igraph object describing a monoplex network. It will be
#' integrated as the fourth layer of the multiplex network. It's optional.
#' @param Layers_Name A vector containing the names of the different layers.
#' This name will be included as an attribute for all the edges of each network.
#' It's optional. See more details below.
#' @param ... Further arguments passed to \code{create.multiplex}
#'
#' @return A Multiplex object. It contains a list of the different graphs
#' integrating the multiplex network, the names and number of its nodes and the
#' number of layers.
#'
#' @seealso \code{\link{create.multiplexHet},\link{is.multiplex}}
#'
#' @author Alberto Valdeolivas Urbelz \email{alvaldeolivas@@gmail.com}
#'
#' @examples
#' m1 <- igraph::graph(c(1,2,1,3,2,3), directed = FALSE)
#' m2 <- igraph::graph(c(1,3,2,3,3,4,1,4), directed = FALSE)
#' multiObject <- create.multiplex(m1,m2)
#'
#'@import igraph
#'@rdname create.multiplex
#'@export
create.multiplex <- function(...){
    UseMethod("create.multiplex")
}

#' @rdname create.multiplex
#' @export
create.multiplex.default <- function(L1, L2=NULL, L3=NULL, L4=NULL,
    Layers_Name,...)
{

    ## We check that all the input layers are igraph objects.
    if (!is_igraph(L1)) {
        stop("Layer 1 is not a igraph object")
    }   else {
        Layer_1 <- simplify.layers(L1)
        if (is.null(V(L1)$name)) {
            Layer_1 <- set_vertex_attr(L1,"name",value=seq(1,vcount(L1),by=1))
        }
    }

    if (!is.null(L2)) {
        if (!is_igraph(L2)) {
            stop("Layer 2 is not a igraph object")
        } else {
            Layer_2 <- simplify.layers(L2)
            if (is.null(V(L2)$name)) {
                Layer_2 <-
                    set_vertex_attr(L2,"name",value=seq(1,vcount(L2),by=1))
            }
        }
    } else {
        Layer_2 <- make_empty_graph(n = 0, directed = FALSE)
    }

    if (!is.null(L3)) {
        if (!is_igraph(L3)) {
            stop("Layer 3 is not a igraph object")
        } else {
            Layer_3 <- simplify.layers(L3)
            if (is.null(V(L3)$name)) {
                Layer_3 <-
                    set_vertex_attr(L3,"name",value=seq(1,vcount(L3),by=1))
            }
        }
    } else {
        Layer_3 <- make_empty_graph(n = 0, directed = FALSE)
    }
    if (!is.null(L4)) {
        if (!is_igraph(L4)) {
            stop("Layer 4 is not a igraph object")
        } else {
            Layer_4 <- simplify.layers(L4)
            if (is.null(V(L4)$name)) {
                Layer_4 <-
                    set_vertex_attr(L4,"name",value=seq(1,vcount(L4),by=1))
            }
        }
    } else {
        Layer_4 <- make_empty_graph(n = 0, directed = FALSE)
    }

    ## We get a pool of nodes (Nodes in any of the layers.)
    Pool_of_Nodes <- sort(unique(c(V(Layer_1)$name,V(Layer_2)$name,
                                    V(Layer_3)$name,V(Layer_4)$name)))
    Number_of_Nodes <- length(Pool_of_Nodes)
    Number_of_Layers <- 1 + as.numeric(!is.null(L2)) +
        as.numeric(!is.null(L3)) +  as.numeric(!is.null(L4))

    if(missing(Layers_Name)){
        Layers_Name <-character()
        for (i in 1:Number_of_Layers){
            Layers_Name <- c(Layers_Name,paste("Layer_", i,":",sep="",
                                                collapse = ""))
        }
    } else {
        if (!is.character(Layers_Name)) {stop("The name of the layer should be
                                                a vector of characters")}
        if (length(Layers_Name) != Number_of_Layers) {stop("The length of the
                    Layer Name vector must be equal to the Number of Layers")}
    }

    Layer_List <- list(Layer_1,Layer_2,Layer_3,Layer_4)
    Layer_List <- add.missing.nodes(Layer_List, Number_of_Layers,Pool_of_Nodes)

    # We set the attributes of the layer
    for (i in 1:Number_of_Layers){
        Layer_List[[i]] <- set_edge_attr(Layer_List[[i]], "type",
                                            E(Layer_List[[i]]),
                                            value = Layers_Name[i])
    }

    names(Layer_List) <- Layers_Name

    MultiplexObject <- c(Layer_List,list(Pool_of_Nodes=Pool_of_Nodes,
                                Number_of_Nodes_Multiplex=Number_of_Nodes,
                                Number_of_Layers=Number_of_Layers))
    class(MultiplexObject) <- "Multiplex"
    return(MultiplexObject)

}

#' @method print Multiplex
#' @export
print.Multiplex <- function(x,...)
{
    cat("Number of Layers:\n")
    print(x$Number_of_Layers)
    cat("\nNumber of Nodes:\n")
    print(x$Number_of_Nodes)
    for (i in 1:x$Number_of_Layers){
        cat("\n")
        print(x[[i]])
    }
}

## Roxy Documentaiton comments
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
    if (!is.multiplex(x) & !is.multiplex.het(x)) {
        stop("Not a Multiplex or Multiplex Heterogeneous object")
    }
    if (delta > 1 || delta <= 0) {
        stop("Delta should be between 0 and 1")
    }

    N <- x$Number_of_Nodes_Multiplex
    L <- x$Number_of_Layers
    ## IDEM_MATRIX.
    Idem_Matrix <- Matrix::Diagonal(N, x = 1)

    SupraAdjacencyMatrix <- Matrix::Matrix(0,ncol=N*L,nrow=N*L,sparse = TRUE)

    Col_Node_Names <- character()
    Row_Node_Names <- character()

    for (i in 1:L){
        Adjacency_Layer <-  as_adjacency_matrix(x[[i]],sparse = TRUE)

        ## We order the matrix by the node name. This way all the matrix will
        ## have the same. Additionally we include a label with the layer number
        ## for each node name.
        Adjacency_Layer <- Adjacency_Layer[order(rownames(Adjacency_Layer)),
                                            order(colnames(Adjacency_Layer))]
        Layer_Col_Names <- paste(colnames(Adjacency_Layer),i,sep="_")
        Layer_Row_Names <- paste(rownames(Adjacency_Layer),i,sep="_")
        Col_Node_Names <- c(Col_Node_Names,Layer_Col_Names)
        Row_Node_Names <- c(Row_Node_Names,Layer_Row_Names)


        ## We fill the diagonal blocks with the transition probability
        ## within a layer
        Position_ini_row <- 1 + (i-1)*N
        Position_end_row <- N + (i-1)*N
        SupraAdjacencyMatrix[(Position_ini_row:Position_end_row),
                                (Position_ini_row:Position_end_row)] <-
                                (1-delta)*(Adjacency_Layer)

        ## We fill the off-diagonal blocks with the transition probability
        ## among layers.
        for (j in 1:L){
            Position_ini_col <- 1 + (j-1)*N
            Position_end_col <- N + (j-1)*N
            if (j != i){
                SupraAdjacencyMatrix[(Position_ini_row:Position_end_row),
                                    (Position_ini_col:Position_end_col)] <-
                                    (delta/(L-1))*Idem_Matrix
            }
        }
    }

    rownames(SupraAdjacencyMatrix) <- Row_Node_Names
    colnames(SupraAdjacencyMatrix) <- Col_Node_Names

    SupraAdjacencyMatrix <- as(SupraAdjacencyMatrix, "dgCMatrix")
    return(SupraAdjacencyMatrix)
}

get.seed.scores.multiplex <- function(Seeds,Number_Layers,tau) {

    Seeds_Seeds_Scores <- numeric(length = length(Seeds)*Number_Layers)
    Seed_Seeds_Layer_Labeled <- character(length = length(Seeds)*Number_Layers)

    Current_Seed_Labeled <- character()

    for (k in 1:Number_Layers){
        Current_Seed_Labeled <- c(Current_Seed_Labeled,
                                    paste(Seeds[k],k,sep="_",collapse = ""))
        for (j in 1:length(Seeds)){
            Seed_Seeds_Layer_Labeled[((k-1)*length(Seeds))+ j] <-
                                paste(Seeds[j],k,sep="_",collapse = "")
            Seeds_Seeds_Scores[((k-1)*length(Seeds))+ j] <-
                                tau[k]/length(Seeds)
        }
    }

    Seeds_Score <- data.frame(Seeds_ID = Seed_Seeds_Layer_Labeled,
                                Score = Seeds_Seeds_Scores,
                                stringsAsFactors = FALSE)

    return(Seeds_Score)
}

## Roxy Documentaiton comments
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
    if (!"dgCMatrix" %in% class(x)){
        stop("Not a dgCMatrix object of package Matrix")
    }

    Adj_Matrix_Norm <- t(t(x)/(Matrix::colSums(x, na.rm = FALSE, dims = 1,
                                                sparseResult = FALSE)))

    return(Adj_Matrix_Norm)
}


## Roxy Documentaiton comments
#' Performs Random Walk with Restart on a Multiplex Network
#'
#' \code{Random.Walk.Restart.Multiplex} is a function to perform a Random Walk
#' with Restart on a Multiplex network (on a \code{Multiplex} object). See
#' more details about the algorithm below.
#'
#' @usage Random.Walk.Restart.Multiplex(...)
#'
#' @details Random Walk with Restart simulates an imaginary particle which
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
#' should belong to any of the layers The lenght of this vector should be
#' smaller than the total number of nodes in the multiplex network.
#' \item \code{r}: A numeric value representing the restart probability on the
#' seeds for the random walker. It must be between 0 and 1. It is set by default
#' to 0.7, which is the most used value in this kind of approaches. It means
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
#' \link{is.RWRM.Results}, \link{Random.Walk.Restart.MultiplexHet}}
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
Random.Walk.Restart.Multiplex.default <- function(x, MultiplexObject, Seeds,
                                                    r=0.7,tau,...){

    ### We control the different values.
    if (!"dgCMatrix" %in% class(x)){
        stop("Not a dgCMatrix object of package Matrix")
    }

    if (!is.multiplex(MultiplexObject)) {
        stop("Not a Multiplex object")
    }

    L <- MultiplexObject$Number_of_Layers
    N <- MultiplexObject$Number_of_Nodes

    Seeds <- as.character(Seeds)
    if (length(Seeds) < 1 | length(Seeds) >= N){
        stop("The lenght of the vector containing the seed
            nodes is not correct")
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
        if (sum(tau)/L != 1) {stop("The sum tau components divided by
                                    the number of layers must be 1")}
    }

    ## We define the threshold and the number maximum of iterations for
    ## the random walker.
    Threeshold <- 1e-10
    NetworkSize <- ncol(x)

    ## We initialize the variables to control the flux in the RW algo.
    residue <- 1
    iter <- 1

    ## We compute the scores for the different seeds.
    Seeds_Score <- get.seed.scores.multiplex(Seeds,L,tau)

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
    rank_global$NodeNames <- gsub("_1", "", row.names(prox_vector)[1:N])
    rank_global$Score <- geometric.mean(as.vector(prox_vector[,1]),L,N)

    ## We sort the genes according to their score.
    Global_results <- rank_global[with(rank_global, order(-Score, NodeNames)), ]

    ### We remove the seed genes from the Ranking and we write the results.
    Global_results_NoSeeds <- Global_results[which(!Global_results$NodeNames
                                                    %in% Seeds),]

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

## Functions to perform Random Walk with Restart on Multiplex-Heterogeneous
## Networks.

get.bipartite.graph <- function(Names_Mul, Names_Het, Nodes_relation,
                                Number_Nodes_1,Number_Nodes_2){

    Bipartite_matrix <- Matrix(data=0, nrow=Number_Nodes_1, ncol=Number_Nodes_2)
    Names_Mul_order <- sort(Names_Mul)
    Names_Het_order <- sort(Names_Het)
    rownames(Bipartite_matrix) <- Names_Mul_order
    colnames(Bipartite_matrix) <- Names_Het_order

    for (i in 1:Number_Nodes_1){
        current_node1 <- Names_Mul_order[i]
        current_node2 <-
            Nodes_relation[which(Nodes_relation[,1] == current_node1),2]


        if (length(current_node2) > 0){
            for (j in 1:length(current_node2)){
                if (!is.na(current_node2[j])){
                    ## We need to identify the position on the matrix.
                    index <-
                        which(colnames(Bipartite_matrix) %in%  current_node2[j])
                    Bipartite_matrix[i,index] <- 1
                }
            }
        }
    }
    return(Bipartite_matrix)
}

expand.bipartite.graph <- function(Number_Nodes_1,Number_Layers,Number_Nodes_2,
                                    Bipartite_matrix){

    Supra_Bipartite_Matrix <- Matrix(0,nrow=Number_Nodes_1*Number_Layers,
                                    ncol=Number_Nodes_2,sparse = TRUE)
    Row_Node_Names <- character()

    for (i in 1:Number_Layers){
        Layer_Row_Names <- paste(rownames(Bipartite_matrix),i,sep="_")
        Row_Node_Names <- c(Row_Node_Names,Layer_Row_Names)
        Position_ini_row <- 1 + (i-1)*Number_Nodes_1
        Position_end_row <- Number_Nodes_1 + (i-1)*Number_Nodes_1
        Supra_Bipartite_Matrix[(Position_ini_row:Position_end_row),] <-
                                                                Bipartite_matrix
    }

    rownames(Supra_Bipartite_Matrix) <- Row_Node_Names
    colnames(Supra_Bipartite_Matrix) <- colnames(Bipartite_matrix)

    return(Supra_Bipartite_Matrix)
}

## Roxy Documentaiton comments
#' Create multiplex heterogeneous graphs from individual networks
#'
#' \code{create.multiplexHet} is a function to create a multiplex
#' and heterogeneous network (\code{MultiplexHet} object). It combines a
#' multiplex network composed from 1 (monoplex case) up to 4 layers with another
#' single network whose nodes are of different nature. See more details below.
#'
#' @usage create.multiplexHet(...)
#'
#' @details A multiplex network is a collection of layers (monoplex networks)
#' sharing the same nodes, but in which the edges represent relationships of
#' different nature. A heterogeneous network is composed of two single networks
#' where the nodes are of different nature.
#'
#' @param Multiplex_object A Multiplex network (\code{Multiplex} object)
#' generated by the function \code{create.multiplex}. This multiplex network
#' will be integrated as the first network of the heterogeneous network.
#' @param Het_graph An igraph object describing a monoplex network. It will be
#' integrated as the second network of the heterogeneous network.
#' @param Nodes_relations A data frame containing the relationships (bipartite
#' interactions) between the nodes of the multiplex network and the nodes of
#' the second network of the heterogeneous system. The data frame should contain
#' two columns: the first one with the nodes of the multiplex network; the
#' second one with the nodes of the second network. Every node should be present
#' in their corresponing network.
#' @param SecondNetworkName A vector containing the name for the second network
#' to be integrated as part of the heterogeneous network. This name will be
#' included as an attribute for all the edges of this network. It's optional and
#' its default value is SecondNetwork.
#' @param ... Further arguments passed to \code{create.multiplexHet}
#'
#' @return A Multiplex Heterogeneous object. It contains a list of the different
#' graphs integrating the multiplex network, the names and number of its nodes
#' and the number of layers. In addition, it contains the graph of the second
#' network integrating the heterogeneous network along with its number of
#' seeds. Finally, it contains a expanded bipartite adjacency matrix
#' describing the relations of the nodes in every layer of the multiplex network
#' with the nodes of the second network.
#'
#' @seealso \code{\link{create.multiplex},\link{is.multiplex.het}}
#'
#' @author Alberto Valdeolivas Urbelz \email{alvaldeolivas@@gmail.com}
#'
#' @examples
#' m1 <- igraph::graph(c(1,2,1,3,2,3), directed = FALSE)
#' m2 <- igraph::graph(c(1,3,2,3,3,4,1,4), directed = FALSE)
#' multiObject <- create.multiplex(m1,m2)
#' h1 <- igraph::graph(c("A","C","B","E","E","D","E","C"), directed = FALSE)
#' bipartite_relations <- data.frame(m=c(1,3),h=c("A","E"))
#' create.multiplexHet(multiObject,h1,bipartite_relations)
#'
#'@import igraph
#'@import Matrix
#'@rdname create.multiplexHet
#'@export
create.multiplexHet <- function(...) {
                                    UseMethod("create.multiplexHet")
}

#'@rdname create.multiplexHet
#'@export
create.multiplexHet.default <- function(Multiplex_object, Het_graph,
    Nodes_relations, SecondNetworkName,...)
{

    ## We check that all the arguments are correct
    print("checking input arguments...")
    if (!is.multiplex(Multiplex_object)) {
        stop("First element should be a multiplex object")
    }


    if (!is_igraph(Het_graph)) {
        stop("Second element should be an igraph object")
    }

    Het_nodes <- sort(V(Het_graph)$name)

    if (!is.data.frame(Nodes_relations)) {
        stop("Third element should be a data frame")
    } else {
        if (ncol(Nodes_relations) != 2) {
            stop("The data frame should contain two columns")
        } else {
            if (nrow(Nodes_relations) == 0) {
                stop("The data frame should contain any bipartite interaction")
            } else {
                names_1 <- unique(c(as.character(Nodes_relations[, 1])))
                names_2 <- unique(c(as.character(Nodes_relations[, 2])))
                if (!all(names_1 %in% Multiplex_object$Pool_of_Nodes)){
                    stop("Some of the nodes in the first column of the data
                        frame are not nodes of the multiplex network")
                } else {
                    if (!all(names_2 %in% Het_nodes)){
                        stop("Some of the nodes in the second column of the data
                            frame are not nodes of the second network")
                    }
                }
            }
        }
    }

    if(missing(SecondNetworkName)){
        SecondNetworkName <-c("SecondNetwork")
    } else {
        if (!is.character(SecondNetworkName)) {stop("The name of the layer
                                        should be a vector of characters")}
        if (length(SecondNetworkName) != 1) {stop("The name of the Second
                                Network name should be a vector of length 1")}
    }

    ## Het graph
    Het_graph <-  simplify.layers(Het_graph)
    Het_graph<- set_edge_attr(Het_graph, "type", E(Het_graph),
                                value = SecondNetworkName)

    M <- vcount(Het_graph)

    ## Multiplex graph
    Nodes_Multiplex <- Multiplex_object$Pool_of_Nodes
    N <- Multiplex_object$Number_of_Nodes
    L <- Multiplex_object$Number_of_Layers

    print("Generating bipartite matrix...")
    Bipartite_Matrix <- get.bipartite.graph(Nodes_Multiplex,
                                            Het_nodes,Nodes_relations,N, M)

    print("Expanding bipartite matrix to fit the multiplex network...")
    Supra_Bipartite_Matrix <- expand.bipartite.graph(N,L,M,Bipartite_Matrix)

    Multiplex_HetObject <- c(Multiplex_object,list(Second_Network=Het_graph,
                            Number_of_Nodes_Second_Network= M,
                            Bipartite_Matrix_Expanded =Supra_Bipartite_Matrix))

    class(Multiplex_HetObject) <- "MultiplexHet"
    return(Multiplex_HetObject)

}

#' @method print MultiplexHet
#' @export
print.MultiplexHet <- function(x,...)
{
    cat("Number of Layers:\n")
    print(x$Number_of_Layers)
    cat("\nNumber of Nodes Multiplex:\n")
    print(x$Number_of_Nodes_Multiplex)
    for (i in 1:x$Number_of_Layers){
        cat("\n")
        print(x[[i]])
    }
    cat("\nNumber of Nodes of the second network:\n")
    print(x$Number_of_Nodes_Second_Network)
    cat("\nSecond Network\n")
    print(x$Second_Network)
}

get.transition.multiplex.secondNet <- function(Number_Nodes_Multiplex,
            Number_Layers,Number_Nodes_secondNet,SupraBipartiteMatrix,lambda){

    Transition_Multiplex_SecondNet <- Matrix(0,
                    nrow=Number_Nodes_Multiplex*Number_Layers,
                    ncol=Number_Nodes_secondNet,sparse = TRUE)
    colnames(Transition_Multiplex_SecondNet) <- colnames(SupraBipartiteMatrix)
    rownames(Transition_Multiplex_SecondNet) <- rownames(SupraBipartiteMatrix)

    Col_Sum_Bipartite <- Matrix::colSums (SupraBipartiteMatrix, na.rm = FALSE,
                                            dims = 1,sparseResult = FALSE)

    for (j in 1:Number_Nodes_secondNet){
        if (Col_Sum_Bipartite[j] != 0){
            Transition_Multiplex_SecondNet[,j] <-
                (lambda*SupraBipartiteMatrix[,j]) /Col_Sum_Bipartite[j]
        }
    }

    return(Transition_Multiplex_SecondNet)
}

get.transition.secondNet.multiplex <- function(Number_Nodes_Multiplex,
            Number_Layers,Number_Nodes_secondNet,SupraBipartiteMatrix,lambda){

    Transition_SecondNet_Multiplex <- Matrix(0,nrow=Number_Nodes_secondNet,
                        ncol=Number_Nodes_Multiplex*Number_Layers,sparse = TRUE)

    colnames(Transition_SecondNet_Multiplex) <- rownames(SupraBipartiteMatrix)
    rownames(Transition_SecondNet_Multiplex) <- colnames(SupraBipartiteMatrix)

    Row_Sum_Bipartite <- Matrix::rowSums (SupraBipartiteMatrix, na.rm = FALSE,
                                            dims = 1,sparseResult = FALSE)

    for (i in 1:(Number_Nodes_Multiplex*Number_Layers)){
        if (Row_Sum_Bipartite[i] != 0){
            Transition_SecondNet_Multiplex[,i] <-
                (lambda*SupraBipartiteMatrix[i,])/Row_Sum_Bipartite[i]
        }
    }

    return(Transition_SecondNet_Multiplex)
}

get.transition.multiplex <- function(Number_Nodes_Multiplex,Number_Layers,
                            lambda,SupraAdjacencyMatrix,SupraBipartiteMatrix){

    Transition_Multiplex_Network <- Matrix(0,
                        nrow=Number_Nodes_Multiplex*Number_Layers,
                        ncol=Number_Nodes_Multiplex*Number_Layers,sparse = TRUE)

    rownames(Transition_Multiplex_Network) <- rownames(SupraAdjacencyMatrix)
    colnames(Transition_Multiplex_Network) <- colnames(SupraAdjacencyMatrix)

    Col_Sum_Multiplex <- Matrix::colSums(SupraAdjacencyMatrix,na.rm=FALSE,
                                        dims=1, sparseResult=FALSE)
    Row_Sum_Bipartite <- Matrix::rowSums (SupraBipartiteMatrix, na.rm = FALSE,
                                        dims = 1,sparseResult = FALSE)

    for (j in 1:(Number_Nodes_Multiplex*Number_Layers)){
        if(Row_Sum_Bipartite[j] != 0){
            Transition_Multiplex_Network[,j] <-
                ((1-lambda)*SupraAdjacencyMatrix[,j]) /Col_Sum_Multiplex[j]
        } else {
            Transition_Multiplex_Network[,j] <-
                SupraAdjacencyMatrix[,j] /Col_Sum_Multiplex[j]
        }
    }

    return(Transition_Multiplex_Network)
}


get.transition.secondNet <- function(Number_Nodes_secondNet,lambda,
                                    AdjMatrix,SupraBipartiteMatrix){

    Transition_Second_Network <- Matrix(0,nrow=Number_Nodes_secondNet,
                                ncol=Number_Nodes_secondNet,sparse = TRUE)

    rownames(Transition_Second_Network) <- rownames(AdjMatrix)
    colnames(Transition_Second_Network) <- colnames(AdjMatrix)

    Col_Sum_SecondNet <- Matrix::colSums (AdjMatrix,na.rm=FALSE,
                                        dims=1,sparseResult=FALSE)
    Col_Sum_Bipartite <- Matrix::colSums (SupraBipartiteMatrix, na.rm = FALSE,
                                        dims = 1,sparseResult = FALSE)

    for (j in 1:Number_Nodes_secondNet){
        if(Col_Sum_Bipartite[j] != 0){
            Transition_Second_Network[,j] <-
                ((1-lambda)*AdjMatrix[,j]) /Col_Sum_SecondNet[j]
        } else {
            Transition_Second_Network[,j] <- AdjMatrix[,j] /Col_Sum_SecondNet[j]
        }
    }

    return(Transition_Second_Network)
}

## Roxy Documentaiton comments
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
#' multiHetObject <- create.multiplexHet(multiObject,
#'     h1,bipartite_relations)
#' compute.transition.matrix(multiHetObject)
#'
#'@import igraph
#'@import Matrix
#'@export
compute.transition.matrix <- function(x,lambda = 0.5, delta=0.5)
{
    if (!is.multiplex.het(x)) {
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

    print("Computing adjacency matrix of the Multiplex network...")
    AdjMatrix_MulNet <- compute.adjacency.matrix(x,delta)

    print("Computing adjacency matrix of the second network component...")
    ## We have to sort the adjacency matrix
    AdjMatrix_HetNet <- as_adjacency_matrix(x$Second_Network)
    AdjMatrix_HetNet <- AdjMatrix_HetNet[order(rownames(AdjMatrix_HetNet)),
                                        order(colnames(AdjMatrix_HetNet))]


    ## Transition Matrix for the inter-subnetworks links
    print("Computing inter-subnetworks transitions...")
    Transition_Multiplex_SecondNet <- get.transition.multiplex.secondNet(N,L,M,
                                                    SupraBipartiteMatrix,lambda)
    Transition_SecondNet_Multiplex <- get.transition.secondNet.multiplex(N,L,M,
                                                    SupraBipartiteMatrix,lambda)

    ## Transition Matrix for the intra-subnetworks links
    print("Computing intra-subnetworks transitions...")
    Transition_Multiplex_Network <- get.transition.multiplex(N,L,
                                lambda,AdjMatrix_MulNet,SupraBipartiteMatrix)
    Transition_SecondNet_Network <- get.transition.secondNet(M,lambda,
                                AdjMatrix_HetNet,SupraBipartiteMatrix)

    ## We generate the global transiction matrix and we return it.
    print("Combining inter e intra layer into the global Transition Matix")
    Transition_Multiplex_Heterogeneous_Matrix_1 <-
        cbind(Transition_Multiplex_Network, Transition_Multiplex_SecondNet)
    Transition_Multiplex_Heterogeneous_Matrix_2 <-
        cbind(Transition_SecondNet_Multiplex, Transition_SecondNet_Network)
    Transition_Multiplex_Heterogeneous_Matrix <-
        rbind(Transition_Multiplex_Heterogeneous_Matrix_1,
                Transition_Multiplex_Heterogeneous_Matrix_2)

    return(Transition_Multiplex_Heterogeneous_Matrix)
}

get.seed.scores.multHet <- function(Multiplex_Seed_Nodes,SecondNet_Seed_Nodes,
                                    eta,L,tau) {


    n <- length(Multiplex_Seed_Nodes)
    m <- length(SecondNet_Seed_Nodes)

    if ((n != 0 && m!= 0)){

        Seed_Multiplex_Layer_Labeled <-character(length = n*L)
        Seeds_Multiplex_Scores <-numeric(length = n*L)

        Current_MultiplexNode_Labeled <- character()

        for (k in 1:L){
            Current_MultiplexNode_Labeled <-
                c(Seed_Multiplex_Layer_Labeled,paste(Multiplex_Seed_Nodes[k],
                                                    k,sep="_",collapse = ""))
            for (j in 1:n){
                Seed_Multiplex_Layer_Labeled[((k-1)*n)+ j] <-
                    paste(Multiplex_Seed_Nodes[j],k,sep="_",collapse = "")
                Seeds_Multiplex_Scores[((k-1)*n)+ j] <- ((1-eta) * tau[k])/n
            }
        }

        SecondNet_Seeds_Score <- eta/m

    } else {
        eta <- 1
        if (n == 0){
            Seed_Multiplex_Layer_Labeled <- character()
            Seeds_Multiplex_Scores <- numeric()
            SecondNet_Seeds_Score <- eta/m
        } else {

            Seed_Multiplex_Layer_Labeled <- character(length = n*L)
            Seeds_Multiplex_Scores <- numeric(length = n*L)

            Current_MultiplexNode_Labeled <- character()

            for (k in 1:L){
                Current_MultiplexNode_Labeled <-
                    c(Current_MultiplexNode_Labeled,
                    paste(Multiplex_Seed_Nodes[k],k,sep="_",collapse = ""))
                for (j in 1:n){
                    Seed_Multiplex_Layer_Labeled[((k-1)*n)+ j] <-
                        paste(Multiplex_Seed_Nodes[j],k,sep="_",collapse = "")
                    Seeds_Multiplex_Scores[((k-1)*n)+ j] <- tau[k]/n
                    SecondNet_Seeds_Score <- numeric()
                }
            }
        }
    }

    ## We prepare a data frame with the seeds.
    Seeds_Score <- data.frame(Seeds_ID = c(Seed_Multiplex_Layer_Labeled,
                                            SecondNet_Seed_Nodes),
                            Score = c(Seeds_Multiplex_Scores,
                                        rep(SecondNet_Seeds_Score,m)),
                            stringsAsFactors = FALSE)

    return(Seeds_Score)
}

rank.nodes.multiplex <- function(N, L, Results,Seeds){

    ## We sort the score to obtain the ranking of multiplex nodes and
    ## seconde network nodes.
    NodeNames <- character(length = N)
    Score <- numeric(length = N)

    nodes_multiplex_rank <- data.frame(NodeNames = NodeNames, Score = Score)
    nodes_multiplex_rank$NodeNames <- gsub("_1", "", row.names(Results)[1:N])

    ## We calculate the Geometric Mean among the proteins in the
    ## different layers.
    nodes_multiplex_rank$Score <- geometric.mean(as.vector(Results[,1]),L,N)

    nodes_multiplex_sort <- nodes_multiplex_rank[with(nodes_multiplex_rank,
                                                order(-Score, NodeNames)), ]

    ## We remove the seed genes from the Ranking
    nodes_multiplex_sort_NoSeeds <-
        nodes_multiplex_sort[which(!nodes_multiplex_sort$NodeNames %in% Seeds),]

    return(nodes_multiplex_sort_NoSeeds)
}

rank.nodes.secondNet <- function(N,L,M,Results,Seeds){

    ## rank_diseases
    SecondNet_node <- character(length = M)
    Score <- character(length = M)
    SecondNet_rank <- data.frame(SecondNet_node = SecondNet_node, Score = Score)
    SecondNet_rank$SecondNet_node <- row.names(Results)[((N*L)+1):nrow(Results)]
    SecondNet_rank$Score <- Results[((N*L)+1):nrow(Results),1]

    SecondNet_rank_sort <- SecondNet_rank[with(SecondNet_rank,
                                            order(-Score, SecondNet_node)), ]
    SecondNet_rank_sort_NoSeeds <-
    SecondNet_rank_sort[which(!SecondNet_rank_sort$SecondNet_node %in% Seeds),]

    return(SecondNet_rank_sort_NoSeeds)
}

## Roxy Documentaiton comments
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
#' should belong to any of the layers of the multiplex network. The lenght of
#' this vector should be smaller than the total number of nodes in the multiplex
#' network.
#' \item \code{SecondNet_Seed_Nodes}: A vector containing the name of the
#' different seed node(s) of the second network. It's mandatory to provide at
#' least one seed (taking in account both types of seeds) The seed(s) node(s)
#' should belong to the second network. The lenght of this vector should be
#' smaller than the total number of nodes in the second network.
#' \item \code{r}: A numeric value representing the restart probability on the
#' seeds for the random walker. It must be between 0 and 1. It is set by default
#' to 0.7, which is the most used value in this kind of approaches. It means
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
#' \link{is.RWRMH.Results}}
#'
#' @author Alberto Valdeolivas Urbelz \email{alvaldeolivas@@gmail.com}
#'
#' @examples
#' m1 <- igraph::graph(c(1,2,1,3,2,3), directed = FALSE)
#' m2 <- igraph::graph(c(1,3,2,3,3,4,1,4), directed = FALSE)
#' multiObject <- create.multiplex(m1,m2)
#' h1 <- igraph::graph(c("A","C","B","E","E","D","E","C"), directed = FALSE)
#' bipartite_relations <- data.frame(m=c(1,3),h=c("A","E"))
#' multiHetObject <- create.multiplexHet(multiObject,
#'     h1,bipartite_relations)
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
Random.Walk.Restart.MultiplexHet.default <- function(x, MultiplexHet_Object,
                Multiplex_Seed_Nodes,SecondNet_Seed_Nodes,r=0.7,tau,
                eta=0.5,...){

    ## We control the different values.
    if (!"dgCMatrix" %in% class(x)){
        stop("Not a dgCMatrix object of package Matrix")
    }

    if (!is.multiplex.het(MultiplexHet_Object)) {
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
            stop("The lenght of some of the vectors containing
                the seed nodes is not correct")
        }  else {
            if (!all(MultiplexSeeds %in% All_nodes_Multiplex)){
                stop("Some of the multiplex seeds are not
                    nodes of the multiplex network")
            } else {
                if (!all(SecondNet_Seed_Nodes %in% All_nodes_SecondNet)){
                    stop("Some of the second network seeds
                        are not nodes of this network")
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
        if (sum(tau)/L != 1) {stop("The sum tau components divided
                                    by the number of layers must be 1")}
    }


    ## We define the threshold and the number maximum of iterations
    ## for the random walker.
    Threeshold <- 1e-10
    NetworkSize <- ncol(x)

    ## We initialize the variables to control the flux in the RW algo.
    residue <- 1
    iter <- 1

    ## We compute the scores for the different seeds.
    Seeds_Score <- get.seed.scores.multHet(Multiplex_Seed_Nodes,
                                            SecondNet_Seed_Nodes,eta,L,tau)

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
    final_rank_multiplex_nodes <- rank.nodes.multiplex(N,L,prox_vector,
                                                        Multiplex_Seed_Nodes)
    final_rank_secondNet_nodes <- rank.nodes.secondNet(N,L,M,prox_vector,
                                                        SecondNet_Seed_Nodes)

    RWRMH_ranking <- list(RWRMH_Results_MultiplexNodes =
                                final_rank_multiplex_nodes,
                            Multiplex_Seed_Nodes = Multiplex_Seed_Nodes,
                            RWRMH_Results_SecondNetNodes =
                                final_rank_secondNet_nodes,
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


## Functions to generate a Network with the Top results of RWR-M and
## RWR-MH results.


## Roxy Documentaiton comments
#' Creates a Network with the top results of the Random Walk with restart on
#' a Multiplex Network
#'
#' \code{create.multiplexNetwork.topResults} is a function to create a network
#' from the top results of the Random Walk with Restart on Multiplex networks
#' algorithm (a \code{RWRM_Results} object).
#'
#' @usage create.multiplexNetwork.topResults(RWRM_Result_Object,
#'   MultiplexObject,k=25)
#'
#' @param RWRM_Result_Object A \code{RWRM_Results} object generated by the
#' function \code{Random.Walk.Restart.Multiplex} representing the results
#' of the Random Ralk with restart on the multiplex network described in the
#' following argument.
#' @param MultiplexObject A \code{Multiplex} object generated by the
#' function \code{create.multiplex} representing a multiplex network.
#' @param k A numeric value between 1 and 200. It is the number of top ranked
#' nodes to be included in the resulting multiplex network.
#'
#' @return An \code{igraph} object containing the top \code{k} ranked
#' multiplex nodes in the Random Walk with Restart on a Multiplex network
#' algorithm. We include all the possible types of interactions between pairs of
#' nodes according to the different layers of the multiplex network.
#'
#' @seealso \code{\link{create.multiplex}, \link{Random.Walk.Restart.Multiplex}
#' \link{is.RWRM.Results}, \link{create.multiplexHetNetwork.topResults}}
#'
#' @author Alberto Valdeolivas Urbelz \email{alvaldeolivas@@gmail.com}
#'
#' @examples
#' m1 <- igraph::graph(c(1,2,1,3,2,3), directed = FALSE)
#' m2 <- igraph::graph(c(1,3,2,3,3,4,1,4), directed = FALSE)
#' multiObject <- create.multiplex(m1,m2)
#' AdjMatrix <- compute.adjacency.matrix(multiObject)
#' AdjMatrixNorm <- normalize.multiplex.adjacency(AdjMatrix)
#' Seed <- c(1)
#' RWR_MultiResults <- Random.Walk.Restart.Multiplex(AdjMatrixNorm, multiObject,
#'     Seed)
#' create.multiplexNetwork.topResults(RWR_MultiResults,multiObject)
#'
#'@import igraph
#'@importFrom dnet dNetInduce
#'@export
create.multiplexNetwork.topResults <- function(RWRM_Result_Object,
                                                MultiplexObject,k=25) {

    if (!is.multiplex(MultiplexObject)) {
        stop("Not a Multiplex object")
    }

    if (!is.RWRM.Results(RWRM_Result_Object)){
        stop("Not Results for RWR-M")
    }

    k <- as.numeric(k)
    if (k  <= 0 || k >= 200) {
        stop("K should be between 0 and 200")
    }

    Multiplex_df <- as_data_frame(MultiplexObject[[1]])
    if (MultiplexObject$Number_of_Layers > 1){
        for (i in 2:MultiplexObject$Number_of_Layers){
            Multiplex_df <- rbind(Multiplex_df,
                                    as_data_frame(MultiplexObject[[i]]))
        }
    }

    Multiplex_Network <- graph.data.frame(Multiplex_df,directed=FALSE)
    Top_Results_Nodes <- RWRM_Result_Object$RWRM_Results$NodeNames[1:k]
    Query_Nodes <- c(RWRM_Result_Object$Seed_Nodes,Top_Results_Nodes)

    Induced_Network <- dNetInduce(g=Multiplex_Network, nodes_query=Query_Nodes,
                                knn=0, remove.loops=FALSE, largest.comp=FALSE)
    return(Induced_Network)
}


## Roxy Documentaiton comments
#' Creates a Network with the top results of the Random Walk with restart on
#' a Multiplex and Heterogeneous Network
#'
#' \code{create.multiplexHetNetwork.topResults} is a function to create a
#' network from the top results of the Random Walk with Restart on Multiplex and
#' Heterogeneous networks algorithm (a \code{RWRMH_Results} object).
#'
#' @usage create.multiplexHetNetwork.topResults(RWRMH_Results_Object,
#'   MultiplexHetObject, bipartite_relations, bipartite_name, k=25)
#'
#' @param RWRMH_Results_Object A \code{RWRMH_Results} object generated by the
#' function \code{Random.Walk.Restart.MultiplexHet} representing the results
#' of the Random Ralk with restart on the multiplex and heterogeneous network
#' described in the following argument.
#' @param MultiplexHetObject A \code{MultiplexHet} object generated by the
#' function \code{create.multiplexHet} representing a multiplex and
#' heterogeneous network.
#' @param bipartite_relations A data frame containing the relationships between
#' the nodes of the multiplex network and the nodes of the second network of the
#' heterogeneous system. The data frame should contain two columns: the first
#' one with the nodes of the multiplex network; the second one with the nodes of
#' the second network. Every node should be present in their corresponing
#' network.
#' @param bipartite_name A vector containing the name for the bipartite
#' relations to be integrated as part of the resulting network. It is included
#' as an attribute for all the bipartite edges of the resulting network. It's
#' optional and its default value is "bipartiteRelations".
#' @param k A numeric value between 1 and 200. It is the number of top ranked
#' nodes to be included in the resulting multiplex network.
#'
#' @return An \code{igraph} object containing the top \code{k} ranked
#' multiplex nodes and the top \code{k} ranked second network nodes in the
#' Random Walk with Restart on a Multiplex and Heterogeneous network algorithm.
#' We include all the possible types of interactions between pairs of
#' nodes according to the different layers of the multiplex network, the
#' bipartite interactions and the second network type of interactions.
#'
#' @seealso \code{\link{create.multiplexHet},
#' \link{is.RWRMH.Results}, \link{Random.Walk.Restart.MultiplexHet}
#' \link{create.multiplexNetwork.topResults}}
#'
#' @author Alberto Valdeolivas Urbelz \email{alvaldeolivas@@gmail.com}
#'
#' @examples
#' m1 <- igraph::graph(c(1,2,1,3,2,3), directed = FALSE)
#' m2 <- igraph::graph(c(1,3,2,3,3,4,1,4), directed = FALSE)
#' multiObject <- create.multiplex(m1,m2)
#' h1 <- igraph::graph(c("A","C","B","E","E","D","E","C"), directed = FALSE)
#' bipartite_relations <- data.frame(m=c(1,3),h=c("A","E"))
#' multiHetObject <- create.multiplexHet(multiObject,
#'     h1,bipartite_relations)
#' MultiHetTranMatrix <- compute.transition.matrix(multiHetObject)
#' Multiplex_Seeds <- c(1)
#' SecondNet_Seeds <- c("E")
#' RWR_MultiHetResults <- Random.Walk.Restart.MultiplexHet(MultiHetTranMatrix,
#'     multiHetObject,Multiplex_Seeds,SecondNet_Seeds)
#' create.multiplexHetNetwork.topResults(RWR_MultiHetResults,multiHetObject,
#'     bipartite_relations)
#'
#'@import igraph
#'@importFrom dnet dNetInduce
#'@export

create.multiplexHetNetwork.topResults <- function(RWRMH_Results_Object,
    MultiplexHetObject,bipartite_relations,bipartite_name, k=25) {

    if (!is.multiplex.het(MultiplexHetObject)) {
        stop("Not a Multiplex Heterogeneous object")
    }

    if (!is.RWRMH.Results(RWRMH_Results_Object)){
        stop("Not Results for RWR-MH")
    }

    if (!is.data.frame(bipartite_relations)) {
        stop("Third element should be a data frame")
    } else {
        if (ncol(bipartite_relations) != 2) {
            stop("The data frame should contain two columns")
        } else {
            if (nrow(bipartite_relations) == 0) {
                stop("The data frame should contain any bipartite interaction")
            } else {
                names_1 <- unique(c(as.character(bipartite_relations[, 1])))
                names_2 <- unique(c(as.character(bipartite_relations[, 2])))
                if (!all(names_1 %in% MultiplexHetObject$Pool_of_Nodes)){
                    stop("Some of the nodes in the first column of the data
                        frame are not nodes of the multiplex network")
                } else {
                    if (!all(names_2 %in%
                            V(MultiplexHetObject$Second_Network)$name)){
                        stop("Some of the nodes in the second column of the
                            data frame are not nodes of the second network")
                    }
                }
            }
        }
    }

    if(missing(bipartite_name)){
        bipartite_name <-c("bipartiteRelations")
    } else {
        if (!is.character(bipartite_name)) {stop("The name of the bipartite
                                                relations should be a vector
                                                of characters")}
        if (length(bipartite_name) != 1) {stop("The name of the bipartite
                                                relations should be a vector
                                                of length 1")}
    }

    k <- as.numeric(k)
    if (k  <= 0 || k >= 200) {
        stop("K should be between 0 and 200")
    }

    bipartite_relations$type <- bipartite_name
    colnames(bipartite_relations) <-c("from","to","type")

    Multiplex_df <- as_data_frame(MultiplexHetObject[[1]])
    if (MultiplexHetObject$Number_of_Layers > 1){
        for (i in 2:MultiplexHetObject$Number_of_Layers){
            Multiplex_df <- rbind(Multiplex_df,
                                    as_data_frame(MultiplexHetObject[[i]]))
        }
    }

    Second_Network_df <- as_data_frame(MultiplexHetObject$Second_Network)

    Multiplex_Heterogeneous_df <- rbind(Multiplex_df,Second_Network_df,
                                        bipartite_relations)


    Multiplex_Heterogeneous_Network <-
        graph.data.frame(Multiplex_Heterogeneous_df,directed=FALSE)

    Top_Results_MultiNodes <-
        RWRMH_Results_Object$RWRMH_Results_MultiplexNodes$NodeNames[1:k]
    Top_Results_SecondNetNodes <-
        RWRMH_Results_Object$RWRMH_Results_SecondNetNodes$SecondNet_node[1:k]

    Query_Nodes <- c(c(RWRMH_Results_Object$Multiplex_Seed_Nodes,
                        RWRMH_Results_Object$SecondNet_Seed_Nodes),
                        Top_Results_MultiNodes,
                        Top_Results_SecondNetNodes)

    Induced_Network <- dNetInduce(g=Multiplex_Heterogeneous_Network,
                                    nodes_query=Query_Nodes, knn=0,
                                    remove.loops=FALSE, largest.comp=FALSE)
    return(Induced_Network)
}
