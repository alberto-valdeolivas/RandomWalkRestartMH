## R Code for the Random Walk with Restart Package (RandomWalkRestartMH).

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Functions to create Multiplex and Multiplex-Heterogeneous objects.
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

## Roxy Documentation comments
#' Create multiplex graphs from individual networks
#'
#' \code{create.multiplex} is a function to create a multiplex network
#' (\code{Multiplex} object) from up to 6 individual networks defined as
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
#' network) up to 6. Therefore, only the first layer is mandatory. We have
#' limited the number of layers to 6 in order to reduce computation times for
#' very large networks.
#' \itemize{
#'     \item \code{Layers_Name}: A vector contianing the name of the different
#'     layers. It's optional, but if provided the number of layers should match
#'     the length of this vector. Its elements should be in the same order than
#'     the graphs (\code{Layers_Name = c(L1_name, L2_name, L3_name, L4_name,
#'     L5_name, L6_name)}}
#'
#' @param L1 An igraph object describing a monoplex network. It will be
#' integrated as the first layer of the multiplex network.
#' @param L2 An igraph object describing a monoplex network. It will be
#' integrated as the second layer of the multiplex network. It's optional.
#' @param L3 An igraph object describing a monoplex network. It will be
#' integrated as the third layer of the multiplex network. It's optional.
#' @param L4 An igraph object describing a monoplex network. It will be
#' integrated as the fourth layer of the multiplex network. It's optional.
#' @param L5 An igraph object describing a monoplex network. It will be
#' integrated as the fifth layer of the multiplex network. It's optional.
#' @param L6 An igraph object describing a monoplex network. It will be
#' integrated as the sixth layer of the multiplex network. It's optional.
#' @param Layers_Name A vector containing the names of the different layers.
#' This name will be included as an attribute for all the edges of each network.
#' It's optional. See more details below.
#' @param ... Further arguments passed to \code{create.multiplex}
#'
#' @return A Multiplex object. It contains a list of the different graphs
#' integrating the multiplex network, the names and number of its nodes and the
#' number of layers.
#'
#' @seealso \code{\link{create.multiplexHet},\link{isMultiplex}}
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
create.multiplex.default <- function(L1, L2=NULL, L3=NULL, L4=NULL,L5=NULL,
    L6=NULL, Layers_Name,...)
{
    
    ## We check that all the input layers are igraph objects.
    if (!is_igraph(L1)) {
        stop("Layer 1 is not an igraph object")
    }   else {
        Layer_1 <- simplify.layers(L1)
        if (is.null(V(L1)$name)) {
            Layer_1 <- set_vertex_attr(L1,"name",value=seq(1,vcount(L1),by=1))
        }
    }
    
    if (!is.null(L2)) {
        if (!is_igraph(L2)) {
            stop("Layer 2 is not an igraph object")
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
            stop("Layer 3 is not an igraph object")
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
            stop("Layer 4 is not an igraph object")
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
    
    if (!is.null(L5)) {
        if (!is_igraph(L5)) {
            stop("Layer 5 is not an igraph object")
        } else {
            Layer_5 <- simplify.layers(L5)
            if (is.null(V(L5)$name)) {
                Layer_5 <-
                    set_vertex_attr(L5,"name",value=seq(1,vcount(L5),by=1))
            }
        }
    } else {
        Layer_5 <- make_empty_graph(n = 0, directed = FALSE)
    }
    
    if (!is.null(L6)) {
        if (!is_igraph(L6)) {
            stop("Layer 6 is not an igraph object")
        } else {
            Layer_6 <- simplify.layers(L6)
            if (is.null(V(L6)$name)) {
                Layer_6 <-
                    set_vertex_attr(L6,"name",value=seq(1,vcount(L6),by=1))
            }
        }
    } else {
        Layer_6 <- make_empty_graph(n = 0, directed = FALSE)
    }
    
    ## We get a pool of nodes (Nodes in any of the layers.)
    Pool_of_Nodes <- sort(unique(c(V(Layer_1)$name,V(Layer_2)$name,
        V(Layer_3)$name,V(Layer_4)$name, V(Layer_5)$name, V(Layer_6)$name)))
    Number_of_Nodes <- length(Pool_of_Nodes)
    Number_of_Layers <- 1 + as.numeric(!is.null(L2)) +
        as.numeric(!is.null(L3)) +  as.numeric(!is.null(L4)) + 
        as.numeric(!is.null(L5)) + as.numeric(!is.null(L6))
    
    if(missing(Layers_Name)){
        Layers_Name <-Layers_Name <- paste0("Layer_", seq_len(Number_of_Layers))
    } else {
        if (!is.character(Layers_Name)) {
            stop("The names of the layers should be provided as a vector 
                 of characters")}
        if (length(Layers_Name) != Number_of_Layers) {
            stop("The length of the Layers Names vector should be equal to the 
                Number of Layers")}
        }
    
    Layer_List <- list(Layer_1,Layer_2,Layer_3,Layer_4,Layer_5,Layer_6)
    Layer_List <-
        lapply(Layer_List, add.missing.nodes,Number_of_Layers,Pool_of_Nodes)
    
    # We set the attributes of the layer
    counter <- 0 
    Layer_List <- lapply(Layer_List, function(x) { 
        counter <<- counter + 1; 
        set_edge_attr(x,"type",E(x), value = Layers_Name[counter])
    })
    
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
    for (i in seq_len(x$Number_of_Layers)){
        cat("\n")
        print(x[[i]])
    }
}

## Roxy Documentation comments
#' Create multiplex heterogeneous graphs from individual networks
#'
#' \code{create.multiplexHet} is a function to create a multiplex
#' and heterogeneous network (\code{MultiplexHet} object). It combines a
#' multiplex network composed from 1 (monoplex case) up to 6 layers with another
#' single network whose nodes are of different nature. See more details below.
#'
#' @usage create.multiplexHet(...)
#'
#' @details A multiplex network is a collection of layers (monoplex networks)
#' sharing the same nodes, but in which the edges represent relationships of
#' different nature. A heterogeneous network is composed of two single networks
#' where the nodes are of different nature. These nodes of different nature
#' are linked through bipartite interactions.
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
#' nodes Finally, it contains a expanded bipartite adjacency matrix
#' describing the relations of the nodes in every layer of the multiplex network
#' with the nodes of the second network.
#'
#' @seealso \code{\link{create.multiplex},\link{isMultiplexHet}}
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
    message("checking input arguments...")
    if (!isMultiplex(Multiplex_object)) {
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
                        frame are not present on the multiplex network")
                } else {
                    if (!all(names_2 %in% Het_nodes)){
                        stop("Some of the nodes in the second column of the data
                            frame are not present on the second network")
                    }
                }
            }
        }
    }
    
    if(missing(SecondNetworkName)){
        SecondNetworkName <-c("SecondNetwork")
    } else {
        if (!is.character(SecondNetworkName)) {
            stop("The names of the layers should be a vector of characters")}
        if (length(SecondNetworkName) != 1) {
            stop("The name of the Second Network name should be a vector of 
                length 1")}
        }
    
    ## Het graph
    Het_graph <-  simplify.layers(Het_graph)
    Het_graph<- 
        set_edge_attr(Het_graph,"type",E(Het_graph),value = SecondNetworkName)
    
    M <- vcount(Het_graph)
    
    ## Multiplex graph
    Nodes_Multiplex <- Multiplex_object$Pool_of_Nodes
    N <- Multiplex_object$Number_of_Nodes
    L <- Multiplex_object$Number_of_Layers
    
    message("Generating bipartite matrix...")
    Bipartite_Matrix <- 
        get.bipartite.graph(Nodes_Multiplex,Het_nodes,Nodes_relations,N, M)
    
    message("Expanding bipartite matrix to fit the multiplex network...")
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
    for (i in seq_len(x$Number_of_Layers)){
        cat("\n")
        print(x[[i]])
    }
    cat("\nNumber of Nodes of the second network:\n")
    print(x$Number_of_Nodes_Second_Network)
    cat("\nSecond Network\n")
    print(x$Second_Network)
}
