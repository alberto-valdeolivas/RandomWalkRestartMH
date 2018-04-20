## R Code for the Random Walk with Restart Package (RandomWalkRestartMH).

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Functions for internal use
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

## Geometric mean: R computation.

geometric.mean <- function(Scores, L, N) {
    
    FinalScore <- numeric(length = N)
    
    for (i in seq_len(N)){
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
    Layer <- igraph::simplify(Layer,remove.multiple = TRUE,remove.loops = TRUE)
    
    return(Layer)
}


## Add missing nodes in some of the layers.
add.missing.nodes <- function (Layers,Nr_Layers,NodeNames) {
    
    add_vertices(Layers,
        length(NodeNames[which(!NodeNames %in% V(Layers)$name)]),
             name=NodeNames[which(!NodeNames %in%  V(Layers)$name)])
}

## Scores for the seeds of the multiplex network.
get.seed.scoresMultiplex <- function(Seeds,Number_Layers,tau) {
    
    Nr_Seeds <- length(Seeds)
    
    Seeds_Seeds_Scores <- rep(tau/Nr_Seeds,Nr_Seeds)
    Seed_Seeds_Layer_Labeled <- 
        paste0(rep(Seeds,Number_Layers),sep="_",rep(seq(Number_Layers), 
            length.out = Nr_Seeds*Number_Layers,each=Nr_Seeds))
    
    Seeds_Score <- data.frame(Seeds_ID = Seed_Seeds_Layer_Labeled,
        Score = Seeds_Seeds_Scores, stringsAsFactors = FALSE)
    
    return(Seeds_Score)
}

## Bipartite graph construction.
get.bipartite.graph <- 
    function(Names_Mul, Names_Het, Nodes_relation, Number_Nodes_1,
        Number_Nodes_2){
        
    Bipartite_matrix <- Matrix(data=0, nrow=Number_Nodes_1, ncol=Number_Nodes_2)
    Names_Mul_order <- sort(Names_Mul)
    Names_Het_order <- sort(Names_Het)
    rownames(Bipartite_matrix) <- Names_Mul_order
    colnames(Bipartite_matrix) <- Names_Het_order
        
    for (i in seq_len(Number_Nodes_1)){
        current_node1 <- Names_Mul_order[i]
        current_node2 <-
            Nodes_relation[which(Nodes_relation[,1] == current_node1),2]
            
            
        if (length(current_node2) > 0){
            for (j in seq_len(length(current_node2))){
                if (!is.na(current_node2[j])){
                    ## We need to identify the position on the matrix.
                    index <- which(colnames(Bipartite_matrix) %in% 
                        current_node2[j]) 
                                
                    Bipartite_matrix[i,index] <- 1
                }
            }
        }
    }
    return(Bipartite_matrix)
}

## Fitting the bipartite graph to the multiplex network.
expand.bipartite.graph <- 
    function(Number_Nodes_1,Number_Layers,Number_Nodes_2,Bipartite_matrix){
        
    Supra_Bipartite_Matrix <- 
        do.call(rbind, replicate(Number_Layers,Bipartite_matrix,simplify=FALSE))
        
    rownames(Supra_Bipartite_Matrix) <- 
        paste0(rownames(Bipartite_matrix), sep="_",rep(seq(Number_Layers),
            each=Number_Nodes_1))
        
    colnames(Supra_Bipartite_Matrix) <- colnames(Bipartite_matrix)
        
    return(Supra_Bipartite_Matrix)
}


## Transitions for the computation of the Multiplex Heterogeneous transition
## Matrix.

get.transition.multiplex.secondNet <- 
    function(Number_Nodes_Multiplex, Number_Layers,Number_Nodes_secondNet,
        SupraBipartiteMatrix,lambda){
        
    Transition_Multiplex_SecondNet <- 
        Matrix(0, nrow=Number_Nodes_Multiplex*Number_Layers,
            ncol=Number_Nodes_secondNet,sparse = TRUE)
    colnames(Transition_Multiplex_SecondNet) <- colnames(SupraBipartiteMatrix)
    rownames(Transition_Multiplex_SecondNet) <- rownames(SupraBipartiteMatrix)
        
    Col_Sum_Bipartite <- 
        Matrix::colSums (SupraBipartiteMatrix, na.rm = FALSE, dims = 1,
            sparseResult = FALSE)
        
    m <- lambda * t(t(SupraBipartiteMatrix) / Col_Sum_Bipartite)
    idx <- Col_Sum_Bipartite != 0
    Transition_Multiplex_SecondNet[,idx] = m[,idx]
        
    return(Transition_Multiplex_SecondNet)
}

## Transitions for the computation of the Multiplex Heterogeneous transition
## Matrix.
get.transition.secondNet.multiplex <- 
    function(Number_Nodes_Multiplex, Number_Layers,Number_Nodes_secondNet,
        SupraBipartiteMatrix,lambda){
        
    Transition_SecondNet_Multiplex <- 
        Matrix(0,nrow=Number_Nodes_secondNet, 
            ncol=Number_Nodes_Multiplex*Number_Layers,sparse = TRUE)
        
    colnames(Transition_SecondNet_Multiplex) <- rownames(SupraBipartiteMatrix)
    rownames(Transition_SecondNet_Multiplex) <- colnames(SupraBipartiteMatrix)
        
    Row_Sum_Bipartite <- 
        Matrix::rowSums (SupraBipartiteMatrix, na.rm = FALSE, dims = 1,
            sparseResult = FALSE)
        
    m <- lambda * t((SupraBipartiteMatrix) / Row_Sum_Bipartite)
    idx <- Row_Sum_Bipartite != 0
    Transition_SecondNet_Multiplex[,idx] = m[,idx]
        
    return(Transition_SecondNet_Multiplex)
}

## Transitions for the computation of the Multiplex Heterogeneous transition
## Matrix.
get.transition.multiplex <- 
    function(Number_Nodes_Multiplex,Number_Layers, lambda,SupraAdjacencyMatrix,
        SupraBipartiteMatrix){
        
    Transition_Multiplex_Network <- 
        Matrix(0, nrow=Number_Nodes_Multiplex*Number_Layers,
            ncol=Number_Nodes_Multiplex*Number_Layers,sparse = TRUE)
        
    rownames(Transition_Multiplex_Network) <- rownames(SupraAdjacencyMatrix)
    colnames(Transition_Multiplex_Network) <- colnames(SupraAdjacencyMatrix)
        
    Col_Sum_Multiplex <- 
        Matrix::colSums(SupraAdjacencyMatrix,na.rm=FALSE, dims=1, 
            sparseResult=FALSE)
    Row_Sum_Bipartite <- 
        Matrix::rowSums (SupraBipartiteMatrix, na.rm = FALSE, dims = 1,
            sparseResult = FALSE)
        
    idx <- Row_Sum_Bipartite != 0
    Transition_Multiplex_Network[,idx] <- 
        ((1-lambda)*t(t(SupraAdjacencyMatrix[,idx])/Col_Sum_Multiplex[idx]))
        
    Transition_Multiplex_Network[,!idx] <-
        t(t(SupraAdjacencyMatrix[,!idx]) / Col_Sum_Multiplex[!idx])
        
    return(Transition_Multiplex_Network)
}


## Transitions for the computation of the Multiplex Heterogeneous transition
## Matrix.
get.transition.secondNet <- 
    function(Number_Nodes_secondNet,lambda, AdjMatrix,SupraBipartiteMatrix){
        
    Transition_Second_Network <- 
        Matrix(0,nrow=Number_Nodes_secondNet, ncol=Number_Nodes_secondNet,
            sparse = TRUE)
        
    rownames(Transition_Second_Network) <- rownames(AdjMatrix)
    colnames(Transition_Second_Network) <- colnames(AdjMatrix)
        
    Col_Sum_SecondNet <- 
        Matrix::colSums (AdjMatrix,na.rm=FALSE, dims=1,sparseResult=FALSE)
    Col_Sum_Bipartite <- 
        Matrix::colSums (SupraBipartiteMatrix, na.rm = FALSE, dims = 1,
            sparseResult = FALSE)
        
    idx <- Col_Sum_Bipartite != 0
    Transition_Second_Network[,idx] <- 
        ((1-lambda)*t(t(AdjMatrix[,idx]) / Col_Sum_SecondNet[idx]))
        
    Transition_Second_Network[,!idx] <-
        t(t(AdjMatrix[,!idx]) / Col_Sum_SecondNet[!idx])
        
    return(Transition_Second_Network)
}


## Compute the scores for the seeds on the RWRMH approach
get.seed.scores.multHet <- 
    function(Multiplex_Seed_Nodes,SecondNet_Seed_Nodes,eta,L,tau) {
        
    n <- length(Multiplex_Seed_Nodes)
    m <- length(SecondNet_Seed_Nodes)
        
    if ((n != 0 && m!= 0)){
            
        Seed_Multiplex_Layer_Labeled <- paste0(rep(Multiplex_Seed_Nodes,L),
            sep="_",rep(seq(L), length.out = n*L,each=n))
            
        Seeds_Multiplex_Scores <- rep(((1-eta) * tau)/n,n)
        SecondNet_Seeds_Score <- eta/m
            
    } else {
        eta <- 1
        if (n == 0){
            Seed_Multiplex_Layer_Labeled <- character()
            Seeds_Multiplex_Scores <- numeric()
            SecondNet_Seeds_Score <- eta/m
        } else {
            Seed_Multiplex_Layer_Labeled <- 
                paste0(rep(Multiplex_Seed_Nodes,L), sep="_",rep(seq(L),
                    length.out = n*L,each=n))
                
            Seeds_Multiplex_Scores <- rep(tau/n,n)
            SecondNet_Seeds_Score <- numeric()
        }
    }
        
    ## We prepare a data frame with the seeds.
    Seeds_Score <- data.frame(Seeds_ID = 
        c(Seed_Multiplex_Layer_Labeled,SecondNet_Seed_Nodes),
            Score = c(Seeds_Multiplex_Scores, rep(SecondNet_Seeds_Score,m)),
                stringsAsFactors = FALSE)
        
    return(Seeds_Score)
}


## Ranking for the multiplex nodes
rank.nodes.multiplex <- function(N, L, Results,Seeds){
    
    ## We sort the score to obtain the ranking of multiplex nodes and
    ## seconde network nodes.
    NodeNames <- character(length = N)
    Score <- numeric(length = N)
    
    nodes_multiplex_rank <- data.frame(NodeNames = NodeNames, Score = Score)
    nodes_multiplex_rank$NodeNames <- 
        gsub("_1", "",row.names(Results)[seq_len(N)])
    
    ## We calculate the Geometric Mean among the nodes in the
    ## different layers.
    nodes_multiplex_rank$Score <- geometric.mean(as.vector(Results[,1]),L,N)
    
    nodes_multiplex_sort <- 
        nodes_multiplex_rank[with(nodes_multiplex_rank, 
                                  order(-Score, NodeNames)), ]
    
    ## We remove the seed nodes from the Ranking
    nodes_multiplex_sort_NoSeeds <-
        nodes_multiplex_sort[which(!nodes_multiplex_sort$NodeNames %in% Seeds),]
    
    return(nodes_multiplex_sort_NoSeeds)
}


## Ranking for the Second Network nodes
rank.nodes.secondNet <- function(N,L,M,Results,Seeds){
    
    SecondNet_node <- character(length = M)
    Score <- character(length = M)
    SecondNet_rank <- data.frame(SecondNet_node = SecondNet_node, Score = Score)
    SecondNet_rank$SecondNet_node <- row.names(Results)[((N*L)+1):nrow(Results)]
    SecondNet_rank$Score <- Results[((N*L)+1):nrow(Results),1]
    
    SecondNet_rank_sort <- 
        SecondNet_rank[with(SecondNet_rank,order(-Score, SecondNet_node)), ]
    SecondNet_rank_sort_NoSeeds <-
        SecondNet_rank_sort[which(!SecondNet_rank_sort$SecondNet_node 
                                  %in% Seeds),]
    
    return(SecondNet_rank_sort_NoSeeds)
}
