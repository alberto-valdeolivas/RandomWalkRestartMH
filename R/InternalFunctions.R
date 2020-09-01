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
    
    ## Unweighted or Weigthed Graphs
    if (is_weighted(Layer)){
        b <- 1
        weigths_layer <- E(Layer)$weight
        if (min(weigths_layer) != max(weigths_layer)){
            a <- min(weigths_layer)/max(weigths_layer)
            range01 <- (b-a)*(weigths_layer-min(weigths_layer))/
                (max(weigths_layer)-min(weigths_layer)) + a
            E(Layer)$weight <- range01
        } else {
            E(Layer)$weight <- rep(1, length(weigths_layer))
        }
    } else {
        E(Layer)$weight <- rep(1, ecount(Layer))
    }
    
    ## Simple Graphs
    Layer <- 
        igraph::simplify(Layer,remove.multiple = TRUE,remove.loops = TRUE, 
            edge.attr.comb=mean)
    
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
    function(Names_Mul_1, Names_Mul_2, Nodes_relation, Number_Nodes_1,
        Number_Nodes_2){
        
    Bipartite_matrix <- Matrix(data=0, nrow=Number_Nodes_1, ncol=Number_Nodes_2)
    Names_Mul1_order <- sort(Names_Mul_1)
    Names_Mul2_order <- sort(Names_Mul_2)
    rownames(Bipartite_matrix) <- Names_Mul1_order
    colnames(Bipartite_matrix) <- Names_Mul2_order
        
    for (i in seq_len(Number_Nodes_1)){
        current_node1 <- Names_Mul_1[i]
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
    function(Number_Nodes_1, Number_Layers_1, Number_Nodes_2,
        Number_Layers_2, Bipartite_matrix){
        
        Supra_Bipartite_Matrix <- 
            do.call(rbind, replicate(Number_Layers_1,Bipartite_matrix,
                simplify=FALSE))
        
        rownames(Supra_Bipartite_Matrix) <- 
            paste0(rownames(Bipartite_matrix), sep="_",rep(seq(Number_Layers_1),
                each=Number_Nodes_1))
        
        
        Supra_Bipartite_Matrix <- 
            do.call(cbind, replicate(Number_Layers_2,Supra_Bipartite_Matrix,
                simplify=FALSE))
        
        colnames(Supra_Bipartite_Matrix) <- 
            paste0(colnames(Bipartite_matrix), sep="_",rep(seq(Number_Layers_2),
                each=Number_Nodes_2))
        
        return(Supra_Bipartite_Matrix)
}


## Transitions for the computation of the Multiplex Heterogeneous transition
## Matrix.

get.transition.multiplex1.multiplex2 <- function(Number_Nodes_Multiplex1, 
    Number_Layers1,Number_Nodes_Multiplex2, Number_Layers2, 
    SupraBipartiteMatrix,lambda){
        
    TransitionMat_Multiplex1_Multiplex2 <- 
        Matrix(0, nrow=Number_Nodes_Multiplex1*Number_Layers1,
            ncol=Number_Nodes_Multiplex2*Number_Layers2,sparse = TRUE)
        
    colnames(TransitionMat_Multiplex1_Multiplex2) <- 
        colnames(SupraBipartiteMatrix)
    rownames(TransitionMat_Multiplex1_Multiplex2) <- 
        rownames(SupraBipartiteMatrix)
        
    Col_Sum_Bipartite <- Matrix::colSums (SupraBipartiteMatrix, na.rm = FALSE, 
        dims = 1, sparseResult = FALSE)
        
    m <- lambda * t(t(SupraBipartiteMatrix) / Col_Sum_Bipartite)
    idx <- Col_Sum_Bipartite != 0
    if (length(idx) > 0){
        TransitionMat_Multiplex1_Multiplex2[,idx] = m[,idx]    
    }
    
    return(TransitionMat_Multiplex1_Multiplex2)
}

## Transitions for the computation of the Multiplex Heterogeneous transition
## Matrix.
get.transition.multiplex2.multiplex1 <- function(Number_Nodes_Multiplex1, 
    Number_Layers1,Number_Nodes_Multiplex2, Number_Layers2,SupraBipartiteMatrix,
    lambda){
        
    TransitionMat_Multiplex2_Multiplex1 <- 
        Matrix(0,nrow=Number_Nodes_Multiplex2*Number_Layers2, 
            ncol=Number_Nodes_Multiplex1*Number_Layers1,sparse = TRUE)
        
    colnames(TransitionMat_Multiplex2_Multiplex1) <- 
        rownames(SupraBipartiteMatrix)
    rownames(TransitionMat_Multiplex2_Multiplex1) <- 
        colnames(SupraBipartiteMatrix)
        
    Row_Sum_Bipartite <- 
        Matrix::rowSums (SupraBipartiteMatrix, na.rm = FALSE, dims = 1,
            sparseResult = FALSE)
        
    m <- lambda * t((SupraBipartiteMatrix) / Row_Sum_Bipartite)
    idx <- Row_Sum_Bipartite != 0
    if (length(idx) > 0) {
        TransitionMat_Multiplex2_Multiplex1[,idx] = m[,idx]
    }   
    
    return(TransitionMat_Multiplex2_Multiplex1)
}

## Transitions for the computation of the Multiplex Heterogeneous transition
## Matrix.
get.transition.multiplex <- 
    function(Number_Nodes,Number_Layers, lambda,SupraAdjacencyMatrix,
        SupraBipartiteMatrix) {
        
    Transition_Multiplex_Network <- Matrix(0, nrow=Number_Nodes*Number_Layers,
        ncol=Number_Nodes*Number_Layers,sparse = TRUE)
        
    rownames(Transition_Multiplex_Network) <- rownames(SupraAdjacencyMatrix)
    colnames(Transition_Multiplex_Network) <- colnames(SupraAdjacencyMatrix)
        
    Col_Sum_Multiplex <- Matrix::colSums(SupraAdjacencyMatrix,na.rm=FALSE, 
        dims=1, sparseResult=FALSE)
    Row_Sum_Bipartite <- Matrix::rowSums (SupraBipartiteMatrix, na.rm = FALSE, 
        dims = 1, sparseResult = FALSE)
        
    idx <- Row_Sum_Bipartite != 0
    if (length(idx) > 0){
        Transition_Multiplex_Network[,idx] <- 
            ((1-lambda)*t(t(SupraAdjacencyMatrix[,idx])/Col_Sum_Multiplex[idx]))
        
        Transition_Multiplex_Network[,!idx] <-
            t(t(SupraAdjacencyMatrix[,!idx]) / Col_Sum_Multiplex[!idx])
    }    
        
    return(Transition_Multiplex_Network)
}

## Transitions for the computation of the Multiplex Heterogeneous transition
## Matrix.
#get.transition.secondNet <- 
#    function(Number_Nodes_secondNet,lambda, AdjMatrix,SupraBipartiteMatrix){
#        
#    Transition_Second_Network <- 
#        Matrix(0,nrow=Number_Nodes_secondNet, ncol=Number_Nodes_secondNet,
#            sparse = TRUE)
#        
#    rownames(Transition_Second_Network) <- rownames(AdjMatrix)
#    colnames(Transition_Second_Network) <- colnames(AdjMatrix)
#        
#    Col_Sum_SecondNet <- 
#        Matrix::colSums (AdjMatrix,na.rm=FALSE, dims=1,sparseResult=FALSE)
#    Col_Sum_Bipartite <- 
#        Matrix::colSums (SupraBipartiteMatrix, na.rm = FALSE, dims = 1,
#            sparseResult = FALSE)
#        
#    idx <- Col_Sum_Bipartite != 0
#    Transition_Second_Network[,idx] <- 
#        ((1-lambda)*t(t(AdjMatrix[,idx]) / Col_Sum_SecondNet[idx]))
#        
#    Transition_Second_Network[,!idx] <-
#        t(t(AdjMatrix[,!idx]) / Col_Sum_SecondNet[!idx])
#        
#    return(Transition_Second_Network)
#}

#### Computing the scores for the different seed nodes

get.seed.scores.multHet <- 
    function(Multiplex1_Seed_Nodes,Multiplex2_Seed_Nodes,eta,L1,L2,tau1,tau2) {
        
    n <- length(Multiplex1_Seed_Nodes)
    m <- length(Multiplex2_Seed_Nodes)
        
    if ((n != 0 && m!= 0)){
            
        Seed_Multiplex1_Layer_Labeled <- paste0(rep(Multiplex1_Seed_Nodes,L1),
            sep="_",rep(seq(L1), length.out = n*L1,each=n))
            
        Seed_Multiplex2_Layer_Labeled <- paste0(rep(Multiplex2_Seed_Nodes,L2),
            sep="_",rep(seq(L2), length.out = m*L2,each=m))
            
        Seeds_Multiplex1_Scores <- rep(((1-eta) * tau1)/n,n)
        Seeds_Multiplex2_Scores <- rep((eta*tau2)/m,m)
            
    } else {
        eta <- 1
        if (n == 0){
            Seed_Multiplex1_Layer_Labeled <- character()
            Seeds_Multiplex1_Scores <- numeric()
            Seed_Multiplex2_Layer_Labeled <- 
                paste0(rep(Multiplex2_Seed_Nodes,L2), sep="_",rep(seq(L2),
                    length.out = m*L2,each=m))
            Seeds_Multiplex2_Scores <- rep(tau2/m,m)
            } else {
                Seed_Multiplex1_Layer_Labeled <- 
                    paste0(rep(Multiplex1_Seed_Nodes,L1), sep="_",rep(seq(L1),
                        length.out = n*L1,each=n))
                Seeds_Multiplex1_Scores <- rep(tau1/n,n)
                Seed_Multiplex2_Layer_Labeled <- character()
                Seeds_Multiplex2_Scores <- numeric()
            }
        }
        
    ## We prepare a data frame with the seeds.
    Seeds_Score <- data.frame(Seeds_ID = 
         c(Seed_Multiplex1_Layer_Labeled,Seed_Multiplex2_Layer_Labeled),
         Score = c(Seeds_Multiplex1_Scores, Seeds_Multiplex2_Scores),
         stringsAsFactors = FALSE)
        
    return(Seeds_Score)
}

geometric.mean <- function(Scores, L, N) {
    
    FinalScore <- numeric(length = N)
    
    for (i in seq_len(N)){
        FinalScore[i] <- prod(Scores[seq(from = i, to = N*L, by=N)])^(1/L)
    }
    
    return(FinalScore)
}

regular.mean <- function(Scores, L, N) {
    
    FinalScore <- numeric(length = N)
    
    for (i in seq_len(N)){
        FinalScore[i] <- mean(Scores[seq(from = i, to = N*L, by=N)])
    }
    
    return(FinalScore)
}

sumValues <- function(Scores, L, N) {
    
    FinalScore <- numeric(length = N)
    
    for (i in seq_len(N)){
        FinalScore[i] <- sum(Scores[seq(from = i, to = N*L, by=N)])
    }
    
    return(FinalScore)
}

## Ranking for the multiplex nodes
#rank.nodes.multiplex <- function(N, L, Results,Seeds){
#    
#    ## We sort the score to obtain the ranking of multiplex nodes and
#    ## seconde network nodes.
#    NodeNames <- character(length = N)
#    Score <- numeric(length = N)
#    
#    nodes_multiplex_rank <- data.frame(NodeNames = NodeNames, Score = Score)
#    nodes_multiplex_rank$NodeNames <- 
#        gsub("_1", "",row.names(Results)[seq_len(N)])
#    
#    ## We calculate the Geometric Mean among the nodes in the
#    ## different layers.
#    nodes_multiplex_rank$Score <- geometric.mean(as.vector(Results[,1]),L,N)
#    
#    nodes_multiplex_sort <- 
#        nodes_multiplex_rank[with(nodes_multiplex_rank, 
#                                  order(-Score, NodeNames)), ]
#    
#    ## We remove the seed nodes from the Ranking
#    nodes_multiplex_sort_NoSeeds <-
#        nodes_multiplex_sort[which(!nodes_multiplex_sort$NodeNames %in% Seeds),]
#    
#    return(nodes_multiplex_sort_NoSeeds)
#}


## Ranking for the Second Network nodes
#rank.nodes.secondNet <- function(N,L,M,Results,Seeds){
#    
#    SecondNet_node <- character(length = M)
#    Score <- character(length = M)
#    SecondNet_rank <- data.frame(SecondNet_node = SecondNet_node, Score = Score)
#    SecondNet_rank$SecondNet_node <- row.names(Results)[((N*L)+1):nrow(Results)]
#    SecondNet_rank$Score <- Results[((N*L)+1):nrow(Results),1]
#    
#    SecondNet_rank_sort <- 
#        SecondNet_rank[with(SecondNet_rank,order(-Score, SecondNet_node)), ]
#    SecondNet_rank_sort_NoSeeds <-
#        SecondNet_rank_sort[which(!SecondNet_rank_sort$SecondNet_node 
#                                  %in% Seeds),]
#    
#    return(SecondNet_rank_sort_NoSeeds)
#}
