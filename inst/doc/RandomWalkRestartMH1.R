### R code from vignette source 'RandomWalkRestartMH1.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: style
###################################################
BiocStyle::latex()


###################################################
### code chunk number 2: installation (eval = FALSE)
###################################################
## source("http://bioconductor.org/biocLite.R")
## biocLite("RandomWalkRestartMH")


###################################################
### code chunk number 3: RandomWalkRestartMH1.Rnw:144-152
###################################################
library(RandomWalkRestartMH)
library(igraph)
data(PPI_Network) # We load the PPI_Network

## We create a Multiplex object of 1 layer(Monoplex) and we display how it
## looks like
PPI_MultiplexObject <- create.multiplex(PPI_Network,Layers_Name=c("PPI"))
PPI_MultiplexObject


###################################################
### code chunk number 4: RandomWalkRestartMH1.Rnw:158-160
###################################################
AdjMatrix_PPI <- compute.adjacency.matrix(PPI_MultiplexObject)
AdjMatrixNorm_PPI <- normalize.multiplex.adjacency(AdjMatrix_PPI)


###################################################
### code chunk number 5: RandomWalkRestartMH1.Rnw:167-173
###################################################
SeedGene <- c("PIK3R1")
## We launch the algorithm with the default parameters (See details on manual)
RWR_PPI_Results <- Random.Walk.Restart.Multiplex(AdjMatrixNorm_PPI,
                        PPI_MultiplexObject,SeedGene)
# We display the results
RWR_PPI_Results


###################################################
### code chunk number 6: RandomWalkRestartMH1.Rnw:181-185
###################################################
## In this case we select to induce a network with the Top 15 genes.
TopResults_PPI <-
    create.multiplexNetwork.topResults(RWR_PPI_Results,PPI_MultiplexObject,
                                       k=15)


###################################################
### code chunk number 7: fig1
###################################################
par(mar=c(0.1,0.1,0.1,0.1))
plot(TopResults_PPI, vertex.label.color="black",vertex.frame.color="#ffffff",
    vertex.size= 20, edge.curved=.2,
    vertex.color = ifelse(igraph::V(TopResults_PPI)$name == "PIK3R1","yellow",
    "#00CCFF"), edge.color="blue",edge.width=0.8)


###################################################
### code chunk number 8: RandomWalkRestartMH1.Rnw:221-239
###################################################
data(Disease_Network) # We load our disease Network

## We load a data frame containing the gene-disease associations.
## See ?create.multiplexHet for details about its format
data(GeneDiseaseRelations)

## We keep gene-diseases associations where genes are present in the PPI
## network
GeneDiseaseRelations_PPI <-
    GeneDiseaseRelations[which(GeneDiseaseRelations$hgnc_symbol %in%
    PPI_MultiplexObject$Pool_of_Nodes),]

## We create the MultiplexHet object.
PPI_Disease_Net <- create.multiplexHet(PPI_MultiplexObject,
    Disease_Network, GeneDiseaseRelations_PPI, c("Disease"))

## The results look like that
PPI_Disease_Net


###################################################
### code chunk number 9: RandomWalkRestartMH1.Rnw:246-247
###################################################
PPIHetTranMatrix <- compute.transition.matrix(PPI_Disease_Net)


###################################################
### code chunk number 10: RandomWalkRestartMH1.Rnw:254-263
###################################################
SeedDisease <- c("269880")

## We launch the algorithm with the default parameters (See details on manual)
RWRH_PPI_Disease_Results <-
    Random.Walk.Restart.MultiplexHet(PPIHetTranMatrix,
    PPI_Disease_Net,SeedGene,SeedDisease)

# We display the results
RWRH_PPI_Disease_Results


###################################################
### code chunk number 11: RandomWalkRestartMH1.Rnw:270-275
###################################################
## In this case we select to induce a network with the Top 10 genes
## and the Top 10 diseases.
TopResults_PPI_Disease <-
    create.multiplexHetNetwork.topResults(RWRH_PPI_Disease_Results,
    PPI_Disease_Net, GeneDiseaseRelations_PPI, k=10)


###################################################
### code chunk number 12: fig2
###################################################
par(mar=c(0.1,0.1,0.1,0.1))
plot(TopResults_PPI_Disease, vertex.label.color="black",
    vertex.frame.color="#ffffff",
    vertex.size= 20, edge.curved=.2,
    vertex.color = ifelse(V(TopResults_PPI_Disease)$name == "PIK3R1"
    | V(TopResults_PPI_Disease)$name == "269880","yellow",
    ifelse(V(TopResults_PPI_Disease)$name %in% PPI_Disease_Net$Pool_of_Nodes,
    "#00CCFF","Grey75")),
    edge.color=ifelse(E(TopResults_PPI_Disease)$type == "PPI","blue",
        ifelse(E(TopResults_PPI_Disease)$type == "Disease","black","grey50")),
    edge.width=0.8,
    edge.lty=ifelse(E(TopResults_PPI_Disease)$type == "bipartiteRelations",
        2,1),
    vertex.shape= ifelse(V(TopResults_PPI_Disease)$name %in%
        PPI_Disease_Net$Pool_of_Nodes,"circle","rectangle"))


###################################################
### code chunk number 13: RandomWalkRestartMH1.Rnw:326-332
###################################################
data(Pathway_Network) # We load the Pathway Network

## We create a 2-layers Multiplex object
PPI_PATH_Multiplex <- create.multiplex(PPI_Network,Pathway_Network,
                        Layers_Name=c("PPI","PATH"))
PPI_PATH_Multiplex


###################################################
### code chunk number 14: RandomWalkRestartMH1.Rnw:338-340
###################################################
AdjMatrix_PPI_PATH <- compute.adjacency.matrix(PPI_PATH_Multiplex)
AdjMatrixNorm_PPI_PATH <- normalize.multiplex.adjacency(AdjMatrix_PPI_PATH)


###################################################
### code chunk number 15: RandomWalkRestartMH1.Rnw:346-351
###################################################
## We launch the algorithm with the default parameters (See details on manual)
RWR_PPI_PATH_Results <- Random.Walk.Restart.Multiplex(AdjMatrixNorm_PPI_PATH,
                        PPI_PATH_Multiplex,SeedGene)
# We display the results
RWR_PPI_PATH_Results


###################################################
### code chunk number 16: RandomWalkRestartMH1.Rnw:357-361
###################################################
## In this case we select to induce a multiplex network with the Top 15 genes.
TopResults_PPI_PATH <-
    create.multiplexNetwork.topResults(RWR_PPI_PATH_Results,
                                        PPI_PATH_Multiplex, k=15)


###################################################
### code chunk number 17: fig3
###################################################
par(mar=c(0.1,0.1,0.1,0.1))
plot(TopResults_PPI_PATH, vertex.label.color="black",
    vertex.frame.color="#ffffff", vertex.size= 20,
    edge.curved= ifelse(E(TopResults_PPI_PATH)$type == "PPI",
                    0.4,0),
    vertex.color = ifelse(igraph::V(TopResults_PPI_PATH)$name == "PIK3R1",
                    "yellow","#00CCFF"),edge.width=0.8,
    edge.color=ifelse(E(TopResults_PPI_PATH)$type == "PPI",
                      "blue","red"))


###################################################
### code chunk number 18: RandomWalkRestartMH1.Rnw:402-414
###################################################
## We keep gene-diseases associations where genes are present in the PPI
## or in the pathway network
GeneDiseaseRelations_PPI_PATH <-
    GeneDiseaseRelations[which(GeneDiseaseRelations$hgnc_symbol %in%
    PPI_PATH_Multiplex$Pool_of_Nodes),]

## We create the MultiplexHet object.
PPI_PATH_Disease_Net <- create.multiplexHet(PPI_PATH_Multiplex,
    Disease_Network, GeneDiseaseRelations_PPI_PATH, c("Disease"))

## The results look like that
PPI_PATH_Disease_Net


###################################################
### code chunk number 19: RandomWalkRestartMH1.Rnw:421-422
###################################################
PPI_PATH_HetTranMatrix <- compute.transition.matrix(PPI_PATH_Disease_Net)


###################################################
### code chunk number 20: RandomWalkRestartMH1.Rnw:428-435
###################################################
## We launch the algorithm with the default parameters (See details on manual)
RWRH_PPI_PATH_Disease_Results <-
    Random.Walk.Restart.MultiplexHet(PPI_PATH_HetTranMatrix,
    PPI_PATH_Disease_Net,SeedGene,SeedDisease)

# We display the results
RWRH_PPI_PATH_Disease_Results


###################################################
### code chunk number 21: RandomWalkRestartMH1.Rnw:442-447
###################################################
## In this case we select to induce a network with the Top 10 genes.
## and the Top 10 diseases.
TopResults_PPI_PATH_Disease <-
    create.multiplexHetNetwork.topResults(RWRH_PPI_PATH_Disease_Results,
    PPI_PATH_Disease_Net, GeneDiseaseRelations_PPI_PATH, k=10)


###################################################
### code chunk number 22: fig4
###################################################
par(mar=c(0.1,0.1,0.1,0.1))
plot(TopResults_PPI_PATH_Disease, vertex.label.color="black",
    vertex.frame.color="#ffffff",
    vertex.size= 20,
    edge.curved=ifelse(E(TopResults_PPI_PATH_Disease)$type == "PATH",
                    0,0.3),
    vertex.color = ifelse(V(TopResults_PPI_PATH_Disease)$name == "PIK3R1"
    | V(TopResults_PPI_Disease)$name == "269880","yellow",
    ifelse(V(TopResults_PPI_PATH_Disease)$name %in%
               PPI_PATH_Disease_Net$Pool_of_Nodes,
    "#00CCFF","Grey75")),
    edge.color=ifelse(E(TopResults_PPI_PATH_Disease)$type == "PPI","blue",
    ifelse(E(TopResults_PPI_PATH_Disease)$type == "PATH","red",
    ifelse(E(TopResults_PPI_PATH_Disease)$type == "Disease","black","grey50"))),
    edge.width=0.8,
    edge.lty=ifelse(E(TopResults_PPI_PATH_Disease)$type ==
        "bipartiteRelations", 2,1),
    vertex.shape= ifelse(V(TopResults_PPI_PATH_Disease)$name %in%
        PPI_PATH_Disease_Net$Pool_of_Nodes,"circle","rectangle"))


###################################################
### code chunk number 23: sessionInfo
###################################################
toLatex(sessionInfo())


