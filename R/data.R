#' A protein-protein physical interaction network (PPI network)
#'
#' An igraph object containing a protein-protein physical interaction network
#' (PPI network). The network is obtained as described in the article 
#' cited in the source section. However, it was reduced in such a way
#' that only contains genes/proteins expressed in the adipose tissue. 
#'
#' @format An igraph object containing 18062 binary interactions between 4317
#' proteins
#'
#' @docType data
#' @keywords datasets
#' @name PPI_Network
#' @usage data(PPI_Network)
#' @source Valdeolivas, A., Tichit, L., Navarro, C., Perrin, S., Odelin, G., 
#' Levy, N., … Baudot, A. (2017). Random Walk With Restart On Multiplex And 
#' Heterogeneous Biological Networks. bioRxiv, 1–31. 
#' https://doi.org/10.1101/134734 
#' \url{https://www.biorxiv.org/content/early/2017/08/30/134734}
NULL

#' A pathway network (Pathway network)
#'
#' An igraph object containing a Pathway network. 
#' The network is obtained as described in the article 
#' cited in the source section. However, it was reduced in such a way
#' that only contains genes/proteins expressed in the adipose tissue.

#' @format An igraph object containing 62602 binary interactions between 3533
#' proteins
#'
#' @docType data
#' @keywords datasets
#' @name Pathway_Network
#' @usage data(Pathway_Network)
#' @source Valdeolivas, A., Tichit, L., Navarro, C., Perrin, S., Odelin, G.,
#' Levy, N., … Baudot, A. (2017). Random Walk With Restart On Multiplex And
#' Heterogeneous Biological Networks. bioRxiv, 1–31.
#' https://doi.org/10.1101/134734
#' \url{https://www.biorxiv.org/content/early/2017/08/30/134734}
NULL

#' A disease-disease similarity network.
#'
#' An igraph object containing a disease-disease similarity network. 
#' The network is obtained as described in the article
#' cited in the source section.
#'
#' @format An igraph object containing 28246 binary relationships between 6947
#' diseases.
#'
#' @docType data
#' @keywords datasets
#' @name Disease_Network
#' @usage data(Disease_Network)
#' @source Valdeolivas, A., Tichit, L., Navarro, C., Perrin, S., Odelin, G.,
#' Levy, N., … Baudot, A. (2017). Random Walk With Restart On Multiplex And
#' Heterogeneous Biological Networks. bioRxiv, 1–31.
#' https://doi.org/10.1101/134734
#' \url{https://www.biorxiv.org/content/early/2017/08/30/134734}
NULL

#' Diseases and their causative genes
#'
#' A dataset containing some diseases and their causative genes.
#' The dataset is obtained as described in the article
#' cited in the source section.
#'
#' @format A data frame with 4496 rows and 2 variables:
#' \describe{
#'     \item{hgnc_symbol}{Gene name, in HGNC format}
#'     \item{mim_morbid}{Disease id, in mim code}
#' }
#'
#' @docType data
#' @keywords datasets
#' @name GeneDiseaseRelations
#' @usage data(GeneDiseaseRelations)
#' @source Valdeolivas, A., Tichit, L., Navarro, C., Perrin, S., Odelin, G.,
#' Levy, N., … Baudot, A. (2017). Random Walk With Restart On Multiplex And
#' Heterogeneous Biological Networks. bioRxiv, 1–31.
#' https://doi.org/10.1101/134734
#' \url{https://www.biorxiv.org/content/early/2017/08/30/134734}
NULL



