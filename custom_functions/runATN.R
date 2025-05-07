runATN <- function(graph,
                   initialBiomasses = NULL,
                   regexForC4Nodes = "^C4-",
                   regexForC3Nodes = "^C3-",
                   customK=NULL,
                   C3IntraspecificComp = 1,
                   C4IntraspecificComp = 1,
                   carnivoreIntraspecificCompetition = 1,
                   herbivoreIntraspecificCompetition = 1
                   )
  {
  
  #'graph' specifies the igraph object which must contain certain vertex and edge attributes
  # most likely, this graph was produced by the custom makeGraphForATN function
  # argument regexForC3Nodes and regexForC4Nodes provides regular expressions which identify
  # producer nodes based on their names in the vertex attributes.
  # C3IntraspecificComp and C4IntraspecificComp specifiy the value of the diagonals
  # of the alpha matrix for C3 and C4 producer nodes, respectively. Interspecific competition is set to 0 by default
  # Lower values for intraspecific competition INCREASE effective carrying capacity for nodes of that type
  # return value is the solved matrix of biomasses, plus the original igraph object
  # for downstream plotting or analysis
  
  if(!"igraph" %in% class(graph)) stop("runATN requires an object of class igraph")
  
  require(igraph)
  require(ATNr)
  
  vertexAttributes <- as.data.frame(vertex_attr(graph))
  edgeAttributes <- as.data.frame(edge_attr(graph))
  
  if(!"BM" %in% names(vertexAttributes)) stop("can't find body masses (variable name 'BM') in vertex attributes")
  if(!"guild" %in% names(vertexAttributes)) stop("can't find guilds (variable name 'guild') in vertex attributes")
  if(!"linkProb" %in% names(edgeAttributes)) stop("can't find link probabilites (variable name 'linkProb') in vertex attributes")
  
  gAdjMat <- as_adjacency_matrix(graph, sparse = FALSE)
  
  isC4 <- grepl(regexForC4Nodes, vertexAttributes$name)
  isC3 <- grepl(regexForC3Nodes, vertexAttributes$name)
  
  mod <- create_model_Scaled(nb_s = dim(gAdjMat)[1], 
                             nb_b = sum(c(isC4, isC3)), 
                             BM   = vertexAttributes$BM, 
                             fw   = gAdjMat)
  mod <- initialise_default_Scaled(mod)
  
  #don't assume we know whether C4 or C3 producer nodes come first in the vertex attributes
  #check to see.  We do assume that all C4 and C3 producer nodes are specified in blocks, not interleaved
  if(which(isC4)[1]==1) diag(mod$alpha)  <- c(rep(C4IntraspecificComp, sum(isC4)), rep(C3IntraspecificComp, sum(isC3)))
  if(!which(isC4)[1]==1) diag(mod$alpha) <- c(rep(C3IntraspecificComp, sum(isC3)), rep(C4IntraspecificComp, sum(isC4)))
  
  if(!is.null(customK)) mod$K <- customK

  mod$c[which(vertex_attr(graph)$guild == "carnivore") - sum(c(isC4, isC3))] <- carnivoreIntraspecificCompetition
  mod$c[which(vertex_attr(graph)$guild %in% c("browser", "mixed feeder", "grazer")) - sum(c(isC4, isC3))] <- herbivoreIntraspecificCompetition
  
  
  
  #modify w matrix
  #convert the graph to an edgelist with two cols resource and consumer
  edgelist <- data.frame(as_edgelist(graph, names = TRUE))
  names(edgelist) <- c("resource", "consumer")
  #add in the link probabilities from the edge attributes
  edgeDF <- data.frame(edgelist, prob=edgeAttributes$linkProb)
  
  w <- matrix(nrow=mod$nb_s, ncol=mod$nb_s - mod$nb_b)
  
  consumerTaxaNames <- vertexAttributes$name[!isC4 & !isC3]

  ##loop through the indices of the resource taxa names (i) and the indices
  ## of the consumer taxa names (j) and fill out the w matrix with the link weight
  for(i in 1:length(vertexAttributes$name)) {
    for(j in 1:length(consumerTaxaNames)) {
      probVal <- edgeDF[edgeDF$resource == vertexAttributes$name[i] & edgeDF$consumer == consumerTaxaNames[j],"prob"]
      if(length(probVal) > 0) {
        w[i,j] <- probVal
      }
    }
  }
  #convert NAs to zero
  w[is.na(w)] <- 0

  #scale by the column sums making each consumer's total feeding weight = 1
  mod$w <- scale(w, center=FALSE, scale=colSums(w))
  
  if(!is.null(initialBiomasses))  biomasses <- initialBiomasses #if biomasses are provided use them
  if( is.null(initialBiomasses))  biomasses <- runif(dim(gAdjMat)[1], 2, 3) #if not make them up

  times <- seq(0, 4000, 5)
  solution <- lsoda_wrapper(times, biomasses, mod, verbose = FALSE)
  
  return(list(solution=solution, graph=graph))
}
