runATN <- function(graph,
                   initialBiomasses = NULL,
                   regexForProducerNodes = "^C4-|^C3-"
                   )
  {
  
  #'graph' specifies the igraph object which must contain certain vertex and edge attributes
  # most likely, this graph was produced by the custom makeGraphForATN function
  # argument regexForProducerNodes provides a regular expression which identifies 
  # producer nodes based on their names in the vertex attributes.  Default argument value matches either
  # C3- or C4- at the beginning of the name.
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
  
  mod <- create_model_Scaled(nb_s = dim(gAdjMat)[1], 
                             nb_b = sum(grepl(regexForProducerNodes, vertexAttributes$name)), 
                             BM   = vertexAttributes$BM, 
                             fw   = gAdjMat)
  mod <- initialise_default_Scaled(mod)
  
  ##TODO: modify default argument values here
  ## we intend to modify the w matrix (relative consumption weights)
  ## note the w matrix located at mod$w consists of a matrix of dimensions nb_s * (nb_s - nb_b)
  ## meaning there is a row for each species (producer/basal nodes + consumer nodes, and a column for each consumer node only)
  ## the column sums must be 1, and by default the weight is spread evenly over each of the rows (resource nodes)
  
  if(!is.null(initialBiomasses))  biomasses <- initialBiomasses #if biomasses are provided use them
  if( is.null(initialBiomasses))  biomasses <- runif(dim(gAdjMat)[1], 2, 3) #if not make them up

  times <- seq(0, 4000, 5)
  solution <- lsoda_wrapper(times, biomasses, mod, verbose = FALSE)
  
  return(list(solution=solution, graph=graph))
}
