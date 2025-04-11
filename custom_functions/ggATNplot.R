ggATNplot <- function(modelsolution,
                      graph,
                      facetOrder = c("C3","C4", "carnivore", "browser", "mixed feeder", "grazer")
                      ) {
  ##function to plot a solved ATN model provided using the modelsolution argument
  ##also requires the igraph object intially used to run the ATN model
  ##both this solved model and the graph are returned by custom runATN.R function writen by WAB
  require(igraph)
  require(dplyr)
  require(ggplot2)
  require(tidyr)
  
  vertexAttributes <- as.data.frame(vertex_attr(graph))
  names(vertexAttributes)[which(names(vertexAttributes)=="name")] <- "taxon"
  
  if(!all(facetOrder %in% unique(vertexAttributes$guild))) stop("all strings in facetOrder must be found among the guilds in the vertex attributes of g")
  if(!"igraph" %in% class(graph)) stop("ggATNplot requires an object of class igraph")
  if(!"guild" %in% names(vertexAttributes)) stop("can't find guilds (variable name 'guild') in vertex attributes")
  if(!"deSolve" %in% class(modelsolution))  stop("modelSolution must be an object of class deSolve")


  forGG <- as.data.frame(modelsolution)
  names(forGG)[2:(length(vertexAttributes$taxon) + 1)] <- vertexAttributes$taxon
  forGG <- pivot_longer(forGG, names_to="taxon", values_to="biomass", -time)
  forGG <- left_join(forGG, vertexAttributes, by = join_by(taxon))
  
  forGG$guild <- factor(forGG$guild, levels=facetOrder, ordered = TRUE)
  
  annotations <- 
    filter(forGG, time==max(time)) %>%
    group_by(guild) %>%
    summarize(nExtinct = sum(biomass==0),
              guildBiomass=sum(biomass))
  
  thePlot <- 
    ggplot(forGG, aes(x=time, y=biomass, group=taxon)) + 
    geom_line(alpha=0.6, aes(color=guild)) + 
    geom_text(data=annotations, 
              aes(label=sprintf("Guild biomass = %0.1f\nnExtinct = %d", round(guildBiomass,1), nExtinct)),
              x=median(forGG$time), 
              y=max(forGG$biomass)*0.9, 
              size=4,
              inherit.aes = FALSE)  + 
    facet_wrap(~guild) + 
    guides(color='none') + 
    theme_bw(12) 
  
  print(thePlot)
  
}

#ggATNplot(solved$solution, solved$graph)