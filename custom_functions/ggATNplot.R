ggATNplot <- function(modelsolution,
                      graph,
                      facetOrder = c("C3","C4", "carnivore", "browser", "mixed feeder", "grazer")
                      ) {
  ##function to plot a solved ATN model provided using the modelsolution argument
  ##also requires the igraph object intially used to run the ATN model
  ##both this solved model and the graph are returned by custom runATN.R function writen by WAB
  ##this function plots the time series of taxon biomasses,  computes average guild biomass for the last
  ##5% of model time steps, and counts extinct taxa (those with biomass < 0.0001)
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
  
  # find the 95 percentile for time ticks, because need to summarize biomasses
  # over multiple cycles to figure out average ending biomass
  uniqueTimes <- unique(forGG$time)
  time95Percentile <- round(quantile(uniqueTimes, 0.95)) 
  nYearsForComputingEndingBiomass <- length(uniqueTimes[which(uniqueTimes >= time95Percentile)])
  
  annotations <- 
    filter(forGG, time >= time95Percentile) %>%
    group_by(taxon) %>%
    mutate(endingAverageTaxonBiomass=sum(biomass)/nYearsForComputingEndingBiomass) %>%
    select(taxon, guild, endingAverageTaxonBiomass) %>%
    distinct() %>%
    ungroup() %>%
    group_by(guild) %>%
    summarize(nExtinct = sum(endingAverageTaxonBiomass < 0.0001), 
              guildBiomass = sum(endingAverageTaxonBiomass))
  
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
  
  return(thePlot)
  
}

#ggATNplot(solved$solution, solved$graph)