plotGraphWithGuilds <- function(g){
  #color order is browser C3 C4 carnivore grazer mixed feeder, the factor level order 
  guildColors <- c("#5FAD91","#411354", "#444585","#42778C","#F9E855", "#91CF65")
  
  plot(g, 
       vertex.color=guildColors[factor(vertex_attr(g)$guild)],
       vertex.label.cex=0.5,
       edge.arrow.size=0.3,
       edge.arrow.width=0.3
  )
}
