## this is an example of the output file from Alyssa's random forest models
## predicting link probabilities for Upper Burgi taxa
##UBurgiLinks <- read.csv("./upperBurgi_species.csv")

makeGraphForATN <- function(linkDF, 
                      preyColName="scientificNamePrey", 
                      carnivoreColName="scientificNameCarni",
                      guildColName="Dietary.Guild",
                      linkProbColName="RelationshipProbability",
                      herbBMcolName="herbivore_BM.kg",
                      carniBMcolName="carnivore_BM.kg",
                      nC4=5,
                      nC3=5,
                      C4bodymass=0.01,
                      C3bodymass=0.1,
                      linkProbThreshold = 0.05
                      ) {
  ## function which takes an edge list (linkDF), adds a user-specified number of producer nodes
  ## and converts the resulting edgelist to a igraph object, preserving link probabilities 
  ## as edge weights in the output, retrievable using edge_attr(output)
  ## vertex attributes of output include the guild assignment of each taxon, BodyMass, and the taxon name
  ## linkDF should contain resource -> consumer links and
  ## has to have, at a minimum, columns names matching the strings given
  ## in arguments carnivoreColName, preyColName, guildColName, herbBMcolName, carniBMcolName, and linkProbColName
  ## nC4 and nC3 arguments specifies how many C4 and C3 producer nodes should be created
  ## C4bodymass and C3bodymass give the body mass (in kg) of individual producers

  if(!preyColName %in% names(linkDF))      stop('value of the preyColName argument does not match a column name in linkDF')
  if(!carnivoreColName %in% names(linkDF)) stop('value of the carnivoreColName argument does not match a column name in linkDF')
  if(!guildColName %in% names(linkDF))     stop('value of the guildColName argument does not match a column name in linkDF')
  if(!linkProbColName %in% names(linkDF))  stop('value of the linkProbColName argument does not match a column name in linkDF')
  
  require(igraph)
  require(dplyr)
  
  ## note the use of {{ }}, which allow us to use arguments to the makeGraph function
  ## with dplyr functions as though these were column names in linkDF
  ## https://dplyr.tidyverse.org/articles/programming.html
  
  distinctPrey <- select(linkDF, 
                         taxon = {{ preyColName }}, 
                         guild = {{ guildColName }},
                         BM = {{ herbBMcolName }}
                         ) %>%
    distinct()

  
  
  ##select relevant columns from input dataframe and rename them resourceTaxon, consumerTaxon, guild
  existingLinks <- select(linkDF, 
                          resourceTaxon={{ preyColName }}, 
                          consumerTaxon={{ carnivoreColName }},  
                          linkProb = {{ linkProbColName }}
                          ) %>%
    filter(linkProb > linkProbThreshold)
  
  grazers <- filter(distinctPrey, guild == "grazer")
  C4taxa <- paste("C4",1:nC4, sep="-")
  grazerLinks <- data.frame(resourceTaxon=rep(C4taxa, each=nrow(grazers)),
                        consumerTaxon=rep(grazers$taxon, times=nC4),
                        linkProb=1)
  
  browsers <- filter(distinctPrey, guild == "browser")
  C3taxa <- paste("C3",1:nC3, sep="-")
  browserLinks <- data.frame(resourceTaxon=rep(C3taxa, each=nrow(browsers)),
                            consumerTaxon=rep(browsers$taxon, times=nC3),
                            linkProb = 1)
  
  MFs <- filter(distinctPrey, guild == "mixed feeder")
  MFC4links <- data.frame(resourceTaxon=rep(C4taxa, each=nrow(MFs)),
                            consumerTaxon=rep(MFs$taxon, times=nC4),
                            linkProb=1
                          )
  
  MFC3links <- data.frame(resourceTaxon=rep(C3taxa, each=nrow(MFs)),
                             consumerTaxon=rep(MFs$taxon, times=nC3),
                             linkProb=1
                          )
  
  finalLinkDF <- do.call(rbind, list(
    ## note because they are listed with existing links at the end, the 
    grazerLinks, browserLinks, MFC4links, MFC3links,existingLinks
  ))
  
  g <- graph_from_data_frame(finalLinkDF)
  
  
  distinctPreds <- select(linkDF, 
                          taxon = {{ carnivoreColName }},
                          BM = {{ carniBMcolName }}
                          ) %>%
    distinct()
  distinctPreds$guild <- "carnivore"
  
  distinctProducerNodes <- data.frame(
    taxon = c(C3taxa, C4taxa), 
    guild = c(rep("C3", nC3), rep("C4", nC4)),
    BM =    c(rep(C3bodymass, nC3), rep(C4bodymass, nC4))
  )
  
  
  
  ##pulling the vertex attributes from the igraph object makes sure vertices are in the correct order
  vertexAttrDF <- data.frame(name=vertex_attr(g)$name)
  
  vertexAttrDF <- left_join(vertexAttrDF, 
                            rbind(distinctPreds, distinctPrey, distinctProducerNodes), 
                            by=c("name" = "taxon"))
  
  vertex_attr(g, name="guild") <- vertexAttrDF$guild
  vertex_attr(g, name="BM") <- vertexAttrDF$BM
  return(g)
  
}