# helper function to add abitrary weights to the network ----
edge.weights <- function(community, network, weight.within = 100, weight.between = 1){
  newNet = network
  
  # add weights to newNet
  bridges <- crossing(communities = community, graph = network)
  weights <- ifelse(test = bridges, yes = weight.between, no = weight.within)
  
  E(newNet)$weight = weights  
  
  # add more edges to newNet
  for(i in unique(membership(community))){
    GroupV = which(membership(community) == i)
    LC_Grouped = add_edges(newNet, combn(GroupV, 2), attr = list(weight = weight.within))
  } 
  
  return(newNet) 
}

# cluster network using Louvain
clustNet = function(igraphNet, minSize = 20, seed = 123){
  set.seed(123)
  cluster = cluster_louvain(igraphNet)
  
  smallClusters = which(table(cluster$membership) < 20)
  keep = V(igraphNet)[!(cluster$membership %in% smallClusters)]
  
  igraphNetPruned = induced_subgraph(igraphNet, keep)
  
  clusterPruned = cluster_louvain(igraphNetPruned)
  
  igraphNetClust = edge.weights(clusterPruned, igraphNetPruned, weight.within = 1e20, weight.between = 1)
  
  V(igraphNetClust)$cluster = clusterPruned$membership
  
  return(igraphNetClust)
}

# add geneList proportion to network cluster ----
propClustNet = function(igraphNetClust, geneList, seed = 123){
  igraphNetSuper = igraph::simplify(contract(igraphNetClust, V(igraphNetClust)$cluster))
  D = unname(degree(igraphNetSuper))
  
  coverage = c()
  size = c()
  for(i in unique(V(igraphNetClust)$cluster)){
    community = V(igraphNetSuper)[i]$name[[1]]
    
    coverage[i] = sum(community %in% geneList) / length(community)
    # coverage[i] = sum(community %in% geneList) / length(geneList)
    
    size[i] = length(community)
  }
  
  V(igraphNetSuper)$coverage = coverage
  V(igraphNetSuper)$ComSize = size
  
  V(igraphNetSuper)$name = c(1:length(V(igraphNetSuper)))
  
  set.seed(seed)
  ggNetSuper = ggnetwork(igraphNetSuper, layout = layout_with_fr(igraphNetSuper, niter = 50000))
  # ggNetSuper = ggnetwork(igraphNetSuper, layout = layout_with_mds(igraphNetSuper))
  
  
  return(ggNetSuper)
}

# plot clusters
plotClusters = function(net, clust, path, prefix, seed = 123){
  # plot all dominant clusters
  for (i in 1:max(V(clust)$cluster)){
    cat('Plot cluster: ', i, '\n')
    
    net_clust = induced_subgraph(clust, V(clust)$cluster == i)
    
    set.seed(seed)
    net_clust = ggnetwork(net_clust, layout = layout_with_graphopt(net_clust))
    
    ggplot(net_clust, aes(x = x, y = y, xend = xend, yend = yend)) + 
      geom_edges(colour = 'grey') + 
      geom_nodes(colour = 'steelblue') + 
      geom_nodetext_repel(aes(label = name)) + 
      theme_blank()
    
    ggsave(paste0(path, 'clust_', i, '_', prefix, '.jpeg'))
  }
  
  # plot cluster 0
  cat('Plot cluster: 0\n')
  
  net_clust = induced_subgraph(net, !(V(net)$name %in% V(clust)$name))
  
  set.seed(seed)
  net_clust = ggnetwork(net_clust, layout = layout_with_graphopt(net_clust))
  
  ggplot(net_clust, aes(x = x, y = y, xend = xend, yend = yend)) + 
    geom_edges(colour = 'grey') + 
    geom_nodes(colour = 'steelblue') + 
    geom_nodetext_repel(aes(label = name)) + 
    theme_blank()
  
  ggsave(paste0(path, 'clust_0', '_', prefix, '.jpeg'))
  
  return('Done')
}

# plot word cloud
plotWordCloudsClust = function(net, clust, path, prefix, seed = 123){
  # gene set analysis ----
  cat('Gene set analysis:\n')
  
  clust_entrez = list()
  clust_ego = list()
  
  for (i in 1:max(V(clust)$cluster)){
    cat('\tCluster: ', i, '\n')
    clust_entrez[[i]] = bitr(V(clust)[V(clust)$cluster == i]$name, 
                             fromType = "SYMBOL", 
                             toType = c("ENTREZID"), 
                             OrgDb = org.Hs.eg.db)
    
    clust_ego[[i]] = enrichGO(gene = clust_entrez[[i]]$ENTREZID, 
                              # universe = zImpute$gene, 
                              OrgDb = org.Hs.eg.db, 
                              ont = 'BP', 
                              pAdjustMethod = 'BH', 
                              pvalueCutoff = 0.001, 
                              qvalueCutoff = 0.001)
    
    clust_ego[[i]] = summary(clust_ego[[i]])
  }
  
  # Word cloud ----
  # Generate world cloud from the GO enrichment analysis using text mining
  cat('Generate word clouds:\n')
  
  for (i in 1:max(V(clust)$cluster)){
    cat('\tCluster: ', i, '\n')
    
    if(dim(clust_ego[[i]])[1] < 1) next
    
    ## create text object ----
    docs = Corpus(VectorSource(clust_ego[[i]]$Description))
    
    inspect(docs)
    
    ## preprocessing ----
    ### text transformation ----
    toSpace <- content_transformer(function (x , pattern ) gsub(pattern, " ", x))
    docs <- tm_map(docs, toSpace, "/")
    docs <- tm_map(docs, toSpace, "@")
    docs <- tm_map(docs, toSpace, "\\|")
    docs <- tm_map(docs, toSpace, "'")
    docs <- tm_map(docs, toSpace, "-")
    docs <- tm_map(docs, toSpace, "\\(")
    docs <- tm_map(docs, toSpace, "\\)")
    
    ### text cleaning ----
    # Convert the text to lower case
    docs <- tm_map(docs, content_transformer(tolower))
    # Remove numbers
    docs <- tm_map(docs, removeNumbers)
    # Remove english common stopwords
    docs <- tm_map(docs, removeWords, stopwords("english"))
    # Remove punctuations
    docs <- tm_map(docs, removePunctuation)
    # Eliminate extra white spaces
    docs <- tm_map(docs, stripWhitespace)
    # Text stemming
    # docs <- tm_map(docs, stemDocument)
    
    ## Build a term-document matrix ----
    dtm <- TermDocumentMatrix(docs)
    m <- as.matrix(dtm)
    v <- sort(rowSums(m),decreasing=TRUE)
    d <- data.frame(word = names(v),freq=v)
    head(d, 10)
    
    ## Generate a word cloud ----
    set.seed(seed)
    setEPS()
    postscript(paste0(path, 'wordCloud_clust_', i, '_', prefix, '.eps'))
    wordcloud(words = d$word, freq = d$freq, min.freq = 1,
              max.words=200, random.order=FALSE, rot.per=0.35,
              colors=brewer.pal(8, "Dark2"))
    dev.off()
  }
}
