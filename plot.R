# # generate network from list of queried genes
getNet = function(netTab, geneList, input, seed = 123){
  # input = list()
  # input$genes = c('ATM', 'ATR')
  # input$impt = 0.08
  # input$freq = 4
  # input$recr = 2
  # 
  # get queried genes
  queriedGenes = input$genes

  # get interaction table based on user's inputs
  tab = netTab %>% 
    filter(maxIptDifNormSum >= input$impt) %>% 
    filter(maxCvFreq >= input$freq)

  # remove genes not in GO:DDR
  if('GO:DDR' %in% input$plotOpts){
    tab = tab %>% 
      filter(gene_1 %in% geneList[[1]]) %>% 
      filter(gene_2 %in% geneList[[1]])
  }
  
  # remove genes not in Pendragon
  if('Pendragon' %in% input$plotOpts){
    tab = tab %>% 
      filter(gene_1 %in% geneList[[11]]) %>% 
      filter(gene_2 %in% geneList[[11]])
  }
  
  # recursively expand the list of queried genes
  if(input$recr > 1){
    for(i in 2:input$recr){
      tmp = tab %>%
        filter(gene_1 %in% queriedGenes | gene_2 %in% queriedGenes)

      queriedGenes = unique(c(queriedGenes, tmp$gene_1, tmp$gene_2))
    }
  }
  
  # get final network table
  net = tab %>%
    filter(gene_1 %in% queriedGenes | gene_2 %in% queriedGenes)

  # generate igraph object
  net = net %>% 
    select(gene_1, gene_2) %>% 
    clusterProfiler::rename(from = gene_1, to = gene_2) %>% 
    graph_from_data_frame(directed = F)
  
  # generate ggnotwork object
  set.seed(seed)
  # net = ggnetwork(net, layout = layout_with_fr(net, niter = 200))
  net = ggnetwork(net, layout = layout_with_graphopt(net))
  # net = ggnetwork(net, layout = layout_with_dh(net, maxiter = 200))
  
  # add GO:DDR and Pendragon list to the ggnetwork object
  net = net %>% 
    mutate(DDR = ifelse(name %in% geneList[[1]], 
                        'in GO:DDR', 'not in GO:DDR')) %>% 
    mutate(Pendragon = ifelse(name %in% geneList[[11]], 
                        'in Pendragon', 'not in Pendragon'))

  return(net)
}

getClust = function(clustFreq4List, input, seed = 123){
  # get cluster ID from user's input
  clustId = as.numeric(input$cluster)
  
  # get network of genes in the chosen cluster
  net = clustFreq4List[[clustId]]
  
  # filter out edges with maxCvFreq lower than user's input
  net = subgraph.edges(net, which(E(net)$maxCvFreq >= input$freqClust), 
                       delete.vertices = F)
  
  # filter out edges with maxIptDifNormSum lower than user's input
  net = subgraph.edges(net, which(E(net)$maxIptDifNormSum >= input$imptClust), 
                       delete.vertices = F)
  
  # remove genes not in GO:DDR
  if('GO:DDR' %in% input$plotOptsClust){
    net = induced_subgraph(net, (V(net)$name %in% geneList[[1]]))
  }
  
  # generate ggnetwork object
  set.seed(seed)
  # net = ggnetwork(net, layout = layout_with_graphopt(net))
  # net = ggnetwork(net, layout = layout_with_fr(net, niter = 10000))
  net = ggnetwork(net, layout = layout_with_dh(net, maxiter = 200))
  
  
  # add GO:DDR and Pendragon list to the ggnetwork object
  net = net %>% 
    mutate(DDR = ifelse(name %in% geneList[[1]], 
                        'in GO:DDR', 'not in GO:DDR')) %>% 
    mutate(Pendragon = ifelse(name %in% geneList[[11]], 
                              'in Pendragon', 'not in Pendragon'))
  
  return(net)
}

# plot network of queried genes
plotNet = function(net, input){
  p = ggplot(net, aes(x = x, y = y, xend = xend, yend = yend)) + 
    geom_edges(colour = 'grey85') + 
    geom_nodes(data = net %>% filter(name %in% input$genes), # highlight queried genes --> point border colour
               aes(shape = DDR), colour = 'grey40', size = 8) +
    geom_nodes(data = net %>% filter(name %in% input$genes), # highlight queried genes ---> point fill colour
               aes(colour = Pendragon, shape = DDR), size = 6) +
    geom_nodes(aes(colour = Pendragon, shape = DDR), size = 3) + # plot all the nodes
    geom_nodetext_repel(aes(label = name), max.overlaps = 100) +
    theme_blank() +
    theme(legend.title = element_blank(),
          legend.position = 'bottom',
          legend.text = element_text(size = 16)) +
    scale_color_d3() +
    scale_fill_npg()

  return(p)
}

# plot network of genes from a cluster
plotClust = function(net, input){
  p = ggplot(net, aes(x = x, y = y, xend = xend, yend = yend)) +
    geom_edges(colour = 'grey85') + 
    geom_nodes(aes(colour = subcluster, shape = DDR), size = 3) + # plot all the nodes
    geom_nodetext_repel(aes(label = name), max.overlaps = 100, colour = '#3c1518') +
    theme_blank() +
    theme(legend.title = element_blank(),
          legend.position = 'bottom',
          legend.text = element_text(size = 16)) +
    scale_color_d3(palette = "category20")
  
  return(p)
}