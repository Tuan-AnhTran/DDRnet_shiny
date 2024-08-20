# # generate network from list of queried genes
getNet = function(netTab, geneList, input){
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
  set.seed(123)
  # net = ggnetwork(net, layout = layout_with_fr(net, niter = 10000))
  net = ggnetwork(net, layout = layout_with_graphopt(net))
  
  # add GO:DDR and Pendragon list to the ggnetwork object
  net = net %>% 
    mutate(DDR = ifelse(name %in% geneList[[1]], 
                        'in GO:DDR', 'not in GO:DDR')) %>% 
    mutate(Pendragon = ifelse(name %in% geneList[[11]], 
                        'in Pendragon', 'not in Pendragon'))

  # net = ggnetwork(net, layout = layout_as_tree(delete_edges(net,
  #                                                           which(E(net)$edge == 'Indirect')),
  #                                              circular = T))

  # net = ggnetwork(net, layout = layout_with_drl(net))

  # net = ggnetwork(net, layout = layout_as_backbone(net)$xy)

  return(net)
}

# get gene cluster from user's input
# getClust = function(netTab, clustFreq4List, geneList, input){
#   # get cluster ID from user's input
#   clustId = as.numeric(input$cluster)
#   
#   # get network of genes in the chosen cluster
#   net = netTab %>% 
#     filter(gene_1 %in% clustFreq4List[[clustId]]) %>% 
#     filter(gene_2 %in% clustFreq4List[[clustId]]) %>% 
#     select(gene_1, gene_2, maxCvFreq, maxIptDifNormSum) %>% 
#     clusterProfiler::rename(from = gene_1, to = gene_2) %>% 
#     graph_from_data_frame(directed = F)
#   
#   # filter out edges with maxCvFreq lower than user's input
#   net = delete_edges(net, which(E(net)$maxCvFreq < input$freqClust))
#   
#   # filter out edges with maxIptDifNormSum lower than user's input
#   net = delete_edges(net, which(E(net)$maxCvFreq < input$imptClust))
#   
#   # remove genes not in GO:DDR
#   if('GO:DDR' %in% input$plotOptsClust){
#     net = induced_subgraph(net, (V(net)$name %in% geneList[[1]]))
#   }
#   
#   # generate ggnetwork object
#   set.seed(123)
#   net = ggnetwork(net, layout = layout_with_graphopt(net))
#   # net = ggnetwork(net, layout = layout_with_fr(net, niter = 10000))
#   
#   # add GO:DDR and Pendragon list to the ggnetwork object
#   net = net %>% 
#     mutate(DDR = ifelse(name %in% geneList[[1]], 
#                         'in GO:DDR', 'not in GO:DDR')) %>% 
#     mutate(Pendragon = ifelse(name %in% geneList[[11]], 
#                               'in Pendragon', 'not in Pendragon'))
#   
#   return(net)
# }

getClust = function(clustFreq4, geneList, input){
  # input = list()
  # input$cluster = 25
  # input$freqClust = 5
  # input$imptClust = 0.03

  
  # get cluster ID from user's input
  clustId = as.numeric(input$cluster)
  
  # get network of genes in the chosen cluster
  net = induced_subgraph(clustFreq4, V(clustFreq4)$cluster == clustId, 
                         impl = 'create_from_scratch')
  
  # filter out edges with maxCvFreq lower than user's input
  # net = delete_edges(net, which(E(net)$maxCvFreq < input$freqClust))
  net = subgraph.edges(net, which(E(net)$maxCvFreq >= input$freqClust), 
                       delete.vertices = F)
  
  # filter out edges with maxIptDifNormSum lower than user's input
  # net = delete_edges(net, which(E(net)$maxIptDifNormSum < input$imptClust))
  net = subgraph.edges(net, which(E(net)$maxIptDifNormSum >= input$imptClust), 
                       delete.vertices = F)
  
  # remove genes not in GO:DDR
  if('GO:DDR' %in% input$plotOptsClust){
    net = induced_subgraph(net, (V(net)$name %in% geneList[[1]]))
  }
  
  # generate ggnetwork object
  set.seed(123)
  net = ggnetwork(net, layout = layout_with_graphopt(net))
  # net = ggnetwork(net, layout = layout_with_fr(net, niter = 10000))
  
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
    geom_nodes(aes(colour = Pendragon, shape = DDR), size = 3) + # plot all the nodes
    geom_nodetext_repel(aes(label = name), max.overlaps = 100, colour = '#3c1518') +
    theme_blank() +
    theme(legend.title = element_blank(),
          legend.position = 'bottom',
          legend.text = element_text(size = 16)) +
    scale_color_d3() +
    scale_fill_npg()
  
  return(p)
}