# # generate network from list of queried genes
getNet = function(igraphFreq, netAll, input){
  # get queried genes
  queriedGenes = input$genes

  # get interaction table based on user's frequency input
  tab = netAll[[input$freq]]

  # recursively expand the list of queried genes
  if(input$recr > 1){
    for(i in 2:input$recr){
      tmp = tab %>%
        filter(gene_1 %in% queriedGenes | gene_2 %in% queriedGenes) %>% 
        filter(edge == 'Direct') %>% 
        filter(maxImptNorm > input$impt)

      queriedGenes = unique(c(queriedGenes, tmp$gene_1, tmp$gene_2))
    }
  }

  # get all the genes connected to the queried genes
  genes = tab %>%
    filter(gene_1 %in% queriedGenes | gene_2 %in% queriedGenes) %>%
    filter(edge == 'Direct') %>% 
    filter(maxImptNorm > input$impt)

  genes = unique(c(genes$gene_1, genes$gene_2))

  # get net with user's frequency input
  net = igraphFreq[[input$freq]]

  # select edges and nodes related to the queried genes
  net = delete_vertices(net, !(V(net)$name %in% genes))
  
  # remove indirect edges
  if('Direct' %in% input$plotOpts){
    net = delete_edges(net, which(E(net)$edge == 'Indirect'))
  }
  
  # remove genes not in GO:DDR
  if('GO:DDR' %in% input$plotOpts){
    net = delete_vertices(net, V(net)$DDR == 'not in GO:DDR')
  }
  
  # remove genes not in Pendragon
  if('Pendragon' %in% input$plotOpts){
    net = delete_vertices(net, V(net)$Pendragon == 'not in Pendragon')
  }

  # generate ggnetwork object
  set.seed(123)
  net = ggnetwork(net, layout = layout_with_fr(net, niter = 10000))

  # net = ggnetwork(net, layout = layout_as_tree(delete_edges(net,
  #                                                           which(E(net)$edge == 'Indirect')),
  #                                              circular = T))

  # net = ggnetwork(net, layout = layout_with_drl(net))

  # net = ggnetwork(net, layout = layout_as_backbone(net)$xy)

  return(net)
}

# get gene cluster from user's input
getClust = function(igraphFreq, input){
  # get cluster ID from user's input
  clustId = as.numeric(input$cluster)
  
  # # get clusters from freq
  # clust = clustFreq[[input$freq]]
  
  # get network corresponding to user's input of freq
  net = igraphFreq[[input$freq]]
  
  # get cluster from cluster ID
  clust = delete_vertices(net, V(net)$cluster != clustId)
  
  # filter out indirect edges
  if('Direct' %in% input$plotOpts){
    clust = delete_edges(clust, which(E(clust)$edge == 'Indirect'))
  }
  
  # filter out genes not in GO:DDR
  if('GO:DDR' %in% input$plotOpts){
    clust = delete_vertices(clust, V(clust)$DDR != 'in GO:DDR')
  }
  
  # generate ggnetwork object
  set.seed(123)
  clust = ggnetwork(clust, layout = layout_with_fr(clust, niter = 10000))
  
  # clust = ggnetwork(clust, layout = layout_with_sugiyama(clust, maxiter = 100))
  # clust = ggnetwork(clust, layout = layout_nicely(clust))
  # clust = ggnetwork(clust, layout = layout_in_circle(clust, V(clust)))
  # clust = ggnetwork(clust, layout = layout_as_tree(clust, circular = T))
  # clust = ggnetwork(clust, layout = layout_as_tree(delete_edges(clust, 
  #                                                               which(E(clust)$edge == 'Indirect')), 
  #                                                  circular = T))
  # clust = ggnetwork(clust, layout = layout_with_dh(clust, maxiter = 100))
  # clust = ggnetwork(clust, layout = layout_as_star(clust))
  # clust = ggnetwork(clust, layout = layout_with_multilvever())
  
  return(clust)
}

# plot network of queried genes
plotNet = function(net, input){
  p = ggplot(net, aes(x = x, y = y, xend = xend, yend = yend)) +
    geom_edges(data = net %>% filter(edge == 'Indirect'), # plot indirect edges
               aes(linetype = edge),
               colour = 'grey85') +
    geom_edges(data = net %>% filter(edge == 'Direct'), # plot direct edges
               aes(linetype = edge),
               colour = 'grey50') +
    geom_nodes(data = net %>% filter(name %in% input$genes), # highlight queried genes --> point boder colour
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
plotClust = function(clust, input){
  p = ggplot(clust, aes(x = x, y = y, xend = xend, yend = yend)) +
    geom_edges(data = clust %>% filter(edge == 'Indirect'), # plot indirect edges
               aes(linetype = edge),
               colour = 'grey85') +
    geom_edges(data = clust %>% filter(edge == 'Direct'), # plot direct edges
               aes(linetype = edge),
               colour = 'grey50') +
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