#' given one set of genes, find the sensitivity and specificity of those genes 
#' against the target gene sets. 
#' Specificity is heavily influenced by the size of the queried gene set, i.e. 
#' bigger queried set leading to smaller specificity. Need to correct for the 
#' size of the queried set. Introduce ratioRelDif as a relative difference 
#' between specificity and the size of target sets normalised by the size of 
#' queried set. Lower ratioRelDif means better specificity. 
evalGeneSet = function(queriedSet, targetSets, targetSetNames){
  results = data.frame(targteSet = c(), sen = c(), spe = c(), meanSenSpe = c(), 
                       ratio = c(), ratioDif = c(), ratioRelDif = c())
  
  for (i in 1:length(targetSets)){
    sen = sum(queriedSet %in% targetSets[[i]]) / length(targetSets[[i]])
    spe = sum(queriedSet %in% targetSets[[i]]) / length(queriedSet)
    meanSenSpe = sqrt(sen * spe)
    ratio = length(targetSets[[i]]) / length(queriedSet)
    ratioDif = ratio - spe
    ratioRelDif = ratioDif / ratio
    
    results = rbind(results, 
                    data.frame(targetSet = targetSetNames[i], 
                               sen = sen, 
                               spe = spe, 
                               meanSenSpe = meanSenSpe, 
                               ratio = ratio, 
                               ratioDif = ratioDif, 
                               ratioRelDif = ratioRelDif))
  }
  
  return(results)
}

# given multiple gene sets, find the sensitivity and specificity of those genes 
# against the target gene sets
evalGeneSets = function(queriedSets, queriedSetNames, targetSets, targetSetNames){
  results = data.frame(queriedSet = c(), targteSet = c(), 
                       sen = c(), spe = c(), meanSenSpe = c(), 
                       ratio = c(), ratioDif = c(), ratioRelDif = c())
  
  for (i in 1:length(queriedSetNames)){
    tmp = evalGeneSet(queriedSets[[i]], targetSets, targetSetNames)
    tmp = cbind(queriedSet = queriedSetNames[i], tmp)
    
    results = rbind(results, tmp)
  }
  
  return(results)
}

# find all genes connected to a set of queried genes
getConnection2GeneSet = function(net, queriedSet, freqCrit, imptCrit){
  genes = net %>% 
    filter(maxCvFreq >= freqCrit & maxIptDifNormSum >= imptCrit) %>% 
    filter(gene_1 %in% queriedSet | gene_2 %in% queriedSet) %>% 
    pull(gene_1, gene_2) %>% 
    unique()
  
  return(genes)
}

# get all genes connected to sets of queried genes
getConnection2GeneSets = function(net, queriedSets, freqCrit, imptCrit){
  genes = list()
  
  for (i in 1:length(queriedSets)){
    genes[[i]] = getConnection2GeneSet(net, queriedSets[[i]], freqCrit, imptCrit)
  }
  
  return(genes)
}
