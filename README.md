DDRnet utilises the xgboost model to explore genetic interactions with 
and amongst DDR genes, offering a complementary approach to CRISPR screen 
analyses. The xgboost model identifies interactions with a target gene, such as 
ATM, by predicting its z-scores based on those of other genes across all screens
within the DDRcs. Through a training process, the xgboost model autonomously 
constructs an ensemble of regression trees with a select few important genes, 
which are considered to interact with the target gene (e.g. ATM). DDRnet applies
this procedure to all genes in the human genome to infer genetic interactions 
with DDR genes.<br><br>

The DDRnet web app offers two modes of network visualization: (1) gene query and
(2) cluster. In gene query mode, users select a list of genes of interest and 
query all the genes interacting with these input genes. Users have the 
flexibility to adjust the network confidence by setting frequency and importance
threshold. Higher frequency and importance yield more stringent analyses but 
result in fewer interacting genes. Additionally, users can choose from various 
display options to customize the network visualisation. In cluster mode, users 
can explore subsets of Pendragon genes, which have been pre-clustered offline 
(frequency \u2265 3 and importance \u2265 0.08). Adjusting the frequency and importance 
threshold results in different display of gene interactions within the clusters.
<br><br>

Since we filtered out some genes to maintain high consistency of the data across
multiple CRISPR screens, DDRnet does not encompass a complete genetic network of
all gene pairs. If you cannot locate your genes of interest in the networks or 
if they only interact with few others, it is likely that these genes were not 
included in our analysis. The list of genes utilised in the DDRnet can be found 
in the “Genes” tab.