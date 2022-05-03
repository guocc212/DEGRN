setwd("D://workshop//SJTU//20211205_scrna_deep//DEGRN//")

library(data.table)
library(dplyr)

select_number_of_data = function(dat, num=NULL, FULL=FALSE){
  if(FULL == FALSE){
    if(nrow(dat) > num){
      dat = dat[1:num, ]
    }
  }
  return(dat)
}

transform_gene_symbol = function(gene){
  symbol = read.delim("data//Ath_total_gene_symbol.txt", header = F, check.names = F, sep = "\t")
  names(symbol) = c("gene", "symbol", "Description")
  symbol$symbol2 = lapply(symbol$symbol, function(x){
                       a = strsplit(x, split = "[||]")[[1]]
                       return(a[1])
                    }) %>% unlist()
  gene.symbol = symbol$symbol2[match(gene, symbol$gene)]
  return(gene.symbol)
}

get_target_data = function(dat, target.number=NULL, FULL=FALSE){
  target.list = lapply(dat$geneID, function(x){
    a = strsplit(x, split = "/")[[1]]
    return(a)
  }) %>% unlist()
  target.data = data.frame(table(target.list))
  names(target.data) = c("target", "Freq")
  target.data$target = as.character(target.data$target)
  
  target.data$symbol = transform_gene_symbol(target.data$target)
  target.data = target.data[order(target.data$Freq, decreasing = T), ]
  target.data = select_number_of_data(target.data, target.number, FULL=FULL)
  return(target.data)
}

network_plot = function(net=net, Protname=Protname, main = "gene", delete.temp = TRUE, show.legend = TRUE, reverse=FALSE,
                        nodes_size = 8, nodes_border_col = "white", pt.cex=2.5,
                        edges_col = "gray", nodes_text_size = 0.5,nodes_text_color = "black" ){
  library(igraph)
  set.seed(123)
  net.graph = graph.adjacency(adjmatrix=net, mode="undirected", weighted=TRUE)
  V(net.graph)$name = Protname$label
  temp.degree = degree(net.graph)
  pname = Protname$label
  pname[temp.degree==0]=NA
  V(net.graph)$label=pname
  weights = E(net.graph)$weight
  E(net.graph)$color = E(net.graph)$weight
  
  if(delete.temp == TRUE){   
    net.graph <- delete.vertices(net.graph, which(degree(net.graph)==0)) 
  }
  
  plot(net.graph, main=main,
       layout=layout_with_fr(net.graph),  
       vertex.size=nodes_size, 
       vertex.frame.color=nodes_border_col,  
       edge.color= edges_col,
       vertex.color=Protname$cols,
       vertex.label.font=1,
       vertex.label.cex=nodes_text_size,
       vertex.label.color=nodes_text_color)
  
  if(show.legend == TRUE){
    lgd.txt = unique(Protname$type)
    lgd.col = unique(Protname$cols)
    if( reverse == TRUE){
      lgd.txt = rev(lgd.txt)
      lgd.col = rev(lgd.col)
    }
    legend("bottomleft", legend = lgd.txt, col = lgd.col,pch=16, pt.cex=pt.cex, bty="n")
  }
}

DEGRNc_tf_analysis = function(DEGRN, gene = gene, go.number = go.number, target.number = target.number, 
                              out.prefix = NULL, cols = NULL, FULL=FALSE, show.legend = TRUE){
  library(dplyr)
  library(igraph)
  # setting the color for nodes and the output prefix
  if(is.null(cols)){
    cols = c("red","#FBC723","light blue")
  }
  if(is.null(out.prefix)){
    out.prefix = gene
  }
  
  # 2. select the interested TF datasets
  gene.data = subset(DEGRN, TF == gene) %>%
    filter(pvalue < 0.05) %>%
    select_number_of_data(., go.number, FULL=FULL)
  if( nrow(gene.data)==0 ){
    message(paste0("[ERROR] There's not enough data for gene", gene, " \n. You can try the other TF genes.\n"), )
    exit(1)
  }
  write.table(gene.data, file=paste0(out.prefix, ".predicted_GOs.xls"), quote = FALSE, sep = "\t", row.names = F)
  
  # check the number of GO
  go.number = nrow(gene.data)
  
  # transform the gene symbol for TF
  gene.symbol = strsplit(gene.data$TF_Symbol, split = ",")[[1]][1]
  # obtain the target genes 
  target.data = get_target_data(gene.data, target.number, FULL = FULL)
  target.number = nrow(target.data)
  
  Protname = data.frame(name = c(gene, gene.data$Description, target.data$target),
                        label = c(gene.symbol, gene.data$Description, target.data$symbol),
                        type = c(gene.symbol, rep("GO Terms", go.number), rep("Target genes", target.number)),
                        cols = c(cols[1], rep(cols[2], go.number), rep(cols[3], target.number)))
  
  tn = 1 + go.number + target.number
  net = matrix(0, nrow = tn, ncol = tn, dimnames = list(Protname$name, Protname$name))
  net[1, 2:(go.number+1)] = 1
  for(i in 1:go.number){
    target.data.i = get_target_data(gene.data[i, ], target.number, FULL = TRUE)
    for(j in target.data.i$target){
      if( j %in% target.data$target){
        net[gene.data$Description[i], j] = 1
      }
    }
  }
  
  network_plot(net=net, Protname=Protname, main = gene.symbol, nodes_size = 8, delete.temp = TRUE, reverse=FALSE)
}
  
DEGRNc_GO_analysis = function(DEGRN, go = go, tf.number = tf.number, target.number = target.number, 
                                out.prefix = NULL, cols = NULL, FULL=FALSE, show.legend = TRUE){
    library(dplyr)
    library(igraph)
    # setting the color for nodes and the output prefix
    if(is.null(cols)){
      cols = c("red","#FBC723","light blue")
    }
    if(is.null(out.prefix)){
      out.prefix = paste0(strsplit(go, split = " ")[[1]], collapse = "_")
    }
    
    # 2. select the interested GO functions
    go.data = subset(DEGRN, Description == go)
    if( nrow(go.data) == 0){
      message(paste0("[ERROR] There's not enough  for the GO ",go, " \n. You can try the other functions.\n" ))
      exit(1)
    }
    go.data$symbol2 = transform_gene_symbol(go.data$TF)
    go.KNOWN = go.data %>%
                 filter(pvalue < 0.05) %>%
                   subset(Status == "KNOWN")
    go.NOVEL = go.data %>%
                 filter(pvalue < 0.05) %>%
                   subset(Status == "NOVEL") %>%
                     select_number_of_data(., tf.number, FULL=FULL)
    write.table(go.data, file=paste0(out.prefix, ".predicted_TFs.xls"), quote = FALSE, sep = "\t", row.names = F)
    
    # check the number of GO
    tf.number = nrow(go.NOVEL)
    known.number = nrow(go.KNOWN)
    
    # obtain the target genes 
    target.data = get_target_data(go.data, target.number, FULL = FULL)
    target.number = nrow(target.data)
    
    # 3. construct the group information for nodes in the network
    Protname = data.frame(name = c(target.data$target,  go.NOVEL$TF, go.KNOWN$TF),
                          label = c(target.data$symbol, go.NOVEL$symbol2, go.KNOWN$symbol2),
                          type = c(rep("Target genes", target.number), rep("NOVEL TFs", tf.number), rep("KNOWN TFs", known.number)),
                          cols = c(rep(cols[3], target.number), rep(cols[2], tf.number), rep(cols[1], known.number)))
    
    tn = known.number+ tf.number + target.number
    net = matrix(0, nrow = tn, ncol = tn, dimnames = list(Protname$name, Protname$name))
    # 4. calculate the matrix for the links between TFs and targets
    total = rbind(go.KNOWN, go.NOVEL)
    for(i in 1:nrow(total)){
      target.data.i = get_target_data(total[i, ], target.number, FULL = TRUE)
      for(j in target.data.i$target){
        if( j %in% target.data$target){
          net[total$TF[i], j] = 1
        }
      }
    }
    
    network_plot(net=net, Protname=Protname, main = go, nodes_size = 8, delete.temp = TRUE, reverse=TRUE)

}

DEGRN_plot_pipeline = function(modes = "go", go = NULL, gene = NULL, out.prefix=nULL,
                               go.number = NULL, tf.number = NULL, target.number = NULL,
                               cols = NULL){
  library(data.table)
  # 1. read the results of DEGRN
  DEGRN = fread("data//Network_by_DEGRN.xls", check.names = F, header = T, sep = "\t")

  # 2. choose the mode for users
  if( modes == "go"){
    DEGRNc_GO_analysis(DEGRN, go = go, tf.number = tf.number, target.number = target.number)
  }else if( modes == "gene"){
    DEGRNc_tf_analysis(DEGRN, gene = gene, go.number = go.number, target.number = target.number)
  }else{
    message("[ERROR] \n There are two mode: go or gene.")
    message("[ERROR] \n 1. mode = 'go':       You can input the interested BP of GO for the TFs which was involved in this GO, including the KNOWN TFs and NOVEL TFs. ")
    message("[ERROR] \n 2. mode = 'gene':       You can input the interested TFs for the BP which was predicted by DEGRN, including the KNOWN functions and NOVEL functions. ")
    
  }
}
