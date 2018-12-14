#Random effects to networks

#' re.net function
#'
#' This function will take the random effect matrix from a mixed model and create networks.
#' @param df.int Dataframe of the random effects where each row is a sample for each dyadic relationship. The column names of such a dataframe should have both individual names seperated by a "_".
#' @param summary If summary is F (Default) then a list of networks is returned, one for each sample. If summary is T then one network is returned with mean edge weights.
#' @param directed Whether the relationships are considered as directed or not.
#' @export
#' @import igraph
#' @importFrom utils setTxtProgressBar txtProgressBar
#'
re.net<-function(df.int, summary=FALSE, directed = TRUE){

  #create list of node names and give them each a node number (needed for adding edges)
  col.names <- colnames(df.int)
  node.names <- vector()
  for(i in 1:length(col.names)){
    split.name<-strsplit(col.names[i], '_')
    if(length(split.name[[1]])==2){
      node.names[length(node.names)+1]<-split.name[[1]][1]
      node.names[length(node.names)+1]<-split.name[[1]][2]
    } else {
      print("error colnames not read properly")
    }
  }
  unique.node.names<-unique(node.names)
  df.names<-data.frame(names=unique.node.names, numb=seq(1,length(unique.node.names),by=1))


  if(summary==FALSE){

    #create a list of graphs
    graph.list <- list()

    #create progress bar
    pb <- txtProgressBar(min = 0, max = nrow(df.int), style = 3)

    #for each row i
    for (i in 1:nrow(df.int)){

      #create empty graph and start adding edges
      if(directed == FALSE) {
        graph.temp<-as.undirected(make_empty_graph(n=length(unique.node.names) ))
      } else {
        graph.temp<-as.directed(make_empty_graph(n=length(unique.node.names) ))
      }


      #take each column j
      for (j in 1:ncol(df.int)){

        #get names of nodes and their respective numbers
        split.name<-strsplit(colnames(df.int)[j], '_')
        node1 <- df.names[which(df.names$names==split.name[[1]][1]),]$numb
        node2 <- df.names[which(df.names$names==split.name[[1]][2]),]$numb

        #get value of
        weight <- df.int[i,j]

        #add to graph
        graph.temp <- graph.temp %>% add_edges(c(node1,node2), weight=df.int[i,j])

        # update progress bar
        setTxtProgressBar(pb, i)
      }


      #add the graph to the list of graphs
      V(graph.temp)$label <- as.character(df.names$names)
      graph.list[[i]] <- graph.temp

    }
    return(graph.list)
  } else {

    #Calculate mean network

    #create empty graph and start adding edges
    if(directed == FALSE) {
      graph.mean<-as.undirected(make_empty_graph(n=length(unique.node.names) ))
    } else {
      graph.mean<-as.directed(make_empty_graph(n=length(unique.node.names) ))
    }

    #take each column j
    for (j in 1:5){ #ncol(df.int)

      #get names of nodes and their respective numbers
      split.name<-strsplit(colnames(df.int)[j], '_')
      node1 <- df.names[which(df.names$names==split.name[[1]][1]),]$numb
      node2 <- df.names[which(df.names$names==split.name[[1]][2]),]$numb

      #get value of
      weight <- mean(df.int[,j])

      #add to graph
      graph.mean <- graph.mean %>% add_edges(c(node1,node2), weight=df.int[i,j])
    }

    V(graph.mean)$label <- as.character(df.names$names)

    return(graph.mean)
  }

}
