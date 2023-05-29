#' Estimate the biological meaningfulness of cluster splits
#'
#' @description Estimate the support of one cluster grouped into several subclusters by the external dataset
#' @param obj a Seurat object with clustering at multiple resolutions
#' @param external a external dataset; bulk or one Seurat object for single cell RNA-seq; if bulk, row is the gene, column is the sample. 
#' @param prefix the prefix of Seurat clustering resolution names (default: RNA_snn_res.)
#' @param logfc.threshold the cutoff of log2fold change to identify markers (default: 0.1)
#' @param min.pct the percent difference to identify markers (default: 0.1)
#' @param adjp the cutoff of adjusted p-value to identify markers (default: 0.01)
#' @param topn the number of marker genes (default: 20)
#' @param bulk indicate the external is from bulk RNAseq (TRUE) or single-cell RNA-seq (FALSE) (default:TRUE)
#' @return a Seurat object with the score and p-value associated with each subcluster stored in misc$BioCluster
#' @export


BioCluster<-function(obj,external,prefix="RNA_snn_res.",logfc.threshold=0.1,min.pct=0.1, adjp=0.01,topn=20, bulk=TRUE){
  
 
  resolution<-gsub(prefix,"",colnames(obj@meta.data)[grep(prefix,colnames(obj@meta.data))])
  
  resolution<-sort(as.numeric(resolution))
  if(length(resolution)<=1) {
    stop("the object doesnot have multiple resolutions or the prefix is incorrect")
  }
  
  result<-NULL
  if (!bulk) {
    transferanchors <- FindTransferAnchors(reference =obj, query = external,dims = 1:30, reference.reduction = "pca")
    external_data<-GetAssayData(external,slot="data")
    genenames<-rownames(external_data)
  } else {genenames<-rownames(external)}
  
  if (length(intersect(genenames,rownames(GetAssayData(obj,assay="RNA",slot="data"))))<1000){
    stop("the genes from the Seurat object and the external have few overlaps (<1000)")
  }
  
  # if the external is single-cell RNAseq, find transfer anchors
  
  
  Idents(obj) <-paste0(prefix,resolution[1])
  if (length(levels(obj))>1) {
      markers<-FindAllMarkers(obj)
      topmarkers<-markers %>% group_by(cluster) %>% top_n(n=topn,wt=avg_log2FC)
     
      geneind<-match(topmarkers$gene,genenames)
      if (bulk) {external_sub<-external[geneind[!is.na(geneind)],]
                 cluster<-topmarkers$cluster[!is.na(geneind)]
         }else {
           
           external_sub<-as.matrix(external_data[geneind[!is.na(geneind)],])
           cluster<-topmarkers$cluster[!is.na(geneind)]
           genevar<-apply(external_sub,1,var)
           external_sub<-external_sub[genevar!=0,]
           cluster<-cluster[genevar!=0]
           
      }
      
      cor_external<-cor(t(external_sub),method="spearman")
      obs_coexp<-as.vector(CalScore(cor_external,cluster))
      rand_coexp<-replicate(1000,CalScore(cor_external,cluster[sample(1:length(cluster),length(cluster))])) 
      pval_coexp<-apply((rand_coexp)>obs_coexp,1,sum)/1000
      result<-rbind(result,cbind(from_clust=0,to_clust=levels(cluster),from_res=0,to_res=resolution[1],score=obs_coexp,pval=pval_coexp))
      
      }
  
    for (resind in 1:(length(resolution)-1)){
        res1<-paste0(prefix,resolution[resind])
        res2<-paste0(prefix,resolution[resind+1])
        tableind<-table(obj[[res1]][,1],obj[[res2]][,1])
        tableprop<-proportions(tableind,margin=2)
        splitind<-apply(tableprop,1,function(x){names(x)[which(x>0.5)]},simplify=FALSE)
        splitlength<-sapply(splitind,length)
        indgreater2<-which(splitlength>=2)
        Idents(obj)<-res2
        if (length(indgreater2)>0){
           for (indi in 1:length(indgreater2)) {
             clustin<-indgreater2[indi]
              
          
          
                 markers<-FindAllMarkers(subset(obj,idents=splitind[[clustin]]),pseudocount.use = 1,min.pct = min.pct,logfc.threshold = logfc.threshold)
                  if(nrow(markers)>=1) {
                      topmarkers<-markers %>% group_by(cluster) %>% filter(p_val_adj<adjp & pct.1-pct.2>=min.pct & avg_log2FC>logfc.threshold) %>% top_n(n=topn,wt=avg_log2FC)} else {topmarkers<-NULL}
                  if (!is.null(topmarkers)) {
                      if (nrow(topmarkers)>2) {
                          geneind<-match(topmarkers$gene,genenames)
                          if (bulk) {
                               external_sub<-external[geneind[!is.na(geneind)],]
                               cluster<-topmarkers$cluster[!is.na(geneind)]
                          }else {
                            predictions <- TransferData(anchorset = transferanchors, refdata = obj@active.ident, dims = 1:30)
                            external_sub<-external_data[,predictions$predicted.id %in% splitind[[clustin]]]
                            external_sub<-as.matrix(external_sub[geneind[!is.na(geneind)],])
                            cluster<-topmarkers$cluster[!is.na(geneind)]
                            genevar<-apply(external_sub,1,var)
                            external_sub<-external_sub[genevar!=0,]
                            cluster<-cluster[genevar!=0]
                          }
                          if(ncol(external_sub)>=20){
                              cor_external<-cor(t(external_sub),method="spearman")
                              obs_coexp<-as.vector(CalScore(cor_external,cluster))
                              clustersum<-table(cluster)
                              if (prod(clustersum[clustersum>0]>=3)>0){
                                  rand_coexp<-replicate(1000,CalScore(cor_external,cluster[sample(1:length(cluster),length(cluster))])) 
                                  pval_coexp<-apply((rand_coexp)>obs_coexp,1,sum)/1000
                              }else pval_coexp<-NA
                             result<-rbind(result,cbind(from_clust=names(splitind[clustin]),to_clust=splitind[[clustin]],from_res=resolution[resind],to_res=resolution[resind+1],score=obs_coexp,pval=pval_coexp))
                          }
                      }  
                 }
              
           }
        }
       }
  colnames(result)<-c("from_clust","to_clust",paste0("from_",prefix),paste0("to_",prefix),"score","pval")
  result<-data.frame(result)
  result[,3:6]<-apply(result[,3:6],2,as.double)
  obj@misc$BioCluster<-result
  return(obj)
}



#' the BioCluster result Visualization 
#'
#' @description Visualizing how clusters at a low resolution are grouped into several subclusters labeled with scores and p-values, suggesting the support of the partition by the external dataset.
#' @param obj Seurat object
#' @param node.size node size; either a numeric value or NULL ;  if NULL, the size of each node is determined by the size of the cluster (default:NULL) 
#' @param node.color node color; either a color or  NULL ; if NULL, the color of each node is determined by the resolution (default:NULL)
#' @param node.text  node text; either "" or NULL; if NULL, the text is the cluster ID; if "", node text is not shown(default:NULL)
#' @param return specifying what to return, either "plot" (a ggplot object) or "graph" (a tbl_graph object) (default:"plot")
#' @return a ggplot object or a tbl_graph object
#' @export

PlotBioCluster<-function(obj,node.size=NULL, node.color=NULL, node.text=NULL,return="plot") {
  
  bioresult<-obj@misc$BioCluster
  if (is.null(bioresult)){
    stop("BioCluster not computed for the object")
  } 
  prefix=gsub("from_","",colnames(bioresult)[3])
  resolution<-gsub(prefix,"",colnames(obj@meta.data)[grep(prefix,colnames(obj@meta.data))])
  
  resolution<-sort(as.numeric(resolution))
  


  Idents(obj) <-paste0(prefix,resolution[1])
  
  #add root node
  if (length(levels(obj))>1) {
    obj[[paste0(prefix,"0")]]<-0
  }
  
  graph<-clustree::clustree(obj,prefix=prefix,return="graph")
  graph <- graph %>% activate("edges") %>% filter(in_prop>0.5) %>% filter(is_core)
 
  ##remove isolate nodes
  graph<-graph %>% activate("nodes") %>% mutate(isolate=node_is_isolated()) %>% filter(!isolate)
  
  # combine the BioCluster result 
  bioresult<-data.frame(bioresult)
  bioresult[,3:6]<-apply(bioresult[,3:6],2,as.double)
  graph<- graph %>% activate("edges") %>% left_join(bioresult,by=colnames(bioresult)[1:4])
  
  ## 
  if (return=="graph") {return(graph)}
  
  node.size<-ifelse(is.null(node.size),as.name("size"),node.size)
  node.color<-ifelse(is.null(node.color),as.name(prefix),node.color)
  
  graphplot<-ggraph(graph,layout="tree")+geom_edge_diagonal()+
    geom_edge_diagonal(angle_calc = 'along',label_dodge = unit(-1.5, 'mm'),label_size=3,aes(filter=!is.na(score),label=paste0("p=",round(pval,2))))+
    geom_edge_diagonal(angle_calc = 'along',label_dodge = unit(1.5, 'mm'),label_size=3,aes(filter=!is.na(score),label=paste0("s=",round(score,2))))+
    scale_size(range=c(4,10))+theme(panel.background = element_blank())+
    geom_node_point(aes_(size=node.size,color=node.color))
 
  if (is.null(node.text))
  {graphplot<-graphplot+geom_node_text(aes_(label=as.name("cluster")))} 
  if (return=="plot") {return(graphplot)}
  
  
}

#' Pruning the clustering tree by removing those splits not supported by external dataset
#'
#' @description Visualizing the pruned result
#' @param obj Seurat object
#' @param score.cutoff the cutoff of score to keep the split: (default:0.2)
#' @param pval.cutoff the cutoff of p-value to keep the split; (default:0.05)
#' @param node.size node size; either a numeric value or NULL ;  if NULL, the size of each node is determined by the size of the cluster (default:NULL) 
#' @param node.color node color; either a color or  NULL ; if NULL, the color of each node is determined by the resolution (default:NULL)
#' @param node.text  node text; either NULL, "cluster" or others; if NULL, the node text is the new cluster id after pruning; if "cluster", node text is the cluster ID at each resolution; if others, node text is not shown(default:NULL)
#' @return a list containing a Seurat object with the final clustering result after pruning and a ggplot object 
#' @export

BioCluster_prune<-function(obj,score.cutoff=0.2,pval.cutoff=0.05,node.size=NULL, node.color=NULL, node.text=NULL) {
  
  bioresult<-data.frame(obj@misc$BioCluster)
 
  
  if (is.null(bioresult)){
    stop("BioCluster not computed for the object")
  }
  
  prefix=gsub("from_","",colnames(bioresult)[3])
  resolution<-gsub(prefix,"",colnames(obj@meta.data)[grep(prefix,colnames(obj@meta.data))])
  resolution<-sort(as.numeric(resolution))
  
  
  
  graph<-PlotBioCluster(obj,return="graph")
  
  ##remove nodes with no children nodes (like resolutin 0.2 cluster 9 in lung data)
  graph <- graph %>% activate("nodes") %>% mutate(node_nochild=node_is_leaf() &  eval(as.symbol(prefix)) != max(resolution)) %>% filter(!node_nochild)
  nodes <-as.list(graph)$nodes
  bioresult<-bioresult[paste0(prefix,bioresult[,4],"C",bioresult[,2]) %in% nodes$node,]
  
  
  nodesinsplit<-paste0(prefix,bioresult[,4],"C",bioresult[,2])
  
  keepind <- bioresult %>% filter((pval<pval.cutoff | score>score.cutoff) & !is.na(score) & score!="NaN" ) %>% group_by(from_clust,eval(as.symbol(colnames(bioresult)[3]))) %>% tally
  colnames(keepind)[1:2]<-colnames(bioresult)[c(1,3)]
  bioresult_keep<-bioresult %>% inner_join(keepind,by=colnames(bioresult)[c(1,3)])
  
  
  ##keep the root node to the first level
  bioresult_keep<-rbind(bioresult[bioresult$from_integrated_snn_res.==bioresult$from_integrated_snn_res.[1],],bioresult_keep[,1:6])
  
  keepv<-paste0(prefix, bioresult_keep[,4],"C",bioresult_keep[,2])
  removev<-nodesinsplit[!nodesinsplit %in% keepv]
  
  ## remove nodes that link to the removev
  dist_keep_remove<-apply(distances(graph,v=which(nodes$node %in% keepv),to=which(nodes$node %in% removev),mode="in"),1,min)
  keepv<-nodes$node[nodes$node %in% keepv][!is.finite(dist_keep_remove)]
  
  #keep the root node
  
  keepv <-unique(c(keepv,nodes$node[1]))
  distval <-apply(distances(graph,to=which(nodes$node %in% keepv),mode="out"),1,min) 
  pruned_graph<-graph %>% activate("nodes") %>% mutate(distval=distval)  %>% filter(is.finite(distval)) 
  
 
  ##generate a final clustering result based on the score and pvalues
  highind<-nodes[[prefix]]==max(resolution)
  distval<-distances(graph,v=which(highind),to=which(nodes$node %in% keepv),mode="in")
  rownames(distval)<-nodes$cluster[highind]
  
 
  ###
  id_min<-apply(distval,1,which.min)
  mapdata<-data.frame(groupid=names(id_min),from=id_min)
  mapdata<-split.data.frame(mapdata,mapdata$from)
  combinenode<-data.frame(node=nodes$node[nodes$node %in% keepv][as.numeric(names(mapdata))],newname=paste0("C",1:length(mapdata)))
  #combineresult<-list()
  biocluster<-unlist(obj[[paste0(prefix,max(resolution))]])
  levels(biocluster)<-c(levels(biocluster),paste0("C",1:length(mapdata)))
    
  for (i in 1:length(mapdata)){
 
  biocluster[biocluster %in% mapdata[[i]]$groupid] <- paste0("C",i)
  #combineresult[[i]]<-combineclusters
  #names(combineresult)[i]<-paste0("C",i)
  }    
  biocluster<-droplevels(biocluster)
  Idents(obj)<-biocluster
  
  node.size<-ifelse(is.null(node.size),as.name("size"),node.size)
  node.color<-ifelse(is.null(node.color),as.name(prefix),node.color)
  
  pruned_graph<-pruned_graph %>% activate("nodes") %>% left_join(combinenode,by="node")
  
  graphplot<-ggraph(pruned_graph,layout="tree")+geom_edge_diagonal()+
    geom_edge_diagonal(angle_calc = 'along',label_dodge = unit(-1.5, 'mm'),label_size=3,aes(filter=!is.na(score),label=paste0("p=",round(pval,2))))+
    geom_edge_diagonal(angle_calc = 'along',label_dodge = unit(1.5, 'mm'),label_size=3,aes(filter=!is.na(score),label=paste0("s=",round(score,2))))+
    scale_size(range=c(4,10))+theme(panel.background = element_blank())+
    geom_node_point(aes_(size=node.size,color=node.color))
  
  if (is.null(node.text))
  {graphplot<-graphplot+geom_node_text(aes(label=newname))} else if(node.text=="cluster")
  {graphplot<-graphplot+geom_node_text(aes(label=cluster))}
  
  
  return(list(obj=obj,ggplot=graphplot))
}


CalScore<-function(cor_bulk,cluster) {
  
  sw<-rep(0,ncol(cor_bulk))
  for (i in 1:ncol(cor_bulk)) {
    incluster<-cluster == cluster[i]
    corval<-cor_bulk[-i,i]
    incluster<-incluster[-i]
    sw[i]<-mean(corval[incluster])-mean(corval[!incluster])
    
    
  }
  
  sw_mean<-tapply(sw,cluster,mean)
  
  
  return(sw_mean)
}

