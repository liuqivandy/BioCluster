% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BioCluster.R
\name{BioCluster_prune}
\alias{BioCluster_prune}
\title{Pruning the clustering tree by removing those splits not supported by external dataset}
\usage{
BioCluster_prune(
  obj,
  score.cutoff = 0.2,
  pval.cutoff = 0.05,
  node.size = NULL,
  node.color = NULL,
  node.text = NULL
)
}
\arguments{
\item{obj}{Seurat object}

\item{score.cutoff}{the cutoff of score to keep the split: (default:0.2)}

\item{pval.cutoff}{the cutoff of p-value to keep the split; (default:0.05)}

\item{node.size}{node size; either a numeric value or NULL ;  if NULL, the size of each node is determined by the size of the cluster (default:NULL)}

\item{node.color}{node color; either a color or  NULL ; if NULL, the color of each node is determined by the resolution (default:NULL)}

\item{node.text}{node text; either NULL, "cluster" or others; if NULL, the node text is the new cluster id after pruning; if "cluster", node text is the cluster ID at each resolution; if others, node text is not shown(default:NULL)}
}
\value{
a list containing a Seurat object with the final clustering result after pruning and a ggplot object
}
\description{
Visualizing the pruned result
}
