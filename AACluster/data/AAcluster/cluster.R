args<-commandArgs(trailingOnly = T)
args
clust <- hclust(as.dist(read.table(args, head=TRUE)), method="complete")
plot(clust)
