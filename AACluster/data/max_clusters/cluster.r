library(cluster)



#read data
tbl<-read.table(commandArgs()[4], head=TRUE)
#create distance matrix
d<-as.dist(tbl)
#cluster
clust<-hclust(d,method="complete")
#plot cluster dendrogramm
plot(clust)
rect.hclust(clust,k=6,border="green")
rect.hclust(clust,k=4,border="red")

#for 1 to 18 extract clusters from the dendrogramm
#seems like either 4 clusters (manually selected k) or 6 clusters (height smaller 60) is reasonable
#for(i in 3:8) {
#write(cutree(clust,4),"tmp.txt",ncolumns=20)
plot(silhouette(cutree(clust,4),d))
plot(silhouette(cutree(clust,6),d))
#}

