library(cluster)



#read data
tbl<-read.table(commandArgs()[4], head=TRUE)
#create distance matrix
d<-as.dist(tbl)
#cluster
clust<-hclust(d,method="complete")
#plot cluster dendrogramm
#plot(clust)
#rect.hclust(clust,k=4,border="red")
#rect.hclust(clust,k=6,border="green")

paint<-rainbow(17)

#for 1 to 18 extract clusters from the dendrogramm
#seems like either 4 clusters (manually selected k) or 6 clusters (height smaller 60) is reasonable
for(i in 2:17) {
#write(cutree(clust,4),"tmp.txt",ncolumns=20)
plot(clust, main=sprintf("Cluster Dendrogram (%d clusters)", i))
#rect.hclust(clust,k=4,border="red")
rect.hclust(clust,k=i,border="red")
plot(silhouette(cutree(clust,i),d),main=sprintf("Silhouette plot of %d clusters", i))
}

