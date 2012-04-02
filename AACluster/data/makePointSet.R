max <- as.dist(read.table("max.txt", head=TRUE))

# Generate 10 random points in 2 dimensions 
nPoints = 20; 
nDim = 2; 

set.seed(10); 
points = matrix(runif(nPoints * nDim), nPoints, nDim); 

# Their distance: 
dst = dist(points) 

# Classical multidimensional scaling 
mds = cmdscale(dst); 

# Distance of the points calculated by mds 
dst2 = dist(mds); 

# The two distances are equal 
all.equal(as.vector(dst), as.vector(dst2)) 
