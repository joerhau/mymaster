args<-commandArgs(trailingOnly = T)
args

rand <- read.delim("RAND.csv")
plot(rand[,],type="l")

