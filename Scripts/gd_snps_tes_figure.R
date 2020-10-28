### Making the GD plot of SNPs vs TEs for Christine's paper ###
## Clearing the global environment
rm(list=ls(all=TRUE))

## Setting up the working directory
getwd()
setwd("/Users/jonathanrenk/Desktop/")

## Loading in the data
data <- read.table("te_dist.txt", header=F, sep = "\t")
str(data)
## Removing taxa name column
data <- data[,c(-1)]
colnames(data) <- NULL # removing column headers
data <- as.matrix(data) # converting to matrix 
data[upper.tri(data)] <- NA # make the upper half of the matrix as NAs
diag(data)=NA # make the diagonal NAs since they are zero
data <- cbind(which(!is.na(data), arr.ind = TRUE), na.omit(as.vector(data))) # getting into a vector of values for plotting
data <- as.data.frame(data)

## Loading in the data
data2 <- read.table("snp_dist.txt", header=F, sep = "\t")
str(data2)
## Removing taxa name column
data2 <- data2[,c(-1)]
colnames(data2) <- NULL # removing column headers
data2 <- as.matrix(data2) # converting to matrix 
data2[upper.tri(data2)] <- NA # make the upper half of the matrix as NAs
diag(data2)=NA # make the diagonal NAs since they are zero
data2 <- cbind(which(!is.na(data2), arr.ind = TRUE), na.omit(as.vector(data2))) # getting into a vector of values for plotting
data2 <- as.data.frame(data2)

### Making a new combined data frame ###
combined_data <- data[,c(1:3)]
combined_data$SNP <- data2[,c(3)]
names(combined_data)[3] <- "TE"
str(combined_data)

plot(combined_data$SNP, combined_data$TE, xlab = "SNP Genetic Distance",ylab = "TE Genetic Distance", col=rgb(0.52,0.52,0.52,0.6), pch=16 , cex=1.3)
