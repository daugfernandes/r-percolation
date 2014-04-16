##############################################################################
# lattice constructor
#  size - side of the squared space
#  p - probability of a site be occupied
#-----------------------------------------------------------------------------
percolation <- function (size, p) {
  return (list("size"=size, 
               "p"=p, 
               "lattice"=ifelse(matrix(runif(size^2,0,1),nrow=size)<p,1,0)))
}

##############################################################################
# find cluster's information
#  space - 
#-----------------------------------------------------------------------------
percolation.clusters <- function(space) {
  
  size = length(space$lattice[,1])
  
  # recursive cluster detection
  set <- function(ol, cl, col, row, value) {

    ol[col, row] <- FALSE
    cl[col, row] <- value
    
    if(col>1) {
      if(ol[col-1, row]==1 && cl[col-1, row]==0)
      {
        ret <- set(ol, cl, col-1, row, value)
        ol <- ret$ol
        cl <- ret$cl
      }
    }
    
    if(col<size) {
      if(ol[col+1, row]==1 && cl[col+1, row]==0)
      {
        ret <- set(ol, cl, col+1, row, value)
        ol <- ret$ol
        cl <- ret$cl
      }
    }
    
    if(row>1) {
      if(ol[col, row-1]==1 && cl[col, row-1]==0)
      {
        ret <- set(ol, cl, col, row-1, value)
        ol <- ret$ol
        cl <- ret$cl
      }
    }
    
    if(row<size) {
      if(ol[col, row+1]==1 && cl[col, row+1]==0)
      {
        ret <- set(ol, cl, col, row+1, value)
        ol <- ret$ol
        cl <- ret$cl
      }
    }

    return(list("ol"=ol,"cl"=cl))
    
  }
  
  originalLattice <- space$lattice
  clustersLattice <- matrix(replicate(space$size^2, 0), nrow=space$size)
  clusterNumber <- 0

  repeat{

    if(length(which(originalLattice==1))==0) break;
    
    # first occurence of an occupied site (if there is one)
    rc <- which(originalLattice==1, arr.ind=TRUE)[1,]
    clusterNumber <- clusterNumber + 1
    ret = set(originalLattice, clustersLattice, rc["row"], rc["col"], clusterNumber)
    originalLattice <- ret$ol
    clustersLattice <- ret$cl
    if(sum(ret$ol)==0) break;
  }
  
  # counting
  summary <- table(clustersLattice)

  # indexes of meaningfull clusters
  indexes <- which(summary[2:length(summary)] > 1, arr.ind = T)
  
  # find spanningClusters
  firstRow <- clustersLattice[1,]
  lastRow <- clustersLattice[size,]
  firstCol <- clustersLattice[,1]
  lastCol <- clustersLattice[,size]

  # horizontal spanning clusters  --
  hsc <- unique(firstCol[which(firstCol[firstCol %in% lastCol]>0)])
  
  # vertical spanning clusters    |
  vsc <- unique(firstRow[which(firstRow[firstRow %in% lastRow]>0)])
  
  h <- summary[which(firstCol[firstCol %in% lastCol]>0)]
  v <- summary[which(firstRow[firstRow %in% lastRow]>0)]
  w <- summary[which(firstRow[hsc %in% vsc]>0)]
  
  return(list("lattice" = clustersLattice, 
              "size" = size,
              "clusters" = ifelse(length(summary)==1,NA,unname(summary[indexes+1])),
              "h.spanning.clusters" = unname(h[which(!is.na(h))]),
              "v.spanning.clusters" = unname(v[which(!is.na(v))]),
              "w.spanning.clusters" = unname(w[which(!is.na(w))])))
  
}

##############################################################################
# general plot function
#  space - space or clusters information
#-----------------------------------------------------------------------------
percolation.plot <- function(space) 
{
  
  # is a percolation
  if(!is.null(space$p)) 
  {
    plot(x= c(0,space$size*10), y=c(0,space$size*10), xlab="", ylab="", xaxt='n',yaxt='n')

    # build graphic representation
    foo <- mapply(function (x)
      mapply(function(y) 
        rect(xleft = (y - 1) * 10,
             ybottom = (x - 1) * 10,
             xright = y * 10,
             ytop = x * 10,
             col = ifelse(space$lattice[x, y], 2, 0)),
             1:space$size), 
           1:space$size)
  } 
  else 
  {
    lattice <- space$lattice
    
    # is a cluster representation
    plot(x = c(0,length(lattice[1,]) * 10), 
         y = c(0,length(lattice[1,]) * 10), 
         xlab = "", 
         ylab = "", 
         xaxt = 'n', 
         yaxt = 'n')

    # array frequency of sites for each cluster number
    # 0  1  2  3  4  5  6  7  8  9 10 11   <- cluster number 
    #64  3  9  6  4  4  1  4  2  1  1  1   <- number of sites
    summary <- table(lattice)

    # indexes of meaningfull clusters
    indexes <- which(summary[2:length(summary)] > 1, arr.ind = T)
    
    # indexes of single-site "clusters"
    singleSites <- which(summary[2:length(summary)] == 1, arr.ind = T)

    # build graphic representation
    foo <- mapply(function (x)
      mapply(function(y) 
        rect(xleft = (y - 1) * 10,
             ybottom = (x - 1) * 10,
             xright = y * 10,
             ytop = x * 10,
             # color of the square
             col=ifelse(lattice[x, y] %in% indexes, 
                        lattice[x, y] + 1,
                        # if we need to color single-site cluster's
                        ifelse(lattice[x, y] %in% singleSites,
                               0,
                               0))),
             1:length(lattice[1,])), 
                1:length(lattice[1,]))
  }
}
