#' Null model 1: permuting link weights among nodes
#'
#' @description Permutation of link weights among nodes in a edgelist
#' @param edgelist 3-column matrix object
#' @param iter Number of random edgelists to be created
#' @details Shuffles the link weights among nodes in a network
#' @return list of randomized edgelists
#' @references "Opsahl et al 2008 Physical review letters 101, 168702"; "Bascompte et al. 2005. PNAS"
#' @author Mauricio Cantor
#' @export
null1 <- function(edgelist, iter=1000){
  
  edgelist <- as.matrix(edgelist)
  el <- list()
  
  for(i in 1: iter){
    rand <- edgelist
    # shuffling link weights among pairs of nodes
    rand[,3] <- sample(edgelist[,3])
    el[[i]] <- rand
  }
  
  el
}




#' Null model 2: Randomize binary matrices based on rows and columns totals
#' @description This null model created random binary matrices by resorting the 1's among the matrix cells according to marginal totals of rows and columns.
#' @param mat Binary matrix
#' @param iter Number of iterations
#' @return mat.t array of randomized matrices
#' @details Each cell has a probability of being filled that is proportional to the number of occurrences of individuals in sites: cij = 1/2*(Pi/C + Pj/R) where Pi= row sums; Pj = column sums; C = number of columns; and R = number of rows.
#' @references "Bascompte J, Jordano P, Melian CJ and Olesen JM (2003) The nested assembly of 
#'              plant-animal mutualistic networks. PNAS, 100(16):9383-9387"
#' @author Mathias Pires with modifications by Diego Barneche and Mauricio Cantor
#' @export
null2 <- function(mat, iter=1000){
  nR<-nrow(mat);    nC<-ncol(mat)
  mR<-rowSums(mat); mC<-colSums(mat)
  
  if(identical(as.vector(mat),as.numeric(as.logical(mat)))==FALSE) stop("this function only takes binary matrices")
  
  #generating a matrix with probablities for each cell
  tR    <-  rep(mR, each=nC); tC<-rep(mC, nR)
  prob  <-  matrix(((tR/nC)+(tC/nR))/2, nR, nC, byrow=TRUE)
  #filling theoretical matrices
  mat.t=array(0,c(nR,nC,iter)) #To store theoretical matrices
  s=1
  
  while (s<=iter){
    rand=matrix(runif(nR*nC),nR,nC)
    aux=mat.t[,,s] #avoid indexing problems
    aux[which(rand<prob)]=1 #filling empty matrix
    #avoid zeroed rows or columns 
    #rows
    rm.aux=rowSums(aux)
    cols=sample(1:nC,sum(rm.aux==0), replace=TRUE) #randomly selecting columns
    for (i in 1:sum(rm.aux==0)){
      aux[which(rm.aux==0)[i],cols[i]]=1
    }
    #columns
    cm.aux=colSums(aux)
    rows=sample(1:nR,sum(cm.aux==0), replace=TRUE) #randomly selecting rows
    for (i in 1:sum(cm.aux==0)){
      aux[rows[i],which(cm.aux==0)[i]]=1
    }
    #storing matrices
    if (sum(aux)==sum(mat)){ #whenever the resulting matrix has a different number of interactions the code runs again
      mat.t[,,s]=aux  #store the matrix within the array
      s=s+1
    }
  }  
  mat.t
}




#' Null model 3: Version of \code{null 2} for quantitative matrices and weighted networks
#' @description This null model followes the same idea of null2 by Bascompte et al. 2003, but can be applied bot to binary and quantitative matrices. It uses a slightly different way of producing a probabilistic matrix for reference, calling the algorithm \code{mgen} (see help(mgen)). 
#' @param mat A binary or quantitative matrix
#' @param iter Number of iteractions, i.e. number of random matrices to be created
#' @param ... Further arguments from \code{mgen}
#' @return a list with iter random matrices
#' @author Mauricio Cantor with modifications by Diego Barneche
#' @export
null3 <- function(mat, iter, ...) {
  # creating the reference probability matrix on which the randomization process will be based
  pmat <- (rowSums(mat)) %*% t(colSums(mat))/(sum(mat))^2 
  null <- list()
  for(i in 1:iter){
    # using the 'mgen' algorithm by Vazquez et al. 2009. see below
    null[[i]] <- mgen(pmat, n=sum(mat), ...) 
  }
  null  #all random matrices in a list
}





#' Null models 4 to 7: 
#' @description shuffles XXX among XXX. It calls algorithms from package \code{vegan} to generate random matrices.
#' @param mat A quantitative matrix
#' @param iter Number of random matrices to be created 
#' @param model Function to be chosen. See \code{Details}.
#' @param ... Further arguments from \code{permatswap} or \code{permatfull}
#' @return a list with \code{iter} random matrices
#' @details If \code{model} = 4, an unrestricted null model is called. All cell values are randomly permuted. Only total total sum is constrained; row sums, column sums, and matrix fill are allow to change.
#' 
#' If \code{model} = 5, cell values are permuted restricting only column sums and total sum.
#' 
#' If \code{model} = 6, cell values are permuted restricting only row sums and total sum.
#' 
#' If \code{model} = 7, totally restricted null model is called. Cell values are permuted restricting all features of the original matrix: column sums, row sums, matrix fill and total sum.
#' @author Mauricio Cantor and Diego Barneche, using functions available in \code{vegan}
#' @references \code{citation("vegan")}
#' @seealso \code{library(vegan)}
#' @export
null4_7 <- function (mat, iter, model, ...){
  mat <- round(mat) #only works with integers...
  
  if(model %in% 4:7 == FALSE)
    stop("Available models are numbered from 4 to 7, see Details in help(null4_7)")
  
  if(model == 4){ 
    aux1 <- permatfull(mat, times=iter, fixedmar="none", shuffle="ind", mtype="count")
    return(aux1$perm)
  }
  
  if(model == 5){ 
    aux2 <- permatfull(mat, times=iter, fixedmar="column", shuffle="ind", mtype="count")
    return(aux2$perm)
  }
  
  if(model == 6){ 
    aux3 <- permatfull(mat, times=iter, fixedmar="row", shuffle="ind", mtype="count")
    return(aux3$perm)
  }
  
  if(model == 7) { 
    aux4 <- permatswap(mat, times=iter, method="quasiswap", fixedmar="both", shuffle="both", mtype="count")
    return(aux4$perm)
  }
}




#' mgen algorithm
#' @description This is a generic function to build null models based on probability matrices. 
#' @param m quantitative matrix
#' @param n total number of occurrences in the matrix (if individuals, default is the total sum of the matrix)
#' @param zs logical. See \code{Details}
#' @param rep.cell logical. See \code{Details}
#' @return a randomized matrix
#' @details \code{Mgen} is general in the sense that it allows any type of probability matrix 
#'          to be used for constructing the simulated matrices. It does not, however, constrain rows 
#'          and column totals, nor does it constrain connectance. If rep.cell=TRUE, repeated interactions
#'          are added, thus generating a quantitative matrix with cell values as positive integers. 
#'          If rep.cell=FALSE, no repeated assignment of interactions is allowed, thus generating a binary 
#'          matrix of zeros and ones. Please note that when rep.cell=FALSE the number of interactions to be 
#'          assinged must be equal or lower than the number of cells in the matrix. if \code{zs} is FALSE 
#'          it does not allow any column or row to be empty, i.e. sum up zero (default);
#'          if TRUE, it does. If \code{rep.cell} is TRUE it returns a quantitative matrix (default);
#'          if FALSE, it returns a binary matrix.
#' @references "Vazquez DP, Chacoff N and Cagnolo L (2009) Evaluating multiple determinants of the structure of mutualistic networks. Ecology, 90:2039-2046"
#' @author Diego Vazquez, with modifications by Diego Barneche
#' @export
mgen  <-  function (m, n=sum(m), zs = FALSE, rep.cell = TRUE) {
  if (rep.cell == FALSE & n > (nrow(m) * ncol(m))) {
    message("Argument n should be smaller than the number of cells in matrix")
  } else {
    mac  <- cumsum(m)
    mint <- matrix(0, nrow(m), ncol(m))
    if (zs == FALSE) {
      c1 <- sample(ncol(m), nrow(m), replace = TRUE, prob = colSums(m))
      for (i in 1:nrow(m)) {
        mint[i, c1[i]] <- 1
      }
      r1 <- sample(nrow(m), ncol(m), replace = TRUE, prob = rowSums(m))
      for (i in 1:ncol(m)) {
        mint[r1[i], i] <- 1
      }
    }
    rand  <-  runif((n-sum(mint)), min(mac), 1)
    for(i in rand){
      ri  <-  min(which(mac >= i))
      if (rep.cell == TRUE) 
        mint[ri] <- mint[ri] + 1
      if (rep.cell == FALSE) 
        mint[ri] <- 1
    }
    mint
  }
}









#' Transforms edgelist into matrix
#'
#'@param edgelist 3-column matrix with node in, node out, weight of interaction
#'@param type \code{binary} returns a binary matrix and \code{weighted} returns a weighted matrix 
#'@author Mauricio Cantor
#'@export
el_mat <- function(edgelist, type='weighted'){
  
  edgelist=as.matrix(edgelist)
  mat=matrix(0, max(edgelist[,1]),max(edgelist[,2]))
  mat.w=mat
  
  for (k in 1:nrow(edgelist)) {
    aux1=edgelist[k,1]
    aux2=edgelist[k,2]
    mat[aux1,aux2]=1 #binary
    mat.w[aux1,aux2]=edgelist[k,3] #weighted
  }
  
  if(type=="binary"){
    return(mat)
  } else {
    return(mat.w)
  }
}







#' Transforms matrix into edgelist
#'
#'@param mat matrix of interactions
#'@param type \code{binary} returns a binary edgelist and \code{weighted} returns a weighted edgelist
#'@author Mauricio Cantor 
#'@export
mat_el <- function(mat, type='weighted'){
  
  el <- reshape2::melt(mat)
  el <- el[which(el[,3]!=0),]
  el <- as.matrix(el)
  colnames(el) = rownames(el) <- NULL
  
  if(type=="binary"){
    return(el[,1:2])
  } else {
    return(el)
  }
}