dyn.load("GK3WN.so")
dyn.load("GK3W.so")
dyn.load("K3.so")
computingKernel <- function(string.set, K1, beta, ni = "GK3WN"){
  len <- length(string.set)
  K3 <- .C("GK3WN", 
    peptides = as.character(string.set), 
    Length = as.integer(len), 
    K3 = as.double(matrix(0.0, nrow = len, ncol = len)), 
    K1 = K1 ^ beta,
    weights = as.double(rep(1, max(nchar(string.set)))), 
    DUP = TRUE, #because of the characters
    PACKAGE = ni
   )$K3
   dim(K3) <- c(len, len);gc()
   colnames(K3) <- rownames(K3) <- names(string.set)
   K3
}
computingKernelLocal <- function(str1, str2, K1, beta, mi = "pair", normalize = TRUE){
  #mode: "pair" and "rectangle". Let "*" denotes the kernel values.
  #in pair mode, compute kernel value in this way:
  #	str1[1]  str2[1]  *
  #	str1[2]  str2[2]  *
  #	str1[3]  str2[3]  *
  #	...................
  #	str1[n]  str2[n]  *
  #where n = min(length(str1), length(str2))
  #in rectangle mode, compute kernel value in this way:
  #		str2[1]  str2[2]  str2[3]  .......  str2[len2]
  #str1[1] 	   *        *        *     .......       *
  #str1[2] 	   *        *        *     .......       *
  #............................................................
  #str1[len1]	   *        *        *     .......       *
  #where len1 = length(str1), len2 = length(str2)
  mode <- 0
  if(!normalize) mode <- mode + 1
  len1 <- length(str1)
  len2 <- length(str2)
  maxNchar <- max(c(nchar(str1), nchar(str2)))
  if(mi == "pair"){
    mode <- mode + 6
    minLen <- min(len1, len2)
    K3 <- .C("K3",
      peptides1 = as.character(str1), 
      peptides2 = as.character(str2), 
      Length1 = as.integer(minLen),
      Length2 = as.integer(minLen),
      maxNchar = maxNchar,
      K3 = double(length = minLen), 
      K1 = K1 ^ beta, #the beta must be handled in R!!!
      weights = as.double(rep(1, maxNchar)),
      computing.mode = as.integer(mode),
      DUP = TRUE, #because of the characters
      NAOK = FALSE,
      PACKAGE = "K3")$K3; gc()
  }else if(mi == "rectangle"){
    mode <- mode + 4
    K3 <- .C("K3",
      peptides1 = as.character(str1), 
      peptides2 = as.character(str2), 
      Length1 = as.integer(len1),
      Length2 = as.integer(len2),
      maxNchar = maxNchar,
      K3 = double(length = len1 * len2), 
      K1 = K1 ^ beta, #the beta must be handled in R!!!
      weights = as.double(rep(1, maxNchar)),
      computing.mode = as.integer(mode),
      DUP = TRUE, #because of the characters
      NAOK = FALSE,
      PACKAGE = "K3")$K3; gc()
    dim(K3) <- c(len1, len2)
    rownames(K3) <- str1
    colnames(K3) <- str2
  }else{
    warning("unknown mode")
  }
  K3
}
