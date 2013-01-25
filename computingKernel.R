getAAK3 <- function(proteins, maxLen, AAK1){
  #input: proteins as a list, each entry is a list, containing:
  #  proteins[[i]]$aa, a string of amino acid sequence
  if(!require("multicore")){
  	mclapply <- lapply
	print("do not worry, it is 100% ok to run without multicore!")
  }
  dyn.load("STR.GEN.K3.so")
  maxLen <- as.integer(maxLen)
  AAK1 <- as.double(AAK1)
  N <- length(proteins)
  K3 <- matrix(0, nrow = N, ncol = N)
  preKernel <- mclapply(1 : N, FUN = function(i){
    g <- proteins[[i]]
    glen <- as.integer(nchar(g$aa))
    unlist(lapply(proteins[i:N], FUN = function(f){
    .C("AaGenK3",
      fa = as.character(f$aa),
      ga = as.character(g$aa),
      flen = as.integer(nchar(f$aa)),
      glen = glen,
      maxLen = maxLen,
      AAK1 = AAK1,
      K3 = as.double(0),
      DUP = T,
      PACKAGE = "STR.GEN.K3"
    )$K3
    }))
  })
  for(i in 1:N) K3[i:N, i] <- K3[i, i:N] <- preKernel[[i]]
  rownames(K3) <- colnames(K3) <- names(proteins)
  K3
}
nmkernel <- function(K){
	if (dim(K)[1]==1) return(1)
	dgn <- 1/sqrt(diag(K))
	K <- sweep(x = K, MARGIN = 1, STATS = dgn, FUN = "*")
	K <- sweep(x = K, MARGIN = 2, STATS = dgn, FUN = "*")
	K
}
computingKernel <- function(string.set, K1, beta, ni = "GK3WN", maxLen = NULL){
	stopifnot(ni %in% c("GK3WN", "GK3W"))
	stopifnot(sum(!(unique(unlist(strsplit(string.set,
		split = "", fixed = T))) %in% rownames(K1))) == 0)
	proteins <- list()
	if(is.null(names(string.set))) names(string.set) <- string.set
	for(a in names(string.set)) proteins[[a]] <- list(aa = string.set[a])
	if(is.null(maxLen)) maxLen <- max(nchar(string.set))
	if(ni == "GK3W") return(getAAK3(proteins, maxLen, AAK1 = K1^beta))
	return(nmkernel(getAAK3(proteins, maxLen, AAK1 = K1^beta)))
}
