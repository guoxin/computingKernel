dyn.load("GK3WN.so")
dyn.load("GK3W.so")
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
