#### Stage 3: analise

## Clustering coefficient
clustering_coeficient <- function(H){
  V = nrow(H)
  exist_tri <- rep(NA, V)
  pos_tri <- rep(NA, V)
  m_v <- 0
  for(v in 1:V){
    exist_tri[v] = sum(H[H[v,] == 1, H[v,] == 1]) / 2
    m_v <- sum(H[,v])
    pos_tri[v] = m_v * (m_v - 1) / 2
  }
  clst_vect = exist_tri / pos_tri
  cat(exist_tri, "\n")
  cat(pos_tri, "\n")

  return(clst_vect)
}

## Hamming distance

hamming_dist <- function(A, B, normalized = TRUE, undirected = TRUE){
  n <- nrow(A)
  hamm_mt <- matrix(0, nrow=n, ncol=n)
  if(undirected == TRUE){
    cll <- upper.tri(A, diag = FALSE)
  }
  else{
    cll <- matrix(TRUE, nrow=n, ncol=n)
  }
  hamm_mt[cll][A[cll] != B[cll]] <- 1
  dist <- sum(hamm_mt)
  if(normalized == TRUE){
    dist <- 2 * dist / (n * (n-1))
  }
  return(dist)
}
