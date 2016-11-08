#Funcao de ajuste: recebe o limiar_base(p), o valor de ajuste(eta) e um valor
#de normalizacao(k)
source("myutils.r")


limiarize_ctl <- function(A, maxA, p){
  #verificar tipos#
  limiar <- maxA * p
  B <- matrix(0, ncol = ncol(A), nrow = nrow(A))
  cll <- upper.tri(A, diag = FALSE) 
  B[cll][A[cll] <= limiar] = 1
  cll <- lower.tri(B, diag = FALSE)
  B[cll] = t(B)[cll]
  return(B)
}

##Tentativa 4
# Cuidado! Pode retornar Inf 
normVal_calc <- function(H, alpha){
  cll <- c(1:ncol(H))
  n <- nrow(H)
  temp <- 1/H[,cll]
  beta <- sum(temp[upper.tri(temp, diag = FALSE)])
  beta <- beta*2 / (n*(n-1))
  k <- 1 / (alpha + beta) 
  if(k == Inf){
    cat("NORMVAL_CALC (ERROR) : Calculated k = Inf!\n")
  }
  return(k)
}

my_create_P <- function(p, H, k, alpha){
  cll <- c(1:ncol(H))
  P <- ((1 / H[,cll]) + alpha) * k * p
  return(P)
}

limiarize <- function(A, max_A, p, H, create_P, alpha){
  #verificar tipos#
  B <- matrix(NA, ncol = ncol(A), nrow = nrow(A))
  k <- normVal_calc(H, alpha)
  if(k == Inf){
    print("LIMIARIZE (ERROR) : Received k = Inf. B with NA's returned!")
  }
  else{
    p_max <- 1 / ((alpha + (1/min(H[upper.tri(H, diag = FALSE)]))) * k)
    myprint(paste0("p_max calculated ", p_max, "\n"))
    if( p >= p_max){
      cat("WARNING: p >= 1 / ((alpha + (1/min(H)))*k)\n")
    }
    myprint(paste0("k = ",k, "\n"))
    P <- create_P(p * max_A, H, k, alpha)
    B <- matrix(0, ncol = ncol(A), nrow = nrow(A))
    cll <- upper.tri(A, diag = FALSE) 
    B[cll][A[cll] >= P[cll]] = 1
    cll <- lower.tri(B, diag = FALSE)
    B[cll] = t(B)[cll]
    myprint("B ............ calculated\n")
  } 
  return(B)
}

##Tentativa 5
#Exemple of pertubation function f: R^(nxn) --> R^(nxn)
pertubation1 <-function(H){
  alpha <- 2
  M <- H**alpha
  return(M)
}

normVal_calc5 <- function(M, n){
  sumF <- sum(M[upper.tri(M, diag = FALSE)])
  k <- n*(n-1)/(2*sumF)
  return(k)
}

adjust <- function(k, p, M){
  P <- k*p*M
  return(P)
}

limiarize5 <- function(A, maxA, P){
  B <- matrix(0, ncol = ncol(A), nrow = nrow(A))
  cll <- upper.tri(A, diag = FALSE) 
  B[cll][A[cll] <= (P[cll]*maxA)] = 1
  cll <- lower.tri(B, diag = FALSE)
  B[cll] = t(B)[cll]
  myprint("B ............ calculated\n")
  return(B)
}

limiarization <- function(A, maxA, p, H, pertubation){
  M <- pertubation(H)
  k <- normVal_calc5(M, ncol(A))
  P <- adjust(k, p, M)

  ### testing p ###
  tri <- M[upper.tri(M, diag = FALSE)]
  n <- ncol(A)
  p_max <- 2 * sum(tri) / (max(tri)*n*(n-1))
  myprint(paste0("   p_max= ",p_max, "\n"))
  if( p >= p_max){
    cat("WARNING: p >= p_max! p = ",p," e p_max = ", p_max, "\n")
  }
  ###

  B <- limiarize5(A, maxA, P)
  return(B)
}

