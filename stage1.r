###Distributions
un_0_to_1 <- function(n){ 
  vals <- runif(n, 0, 1)
  return(vals)
}

my_ones <- function(n){ 
             vals <- rep(1, n);
             return(vals)
           }

###Makes
set.seed(3221)

case_make_matrix <- function(name, n1, n2){
  n <- n1*n2
    switch(name,
        "path"   = {
                    cat("Making  Raw_path..")
                    Ret <- make_path_adj(n, my_ones)
                    cat(".......ok\n") 
                   },
        "cycle"  = {
                    cat("Making  Raw_cycle..")
                    Ret <- make_cycle_adj(n, my_ones)
                    cat(".......ok\n")
                   },
        "star"   = {
                    cat("Making  Raw_star..")
                    Ret <- make_star_adj(n, my_ones)
        "star"   = cat(".......ok\n")
                   },
        "board" = {
                    cat("Making  Raw_board..")
                    Ret <- make_board_adj(n1, n2, my_ones)
                    cat(".......ok\n")
                   },
        "torBoard" = {
                    cat("Making  Raw_torBoard..")
                    Ret <- make_torBoard_adj(n1, n2, my_ones)
                    cat(".......ok\n")
                   }
      )
  return(Ret)
}

make_matrix <- function(n, func_dist){
  #verificar tipos#
  vals <- func_dist(n*(n-1)/2)
  A <- matrix(0, nrow = n, ncol = n, byrow = TRUE)
  cll <- upper.tri(A, diag = FALSE) 
  A[cll] = vals 
  cll <- lower.tri(A, diag = FALSE)
  A[cll] = t(A)[cll]

  return(A)
}

make_path_adj <- function(n, func_dist){
  #verificar tipos#
  A <- matrix(0, nrow = n, ncol = n)
  for(i in 1:(n - 1)){
    A[i,i + 1] <- func_dist(1)[1]
  }
  return(A)
}

make_cycle_adj <- function(n, func_dist){
  #verificar tipos#
  A <- matrix(0, nrow = n, ncol = n)
  A[n,1] <- func_dist(1)[1]
  for(i in 1:(n - 1)){
    A[i,i + 1] <- func_dist(1)[1]
  }
  return(A)
}

make_star_adj <- function(n, func_dist){
  #verificar tipos#
  A <- matrix(0, nrow = n, ncol = n)
  for(i in 1:(n - 1)){
    A[1,i + 1] <- func_dist(1)[1]
  }
  return(A)
}

make_board_adj <- function(n1, n2, func_dist){
  #verificar tipos#
  n <- n1*n2
  A <- matrix(0, nrow = n, ncol = n)
  for(i in 1:n){
    if((i %% n2) != 0){
      A[i,i + 1] <- func_dist(1)[1]
    }
    if(((i - 1) %/% n2) < n1 - 1){
      A[i,i + n2] <- func_dist(1)[1]
    }
  }
  return(A)
}

make_torBoard_adj <- function(n1, n2, func_dist){
  #verificar tipos#
  n <- n1*n2
  A <- matrix(0, nrow = n, ncol = n)
  for(i in 1:n){
    if(i %% n2 != 0){
      A[i,i + 1] <- func_dist(1)[1]
    }
    else{
      A[i, i - n2 + 1] <- func_dist(1)[1]
    }
    if(((i - 1) %/% n2) < n1 - 1){
      A[i,i + n2] <- func_dist(1)[1]
    }
    else{
      A[i,(i - 1) %% n2 + 1] <- func_dist(1)[1]
    }
  }
  return(A)
}

make_raw_data <- function(num_matrices, din_matrices, func_dist){
  #verificar tipos#
  data <- vector("list", num_matrices)
  for (i in 1:num_matrices) {
    data[[i]] <- make_matrix(din_matrices, func_dist) 
  }
  return(data)
}
