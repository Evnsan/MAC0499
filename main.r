#!/usr/bin/env Rscript
source("stage1.r")
source("stage2.r")
source("stage3.r")
source("myutils.r")
library("igraph")

###images path
#img_path <- "~/MAC0499/R_module/img"
img_path <- "~/git/MAC0499/R_module/img"

if(!file.exists(img_path)){
  dir.create(img_path)
}
###
###results path
#results_path <- "~/MAC0499/R_module/results"
results_path <- "~/git/MAC0499/R_module/results"

if(!file.exists(results_path)){
  dir.create(results_path)
}

### Script settings
## times to run the script
tries <- 100 
## number of vertices
#n1 <- 80 
#n2 <- 80 
n1 <- 15 
n2 <- 15
n <- n1 * n2
## base limiar
p_stages <- c(0.25, 0.5, 0.75)
## alpha => adjust power
alpha <- 0.2 
##pertubation function
my_pertubation <- function(H){
  M <- H**alpha
  return(M)
}
## settings for matrix A 
A_unif_min <- 0
A_unif_max <- 1 
my_dist <- function(n){ 
             vals <- runif(n, A_unif_min, A_unif_max);
             return(vals)
           }
## base graphs
Raw_names = c("path", "cycle", "star", "board", "torBoard")
###

### Flux controls
## obs. to see matrices n must be less than or equal 316
#see_matrices <- TRUE 
see_matrices <- FALSE 
#plot_G <- TRUE
plot_G <- FALSE
#save_results <- TRUE
save_results <- FALSE 
#show_results <- TRUE
show_results <- FALSE

###

### Internals

## control results adjacence matrices
ctl_adj_mt <- array(NA, c(n, n, length(p_stages)) )

## result matrices
result_mt <- list()
result_mt_names <- c("DEGREE AVERAGE",
                     "CLOSENESS",
                     "CORENESS",
                     "CLUSTERING COEFFICIENT",
                     "TRANSITIVITY",
                     "HAMMING DISTANCE",
                     "MODULARITY IN COMMUNITIES BY CLUSTER FAST GREEDY",
                     "VAR OF SIZES IN COMMUNITIES BY CLUSTER FAST GREEDY"
                    )
stats_number <- length(result_mt_names)

label <- rep("",stats_number)
for(result_name_index in 1:stats_number){
    label[result_name_index] <- paste0("[[", result_name_index,"]]")
}
names(result_mt_names) <- label 

result_mt_desc <- c("sum(deg(v)) / n for all v in V",

                    "n / sum(d(v,i)), for all i not v, in V",

                    "largest k-core for every vertex",

                    "the ratio of the triangles connected
                     to the vertex and the triples centered
                     on the vertex",

                    "the ratio of the triangles and the
                     connected triples in the graph",
                    
                    "the number of edges in one graph but
                     not in the other (the base to compare
                     is the control graph)",

                    "modularity in communities of graphs
                     via directly optimizing a modularity
                     score",

                     "square root of the sum of sizes in
                      communities of graphs computed via
                      directly optimizing a modularity
                      score"
                   )
## mean results
result_mt_mean <- array(NA, c(
                              length(Raw_names) + 1,
                              length(p_stages),
                              tries,
                              2*stats_number
                             )
                       )

tries_names <- rep("", tries)
for(try in 1:tries){
    tries_names[try] <- paste0("Try ", try)
}

rate_names <- rep("",stats_number)
for(name in 1:stats_number){
    rate_names[name] <- paste0(result_mt_names[name], " (RATE)")
}
mean_mt_names <- c(result_mt_names, rate_names)

dimnames(result_mt_mean) <- list(  c(c("control"), Raw_names),
                                   p_stages,
                                   tries_names,
                                   mean_mt_names
                                )
## degreee avarage
deg <- matrix(NA, ncol = length(p_stages), nrow = length(Raw_names) + 1)
rownames(deg) <- c(c("control"), Raw_names)
colnames(deg) <- p_stages
deg_ctl <- rep(NA,length(p_stages))
result_mt[[1]] <- deg
##closeness
closeness_values <- array(NA, c(length(Raw_names) + 1, length(p_stages), 2))
dimnames(closeness_values) <- list(c(c("control"), Raw_names),
                                   p_stages,
                                   c("mean" , "standard deviation")
                                  )
result_mt[[2]] <- closeness_values
##coreness
coreness_values <- array(NA, c(length(Raw_names) + 1, length(p_stages), 2))
dimnames(coreness_values) <- list(c(c("control"), Raw_names),
                                   p_stages,
                                   c("mean" , "standard deviation")
                                  )
result_mt[[3]] <- coreness_values 
###
##clustering_coefficient
ctl_coeff_values <- array(NA, c(length(Raw_names) + 1, length(p_stages), 2))
dimnames(ctl_coeff_values) <- list(c(c("control"), Raw_names),
                                   p_stages,
                                   c("mean" , "standard deviation")
                                  )
result_mt[[4]] <- ctl_coeff_values 
###
##transitivity
transitivity_values <- matrix(NA, ncol = length(p_stages),
                               nrow = length(Raw_names) + 1
                             )
rownames(transitivity_values) <- c(c("control"), Raw_names)
colnames(transitivity_values) <- p_stages
result_mt[[5]] <- transitivity_values 
###
##Hamming distance
hamm_dist <- matrix(NA, ncol = length(p_stages),
                               nrow = length(Raw_names) + 1
                             )
rownames(hamm_dist) <- c(c("control"), Raw_names)
colnames(hamm_dist) <- p_stages
result_mt[[6]] <- hamm_dist
###
##Cluster fast and greed
modularity <- matrix(NA, ncol = length(p_stages),
                               nrow = length(Raw_names) + 1
                             )
rownames(modularity) <- c(c("control"), Raw_names)
colnames(modularity) <- p_stages
result_mt[[7]] <- modularity 
###
##Cluster fast and greed
var_sizes <- matrix(NA, ncol = length(p_stages),
                               nrow = length(Raw_names) + 1
                             )
rownames(var_sizes) <- c(c("control"), Raw_names)
colnames(var_sizes) <- p_stages
result_mt[[8]] <- var_sizes 
###

### Printing parameters
cat("Parameters:\n")
cat("tries = ", tries, "\n")
cat("n1 = ", n1, "\n")
cat("n2 = ", n2, "\n")
cat("n = n1 * n2 = ", n, "\n")
cat("alpha = ", alpha, "\n")
cat("p_satges = ", p_stages, "\n")
cat("p_list = ", p_stages * A_unif_max, "\n")
cat("\nlimits for entries of A:\n")
cat("min = ", A_unif_min, "\n")
cat("max = ", A_unif_max, "\n")
###

###Script beginning##
cat("\n\nStarting script...\n\n")

for(try in 1:tries){
  cat("\n\n################################################################\n")
  cat("TRY ", try, ":\n")
  #creating A#
  myprint("\n\n#########################################\n")
  myprint("Making  A..")
  A <- make_matrix(n, my_dist)
  myprint("..............ok\n")
  if(see_matrices){
    cat("########### A ##########\n")
    print(A)
    cat("########################\n")
  }
  myprint("\n\n")
  myprint("  ------------------------------  \n")
  myprint("   BCONTROL\n")
  myprint("  ------------------------------  \n")

  for(p_index in 1:length(p_stages)){
    myprint("=========================\n")
    myprint(paste0("   using p = ", p_stages[p_index],"\n"))

    ctl_adj_mt[,,p_index] <- limiarize_ctl(A,A_unif_max, p_stages[p_index])

    if(see_matrices){
      cat("----------- B_CONTROL ----------\n")
      print(ctl_adj_mt[,,p_index])
      cat("--------------------------------\n")
    }
    sum_B_ctl = sum(ctl_adj_mt[,,p_index])
    deg_ctl = sum_B_ctl / n
    myprint(paste0("Average_degree = ", deg_ctl[p_index], "\n"))
    myprint(paste0("Density = ", sum_B_ctl / (n * (n-1)), "\n"))

    ### graph analise ###
    G <- graph_from_adjacency_matrix(ctl_adj_mt[,,p_index], mode = c("undirected"),
                                    weighted = NULL 
                                  )
    clos_vect <- closeness(G, v = V(G),mode = c("all"), weights = NULL,
                            normalized = TRUE 
                          )
    core_vect <- coreness(G, mode = c("all"))

    ctl_coeff_vect <- transitivity(G, type = c("localundirected"), vids = NULL,
                                   weights = NULL, isolates = c("NaN")
                                  )

    community <- cluster_fast_greedy(G, weights = NULL)

    ## degree_average ##
    result_mt[[1]][1,p_index] <- deg_ctl 

    result_mt_mean[1,p_index,try,1] <- result_mt[[1]][1,p_index] 
    result_mt_mean[1,p_index,try, stats_number + 1] <- 1 
    ## closeness ##
    result_mt[[2]][1,p_index,1] <- mean(clos_vect, na.rm = TRUE)
    result_mt[[2]][1,p_index,2] <- sd(clos_vect, na.rm = TRUE)
    
    result_mt_mean[1,p_index,try,2] <- result_mt[[2]][1,p_index,1] 
    result_mt_mean[1,p_index,try, stats_number + 2] <- 1 
    ## coreness ##
    result_mt[[3]][1,p_index,1] <- mean(core_vect, na.rm = TRUE)
    result_mt[[3]][1,p_index,2] <- sd(core_vect, na.rm = TRUE)
    
    result_mt_mean[1,p_index,try,3] <- result_mt[[3]][1,p_index,1] 
    result_mt_mean[1,p_index,try, stats_number + 3] <- 1 
    ##clustering_coefficient
    result_mt[[4]][1,p_index,1] <- mean(ctl_coeff_vect, na.rm = TRUE)
    result_mt[[4]][1,p_index,2] <- sd(ctl_coeff_vect, na.rm = TRUE)
    
    result_mt_mean[1,p_index,try,4] <- result_mt[[4]][1,p_index,1] 
    result_mt_mean[1,p_index,try, stats_number + 4] <- 1 
    ##transitivity
    result_mt[[5]][1,p_index] <- transitivity(G, type = c("globalundirected"),
                                              vids = NULL,
                                              weights = NULL,
                                              isolates = c("NaN")
                                             ) 
    result_mt_mean[1,p_index,try,5] <- result_mt[[5]][1,p_index] 
    result_mt_mean[1,p_index,try, stats_number + 5] <- 1 
    ##hamming distance
    result_mt[[6]][1,p_index] <- 0
    result_mt_mean[1,p_index, try, 6] <- 0
    result_mt_mean[1,p_index,try, stats_number + 6] <- 1 
    ##modularity in communities  by cluster_fg
    result_mt[[7]][1,p_index] <- modularity(community) 
    result_mt_mean[1,p_index, try, 7] <- result_mt[[7]][1,p_index] 
    result_mt_mean[1,p_index,try, stats_number + 7] <- 1 
    ##var of sizes in communities  by cluster_fg
    communities_sizes <- sizes(community)
    result_mt[[8]][1,p_index] <- ( sqrt(sum(communities_sizes^2))
                                   /sum(communities_sizes)
                                 )
    result_mt_mean[1,p_index, try, 8] <- result_mt[[8]][1,p_index] 
    result_mt_mean[1,p_index,try, stats_number + 8] <- 1 
    ###################


    ###### memory cleaning ######
    gc()
    #############################
  }
  
  #looping over H's#
  myprint("\n\n#########################################\n\n\n")

  for(name_index in 1:length(Raw_names)){
    myprint("\n\n#########################################\n")
    myprint("  ------------------------------  \n")
    myprint(paste0("   B", "_", Raw_names[name_index],"\n", sep=""))
    myprint("  ------------------------------  \n")
    
    Raw <- case_make_matrix(Raw_names[name_index],n1,n2)

    myprint(paste0("Making  G", "_", Raw_names[name_index],"..", sep=""))
    G <- graph_from_adjacency_matrix(Raw, mode = c("undirected"),
                                    weighted = TRUE
                                  )
    myprint(".........ok\n")

    if(plot_G){
      myprint(paste0("Ploting G", "_", Raw_names[name_index],"..", sep=""))
      png(filename = file.path(img_path,
                               paste(Raw_names[name_index], ".png", sep = "")
                              ),
          width=1200,
          height=700
         )
      plot(G, edge.label = round(E(G)$weight, 3))
      invisible(dev.off())
      myprint(".........ok\n")
    }

    myprint(paste0("Making  H", "_", Raw_names[name_index],"..", sep=""))
    H <- distances(G, v = V(G), to = V(G), mode = c("all"),
      weights = E(G)$weight, algorithm = c("bellman-ford")
    )
    myprint(".........ok\n")

    if(see_matrices){
      cat("\n")
      cat("########### H", "_", Raw_names[name_index]," ##########\n", sep="")
      print(H)
      cat("########################\n")
      cat("\n")
    }
    
    ###### memory cleaning ######
    rm(Raw)
    rm(G)
    gc()
    #############################
    
    myprint("\n")
    myprint("####### ADJUSTING ########\n")
    for( p_index in 1:length(p_stages)){
    
      myprint("=========================\n")
      myprint(paste0("   using p = ", p_stages[p_index],"\n"))
      B <- limiarization(A, A_unif_max, p_stages[p_index], H, my_pertubation)
      if(see_matrices){
        myprint(
          paste0(
            "----------- B",
            "_",
            Raw_names[name_index]," ----------\n",
            sep=""
          )
        )
        print(B)
        myprint("------------------------\n")
      }
      G <- graph_from_adjacency_matrix(B, mode = c("undirected"),
        weighted = NULL 
      )
      ### data collect ###
      sum_B_adj = sum(B)
      deg_adj = sum_B_adj / n
      myprint(paste0("Average degree = ", deg_adj, "\n"))
      myprint(paste0("Density = ", sum_B_adj / (n * (n-1)), "\n"))
      clos_vect <- closeness( G, v = V(G),mode = c("all"), weights = NULL,
        normalized = TRUE 
      )
      core_vect <- coreness(G, mode = c("all"))
      ctl_coeff_vect <- transitivity(G, type = c("localundirected"),
        vids = NULL, weights = NULL,
        isolates = c("zero")
      ) 
      community <- cluster_fast_greedy(G, weights = NULL)
      ## degree_average ##
      result_mt[[1]][name_index + 1, p_index] <- deg_adj

      result_mt_mean[name_index + 1,p_index,try,1] <- (
        result_mt[[1]][name_index + 1,p_index]
      ) 
      result_mt_mean[name_index + 1,p_index,try, stats_number + 1] <- (
        result_mt[[1]][name_index + 1,p_index] /
        result_mt[[1]][1,p_index] 
      )
      ## closeness ##
      result_mt[[2]][name_index + 1,p_index,1] <- mean(clos_vect, na.rm = TRUE)
      result_mt[[2]][name_index + 1,p_index,2] <- sd(clos_vect, na.rm = TRUE)

      result_mt_mean[name_index + 1,p_index,try,2] <- (
        result_mt[[2]][name_index + 1,p_index,1]
      ) 
      result_mt_mean[name_index + 1,p_index,try, stats_number + 2] <- (
        result_mt[[2]][name_index + 1,p_index,1] /
        result_mt[[2]][1,p_index,1] 
      )
      ## coreness ##
      result_mt[[3]][name_index + 1,p_index,1] <- mean(core_vect, na.rm = TRUE)
      result_mt[[3]][name_index + 1,p_index,2] <- sd(core_vect, na.rm = TRUE)

      result_mt_mean[name_index + 1,p_index,try,3] <- (
        result_mt[[3]][name_index + 1,p_index,1]
      ) 
      result_mt_mean[name_index + 1,p_index,try, stats_number + 3] <- (
        result_mt[[3]][name_index + 1,p_index,1] /
        result_mt[[3]][1,p_index,1] 
      )
      ##clustering_coefficient
      result_mt[[4]][name_index + 1,p_index,1] <- mean(ctl_coeff_vect,
        na.rm = TRUE
      )
      result_mt[[4]][name_index + 1,p_index,2] <- sd(ctl_coeff_vect,
        na.rm = TRUE
      )

      result_mt_mean[name_index + 1,p_index,try,4] <- (
        result_mt[[4]][name_index + 1,p_index,1]
      ) 
      result_mt_mean[name_index + 1,p_index,try, stats_number + 4] <- (
        result_mt[[4]][name_index + 1,p_index,1] /
        result_mt[[4]][1,p_index,1] 
      )
      ##transitivity
      result_mt[[5]][name_index + 1,p_index] <- transitivity(
        G,
        type = c("globalundirected"),
        vids = NULL,
        weights = NULL,
        isolates = c("zero")
      )

      result_mt_mean[name_index + 1,p_index,try,5] <- (
        result_mt[[5]][name_index + 1,p_index]
      ) 
      result_mt_mean[name_index + 1,p_index,try, stats_number + 5] <- (
        result_mt[[5]][name_index + 1,p_index] /
        result_mt[[5]][1,p_index] 
      )
      ##hamming distance
      result_mt[[6]][name_index + 1,p_index] <- hamming_dist(
        B, ctl_adj_mt[,,p_index], normalized=FALSE
      )
      result_mt_mean[name_index +1,p_index, try, 6] <- (
        result_mt[[6]][name_index + 1, p_index]
      )
      result_mt_mean[name_index + 1,p_index,try, stats_number + 6] <- (
        result_mt[[6]][name_index + 1,p_index] /
        result_mt[[6]][1,p_index] 
      )
      ##modularity in communities  by cluster_fg
      result_mt[[7]][name_index + 1,p_index] <- modularity(community) 
      result_mt_mean[name_index +1,p_index, try, 7] <-( 
        result_mt[[7]][name_index + 1,p_index]
      ) 
      result_mt_mean[name_index + 1,p_index,try, stats_number + 7] <- (
        result_mt[[7]][name_index + 1,p_index] /
        result_mt[[7]][1,p_index] 
      )
      ##var of sizes in communities  by cluster_fg
      communities_sizes <- sizes(community)
      result_mt[[8]][name_index + 1,p_index] <- (
        sqrt(sum(communities_sizes^2))
        /sum(communities_sizes)
      )
      result_mt_mean[name_index + 1,p_index, try, 8] <- (
        result_mt[[8]][name_index + 1,p_index]
      )
      result_mt_mean[name_index + 1,p_index,try, stats_number + 8] <- (
        result_mt[[8]][name_index + 1,p_index] /
        result_mt[[8]][1,p_index] 
      )
      ####################
      
      myprint("=========================\n")
    
      ###### memory cleaning ######
      rm(B)
      rm(G)
      gc()
      #############################
    }
    myprint("##########################\n")
    myprint("\n#########################################\n")
  
    ###### memory cleaning ######
    rm(H)
    gc()
    #############################
   
  }  
  
  if(show_results){
    cat("\n\n\n############# SHOWING RESULTS ###########\n")
    cat("TRY", try, "\n")
  
    for(i in 1:stats_number){
      cat("\n")
      cat(result_mt_names[i], "MATRIX:\n")
      cat("Contains values of (",
          result_mt_desc[i],
          ")\nfor different H's and values of p\n", sep = "")
      cat("\ncols => percentage to calculate a value for p\n")
      cat("rows => base graph for H\n\n")
      print(result_mt[[i]])
      cat("-------------------------------------\n")
    }
    cat("#########################################\n")
  }
  if(save_results){
    cat("\nsaving...\n")
    capture.output(
                   cat("\ncols => percentage to calculate a value for p\n"),
                   cat("rows => base graph for H\n\n"),
                   file = file.path(results_path,
                                    paste0("result_try", try, ".txt")
                                   )
                  )
    for(i in 1:stats_number){
      capture.output(cat(result_mt_names[i], "MATRIX:\n"),
                     file = file.path(results_path,
                                      paste0("result_try", try, ".txt")
                                     ),
                     append = TRUE
                    )
      capture.output(print(result_mt[[i]]),
                     cat("-------------------------------------\n"),
                     file = file.path(results_path,
                                      paste0("result_try", try, ".txt")
                                     ),
                     append = TRUE
                    )
    }
    cat("saved!\n")
  }
  cat("\n\n################################################################\n")
}

#### Printing MEANS BETWEANS TRIES ####
res_tmp <- array(NA, c(nrow=length(Raw_names) + 1, ncol=length(p_stages), 2))
dimnames(res_tmp) <- list( c(c("control"), Raw_names),
                           p_stages,
                           c("mean" , "standard deviation")
                         )

cat("########## MEANS BETWEEN TRIES ##########\n")
for(atr in 1:(2*stats_number)){

  cat("\n")
  cat(mean_mt_names[atr], "MATRIX:\n")
  cat("Contains values of (",
      result_mt_desc[((atr - 1)%%stats_number) + 1],
      ")\nfor different H's and values of p\n", sep = "")
  cat("\ncols => percentage to calculate a value for p\n")
  cat("rows => base graph for H\n\n")
  
  for(i in 1:(length(Raw_names) + 1)){
    for(j in 1:length(p_stages)){
      res_tmp[i,j,1] <- mean(result_mt_mean[i,j,,atr])
      res_tmp[i,j,2] <- sd(result_mt_mean[i,j,,atr])
    }
  }
  print(res_tmp)
  cat("-------------------------------------\n")
}
cat("#########################################\n")

cat("Hamm / (p * n2a2)\n")
tmp = result_mt[[6]]
#print(tmp)
for (p_index in 1:length(p_stages)){
  tmp[,p_index] = tmp[,p_index]/(p_stages[p_index]*n*(n-1)/2)
}
print(tmp)

cat("\nScript done!\n")
###Script ending##
