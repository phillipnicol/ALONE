##### main.R
##### Phillip Nicol
##### November 19, 2019

require(igraph)
require(Matrix)
require(Rcpp)
require(inline)


##### Function GA performs a genetic algorithm to approximate MLE
##### Data - Binary Data frame, as generated from function Generate_Data
##### N - number of individuals per generation
##### MAX_GEN - number of generations
##### p_mat_c - probability of matrix crossover between two individuals
##### p_cx - probability of cycle crossover between two individuals
##### p_perm_mut - probability of permutation mutation between two individuals
##### e - false positive rate
##### d - false negative rate
##### suppress -  print to console?
##### method - "SA" (Spontaneous activation) or "ME" (Measurement error)
##### penalty - "hard" (3 nodes per parent), "BIC", or "none"
                

GA <- function(Data, N=100, MAX_GEN=100, p_mat_c=0.5, p_cx=0.5, p_ba=0.01, p_ea=0.05, p_perm_mut=0.01, suppress = FALSE,
               method = "SA", penalty = "hard")
{
  
  if(colMeans(Data)[1] != 1) {
    warning("The first column of data should be all 1's")
    Data <- cbind(1,Data)
  }
  
  leaky <- rep(min(colMeans(Data))/2, ncol(Data))
  
  n <- ncol(Data)
  num_samples <- nrow(Data)
  
  if(penalty != "hard" && penalty != "BIC")
  {
    penalty = "none"
  }
  
  if(method == "ME")
  {
    e <- leaky[1]
    d <- leaky[2]
  }
  
  if(!suppress)
  {
    cat("Preparing Data... ...\n")
  }
  Data <- as.matrix(Data)
  Unique_Data <- matrix(0, nrow = 1, ncol = n)
  Unique_Data[1,] = Data[1,]
  counts <- 1
  ###Given Data, Find Unique Data Values. Extract counts
  for(i in 2:nrow(Data))
  {
    for(j in 1:nrow(Unique_Data))
    {
      if(identical(Data[i,], Unique_Data[j,]) == TRUE)
      {
        counts[j] = counts[j] + 1
        break
      }
      if(j == nrow(Unique_Data))
      {
        Unique_Data = rbind(Unique_Data, Data[i,])
        counts[j+1] = 1
      }
    }
  }
  
  if(method == "ME")
  {
    ### Create all binary vectors of size n such that 1 is in the first slot
    
    dec_vecs = c(0 : (2^{n-1} - 1))
    Bin_vec = matrix(0, nrow = 2^{n-1}, ncol = n-1)
    for(i in 1:nrow(Bin_vec))
    {
      Bin_vec[i,] = number2binary(dec_vecs[i], n-1)
    }
    Bin_vec = cbind(1, Bin_vec)
    
    ### Create a table of P(Y=y|X=x)
    ### Rows of CPT are unique Data points
    CPT <- matrix(0, nrow = nrow(Unique_Data), ncol = nrow(Bin_vec))
    
    for(i in 1:nrow(Unique_Data))
    {
      for(j in 1:nrow(Bin_vec))
      {
        CPT[i,j] = Y_Given_X(Bin_vec[j,-1], Unique_Data[i,-1], e, d)
      }
    }
  }
  
  mode(Unique_Data) <- "integer"
  
  Node_counts <- list()
  for(i in 2:n)
  {
    Node_counts[[i]] = which(Data[,i] == 1)
  }
  
  ##### INITIALIZE POPULATION ######
  
  if(!suppress)
  {
    cat("Initializing Population ... ... \n")   
  }
  
  ### Initialize Binary matrix chromosomes
  Mats <- list()
  for(i in 1:N)
  {
    Mats[[i]] <- Generate_Tree_DAG(n)
  }
  
  ### The other chromosome will be random permutations
  Perms <- list()
  for(i in 1:N)
  {
    Perms[[i]] = c(1, sample(2:n, size = n-1, replace = FALSE))
    Perms[[i]] = Perm_Repair(Mats[[i]], Perms[[i]])
  }
  
  ### For each Graph, infer theta
  ### Calculate fitness for each member of the population
  fitness <- c(1:N)
  for(i in 1:N)
  {
    p = invPerm(Perms[[i]])
    Temp = Mats[[i]][p,]
    theta = Infer_Theta(Temp[,p], Node_counts, counts)
    if(method == "ME")
    {
      List = EM(Temp[,p], theta, CPT, Unique_Data, counts, 0.01)
      fitness[i] = List[[1]]
    }
    else
    {
      fitness[i] = Fitness(Mats[[i]], Perms[[i]], theta, leaky, Unique_Data, counts)
    }
    
    if(penalty == "BIC")
    {
      fitness[i] = fitness[i] -  log(num_samples)*log(n)*sum(Mats[[i]])
    }
  }
  
  ###Save the maximum value of the graph
  key = which(fitness == max(fitness))
  Max_fitness = max(fitness)
  Max_perm = Perms[[key[1]]]
  Max_edge = Mats[[key[1]]]
  max_gen = 1
  
  ###Initialize Data Frame to store information
  Results <- matrix(0, nrow = MAX_GEN, ncol = 2)
  Results <- as.data.frame(Results)
  
  if(!suppress)
  {
    cat("Done! Starting GA ... ... \n")    
  }
  
  ##### Genetic Algorithm #####
  generation = 1
  while(generation <= MAX_GEN)
  {
    if(generation %% 10 == 0)
    {
      if(!suppress)
      {
        cat("Generation ", generation, " ... ... \n")
      }
    }
    
    ###Create copies of the lists we currently have
    Mats_c = Mats
    Perms_c = Perms
    fitness_c = fitness
    
    dart_board = Create_Dart_Board(fitness, N)
    
    ### Selection
    for(i in 1:(N/2))
    {
      selected = Selection(dart_board)
      
      Children = DAGrR(Mats[[selected[1]]], Mats[[selected[2]]], Perms[[selected[1]]],
                       Perms[[selected[2]]], fitness[selected[1]], fitness[selected[2]],
                       p_mat_c, p_cx, p_ba, p_ea, p_perm_mut, penalty)
      
      Mats_c[[2*i - 1]] = Children[[1]]
      Mats_c[[2*i]] = Children[[3]]
      Perms_c[[2*i - 1]] = Children[[2]]
      Perms_c[[2*i]] = Children[[4]]
      
      if(identical(Children[[1]], Mats[[selected[1]]]) & 
         identical(Children[[2]], Perms[[selected[1]]]))
      {
        fitness_c[2*i - 1] = fitness[selected[1]]
      }
      else
      {
        p = invPerm(Children[[2]])
        Temp = Children[[1]][p,]
        theta = Infer_Theta(Temp[,p], Node_counts, counts)
        if(method == "ME")
        {
          List = EM(Temp[,p], theta, CPT, Unique_Data, counts, 0.01)
          fitness_c[2*i - 1] = List[[1]]
        }
        else
        {
          fitness_c[2*i - 1] = Fitness(Children[[1]], Children[[2]], theta, leaky, Unique_Data,
                                     counts)
        }
        
        if(penalty == "BIC")
        {
          fitness_c[2*i - 1] = fitness_c[2*i - 1] - log(n)*log(num_samples)*sum(Temp)
        }
      }
      
      if(identical(Children[[3]], Mats[[selected[2]]]) & 
         identical(Children[[4]], Perms[[selected[2]]]))
      {
        fitness_c[2*i] = fitness[selected[2]]
      }
      else
      {
        p = invPerm(Children[[4]])
        Temp = Children[[3]][p,]
        theta = Infer_Theta(Temp[,p], Node_counts, counts)
        if(method == "ME")
        {
          List = EM(Temp[,p], theta, CPT, Unique_Data, counts, 0.01)
          fitness_c[2*i] = List[[1]]
        }
        else
        {
          fitness_c[2*i] = Fitness(Children[[3]], Children[[4]], theta, leaky, Unique_Data, counts)
        }
        
        if(penalty  == "BIC")
        {
          fitness_c[2*i] = fitness_c[2*i] - log(n)*log(num_samples)*sum(Temp)
        }
      }
    }
    
    Mats = Mats_c
    Perms = Perms_c
    fitness = fitness_c
    
    Results[generation, 1] = max(fitness)
    Results[generation, 2] = mean(fitness)
    
    if(max(fitness) > Max_fitness)
    {
      key = which(fitness == max(fitness))
      Max_fitness = max(fitness)
      Max_perm = Perms[[key[1]]]
      Max_edge = Mats[[key[1]]]
      max_gen = generation
    }
    generation = generation + 1
  }

  if(!suppress)
  {
    cat("Done! \n")
  }
  p = invPerm(Max_perm)
  Temp = Max_edge[p,]
  theta <- Infer_Theta(Temp[,p], Node_counts, counts)
  
  G <- graph_from_adjacency_matrix(t(Temp[,p]), mode = "directed")
  
  Return_List <- list()
  
  cat("Call summary: \n Method: ", method, "\n Penalty:", penalty, "\n")
  
  Return_List$most_fit = Results$V1
  Return_List$mean_fit = Results$V2
  Return_List$adj_mat = Temp[,p]
  Return_List$igraph_obj = G
  Return_List$score = Max_fitness
  Return_List$theta = theta
  return(Return_List)
}


