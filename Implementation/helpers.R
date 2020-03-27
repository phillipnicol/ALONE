###Helpers.R
###November 19, 2019
###Phillip Nicol


##### Genetic Operators

### Function Matrix Crossover
### Takes two topologically sorted DAGs and produces 2 baby DAGs

Matrix_Crossover <- function(M1, M2)
{
  n <- nrow(M1)
  ###Define children
  C1 = M1
  C2 = M2
  
  for(i in 3:n)
  {
    if(i %% 2 == 0)
    {
      C1[i,] = 0
      x <- which(M1[i,] != 0)
      
      C1[i,x[1]] = 1
      
      x <- x[-1]
      for(j in x)
      {
        if(Redundant(j, i, C1) == FALSE)
        {
          C1[i, j] = 1
        }
      }
      
      C2[i,] = 0
      x <- which(M2[i,] != 0)
      
      C2[i,x[1]] = 1
      
      x <- x[-1]
      for(j in x)
      {
        if(Redundant(j, i, C2) == FALSE)
        {
          C2[i, j] = 1
        }
      }
    }
    else
    {
      C1[i,] = 0
      x <- which(M2[i,] != 0)
      
      C1[i,x[1]] = 1
      
      x <- x[-1]
      for(j in x)
      {
        if(Redundant(j, i, C1) == FALSE)
        {
          C1[i, j] = 1
        }
      }
      
      C2[i,] = 0
      x <- which(M1[i,] != 0)
      
      C2[i,x[1]] = 1
      
      x <- x[-1]
      for(j in x)
      {
        if(Redundant(j, i, C2) == FALSE)
        {
          C2[i, j] = 1
        }
      }
    }
  }
  
  Children = list()
  Children[[1]] = C1
  Children[[2]] = C2
  return(Children)
}

Edge_Mutation <- function(M, hard_penalty)
{
  Candidate_Changes <- list()
  count <- 1
  n <- nrow(M)
  
  ### Consider all possible edge additions
  for(i in 3:n)
  {
    if(sum(M[i,]) < 3 & hard_penalty)
    {
      for(j in 2:(i-1))
      {
        if(M[i, j] == 0 & Redundant(j, i, M) == FALSE & Overpower(j, i, M) == FALSE)
        {
          M_new = M
          M_new[i, j] = 1
          Candidate_Changes[[count]] = M_new
          count = count + 1
        }
      }
    }
    else
    {
      for(j in 2:(i-1))
      {
        if(M[i, j] == 0 & Redundant(j, i, M) == FALSE & Overpower(j, i, M) == FALSE)
        {
          M_new = M
          M_new[i, j] = 1
          Candidate_Changes[[count]] = M_new
          count = count + 1
        }
      }     
    }
  }
  
  if(rbinom(1, 1, 0.5) == 1 & length(Candidate_Changes) != 0)
  {
    key = sample(1:length(Candidate_Changes), 1, replace = FALSE)
    return(Candidate_Changes[[key]])
  }
  
  
  ### Consider all possible edge removals
  for(i in 2:n)
  {
    parents = which(M[i,] != 0)
    
    if(length(parents) > 1)
    {
      for(j in parents)
      {
        M_new = M
        M_new[i, j] = 0
        Candidate_Changes[[count]] <- M_new
        count = count + 1
      }
    }
  }
  
  if(length(Candidate_Changes) == 0)
  {
    return(M)
  }
  
  key = sample(1:length(Candidate_Changes), 1, replace = FALSE)
  return(Candidate_Changes[[key]])
}

Branch_Adjust <- function(M, node)
{
  ### Remove ties with other branches
  children <- c(node, which(M[,node] != 0))
  
  repeat
  {
    children_orig = children
    
    for(i in children)
    {
      children <- c(children, which(M[,i] != 0))
    }
    
    children = unique(children)
    
    if(length(children_orig) == length(children))
    {
      break
    }
  }
  
  set <- c(1:(node - 1))
  parents = which(M[node,] != 0)
  set = set[-parents]
  if(length(set) == 0)
  {
    return(M)
  }
  
  for(i in children)
  {
    for(j in 1:nrow(M))
    {
      if(M[i,j] == 1 & (j %in% children == FALSE))
      {
        M[i,j] = 0
      }
    }
  }
  
  ### Move the branch somewhere
  new_node = sample(set, 1, replace = FALSE)
  M[node, new_node] = 1
  
  return(M)
}


Redundant <- function(parent_node, child_node, M)
{
  G = graph_from_adjacency_matrix(t(M), mode = "directed")
  Paths = all_simple_paths(G, 1, parent_node)
  
  if(length(Paths) == 0)
  {
    return(TRUE)
  }
  for(i in 1:length(Paths))
  {
    Paths[[i]] = as_ids(Paths[[i]])
  }
  
  parents <- which(M[child_node,] != 0)
  
  count <- 0
  for(i in 1:length(Paths))
  {
    if(sum(parents %in% Paths[[i]]) > 0)
    {
      count = count + 1
    }
  }
  
  if(count == length(Paths))
  {
    return(TRUE)
  }
  else
  {
    return(FALSE)
  }
}

Overpower <- function(parent_node, child_node, M)
{
  vec = child_node
  
  repeat
  {
    vec_copy = vec
    for(i in vec)
    {
      vec = c(vec, which(M[i,] != 0))
    }
    vec = vec[-c(1:length(vec_copy))]
    vec = unique(vec)
    
    if(length(vec) == 1 & vec[1] == parent_node)
    {
      return(TRUE)
    }
    else if(length(vec) == 1 & vec[1] == 1)
    {
      break
    }
  }
  
  return(FALSE)
}



CX <- function(p1, p2, n)
{
  Children = matrix(1, nrow = 2, ncol = n)
  p1 = p1[-1]
  p2 = p2[-1]
  p1 = p1 - 1
  p2 = p2 - 1
  
  ###Identify cycles in first parent
  cycle = p1[1]
  cycle_over = FALSE
  
  while(cycle_over == FALSE)
  {
    key = which(p1 == cycle[length(cycle)])
    if(p2[key] != p1[1])
    {
      cycle = c(cycle, p2[key])
    }
    else
    {
      cycle_over = TRUE
    }
  }
  for(i in 1:(n-1))
  {
    Children[1, i+1] = ifelse(p1[i] %in% cycle, p1[i], p2[i])
  }
  
  Children[1, 2:n] = Children[1,2:n] + 1
  
  ###Rinse and repeat
  
  cycle = p2[1]
  cycle_over = FALSE
  
  while(cycle_over == FALSE)
  {
    key = which(p2 == cycle[length(cycle)])
    if(p1[key] != p2[1])
    {
      cycle = c(cycle, p1[key])
    }
    else
    {
      cycle_over = TRUE
    }
  }
  
  for(i in 1:(n-1))
  {
    Children[2, i+1] = ifelse(p2[i] %in% cycle, p2[i], p1[i])
  }
  
  Children[2, 2:n] = Children[2, 2:n] + 1
  
  return(Children)
}



Permutation_Mutation <- function(perm_chromosome, n)
{
  keys = sample(2:n, 2, replace = FALSE)
  
  perm_chromosome <- replace(perm_chromosome, keys, perm_chromosome[c(keys[2], keys[1])])
  
  return(perm_chromosome)
}


Perm_Repair <- function(G, p)
{
  n <- nrow(G)
  set <- c(2:n)
  
  while(length(set) != 0)
  {
    row = set[1]
    eq_rows = which(apply(G, 1, function(x) return(all(x == G[row,]))))
    
    if(length(eq_rows) != 1)
    {
      eq_cols = which(apply(G, 2, function(x) return(all(x == G[,row]))))
      eq_cols = intersect(eq_cols, eq_rows)
      eq_rows = eq_cols
      
      if(length(eq_rows) > 0)
      {
        ###Sort
        extract = p[eq_rows]
        extract = extract[order(extract, decreasing = FALSE)]
        p[eq_rows] = extract
      }
    }
    
    set = setdiff(set, eq_rows)
  }
  
  return(p)
}



DAGrR <- function(m1, m2, p1, p2, f1, f2, p_mc, p_cx, p_ba, p_ea, p_perm_mut, penalty)
{
  hard <- FALSE
  if(penalty == "hard")
  {
    hard = TRUE
  }
  
  ### Define return
  Children <- list()
  n <- nrow(m1)
  
  recomb = TRUE
  
  ### Check for equivalence
  ### THIS IS NOT CURRENTLY ACCURATE
  ### TO DO: IMPROVE THIS CHECK
  if(f1 == f2)
  {
    recomb = FALSE
  }
  
  child_m1 = m1
  child_m2 = m2
  child_p1 = p1
  child_p2 = p2
  
  if(rbinom(1, 1, p_mc) == 1 & recomb == TRUE)
  {
    Children = Matrix_Crossover(child_m1, child_m2)
    child_m1 = Children[[1]]
    child_m2 = Children[[2]]
  }
  
  if(rbinom(1, 1, p_cx) == 1 & recomb == TRUE)
  {
    child_p = CX(p1, p2, nrow(m1))
    child_p1 = child_p[1,]
    child_p2 = child_p[2,]
  }
  
  if(rbinom(1, 1, p_ba) == 1)
  {
    child_m1 = Branch_Adjust(child_m1, sample(3:n, 1, replace = FALSE))
  }
  
  if(rbinom(1, 1, p_ba) == 1)
  {
    child_m2 = Branch_Adjust(child_m1, sample(3:n, 1, replace = FALSE))
  }
  
  if(rbinom(1, 1, p_ea) == 1)
  {
    child_m1 = Edge_Mutation(child_m1, hard)
  }
  
  if(rbinom(1, 1, p_ea) == 1)
  {
    child_m2 = Edge_Mutation(child_m2, hard)
  }
  
  if(rbinom(1, 1, p_perm_mut) == 1)
  {
    child_p1 = Permutation_Mutation(child_p1, nrow(m1))
  }
  
  if(rbinom(1, 1, p_perm_mut) == 1)
  {
    child_p2 = Permutation_Mutation(child_p2, nrow(m1))
  }
  
  child_p1 = Perm_Repair(child_m1, child_p1)
  child_p2 = Perm_Repair(child_m2, child_p2)
  
  Children[[1]] = child_m1
  Children[[2]] = child_p1
  Children[[3]] = child_m2
  Children[[4]] = child_p2
  return(Children)
}




Fitness <- function(m, p, theta, leaky, Unique_Data, counts)
{
  ###Permute the matrix
  p = invPerm(p)
  Temp = m[p,]
  m = Temp[,p]
  
  val = 0
  
  for(i in 1:nrow(Unique_Data))
  {
    val = val + counts[i]*log(Likelihood(Unique_Data[i,], m, theta, leaky))
  }
  
  return(val)
}

Likelihood <- function(x, G, theta, leaky)
{
  N <- nrow(G)
  val <- 1
  for(i in 2:N)
  {
    G_t = G[i,]
    if(min(1, sum(G_t * x)) == 1)
    {
      if(x[i] == 1)
      {
        val = val*theta[i]
      }
      else
      {
        val = val*(1 - theta[i])
      }
    }
    else
    {
      if(x[i] == 1)
      {
        val = val*leaky[i]
      }
      else
      {
        val = val*(1 - leaky[i])
      }
    }
  }
  return(val)
}

Create_Dart_Board <- function(fitness, N)
{
  fitness = fitness - min(fitness)
  dart_board = c(1:N)
  sum_vec = sum(fitness)
  
  if(sum_vec == 0)
  {
    dart_board = seq(0, 1, 1/N)
    dart_board = dart_board[-1]
    return(dart_board)
  }
  
  for(i in 1:length(fitness))
  {
    dart_board[i] = (sum(fitness[1:i])) / sum_vec
  }
  
  return(dart_board)
}

Selection <- function(dart_board)
{
  indices <- c(1,2)
  
  dart_1 <- runif(1, 0, 1)
  
  for(i in 2:length(dart_board))
  {
    if(dart_1 < dart_board[i])
    {
      indices[1] = i
      break
    }
  }
  
  dart_2 <- runif(1, 0, 1)
  
  for(i in 2:length(dart_board))
  {
    if(dart_2 < dart_board[i])
    {
      indices[2] = i
      break
    }
  }
  
  return(indices)
}


Generate_Tree_DAG <- function(n)
{
  M <- matrix(0, nrow = n, ncol = n)
  M[2, 1] = 1
  
  for(i in 3:n)
  {
    key = sample(1:(i - 1), 1, replace = FALSE)
    M[i, key] = 1
  }
  
  return(M)
}


Infer_Theta <- function(G, Node_counts, counts)
{
  n <- nrow(G)
  theta <- c(1:n)
  
  for(i in 2:n)
  {
    parents = Parents(G, i)
    if(1 %in% parents)
    {
      denom = sum(counts)
      num = length(Node_counts[[i]])
    }
    else
    {
      set = Node_counts[[parents[1]]]
      for(j in parents)
      {
        set = union(set, Node_counts[[j]])
      }
      denom = length(set)
      set = intersect(set, Node_counts[[i]])
      num = length(set)
    }
    
    theta[i] = num/denom
    
    if(num == 0)
    {
      theta[i] = 1/sum(counts)
    }
    if(denom == 0)
    {
      theta[i] = 1/sum(counts)
    }
    if(theta[i] == 1)
    {
      theta[i] = (sum(counts) - 1)/sum(counts)
    }
  }
  return(theta)
}



Infer_Spontaneous <- function(G, Node_counts, counts)
{
  n <- nrow(G)
  spont <- c(1:n)
  spont[1] = 0
  
  for(i in 2:n)
  {
    parents = Parents(G, i)
    if(1 %in% parents)
    {
      spont[i] = 0
    }
    else
    {
      set = setdiff(c(1:sum(counts)), Node_counts[[parents[1]]])
      for(j in parents)
      {
        set = intersect(set, setdiff(c(1:sum(counts)), Node_counts[[j]]))
      }
      set = intersect(set, Node_counts[[i]])
      num = length(set)
      denom = length(Node_counts[[i]])
      
      if(denom == 0 | num == 0)
      {
        spont[i] = 0.01
      }
      else
      {
        spont[i] = num/denom
      }
      
      if(spont[i] == 1)
      {
        spont[i] = 0.99
      }
    }
  }
  return(spont)
}

Parents <- function(graph, node)
{
  return(which(graph[node,] != 0))
}





#### CODE SPECIFICALLY FOR ME MODEL
EM <- function(G, theta_init, CPT, Unique_Data, counts, CONVERGENCE)
{
  
  n <- length(theta_init)
  difference <- replicate(n, 1)
  counter <- 1
  vec = convert_mat(as.integer(n), as.integer(G))
  theta <- theta_init
  
  while(max(difference) > CONVERGENCE)
  {
    likelihood_vec = likelihood(theta, vec)
    theta_prior = theta
    P_y = as.vector(t(CPT %*% likelihood_vec))
    
    for(j in 2:n)
    {
      parents = Parents(G, j)
      if(1 %in% parents)
      {
        denom = sum(counts)
        set_num = indices_denom(as.integer(n), as.integer(j))
        a <- as.vector(t(CPT[,set_num] %*% likelihood_vec[set_num]))
        num = sum(counts*(P_y)^{-1}*a)
      }
      else
      {
        set_denom = indices_denom(as.integer(n), as.integer(parents))
        set_num = indices_num(as.integer(n), as.integer(parents), as.integer(j))
        set_denom = unique(set_denom)
        set_num = unique(set_num)
        a <- as.vector(t(CPT[,set_num] %*% likelihood_vec[set_num]))
        b <- as.vector(t(CPT[,set_denom] %*% likelihood_vec[set_denom]))
        num = sum(counts*(P_y)^{-1}*a)
        denom = sum(counts*(P_y)^{-1}*b)
      }
      
      theta[j] = num/denom
    }
    
    difference = abs(theta - theta_prior)
    counter = counter+1
  }
  
  likelihood_vec = likelihood(theta, vec)
  P_y = as.vector(t(CPT %*% likelihood_vec))
  Return_list <- list()
  Return_list[[1]] = sum(counts*log(P_y))
  Return_list[[2]] = theta
  return(Return_list)
}

cppFunction('IntegerVector indices_denom(int n, IntegerVector parents) 
            {
            IntegerVector out(0);
            
            for(int i = 0; i < (1 << (n-1)); i++)
            {
            for(int j = 0; j < parents.size(); j++)
            {
            if((i & (1 << (n-parents[j]))) != 0)
            {
            out.push_back(i+1);
            }
            }
            }
            
            return out;
            }')


cppFunction('IntegerVector indices_num(int n, IntegerVector parents, int node) 
            {
            IntegerVector out(0);
            
            for(int i = 0; i < (1 << (n-1)); i++)
            {
            for(int j = 0; j < parents.size(); j++)
            {
            if(((i & (1 << (n-parents[j]))) != 0) & ((i & (1 << (n - node))) != 0))
            {
            out.push_back(i+1);
            }
            }
            }
            
            return out;
            }')


code_likelihood <- "
double *ptheta = REAL(Theta);
double *pout;

int *G = INTEGER(G_);

int n = length(Theta);
int m = 1 << (n-1);
int data = m;
int test = 2;
double theta;
int row;
SEXP out = PROTECT(allocVector(REALSXP, m));
pout = REAL(out);

for(int i = 0; i < m; i++)
{
  pout[i] = 1;
  for(int j = 1; j < n; j++)
  {
  test = 1 << (n-(j + 1));
  row = G[j];
  if((G[j] & data) == 0)
  {
  if((data & test) != 0)
  {
  pout[i] = 0;
  break;
  }
  }
  else
  {
  theta = ptheta[j];
  if((data & test) != 0)
  {
  pout[i] *= theta;
  }
  else
  {
  pout[i] *= (1 - theta);
  }
  }
  }
  data++;
}

UNPROTECT(1);
return out;
"

likelihood <- cfunction(c(Theta = "numeric", G_ = "integer"), body = code_likelihood, language = c("C"))



code_convert_mat <- "
int n = asInteger(N);
int *G = INTEGER(G_);
int *pout;
int add;

SEXP out = PROTECT(allocVector(INTSXP, n));
pout = INTEGER(out);

for(int i = 0; i < n; ++i)
{
  pout[i] = 0;
}

for(int i = 0; i < n; ++i)
{
  for(int j = 0; j < n; ++j)
  {
  if(G[n*i + j] == 1)
  {
  add = 1 << (n - (i+1));
  pout[j] += add;
  }
  }
}

UNPROTECT(1);
return out;
"

convert_mat <- cfunction(c(N = "integer", G_ = "integer"), body = code_convert_mat, language = c("C"))


### Code copied from
###https://stackoverflow.com/questions/12088080/how-to-convert-integer-number-into-binary-vector

number2binary = function(number, noBits) {
  binary_vector = rev(as.numeric(intToBits(number)))
  if(missing(noBits)) {
    return(binary_vector)
  } else {
    binary_vector[-(1:(length(binary_vector) - noBits))]
  }
}


Dist <- function(A, B)
{
  if(length(A) == 0)
  {
    return(0)
  }
  
  return(length(which(A != B)))
}

Y_Given_X <- function(hidden_data, true_data, e, d)
{
  hidden_data = hidden_data[-1]
  true_data = true_data[-1]
  
  hidden_data_1 <- hidden_data[which(hidden_data == 1)]
  hidden_data_0 <- hidden_data[which(hidden_data == 0)]
  
  distance_1 <- Dist(hidden_data_1, true_data[which(hidden_data == 1)])
  distance_0 <- Dist(hidden_data_0, true_data[which(hidden_data == 0)])
  
  prob = e^{distance_0}*(1-d)^{length(hidden_data_0) - distance_0}
  prob = prob*d^{distance_1}*(1-e)^{length(hidden_data_1) - distance_1}
  
  return(prob)        
}
