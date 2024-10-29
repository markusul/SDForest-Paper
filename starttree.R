predict_outsample <- function(tree, X){
  # predict for every observation in X f(x)
  # using the splitting rules from the tree
  if(is.null(dim(X))){
    return(traverse_tree(tree, X))
  }
  apply(X, 1, function(x)traverse_tree(tree, x))
}

#helper functions to label nodes for plotting

split_names <- function(node, var_names = NULL){
  if(is.null(var_names)){
    node$label <- paste('X', node$j, ' <= ', round(node$s, 2), sep = '')
  }else{
    node$label <- paste(var_names[node$j], ' <= ', round(node$s, 2), sep = '')
  }
}

leave_names <- function(node){
  new_name <- as.character(round(node$value, 1))
  if(new_name %in% node$Get('name', filterFun = data.tree::isLeaf)){
    new_name <- paste(new_name, '')
  }
  node$label <- new_name
}

# finds all the reasonable spliting points in a data matrix
find_s <- function(X, max_candidates = 100){
  p <- ncol(X)
  if(p == 1){
    X <- matrix(X, ncol = 1)
  }
  n <- nrow(X)
  
  X_sort <- apply(X, 2, sort, method = 'quick')
  
  if(is.null(dim(X_sort))){
    X_sort <- matrix(X_sort, ncol = p)
  }
  
  # find middle points between observed x values
  s <- X_sort[-nrow(X_sort), ] + diff(X_sort)/2
  
  # for runtime reasons limit split candidates
  if(nrow(s) > max_candidates){
    s <- s[unique(floor(seq(1, dim(s)[1], length.out = max_candidates))), ]
  }
  
  if(is.null(dim(s))){
    matrix(s, ncol = p)
  }else{
    s
  }
}

traverse_tree <- function(tree, x){
  # traverse the tree using the splitting rules and 
  # returns point estimate for f(x)
  if(tree$isLeaf){
    return(tree$value)
  }
  if(x[tree$j] <= tree$s){
    traverse_tree(tree$children[[1]], x)
  }else {
    traverse_tree(tree$children[[2]], x)
  }
}

loss <- function(Y, f_X){
  as.numeric(sum((Y - f_X)^2) / length(Y))
}

pruned_loss <- function(tree, X_val, Y_val, Q_val, t){
  # function to prune tree using the minimum loss decrease t
  # and return spectral loss on the validation set
  
  tree_t <- data.tree::Clone(tree)
  
  # prune tree
  data.tree::Prune(tree_t, function(x) x$dloss > t)
  
  # predict on test set
  f_X_hat_val <- predict_outsample(tree_t, X_val)
  
  # return spectral loss
  sum((Q_val %*% Y_val - Q_val %*% f_X_hat_val) ** 2) / length(Y_val)
}

randomTree <- function(X, m = nrow(X) / 10, js = 1:ncol(X), min_sample = 5, 
                       outInter = c(-50, 50), make_tree = T){
  
  n <- nrow(X)
  p <- ncol(X)
  cp_max <- 1 
  # generate tree
  if(make_tree){
    tree <- data.tree::Node$new(name = '1', 
                                value = runif(1, outInter[1], outInter[2]), 
                                cp_max = 10)
  }
  
  var_imp <- rep(0, p)
  var_names <- colnames(data.frame(X))
  names(var_imp) <- var_names
  
  # partitions of observations
  index <- list(1:n)
  i <- 1
  while (i <= m){
    # get number of observations in each partition
    samples_per_part <- unlist(lapply(index, function(x)length(x)))
    
    # get potential splits (partitions with enough observations)
    potential_splitts <- which(samples_per_part >= 2 * min_sample)
    if(length(potential_splitts) == 0){
      break
    }
    
    # sample potential partition to split
    # probability of each partition is proportional to number of observations
    if(length(potential_splitts) == 1){
      branch <- potential_splitts
    }else {
      prob_part <- samples_per_part[potential_splitts]/
        sum(samples_per_part[potential_splitts])
      branch <- sample(potential_splitts, 1, prob = prob_part)
    }
    
    # sample covariate to split on
    j <- sample(js, 1)
    var_imp[j] <- var_imp[j] + 1
    
    branch_index <- index[[branch]]
    
    # sample split point
    X_b_j <- X[branch_index, j]
    potential_s <- sort(X_b_j)
    low_s <- potential_s[min_sample]
    top_s <- potential_s[length(potential_s) - min_sample + 1]
    s <- rbeta(1, 2, 2) * (top_s - low_s) + low_s
    
    # split partition
    index <- append(index, list(branch_index[X_b_j > s]))
    index[[branch]] <- branch_index[X_b_j <= s]
    
    # add split to tree
    if(make_tree){
      if(tree$height == 1){
        leave <- tree
      }else{
        leaves <- tree$leaves
        leave <- leaves[[which(tree$Get('name', filterFun = data.tree::isLeaf) == branch)]]
      }
      leave$j <- j
      leave$s <- s
      leave$res_dloss <- 1
      cp_max <- cp_max - 0.001
      leave$AddChild(branch, value = runif(1, outInter[1], outInter[2]), cp_max = cp_max)
      leave$AddChild(i + 1, value = runif(1, outInter[1], outInter[2]), cp_max = cp_max)
    }
    i <- i + 1
  }
  
  # sample means per partition
  if(make_tree){
    f_X <- predict_outsample(tree, X)
  }else {
    # sample means per partition
    f_X_means <- runif(length(index), outInter[1], outInter[2])
    
    # generate f_X
    f_X <- rep(0, n)
    for (i in 1:length(index)){
      f_X[index[[i]]] <- f_X_means[i]
    }
  }
  
  res <- list(predictions = f_X, var_names = var_names, 
              var_importance = var_imp)
  
  # return data
  if(make_tree){
    # labels for the nodes
    tree$Do(split_names, filterFun = data.tree::isNotLeaf, var_names = var_names)
    tree$Do(leave_names, filterFun = data.tree::isLeaf)    

    res$tree <- tree
    class(res) <- 'SDTree'
    return(res)
  }
  
  return(res)
}

randomForest <- function(X, nTree = 100, js = 1:ncol(X), min_sample = 5,
                         outInter = c(-50, 50), make_model = T, 
                         m = nrow(X) / 10){
  p <- ncol(X)
  n <- nrow(X) 
  forest <- pbapply::pblapply(1:nTree, function(tree) randomTree(X, js = js, 
                                                  min_sample = min_sample, 
                                                  outInter = outInter, 
                                                  make_tree = make_model))
  f_X <- rowMeans(sapply(forest, function(tree) tree$predictions))
  
  if(!make_model)
    return(f_X)
  
  # variable importance
  var_imp <- sapply(forest, function(x){x$var_importance})
  if(p > 1){
    var_imp <- rowMeans(var_imp)
  }else {
    var_imp <- mean(var_imp)
  }
  
  oob_ind <- lapply(1:n, function(i) 1:nTree)
  
  output <- list(predictions = f_X, 
                 forest = forest, 
                 var_names = colnames(data.frame(X)), 
                 oob_loss = 0, 
                 oob_SDloss = 0, 
                 var_importance = var_imp, 
                 oob_ind = oob_ind, 
                 oob_predictions = f_X)
  class(output) <- 'SDForest'
  output
}


X <- matrix(rnorm(1000 * 2), nrow = 1000)

tree <- randomTree(X, make_tree = T)

tree
library(SDForest)

plot(tree)
predict(tree, data.frame(X)) == tree$predictions
prune(tree, 0.9999)
print(tree$tree, 'cp_max')

tree
tree$var_importance
tree2 <- SDTree(y = tree$predictions, x = X, Q_type = 'no_deconfounding', cp = 0.01)
plot(tree2)
tree$var_names
predict_outsample(tree$tree, data.frame(X))
predict(tree, data.frame(X))

plot(tree$predictions, tree2$predictions)

dep <- partDependence(tree, 1, X)
plot(dep)


forest <- randomForest(X, nTree = 2, make_model = T, js = c(1, 2), m = 5)
forest$predictions == predict(forest, data.frame(X))

dep <- partDependence(forest, 1, X, subSample = 100)
plot(dep)
