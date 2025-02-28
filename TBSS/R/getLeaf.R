
getLeafGeneric = function(){
  result = vector("character")
  depth  = 0
  function(node,tree,leafs)
  {  
    if(node %in% leafs)
    {
      res <- node
      result[length(result) +1] <<- res
    }
    children = unique(tree[tree[,2] == node,1])
    if(length(children)>0)
    {
      for(ch in children)
      {
         res <- Recall(ch,tree,leafs)
      }
    }
    else
    {
      if(!(node %in% leafs))
      {
        res <- node
        result[length(result) +1] <<- res
      }
            depth <<- depth +1
	    return(result)
    }
    return(res)
  }
}


#'  Retrieve leaves linked to a node
#'@description Internal function that list the leveas connected to a node  
#'@param `node` a character vector with the id of a node
#'@param `tree`: a data.frame with two columns. The first column is the id of each node of a tree while the second is its direct ancestor. Root node as no ancestor that is represented with "".
#'@param `leafs`:  a list of the terminal nodes, i.e. the leaves

getLeaf = function(node,tree,leafs)
{
  leafs_fun = getLeafGeneric()
  leafs_fun(node,tree,leafs)
}
########
