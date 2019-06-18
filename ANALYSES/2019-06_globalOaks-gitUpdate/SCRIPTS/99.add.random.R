add.random.startNode <- function (tree, n = NULL,
                                  startNodeGrep = NA,
                                  excludeNodeGrep = NA,
                                  tips = NULL, edge.length = NULL, order = c("random",
    "input"))
{
    if (!inherits(tree, "phylo"))
        stop("tree should be an object of class \"phylo\".")
    randomPosn <- function(tree, node, excludeNode = NA) {
        nodesToDo <- c(node, descendants(tree, node, type = 'internal'))
        if(!is.na(excludeNode)) {
          nodesToDrop <- c(excludeNode, descendants(tree, excludeNode, type = 'internal'))
          nodesToDo <- setdiff(nodesToDo, nodesToDrop)
        }
        edgesToDo <- which(tree$edge[, 1] %in% nodesToDo)
        cum.edge <- cumsum(tree$edge.length[edgesToDo])
        index <- tree$edge[edgesToDo, 2]
        pos <- runif(1) * sum(tree$edge.length[edgesToDo])
        edge <- 1
        while (pos > cum.edge[edge]) edge <- edge + 1
        return(list(node = index[edge], posn = cum.edge[edge] -
            pos))
    }
#    debug(randomPosn)
    if (is.ultrametric(tree))
        um <- TRUE
    else um <- FALSE
    if (is.null(tips)) {
        if (is.null(n))
            n <- 1
        tips <- paste("t", length(tree$tip) + 1:n, sep = "")
    }
    else n <- length(tips)
    if (is.null(edge.length))
        if (!um)
            edge.length <- runif(n = n, min = min(tree$edge.length),
                max = max(tree$edge.length))
    if (order[1] == "random") {
        o <- sample(1:n)
        tips <- tips[o]
        if (!is.null(edge.length))
            edge.length <- edge.length[o]
    }
    for (i in 1:n) {
        treeNode <- findMRCA(tree, grep(startNodeGrep,tree$tip.label))
        if(is.na(excludeNodeGrep[1])) {
          where <- randomPosn(tree, treeNode)
          } else {
            excludeNode <- findMRCA(tree, grep(excludeNodeGrep, tree$tip.label))
            where <- randomPosn(tree, treeNode, excludeNode)
          } # close else
        if (is.null(edge.length))
            tree <- bind.tip(tree, tips[i], where = where$node,
                position = where$posn)
        else tree <- bind.tip(tree, tips[i], where = where$node,
            position = where$posn, edge.length = edge.length[i])
    } # close for(i in 1:n)
    return(tree)
}
