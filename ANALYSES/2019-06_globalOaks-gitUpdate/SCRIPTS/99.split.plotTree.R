split.plotTree <- function (tree, splits = NULL, file = NULL, nodeLabs = TRUE, ...)
{
    ef <- 0.037037037037
    if (!is.null(file))
        pdf(file, width = 8.5, height = 11)
    if (is.null(splits))
        splits <- (floor(0.5 * Ntip(tree)) + 0.5)/Ntip(tree)
    S <- matrix(c(0, splits, splits, 1 + 1/Ntip(tree)), length(splits) +
        1, 2)
    S <- cbind(S[, 1] + ef * (S[, 2] - S[, 1]), S[, 2] - ef *
        (S[, 2] - S[, 1]))
    for (i in nrow(S):1) {
        if (is.null(file) && i < nrow(S))
            par(ask = TRUE)
        plotTree(tree, ylim = Ntip(tree) * S[i, ], ...)
        if(nodeLabs) nodelabels(text = tree$node.label, adj=c(1,-0.5),
                    frame = 'none', cex = 0.5)
    }
    if (!is.null(file))
        oo <- dev.off()
}
