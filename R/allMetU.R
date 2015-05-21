allMetU <-
function(A) {
    if ((sum(A == 0) + sum(A == 1)) != length(A)) {
      fw <- ifelse(A > 0, 1, 0)
    } else {
      fw <- A
    }
    A.graph <- graph.adjacency(A, mode="directed", diag=FALSE)
    A.graph2 <- graph.adjacency(A, mode="undirected", diag=FALSE)
    diag(fw) <- 0
    S <- ncol(fw)
    L <- sum(fw)
    # connectance
    connect <- L / (S * S)
    # no. basal
    basal <- sum(apply(fw, 2, sum, na.rm=TRUE) == 0, na.rm=TRUE)
    # no. int.
    int <- sum((apply(fw, 1, sum, na.rm=TRUE) != 0) &
               (apply(fw, 2, sum, na.rm=TRUE) != 0))
    # no. top
    top <- sum(apply(fw, 1, sum, na.rm=TRUE) == 0, na.rm=TRUE)
    # generality
    gen <- L / (top + int)
    # vulnerability
    vul <- L / (basal + int)
    # standard deviation of generality
    genk <- (S / L) * apply(fw, 2, sum, na.rm=TRUE)
    sdGen <- sd(genk, na.rm=TRUE)
    # standard deviation of vulnerability
    vulk <- (S / L) * apply(fw, 1, sum, na.rm=TRUE)
    sdVul <- sd(vulk, na.rm=TRUE)
    # mean path length
    mean.pl <- average.path.length(A.graph)
    # mean trophic level
    mean.tl <- mean(TrophInd(fw)$TL, na.rm=TRUE)
    # mean omnivory
    mean.oi <- mean(TrophInd(fw)$OI, na.rm=TRUE)
    # modularity
    fg.mod <- modularity(fastgreedy.community(A.graph2))
    sg.mod <- 0#modularity(spinglass.community(A.graph))
    le.mod <- modularity(leading.eigenvector.community(A.graph2))
    rw.mod <- modularity(walktrap.community(A.graph))
    # clustering
    clust <- clustering_w(fw, measure="bi")
    # transitivity
    trans <- gtrans(fw, mode="digraph")
	return(list(connect=connect, basal=basal, int=int, top=top, gen=gen, 
	            vul=vul, sdGen=sdGen, sdVul=sdVul, mean.pl=mean.pl,
	            mean.tl=mean.tl, mean.oi=mean.oi, fg.mod=fg.mod, sg.mod=sg.mod,
	            le.mod=le.mod, rw.mod=rw.mod, clust=clust, trans=trans))
}
