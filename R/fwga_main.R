fwga <-
function(n.pop=200, n.iter=100, fits=list(fit_tst, fit_robust), suggest=NULL, k=5, kp=0.5,
                        n.keep=round(n.pop / 5), damp.fact=1) {
  # check and format starting food web
  if (is.null(suggest)) {
    stop("An initial food web needs to be provided...", call.=FALSE)
  }
  suggest <- as.matrix(suggest)
  # calculate summary stats and create initial population
  S <- ncol(suggest)
  C <- sum(ifelse(suggest > 0, 1, 0)) / (S * S)
  start.pop <- vector("list", length=n.pop)
  if ((sum(suggest == 0) + sum(suggest == 1)) == length(suggest)) {
    start.pop[[1]] <- suggest * runif(n = (S * S), min=0, max=1)
  } else {
    start.pop[[1]] <- suggest
  }
  for (i in 2:n.pop) {
  	if ((sum(suggest == 0) + sum(suggest == 1)) == length(suggest)) {
      start.pop[[i]] <- runif(n = (S * S), min=0, max=1) * null2(suggest, iter=1)[, , 1]
    } else {
      start.pop[[i]] <- null3(suggest, iter=1)[[1]]
    }
  }
  names(start.pop) <- NULL
  # calculate fitness for starting population
  n.fit <- length(fits)
  pops.temp <- multiGA(start.pop, n.fit=n.fit, n.pop=n.pop, fits=fits, S=S)
  pops <- pops.temp$save.pop
  doms <- pops.temp$dom
  fitness <- vector("list", length=n.fit)
  for (i in 1:n.fit) {
    fitness[[i]] <- pops.temp$fitness[[i]][1:n.keep]
  }
  # selection and mutation step
  damper <- 1
  for (i in 2:n.iter) {
    pops.temp <- NULL
    start.pop <- vector("list", length=n.pop)
    start.pop[1:n.keep] <- pops[1:n.keep]
    k.prob <- vector("numeric", length=k)
    for (j in 1:k) {
      k.prob[j] <- kp * ((1 - kp) ^ (j - 1))
    }
    counter <- n.keep
    while (counter < n.pop) {
      # start k tournament selection to choose cand.pop
      pop.choice <- sample(1:counter, size=k, replace=FALSE)
      pop.choice <- pop.choice[order(doms[pop.choice])]
      pop.choice <- pop.choice[sample(1:k, size=1, replace=FALSE, prob=k.prob)]
      cand.pop <- start.pop[[pop.choice]]
      # start mutation, constraining row and col sums
      num.zero <- length(which(cand.pop == 0))
      num.pos <- length(which(cand.pop > 0))
      new.pop <- GA_mutate_select2(pass.pop=cand.pop, n_species=S, n_ones=num.pos,
                                   n_zero=num.zero, damp_val=damper)
      start.pop[[counter + 1]] <- matrix(new.pop, ncol=S)
      counter <- counter + 1
    }
    damper <- damper * damp.fact
    # calculate fitness for new population
    pops.temp <- multiGA(start.pop, n.fit=n.fit, n.pop=n.pop, fits=fits, S=S)
    pops <- pops.temp$save.pop
    doms <- pops.temp$dom
    for (i in 1:n.fit) {
      fitness[[i]] <- rbind(fitness[[i]], pops.temp$fitness[[i]][1:n.keep])
    }
  }
  return(list(pops=pops[1:n.keep], fitness=fitness, doms=doms[1:n.keep]))
}
