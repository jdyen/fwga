multiGA <-
function(start.pop, n.fit, n.pop, fits, S=S) {
  fitness <- vector("numeric", length=(n.pop * n.fit))
  for (i in 1:n.fit) {
    fitness[((i - 1) * n.pop + 1):(i * n.pop)] <- sapply(start.pop, fits[[i]], S=S)
  }
  dom <- vector(mode='integer', length=n.pop)
  GA_out <- dom_calc(fitness=fitness, n_fit=n.fit, n_pop=n.pop, dom=dom)
  dom <- GA_out$dom
  save.pop <- start.pop[order(dom)]
  fitness.store <- GA_out$fitness
  fitness <- vector("list", length=n.fit)
  for (i in 1:n.fit) {
    fitness[[i]] <- fitness.store[((i - 1) * n.pop + 1):(i * n.pop)][order(dom)]
  }
  return(list(save.pop=save.pop, dom=dom, fitness=fitness))
}
