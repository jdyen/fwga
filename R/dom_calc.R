dom_calc <-
function(fitness, n_fit, n_pop, dom){
  tmp <- .C("dom_calc", fitness=as.double(fitness), n_fit=as.integer(n_fit),
            n_pop=as.integer(n_pop), dom=as.integer(dom))
  return(list(fitness=tmp$fitness, dom=tmp$dom))
}
