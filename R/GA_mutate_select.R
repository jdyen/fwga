GA_mutate_select <-
function(pass.pop, dom, n_pop, n_keep, n_species, k_val, k_prob) {
  pops <- unlist(pass.pop)
  pops <- c(pops, rep(0, ((n_pop - n_keep) * (n_species * n_species))))
  tmp <- .C("GA_mutate_select", pops=as.double(pops), dom=as.integer(dom),
            n_pop=as.integer(n_pop), n_keep=as.integer(n_keep),
            n_species=as.integer(n_species), k_val=as.integer(k_val),
            k_prob=as.double(k_prob))
  new.pop <- vector("list", length=n_pop)
  for (i in seq(along=new.pop)) {
    new.pop[[i]] <- matrix(tmp$pops[((i - 1) * (n_species * n_species) + 1):(i *
                      (n_species * n_species))], ncol=n_species)
  }
  return(list(pops=new.pop))
}
