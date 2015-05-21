GA_mutate_select2 <-
function(pass.pop, n_species, n_ones, n_zero, damp_val) {
  pops <- c(pass.pop)
  tmp <- .C("GA_mutate_select2", pops=as.double(pops), n_species=as.integer(n_species),
            n_ones=as.integer(n_ones), n_zero=as.integer(n_zero), damp_val=as.double(damp_val)) 
  return(tmp$pops)
}
