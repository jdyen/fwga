extract.metrics <-
function(obj) {
  connect <- sapply(obj, function(x) x$connect)
  basal <- sapply(obj, function(x) x$basal)
  int <- sapply(obj, function(x) x$int)
  top <- sapply(obj, function(x) x$top)
  gen <- sapply(obj, function(x) x$gen)
  vul <- sapply(obj, function(x) x$vul)
  sdGen <- sapply(obj, function(x) x$sdGen)
  sdVul <- sapply(obj, function(x) x$sdVul)
  mean.pl <- sapply(obj, function(x) x$mean.pl)
  mean.tl <- sapply(obj, function(x) x$mean.tl)
  mean.oi <- sapply(obj, function(x) x$mean.oi)
  fg.mod <- sapply(obj, function(x) x$fg.mod)
  sg.mod <- sapply(obj, function(x) x$sg.mod)
  le.mod <- sapply(obj, function(x) x$le.mod)
  rw.mod <- sapply(obj, function(x) x$rw.mod)
  clust <- sapply(obj, function(x) x$clust)
  trans <- sapply(obj, function(x) x$trans)
  return(list(connect=connect, basal=basal, int=int, top=top, gen=gen, vul=vul, sdGen=sdGen,
              sdVul=sdVul, mean.pl=mean.pl, mean.tl=mean.tl, mean.oi=mean.oi, fg.mod=fg.mod,
              sg.mod=sg.mod, le.mod=le.mod, rw.mod=rw.mod, clust=clust, trans=trans))
}
