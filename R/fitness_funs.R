fit_rand <-
function(x, S) { # pseudo-random
  return(5)
}


fit_tst <-
function(x, S) { # total throughput
  return(-sum(x))
}


fit_robust <-
function(x, S) { # robustness
  diag(x) <- 0
  sec.ext <- 0
  prim.prod <- which(apply(x, 2, sum) == 0)
  robust50 <- 0
  i <- 1
  while(sec.ext < (S / 2)) {
    x[i, ] <- 0
    sec.ext <- sec.ext + sum(apply(x, 2, sum)[-c(i, prim.prod)] == 0) + 1
    robust50 <- robust50 + 1
    i <- i + 1
  }
  robust50 <- -robust50 / S
  return(robust50)
}


fit_se2  <-
function(x, S) { # secondary extinctions
  diag(x) <- 0
  sec.ext <- 0
  prim.prod <- which(apply(x, 2, sum) == 0)
  for (i in 1:S) {
    x.store <- x
    x[i, ] <- 0
    sec.ext <- sec.ext + sum(apply(x, 2, sum)[-c(i, prim.prod)] == 0)
    x <- x.store
  }
  sec.ext <- sec.ext / S
  return(sec.ext)
}


fit_se <-
function(x, S) {
  fw <- ifelse(x > 0, 1, 0)
  n_sp <- S
  sec_e <- 0
  fw_test <- fw
  diag(fw_test) <- 0
  n_prim_prod <- sum(apply(fw_test, 2, sum) == 0)
  fw_in <- vector("numeric", length=(n_sp * n_sp))
  for (i in 1:n_sp) {
    fw_in[((i - 1) * n_sp + 1):(i * n_sp)] <- fw[i, ]
  }
  tmp <- .C("fitness_se", fw_in=as.double(fw_in), n_sp=as.integer(n_sp),
            sec_e=as.double(sec_e))
  for (i in 1:n_sp) {
    fw[i, ] <- fw_in[((i - 1) * n_sp + 1):(i * n_sp)]
  }
  sec_ext <- (tmp$sec_e - ((n_sp - 1) * n_prim_prod)) / n_sp
  return(sec_ext)
}

fit_wtst <- function(x, S) {
  fw <- ifelse(x > 0, 1, 0)
  diag(fw) <- 0
  TL <- TrophInd(Tij=fw)$TL
  wFW <- sweep(fw, 1, TL, function(x, y) x * (10 ^ y))
  if (is.nan(sum(wFW)) | (sum(wFW, na.rm=TRUE) == 0) | !is.finite(sum(wFW, na.rm=TRUE))) {
    out <- 0
  } else {
    out <- -log(sum(wFW, na.rm=TRUE))
  }
  return(out)
}