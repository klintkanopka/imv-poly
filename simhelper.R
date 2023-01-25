sim.fun <- function(b, b_fixed=0, dgm='graded', dam=NULL, n_i=25, n_p=500){
  require(mirt)

  if (is.null(dam)){
    dam <- dgm
  }

  make.resp <- function(b=b, b_grid=b_fixed, type=dgm, n_items=n_i, n_resp=n_p){

    thetas <- rnorm(n_resp)

    if(type == 'graded'){
      pars <- sort(c(b, b_grid), decreasing=TRUE)
    } else if(type == 'gpcm'){
      pars <- sort(c(b, b_grid), decreasing=TRUE)
    }

    par_mat <- c(runif(n_items, min=-1.5, max=1.5), pars[1])

    for (i in 2:length(pars)){
      par_mat <- c(par_mat, rep(NA, n_items), pars[i])
    }

    par_mat <- matrix(par_mat, ncol=length(pars))

    resp <- simdata(a=rep(1, n_items+1),
                    d=par_mat,
                    itemtype = c(rep('dich', n_items), type),
                    Theta=thetas)

    resp <- data.frame(resp)
    resp[[paste0('Item_', n_items+2)]] <- ifelse(resp[,ncol(resp)] == 0, 0, 1)
    return(resp)
  }

  mm <- mirt.model(paste0("F = 1-", n_i+2,
                          "\nFIXED = (", n_i+1,
                          ", a1)\nSTART = (", n_i+1,
                          ", a1, 1.0)"))
  resp <- make.resp(b)
  m <- mirt(resp, model=mm, itemtype = c(rep('Rasch', n_i), dam, 'Rasch'))
  th <- fscores(m)
  prob <- data.frame(probtrace(m, th)) # creates matrix of estimated probabilities

  # omega_0

  r_0 <- data.frame(resp=resp[[paste0('Item_', n_i+2)]])
  r_0$ctt <- rep(mean(r_0$resp), n_p)
  r_0$irt <- prob[[paste0("Item_", n_i+2, ".P.1")]]

  y <- data.frame(resp = resp[[paste0('Item_', n_i+1)]])

  pctt.tab <- c()
  for (i in 0:(length(unique(y$resp))-1)){
    pctt.tab <- c(pctt.tab, sum(y$resp == i)/length(y$resp))
  }

  for (i in 1:length(unique(y$resp))){
    y[[paste0('ctt', i-1)]] <- pctt.tab[i]
    y[[paste0('irt', i-1)]] <- prob[[paste0('Item_', n_i+1, '.P.', i)]]
  }

  out <- c(imv(r_0, 'ctt', 'irt'),
           imv_c(y, pctt.tab, 'ctt', 'irt'),
           imv_t(y, pctt.tab, 'ctt', 'irt'))

  return(out)
}


sim.fun.2 <- function(b1, b2, b_fixed=0, dgm='graded', dam=NULL, n_i=25, n_p=500){
  require(mirt)

  if (is.null(dam)){
    dam <- dgm
  }

  make.resp <- function(b1=b1, b2=b2, b_grid=b_fixed, type=dgm, n_items=n_i, n_resp=n_p){

    thetas <- rnorm(n_resp)

    if(type == 'graded'){
      pars <- sort(c(b1, b2, b_grid), decreasing=TRUE)
    } else if(type == 'gpcm'){
      pars <- sort(c(b1, b2, b_grid), decreasing=TRUE)
    }

    par_mat <- c(runif(n_items, min=-1.5, max=1.5), pars[1])

    for (i in 2:length(pars)){
      par_mat <- c(par_mat, rep(NA, n_items), pars[i])
    }

    par_mat <- matrix(par_mat, ncol=length(pars))

    resp <- simdata(a=rep(1, n_items+1),
                    d=par_mat,
                    itemtype = c(rep('dich', n_items), type),
                    Theta=thetas)

    resp <- data.frame(resp)
    resp[[paste0('Item_', n_items+2)]] <- ifelse(resp[,ncol(resp)] == 0, 0, 1)
    return(resp)
  }

  mm <- mirt.model(paste0("F = 1-", n_i+2,
                          "\nFIXED = (", n_i+1,
                          ", a1)\nSTART = (", n_i+1,
                          ", a1, 1.0)"))
  resp <- make.resp(b)
  m <- mirt(resp, model=mm, itemtype = c(rep('Rasch', n_i), dam, 'Rasch'))
  th <- fscores(m)
  prob <- data.frame(probtrace(m, th)) # creates matrix of estimated probabilities

  # omega_0

  r_0 <- data.frame(resp=resp[[paste0('Item_', n_i+2)]])
  r_0$ctt <- rep(mean(r_0$resp), n_p)
  r_0$irt <- prob[[paste0("Item_", n_i+2, ".P.1")]]

  y <- data.frame(resp = resp[[paste0('Item_', n_i+1)]])

  pctt.tab <- c()
  for (i in 0:(length(unique(y$resp))-1)){
    pctt.tab <- c(pctt.tab, sum(y$resp == i)/length(y$resp))
  }

  for (i in 1:length(unique(y$resp))){
    y[[paste0('ctt', i-1)]] <- pctt.tab[i]
    y[[paste0('irt', i-1)]] <- prob[[paste0('Item_', n_i+1, '.P.', i)]]
  }

  out <- c(imv(r_0, 'ctt', 'irt'),
           imv_c(y, pctt.tab, 'ctt', 'irt'),
           imv_t(y, pctt.tab, 'ctt', 'irt'))

  return(out)
}

prep.results <- function(control, output){
  control$omega_0 <- sim_out[1,]
  control$omega_c <- sim_out[2,]
  control$omega_t <- sim_out[3,]
  return(control)
}
