#' Original RobNorm Function
#'
#' To robustly normalize expression data (Author: Meng Wang, https://github.com/mwgrassgreen/RobNorm).
#'
#' @param X.0 The expression matrix in log scale.
#' @param gamma.0 The density exponent parameter gamma, in practice, taking gamma.0 = 0.5 or 1.
#' @param tol The tolerance for interations (default: 10^(-4)).
#' @param step The step limit (default: 50).
#'
#' @return Normalized expression data
#' 
#' @keywords internal
#'
RobNorm <- function (X.0, gamma.0=0.5, tol=10^(-4), step=200) {

  # to select the rows that are nonmissing in at least half of the samples
  id.norm <- rownames(X.0)[rowSums(!is.na(X.0)) >= (ncol(X.0)/2)]
  X.1 <- X.0[id.norm,]

  # to construct the standard samples from row medians
  x.stand <-  apply(X.1, 1, stats::median)
  names(x.stand) <- rownames(X.1)

  X <- cbind(x.stand, X.1)
  I <- nrow(X)
  J <- ncol(X)

  # to initialize parameters
  Nu.0 <- apply(X, 2, function(x) stats::median(x - x.stand, na.rm=TRUE))
  Mu.0 <- apply(X-outer(rep(1, I), Nu.0), 1, mean, na.rm=TRUE)
  sigma2.0 <- apply( (X-outer(rep(1, I), Nu.0)-outer(Mu.0, rep(1, J)))^2, 1, mean, na.rm=TRUE)

  mean.0.mx <- Mu.0 %*% t(rep(1, J)) + rep(1, I) %*% t(Nu.0)
  var.0.mx <- sigma2.0 %*% t(rep(1, J))
  den.0 <- stats::dnorm(X, mean.0.mx, sqrt(var.0.mx))
  dum <- den.0^gamma.0
  dum[is.na(den.0)] <- NA
  w.mx <- dum/(rowSums(dum, na.rm=TRUE) %*% t(rep(1,J)))
  w.size <- rowSums(dum, na.rm=TRUE)

  # to obtain the initial divergence
  divergence.0 <- sum( w.size * ((gamma.0+1) * rowSums(w.mx * (X - Mu.0 %*% t(rep(1, J)) - rep(1, I) %*% t(Nu.0))^2, na.rm=TRUE)  /sigma2.0 + log(2*pi*sigma2.0/(gamma.0+1))  ) )

  # iterations
  para.diff.int <- c()
  divergence.int <- divergence.0
  int <- 0
  flag <- FALSE
  while ( flag == FALSE ) {

    # to update mu and sigma-squared
    Mu.new <- rowSums((X - rep(1,I) %*% t(Nu.0)) * w.mx, na.rm=TRUE)
    sigma2.new <- (rowSums((X - rep(1,I) %*% t(Nu.0) - Mu.new %*% t(rep(1, J)))^2 * w.mx, na.rm=TRUE))*(1+gamma.0)

    # to break the iterations if some sigma-squared is too small (locally trapped)
    if(min(sigma2.new, na.rm=TRUE) < 10^(-10)) {
      Mu.0 <- NULL
      Nu.0 <- NULL
      sigma2.0 <-NULL
      divergence.int <- NULL
      para.diff.int <- NULL
      stop("the pre-chosen gamma is large, to choose a smaller nonnegative gamma")
    }

    # to udpate nu
    weight.nu <- w.mx * ((w.size /(sigma2.new/(1+gamma.0))) %*% t(rep(1, J)))
    Nu.new <- colSums(weight.nu * (X - Mu.new %*% t(rep(1,J))), na.rm=TRUE)/colSums(weight.nu, na.rm=TRUE)
    Mu.new <- Mu.new + Nu.new[1]
    Nu.new <- Nu.new - Nu.new[1]

    para.diff.new <- sum(abs(Mu.new - Mu.0), na.rm=TRUE) + sum(abs(Nu.new - Nu.0), na.rm=TRUE) + sum(abs(sigma2.new - sigma2.0), na.rm=TRUE)
    para.diff.int <- c(para.diff.int, para.diff.new)

    # to update the divergence
    divergence.new <- sum( w.size * ( (gamma.0+1) * rowSums(w.mx * (X - Mu.new %*% t(rep(1, J)) - rep(1, I) %*% t(Nu.new))^2, na.rm=TRUE)  /sigma2.new + log(2*pi*sigma2.new/(gamma.0+1))  ) , na.rm=TRUE)
    divergence.int <- c(divergence.int, divergence.new)

    int <- int + 1
    flag <- (para.diff.new < tol) | (int >= step)

    if (flag == TRUE) break;

    Mu.0 <- Mu.new
    Nu.0 <- Nu.new
    sigma2.0 <- sigma2.new
    mean.0.mx <- Mu.0 %*% t(rep(1, J)) + rep(1, I) %*% t(Nu.0)
    var.0.mx <- sigma2.0 %*% t(rep(1, J))
    den.0 <- stats::dnorm(X, mean.0.mx, sqrt(var.0.mx))

    dum <- den.0^gamma.0
    dum[is.na(den.0)]  <- NA
    w.mx <- dum/ (rowSums(dum, na.rm=TRUE) %*% t(rep(1,J)))
    w.size <- rowSums(dum, na.rm=TRUE)

  }

  nu.rob.est <- Nu.0[-1]
  X.0.rob <- X.0 - outer(rep(1, nrow(X.0)), nu.rob.est)

  return(list(norm.data=X.0.rob, id.norm=id.norm, stand.s=x.stand, nu.est=nu.rob.est, mu.est=Mu.0, sigma2.est=sigma2.0, w.mx=w.mx, divergence=divergence.int, para.diff=para.diff.int))
}
