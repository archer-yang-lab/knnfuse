#' Estimate the Number of Components in a Multivariate Normal Location Mixture Model.
#'
#' Estimate the order of a finite mixture of multivariate normal densities
#' with respect to the mean parameter, whose variance-covariance matrices
#' are common but potentially unknown. 
#'
#' @param y n by D matrix consisting of the data, where n is the sample size 
#'          and D is the dimension.
#' @param K Upper bound on the true number of components. 
#'          If \code{K} is \code{NULL}, at least one of \code{theta} and
#'          \code{pii} must be non-\code{NULL}, and K is inferred from their
#'          number of columns.
#' @param lambdas Vector of tuning parameter values. 
#' @param sigma   D by D matrix, which is the starting value for the common
#'                variance-covariance matrix. If \code{NULL}, \code{arbSigma}
#'                must be \code{TRUE}, and in this case the starting value is
#'                set to be the sample variance-covariance of the data.
#' @param arbSigma Equals \code{TRUE} if the common variance-covariance matrix should
#'                 be estimated, and FALSE if it should stay fixed, and equal to
#'                 \code{sigma}.
#' @param ... Additional control parameters. See the \strong{Details} section.
#'
#' @return An object with S3 classes \code{gsf} and \code{normalLocGsf}, 
#'         consisting of a list with the estimates produced for every tuning
#'         parameter in \code{lambdas}. 
#'
#' @details The following is a list of additional control parameters.
#'
#'  \describe{
#'   \item{\code{mu}}{D by K matrix of starting values where each column is the mean
#'                      vector for one component. If \code{theta=NULL}, the starting 
#'                      values are chosen using the procedure of Benaglia et al. (2009).}
#'   \item{\code{pii}}{Vector of size K whose elements must sum to 1, consisting of 
#'                      the starting values for the mixing proportions.  
#'                      If \code{NULL}, it will be set to a discrete 
#'                      uniform distribution with K support points.}
#'   \item{\code{penalty}}{Choice of penalty, which may be "\code{SCAD}", "\code{MCP}", 
#'                          "\code{SCAD-LLA}", "\code{MCP-LLA}" or "\code{ADAPTIVE-LASSO}". 
#'                          Default is "\code{SCAD}".}
#'   \item{\code{C}}{Tuning parameter for penalizing the mixing proportions.}
#'   \item{\code{a}}{Tuning parameter for the SCAD or MCP penalty. Default is \code{3.7}.}
#'   \item{\code{convMem}}{Convergence criterion for the modified EM algorithm.}
#'   \item{\code{maxMem}}{Maximum number of iterations of the modified EM algorithm.}
#'   \item{\code{verbose}}{If \code{TRUE}, print updates while the function is running.}}
#' 
#'
#' @references 
#'  Manole, T., Khalili, A. 2019. "Estimating the Number of Components in Finite Mixture Models 
#'  via the Group-Sort-Fuse Procedure".
#'
#'  Benaglia, T., Chauveau, D., Hunter, D., Young, D. 2009. "mixtools: An R package for analyzing
#'  finite mixture models". Journal of Statistical Software. 32(6): 1-29.
#'
#' @examples
#'  # Example 1: Seeds Data.
#'    data(seeds) 
#'    y <- cbind(seeds[,2], seeds[,6])
#'    set.seed(1)
#'    out <- normalLocOrder(y, K=12, lambdas=seq(0.1, 1.7, 0.3), maxPgd=200, maxMem=500)
#'    tuning <- bicTuning(y, out)
#'    plot(out, gg=FALSE, eta=TRUE, vlines=TRUE, opt=tuning$result$lambda)
#'   
#'  # Example 2: Old Faithful Data.
#'    data(faithful)
#'    set.seed(1)
#'    out <- normalLocOrder(faithful, K=10, 
#'              lambdas=c(0.1, 0.25, 0.5, 0.75, 1.0, 2), penalty="MCP-LLA", a=2, maxPgd=200, maxMem=500)
#' 
#'    # Requires ggplot2.
#'    plot(out, gg=TRUE, eta=FALSE)
normalLocOrder <- function (y, m, lambdas, graphtype = "MNN", K = NULL, sigma = NULL, arbSigma = TRUE, ...) {
  input <- list(...)
  
  if (K==1) stop("Error: K must be at least 2")
  
  mu  <- input$mu
  pii <- input$pii
  
  if(graphtype == "MNN") graphtype = 1
  else if (graphtype == "GSF") graphtype = 2
  else if (graphtype == "MST") graphtype = 3
  else if (graphtype == "SE") graphtype = 4
  
  if (is.null(input[["penalty"]])) penalty <- "SCAD"
  else penalty <- input[["penalty"]]
  
  if (is.null(input[["C"]])) C <- 3
  else C <- input[["C"]]
  
  if (is.null(input[["a"]])) a <- 3.7
  else a <- input[["a"]]
  
  if (is.null(input[["convADMM"]])) delta <- 1e-5
  else delta <- input[["convADMM"]]
  
  if (is.null(input[["convMem"]])) epsilon <- 1e-8
  else epsilon <- input[["convMem"]]
  
  if (is.null(input[["maxadmm"]])) maxadmm <- 1000
  else maxMem <- input[["maxadmm"]]
  
  if (is.null(input[["maxNR"]])) maxNR <- 100
  else maxMem <- input[["maxNR"]]
  
  if (is.null(input[["maxMem"]])) maxMem <- 2500
  else maxMem <- input[["maxMem"]]
  
  if (is.null(input[["verbose"]])) verbose <- T
  else verbose <- input[["verbose"]]
  
  if (is.null(input[["dynamic"]])) dynamic <- F
  else dynamic <- input[["dynamic"]]
  
  if (is.data.frame(y)) {
    y <- as.matrix(y)
  }
  
  # .validateParams(y, mu, "mu", pii, K)
  # .validateOther(maxMem, maxPgd, lambdas, penalty, a, C, epsilon, delta)
  # .validateNormalLoc(y, mu, sigma, arbSigma)
  
  lambdas <- sort(lambdas)
  
  if (!is.null(mu) && !is.null(pii)) { #when the starting values are pre-specified
    if (arbSigma) sigma <- cov(y)
    
    out <- .myEm(y, graphtype, m, mu, sigma, pii, arbSigma, -1, 1, maxadmm, maxNR, C, a,
                 .penCode(penalty), lambdas, epsilon, delta, maxMem, verbose)
    
  } else if (!is.null(mu) && is.null(pii)) {
    K <- ncol(mu)
    if (arbSigma) sigma <- cov(y)
    
    out <- .myEm(y, graphtype, m, mu, sigma, rep(1.0/K, K), arbSigma,-1, 1, maxadmm, maxNR, C, a,
                 .penCode(penalty), lambdas, epsilon, delta, maxMem,  verbose)
    
  } else {
    if (is.null(pii)) {
      pii <- rep(1.0/K, K)
      
    } else {
      K <- length(pii)
    }
    
    if (arbSigma) sigma <- cov(y)
    
    # Compute starting values.
    means <- apply(y, 1, mean)
    yOrdered <- y[order(means), ]
    yBinned  <- list()
    
    n <- nrow(y)
    
    for (j in 1:K) {
      yBinned[[j]] <- yOrdered[max(1, floor((j - 1) * n/K)):ceiling(j * n/K), ]
    }
    
    hypTheta <- lapply(1:K, function(i) apply(yBinned[[i]], 2, mean))
    
    theta <- matrix(NA, ncol(y), K)
    for(k in 1:K) {
      theta[,k] <- as.vector(rmvnorm(1, mu = as.vector(hypTheta[[k]]), Sigma = sigma))
    }
    
    out <- .myEm(y, graphtype, m, theta, sigma, pii, arbSigma, -1, 1, maxadmm, maxNR, C, a,
                 .penCode(penalty), lambdas, epsilon, delta, maxMem,  verbose, dynamic)
  }
  class(out) <- c("normalLocGsf", "gsf")
  out
}

#' Estimate the Number of Components in a Multivariate Normal Location Mixture Model.
#'
#' Estimate the order of a finite mixture of multivariate normal densities
#' with respect to the mean parameter, whose variance-covariance matrices
#' are common but potentially unknown. 
#'
#' @param y n by D matrix consisting of the data, where n is the sample size 
#'          and D is the dimension.
#' @param K Upper bound on the true number of components. 
#'          If \code{K} is \code{NULL}, at least one of \code{theta} and
#'          \code{pii} must be non-\code{NULL}, and K is inferred from their
#'          number of columns.
#' @param lambdas Vector of tuning parameter values. 
#' @param sigma   D by D matrix, which is the starting value for the common
#'                variance-covariance matrix. If \code{NULL}, \code{arbSigma}
#'                must be \code{TRUE}, and in this case the starting value is
#'                set to be the sample variance-covariance of the data.
#' @param arbSigma Equals \code{TRUE} if the common variance-covariance matrix should
#'                 be estimated, and FALSE if it should stay fixed, and equal to
#'                 \code{sigma}.
#' @param ... Additional control parameters. See the \strong{Details} section.
#'
#' @return An object with S3 classes \code{gsf} and \code{normalLocGsf}, 
#'         consisting of a list with the estimates produced for every tuning
#'         parameter in \code{lambdas}. 
#'
#' @details The following is a list of additional control parameters.
#'
#'  \describe{
#'   \item{\code{mu}}{D by K matrix of starting values where each column is the mean
#'                      vector for one component. If \code{theta=NULL}, the starting 
#'                      values are chosen using the procedure of Benaglia et al. (2009).}
#'   \item{\code{pii}}{Vector of size K whose elements must sum to 1, consisting of 
#'                      the starting values for the mixing proportions.  
#'                      If \code{NULL}, it will be set to a discrete 
#'                      uniform distribution with K support points.}
#'   \item{\code{penalty}}{Choice of penalty, which may be "\code{SCAD}", "\code{MCP}", 
#'                          "\code{SCAD-LLA}", "\code{MCP-LLA}" or "\code{ADAPTIVE-LASSO}". 
#'                          Default is "\code{SCAD}".}
#'   \item{\code{C}}{Tuning parameter for penalizing the mixing proportions.}
#'   \item{\code{a}}{Tuning parameter for the SCAD or MCP penalty. Default is \code{3.7}.}
#'   \item{\code{convMem}}{Convergence criterion for the modified EM algorithm.}
#'   \item{\code{maxMem}}{Maximum number of iterations of the modified EM algorithm.}
#'   \item{\code{verbose}}{If \code{TRUE}, print updates while the function is running.}}
#' 
#'
#' @references 
#'  Manole, T., Khalili, A. 2019. "Estimating the Number of Components in Finite Mixture Models 
#'  via the Group-Sort-Fuse Procedure".
#'
#'  Benaglia, T., Chauveau, D., Hunter, D., Young, D. 2009. "mixtools: An R package for analyzing
#'  finite mixture models". Journal of Statistical Software. 32(6): 1-29.
#'
#' @examples
#'  # Example 1: Seeds Data.
#'    data(seeds) 
#'    y <- cbind(seeds[,2], seeds[,6])
#'    set.seed(1)
#'    out <- normalLocOrder(y, K=12, lambdas=seq(0.1, 1.7, 0.3), maxPgd=200, maxMem=500)
#'    tuning <- bicTuning(y, out)
#'    plot(out, gg=FALSE, eta=TRUE, vlines=TRUE, opt=tuning$result$lambda)
#'   
#'  # Example 2: Old Faithful Data.
#'    data(faithful)
#'    set.seed(1)
#'    out <- normalLocOrder(faithful, K=10, 
#'              lambdas=c(0.1, 0.25, 0.5, 0.75, 1.0, 2), penalty="MCP-LLA", a=2, maxPgd=200, maxMem=500)
#' 
#'    # Requires ggplot2.
#'    plot(out, gg=TRUE, eta=FALSE)
tLocOrder <- function (y, m, lambdas, df = NULL, graphtype = "MNN", K = NULL, sigma = NULL, arbSigma = TRUE, ...) {
  input <- list(...)
  
  if (K==1) stop("Error: K must be at least 2")
  
  mu  <- input$mu
  pii <- input$pii
  
  if(graphtype == "MNN") graphtype = 1
  else if (graphtype == "GSF") graphtype = 2
  else if (graphtype == "MST") graphtype = 3
  else if (graphtype == "SE") graphtype = 4
  
  if (is.null(input[["df"]])) df <- rep(1,K)
  else df <- input[["df"]]
  
  if (is.null(input[["penalty"]])) penalty <- "SCAD"
  else penalty <- input[["penalty"]]
  
  if (is.null(input[["C"]])) C <- 3
  else C <- input[["C"]]
  
  if (is.null(input[["a"]])) a <- 3.7
  else a <- input[["a"]]
  
  if (is.null(input[["convADMM"]])) delta <- 1e-5
  else delta <- input[["convADMM"]]
  
  if (is.null(input[["convMem"]])) epsilon <- 1e-8
  else epsilon <- input[["convMem"]]
  
  if (is.null(input[["maxadmm"]])) maxadmm <- 1000
  else maxMem <- input[["maxadmm"]]
  
  if (is.null(input[["maxNR"]])) maxNR <- 100
  else maxMem <- input[["maxNR"]]
  
  if (is.null(input[["maxMem"]])) maxMem <- 2500
  else maxMem <- input[["maxMem"]]
  
  if (is.null(input[["verbose"]])) verbose <- T
  else verbose <- input[["verbose"]]
  
  if (is.null(input[["dynamic"]])) dynamic <- F
  else dynamic <- input[["dynamic"]]
  
  if (is.data.frame(y)) {
    y <- as.matrix(y)
  }
  
  # .validateParams(y, mu, "mu", pii, K)
  # .validateOther(maxMem, maxPgd, lambdas, penalty, a, C, epsilon, delta)
  # .validateNormalLoc(y, mu, sigma, arbSigma)
  
  lambdas <- sort(lambdas)
  
  if (!is.null(mu) && !is.null(pii)) { #when the starting values are pre-specified
    if (arbSigma) sigma <- cov(y)
    
    out <- .myEm2(y, df, graphtype, m, mu, sigma, pii, arbSigma, -1, 1, maxadmm, maxNR, C, a,
                 .penCode(penalty), lambdas, epsilon, delta, maxMem, verbose)
    
  } else if (!is.null(mu) && is.null(pii)) {
    K <- ncol(mu)
    if (arbSigma) sigma <- cov(y)
    
    out <- .myEm2(y, df, graphtype, m, mu, sigma, rep(1.0/K, K), arbSigma,-1, 1, maxadmm, maxNR, C, a,
                 .penCode(penalty), lambdas, epsilon, delta, maxMem,  verbose)
    
  } else {
    if (is.null(pii)) {
      pii <- rep(1.0/K, K)
      
    } else {
      K <- length(pii)
    }
    
    if (arbSigma) sigma <- cov(y)
    
    # Compute starting values.
    means <- apply(y, 1, mean)
    yOrdered <- y[order(means), ]
    yBinned  <- list()
    
    n <- nrow(y)
    
    for (j in 1:K) {
      yBinned[[j]] <- yOrdered[max(1, floor((j - 1) * n/K)):ceiling(j * n/K), ]
    }
    
    hypTheta <- lapply(1:K, function(i) apply(yBinned[[i]], 2, mean))
    
    theta <- matrix(NA, ncol(y), K)
    for(k in 1:K) {
      theta[,k] <- as.vector(rmvnorm(1, mu = as.vector(hypTheta[[k]]), Sigma = sigma))
    }
    
    out <- .myEm2(y, df, graphtype, m, theta, sigma, pii, arbSigma, -1, 1, maxadmm, maxNR, C, a,
                 .penCode(penalty), lambdas, epsilon, delta, maxMem,  verbose, dynamic)
  }
  class(out) <- c("normalLocGsf", "gsf")
  out
}

#' Estimate the Number of Components in a Multinomial Mixture Model.
#'
#' Estimate the order of a finite mixture of multinomial models with fixed and 
#' known number of trials.
#'
#' @param y n by D matrix consisting of the data, where n is the sample size 
#'          and D is the number of categories.
#'          The rows of \code{y} must sum to a constant value, 
#'          taken to be the number of trials. 
#' @param K Upper bound on the true number of components. 
#'          If \code{K} is \code{NULL}, at least one of the control parameters \code{theta} and
#'          \code{pii} must be non-\code{NULL}, and K is inferred from their
#'          number of columns.
#' @inheritParams normalLocOrder
#'
#' @return An object with S3 classes \code{gsf} and \code{multinomialGsf}, 
#'         consisting of a list with the estimates produced for every tuning
#'         parameter in \code{lambdas}. 
#'
#' @details The following is a list of additional control parameters.
#'
#' \itemize{
#'   \item{\code{theta}}    {D by K matrix of starting values where each column is the vector
#'                      of multinomial probabilities for one mixture component.              
#'                      The columns of \code{theta} should therefore 
#'                      sum to 1. If \code{theta=NULL}, the starting values are
#'                      chosen using the MCMC algorithm described by Grenier (2016).} 
#'   \item{\code{pii}}      {Vector of size K whose elements must sum to 1, consisting of 
#'                      the starting values for the mixing proportions.  
#'                      If \code{NULL}, it will be set to a discrete 
#'                      uniform distribution with K support points.}
#'   \item{\code{penalty}} {Choice of penalty, which may be \code{"SCAD"}, \code{"MCP"}, 
#'                          \code{"SCAD-LLA"}, \code{"MCP-LLA"} or \code{"ADAPTIVE-LASSO"}. 
#'                          Default is \code{"SCAD"}.}
#'   \item{\code{mcmcIter}} {Number of iterations for the starting value algorithm described 
#'                           by Grenier (2016).}
#'   \item{\code{C}}       {Tuning parameter for penalizing the mixing proportions.}
#'   \item{\code{a}}        {Tuning parameter for the SCAD or MCP penalty. Default is \code{3.7}.}
#'   \item{\code{convMem}}  {Convergence criterion for the modified EM algorithm.}
#'   \item{\code{maxMem}}   {Maximum number of iterations of the Modified EM algorithm.}
#'   \item{\code{maxNR}}   {Maximum number of iterations of the Newton-Ralphson algorithm.}
#'   \item{\code{verbose}}  {If \code{TRUE}, print updates while the function is running.}
#' }
#'
#' @references 
#'  Manole, T., Khalili, A. 2019. "Estimating the Number of Components in Finite Mixture Models 
#'  via the Group-Sort-Fuse Procedure".
#'
#'  Grenier, I. (2016) Bayesian Model Selection for Deep Exponential Families. 
#'  M.Sc. dissertation, McGill University Libraries.
#'
#' @examples 
#'  require(MM)
#'  data(pollen)
#'  set.seed(1)
#'  out <- multinomialOrder(pollen, K=12, lambdas=seq(0.1, 1.6, 0.2))
#'  tuning <- bicTuning(pollen, out)
#'  plot(out, eta=TRUE, gg=FALSE, opt=tuning$result$lambda)
multinomialOrder <- function (y, m, lambdas, graphtype = "MNN", K = NULL, ...) {
  input <- list(...)

  if (K==1) stop("Error: K must be at least 2")
  
  theta <- input$theta
  pii <- input$pii

  if(graphtype == "MNN") graphtype = 1
  else if (graphtype == "GSF") graphtype = 2
  else if (graphtype == "MST") graphtype = 3
  else if (graphtype == "SE") graphtype = 4
  
  if (is.null(input[["penalty"]])) penalty <- "SCAD"
  else penalty <- input[["penalty"]]

  if (is.null(input[["uBound"]])) uBound <- 0.1
  else uBound <- input[["uBound"]]
  
  if (is.null(input[["C"]])) C <- 3
  else C <- input[["C"]]

  if (is.null(input[["a"]])) a <- 3.7
  else a <- input[["a"]]

  if (is.null(input[["mcmcIter"]])) mcmcIter <- 50
  else mcmcIter <- input[["mcmcIter"]]

  if (is.null(input[["convPgd"]])) delta <- 1e-5
  else delta <- input[["convPgd"]]

  if (is.null(input[["convMem"]])) epsilon <- 1e-8
  else epsilon <- input[["convMem"]]
  
  if (is.null(input[["maxadmm"]])) maxadmm <- 500
  else maxMem <- input[["maxadmm"]]
  
  if (is.null(input[["maxNR"]])) maxNR <- 10
  else maxMem <- input[["maxNR"]]

  if (is.null(input[["maxMem"]])) maxMem <- 2500
  else maxMem <- input[["maxMem"]]

  if (is.null(input[["verbose"]])) verbose <- T
  else verbose <- input[["verbose"]]
  
  if (is.null(input[["dynamic"]])) dynamic <- F
  else dynamic <- input[["dynamic"]]

  if (is.data.frame(y)) {
    y <- as.matrix(y)
  }

  # .validateParams(y, theta, "theta", pii, K)
  # .validateOther(maxMem, maxPgd, lambdas, penalty, a, C, epsilon, delta)
  # .validateMultinomial(y, theta, mcmcIter)

  lambdas <- sort(lambdas)

  D <- ncol(y) - 1

  if(is.null(K)) {
    out <- .completeMultinomialCols(.myEm(y[, 1:D], graphtype, m, theta[1:D, ], diag(D), pii, FALSE, sum(y[1,]), 3, maxadmm, maxNR, C, a,
                .penCode(penalty), lambdas, epsilon, delta, maxMem, verbose))
  } else {
    n <- nrow(y)

    newPii <- array(1 / K, dim = c(n, K))
    newZ    <- t(apply(newPii, 1, function(x) rmultinom(1,1,x)))
    piiSamp <- matrix(0, mcmcIter, K)

    for (d in 1:mcmcIter) {
      nkWords  <- apply(newZ, 2, function(x) colSums(x*y))
      newTheta <- apply(nkWords, 2, function(x) gtools::rdirichlet(1,1 + x))

      pxGivTheta <- 0
      for (k in 1:K) {
        pxGivTheta <- pxGivTheta + t(apply(y, 1, function(x) dmultinom(x, sum(x), newTheta[, k])))
      }

      for (k in 1:K) {
        newPii[, k] <- t(apply(y, 1, function(x) dmultinom(x, sum(x), newTheta[, k]))) / sum(pxGivTheta)
      }

      newZ <- t(apply(newPii, 1, function(x) rmultinom(1, 1, x)))
      piiSamp[d, ] <- colSums(newZ)/n
    }

    if (is.null(pii)) {
      newPii <- as.matrix(piiSamp[mcmcIter,])
  
    } else {
      newPii <- as.matrix(pii)
    }

    out <- .completeMultinomialCols(.myEm(y[,1:D], graphtype, m, newTheta[1:D, ], diag(D), newPii, FALSE, sum(y[1,]), 3, maxadmm, maxNR, C, a,
                                          .penCode(penalty), lambdas, epsilon, delta, maxMem, verbose, dynamic))
  }
   
  class(out) <- c("multinomialGsf", "gsf")
  out
}


#' Tuning Parameter Selection via the Bayesian Information Criterion 
#'
#' Bayesian Information Criterion (BIC) for tuning parameter selection. 
#' @param y n by D matrix consisting of the data.
#' @param result A \code{gsf} object.
#'
#' @details The BIC selects the best tuning parameter out of the ones used in a
#'  \code{gsf} object by minimizing the following function
#'
#' \deqn{\textrm{BIC}(\lambda) = -2 l_n(\hat{\mathbf{\Psi}}_{\lambda}) + 
#' \textrm{df}(\hat{\mathbf{\Psi}}_{\lambda}) \log n}
#'
#' where \eqn{l_n} is the log-likelihood function, and \eqn{
#' \hat{\mathbf{\Psi}}_{\lambda}} is the set of estimated 
#' parameters \code{theta} and \code{pii} corresponding to the tuning parameter \eqn{\lambda}.
#' \eqn{\textrm{df}(\hat{\mathbf{\Psi}}_{\lambda})} denotes the degrees
#' of freedom of the estimates.
#'
#' @return A list containing the selected tuning parameter and 
#'         corresponding estimates, as well as a summary of the
#'         computed RBIC values, log-likelihood values, and corresponding
#'         orders of the estimates.
#'
#' @examples 
#'  require(MM)
#'  data(pollen)
#'  set.seed(1)
#'  out <- multinomialOrder(pollen, K=12, lambdas=seq(0.1, 1.6, 0.2))
#'  tuning <- bicTuning(pollen, out)
#'  plot(out, eta=TRUE, gg=FALSE, opt=tuning$result$lambda)
bicTuning <- function(y, result) {
  if (is.vector(y)) {
    n <- length(y)

  } else {
    n <- nrow(y)
  }

  if (is.data.frame(y)) {
    y <- as.matrix(y)
  }
 
  if (is.vector(y)) {
    y <- matrix(y, ncol=1)
  }
 
  l <- length(result)
  thetas <- list()

  for (i in 1:l) {
    thetas[[i]] <- switch(class(result)[1], 
                  "normalLocGsf" = result[[i]]$mu,
                  "multinomialGsf" = result[[i]]$theta,
                  "poissonGsf" = result[[i]]$theta,
                  "exponentialGsf" = result[[i]]$theta,
                   stop("Error: Class not recognized."))
  }  

  if (is.vector(thetas[[1]])) {
    D <- length(thetas[[1]])

  } else {
    D <- nrow(thetas[[1]])
  }

  out      <- list()
  estimate <- list() 

  lambdaVals <- c()
  order      <- c()
  df         <- c()
  loglik     <- c()
  bic        <- c()  

  minRbic  <- 9999999999999999 #.Machine$double.max
  minIndex <- -1 
  
  arbSigma <- ("sigma" %in% names(result[[1]])) && (class(result)[1] == "normalLocGsf")

  index <- switch(class(result)[1], 
            "normalLocGsf" = 1,
            "multinomialGsf" = 3,
            "poissonGsf" = 5,
            "exponentialGsf" = 6,
            stop("Error: Class not recognized."))

  for (i in 1:l) {
    lambdaVals[i] <- result[[i]]$lambda

    if(lambdaVals[i] == 0) next
  
    order[i]      <- result[[i]]$order
    df[i] <- switch(class(result)[1], 
            "normalLocGsf" = {
              if(arbSigma) {
                order[i] * (D+1) + 0.5 * D * (D+1) - 1
              
              } else {
                order[i] * (D+1) - 1
              }
            }, 
            "multinomialGsf" = D * order[i] - 1,
            "poissonGsf" = 2 * order[i] - 1,
            "exponentialGsf" = 2 * order[i] - 1,
            stop("Error: Class not recognized."))
  
    if (arbSigma) {
      loglik[i] <- .bicLogLik(y, as.matrix(thetas[[i]], nrow=D), result[[i]]$pii, result[[i]]$sigma, index)

    } else {
      loglik[i] <- .bicLogLik(y, as.matrix(thetas[[i]], nrow=D), result[[i]]$pii, diag(D), index)
    }

    if (!is.finite(loglik[i])) {
      stop("Error: Log-likelihood is infinite.")
    }

    bic[i] <- -2 * loglik[i] + df[i] * log(n)

    if (bic[i] < minRbic) {
      minRbic <- bic[i]
      minIndex  <- i
    } 
  }

  if (index == 1) {
    estimate[["mu"]] <- result[[minIndex]]$mu
    estimate[["sigma"]] <- result[[minIndex]]$sigma

  } else {
    estimate[["theta"]] <- result[[minIndex]]$theta
  }

  estimate[["pii"]]   <- result[[minIndex]]$pii
  estimate[["order"]] <- result[[minIndex]]$order
  estimate[["lambda"]]<- result[[minIndex]]$lambda
  estimate[["bic"]] <- bic[minIndex]

  out[["summary"]] = data.frame(lambdaVals, order, bic, loglik, df)
  out[["result"]] = estimate

  out 
}
 

