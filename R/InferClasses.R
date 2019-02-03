InferStructure <- function(X, type, control) {

  ## Define setting
  if (type == "steady-state") {
    is.directed <- FALSE
  }
  if (type == "time-course") {
    is.directed <- TRUE
  }

  ## Initial Matrix Estimate
  control$penalties <- NULL # enforce NULL penalties  
  res <- simone(X, type = type, control = control)
  adj <- getNetwork(res, selection = control$clusters.crit)$A
  if (sum(adj) == 0) {
    cat("\nEmpty initialization network. Let me try something else...")
    adj <- getNetwork(res)$A
    if (sum(adj) > 0) {
      cat(" That's better, I keep this one.")
    }
  }
  
  ## Remove Dust from analysis
  dust <- (rowSums(adj) + colSums(adj) == 0)
  if (sum(dust) == ncol(X)) {
    stop("Empty initialization network, choose another selection criterion.")
  }

  mySBM_collection <- blockmodels::BM_bernoulli(
    membership_type = ifelse(is.directed, "SBM", "SBM_sym"),
    adj             = adj[!dust, !dust],
    verbosity       = 0,
    plotting        = ""
  )
  mySBM_collection$estimate()
  
  ## Retrieve best classification
  best_id <- which.max(mySBM_collection$ICL)
  Z       <- mySBM_collection$memberships[[best_id]]$Z
  cl      <- apply(Z, 1, which.max)
  PI      <- mySBM_collection$model_parameters[[best_id]]$pi

  clusters <- rep(0,ncol(X))
  clusters[!dust] <- cl

  ## Weight structure
  weights  <- matrix(1000, ncol(X), ncol(X))
  weights[!dust, !dust] <- 1 - Z %*% PI %*% t(Z)

  list(
    clusters = clusters,
    weights  = weights
  )
}
