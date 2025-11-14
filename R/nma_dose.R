nma_dose <- function(TE, seTE, weights, studlab,
                     agent1, agent2, treat1, treat2,
                     narms,
                     Xd, D,
                     n,
                     level, reference) {
  
  m <- length(TE) # number of pairwise comparisons
  a <- length(unique(c(agent1, agent2))) # Number of agents
  #
  df1 <- 2 * sum(1 / narms) # sum of degrees of freedom per study
  
  #
  # Adjusted weights
  #
  W <- diag(weights, nrow = length(weights))
  #
  # Laplacian matrix and pseudoinverse of L
  #
  L <- t(Xd) %*% W %*% Xd
  #
  Lplus <- ginv(L) # Cov matrix of beta
  Lplus[is_zero(Lplus)] <- 0
  colnames(Lplus) <- colnames(L)
  rownames(Lplus) <- rownames(L)
  
  # R resistance distance (variance) matrix (n x n)
  #
  R <- matrix(0, nrow = a, ncol = a)
  for (i in 1:a) {
    for (j in 1:a) {
      R[i, j] <- Lplus[i, i] + Lplus[j, j] - 2 * Lplus[i, j]
    }
  }
  #
  # H matrix
  #
  H <- Xd %*% Lplus %*% t(Xd) %*% W
  
  labels <- colnames(Lplus)
  #
  # beta = dose-response effects for agents
  #
  beta <- as.vector(Lplus %*% t(Xd) %*% W %*% TE)
  names(beta) <- labels
  #
  # Adjust treatment effects
  if (beta[reference] != 0)
    beta <- beta - beta[reference]
  #
  se.beta <- sqrt(diag(Lplus))
  ci.beta <- ci(beta, se.beta, level = level)
  #
  # delta = treatment estimates for observed comparisons
  #
  delta <- as.vector(Xd %*% beta) # = H %*% TE
  se.delta <- unname(sqrt(diag(Xd %*% Lplus %*% t(Xd))))
  #
  # delta.all(.matrix) = all direct and indirect treatment estimates
  #
  B.full <- createB(ncol = n)
  X.all <- B.full %*% D
  colnames(X.all) <- colnames(D)
  #
  delta.all <- as.vector(X.all %*% beta)
  se.delta.all <- sqrt(diag(X.all %*% Lplus %*% t(X.all)))
  names(delta.all) <- names(se.delta.all)
  #
  delta.all.matrix <- se.delta.all.matrix <- matrix(0, ncol = a, nrow = a)
  #
  k <- 0
  #
  for (i in 1:(a - 1)) {
    for (j in (i + 1):a) {
      k <- k + 1
      delta.all.matrix[i, j] <- delta.all[k]
      delta.all.matrix[j, i] <- -delta.all[k]
      se.delta.all.matrix[i, j] <-
        se.delta.all.matrix[j, i] <- se.delta.all[k]
    }
  }
  #
  #
  colnames(delta.all.matrix) <- rownames(delta.all.matrix) <-
    colnames(se.delta.all.matrix) <- rownames(se.delta.all.matrix) <-
    unique(labels)
  
  comparisons <- c(
    list(studlab = studlab, agent1 = agent1, agent2 = agent2),
    ci(delta, se.delta, level = level))
  #
  #
  all.comparisons <- ci(delta.all.matrix, se.delta.all.matrix,
                        level = level)
  #
  # Test of total heterogeneity / inconsistency:
  #
  Q <- as.vector(t(delta - TE) %*% W %*% (delta - TE))
  
  Q <- max(0, Q)
  
  as.vector(t(delta - TE) %*% W %*% (delta - TE))
  as.vector(t(TE - delta) %*% W %*% (TE - delta))
  
  df.Q <- df1 - qr(Xd)$rank
  #
  pval.Q <- pvalQ(Q, df.Q)
  #
  # Heterogeneity variance
  #
  I <- diag(m)
  E <- matrix(0, nrow = m, ncol = m)
  #
  for (i in 1:m) {
    for (j in 1:m) {
      E[i, j] <- as.numeric(studlab[i] == studlab[j])
    }
  }
  #
  if (df.Q == 0) {
    tau2 <- NA
    tau <- NA
    I2 <- NA
    lower.I2 <- NA
    upper.I2 <- NA
  }
  else {
    tau2 <- max(0, (Q - df.Q) /
                  sum(diag((I - H) %*% (Xd %*% t(Xd) * E / 2) %*% W)))
    #
    tau <- sqrt(tau2)
    #
    ci.I2 <- isquared(Q, df.Q, level)
    #
    I2 <- ci.I2$TE
    lower.I2 <- ci.I2$lower
    upper.I2 <- ci.I2$upper
  }
  
  res <- list(
    comparisons = comparisons,
    all.comparisons = all.comparisons,
    #
    delta = delta,
    se.delta = se.delta,
    #
    TE.nma = ci.beta$TE,
    seTE.nma = ci.beta$seTE,
    lower.nma = ci.beta$lower,
    upper.nma = ci.beta$upper,
    statistic.nma = ci.beta$statistic,
    pval.nma = ci.beta$p,
    #
    Q = Q,
    df.Q = df.Q,
    pval.Q = pval.Q,
    H = sqrt(Q / df.Q),
    #
    tau = tau,
    I2 = I2, lower.I2 = lower.I2, upper.I2 = upper.I2,
    #
    L.matrix = L,
    Lplus.matrix = Lplus,
    H.matrix = H,
    #
    k = length(unique(studlab)),
    m = m,
    a = a
  )
  
  res
}

