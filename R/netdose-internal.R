#
# Auxiliary functions
#
# Package: netdose
# Authors: Guido Schwarzer <guido.schwarzer@uniklinik-freiburg.de>,
# Maria Petropoulou <maria.petropoulou@uniklinik-freiburg.de>
# License: GPL (>= 2)
#

allNA <- function(x) {
  all(is.na(x))
}

catch <- function(argname, matchcall, data, encl) {
  eval(matchcall[[match(argname, names(matchcall))]], data, enclos = encl)
}

int2num <- function(x) {
  #
  # Convert integer to numeric
  #
  if (is.integer(x)) {
    res <- as.numeric(x)
  } else {
    res <- x
  }
  #
  res
}

npn <- function(x) {
  #
  # Check for non-positive values in vector
  #
  selNA <- is.na(x)
  res <- selNA
  if (sum(!selNA) > 0) {
    x[!selNA] <- x[!selNA] <= 0
  }
  #
  res
}

replaceNULL <- function(x, replace = NA) {
  if (is.null(x)) {
    return(replace)
  }
  x
}

replaceNA <- function(x, replace = NA) {
  if (is.null(x)) {
    return(x)
  } else {
    x[is.na(x)] <- replace
  }
  x
}

replaceVal <- function(x, old, new) {
  if (is.null(x)) {
    return(x)
  } else {
    x[x == old] <- new
  }
  x
}

extrVar <- function(x, name) {
  x[[name]]
}

calcPercent <- function(x) {
  100 * x / sum(x, na.rm = TRUE)
}

list2vec <- function(x) {
  if (is.list(x)) {
    return(do.call("c", x))
  } else {
    return(x)
  }
}

setsv <- function(x) {
  if (is.null(x)) {
    res <- "desirable"
  } else {
    res <- setchar(x, c("good", "bad"), stop.at.error = FALSE)
    #
    if (!is.null(res)) {
      res <- switch(res,
        good = "desirable",
        bad = "undesirable"
      )
    } else {
      res <- x
    }
  }
  #
  setchar(res, c("desirable", "undesirable"))
}

createXd1 <- function(agent1, dose1, agent2, dose2, studlab, data = NULL,
                      g = dose2dose, param = NULL, seq = NULL) {
  sfsp <- sys.frame(sys.parent())
  mc <- match.call()
  #
  agent1 <- catch("agent1", mc, data, sfsp)
  dose1 <- catch("dose1", mc, data, sfsp)
  agent2 <- catch("agent2", mc, data, sfsp)
  dose2 <- catch("dose2", mc, data, sfsp)
  studlab <- catch("studlab", mc, data, sfsp)
  #
  treat1 <- paste(agent1, dose1)
  treat2 <- paste(agent2, dose2)
  #
  agents <- sort(unique(c(agent1, agent2)))
  trts <- sort(unique(c(treat1, treat2)))
  #
  if (is.null(seq)) {
    seq <- agents
  } else {
    seq <- setseq(seq, agents)
  }
  #
  B.matrix <-
    matrix(0,
      nrow = length(agent1), ncol = length(trts),
      dimnames = list(studlab, trts)
    )
  #
  D <-
    matrix(0,
      nrow = length(trts), ncol = length(agents),
      dimnames = list(trts, seq)
    )

  # Split the vector into two parts
  split_vector <- strsplit(trts, " ")

  # Extract the values and characters
  values <- sapply(split_vector, function(x) as.numeric(x[2]))
  characters <- sapply(split_vector, function(x) x[1])
  #
  for (i in seq_len(nrow(B.matrix))) {
    B.matrix[i, treat1[i]] <- 1
    B.matrix[i, treat2[i]] <- -1
  }
  #
  for (i in seq_len(nrow(D))) {
    D[i, characters[i]] <- do.call(g, list(x = values[i], p = param))
  }
  #
  Xd <- B.matrix %*% D
  #
  attr(Xd, "g") <- deparse(substitute(g))
  attr(Xd, "param") <- param
  attr(Xd, "seq") <- seq
  #
  class(Xd) <- c(class(Xd), "Xd1", "Xd")
  #
  list(Xd = Xd, D = D)
}


createXd2 <- function(agent1, dose1, agent2, dose2, studlab, data = NULL,
                      g1 = dose2dose, g2 = dose2poly,
                      param = NULL, seq = NULL) {
  sfsp <- sys.frame(sys.parent())
  mc <- match.call()
  #
  agent1 <- catch("agent1", mc, data, sfsp)
  dose1 <- catch("dose1", mc, data, sfsp)
  agent2 <- catch("agent2", mc, data, sfsp)
  dose2 <- catch("dose2", mc, data, sfsp)
  studlab <- catch("studlab", mc, data, sfsp)
  #
  treat1 <- paste(agent1, dose1)
  treat2 <- paste(agent2, dose2)
  #
  agents <- sort(unique(c(agent1, agent2)))
  trts <- sort(unique(c(treat1, treat2)))
  #
  if (is.null(seq)) {
    seq <- agents
  } else {
    seq <- setseq(seq, agents)
  }
  #
  B.matrix <- matrix(0,
    nrow = length(agent1), ncol = length(trts),
    dimnames = list(studlab, trts)
  )
  #
  D1 <- D2 <-
    matrix(0,
      nrow = length(trts), ncol = length(agents),
      dimnames = list(trts, seq)
    )
  #
  # Split the vector into two parts
  split_vector <- strsplit(trts, " ")

  # Extract the values and characters
  values <- sapply(split_vector, function(x) as.numeric(x[2]))
  characters <- sapply(split_vector, function(x) x[1])
  #
  for (i in seq_len(nrow(B.matrix))) {
    B.matrix[i, treat1[i]] <- 1
    B.matrix[i, treat2[i]] <- -1
  }
  #
  for (i in seq_len(nrow(D1))) {
    D1[i, characters[i]] <- do.call(g1, list(x = values[i], p = param[1]))
    D2[i, characters[i]] <- do.call(g2, list(x = values[i], p = param[2]))
  }
  #
  Xd <- cbind(B.matrix %*% D1, B.matrix %*% D2)
  D <- cbind(D1, D2)
  #
  #
  attr(Xd, "g1") <- deparse(substitute(g1))
  attr(Xd, "g2") <- deparse(substitute(g2))
  attr(Xd, "param1") <- param[1]
  attr(Xd, "param2") <- param[2]
  attr(Xd, "seq") <- seq
  #
  class(Xd) <- c(class(Xd), "Xd2", "Xd")
  #
  list(Xd = Xd, D = D)
}


createXd_rcs <- function(agent1, dose1, agent2, dose2, studlab, data = NULL,
                         seq = NULL, param = NULL) {
  sfsp <- sys.frame(sys.parent())
  mc <- match.call()
  #
  agent1 <- catch("agent1", mc, data, sfsp)
  dose1 <- catch("dose1", mc, data, sfsp)
  agent2 <- catch("agent2", mc, data, sfsp)
  dose2 <- catch("dose2", mc, data, sfsp)
  studlab <- catch("studlab", mc, data, sfsp)
  #
  treat1 <- paste(agent1, dose1)
  treat2 <- paste(agent2, dose2)
  #
  agents <- sort(unique(c(agent1, agent2)))
  trts <- sort(unique(c(treat1, treat2)))
  #
  if (is.null(seq)) {
    seq <- agents
  } else {
    seq <- setseq(seq, agents)
  }
  #
  # Calculate knots for agents
  #
  knots <- vector("list", length(agents))
  names(knots) <- agents
  #
  for (i in seq_along(agents)) {
    # Harrell's suggestion (fixed sample quantiles)
    if (is.null(param)) {
      param <- c(0.10, 0.50, 0.90)
    }
    dose.i <- c(dose1[agent1 == agents[i]], dose2[agent2 == agents[i]])
    #
    knots.i <-
      quantile(c(min(dose.i, na.rm = TRUE), max(dose.i, na.rm = TRUE)),
        probs = param
      )
    #
    knots[[agents[i]]] <- knots.i
  }
  #
  B.matrix <- matrix(0,
    nrow = length(agent1), ncol = length(trts),
    dimnames = list(studlab, trts)
  )
  #
  D1 <- D2 <-
    matrix(0,
      nrow = length(trts), ncol = length(agents),
      dimnames = list(trts, seq)
    )
  #
  #
  g1 <- dose2dose
  g2 <- dose2rcs
  #
  # Split the vector into two parts
  split_vector <- strsplit(trts, " ")

  # Extract the values and characters
  values <- sapply(split_vector, function(x) as.numeric(x[2]))
  characters <- sapply(split_vector, function(x) x[1])
  #
  for (i in seq_len(nrow(B.matrix))) {
    B.matrix[i, treat1[i]] <- 1
    B.matrix[i, treat2[i]] <- -1
  }
  #

  for (i in seq_len(nrow(D1))) {
    D1[i, characters[i]] <- do.call(g1, list(x = values[i]))
    D2[i, characters[i]] <- do.call(g2, list(x = values[i], p = knots[[characters[i]]]))
  }
  #
  Xd <- cbind(B.matrix %*% D1, B.matrix %*% D2)
  D <- cbind(D1, D2)
  #
  attr(Xd, "g1") <- "dose2dose"
  attr(Xd, "g2") <- "dose2rcs"
  attr(Xd, "param") <- knots
  attr(Xd, "seq") <- seq
  #
  class(Xd) <- c(class(Xd), "Xd_rcs", "Xd")
  #
  list(Xd = Xd, D = D)
}


dose2dose <- function(x, p = NULL) {
  return(x)
}


dose2poly <- function(x, p = NULL) {
  if (is.null(p)) {
    p <- 2
  }
  x^p
}

# Define the fractional polynomial transformation function
dose2fp <- function(x, p = -0.5, epsilon = 0.001) {
  if (p == 0) {
    # For p = 0, apply logarithmic transformation with epsilon, and return exact 0 when x is exactly 0
    return(ifelse(x == 0, 0, log(x + epsilon)))
  } else {
    # For p != 0, return x^p for non-zero values, and exact 0 for x = 0
    return(ifelse(x == 0, 0, x^p))
  }
}

#  TEST
# dose2fp <- function(x, p = NULL) {
#   if (is.null(p)) {
#    p <- -0.5
#   }
#   ifelse(x == 0, 1, x^p)
# }

dose2exp <- function(x, p = NULL) {
  return(1 - exp(-x))
}

dose2rcs <- function(x, p = NULL) {
  if (length(unique(p)) == 1) {
    return(x)
  } else {
    return(rcspline.eval(x, knots = p, inclx = FALSE))
  }
}

sel_coef <- function(x, agent, id = 1) {
  x[names(x) %in% agent][id]
}
