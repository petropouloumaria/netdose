#
# Auxiliary functions
#
# Package: netdose
# Authors: Maria Petropoulou <m.petropoulou.a@gmail.com>,
# Guido Schwarzer <guido.schwarzer@uniklinik-freiburg.de>,
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
  if (is.null(seq))
    seq <- agents
  else
    seq <- setseq(seq, agents)
  #
  B.matrix <-
    matrix(0,
           nrow = length(agent1), ncol = length(trts),
           dimnames = list(studlab, trts))
  #
  D.matrix <- D_obs.matrix <-
    matrix(0,
           nrow = length(trts), ncol = length(agents),
           dimnames = list(trts, seq))
  
  # Use regex to separate the agent (all words except the last) and
  # dose (last word)
  data_list <- strsplit(trts, "(?<=\\D)\\s(?=\\d)", perl = TRUE)
  
  # Extract agents and doses
  agent_names <- sapply(data_list, function(x) {
    # If only one element exists, try to split it manually
    if (length(x) == 1) {
      split_x <- unlist(strsplit(x, " "))  # Split by space
      if (length(split_x) > 1) {
        paste(split_x[-length(split_x)], collapse = " ")  # Take all except last
      }
      else {
        split_x  # If it cannot be split, return as is
      }
    }
    else {
      paste(x[-length(x)], collapse = " ")  # Normal case: take all except last
    }
  })
  
  doses <- sapply(data_list, function(x) {
    if (length(x) == 1) {
      split_x <- unlist(strsplit(x, " "))  # Split manually if needed
      dose <- suppressWarnings(as.numeric(split_x[length(split_x)]))  # Convert to number
    }
    else {
      dose <- suppressWarnings(as.numeric(x[length(x)]))  # Normal case
    }
    ifelse(is.na(dose), NA, dose)  # Return NA if conversion fails
  })
  #
  for (i in seq_len(nrow(B.matrix))) {
    B.matrix[i, treat1[i]] <- 1
    B.matrix[i, treat2[i]] <- -1
  }
  #
  for (i in seq_len(nrow(D.matrix))) {
    D_obs.matrix[i, agent_names[i]] <- doses[i]
    D.matrix[i, agent_names[i]] <- do.call(g, list(x = doses[i], p = param))
  }
  #
  X.matrix <- B.matrix %*% D.matrix
  #
  attr(X.matrix, "g") <- deparse(substitute(g))
  attr(X.matrix, "param") <- param
  attr(X.matrix, "seq") <- seq
  #
  class(X.matrix) <- c(class(X.matrix), "Xd1", "Xd")
  #
  list(X.matrix = X.matrix, B.matrix = B.matrix,
       D.matrix = D.matrix, D_obs.matrix = D_obs.matrix)
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
  if (is.null(seq))
    seq <- agents
  else
    seq <- setseq(seq, agents)
  #
  B.matrix <-
    matrix(0,
           nrow = length(agent1), ncol = length(trts),
           dimnames = list(studlab, trts))
  #
  D1 <- D2 <- D_obs.matrix <-
    matrix(0,
           nrow = length(trts), ncol = length(agents),
           dimnames = list(trts, seq))
  
  # Use regex to separate the agent (all words except the last) and dose (last word)
  data_list <- strsplit(trts, "(?<=\\D)\\s(?=\\d)", perl = TRUE)
  
  # Extract agents and doses
  agent_names <- sapply(data_list, function(x) {
    # If only one element exists, try to split it manually
    if (length(x) == 1) {
      split_x <- unlist(strsplit(x, " "))  # Split by space
      if (length(split_x) > 1) {
        paste(split_x[-length(split_x)], collapse = " ")  # Take all except last
      }
      else {
        split_x  # If it cannot be split, return as is
      }
    }
    else {
      paste(x[-length(x)], collapse = " ")  # Normal case: take all except last
    }
  })
  
  doses <- sapply(data_list, function(x) {
    if (length(x) == 1) {
      split_x <- unlist(strsplit(x, " "))  # Split manually if needed
      dose <- suppressWarnings(as.numeric(split_x[length(split_x)])) # Convert to number
    }
    else {
      dose <- suppressWarnings(as.numeric(x[length(x)]))  # Normal case
    }
    ifelse(is.na(dose), NA, dose)  # Return NA if conversion fails
  })
  #
  for (i in seq_len(nrow(B.matrix))) {
    B.matrix[i, treat1[i]] <- 1
    B.matrix[i, treat2[i]] <- -1
  }
  #
  for (i in seq_len(nrow(D1))) {
    D_obs.matrix[i, agent_names[i]] <- doses[i]
    #
    D1[i, agent_names[i]] <- do.call(g1, list(x = doses[i], p = param[1]))
    D2[i, agent_names[i]] <- do.call(g2, list(x = doses[i], p = param[2]))
  }
  #
  X.matrix <- cbind(B.matrix %*% D1, B.matrix %*% D2)
  D.matrix <- cbind(D1, D2)
  #
  attr(X.matrix, "g1") <- deparse(substitute(g1))
  attr(X.matrix, "g2") <- deparse(substitute(g2))
  attr(X.matrix, "param1") <- param[1]
  attr(X.matrix, "param2") <- param[2]
  attr(X.matrix, "seq") <- seq
  #
  class(X.matrix) <- c(class(X.matrix), "Xd2", "Xd")
  #
  list(X.matrix = X.matrix, B.matrix = B.matrix,
       D.matrix = D.matrix, D_obs.matrix = D_obs.matrix)
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
  if (is.null(seq))
    seq <- agents
  else
    seq <- setseq(seq, agents)
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
    mindose.i <- min(dose.i)
    maxdose.i <- max(dose.i)
    #
    knots.i <- quantile(c(mindose.i:maxdose.i), probs = param)
    #
    knots[[agents[i]]] <- knots.i
  }
  #
  B.matrix <- matrix(0,
                     nrow = length(agent1), ncol = length(trts),
                     dimnames = list(studlab, trts)
  )
  #
  D1 <- D2 <- D_obs.matrix <-
    matrix(0,
           nrow = length(trts), ncol = length(agents),
           dimnames = list(trts, seq)
    )
  #
  #
  g1 <- dose2dose
  g2 <- dose2rcs
  #
  
  # Use regex to separate the agent (all words except the last) and dose (last word)
  data_list <- strsplit(trts, "(?<=\\D)\\s(?=\\d)", perl = TRUE)
  
  # Extract agents and doses
  agent_names <- sapply(data_list, function(x) {
    # If only one element exists, try to split it manually
    if (length(x) == 1) {
      split_x <- unlist(strsplit(x, " "))  # Split by space
      if (length(split_x) > 1) {
        paste(split_x[-length(split_x)], collapse = " ")  # Take all except last
      }
      else {
        split_x  # If it cannot be split, return as is
      }
    }
    else {
      paste(x[-length(x)], collapse = " ")  # Normal case: take all except last
    }
  })
  
  doses <- sapply(data_list, function(x) {
    if (length(x) == 1) {
      split_x <- unlist(strsplit(x, " "))  # Split manually if needed
      dose <- suppressWarnings(as.numeric(split_x[length(split_x)]))  # Convert to number
    } else {
      dose <- suppressWarnings(as.numeric(x[length(x)]))  # Normal case
    }
    ifelse(is.na(dose), NA, dose)  # Return NA if conversion fails
  })
  
  #
  for (i in seq_len(nrow(B.matrix))) {
    B.matrix[i, treat1[i]] <- 1
    B.matrix[i, treat2[i]] <- -1
  }
  #
  for (i in seq_len(nrow(D1))) {
    D_obs.matrix[i, agent_names[i]] <- doses[i]
    #
    knots1 <- as.numeric(knots[[agent_names[[i]]]])
    #
    D1[i, agent_names[i]] <- do.call(g1, list(x = doses[i]))
    D2[i, agent_names[i]] <- do.call(g2, list(x = doses[i], p = knots1))
  }
  #
  X.matrix <- cbind(B.matrix %*% D1, B.matrix %*% D2)
  D.matrix <- cbind(D1, D2)
  #
  attr(X.matrix, "g1") <- "dose2dose"
  attr(X.matrix, "g2") <- "dose2rcs"
  attr(X.matrix, "param") <- knots
  attr(X.matrix, "seq") <- seq
  #
  class(X.matrix) <- c(class(X.matrix), "Xd_rcs", "Xd")
  #
  list(X.matrix = X.matrix, B.matrix = B.matrix,
       D.matrix = D.matrix, D_obs.matrix = D_obs.matrix)
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

dose2exp <- function(x, p = NULL) {
  return(1 - exp(-x))
}


dose2rcs <- function(x, p = NULL) {
  if (length(unique(p)) == 1 || length(unique(p)) == 2) {
     return(rep(0, length(x)))
     }
   else {
       return(rcspline.eval(x, knots = p, inclx = FALSE))
     }
}

sel_coef <- function(x, agent, id = 1) {
  x[names(x) %in% agent][id]
}
