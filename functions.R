# ---- Load Libraries ----
library(glmnet)
library(ggplot2)
library(MASS)
library(ggpubr)
library(robustbase)

# ---- General Functions ----

sigmoid <- function(x) {
  return(1/(1 + exp(-x)))
}

calcprob <- function(x, model) {
  # x     - vector of values of explanatory variables, no intercept
  # model - logistic regression model
  # Return: probability of y = 1 given x under the assumption of the logistic
  # model.
  
  x <- c(1, x)
  prob <- sigmoid(coef(model) %*% x)
}

precisionrecall <- function(outliers, alertrates, scores) {
  # outliers   - vector of indexes of observations that are actual outliers
  # alertrates - vector of alert rate values to be considered for precision and
  #              recall evaluation
  # scores     - vector of outlier scores of all observations
  # Return: a list of two vectors - precision and recall calculated depending
  # on the value of alert rate.
  
  n <- length(scores)
  count <- length(outliers)
  precision <- c()
  recall <- c()
  for (a in alertrates) {
    len <- a * n
    prectmp <- sum(outliers %in% order(-scores)[1:len])/len
    precision <- c(precision, prectmp)
    recalltmp <- sum(outliers %in% order(-scores)[1:len])/count
    recall <- c(recall, recalltmp)
  }
  return(list(prec = precision, rec = recall))
}

plotprecisionrecall <- function(precisionrecall, alertrates) {
  # precisionrecall - a list of two vectors (precision and recall) 
  # alertrates      - a vector of alert rates for which the provided precision
  #                   and recall were calculated
  # Return: a figure containing plots of PAR, RAR, and PR curves.
  
  plotdata <- data.frame(alertrates = alertrates,
                         precision = precisionrecall$prec,
                         recall = precisionrecall$rec)
  plot1 <- ggplot(plotdata) + 
    geom_line(aes(x = alertrates, y = precision), size = 0.7) +
    labs(x = "Alert Rate", y = "Precision") +
    theme_bw() + 
    ylim(0, 1)
  plot2 <- ggplot(plotdata) + 
    geom_line(aes(x = alertrates, y = recall), size = 0.7) +
    labs(x = "Alert Rate", y = "Recall") +
    theme_bw() + 
    ylim(0, 1)
  plot3 <- ggplot(plotdata) + 
    geom_line(aes(x = recall, y = precision), size = 0.7) +
    labs(x = "Recall", y = "Precision") +
    theme_bw() + 
    ylim(0, 1) + 
    xlim(0, 1)
  finalplot <- ggarrange(plot1, plot2, plot3, ncol = 3, nrow = 1)
  return(finalplot)
}

# ---- Binary Relevance ----

binaryrel <- function(X, Y, type = "Penalized", alpha = 0, outliers = NA) {
  # X        - a data frame with p-1 columns and n rows, the constant column
  #            should not be included
  # Y        - a data frame with k columns and n rows
  # type     - specifies the type of underlying logistic regression in the
  #            BR model, the options are "Penalized", "Ordinary", "Pearson",
  #            "Oracle", "BY" (Bianco-Yohai), and "MH" (Mallows-Huber)
  # alpha    - penalty parameter, only relevant when type == "Penalized"
  # outliers - vector of indexes of actual outliers, only relevant when 
  #            type == "Oracle"
  # Return: a list of k models for each of the response variables.
  
  k <- ncol(Y)
  names <- c()
  models <- list()
  
  for (i in 1:k) {
    if (type == "Penalized") {
      if (ncol(X) == 1) X$dummy <- 1
      model <- cv.glmnet(as.matrix(X), 
                         as.factor(Y[, i]), 
                         family = "binomial", 
                         alpha = alpha)
    } else if (type == "Ordinary") {
      model <- glm(as.factor(Y[, i]) ~., family = "binomial", data = X)
    } else if (type == "Pearson") {
      model <- glm(as.factor(Y[, i]) ~., family = "binomial", data = X)
      pred <- predict(model, X, type = "response")
      residpearson <- (Y[, i] - pred) / sqrt((pred * (1 - pred)))
      w <- ifelse(abs(residpearson) > 3, 0, 1)
      model <- 
        glm(as.factor(Y[, i]) ~., family = "binomial", data = X, weights = w)
    } else if (type == "Oracle") {
      w <- rep(1, n)
      w[outliers] <- 0
      model <- 
        glm(as.factor(Y[, i]) ~., family = "binomial", data = X, weights = w)
    } else if (type == "BY") {
      model <- glmrob(Y[, i] ~., data = X, family = "binomial", method = "BY") 
    } else if (type == "MH") {
      model <- glmrob(as.factor(Y[, i]) ~., 
                      data = X, 
                      family = "binomial", 
                      method = "Mqle", 
                      control = glmrobMqle.control(tcc = 1.2))
    }
    
    models[[i]] <- model
    names <- c(names, paste("m", i, sep = ""))
  }
  
  names(models) <- names
  
  return(models) 
}

estprobbr <- function(models, X, Y) {
  # models - a list of k models representing the BR model
  # X      - a data frame with p-1 columns and n rows, the constant column
  #            should not be included
  # Y      - a data frame with k columns and n rows
  # Return: a vector of length n containing the values of the conditional
  # probability P(y|x) for the corresponding observation estimated based on the
  # BR model.
  
  k <- length(models)
  probs <- rep(1, nrow(X))
  
  if (class(models$m1)[1] == "cv.glmnet") {
    for (i in 1:k) {
      probstmp <- 
        predict(models[[i]], as.matrix(X), type = "response", s = "lambda.min")
      probstmp <- probstmp^Y[, i] * (1 - probstmp)^(1 - Y[, i])
      probs <- probs * probstmp
    } 
  }
  else {
    for (i in 1:k) {
      probstmp <- apply(X, 1, calcprob, model = models[[i]])
      probstmp <- probstmp^Y[, i] * (1 - probstmp)^(1 - Y[, i])
      probs <- probs * probstmp
    }
  }
  
  return(probs)
}

# ---- Classifier Chains ----

classchain <- function(X, Y, type = "Penalized", alpha = 0, outliers = NA) {
  # X        - a data frame with p-1 columns and n rows, the constant column
  #            should not be included
  # Y        - a data frame with k columns and n rows
  # type     - specifies the type of underlying logistic regression in the
  #            CC model, the options are "Penalized", "Ordinary", "Pearson",
  #            "Oracle", "BY" (Bianco-Yohai), and "MH" (Mallows-Huber)
  # alpha    - penalty parameter, only relevant when type == "Penalized"
  # outliers - vector of indexes of actual outliers, only relevant when 
  #            type == "Oracle"
  # Return: a list of k models for each of the response variables.
  
  k <- ncol(Y)
  names <- c()
  models <- list()
  
  for (i in 1:k) {
    if (i > 1) {
      X[, paste("Y", (i - 1), sep = "")] <- as.numeric(Y[, (i - 1)])
    } 
    if (type == "Penalized") {
      if (ncol(X) == 1) X$dummy <- 1
      model <- cv.glmnet(as.matrix(X),
                         as.factor(Y[, i]),
                         family = "binomial",
                         alpha = alpha)
    } else if (type == "Ordinary") {
      model <- glm(as.factor(Y[, i]) ~., family = "binomial", data = X)
    } else if (type == "Pearson") {
      model <- glm(as.factor(Y[, i]) ~., family = "binomial", data = X)
      pred <- predict(model, X, type = "response")
      residpearson <- (Y[, i] - pred) / sqrt((pred * (1 - pred)))
      w <- ifelse(abs(residpearson) > 3, 0, 1)
      model <- 
        glm(as.factor(Y[, i]) ~., family = "binomial", data = X, weights = w)
    } else if (type == "Oracle") {
      w <- rep(1, n)
      w[outliers] <- 0
      model <- 
        glm(as.factor(Y[, i]) ~., family = "binomial", data = X, weights = w)
    } else if (type == "BY") {
      model <- glmrob(Y[, i] ~., data = X, family = "binomial", method = "BY") 
    } else if (type == "MH") {
      model <- glmrob(as.factor(Y[, i]) ~., 
                      data = X, 
                      family = "binomial", 
                      method = "Mqle", 
                      control = glmrobMqle.control(tcc = 1.2))
    }
    
    models[[i]] <- model
    names <- c(names, paste("m", i, sep = ""))
  }
  
  names(models) <- names
  
  return(models) 
}

estprobcc <- function(models, X, Y) {
  # models - a list of k models representing the CC model
  # X      - a data frame with p-1 columns and n rows, the constant column
  #            should not be included
  # Y      - a data frame with k columns and n rows
  # Return: a vector of length n containing the values of the conditional
  # probability P(y|x) for the corresponding observation estimated based on the
  # CC model.
  
  k <- length(models)
  
  if (class(models$m1)[1] == "cv.glmnet") {
    probs <- 
      predict(models[[1]], as.matrix(X), type = "response", s = "lambda.min")
    probs <- probs^Y[, 1] * (1 - probs)^(1 - Y[, 1])
    
    for (i in 2:k) {
      X[, paste("Y", (i - 1), sep = "")] <- as.numeric(Y[, (i - 1)])
      probstmp <- 
        predict(models[[i]], as.matrix(X), type = "response", s = "lambda.min")
      probstmp <- probstmp^Y[, i] * (1 - probstmp)^(1 - Y[, i])
      probs <- probs * probstmp
    }
  }
  else {
    probs <- apply(X, 1, calcprob, model = models[[1]])
    probs <- probs^Y[, 1] * (1 - probs)^(1 - Y[, 1])
    
    for (i in 2:k) {
      X[, paste("Y", (i - 1), sep = "")] <- as.numeric(Y[, (i - 1)])
      probstmp <- apply(X, 1, calcprob, model = models[[i]])
      probstmp <- probstmp^Y[, i] * (1 - probstmp)^(1 - Y[, i])
      probs <- probs * probstmp
    }
  }
  
  return(probs)
}

# ---- Ising Model ----

ising <- function(X, Y, type = "Penalized", alpha = 0, outliers = NA) {
  # X        - a data frame with p-1 columns and n rows, the constant column
  #            should not be included
  # Y        - a data frame with k columns and n rows
  # type     - specifies the type of underlying logistic regression in the
  #            Ising model, the options are "Penalized", "Ordinary", "Pearson",
  #            "Oracle", "BY" (Bianco-Yohai), and "MH" (Mallows-Huber)
  # alpha    - penalty parameter, only relevant when type == "Penalized"
  # outliers - vector of indexes of actual outliers, only relevant when 
  #            type == "Oracle"
  # Return: a list of k models for each of the response variables
  
  k <- ncol(Y)
  names <- c()
  models <- list()
  
  for (i in 1:k) {
    Xtmp <- data.frame(X, Y[, -i])
    
    if (type == "Penalized") {
      model <- cv.glmnet(as.matrix(Xtmp),
                         as.factor(Y[, i]),
                         family = "binomial",
                         alpha = alpha)
    } else if (type == "Ordinary") {
      model <- glm(as.factor(Y[, i]) ~., family = "binomial", data = Xtmp)
    } else if (type == "Pearson") {
      model <- glm(as.factor(Y[, i]) ~., family = "binomial", data = Xtmp)
      pred <- predict(model, Xtmp, type = "response")
      residpearson <- (Y[, i] - pred) / sqrt((pred * (1 - pred)))
      w <- ifelse(abs(residpearson) > 3, 0, 1)
      model <- 
        glm(as.factor(Y[, i]) ~., family = "binomial", data = Xtmp, weights = w)
    } else if (type == "Oracle") {
      w <- rep(1, n)
      w[outliers] <- 0
      model <- 
        glm(as.factor(Y[, i]) ~., family = "binomial", data = Xtmp, weights = w)
    } else if (type == "BY") {
      model <- glmrob(Y[, i] ~., data = Xtmp, family = "binomial", method = "BY") 
    } else if (type == "MH") {
      model <- glmrob(as.factor(Y[, i]) ~., 
                      data = Xtmp, 
                      family = "binomial", 
                      method = "Mqle", 
                      control = glmrobMqle.control(tcc = 1.2))
    }
    
    models[[i]] <- model
    names <- c(names, paste("m", i, sep = ""))
  }
  
  names(models) <- names
  
  return(models) 
}

coefising <- function(models) {
  # models - a list of k models representing the Ising model
  # Return: a list containing a matrix and a vector. Each row of the matrix is
  # a beta vector from the first sum specified in 3.3. Each element of the
  # vector is a beta value from the second sum. Altogether the list includes
  # all coefficients of the Ising model. The coefficients are estimated using
  # Ising's model connection to the logistic regression. If the coefficients 
  # are estimated by more than one logistic regression model, they are averaged.
  
  p <- sum(grepl("X", rownames(as.matrix(coef(models$m1))))) + 1
  k <- length(models)
  vec <- 1:k
  
  b.ll <- matrix(nrow = k, ncol = p)
  b.il <- matrix(nrow = k, ncol = k)
  
  if (class(models$m1)[1] == "cv.glmnet") {
    for (l in vec) {
      coeftmp <- t(as.matrix(coef(models[[l]], s = "lambda.min")))
      b.ll[l, ] <- coeftmp[1:p]
      
      coeftmp <- append(coeftmp[-(1:p)], NA, after = l - 1)
      for (i in vec[vec != l]) {
        b.il[i, l] <- coeftmp[i]
      }
    }
  }
  else {
    for (l in vec) {
      coeftmp <- t(as.matrix(coef(models[[l]])))
      b.ll[l, ] <- coeftmp[1:p]
      
      coeftmp <- append(coeftmp[-(1:p)], NA, after = l - 1)
      for (i in vec[vec != l]) {
        b.il[i, l] <- coeftmp[i]
      }
    }
  }
  
  tmp <- matrix(nrow = 2, ncol = k * (k - 1)/2)
  tmp[1, ] <- b.il[lower.tri(b.il)]
  b.il <- t(b.il)
  tmp[2, ] <- b.il[lower.tri(b.il)]
  b.il <- apply(tmp, 2, function(x) mean(x))
  
  return(list(b.ll = b.ll, b.il = b.il))
}

Nising <- function(x, coef) {
  # x    - a vector of length p
  # coef - a list returned by coefising function, coefficients of the Ising
  #        model
  # Return: a number which is the normalizing constant present in the Ising 
  # model.
  
  k <- nrow(coef$b.ll)
  
  x <- t(as.matrix(x))
  X <- t(x[rep(1, k), ])
  
  Yall <- expand.grid(replicate(k, 0:1, simplify = F))
  sum <- 0
  
  for (l in 1:2^k) {
    Y <- Yall[l, ]
    insum <- (coef$b.ll %*% X)[, 1] %*% t(Y)
    
    b.il <- matrix(nrow = k, ncol = k)
    b.il[lower.tri(b.il)] <- coef$b.il
    b.il <- t(b.il)
    
    for (i in 1:(k - 1)) {
      for (j in (i + 1):k) {
        insum <- insum + b.il[i, j] * Y[i] * Y[j]
      }
    }
    sum <- sum + exp(insum)
  }
  
  return(sum)
}

estprobising <- function(coef, X, Y) {
  # coef - a list returned by coefising function, coefficients of the Ising
  #        model
  # X    - a data frame with p-1 columns and n rows, the constant column should
  #        not be included
  # Y    - a data frame with k columns and n rows
  # Return: a vector of length n containing the values of the conditional
  # probability P(y|x) for the corresponding observation estimated based on the
  # Ising model.
  
  n <- nrow(X)
  k <- ncol(Y)
  probs <- matrix(nrow = n, ncol = 1)
  
  b.il <- matrix(nrow = k, ncol = k)
  b.il[lower.tri(b.il)] <- coef$b.il
  b.il <- t(b.il)
  
  for (l in 1:n) {
    Xl <- as.matrix(c(1, as.matrix(X[l, ])))
    Yl <- Y[l, ]
    sum <- (coef$b.ll %*% Xl)[, 1] %*% t(Yl)
    
    for (i in 1:(k - 1)) {
      for (j in (i + 1):k) {
        sum <- sum + b.il[i, j] * Yl[i] * Yl[j]
      }
    }
    
    probs[l, 1] <- as.numeric(exp(sum) / Nising(Xl, coef))
  }
  
  return(probs)
}

# ---- Synthetic Data 1 ----

  generatesd1 <- function(n, p, k, coefX, coefY, X = NA) {
  # n     - the number of observations to be generated
  # p     - the number of explanatory variables, including the intercept
  # k     - the number of response variables
  # coefX - a vector of length p specifying the coefficients used in generation
  # coefY - a vector of length k-1 specifying the coefficients used in 
  #         generation
  # X     - a data frame; if not given, every column is generated from the
  #         standard normal distribution
  # Return: a list of the generated explanatory variables X and response 
  # variables Y, X does not include the constant column. The data is generated
  # according to the SD1 method.
  
  p <- p - 1
  
  if (is.na(X)) {
    X <- data.frame(matrix(nrow = n, ncol = p))
    for (i in 1:p) {
      X[, i] <- rnorm(n)
    } 
  }
  
  Y <- data.frame(matrix(nrow = n, ncol = k))
  for (i in 1:n) {
    Y[i, 1] <- rbinom(1, 1, sigmoid(sum(coefX * c(1, as.numeric(X[i, ])))))
    for (j in 2:k) {
      Y[i, j] <- rbinom(1, 1, sigmoid(sum(coefX * c(1, as.numeric(X[i, ]))) 
                                      + sum(coefY[1:(j-1)] * Y[i, 1:(j-1)])))
    }
  }
  names(Y) <- paste("Y", 1:k, sep = "")
  
  return(list(X = X, Y = Y))
}

probsd1 <- function(x, coefX, coefY) {
  # x     - a vector of explanatory variables, no intercept
  # coefX - a vector of coefficients used in data generation
  # coefY - a vector of coefficients used in data generation
  # Return: a vector of probabilities P(y|x) for all possible y.
  
  k <- length(coefY) + 1
  yall <- expand.grid(replicate(k, 0:1, simplify = F))
  prob <- c()
  
  probsinglesd1 <- function(x, y, coefX, coefY) {
    k <- length(y)
    prob <- sigmoid(sum(coefX * c(1, x)))^y[1] * 
      (1 - sigmoid(sum(coefX * c(1, x))))^(1 - y[1])
    for (i in 2:k) {
      lincomb <- sum(coefX * c(1, x)) + sum(coefY[1:(i-1)] * y[1:(i-1)])
      prob <- prob * sigmoid(lincomb)^y[i] * (1 - sigmoid(lincomb))^(1 - y[i])
    }
    
    return(prob)
  }
  
  for (i in 1:2^k) {
    prob <- c(prob, as.numeric(probsinglesd1(x, yall[i, ], coefX, coefY)))
  }
  
  return(prob)
}

addoutlierssd1 <- function(X, Y, percent, coefX, coefY) {
  # X       - a data frame
  # Y       - a data frame
  # percent - a value, percent of outliers to be added to the data set
  # coefX   - a vector of coefficients used in data generation
  # coefY   - a vector of coefficients used in data generation
  # Return: a list of a vector and a data frame, the vector contains the 
  # indexes of outliers and the data frame contains new values of the response
  # variables.
  
  n <- nrow(Y)
  k <- ncol(Y)
  count <- percent * n
  outliers <- sample.int(n, count)
  Ynew <- Y
  Yall <- expand.grid(replicate(k, 0:1, simplify = F))
  
  for (o in outliers) {
    Ynew[o, ] <- 
      Yall[which.min(probsd1(as.numeric(X[o, ]), coefX, coefY)), ]
  }
  
  return(list(outliers = outliers, Ynew = Ynew))
}

# ---- Synthetic Data 2 ----

generatesd2 <- function(n, mu, sigma, probvec) {
  # n       - the number of observations to be generated
  # mu      - a list of four two-dimensional vectors
  # sigma   - a list of four 2x2 matrices
  # probvec - a vector of length 2
  # Return: a list of the generated explanatory variables X and response 
  # variables Y, X does not include the constant column. The data is generated
  # according to the SD2 method. Each column of Y is generated from Bernoulli
  # distribution with probability of success sourced from probvec.
  
  Y <- data.frame(Y1 = rbinom(n, 1, probvec[1]), Y2 = rbinom(n, 1, probvec[2]))
  X <- data.frame(matrix(nrow = n, ncol = 2))
  
  for (i in 1:n) {
    if (sum(Y[i,] == c(0, 0)) == 2) X[i,] = mvrnorm(1, mu$mu00, sigma$s00)
    else if (sum(Y[i,] == c(0, 1)) == 2) X[i,] = mvrnorm(1, mu$mu01, sigma$s01)
    else if (sum(Y[i,] == c(1, 0)) == 2) X[i,] = mvrnorm(1, mu$mu10, sigma$s10)
    else if (sum(Y[i,] == c(1, 1)) == 2) X[i,] = mvrnorm(1, mu$mu11, sigma$s11)
  }
  
  return(list(X = X, Y = Y))
}

addnoise <- function(X, m) {
  # X - a data frame, clean data generated by generatesd2 function
  # m - the number of noise variables to be added to the data
  # Return: new data frame X with added noise variables.
  
  n <- nrow(X)
  Xnoise <- mvrnorm(n, rep(0, m), 0.25 * sqrt(m + 2) * diag(m))
  Xnew <- data.frame(X, Xnoise)
  names(Xnew) <- paste("X", 1:(m + 2), sep = "")
  return(Xnew)
}

probsd2 <- function(x, mu, sigma, probvec) {
  # x       - a vector of explanatory variables, no intercept
  # mu      - a list of four two-dimensional vectors used in data generation
  # sigma   - a list of four 2x2 matrices used in data generation
  # probvec - a vector of length 2 used in data generation
  # Return: a vector of probabilities P(y|x) for all possible y.
  
  yall <- expand.grid(replicate(2, 0:1, simplify = F))
  prob <- c()
  
  probsinglesd2 <- function(x, y, mu, sigma, probvec) {
    prob00 <- mvdnorm(x, mu$mu00, sigma$s00) * prod(1 - probvec)
    prob01 <- mvdnorm(x, mu$mu01, sigma$s01) * (1 - probvec[1]) * probvec[2]
    prob10 <- mvdnorm(x, mu$mu10, sigma$s10) * probvec[1] * (1 - probvec[2])
    prob11 <- mvdnorm(x, mu$mu11, sigma$s11) * prod(probvec)
    const <- prob00 + prob01 + prob10 + prob11
    if (sum(y == c(0, 0)) == 2) {
      prob <- prob00/const
    } else if (sum(y == c(0, 1)) == 2) {
      prob <- prob01/const
    } else if (sum(y == c(1, 0)) == 2) {
      prob <- prob10/const
    } else if (sum(y == c(1, 1)) == 2) {
      prob <- prob11/const
    }
    
    return(prob)
  }
  
  for (i in 1:2^2) {
    prob <- 
      c(prob, as.numeric(probsinglesd2(x, yall[i, ], mu, sigma, probvec)))
  }
  
  return(prob)
}

addoutlierssd2 <- function(X, Y, percent, mu, sigma, probvec) {
  # X       - a data frame
  # Y       - a data frame
  # percent - a value, percent of outliers to be added to the data set
  # mu      - a list of four two-dimensional vectors used in data generation
  # sigma   - a list of four 2x2 matrices used in data generation
  # probvec - a vector of length 2 used in data generation
  # Return: a list of a vector and a data frame, the vector contains the 
  # indexes of outliers and the data frame contains new values of the response
  # variables.
  
  n <- nrow(Y)
  count <- percent * n
  outliers <- sample.int(n, count)
  Ynew <- Y
  Yall <- expand.grid(replicate(2, 0:1, simplify = F))
  
  for (o in outliers) {
    Ynew[o, ] <- Yall[which.min(probsd2(X[o, ], mu, sigma, probvec)), ]
  }
  
  return(list(outliers = outliers, Ynew = Ynew))
}
