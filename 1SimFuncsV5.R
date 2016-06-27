library(metafor)
library(rpart)
library(partykit)

NumGen <- function(n, k){
  # a function to generate the within-study sample sizes based on the average sample size
  #
  # Arguments:
  #   n: the average sample size, needs to be a scalar
  #   k: the number of sudies, needs to be a scalar
  #
  # Returns:
  # a numeric vector of length k
  n <- n
  k <- k
  samp <- rnorm(k, n, n/3)
  samp <- round(samp, 0)  # make the generated number integer
  samp[samp < 10] <- 10  # avoid negative, 0 or too small sample size that is not realistic
  return(samp)
}

ModsGen.dich <- function(m, K, p = 0.5){
  # m : the number of moderators
  # K : the number of studies
  # p : the probability of 1
  ms <- replicate(m,rbinom(K, 1, p))  # generate potential moderators from a Bernoulli distribution
  ms <- as.data.frame(ms) 
  colnames(ms) <- paste("m", 1:m, sep="")
  ms
}

ModsGen.multi <- function(m, lvls, K){
  # m : the number of moderators
  # lvls : the numberof levels
  # K : the number of studies
  # p : assuming all levels having the same probability
  if (length(lvls) !=  m) stop("The length of levels does not equal to the number of moderators")
  ms <- sapply(1:m, function(x) sample(LETTERS[1:lvls[x]], size = K, replace = T ))
  ms <- as.data.frame(ms)
  colnames(ms) <- paste("m", 1:m, sep="")
  ms
}

ModsGen.ordinal <- function(m, levels, K){
  # m : the number of moderators
  # lvls : the numberof levels, a vector of length m
  # K : the number of studies
  # p : assuming all levels having the same probability
  if (length(lvls) !=  m) stop("The length of levels does not equal to the number of moderators")
  ms <- sapply(1:m, function(x) sample(1:lvls[x], size = K, replace = T ))
  ms <- as.data.frame(ms)
  colnames(ms) <- paste("m", 1:m, sep="")
  ms
}

SimData <- function(beta, mods, formula,
                    K, n, tau){
  # a function to simulate meta-analysis data set with moderators  
  # corrected for non-central t
  #
  # Arguments :
  #      beta : a numercia vector which contains the overall effect size
  #             and other coefficients. The length should equal to the 
  #             column number of model matrix
  #      mods : the dataframe of moderators.
  #   formula : a formula object to specify the relationship between the 
  #             effect size and moderators  
  #         K : a scalar which is the number of studies
  #        n  : a numerical vector of length K, the  within-study sample
  #             size
  #       tau : a scalar, the effect sizes variance
  #     
  #
  # Returns   :
  # A data frame which is simulated by the given parameters
  mods <- mods
  beta <- beta
  formula <- as.formula(formula)
  K <- K
  n <- n
  tau <- tau
  # step 1
  X <- model.matrix(formula, data = mods)
  # step 2
  delta <- X%*%beta  # compute the average effect size corresponding to model and the coefficients
  # step 3
  deltak <- rnorm(rep(1,K), delta, tau)  # sample the true effect size for single study
  cmk <- 1 - 3/(8*n-9)  #  approximation for
  #gamma(n - 1) / (sqrt(n - 1) * gamma((2 * n - 3) / 2))
  efk <- vark <- numeric(K)  # generate the effect size from a non-central t-distribution
  for (k in 1:K) {
    samp <- rt(1, df=2*n[k]-2, ncp=deltak[k]*sqrt(n[k]/2))
    efk[k] <- samp*cmk[k]/sqrt(n[k]/2)
    vark[k] <- cmk[k]^2*(2*n[k]-2)*
      (1+n[k]*deltak[k]^2/2)/((2*n[k]-4)*n[k]/2)-deltak[k]^2  # note the vark is the sampling variance but not the sampling error
  }
  dat   <- data.frame(trail=1:K, efk=efk, vark=vark)
  return(dat)
}

treepruner <- function(tree, c){
  # function that prunes a CART with c-SE prunign rule
  #
  # Argument:
  #  tree: a classification tree or a regression tree fit by rpart
  #     c: the pruning parameter, needs to be a scalar
  # 
  # Returns:
  # a pruned tree
  tree <- tree
  c <- c
  mindex <- which.min(tree$cptable[,4])  # find the row of the minimum x-error
  cp.minse <- tree$cptable[mindex,4] + c*tree$cptable[mindex,5]  # the minimum x-error + c*SE
  cp.row <- min(which(tree$cptable[,4]<= cp.minse))  # find the smallest tree within the minimum x-error + c*SE
  cp.take <- tree$cptable[cp.row, 1]  # get the cp value for the smallest tree
  prune(tree, cp=cp.take)  # prune the tree
}

getp <- function(prunetree, testdat){
  # compute the significance of moderator effects. the second step of meta-CART. 
  #
  # Arguments:
  #  prunetree: the final tree in the first step of meta-CART, needs to be regression tree
  #  testdat  : the test data set that is used for the Q test
  #
  # Returns:
  # the pvalue of the Q-test. If there is no interaction detected in the final tree, then 
  # a NA value will be returned.
  testnodes <- predict(prunetree, testdat)  # make the subgrouping variable
  if (length(unique(testnodes)) > 2) {  # if there is interaction detected in the final tree
    res <- rma(yi = efk, vi = vark, 
               data = testdat, measure = "SMD", 
               mods = ~factor(testnodes), intercept=FALSE)
    return(res$QMp)  # compute the pvalue of the Q -test
  } else {
    return(NA)
  }
}
 
getpC <- function(prunetree, testdat){
  # compute the significance of moderator effects. the second step of meta-CART. 
  #
  # Arguments:
  #  prunetree: the final tree in the first step of meta-CART, needs to be regression tree
  #  testdat  : the test data set that is used for the Q test
  #
  # Returns:
  # the pvalue of the Q-test. If there is no moderator detected in the final tree, then 
  # a NA value will be returned.
  testnodes <- predict(prunetree, testdat, type="prob")
  testnodes <- testnodes[, 1]
  if (length(unique(testnodes)) > 2) {
    res <- rma(yi = efk, vi = vark, 
               data = testdat, measure = "SMD", 
               mods = ~factor(testnodes), intercept=FALSE)
    return(res$QMp)
  } else {
    return(NA)
  }
}
