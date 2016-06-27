source("1SimFuncsV5.R")
# ordinal data with 3 levels
set.seed(1552)
K <- c(40,  80 , 120)  # will be devided into training and test data
n.bar <- c(40, 80, 160)
tau <- c(0,  sqrt(0.025), sqrt(0.05))
mods <- c(5, 10, 20)
tree <- 1:3
cons <- expand.grid(K, n.bar, tau, mods, tree)
colnames(cons) <- c("K", "n.bar", "tau", "mods", "tree")

# Tree stuctures
formula <- c("~m1", "~ m1:m2", "~m1:m2:m3")
beta <- list( b1 = c(0, 0.5), b2 = c(0, 0.5),
             b3 = c(0, 0.5))
beta2 <- list( b1 = c(0, 0.8), b2 = c(0, 0.8),
               b3 = c(0, 0.8))
# beta2 <- list(b1 = 0.3, b2 = c(0.2, 0.6), b3 = c(0.2, 0.6),
#               b4 = c(0.2, 0.6, 0.6, -0.6), b5 = c(0.2, 0.6))
truemods <- list(mods1="m1", mods2 = c("m1", "m2"),
                 mods3 = c("m1", "m2", "m3"))
temp0 <-rpart.control(xval=10,minbucket=5,minsplit=10,cp=0.0001)
c <- 1
datrep <- list()  # data sets with medium interaction
for (i in 1:nrow(cons)) {  
  datcons <- list()
  print(i)
  for (j in 1:50){
    m <- cons[i, ]$mods
    lvls <- rep(3, m)
    K <- cons[i, ]$K
    tree <- cons[i, ]$tree
    tau <- cons[i, ]$tau
    n.bar <- cons[i, ]$n.bar
    ms.data <- ModsGen.ordinal(m, lvls, K)
    ms.formula <- as.data.frame((ms.data > 1 )*1)
    dat <- cbind(SimData(beta = beta[[tree]], mods = ms.formula, formula = as.formula(formula[tree]),
                         K = K, n = NumGen(n.bar, K), tau = tau), ms.data)
    datcons[[j]] <- dat
  }
  datrep[[i]] <- datcons
}

datrep2 <- list()  # data sets with medium interaction
for (i in 1:nrow(cons)) {  
  datcons <- list()
  print(i)
  for (j in 1:50){
    m <- cons[i, ]$mods
    lvls <- rep(3, m)
    K <- cons[i, ]$K
    tree <- cons[i, ]$tree
    tau <- cons[i, ]$tau
    n.bar <- cons[i, ]$n.bar
    ms.data <- ModsGen.ordinal(m, lvls, K)
    ms.formula <- as.data.frame((ms.data > 1 )*1)
    dat <- cbind(SimData(beta = beta2[[tree]], mods = ms.formula, formula = as.formula(formula[tree]),
                         K = K, n = NumGen(n.bar, K), tau = tau), ms.data)
    datcons[[j]] <- dat
  }
  datrep2[[i]] <- datcons
}



# m <- cons[1, ]$mods
# lvls <- rep(3, m)
# K <- cons[1, ]$K
# tree <- cons[1, ]$tree
# tau <- cons[1, ]$tau
# n.bar <- cons[1, ]$n.bar
# ms.data <- ModsGen.ordinal(m, lvls, K)
# ms.formula <- as.data.frame((ms.data > 1 )*1)
# dat <- cbind(SimData(beta = beta[[tree]], mods = ms.formula, formula = as.formula(formula[tree]),
#         K = K, n = NumGen(n.bar, K), tau = tau), ms.data)
# FEres <- rpart(efk~.-trail-vark, data= dat, weights = 1/vark, control = temp0)
# FEtree <- treepruner(FEres, sqrt(mean(1/dat$vark)))

save.image("DataOrd")
