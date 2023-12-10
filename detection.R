#Data generation
models <- function(family = c("gaussian", "binomial", "poisson"), type = c("all", "source", "target"), cov.type = 1, h = 5, K = 5, n.target = 200, n.source = rep(100, K), s = 5, p = 500, Ka = K) {
  family <- match.arg(family)
  target <- NULL
  source <- NULL
  
  type <- match.arg(type)
  sig.strength <- 0.5
  
  if (family == "gaussian" || family == "binomial") {
    if(type == "all" || type == "target") {
      wk <- c(rep(sig.strength, s), rep(0, p-s))
      if (cov.type == 1) {
        Sigma <- outer(1:p, 1:p, function(x,y){
          0.5^(abs(x-y))
        })
        R <- chol(Sigma)
        target <- list(x = NULL, y = NULL)
        
        target$x <- matrix(rnorm(n.target*p), nrow = n.target) %*% R
        
      } else if (cov.type == 2) {
        Sigma <- outer(1:p, 1:p, function(x,y){
          0.9^(abs(x-y))
        })
        R <- chol(Sigma)
        target$x <- matrix(rnorm(n.target*p), nrow = n.target) %*% R
      }
      
      if (family == "gaussian") {
        target$y <- as.numeric(target$x %*% wk + rnorm(n.target))
      } else if (family == "binomial") {
        pr <- 1/(1+exp(-target$x %*% wk))
        target$y <- sapply(1:n.target, function(i){sample(0:1, size = 1, prob = c(1-pr[i], pr[i]))})
      }
    }
    
    if(type == "all" || type == "source") {
      if (cov.type == 1) {
        Sigma <- outer(1:p, 1:p, function(x,y){
          0.5^(abs(x-y))
        })
        eps <- rnorm(p, sd = 0.3)
        Sigma <- Sigma + eps %*% t(eps)
        
        R <- chol(Sigma)
      }
      
      source <- sapply(1:K, function(k){
        if (k <= Ka){
          wk <- c(rep(sig.strength, s), rep(0, p-s)) + h/p*sample(c(-1,1), size = p, replace = TRUE)
        } else {
          sig.index <- c(s+1:s, sample((2*s+1):p, s))
          wk <- rep(0, p)
          wk[sig.index] <- sig.strength
          wk <- wk + 2*h/p*sample(c(-1,1), size = p, replace = TRUE)
        }
        if (cov.type == 1) {
          x <- matrix(rnorm(n.source[k]*p), nrow = n.source[k]) %*% R
        } else if (cov.type == 2) {
          x <- matrix(rt(n.source[k]*p, df = 4), nrow = n.source[k])
        }
        
        if (family == "gaussian") {
          y <- as.numeric(0.5*I(k > Ka) + x %*% wk + rnorm(n.source[k]))
        } else if (family == "binomial") {
          pr <- 1/(1+exp(-0.5*I(k > Ka) -x %*% wk))
          y <- sapply(1:n.source[k], function(i){
            sample(0:1, size = 1, prob = c(1-pr[i], pr[i]))
          })
        }
        list(x = x, y = y)
      }, simplify = FALSE)
    }
    
    
  } else { # model == "poisson
    if(type == "all" || type == "target") {
      if (cov.type == 1) {
        Sigma <- outer(1:p, 1:p, function(x,y){
          0.5^(abs(x-y))
        })
        R <- chol(Sigma)
      } else if (cov.type == 2) {
        Sigma <- outer(1:p, 1:p, function(x,y){
          0.9^(abs(x-y))
        })
        R <- chol(Sigma)
      }
      wk <- c(rep(sig.strength, s), rep(0, p-s))
      
      target$x <- matrix(rnorm(n.target*p), nrow = n.target) %*% R
      
      target$x[target$x > 0.5] <- 0.5
      target$x[target$x < -0.5] <- -0.5
      lambda <- as.numeric(exp(target$x %*% wk))
      target$y <- rpois(n.target, lambda)
    }
    
    if(type == "all" || type == "source") {
      if (cov.type == 1) {
        Sigma <- outer(1:p, 1:p, function(x,y){
          0.5^(abs(x-y))
        })
        eps <- rnorm(p, sd = 0.3)
        Sigma <- Sigma + eps %*% t(eps)
        
        R <- chol(Sigma)
      }
      
      source <- sapply(1:K, function(k){
        if (k <= Ka){
          wk <- c(rep(sig.strength, s), rep(0, p-s)) + h/p*sample(c(-1,1), size = p, replace = TRUE)
        } else {
          sig.index <- c(s+1:s, sample((2*s+1):p, s))
          wk <- rep(0, p)
          wk[sig.index] <- sig.strength
          wk <- wk + 2*h/p*sample(c(-1,1), size = p, replace = TRUE)
        }
        if (cov.type == 1) {
          x <- matrix(rnorm(n.source[k]*p), nrow = n.source[k]) %*% R
        } else if (cov.type == 2) {
          x <- matrix(rt(n.source[k]*p, df = 4), nrow = n.source[k])
        }
        x[x > 0.5] <- 0.5
        x[x < -0.5] <- -0.5
        lambda <- as.numeric(exp(0.5*I(k > Ka) + x %*% wk))
        y <- rpois(n.source[k], lambda)
        
        list(x = x, y = y)
      }, simplify = FALSE)
    }
    
  }
  
  
  if (type == "all") {
    return(list(target = target, source = source))
  } else if (type == "target") {
    return(list(target = target))
  } else {
    return(list(source = source))
  }
  
}
#Algorithm 2
library(caret)
valid.nfolds=3
folds <- createFolds(target$y,valid.nfolds)#进行索引
loss.cv.source <- sapply(1:valid.nfolds, function(i){
  source.loss <- sapply(1:length(source), function(k){
    yc = c(target$y[-folds[[i]]], source[[k]]$y)
    yp = array(yc)
    xc = as.matrix(rbind(target$x[-folds[[i]], , drop = F], source[[k]]$x))
    A <- t(xc)-colMeans(xc)
    B1 <- t(A)
    B1 <- np_array(B1)
    V <- os$cqrpadmm(B1,yp,intercept = FALSE)
    tau1 = c(0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95)
    tau1=array(tau1)
    F <- V$cqrp_admm_smw(tau=tau1,maxit=20000L)
    s2 = F['beta']
    wa <- as.numeric(unlist(s2))
    xv <- as.matrix(target$x[folds[[i]], , drop = F])
    yv <- target$y[folds[[i]]]
    xcenter1 <- t(xv)-colMeans(xv)
    xcenter <- t(xcenter1)
    loss(wa,xcenter,yv,tau)})})
loss.cv.target <- sapply(1:valid.nfolds, function(i){
  yd = target$y[-folds[[i]]]
  yp = array(yd)
  xd = as.matrix(target$x[-folds[[i]], , drop = F])
  A <- t(xd)-colMeans(xd)
  B1 <- t(A)
  B1 <- np_array(B1)
  V <- os$cqrpadmm(B1,yp,intercept = FALSE)
  tau = c(0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95)
  tau2=array(tau)
  F <- V$cqrp_admm_smw(tau=tau2,maxit=20000L)
  s2 = F['beta']
  wa.target <- as.numeric(unlist(s2))
  xm <- as.matrix(target$x[folds[[i]], , drop = F])
  ym <- target$y[folds[[i]]]
  xcenter2 <- t(xm)-colMeans(xm)
  xcentert <- t(xcenter2)
  targrt.loss <- loss(wa.target,xcentert,ym,tau)
}) 
loss.cv.source1 <- t(loss.cv.source)
dim(loss.cv.source1)
dim(loss.cv.target)
loss.cv <- cbind(loss.cv.source1,loss.cv.target)
source.loss <- colMeans(loss.cv)[1:(ncol(loss.cv)-1)]
target.valid.loss <- colMeans(loss.cv)[ncol(loss.cv)]
target.valid.loss.sd <- sd(loss.cv[, ncol(loss.cv)])
C0 <- 0.01
threshold <- (1+C0)*target.valid.loss
transfer.source.id <- which(source.loss <= threshold)
transfer.source.id
#Algorithm 3
selection <- function(target, source = NULL, family = c("gaussian", "binomial", "poisson"),
                      standardize = TRUE,valid.nfolds = valid.nfolds, intercept = FALSE, nfolds = 10)
{folds <- createFolds(target$y, valid.nfolds)
#交叉求损失得到索引结果
loss.cv.source <- sapply(1:valid.nfolds, function(i){
  source.loss <- sapply(1:length(source), function(k){
    yc = c(target$y[-folds[[i]]], source[[k]]$y)
    yp = array(yc)
    xc = as.matrix(rbind(target$x[-folds[[i]], , drop = F], source[[k]]$x))
    A <- t(xc)-colMeans(xc)
    B1 <- t(A)
    B1 <- np_array(B1)
    V <- os$cqrpadmm(B1,yp,intercept = FALSE)
    tau = c(0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95)
    tau1=array(tau)
    F <- V$cqrp_admm_smw(tau=tau1,maxit=20000L)
    s2 = F['beta']
    wa <- as.numeric(unlist(s2))
    xv <- as.matrix(target$x[folds[[i]], , drop = F])
    yv <- target$y[folds[[i]]]
    xcenter1 <- t(xv)-colMeans(xv)
    xcenter <- t(xcenter1)
    loss(wa, xcenter, yv,tau)})})
loss.cv.target <- sapply(1:valid.nfolds, function(i){
  yd = target$y[-folds[[i]]]
  yp = array(yd)
  xd = as.matrix(target$x[-folds[[i]], , drop = F])
  A <- t(xd)-colMeans(xd)
  B1 <- t(A)
  B1 <- np_array(B1)
  V <- os$cqrpadmm(B1,yp,intercept = FALSE)
  tau = c(0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95)
  tau2=array(tau)
  F <- V$cqrp_admm_smw(tau=tau2,maxit=20000L)
  s2 = F['beta']
  wa.target <- as.numeric(unlist(s2))
  xm <- as.matrix(target$x[folds[[i]], , drop = F])
  ym <- target$y[folds[[i]]]
  xcenter2 <- t(xm)-colMeans(xm)
  xcentert <- t(xcenter2)
  targrt.loss <- loss(wa.target,xcentert,ym,tau)
})
loss.cv.source1 <- t(loss.cv.source)
loss.cv <- cbind(loss.cv.source1,loss.cv.target)
source.loss <- colMeans(loss.cv)[1:(ncol(loss.cv)-1)]
target.valid.loss <- colMeans(loss.cv)[ncol(loss.cv)]
if (min(source.loss) >= target.valid.loss){
  paste0("There is no transferable source", "\n")
}else{
  transfer.source.firstid = which.min(source.loss)
}
source_list = c(transfer.source.firstid)
source_loss_in = 0
while (TRUE){
  #source.loss = source_loss_new
  print(c("source.loss",source.loss))
  loss_invalue = c();loss_inkey = c()
  loss_outvalue = c();loss_outkey = c()
  # 计算引入源域的loss
  #source_list = c()
  for (k in 1:length(source)){
    if (!k %in% source_list){
      transfer.source.x = rbind(source[[k]]$x,as.matrix(foreach(source_index=source_list,.combine="rbind") %do% {
        source[[source_index]]$x
      }))
      transfer.source.y = c(source[[k]]$y,as.matrix(foreach(source_index=source_list,.combine="cbind") %do% {
        source[[source_index]]$y
      }))
      #测试挑选出的迁移域的效果
      #source_list=c(1,2,3,4,5,6)
      #transfer.source.x = as.matrix(foreach(source_index=source_list,.combine="rbind") %do% {
      #source[[source_index]]$x
      #})
      #transfer.source.y = c(as.matrix(foreach(source_index=source_list,.combine="cbind") %do% {
      #source[[source_index]]$y
      #}))
      target.valid.id <- sample(1:length(target$y), size = floor(length(target$y)*(1/2)))#分为训练集和测试集2个
      ym = c(target$y[-target.valid.id], transfer.source.y)
      yl = array(ym)
      xm = as.matrix(rbind(target$x[-target.valid.id, , drop = F], transfer.source.x))
      Am <- t(xm)-colMeans(xm)
      Bm <- t(Am)
      Bm <- np_array(Bm)
      Vm <- os$cqrpadmm(Bm,yl,intercept = FALSE)
      tau = c(0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95)
      taum=array(tau)
      F <- Vm$cqrp_admm_smw(tau=taum,maxit=20000L)
      sm = F['beta']
      wa <- as.numeric(unlist(sm))
      xm <- as.matrix(target$x[target.valid.id, , drop = F])
      ym <- target$y[target.valid.id]
      xcenterm <- t(xm)-colMeans(xm)
      xcenter <- t(xcenterm)
      loss(wa, xcenter, ym,tau)
      loss_invalue[k] = loss(wa, xcenter, ym, tau)
    }else{
      loss_invalue[k] = NA
    }
  }
  #print(c("loss_invalue : ",loss_invalue))
  source_loss_in = min(na.omit(loss_invalue))
  #print(c("source loss_in : ",source_loss_in))
  if (source_loss_in < min(source.loss)){
    source_list = append(source_list,which.min(loss_invalue))
    source.loss = source_loss_in
    # 计算已引入的影响（是否剔除）
    for (k_index in source_list[-length(source_list)]){
      transfer.source.x = as.matrix(foreach(source_index = source_list[-which(source_list == k_index)],.combine="rbind") %do% {
        source[[source_index]]$x
      })
      transfer.source.y = as.matrix(foreach(source_index = source_list[-which(source_list == k_index)],.combine="cbind") %do% {
        source[[source_index]]$y
      })
      yn = c(target$y[-target.valid.id], transfer.source.y)
      yl = array(yn)
      xn = as.matrix(rbind(target$x[-target.valid.id, , drop = F], transfer.source.x))
      An <- t(xn)-colMeans(xn)
      Bn <- t(An)
      Bn <- np_array(Bn)
      Vn <- os$cqrpadmm(Bn,yl,intercept = FALSE)
      tau = c(0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95)
      taun=array(tau)
      F <- Vn$cqrp_admm_smw(tau=taun,maxit=20000L)
      sn = F['beta']
      wa <- as.numeric(unlist(sn))
      xn <- as.matrix(target$x[target.valid.id, , drop = F])
      yn <- target$y[target.valid.id]
      xcentern <- t(xn)-colMeans(xn)
      xcenter <- t(xcentern)
      loss_outvalue[k_index] = loss(wa, xcenter, yn, tau)
    }
    #print(c("loss_outvalue : ",loss_outvalue))
    source_loss_out = min(na.omit(loss_outvalue))
    #print(c("source loss_out : ",source_loss_out))
    if (source_loss_out <= source_loss_in){
      #source_list = source_list[-which(loss_outvalue==which.min(loss_outvalue))]
      source_list = source_list[-which.min(loss_outvalue)]
      source.loss = source_loss_out
    }
  }else{
    break
  }
}
return(source_list)
}
selection(target = target, source = source,family = "gaussian",standardize = TRUE, valid.nfolds = 2,intercept = FALSE)
