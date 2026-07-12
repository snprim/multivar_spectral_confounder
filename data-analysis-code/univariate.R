# full model with inverse gamma prior for beta
# lambda1 \sim invgamma(1/2, 1/gamma1) # this is lambdaKHC
# lambda2 \sim invgamma(1/2, 1/gamma2) # this is lambdaEHC
# lambda3 (lambdaRHC) is the same as tauHC
# allow for covariates (Z) and treatments (X)
# save theta (AU)
# univariate

MCMC_model4_9_XZ <- function(Y, X, Z, Qmat, E, n, R, Q, L, K, B, iters = 1000, burn = 50, W, thin = 2){
  tik <- proc.time()
  R <- Q <- 1
  iters_thin <- (iters-burn)%/%thin
  # keepers_beta <- array(0, c(n, E, R, iters_thin))
  # the draws group by r (beta, A, tau2) and by q (U and sig2)
  # draws are stored in 3d arrays of the dimension rows by cols by iters
  # keepers_beta <- array(0, c(n, E, R, iters))
  # keepers_gamma <- array(0, c(L, E, R, iters))
  p <- ncol(Z)
  # keepers_alpha <- array(0, c(p, R, iters))
  # keepers_AU <- array(0, c(n,R,iters))
  att <- acc <- 0
  # priors and initial values
  MH <- 0.1
  # W<- eigen(Qmat)$values
  a_tau <- .5
  b_tau <- .005
  Umat <- matrix(0, nrow = n, ncol = Q)
  keepers_Tl <- array(0, c(L, K, iters_thin))
  keepers_Te <- array(0, c(E, K, iters_thin))
  keepers_U <- array(0, c(n, Q, iters_thin))
  Tl <- matrix(0, nrow = L, ncol = K)
  Te <- matrix(0, nrow = E, ncol = K)
  sig2 <- rep(0.5, Q)
  tau2 <- rep(0.5, R)
  lambda <- 0.5
  lambdaKHC <- 1
  lambdaEHC <- 1
  tauHC <- 1
  tauQHC <- 1
  nuKHC <- nuEHC <- xiHC <- 1
  estimate <- find_estimate3(X, B, Tl, Te)
  alpha <- matrix(0, nrow = p, ncol = R)
  pb <- txtProgressBar(min = 2, max = iters, style = 3)
  for (iter in 2:iters){
    # print(paste0("This is the ", iter, "th iterations."))
    # update tensor margins for beta
    Zalpha <- Z%*%alpha
    for (k in 1:K){
      for (l in 1:L){
        Cl <- B[,l] * (X %*% Te[,k])
        Vl <- sum(colSums(Cl^2)/tau2) + 1/lambdaKHC
        ### Vl <- sum(colSums(Cl^2)/tau2) + 1/10
        Tl_prev <- Tl[l,k]
        Ml <- sum(colSums(Cl*(Y - Umat - Zalpha - estimate + Tl_prev*Cl))/tau2)
        Tl[l,k] <- rnorm(1, Ml/Vl, 1/sqrt(Vl))
        # estimate <- find_estimate(X_tensor, B, Tl, Te, Tr, n, R)
        estimate <- estimate + (Tl[l,k] - Tl_prev)*Cl 
      }
      for (e in 1:E){
        Ce <- X[,e] * (B %*% Tl[,k])
        Ve <- sum(colSums(Ce^2)/tau2) + 1/lambdaEHC
        ### Ve <- sum(colSums(Ce^2)/tau2) + 1/10
        Te_prev <- Te[e,k]
        Me <- sum(colSums(Ce*(Y - Umat - Zalpha - estimate + Te_prev*Ce))/tau2)
        Te[e,k] <- rnorm(1, Me/Ve, 1/sqrt(Ve))
        # estimate <- find_estimate(X_tensor, B, Tl, Te, Tr, n, R)
        estimate <- estimate + (Te[e,k] - Te_prev)*Ce
      }
      
    }
    # update shrinkage prior parameters
    sumbyl <- rowSums(Tl^2/2)
    sumbye <- rowSums(Te^2/2)
    
    lambdaKHC <- rinvgamma(1, (L*K+1)/2, sum(sumbyl) + 1/nuKHC)
    lambdaEHC <- rinvgamma(1, (E*K+1)/2, sum(sumbye)/tauHC + 1/nuEHC)
    tauHC <- rinvgamma(1, (K*E+1)/2, sum(sumbye/lambdaEHC) + 1/xiHC)
    
    nuKHC <- rinvgamma(1, 1, 1/lambdaKHC + 1)
    nuEHC <- rinvgamma(1, 1, 1/(lambdaEHC*tauHC) + 1)
    xiHC <- rinvgamma(1, 1, 1/tauHC + 1)
    
    # tauKHC <- rinvgamma(1, (K+1)/2, sum(1/lambdaKHC) + 1)
    # tauEHC <- rinvgamma(1, (E+1)/2, sum(1/lambdaEHC) + 1)
    
    # update tau2_r
    tau2 <- rinvgamma(1, (n+K)/2 + a_tau, 0.5*sum((Y-Zalpha-estimate-Umat)^2) + b_tau)
    
  
    # update alpha_r
    tZK <- sweep(t(Z), 2, 1/tau2, "*")
    AAZ <- tZK %*% Z + diag(1/100, p)
    BBZ <- tZK %*% (Y - Umat - estimate)
    alpha <- t(spam::rmvnorm.canonical(1, BBZ, AAZ))
    
    
    # update sig2Q
    Kinv_q <- makeKinv(Umat, lambda, W)
    # sig2[j] <- rinvgamma(1, n/2+a_sig2, Kinv_q/2+b_sig2)
    sig2 <- rinvgamma(1, (n+1)/2, Kinv_q/2 + 1/tauQHC)
    # update Uq
    sumR <- (Y-Zalpha-estimate)/tau2
    vark <- 1/tau2 + (1-lambda+lambda*W)/sig2
      # Qindex <- setdiff(1:Q, j)
      # sumR <- rep(0, n)
      # sumArqdivTaur <- 0
      # for (w in 1:R){
      #   sumP <- rep(0, n)
      #   for (qprime in Qindex){
      #     sumP <- sumP + Umat[,qprime] * Amat[w,qprime]
      #   }
      #   sumArqdivTaur <- sumArqdivTaur + ((Amat[w,j])^2/tau2[w])
      #   sumR <- sumR + ((Y[,w]-Zalpha[,w]-estimate[,w]-sumP)*Amat[w,j]/tau2[w])
      # }
      # vark <- sumArqdivTaur+(1-lambda+lambda*W)/sig2[j]
    for (d in 1:n){
      Umat[d] <- rnorm(1, sumR[d]/vark[d], 1/sqrt(vark[d]))
    }
  
    tauQHC <- rinvgamma(1, (Q+1)/2, 1/sig2 + 1)
    # update lambdaU
    att <- att + 1
    can <- pnorm(rnorm(1,qnorm(lambda), MH))
    # calculate current and candidate loglikelihood
    curlp <- 0
    canlp <- 0

    curK <- makeK(sig2, lambda, W)
    canK <- makeK(sig2, can, W)
    curKinv <- makeKinv(Umat, lambda, W)
    canKinv <- makeKinv(Umat, can, W)
    curlp <- curlp + log_like(curK, curKinv, sig2)
    canlp <- canlp + log_like(canK, canKinv, sig2)
  
    Rval <- canlp - curlp + dnorm(qnorm(can), log = TRUE) - dnorm(qnorm(lambda), log = TRUE)
    if (!is.na(Rval) & log(runif(1)) < Rval & lambda < .9999){
      acc <- acc + 1
      lambda <- can
    }
    if(iter < burn){
      if(att > 50){
        if(acc/att < 0.3){MH <- MH*0.8}
        if(acc/att > 0.5){MH <- MH*1.2}
        acc <- att <- 0
      }
    }
    
    if (iter > burn){
      if ((iter-burn)%%thin == 0){
        # keepers_beta[,,,(iter-burn)%/%thin] <- bbygamma(B, find_tol_gamma(Tl, Te, Tr), n)  
        keepers_Tl[,,(iter-burn)%/%thin] <- Tl
        keepers_Te[,,(iter-burn)%/%thin] <- Te
        keepers_U[,,(iter-burn)%/%thin] <- Umat
      }
    }
    
    # Sleep for 0.1 seconds
    Sys.sleep(0.01)
    
    # Print progress
    setTxtProgressBar(pb, iter)
    
    # keepers_gamma[,,,iter] <- find_tol_gamma(Tl, Te, Tr)
    # keepers_beta[,,,iter] <- bbygamma(B, find_tol_gamma(Tl, Te, Tr), n)
    # keepers_beta[,,,iter] <- bbygamma(B, find_tol_gamma(Tl, Te, Tr), n)
    # keepers_alpha[,,iter] <- alpha
    # keepers_AU[,,iter] <- AU
    
  }
  close(pb)
  
  tok <- proc.time()
  out <- list(Tl = keepers_Tl, Te = keepers_Te, U = keepers_U, acc_rate = acc/att, time = tok - tik, MH = MH)
  return(out)
}
