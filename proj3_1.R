#lmm
lmm <- function(form, dat, ref = list()) {
  # Set up model matrices and other initial values using LMMsetup function
  setup <- LMMsetup(form, dat, ref)
  
  # Initialize theta with zeros, with length equal to the number of random effect groups
  theta_init <- rep(0, 1+ length(setup$ref)) # 包括log(sigma)和每个随机效应
  
  # Use optim to find the maximum likelihood estimate (MLE) of theta
  opt <- optim(par = theta_init, fn = LMMprof, setup = setup, method = "BFGS")
  
  # Extract optimized theta and compute corresponding beta
  theta <- opt$par
  nlog_lik <- LMMprof(theta,setup)
  beta <- attr(LMMprof(theta, setup), "beta")
  
  # Return the estimated beta and theta values
  return(list(beta = beta, theta = theta))
}

  #LMMsetup
LMMsetup <- function(form, dat, ref = list()) {
  
  # 将公式转换为公式对象
  form <- as.formula(form)
  
  # 设置固定效应矩阵 X 和响应向量 y
  X <- model.matrix(form, dat)
  y <- model.response(model.frame(form, dat))
    
  # 设置随机效应矩阵 Z
  # 对 ref 中的每个变量名向量，创建相应的模型矩阵
  Z_list <- lapply(ref, function(vars) {
    if (length(vars) == 1) {
      
      #如果仅有一个变量，创建无截距的模型矩阵
      model.matrix(as.formula(paste("~", vars, "-1")), dat)
      
      } else {
        # 如果有多个变量，创建这些变量的交互项模型矩阵（无截距）
        model.matrix(as.formula(paste("~", paste(vars, collapse = ":"), "-1")), dat)
      }
    })
    
  # 将所有随机效应模型矩阵合并为一个矩阵 Z
    if (length(Z_list) > 0){
      Z <- do.call(cbind, Z_list)
    }else{
      Z <- matrix(0,nrow = nrow(X),ncol = 0) # 若无随机效应，则Z为空矩阵
    }
   
    # 对Z进行QR分解
    qrz <- qr(Z)
    Q <- qr.Q(qrz) # 提取 Q 矩阵
    p <- qrz$rank
    R_full <- qr.R(qrz) # 提取 R 矩阵
    R <- R_full[1:p,]    
    
    # 记录样本数n和Z的秩p
    
    n <- nrow(X)
    
    # Return the fixed effect matrix X, response vector y, random effect matrix Z, and ref list
    return(list(X = X, y = y, Z = Z, Zblocks = Z_list, ref = ref, qrz = qrz, R = R, n = n, p = p))
}

  
# LMM prof
LMMprof <- function(theta, setup) {
  # 提取 sigma
  sigma <- exp(theta[1])
  variance <- sigma^2  # 方差
  
  # 提取随机效应的标准差和方差
  psi <- exp(theta[-1])
  psi2 <- psi^2  # 方差
  
  # 构建随机效应的协方差矩阵
  psi_theta <- list()
  for (i in seq_along(setup$Zblocks)) {
    z_i <- setup$Zblocks[[i]]  # 获取第 i 个 Z 块
    z_i_ncol <- ncol(z_i)      # Z 的列数，表示有多少随机效应
    psi_theta[[i]] <- diag(psi2[i], z_i_ncol, z_i_ncol)  # 创建 k*k 的对角矩阵
  }
  
   p <- ncol(setup$R)
  # 使用 bdiag 构建 psi_theta_mat
  library(Matrix)
  if (length(psi_theta) > 0) {
    psi_theta_mat <- as.matrix(bdiag(psi_theta))
    R <- setup$R  # 获取调整后的 R 矩阵
    U <- R %*% psi_theta_mat %*% t(R) + variance * diag(p)
  } else {
    U <- variance * diag(p)
  }
  
  # 对 U 进行 Cholesky 分解，返回下三角矩阵
  L <- chol(U)
  
  # 计算 log|U|
  log_det_U <- 2 * sum(log(diag(L)))
  
  # 计算 Q^T y 和 Q^T X
  y <- setup$y
  x <- setup$X
  qrz <- setup$qrz
  qty <- qr.qty(qrz, y)  # Q^T %*% y
  qtx <- qr.qty(qrz, x)  # Q^T %*% X
  
  # 分割 qty 和 qtx
  p <- setup$p
  qty1 <- qty[1:p]
  qty2 <- qty[(p+1):length(qty)]
  qtx1 <- qtx[1:p, , drop = FALSE]
  qtx2 <- qtx[(p+1):nrow(qtx), , drop = FALSE]
  
  # 求解 wy (v = wy)
  v1 <- forwardsolve(L, qty1)
  
  wy1 <- v1
  wy2 <- qty2 / variance
  
  # 求解 wx
  wx1 <- matrix(0, nrow = nrow(qtx1), ncol = ncol(qtx1))
  for (i in 1:ncol(qtx1)) {
    v2 <- forwardsolve(L, qtx1[, i])
    wx1[, i] <- v2
  }
  wx2 <- qtx2 / variance
  
  # 合并 wy 和 wx
  wy <- c(wy1, wy2)
  wx <- rbind(wx1, wx2)
  
  # 计算 xtwx 和 xtwy
  xtwx <- t(x) %*% wx
  xtwy <- t(x) %*% wy
  
  # 计算 beta_hat
  r <- chol(xtwx)
  rv3 <- forwardsolve(t(r), xtwy)
  beta_hat <- backsolve(r, rv3)
  
  # 计算残差 e = y - x %*% beta_hat
  e <- y - x %*% beta_hat
  
  # 计算 Q^T e
  qte <- qr.qty(qrz, e)
  qte1 <- qte[1:p]
  qte2 <- qte[(p+1):length(qte)]
  
  v4 <- forwardsolve(L, qte1)
  
  etwe <- sum(v4^2) + sum((qte2 / sigma)^2)
  
  # 计算负对数似然值
  n <- setup$n
  nlog_lik <- (etwe) / 2 + (n - p) * log(sigma) + log_det_U / 2
  attr(nlog_lik, "beta") <- beta_hat
  return(nlog_lik)
}


library(nlme);library(lme4)
lmm(score ~ Machine,Machines,list("Worker",c("Worker","Machine")))
lmer(score ~ Machine + (1|Worker) + (1|Worker:Machine),data=Machines,
     REML=FALSE)
  
  
  