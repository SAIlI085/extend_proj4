# Contribution:
# Sai Li s2714242; Bingxin Li s2723681; Yuqi Tan s2703761
# Contributions: Sai Li : Bingxin Li : Yuqi Tan = 50%:50%:0%
# Sai Li & Bingxin Li work together to explored implementation options, completed code writing, debugging and commenting


# Main function LMM
lmm <- function(form, dat, ref = list()) {
  # LMM：Main function to estimate parameters in a linear mixed model
  # Arguments:
    # form: model formula to specify the fixed effect of the model
    # dat: the data frame containing all variables in the model
    # ref: a list of vectors specifying the random effects
  # Returns:
    # A list containing the likelihood estimates for beta (fixed effects) and theta(random effects)
  
  # Set up and processing the data using LMMsetup function
  setup <- LMMsetup(form, dat, ref)     # a list containing model matrix and initial components
  
  # Optimize the negative log-likelihood function to find the maximum likelihood estimate of theta
  opt <- optim(par = setup$theta_init,  # Initial values for theta
               fn = LMMprof,            # Function to compute negative log likelihood
               setup = setup,           # Additional data setup required by LMMprof
               ref = ref,               # Random effects
               method = "BFGS")         # Optimization method
  
  # Compute LMMprof with the optimized theta to generate the estimated beta coefficients
  # final_neg_log_likelihood stores the value of the negative log-likelihood function at the optimized theta
  final_neg_log_likelihood <- LMMprof(opt$par, setup, dat, ref)
  
  # Add the attribute of the final_neg_log_likelihood
  beta <- attr(final_neg_log_likelihood, "beta")
  
  # Return a list containing the likelihood estimates
  list(beta = beta,       # estimated coefficients for fixed effects
       theta = opt$par)   # optimized values of theta, representing variance parameters
}


# LMMsetup
LMMsetup <- function(form, dat, ref) {
  # LMMsetup: to set up the model matrix X,Z and some parameters
    # Arguments:
      # form: model formula to specify the fixed effect of the model
      # dat: the data frame containing all variables in the model
      # ref: a list of vectors specifying the random effects
    # Returns:
      # a list containing X,Z,R,y,initial theta and the status of random effects
  
  # Extract response variable y from data
  y <- dat[[all.vars(form)[1]]]  # the first variable in the formula is y
  
  # Create the fixed effects matrix
  X <- model.matrix(form, dat)

  # Check if it exists random effects
  if (length(ref) == 0) {  # Case with no random effects
    Z <- matrix(0, nrow = nrow(X), ncol = 0) # Set Z to an empty matrix with zero columns
    R <- matrix(0, nrow = 0, ncol = 0) # Set R to an empty matrix
    theta_init <- log(var(y))          # initialize theta as log variance of y
    # Return the list with parameters with no random effects 
    return(list(X = X, qr_decomp = NULL, R = R, y = y, theta_init = theta_init, ref_empty = TRUE))
  }
  
  # If there exists random effects
  Z_list <- lapply(ref, function(variables) {
    # If only exists one parameter with random effects
    if (length(variables) == 1) {
      model.matrix(as.formula(paste("~", variables, "- 1")), dat) # Create model matrix with no intercept
    } else {
      # If more than one parameter, create an interaction term for the model matrix
      model.matrix(as.formula(paste("~", paste(variables, collapse = ":"), "- 1")), dat)
    }
  })
  Z <- do.call(cbind, Z_list) # Combine all Z matrices to full Z matrix for random effects
  
  # Advance processing of Z-matrix using QR decomposition
  qr_decomp <- qr(Z)   # QR decomposition
  R <- qr.R(qr_decomp) # Extract the R matrix from QR decomposition of Z
  
  # Initialize theta
  # The first element is the log variance of y followed by zero for each random effect
  theta_init <- c(log(var(y)), rep(0, length(ref)))
  
  # Return a list containing all the set up components needed for next steps
  list(X = X,  # Fixed effects design matrix
       qr_decomp = qr_decomp,  # QR decomposition of Z
       R = R,  # Upper triangular matrix from QR decomposition of Z
       y = y,  # Response variables
       theta_init = theta_init,  # Initial values for theta
       ref_empty = FALSE) # the current status of random effects
}


#LMMprof
LMMprof <- function(theta, setup, dat = Machines, ref) {
  # LMMprof: Computes the negative log-likelihood for a linear mixed model (LMM)
    # Arguments:
      #   - theta: Numeric vector of model parameters
      #       - theta[1]: Log-transformed standard deviation of the residuals
      #       - theta[2:length(theta)]: Log-transformed standard deviations for random effects
      #   - setup: List of model setup components including:
      #       - X: Design matrix for fixed effects
      #       - y: Response vector
      #       - qr_decomp: QR decomposition of the design matrix
      #       - R: Upper triangular matrix from QR decomposition, aids in efficient computations
      #       - ref_empty: Logical indicating if there are no random effects
      #   - dat (optional): Data frame containing the model data
      #   - ref: List of grouping factors for the random effects structure, used to determine covariance structure
    # Returns:
      #   - Returns a single numeric value representing the negative log-likelihood of the model given parameters in theta.
      #   - Adds the estimated fixed effects coefficients as an attribute to the output.
  
  # Extract the parameters in the setup
  X <- setup$X            # Design matrix for the model
  y <- setup$y            # Response variable
  n <- nrow(X)            # Number of observations
  
  # Extract residual variance parameters from theta
  sigma <- exp(theta[1])   # Standard deviation of the residuals
  sigma2 <- sigma^2        # Variance of the residuals
  
  if (setup$ref_empty) {  # Case with no random effects
    # Directly construct the weight matrix W as a diagonal matrix with 1/sigma2 on the diagonal
    W <- diag(1 / sigma2, n)
    
    # Compute X^TWX and X^TWy 
    xTwx <- crossprod(X, W %*% X)  
    xTwy <- crossprod(X, W %*% y)  
    
    # Solve for beta_hat using Cholesky decomposition
    A <- xTwx
    B <- xTwy
    beta_hat <- backsolve(chol(A), forwardsolve(t(chol(A)), B))
    
    # Calculate the residuals using resid = y - X*beta
    resid <- y - X %*% beta_hat  
    CTwC <- crossprod(resid, W %*% resid)  # Calculated the weighted sum of squared residuals
                                           # resid^T %*% W %*% resid
    
    # Calculate the log determinant of W
    log_likelihood_det <- n * log(sigma2)
    
    # Compute the negative log-likelihood
    negative_log_likelihood <- CTwC / 2 + log_likelihood_det / 2
    
  } else {                                 # Case with random effects
    qr_decomp <- setup$qr_decomp           # QR decomposition of Z for efficient computations
    R <- setup$R                           # Upper triangular matrix from QR decomposition
    
    # Calculate the length of psi for each grouping factor in 'ref'
    length_of_psi <- sapply(ref, function(vars) {
      prod(sapply(vars, function(x) nlevels(dat[[x]])))  # Product of levels in each factor
    })
    # Construct psi_theta, scaling each theta component by its respective factor level count
    psi_theta <- rep(exp(theta[2:length(theta)] * 2), times = length_of_psi)
    psi <- diag(psi_theta, length(psi_theta), length(psi_theta))  # Diagonal matrix with psi values
    
    # Set up matrix dimensions
    p <- nrow(R)                       # Rank of R matrix
    q <- n - p                         # Degrees of freedom in the residual
    Ip <- diag(1, nrow = p, ncol = p)  # Identity matrix of size p
    Iq <- diag(1, nrow = q, ncol = q)  # Identity matrix of size (n - p)
    IQ <- Iq / sigma2                  # Rescaled Iq matrix
    
    # Calculate U = R * psi * R^T + Ip * sigma2
    U <- R %*% psi %*% t(R) + Ip * sigma2
    
    # Cholesky decomposition of U
    LT <- chol(U)
    U_inverse <- forwardsolve(LT, backsolve(t(LT), Ip))  # Inverse of U using Cholesky decomposition
    
    # Compute XtWX (X'WX), part of the weighted least squares solution
    Qt_X <- qr.qty(qr_decomp, X)            # Transformed X using QR decomposition
    x1 <- Qt_X[1:p, , drop = FALSE]         # Top part of transformed X
    x2 <- Qt_X[(p + 1):n, , drop = FALSE]   # Bottom part of transformed X
    # Block matrix multiplication of W
    xTwx <- t(x1) %*% U_inverse %*% x1 + t(x2) %*% IQ %*% x2  # Combined XtWX computation
    
    # Compute XtWy (X'Wy), the weighted response vector
    Qt_y <- qr.qty(qr_decomp, y)            # Transformed y using QR decomposition
    y1 <- Qt_y[1:p]                         # Top part of transformed y
    y2 <- Qt_y[(p + 1):n]                   # Bottom part of transformed y
    
    # Block matrix multiplication of W
    xTwy <- t(x1) %*% U_inverse %*% y1 + t(x2) %*% IQ %*% y2  # Combined XtWy computation
    
    # Compute beta_hat (estimated fixed effects coefficients)
    A <- xTwx                               # XtWX matrix
    B <- xTwy                               # XtWy vector
    beta_hat <- backsolve(chol(A), forwardsolve(t(chol(A)), B))  # A * β = B
    
    # Calculate part1: log(det(W^(-1))) and part2: the weighted sum of squared residuals
      # part1: log(det(W^(-1))
    log_likelihood_det <- 2 * sum(log(diag(chol(U)))) + (n - p) * log(sigma2)  # Log determinant
    
      # part2: the weighted sum of squared residuals
    c <- y - X %*% beta_hat                 # Residual vector
    Qt_c <- qr.qty(qr_decomp, c)            # Transformed residuals using QR decomposition
    c1 <- Qt_c[1:p]                         # Top part of transformed residuals
    c2 <- Qt_c[(p + 1):n]                   # Bottom part of transformed residuals
    
    # Block matrix multiplication of W for the weighted residual sum of squares 
    CTwC <- t(c1) %*% U_inverse %*% c1 + t(c2) %*% IQ %*% c2  # Weighted sum of squared residuals
    
    # Compute the negative log-likelihood as the sum of weighted residuals and log determinant
    negative_log_likelihood <- CTwC / 2 + log_likelihood_det / 2
  }
  
  # Store beta_hat as an attribute of the result for later use or inspection
  attr(negative_log_likelihood, "beta") <- beta_hat
  
  return(negative_log_likelihood)  # Return the negative log-likelihood
}


