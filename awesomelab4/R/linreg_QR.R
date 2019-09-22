#' Multiple Regression Model using QR decomposition
#'
#' Multiple Regression Model works by... #TODO#
#' 
#' @param formula an object of class "\code{formula}" describing the model sctructure. 
#' @param data an object of class "data frame" containing the variables in the model.
#'
#' @return \code{linreg} returns an object of class "\code{}"
#'
#' @examples
#' linreg(formula=Petal.Length~Species, data=iris)
#' 
#' @references #TODO#
#' @importFrom #TODO#
#'
#' @export
#'

linreg <- function(formula, data){
  
  ## 0. Check that the class of the formula argument is correct:
  stopifnot(class(formula) == "formula")
  
  ## 1. Initialization
  x <- model.matrix(formula, data) # X matrix (independent variables)
  name_dep <- all.vars(formula)[1] # Dependent variable/s name/s
  y <- data[, name_dep] # y (dependent variable/s)
  
  ## 2. Calculation of Q and R
  qr_x <- qr(x)
  Q <- qr.Q(qr_x) # orthogonal matrix
  R <- qr.R(qr_x) # triangular matrix
  
  ## 3. Estimations (Computations using ordinary least squares).
  # Regression coefficients:
  beta_hat <- as.vector(backsolve(R, (t(Q) %*% y)))
  
  # The fitted values:
  y_hat <- x %*% beta_hat
  
  # The residuals:
  e_hat <- y - (x %*% beta_hat)
  
  # The degrees of freedom:
  df <- nrow(data) - ncol(x) # number of observations - number of parameters in the model
  
  # The residual variance:
  sigma2_hat <- as.numeric((t(e_hat)%*%e_hat)/df)
  
  # The variance of the regression coefficients:
  var_hat <- diag(solve(t(R)%*%R) * sigma2_hat) # Var(beta_hat)= (R^T * R)^(-1) * sigma_hat^2
  std_error <- t(sqrt(var_hat))    # Std. Errors
  
  # The t-values for each coefficient:
  t_values <- beta_hat / std_error
  
  # The p-values for each coefficient:
  p_values <- 2*abs(pt(t_values, df, log.p = T))
  
  return(cbind(beta_hat,var_hat))
}