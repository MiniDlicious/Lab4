#' Multiple Regression Model
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
#' @importFrom stats median model.matrix pt quantile
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
  
  ## 2. Estimations (Computations using ordinary least squares).
  # Regression coefficients:
  xtx <- t(x) %*% x
  beta_hat <- solve(xtx) %*% t(x) %*% y
  
  # The fitted values:
  y_hat <- x %*% beta_hat
  
  # The residuals:
  e_hat <- y - y_hat
  
  # The degrees of freedom:
  n <- nrow(data) # number of observations
  p <- ncol(x)    # number of parameters in the model
  df <- n-p
  
  # The residual variance:
  sigma2_hat <- (t(e_hat)%*%e_hat)/df
  RSE <- sqrt(sigma2_hat)    # Residual standard error
  
  # The variance of the regression coefficients:
  var_hat <- sigma2_hat %*% diag(solve(xtx))
  std_error <- t(sqrt(var_hat))    # Std. Errors
  
  # The t-values for each coefficient:
  t_values <- beta_hat / std_error
  
  # The p-values for each coefficient:
  p_values <- 2*abs(pt(t_values, df, log.p = T))
  p_values
  
  ## 3. Output:
  residuals <- c(min(e_hat), quantile(e_hat,0.25), median(e_hat),quantile(e_hat,0.75),max(e_hat))
  names(residuals) <- c("Min","1Q","Median","3Q","Max")
  
  coef <- cbind(beta_hat, std_error, t_values, p_values)
  colnames(coef) <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
  
  #text <- cat("Residual standard error:", RSE, "on", df , "degrees of freedom")
  
  model <- list()
  class(model) <- "lm"
  model$coefficients <- coef
  model$call <- paste("lm(formula=", name_dep, "~", all.vars(formula)[-1] , ", data=", "data" , ")")
  return(model)
}
