#' A Reference Class.
#'
#' @field formula TODO
#' @field data TODO
#' @field regression_coefficients TODO
#' @field fitted_values TODO
#' @field residuals TODO
#' @field degrees_of_freedom TODO
#' @field residual_variance TODO
#' @field variance_of_coefficients TODO
#' @field std_error TODO
#' @field t_values TODO
#' @field p_values TODO
#' 
#' @examples
#' data("iris")
#' example <- linreg$new(formula=Petal.Length~Species, data=iris)
#' 
#' @references #TODO#
#' @importFrom methods new
#'
#' @export linreg
#' @exportClass linreg
#'

linreg <- setRefClass ("linreg",
  fields = c (
    formula = "formula",
    data = "ANY",
    regression_coefficients = "ANY", 
    fitted_values = "ANY",
    residuals = "ANY",
    degrees_of_freedom = "ANY",
    residual_variance = "ANY",
    variance_of_coefficients = "ANY",
    std_error = "ANY",
    t_values = "ANY",
    p_values = "ANY"),
  methods = c (
    initialize = function(formula, data){
      "Initialize function calculates all values needed from formula and data."
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
      regression_coefficients <<- as.vector(backsolve(R, (t(Q) %*% y)))
      
      # The fitted values:
      fitted_values <<- x %*% regression_coefficients
      
      # The residuals:
      residuals <<- y - (x %*% regression_coefficients)
      
      # The degrees of freedom:
      degrees_of_freedom <<- nrow(data) - ncol(x) # number of observations - number of parameters in the model
      
      # The residual variance:
      residual_variance <<- as.numeric((t(residuals)%*%residuals)/degrees_of_freedom)
      
      # The variance of the regression coefficients:
      variance_of_coefficients <<- diag(solve(t(R)%*%R) * residual_variance) # Var(beta_hat)= (R^T * R)^(-1) * sigma_hat^2
      std_error <<- t(sqrt(variance_of_coefficients))    # Std. Errors
      
      # The t-values for each coefficient:
      t_values <<- regression_coefficients / std_error
      
      # The p-values for each coefficient:
      p_values <<- 2*abs(pt(t_values, degrees_of_freedom, log.p = T))
      return(NULL)
    },
    show = function() {
      "Modifies the print() function for class."
      coeff <- matrix(c(names(std_error), regression_coefficients), ncol=3)
      colnames(coeff) <- colnames(t_values)
      return(coeff)
    },
    print = function() {
      "Same as show, but allows for linreg$print(). Will break print(linreg) that show allows, bad test case!"
      show()
    },
    resid = function(){
      "Returns the residuals."
      residual_list <- c(min(residuals), quantile(residuals,0.25), median(residuals),quantile(residuals,0.75),max(residuals))
      names(residual_list) <- c("Min","1Q","Median","3Q","Max")
      return(residual_list)
    },
    pred = function(){
      "Returns the predictions."
      return(fitted_values)
    },
    coef = function(){
      "Returns the coefficients."
      return(regression_coefficients)
    },
    summary = function() {
      "Prints a summary for the calculated data."
      coef <- cbind(regression_coefficients, std_error[1,], t_values[1,], p_values[1,])
      colnames(coef) <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
      print(coef)
    }
  )
)

# customer <- setRefClass ("customer",
#    fields = c(
#      money = "ANY"),
#    methods = c(
#      initialize = function(money){
#        money <<- money * 2
#      },
#      add_funds = function(amount){
#        money <<- money + amount
#      },
#      show = function(){
#        print("Money is not infinite")
#      }
#    )
# )