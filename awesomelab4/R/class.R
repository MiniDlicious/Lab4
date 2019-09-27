#' A Reference Class.
#'
#' \code{linreg()} calculates a number of statistics (see section \code{"Fields"}) by default given both arguments: formula and data. 
#' It also stores diverse functions implemented using RC object oriented programming. Find more information regarding these functions on section \code{"Methods"}.
#' 
#'
#' @field formula an object of class "formula".
#' @field data a data frame.
#' @field regression_coefficients the regression coefficients of the model calculated using the QR decomposition.
#' @field fitted_values the fitted values.
#' @field residuals the residuals. It is the difference between the real values and the fitted values.
#' @field degrees_of_freedom the degrees of freedom: the number of observations minus the number of parameters in the model.
#' @field residual_variance the residual variance. 
#' @field variance_of_coefficients the variance of the regression coefficients.
#' @field std_error the standard error of the regression coefficients. It is the square root of the variance.
#' @field t_values the t-values for each coefficient.
#' @field p_values the p-values for each coefficient.
#' 
#' @examples
#' data("iris")
#' example <- linreg$new(formula=Petal.Length~Species, data=iris)
#' 
#' @references Reference Classes: \url{http://adv-r.had.co.nz/R5.html}
#' 
#' QR decomposition: \url{https://en.wikipedia.org/wiki/QR_decomposition}
#' 
#' Vignettes: \url{http://r-pkgs.had.co.nz/vignettes.html}
#' 
#' @importFrom methods new
#' @importFrom ggplot2 ggplot
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
    std_residuals = "ANY",
    degrees_of_freedom = "ANY",
    residual_variance = "ANY",
    variance_of_coefficients = "ANY",
    std_error = "ANY",
    t_values = "ANY",
    p_values = "ANY",
    arguments = "ANY"),
  methods = c (
    initialize = function(formula, data){
      "Initialize function calculates all values needed from formula and data."
      ## 0. Check that the class of the formula argument is correct:
      stopifnot(class(formula) == "formula")
      #base::print(as.list(sys.calls()))
      # Extract arguments
      #arg_string <- as.character(as.list(sys.calls())[[1]]) # only works for user, not for testthat or R CMD Check
      arg_string <- paste(as.list(sys.calls()), collapse='')
      #formula_pattern <- "(?<=(linreg\$new\())([^,=]*)(?=\,.*\))|(?<=(linreg\$new*\((formula)=))([^,]*)(?=\,.*\))"
      #data_pattern <- "(?<=(linreg\\$new.*\\(.*)(data).*=?)(\\w+)(?=\\))|(?<=(linreg\\$new.*\\(.*),)(\\w+)(?=\\))"
      
      # Extract function call
      function_pattern <- "linreg\\$new\\(.*?\\)"
      function_call <- regmatches(arg_string, regexpr(function_pattern, arg_string, perl = TRUE))
      
      # From function call, extract data argument
      data_pattern <- "(\\w+)(?=\\))"
      data_argument <- regmatches(function_call, regexpr(data_pattern, function_call, perl = TRUE))
      
      # Save arguments
      arguments <<- c(formula, data_argument)
      
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
      std_residuals <<- residuals/(sqrt((1/(nrow(data)-1))*sum(residuals^2)))
      
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
      #p_values <<- 2*abs(pt(t_values, degrees_of_freedom, log.p = T))
      p_values <<- pt(-abs(t_values), degrees_of_freedom)
      return(NULL)
    },
    show = function() {
      "Modifies the print() function for class."
      base::print("Use ...$print()")
    },
    print = function() {
      "Same as show, but allows for linreg$print(). Will break print(linreg) that show allows, bad test case!"
      coeff <- matrix(c(names(std_error), regression_coefficients), ncol=3)
      colnames(coeff) <- colnames(t_values)
      
      model <- list()
      model$Call <- paste("linreg(formula = ", arguments[1], ", data = ", arguments[2], ")", sep="")
      model$Coefficients <- coeff
      base::print(model)
    },
    plot = function(){
      liu_theme <- theme_bw() + theme(
        plot.title = element_text(color = "#00b9e7"),
        panel.background = element_rect(fill = "#c7e6e3", colour = "#17c7d2",
                                        size = 2, linetype = "solid"),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                        colour = "white"), 
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                        colour = "white"),
        #scale_color_manual(values=c('#ff6442','#8781d3', '#fcf05f', '#687f91')),
        legend.title = element_text(face = "bold"),
        legend.position = "none",
        axis.text.x = element_text(color="#687f91"),
        axis.text.y = element_text(color="#687f91")
        )
      
      df <- data.frame(residuals, std_residuals, fitted_values)
      plot1 <- ggplot(df, aes(y=residuals, x=fitted_values)) + geom_point() + geom_line(aes(y = mean(residuals), x= fitted_values, colour="red")) + ggtitle("Residuals vs Fitted") + xlab("Fitted values") + ylab("Residuals") + liu_theme
      plot2 <- ggplot(df, aes(y=sqrt(abs(std_residuals)), x=fitted_values)) + geom_point() + geom_line(aes(y = mean(sqrt(abs(std_residuals))), x= fitted_values, colour="red")) + ggtitle("Scale-Location") + xlab("Fitted values") + ylab(expression(sqrt("|Standardized residuals|"))) + liu_theme
      plot_grid(plot1, plot2, labels = "AUTO")
      },
    resid = function(){
      "Returns the residuals."
      return(residuals)
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
      .categorise_significance <- function(x){
        if(x < 0.001){
          return("***")
        }else if(x < 0.01){
          return("**")
        }else if(x < 0.05){
          return("*")
        }else if(x < 0.1){
          return(".")
        }else{
          return(" ")
        }
      }
      significance <- lapply(as.vector(p_values), .categorise_significance)
      coef <- cbind(regression_coefficients, std_error[1,], t_values[1,], p_values[1,], significance)
      colnames(coef) <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)", "Significance")
      base::print(coef)
      
      base::print("Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1")
      rse <- sqrt(sum(residuals**2)/degrees_of_freedom)
      base::print(paste("Residual standard error: ", rse, " on ", degrees_of_freedom, " degrees of freedom", sep=""))
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