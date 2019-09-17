linreg <- setRefClass ("linreg",
  fields = c (
    regression_coefficients = "ANY", 
    fitted_values = "ANY",
    residuals = "ANY",
    degrees_of_freedom = "ANY",
    residual_variance = "ANY",
    variance_of_coefficients = "ANY",
    t_values = "ANY",
    p_values = "ANY"),
  methods = c (
    print = function() {
      print("hey")
    }
  )
)

customer <- setRefClass ("customer",
   fields = c(
     money = "ANY"),
   methods = c(
     add_funds = function(amount){
       money <<- money + amount
     }
   )
)